"""File Sorter plugin for AMPA.

Recursively scans a user-selected folder for FITS files and presents
them in one tab per instrument (the ``INSTRUME`` header, normalised via
aliases - e.g. raw ``ALFOSC_FASU`` groups under **ALFOSC**). Each tab
is a table with one row per file: common columns (``TCSTGT``, ``OBJECT``,
``IMAGETYP``, ``FILTER``, ``OBS_MODE``, ``EXPTIME``, ``DATE-OBS``) plus
the instrument's dedicated configuration columns (ALFOSC: ``FAFLTNM``,
``FBFLTNM``, ``ALGRNM``). Double-click a row (or use **Load Selected**)
to open the file in the AMPA viewer.

Files missing headers land under placeholder tabs (``(no TCSTGT)`` /
``(no INSTRUME)``); files whose header cannot be read (corrupt,
truncated, ...) land on an ``(unreadable)`` tab with an ``Error``
column. The scan itself runs on a background task with a progress bar
and cancel support, so huge data folders never freeze the GUI.

The plugin is distributed via the ``not-tools-ampa-plugin`` git
repository: clone it and add the folder as a Local Plugin Directory, or
install it from the registry manifest (``plugins.json``) in the repo root.
"""

from __future__ import annotations

import os
import warnings
from typing import Any, Dict, List, Optional

from astropy.io import fits
from PySide6 import QtCore, QtWidgets

from ampa.core.apis import settings_api, ui_api
from ampa.core.basemodule import BaseModule
from ampa.core.logging import log

__LOGMODULE__ = "FileSorter"

_SETTINGS_GROUP = "FileSorter"
_KEY_LAST_FOLDER = f"{_SETTINGS_GROUP}/last_folder"

# Recognised FITS file suffixes (lower-case). gz-compressed variants are
# handled transparently by astropy.
_FITS_SUFFIXES = (".fits", ".fit", ".fts", ".fits.gz", ".fit.gz", ".fts.gz")

# Placeholder group labels for missing headers / unreadable files.
_NO_TARGET = "(no TCSTGT)"
_NO_INSTRUMENT = "(no INSTRUME)"
_UNREADABLE = "(unreadable)"

# Shown in table cells when a header key is absent.
_MISSING = "\u2014"  # em dash

# Raw INSTRUME value -> display/tab name. Files from the same
# instrument arrive under different INSTRUME spellings; the alias
# normalises them into one tab.
_INSTRUMENT_ALIASES = {
    "ALFOSC_FASU": "ALFOSC",
}

# Instrument-specific header keys that become extra table columns on
# that instrument's tab (matched on the display name). Extend as
# instruments get dedicated views (e.g. FIES).
_INSTRUMENT_KEYS = {
    "ALFOSC": ("FAFLTNM", "FBFLTNM", "ALGRNM"),
}

# Header keys shown as common columns on every tab (after the File and
# TCSTGT columns): observation identifiers first, then context and
# observation-mode keys.
_COMMON_HEADER_KEYS = ("GROUPID", "BLOCKID", "SEQID",
                       "OBJECT", "IMAGETYP", "FILTER",
                       "OBS_MODE", "EXPTIME", "DATE-OBS")


# ======================================================================
# Pure helpers (module level so they can be tested without the GUI)
# ======================================================================

def find_fits_files(root: str) -> List[str]:
    """Recursively collect FITS files under *root*.

    Returns absolute, deterministic paths (directory walk and file names
    are sorted).
    """
    matches: List[str] = []
    for dirpath, dirnames, filenames in os.walk(root):
        dirnames.sort()
        for name in sorted(filenames):
            if name.lower().endswith(_FITS_SUFFIXES):
                matches.append(os.path.abspath(os.path.join(dirpath, name)))
    return matches


def read_grouping_headers(path: str):
    """Read the primary header of a FITS file.

    Returns ``(tcstgt, instrume, header)`` with empty strings for missing
    grouping keywords. Raises on unreadable files - callers handle that.
    """
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        header = fits.getheader(path, 0)
    tcstgt = str(header.get("TCSTGT", "") or "").strip()
    instrume = str(header.get("INSTRUME", "") or "").strip()
    return tcstgt, instrume, header


def scan_folder(root: str, handle=None) -> Dict[str, Any]:
    """Scan *root* recursively and read the grouping header of every FITS.

    ``handle`` is an optional :class:`~ampa.core.tasks.TaskHandle` used
    for progress reporting and cooperative cancellation. Returns a dict::

        {
            "root": <scanned folder>,
            "records": [{"path", "tcstgt", "instrume", "instrume_raw",
                         "header", "error"}],
        }

    Each record describes one file; ``header`` is the astropy header (or
    ``None`` for unreadable files, with ``error`` carrying the reason).
    ``instrume`` carries the display name (alias-normalised, see
    ``_INSTRUMENT_ALIASES``) used for grouping; ``instrume_raw`` the
    original header value.
    """
    files = find_fits_files(root)
    total = len(files)
    records: List[Dict[str, Any]] = []
    for index, path in enumerate(files):
        if handle is not None:
            handle.check_cancelled()
            if total:
                handle.report(
                    f"Reading headers {index + 1}/{total}",
                    percent=100.0 * (index + 1) / total,
                )
        try:
            tcstgt, instrume, header = read_grouping_headers(path)
        except Exception as exc:
            log(__name__, __LOGMODULE__, "debug",
                f"Could not read header of {path}: {exc}")
            records.append({
                "path": path,
                "tcstgt": _UNREADABLE,
                "instrume": _UNREADABLE,
                "instrume_raw": _UNREADABLE,
                "header": None,
                "error": str(exc),
            })
            continue
        display = _INSTRUMENT_ALIASES.get(instrume.upper(), instrume)
        records.append({
            "path": path,
            "tcstgt": tcstgt or _NO_TARGET,
            "instrume": display or _NO_INSTRUMENT,
            "instrume_raw": instrume,
            "header": header,
            "error": None,
        })
    return {"root": root, "records": records}


def group_by_instrument(records: List[Dict[str, Any]]
                        ) -> Dict[str, List[Dict[str, Any]]]:
    """Group scan records into ``{instrument display name: [record, ...]}``.

    The order within each group follows the (already deterministic)
    record order.
    """
    groups: Dict[str, List[Dict[str, Any]]] = {}
    for record in records:
        groups.setdefault(record["instrume"], []).append(record)
    return groups


def format_value(value) -> str:
    """Generic cell text: trimmed value, ``_MISSING`` when absent."""
    text = str(value or "").strip()
    return text or _MISSING


def format_exptime(value) -> str:
    """EXPTIME cell text: compact number (``300`` for 300.0)."""
    if value in (None, ""):
        return _MISSING
    try:
        num = float(value)
    except (TypeError, ValueError):
        return str(value).strip() or _MISSING
    if num.is_integer():
        return str(int(num))
    return str(num)


def format_date_obs(value) -> str:
    """DATE-OBS cell text: ISO trimmed to ``YYYY-MM-DD HH:MM``."""
    text = str(value or "").strip()
    if not text:
        return _MISSING
    if len(text) >= 16 and text[10] in ("T", " "):
        return f"{text[:10]} {text[11:16]}"
    return text


_KEY_FORMATTERS = {
    "OBS_MODE": format_value,
    "EXPTIME": format_exptime,
    "DATE-OBS": format_date_obs,
}


def format_key(header, key: str) -> str:
    """Format one header-key cell using the key's dedicated formatter."""
    if header is None or key not in header:
        return _MISSING
    return _KEY_FORMATTERS.get(key, format_value)(header[key])


def table_headers(instrument: str) -> List[str]:
    """Column headers for one instrument tab."""
    headers = ["File", "TCSTGT", *_COMMON_HEADER_KEYS]
    headers.extend(_INSTRUMENT_KEYS.get(instrument, ()))
    if instrument == _UNREADABLE:
        headers.append("Error")
    return headers


def table_row_values(record: Dict[str, Any], instrument: str) -> List[str]:
    """Cell texts for one record, aligned with :func:`table_headers`."""
    header = record.get("header")
    row = [
        os.path.basename(record["path"]),
        _MISSING if record["tcstgt"] in (_NO_TARGET, _UNREADABLE)
        else record["tcstgt"],
    ]
    for key in _COMMON_HEADER_KEYS:
        row.append(format_key(header, key))
    for key in _INSTRUMENT_KEYS.get(instrument, ()):
        row.append(format_key(header, key))
    if instrument == _UNREADABLE:
        row.append(format_value(record.get("error")))
    return row


# ======================================================================
# Plugin
# ======================================================================

class _NumericItem(QtWidgets.QTableWidgetItem):
    """Table item that sorts numerically (used for the EXPTIME column).

    Qt compares table items by display text, which orders ``"300"``
    before ``"90"``; this subclass compares the underlying float instead.
    Missing values sort last (``+inf``).
    """

    def __init__(self, text: str, number: float):
        super().__init__(text)
        self._number = number

    def __lt__(self, other):
        if isinstance(other, _NumericItem):
            return self._number < other._number
        return super().__lt__(other)


class FileSorterPlugin(BaseModule):
    """Plugins ▸ NOT Toolkit ▸ File Sorter."""

    def __init__(self):
        super().__init__(
            title="File Sorter",
            category="Plugins",
            section="NOT Toolkit",
        )
        settings_api.define_setting(
            group=_SETTINGS_GROUP,
            key="last_folder",
            default_value="",
            value_type=str,
            description="Last folder scanned by the File Sorter.",
        )
        self._scan_root: Optional[str] = None
        self._records: List[Dict[str, Any]] = []

    # -- lifecycle --------------------------------------------------------

    def show_gui(self):
        """Create the window lazily, then show/raise/focus it.

        Overridden because the base implementation resets
        ``gui_widget`` to ``None`` after ``on_activated()`` returns,
        which would discard a window built there.
        """
        if self.gui_widget is None:
            self._build_window()
            last = settings_api.get_value(_KEY_LAST_FOLDER, "")
            if last and os.path.isdir(last):
                self.folder_input.setText(last)
                self._scan_root = last
        self.gui_widget.show()
        self.gui_widget.raise_()
        self.gui_widget.activateWindow()

    # -- GUI creation -----------------------------------------------------

    def _build_window(self):
        window = QtWidgets.QWidget()
        window.setWindowTitle("File Sorter")
        layout = QtWidgets.QVBoxLayout(window)

        # Folder selection row
        folder_row = QtWidgets.QHBoxLayout()
        folder_row.addWidget(QtWidgets.QLabel("Folder:"))
        self.folder_input = QtWidgets.QLineEdit()
        self.folder_input.setReadOnly(True)
        self.folder_input.setPlaceholderText("No folder selected")
        folder_row.addWidget(self.folder_input, 1)
        self.browse_button = QtWidgets.QPushButton("Browse…")
        self.browse_button.clicked.connect(self._choose_folder)
        folder_row.addWidget(self.browse_button)
        self.rescan_button = QtWidgets.QPushButton("Rescan")
        self.rescan_button.clicked.connect(self._start_scan)
        self.rescan_button.setEnabled(False)
        folder_row.addWidget(self.rescan_button)
        layout.addLayout(folder_row)

        # One tab per instrument; each tab is a table with one row per
        # file and that instrument's columns.
        self.tabs = QtWidgets.QTabWidget()
        self.tabs.setDocumentMode(True)
        layout.addWidget(self.tabs, 1)

        # Bottom row: load button + status label
        bottom_row = QtWidgets.QHBoxLayout()
        self.load_button = QtWidgets.QPushButton("Load Selected")
        self.load_button.setToolTip("Open the selected FITS file in the AMPA viewer")
        self.load_button.setEnabled(False)
        self.load_button.clicked.connect(self._load_selected)
        bottom_row.addWidget(self.load_button)
        self.status_label = QtWidgets.QLabel("No folder scanned yet.")
        bottom_row.addWidget(self.status_label, 1)
        layout.addLayout(bottom_row)

        window.resize(1000, 620)
        self.gui_widget = window

    # -- folder selection / scanning ---------------------------------------

    def _choose_folder(self):
        folder = QtWidgets.QFileDialog.getExistingDirectory(
            self.gui_widget,
            "Select folder to scan for FITS files",
            self.folder_input.text() or "",
            QtWidgets.QFileDialog.ShowDirsOnly
            | QtWidgets.QFileDialog.DontResolveSymlinks,
        )
        if not folder:
            return
        self.folder_input.setText(folder)
        settings_api.set_value(_KEY_LAST_FOLDER, folder)
        self._start_scan()

    def _start_scan(self) -> bool:
        folder = self.folder_input.text().strip()
        if not folder or not os.path.isdir(folder):
            ui_api.show_error_dialog(
                "File Sorter", "Please select a valid folder first.")
            return False
        self._scan_root = folder
        self.tabs.clear()
        self.status_label.setText(f"Scanning {folder} …")
        self.rescan_button.setEnabled(False)
        started = self.run_task(
            self._scan_work,
            on_done=self._on_scan_done,
            on_error=self._on_scan_error,
            on_cancel=self._on_scan_cancelled,
            title="Scanning FITS headers",
            message="Collecting FITS files…",
            cancelable=True,
            widgets=(self.browse_button, self.rescan_button,
                     self.load_button),
        )
        if not started:
            self.rescan_button.setEnabled(True)
            self.status_label.setText("A scan is already running.")
        return started

    def _scan_work(self, handle):
        return scan_folder(self._scan_root, handle=handle)

    def _on_scan_done(self, result):
        self._records = result["records"]
        self._populate_tabs(self._records)
        self.rescan_button.setEnabled(True)
        groups = group_by_instrument(self._records)
        unreadable = len(groups.get(_UNREADABLE, ()))
        msg = (f"{len(groups) - (1 if unreadable else 0)} instrument(s) · "
               f"{len(self._records)} file(s)")
        if unreadable:
            msg += f" · {unreadable} unreadable"
        self.status_label.setText(msg)
        if not self._records:
            self.status_label.setText("No FITS files found in this folder.")
        log(__name__, __LOGMODULE__, "info",
            f"Scan of {result['root']}: {msg}")

    def _on_scan_error(self, message):
        self.rescan_button.setEnabled(True)
        self.status_label.setText("Scan failed.")
        ui_api.show_error_dialog("File Sorter",
                                 f"Scanning failed:\n{message}")

    def _on_scan_cancelled(self):
        self.rescan_button.setEnabled(True)
        self.status_label.setText("Scan cancelled.")

    # -- tab / table population ---------------------------------------------

    def _populate_tabs(self, records):
        """Rebuild the instrument tabs from scan records."""
        self.tabs.clear()
        self.load_button.setEnabled(False)
        groups = group_by_instrument(records)
        for instrument in sorted(groups, key=str.lower):
            table = self._build_table(instrument, groups[instrument])
            self.tabs.addTab(table,
                             f"{instrument} ({len(groups[instrument])})")

    def _build_table(self, instrument: str,
                    records: List[Dict[str, Any]]) -> QtWidgets.QTableWidget:
        headers = table_headers(instrument)
        try:
            exptime_column = headers.index("EXPTIME")
        except ValueError:
            exptime_column = -1

        table = QtWidgets.QTableWidget(len(records), len(headers))
        table.setHorizontalHeaderLabels(headers)
        table.verticalHeader().setVisible(False)
        table.setSelectionBehavior(
            QtWidgets.QAbstractItemView.SelectionBehavior.SelectRows)
        table.setSelectionMode(
            QtWidgets.QAbstractItemView.SelectionMode.SingleSelection)
        table.setEditTriggers(
            QtWidgets.QAbstractItemView.EditTrigger.NoEditTriggers)
        table.setAlternatingRowColors(True)
        header_view = table.horizontalHeader()
        header_view.setSectionResizeMode(
            QtWidgets.QHeaderView.ResizeMode.ResizeToContents)
        header_view.setSectionResizeMode(
            0, QtWidgets.QHeaderView.ResizeMode.Stretch)

        for row, record in enumerate(sorted(records,
                                             key=lambda r: r["path"].lower())):
            for column, text in enumerate(table_row_values(record, instrument)):
                if column == 0:
                    item = QtWidgets.QTableWidgetItem(text)
                    item.setData(QtCore.Qt.ItemDataRole.UserRole, record)
                    item.setToolTip(record["path"])
                elif column == exptime_column:
                    # Numeric item so the EXPTIME column sorts
                    # numerically (text sort would order "9" after "10").
                    fits_header = record.get("header") or {}
                    try:
                        number = float(fits_header["EXPTIME"])
                    except (KeyError, TypeError, ValueError):
                        number = float("inf")
                    item = _NumericItem(text, number)
                else:
                    item = QtWidgets.QTableWidgetItem(text)
                table.setItem(row, column, item)

        # Qt's default sort indicator is *descending*; pin an explicit
        # ascending order so enabling sorting doesn't reverse the rows.
        table.horizontalHeader().setSortIndicator(
            0, QtCore.Qt.SortOrder.AscendingOrder)
        table.setSortingEnabled(True)
        table.itemDoubleClicked.connect(self._on_item_double_clicked)
        table.itemSelectionChanged.connect(self._update_load_button)
        return table

    # -- selection / loading --------------------------------------------------

    def _current_record(self) -> Optional[Dict[str, Any]]:
        """The record on the currently selected row of the current tab."""
        table = self.tabs.currentWidget()
        if not isinstance(table, QtWidgets.QTableWidget):
            return None
        item = table.item(table.currentRow(), 0)
        if item is None:
            return None
        record = item.data(QtCore.Qt.ItemDataRole.UserRole)
        return record if isinstance(record, dict) else None

    def _update_load_button(self):
        self.load_button.setEnabled(self._current_record() is not None)

    def _on_item_double_clicked(self, item):
        record = item.data(QtCore.Qt.ItemDataRole.UserRole)
        if isinstance(record, dict):
            self._load_record(record)

    def _load_selected(self):
        record = self._current_record()
        if record is not None:
            self._load_record(record)

    def _load_record(self, record: Dict[str, Any]):
        path = record["path"]
        if record.get("header") is None:
            if not ui_api.confirm_dialog(
                    "File Sorter",
                    "This file could not be read during the scan.\n"
                    f"Try to open it anyway?\n\n{path}"):
                return
        if ui_api.load_fits_file(path):
            ui_api.write_to_statusbar(
                f"Loaded {os.path.basename(path)}", timeout=4000)
        else:
            ui_api.show_error_dialog(
                "File Sorter", f"Could not load:\n{path}")
