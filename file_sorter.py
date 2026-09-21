"""File Sorter plugin for AMPA.

Recursively scans a user-selected folder for FITS files and presents
them in one tab per instrument (the ``INSTRUME`` header, normalised via
aliases - e.g. raw ``ALFOSC_FASU`` groups under **ALFOSC**). Files that
are not science observations (``IMAGECAT`` other than ``SCIENCE``) or
that carry no target (``TCSTGT``) get their own ``<INSTRUMENT> CALIB``
tab. Each tab is a table with one row per file: common columns
(``TARGET``, ``OBJECT``, ``IMAGETYPE``, ``OBSMODE``, ``EXPTIME``,
``DATE-OBS``, ``TELALT`` (1 decimal), ``AIRMASS`` (2 decimals)) plus
the instrument's dedicated configuration columns (ALFOSC: ``FASU A``,
``FASU B``, ``GRISM``; FIES: ``FIBER``). Column headers use friendlier
labels than the raw FITS keywords where a mapping exists (see
``_KEY_LABELS``). Double-click a row (or use **Load Selected**) to open
the file in the AMPA viewer; select a range of rows and use **New
Sequence From Selected ...** (or the table's context menu) to create a
session-only AMPA sequence from them. The search bar above the tabs
filters every table across all columns (case-insensitive); tabs that
hide rows get a ``*`` appended to their name.

Files missing headers land under placeholder tabs (``(no TARGET)`` /
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
from PySide6 import QtCore, QtGui, QtWidgets

from ampa.core.apis import (Sequence, sequence_api, settings_api,
                            ui_api)
from ampa.core.basemodule import BaseModule
from ampa.core.logging import log

__LOGMODULE__ = "FileSorter"

_SETTINGS_GROUP = "FileSorter"
_KEY_LAST_FOLDER = f"{_SETTINGS_GROUP}/last_folder"

# Recognised FITS file suffixes (lower-case). gz-compressed variants are
# handled transparently by astropy.
_FITS_SUFFIXES = (".fits", ".fit", ".fts", ".fits.gz", ".fit.gz", ".fts.gz")

# Placeholder group labels for missing headers / unreadable files.
_NO_TARGET = "(no TARGET)"
_NO_INSTRUMENT = "(no INSTRUME)"
_UNREADABLE = "(unreadable)"

# Shown in table cells when a header key is absent.
_MISSING = "-"

# Raw INSTRUME value -> display/tab name. Files from the same
# instrument arrive under different INSTRUME spellings; the alias
# normalises them into one tab.
_INSTRUMENT_ALIASES = {
    "ALFOSC_FASU": "ALFOSC",
}

# Instrument-specific header keys that become extra table columns on
# that instrument's tab (matched on the display name; CALIB tabs reuse
# their instrument's columns). Extend as instruments get dedicated
# views.
_INSTRUMENT_KEYS = {
    "ALFOSC": ("FAFLTNM", "FBFLTNM", "ALGRNM"),
    "FIES": ("FIFMSKNM",),
}

# Friendlier display labels for FITS keywords shown as table columns.
# Lookup always uses the raw keyword; only the visible header text
# changes.
_KEY_LABELS = {
    "TCSTGT": "TARGET",
    "IMAGETYP": "IMAGETYPE",
    "OBS_MODE": "OBSMODE",
    "FAFLTNM": "FASU A",
    "FBFLTNM": "FASU B",
    "ALGRNM": "GRISM",
    "FIFMSKNM": "FIBER",
}

# Header keys shown as common columns on every tab (after the File and
# TCSTGT columns): observation identifiers first, then context and
# observation-mode keys.
_COMMON_HEADER_KEYS = ("GROUPID", "BLOCKID", "SEQID",
                        "OBJECT", "IMAGETYP",
                        "OBS_MODE", "EXPTIME", "DATE-OBS",
                        "TELALT", "AIRMASS")


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
                         "is_calib", "header", "error"}],
        }

    Each record describes one file; ``header`` is the astropy header (or
    ``None`` for unreadable files, with ``error`` carrying the reason).
    ``instrume`` carries the display name (alias-normalised, see
    ``_INSTRUMENT_ALIASES``) used for grouping; ``instrume_raw`` the
    original header value. ``is_calib`` flags non-science files
    (``IMAGECAT`` other than ``SCIENCE``, case-insensitive) and files
    without a ``TCSTGT``; both are shown on a ``<INSTRUMENT> CALIB``
    tab.
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
                "is_calib": False,
                "header": None,
                "error": str(exc),
            })
            continue
        display = _INSTRUMENT_ALIASES.get(instrume.upper(), instrume)
        imagecat = str(header.get("IMAGECAT", "") or "").strip().upper()
        records.append({
            "path": path,
            "tcstgt": tcstgt or _NO_TARGET,
            "instrume": display or _NO_INSTRUMENT,
            "instrume_raw": instrume,
            "is_calib": imagecat != "SCIENCE" or not tcstgt,
            "header": header,
            "error": None,
        })
    return {"root": root, "records": records}


def group_by_instrument(records: List[Dict[str, Any]]
                        ) -> Dict[str, List[Dict[str, Any]]]:
    """Group scan records into ``{tab name: [record, ...]}``.

    The order within each group follows the (already deterministic)
    record order. Non-science or targetless files (``is_calib``) are
    separated into a ``<INSTRUMENT> CALIB`` tab instead of the plain
    instrument tab.
    """
    groups: Dict[str, List[Dict[str, Any]]] = {}
    for record in records:
        name = record["instrume"]
        if record.get("is_calib"):
            name = f"{name} CALIB"
        groups.setdefault(name, []).append(record)
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


def _format_number(value, decimals: int) -> str:
    """Numeric cell text rounded to *decimals*; ``_MISSING`` when absent.

    Non-numeric values pass through as trimmed text.
    """
    if value in (None, ""):
        return _MISSING
    try:
        num = float(value)
    except (TypeError, ValueError):
        return str(value).strip() or _MISSING
    return f"{num:.{decimals}f}"


def format_airmass(value) -> str:
    """AIRMASS cell text: number rounded to 2 decimals."""
    return _format_number(value, 2)


def format_telalt(value) -> str:
    """TELALT cell text: number rounded to 1 decimal."""
    return _format_number(value, 1)


# Dedicated formatters per header key (everything else falls back to
# format_value()).
_KEY_FORMATTERS = {
    "OBS_MODE": format_value,
    "EXPTIME": format_exptime,
    "DATE-OBS": format_date_obs,
    "TELALT": format_telalt,
    "AIRMASS": format_airmass,
}

# Row-background hues cycled through per GROUPID within a tab so
# adjacent groups are easy to tell apart. The saturation/lightness are
# derived from the current application palette at call time (see
# :func:`group_palette`): dark themes get deep muted shades that keep
# light text readable, light themes get soft pastels.
_GROUP_HUES = (210, 130, 55, 0, 25, 170, 265, 315)


def group_palette() -> List[QtGui.QColor]:
    """Shading colors adapted to the active AMPA color theme.

    Returns one color per entry of ``_GROUP_HUES``: muted dark tints
    when the application palette is dark (so light text stays
    readable) and soft pastels when it is light.
    """
    app = QtWidgets.QApplication.instance()
    base = (app.palette() if app is not None else QtGui.QPalette()).color(
        QtGui.QPalette.ColorRole.Window)
    dark = base.lightness() < 128
    saturation = 110
    lightness = 95 if dark else 235
    return [QtGui.QColor.fromHsl(hue, saturation, lightness)
            for hue in _GROUP_HUES]


def format_key(header, key: str) -> str:
    """Format one header-key cell using the key's dedicated formatter."""
    if header is None or key not in header:
        return _MISSING
    return _KEY_FORMATTERS.get(key, format_value)(header[key])


def base_instrument(tab_name: str) -> str:
    """The plain instrument name behind a tab/group name.

    Strips the `` CALIB`` suffix from ``<INSTRUMENT> CALIB`` tabs so
    their tables keep the instrument's dedicated columns; other names
    pass through unchanged.
    """
    if tab_name.endswith(" CALIB"):
        return tab_name[:-len(" CALIB")]
    return tab_name


def table_headers(instrument: str) -> List[str]:
    """Column headers for one tab, with friendly labels where mapped."""
    headers = ["File", _KEY_LABELS["TCSTGT"]]
    headers.extend(_KEY_LABELS.get(key, key) for key in _COMMON_HEADER_KEYS)
    headers.extend(_KEY_LABELS.get(key, key)
                   for key in _INSTRUMENT_KEYS.get(
                       base_instrument(instrument), ()))
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
    for key in _INSTRUMENT_KEYS.get(base_instrument(instrument), ()):
        row.append(format_key(header, key))
    if instrument == _UNREADABLE:
        row.append(format_value(record.get("error")))
    return row


def group_id_value(record: Dict[str, Any]) -> str:
    """The record's GROUPID (``""`` when absent or unreadable)."""
    header = record.get("header")
    if header is None:
        return ""
    return str(header.get("GROUPID", "") or "").strip()


def assign_group_colors(records: List[Dict[str, Any]]) -> Dict[str, Any]:
    """Map every GROUPID in *records* to a background color.

    Colors cycle through the palette in order of first appearance, so
    adjacent groups always differ; the same GROUPID maps to the same
    color. The empty GROUPID (missing/unreadable) maps to ``None`` - no
    shading.
    """
    colors: Dict[str, Any] = {}
    palette = group_palette()
    assigned = 0
    for record in records:
        gid = group_id_value(record)
        if not gid or gid in colors:
            continue
        colors[gid] = palette[assigned % len(palette)]
        assigned += 1
    colors[""] = None
    return colors


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
                # A remembered folder is already shown - scan it right
                # away instead of leaving the user with a dead Rescan
                # button (it only enables after a scan).
                self._start_scan()
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
        self.browse_button = QtWidgets.QPushButton("Browse...")
        self.browse_button.clicked.connect(self._choose_folder)
        folder_row.addWidget(self.browse_button)
        self.rescan_button = QtWidgets.QPushButton("Rescan")
        self.rescan_button.clicked.connect(self._start_scan)
        self.rescan_button.setEnabled(False)
        folder_row.addWidget(self.rescan_button)
        layout.addLayout(folder_row)

        # Search bar above the tabs: filters every table at once.
        self.search_input = QtWidgets.QLineEdit()
        self.search_input.setPlaceholderText("Search all columns")
        self.search_input.setClearButtonEnabled(True)
        self.search_input.textChanged.connect(self._apply_search_filter)
        layout.addWidget(self.search_input)

        # One tab per instrument; each tab is a table with one row per
        # file and that instrument's columns.
        self.tabs = QtWidgets.QTabWidget()
        self.tabs.setDocumentMode(True)
        self.tabs.currentChanged.connect(self._update_action_buttons)
        layout.addWidget(self.tabs, 1)

        # Bottom row: sequence + load buttons, status label
        bottom_row = QtWidgets.QHBoxLayout()
        self.sequence_button = QtWidgets.QPushButton(
            "New Sequence From Selected ...")
        self.sequence_button.setToolTip(
            "Create a new AMPA sequence from the selected rows\n"
            "(session-only - save it via the Sequence Manager)")
        self.sequence_button.setEnabled(False)
        self.sequence_button.clicked.connect(self._new_sequence_from_selection)
        bottom_row.addWidget(self.sequence_button)
        self.load_button = QtWidgets.QPushButton("Load Selected")
        self.load_button.setToolTip("Open the selected FITS file in the AMPA viewer")
        self.load_button.setEnabled(False)
        self.load_button.clicked.connect(self._load_selected)
        bottom_row.addWidget(self.load_button)
        self.status_label = QtWidgets.QLabel("No folder scanned yet.")
        bottom_row.addWidget(self.status_label, 1)
        layout.addLayout(bottom_row)

        window.resize(1280, 680)
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
        self.status_label.setText(f"Scanning {folder}...")
        self.rescan_button.setEnabled(False)
        started = self.run_task(
            self._scan_work,
            on_done=self._on_scan_done,
            on_error=self._on_scan_error,
            on_cancel=self._on_scan_cancelled,
            title="Scanning FITS headers",
            message="Collecting FITS files...",
            cancelable=True,
            widgets=(self.browse_button, self.rescan_button,
                     self.load_button, self.sequence_button),
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
        instruments = {base_instrument(name) for name in groups}
        unreadable = len(groups.get(_UNREADABLE, ()))
        msg = (f"{len(instruments) - (1 if unreadable else 0)} instrument(s), "
               f"{len(self._records)} file(s)")
        if unreadable:
            msg += f", {unreadable} unreadable"
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
        self.sequence_button.setEnabled(False)
        groups = group_by_instrument(records)
        for instrument in sorted(groups, key=str.lower):
            table = self._build_table(instrument, groups[instrument])
            index = self.tabs.addTab(table, "")
            table.setProperty("tabName", instrument)
            self._refresh_tab_label(index)
        # A search may already be active (e.g. after a rescan) - re-apply it.
        self._apply_search_filter(self.search_input.text())

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
            QtWidgets.QAbstractItemView.SelectionMode.ExtendedSelection)
        table.setEditTriggers(
            QtWidgets.QAbstractItemView.EditTrigger.NoEditTriggers)
        table.setAlternatingRowColors(True)
        header_view = table.horizontalHeader()
        # Every column sizes to its content so long filenames are never
        # cut off; when the window is narrower than the table a horizontal
        # scrollbar appears instead of squeezing columns. The last column
        # stretches to absorb any spare window width.
        header_view.setSectionResizeMode(
            QtWidgets.QHeaderView.ResizeMode.ResizeToContents)
        header_view.setSectionResizeMode(
            headers.index("DATE-OBS"),
            QtWidgets.QHeaderView.ResizeMode.Stretch)

        ordered = sorted(records, key=lambda r: r["path"].lower())
        group_colors = assign_group_colors(ordered)
        if any(color is not None for color in group_colors.values()):
            # Row shading by GROUPID replaces the alternating stripes.
            table.setAlternatingRowColors(False)

        for row, record in enumerate(ordered):
            color = group_colors.get(group_id_value(record))
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
                if color is not None:
                    item.setBackground(color)
                table.setItem(row, column, item)

        # Qt's default sort indicator is *descending*; pin an explicit
        # ascending order so enabling sorting doesn't reverse the rows.
        table.horizontalHeader().setSortIndicator(
            0, QtCore.Qt.SortOrder.AscendingOrder)
        table.setSortingEnabled(True)
        table.setContextMenuPolicy(
            QtCore.Qt.ContextMenuPolicy.CustomContextMenu)
        table.customContextMenuRequested.connect(
            lambda pos, t=table: self._show_table_context_menu(t, pos))
        table.itemDoubleClicked.connect(self._on_item_double_clicked)
        table.itemSelectionChanged.connect(self._update_action_buttons)
        return table

    # -- search / filtering ---------------------------------------------

    def _apply_search_filter(self, text: str):
        """Filter every table to rows matching *text* across all columns.

        The match is a case-insensitive substring over every visible
        cell of a row. Tabs that hide rows get a ``*`` appended to
        their name and show the visible count; an empty search shows
        everything.
        """
        needle = text.strip().lower()
        for index in range(self.tabs.count()):
            table = self.tabs.widget(index)
            if not isinstance(table, QtWidgets.QTableWidget):
                continue
            for row in range(table.rowCount()):
                table.setRowHidden(
                    row,
                    bool(needle) and not self._row_matches(table, row, needle))
            self._refresh_tab_label(index)

    @staticmethod
    def _row_matches(table: QtWidgets.QTableWidget, row: int,
                     needle: str) -> bool:
        """True when any cell of the row's visible columns matches."""
        for column in range(table.columnCount()):
            item = table.item(row, column)
            if item is not None and needle in item.text().lower():
                return True
        return False

    def _refresh_tab_label(self, index: int):
        """Set a tab's text to ``NAME (visible)`` plus ``*`` when filtered."""
        table = self.tabs.widget(index)
        if not isinstance(table, QtWidgets.QTableWidget):
            return
        name = str(table.property("tabName") or "")
        total = table.rowCount()
        visible = sum(1 for row in range(total)
                      if not table.isRowHidden(row))
        star = "*" if visible < total else ""
        self.tabs.setTabText(index, f"{name}{star} ({visible})")

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

    def _update_action_buttons(self):
        """Keep the Load / New Sequence buttons in sync with the selection."""
        self.load_button.setEnabled(self._current_record() is not None)
        table = self.tabs.currentWidget()
        has_selection = (
            isinstance(table, QtWidgets.QTableWidget)
            and table.selectionModel() is not None
            and table.selectionModel().hasSelection())
        self.sequence_button.setEnabled(has_selection)

    def _selected_records(self,
                          table: QtWidgets.QTableWidget
                          ) -> List[Dict[str, Any]]:
        """Records of the selected rows, in the table's visual order.

        The visual order follows the column the tab is currently sorted
        by, so a sequence built from a selection keeps that order.
        """
        records: List[Dict[str, Any]] = []
        selection_model = table.selectionModel()
        if selection_model is None:
            return records
        for index in selection_model.selectedRows(0):
            item = table.item(index.row(), 0)
            record = (item.data(QtCore.Qt.ItemDataRole.UserRole)
                      if item is not None else None)
            if isinstance(record, dict):
                records.append(record)
        return records

    def _show_table_context_menu(self, table: QtWidgets.QTableWidget, pos):
        records = self._selected_records(table)
        if not records:
            return
        menu = QtWidgets.QMenu(table)
        menu.addAction(
            f"New Sequence from Selection ({len(records)} file(s))")
        if menu.exec(table.viewport().mapToGlobal(pos)) is not None:
            self._new_sequence_from_selection()

    def _new_sequence_from_selection(self):
        """Create a session-only AMPA sequence from the selected rows."""
        table = self.tabs.currentWidget()
        if not isinstance(table, QtWidgets.QTableWidget):
            return
        records = self._selected_records(table)
        if not records:
            ui_api.show_info_dialog(
                "File Sorter",
                "Select one or more rows first "
                "(click, shift-click or ctrl-click).")
            return
        name = ui_api.prompt_user(
            "File Sorter",
            f"Name for the new sequence with {len(records)} file(s):")
        if name is None:
            return
        name = name.strip()
        if not name:
            ui_api.show_error_dialog(
                "File Sorter", "The sequence name is empty.")
            return
        if sequence_api.get_sequence_by_name(name) is not None:
            ui_api.show_error_dialog(
                "File Sorter",
                f"A sequence named '{name}' already exists.")
            return
        sequence = Sequence(
            name=name,
            frames=[record["path"] for record in records],
            temporary=True,
        )
        sequence_api.add_sequence(sequence)
        ui_api.write_to_statusbar(
            f"Created sequence '{name}' with {len(records)} file(s)",
            timeout=4000)
        log(__name__, __LOGMODULE__, "info",
            f"Created session sequence '{name}' with "
            f"{len(records)} file(s) from File Sorter")

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
