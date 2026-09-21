"""File Sorter plugin for AMPA.

Recursively scans a user-selected folder for FITS files and presents them
in a tree grouped first by the ``TCSTGT`` header entry (telescope target)
and then by the ``INSTRUME`` key (instrument). Selecting a file shows its
details (grouping headers plus a few common observation keywords); the
**Load Selected** button opens it in the AMPA viewer.

Files missing ``TCSTGT`` / ``INSTRUME`` are kept under ``(no TCSTGT)`` /
``(no INSTRUME)`` placeholder groups; files whose header cannot be read
(corrupt, truncated, ...) land under ``(unreadable)``. The scan itself
runs on a background task with a progress bar and cancel support, so
huge data folders never freeze the GUI.

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

# Additional header keywords shown in the details pane (first hit wins;
# missing keys are simply skipped).
_DETAIL_KEYS = ("OBJECT", "DATE-OBS", "EXPTIME", "IMAGETYP", "FILTER",
                "TELESCOP", "OBSERVER", "NAXIS1", "NAXIS2")


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
            "records": [{"path", "tcstgt", "instrume", "header", "error"}],
        }

    Each record describes one file; ``header`` is the astropy header (or
    ``None`` for unreadable files, with ``error`` carrying the reason).
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
                "header": None,
                "error": str(exc),
            })
            continue
        records.append({
            "path": path,
            "tcstgt": tcstgt or _NO_TARGET,
            "instrume": instrume or _NO_INSTRUMENT,
            "header": header,
            "error": None,
        })
    return {"root": root, "records": records}


def build_tree(records: List[Dict[str, Any]]) -> Dict[str, Dict[str, List[Dict[str, Any]]]]:
    """Group scan records into ``{tcstgt: {instrume: [record, ...]}}``."""
    tree: Dict[str, Dict[str, List[Dict[str, Any]]]] = {}
    for record in records:
        tree.setdefault(record["tcstgt"], {}) \
            .setdefault(record["instrume"], []).append(record)
    return tree


def format_size(num_bytes: int) -> str:
    """Human-readable byte count."""
    size = float(num_bytes)
    for unit in ("B", "KiB", "MiB", "GiB", "TiB"):
        if size < 1024.0 or unit == "TiB":
            return f"{size:.1f} {unit}" if unit != "B" else f"{int(size)} B"
        size /= 1024.0
    return f"{size:.1f} TiB"


def format_details(record: Dict[str, Any]) -> str:
    """Format the details-pane text for one scan record."""
    path = record["path"]
    lines = [
        f"File:      {os.path.basename(path)}",
        f"Path:      {path}",
        f"TCSTGT:    {record['tcstgt']}",
        f"INSTRUME:  {record['instrume']}",
    ]
    try:
        lines.append(f"Size:      {format_size(os.path.getsize(path))}")
    except OSError:
        pass
    header = record.get("header")
    if header is None:
        lines.append(f"Error:     {record.get('error') or 'header unreadable'}")
    else:
        for key in _DETAIL_KEYS:
            if key in header:
                lines.append(f"{key:<10} {header[key]}")
    return "\n".join(lines)


# ======================================================================
# Plugin
# ======================================================================

class FileSorterPlugin(BaseModule):
    """NOT Toolkit ▸ File Sorter."""

    def __init__(self):
        super().__init__(
            title="File Sorter",
            category="NOT Toolkit",
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

        # Tree: TCSTGT -> INSTRUME -> files
        self.tree = QtWidgets.QTreeWidget()
        self.tree.setColumnCount(2)
        self.tree.setHeaderLabels(["Name", "Files"])
        self.tree.setRootIsDecorated(True)
        self.tree.setAlternatingRowColors(True)
        self.tree.header().setSectionResizeMode(
            0, QtWidgets.QHeaderView.ResizeMode.Stretch)
        self.tree.currentItemChanged.connect(self._on_tree_selection)
        layout.addWidget(self.tree, 1)

        # Details pane
        details_group = QtWidgets.QGroupBox("Details")
        details_layout = QtWidgets.QVBoxLayout(details_group)
        self.details_view = QtWidgets.QPlainTextEdit()
        self.details_view.setReadOnly(True)
        self.details_view.setMaximumHeight(180)
        details_layout.addWidget(self.details_view)
        layout.addWidget(details_group)

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

        window.resize(760, 600)
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
        self.tree.clear()
        self.details_view.clear()
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
        self._populate_tree(build_tree(self._records))
        self.rescan_button.setEnabled(True)
        unreadable = sum(1 for r in self._records if r["header"] is None)
        n_targets = len({r["tcstgt"] for r in self._records})
        n_instr = len({(r["tcstgt"], r["instrume"]) for r in self._records})
        msg = (f"{n_targets} target(s) · {n_instr} instrument group(s) · "
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

    # -- tree population ----------------------------------------------------

    def _populate_tree(self, tree):
        self.tree.clear()
        for tcstgt in sorted(tree, key=str.lower):
            instruments = tree[tcstgt]
            n_target = sum(len(recs) for recs in instruments.values())
            target_item = QtWidgets.QTreeWidgetItem(
                [str(tcstgt), str(n_target)])
            target_item.setData(
                0, QtCore.Qt.ItemDataRole.UserRole, "group")
            for instrume in sorted(instruments, key=str.lower):
                records = instruments[instrume]
                instr_item = QtWidgets.QTreeWidgetItem(
                    [str(instrume), str(len(records))])
                instr_item.setData(
                    0, QtCore.Qt.ItemDataRole.UserRole, "group")
                for record in sorted(records, key=lambda r: r["path"].lower()):
                    leaf = QtWidgets.QTreeWidgetItem(
                        [os.path.basename(record["path"]), ""])
                    leaf.setData(0, QtCore.Qt.ItemDataRole.UserRole, record)
                    leaf.setToolTip(0, record["path"])
                    instr_item.addChild(leaf)
                target_item.addChild(instr_item)
            self.tree.addTopLevelItem(target_item)
        self.tree.expandToDepth(0)

    # -- selection / details -------------------------------------------------

    def _on_tree_selection(self, current, _previous):
        if current is None:
            self.load_button.setEnabled(False)
            self.details_view.clear()
            return
        data = current.data(0, QtCore.Qt.ItemDataRole.UserRole)
        if data == "group":
            self.load_button.setEnabled(False)
            self.details_view.setPlainText(
                f"{current.text(0)} — {current.text(1)} file(s)")
        else:
            record = data
            self.load_button.setEnabled(True)
            self.details_view.setPlainText(format_details(record))

    def _load_selected(self):
        item = self.tree.currentItem()
        if item is None:
            return
        record = item.data(0, QtCore.Qt.ItemDataRole.UserRole)
        if not isinstance(record, dict):
            return
        path = record["path"]
        if record["header"] is None:
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
