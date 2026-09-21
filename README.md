# not-tools-ampa-plugin

Remote plugin repository for [AMPA](https://github.com/samgrund/not-tools-ampa-plugin)
(Astro Multi-Purpose Analyzer). Each `*.py` file at the repo root is an
installable AMPA plugin; `plugins.json` is the registry manifest for
one-click installation via AMPA's Plugin Browser.

## Plugins

| File | Description |
|------|-------------|
| `file_sorter.py` | **File Sorter** — scan a folder recursively for FITS files and browse them grouped by `TCSTGT` (target) → `INSTRUME` (instrument). |

## Install

### Option A — Remote registry (one click)

In AMPA open **Plugins → Plugin Browser → Add Plugins → Remote Registry
URLs** and add:

```
https://github.com/samgrund/not-tools-ampa-plugin/raw/main/plugins.json
```

The **File Sorter** plugin appears in the browser list — install it with
the button, then restart AMPA when prompted.

### Option B — Clone as a local plugin directory

```console
$ git clone https://github.com/samgrund/not-tools-ampa-plugin.git
```

Then add the cloned folder under **Plugins → Plugin Browser → Add
Plugins → Local Plugin Directories** (e.g.
`/home/you/not-tools-ampa-plugin`) and restart AMPA. Update with
`git pull` + restart.

## File Sorter usage

1. **Plugins → NOT Toolkit → File Sorter** opens the plugin window.
2. Pick a folder with **Browse…** — every FITS file below it
   (`.fits`, `.fit`, `.fts`, plus `.gz` variants, case-insensitive) is
   scanned in the background with a progress bar and cancel button.
3. Files appear as a three-level tree:

   ```text
   TCSTGT (target)          INSTRUME (instrument)     files
   ├── M 31 (3)             ├── CCD1 (2)              ├── m31_a.fits
   │                        └── CCD2 (1)              └── ...
   ├── BIAS (1)             └── CCD1 (1)
   ├── (no TCSTGT)          └── (no INSTRUME)         (header missing)
   └── (unreadable)         └── (unreadable)          (header unreadable)
   ```

   Files without `TCSTGT` / `INSTRUME` headers are kept under
   placeholder groups, corrupt or truncated files under
   `(unreadable)` — nothing disappears from the listing.
4. Selecting a file shows its details (path, `TCSTGT`, `INSTRUME`, size
   and common observation headers such as `OBJECT`, `DATE-OBS`,
   `EXPTIME`).
5. **Load Selected** opens the selected file in the AMPA viewer.
   **Rescan** re-runs the scan (the folder is remembered between
   sessions).

The scan runs as a cancellable background task, so even huge data
folders never freeze the GUI. No extra dependencies — only `astropy`
and `PySide6`, which AMPA ships already.

## Development

Edit `file_sorter.py`, restart AMPA (local-directory plugins are
re-imported at every launch). When publishing:

1. Commit the changed `.py` and bump `version` in `plugins.json`.
2. Regenerate the `sha256:` checksum from the exact uploaded file
   (`sha256sum file_sorter.py`) and update `download_url`/`download_size`
   — a stale checksum breaks remote installation.
3. Push; users refresh (or restart AMPA) to pick up the update.
