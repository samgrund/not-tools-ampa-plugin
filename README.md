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
3. The files are shown in **one tab per instrument** — each tab a
   table with one row per file:

   ```text
   [ ALFOSC (23) ] [ CCD1 (7) ] [ (no INSTRUME) ] [ (unreadable) ]
    File     TCSTGT  GROUPID  BLOCKID  SEQID  OBJECT  IMAGETYP  FILTER  OBS_MODE  EXPTIME  DATE-OBS       FAFLTNM  FBFLTNM  ALGRNM
    a.fits   M 31    G-7      BLK-1    SEQ-9  M31     —         —       Imaging   300      2026-09-21 …   B_V      Empty    Grism#4
   ```

   - Common columns on every tab: `File`, `TCSTGT`, `GROUPID`,
     `BLOCKID`, `SEQID`, `OBJECT`, `IMAGETYP`, `FILTER`, `OBS_MODE`,
     `EXPTIME`, `DATE-OBS`; missing values show `—`.
   - Instruments with a dedicated view get extra columns — currently
     **ALFOSC** (`FAFLTNM`, `FBFLTNM`, `ALGRNM`); more (e.g. FIES)
     can be added to the registry in `file_sorter.py`.
   - Instrument names are normalised (raw `ALFOSC_FASU` groups under
     **ALFOSC**).
   - Files missing `INSTRUME` land on a `(no INSTRUME)` tab; corrupt or
     truncated files on an `(unreadable)` tab with an `Error` column —
     nothing disappears from the listing.
   - Click any column header to sort (`EXPTIME` sorts numerically);
     row tooltips show the absolute file path.
   - Rows are **shaded by `GROUPID`** — each group gets a soft pastel
     colour (cycled per tab so adjacent groups always differ; rows
     without a `GROUPID` stay unshaded). Grouped tabs drop the
     alternating-row stripes in favour of the shading.
4. **Double-click a row** (or select it and press **Load Selected**)
   to open the file in the AMPA viewer. Unreadable files ask for
   confirmation first.
5. **Rescan** re-runs the scan (the folder is remembered between
   sessions and scanned automatically when the plugin window opens).

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
