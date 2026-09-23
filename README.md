# not-tools-ampa-plugin

This repository contains plugins for AMPA (Astro Multi-Purpose
Analyzer) specific to data from the Nordic Optical Telescope (NOT).

## Install

### Remote registry

In AMPA open **Plugins -> Plugin Browser -> Add Plugins -> Remote
Registry URLs** and add:

```
https://github.com/samgrund/not-tools-ampa-plugin/raw/main/plugins.json
```

Install the plugin from the browser list, then restart AMPA when
prompted.

### Local plugin directory

Clone the repository:

```console
$ git clone https://github.com/samgrund/not-tools-ampa-plugin.git
```

Add the cloned folder under **Plugins -> Plugin Browser -> Add
Plugins -> Local Plugin Directories**, and restart AMPA. Update with
`git pull` and a restart.

## Plugins

### File Sorter

Recursively scans a folder for FITS files and groups them in tabs by
instrument. Non-science files and files without a target go to their
own per-instrument CALIB tabs. Select a range of rows to create a
session-only AMPA sequence; double-click a file to load it into the
AMPA viewer.

### FIES Spectrum Viewer

Opens 1-D merged FIES echelle spectra (FIEStool products) and plots
flux against a wavelength axis calibrated from the FITS header
(`CRVAL1`/`CDELT1`/`CRPIX1`, air Angstroms). When AMPA already has a
FIES spectrum loaded, the viewer picks it up automatically; otherwise
use **Open...**.

The **Telluric lines** checkbox overlays atmospheric absorption lines
for comparison, colour-coded by absorber (O2 red, H2O blue). Hovering
the plot near a line identifies it: air and vacuum wavelength,
absorbing molecule and band (e.g. O2 A-band), and an indicative depth.
**Labels** annotates the strongest features directly on the plot. The
line table is derived from the DKIST Telluric Atlas (HITRAN2020 /
LBLRTM) and stored in vacuum wavelengths, which the plugin converts to
air with the Morton (1991) dispersion relation to match the FIES
wavelength scale. It can be regenerated with
`tools/generate_telluric_table.py`.
