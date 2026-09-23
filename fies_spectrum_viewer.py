"""FIES Spectrum Viewer plugin for AMPA.

Open* 1-D merged FIES (Nordic Optical Telescope) echelle spectra as
produced by FIEStool and plots flux against wavelength.  The spectrum
curve line width is user-adjustable and persisted across sessions.
The first and last 2% of the wavelength range are trimmed off - the
edges of merged FIES spectra are dominated by noise blow-up.  The
wavelength
axis is calibrated from the FITS header linear solution (``CRVAL1`` /
``CDELT1`` / ``CRPIX1``, ``DC-FLAG`` = 0, ``WAT1_001`` "Wavelength
units=Angstroms"), i.e. air Angstroms in the observed frame - exactly
what the ThAr-based FIES wavelength calibration delivers.

**Telluric overlay.**  A checkbox draws atmospheric absorption lines for
comparison, colour-coded by absorber (O2 red, H2O blue), with a
user-adjustable line width (persisted across sessions).  Hovering the
plot near a line shows its wavelength (air and vacuum), the absorbing
molecule and band, and an indicative depth.  The embedded line table is
stored in *vacuum* wavelengths (the native convention of the HITRAN
database) and converted to air at runtime with the Edlen/Morton
dispersion relation so the lines match the spectrograph's air
wavelength scale.

Sources for the telluric table:

- DKIST Telluric Atlas (Schad), https://github.com/tschad/dkist_telluric_atlas
  (Zenodo DOI 10.5281/zenodo.13922154): synthetic LBLRTM/Py4CAtS
  transmission spectra built on the HITRAN2020 database
  (Gordon et al. 2022, JQSRT 277, 107949), US standard atmosphere above
  a 3 km base height, reference file PWV 3 mm / airmass 1.5.
- Air/vacuum conversion: Morton (1991, ApJS 77, 119):
  ``n - 1 = 2.7357e-4 + 131.4182 / wl^2 + 2.762e8 / wl^4`` (wl in A).

The table can be regenerated with ``tools/generate_telluric_table.py``
in the plugin repository, which re-derives line positions, depths and
species directly from the atlas files.

The plugin is distributed via the ``not-tools-ampa-plugin`` repository:
clone it and add the folder as a Local Plugin Directory, or install it
from the registry manifest (``plugins.json``) in the repo root.
"""

from __future__ import annotations

import os
from typing import List, Optional, Tuple

import numpy as np
from astropy.io import fits
from PySide6 import QtWidgets
import pyqtgraph as pg

from ampa.core.apis import settings_api, ui_api
from ampa.core.basemodule import BaseModule
from ampa.core.logging import log

__LOGMODULE__ = "FiesSpectrumViewer"

_SETTINGS_GROUP = "FIESSpectrumViewer"
_KEY_LAST_FOLDER = f"{_SETTINGS_GROUP}/last_folder"
_KEY_LINE_WIDTH = f"{_SETTINGS_GROUP}/line_width"
_KEY_SPECTRUM_WIDTH = f"{_SETTINGS_GROUP}/spectrum_width"

# Colours for the telluric overlay, per absorbing molecule.
_SPECIES_COLORS = {
    "O2": "#ff6b6b",
    "H2O": "#6bb8ff",
}

# How many strongest H2O lines receive text labels (band anchors for O2
# are always labelled); lines closer than this many Angstroms to an
# already-labelled line are skipped to avoid overlap.
_LABEL_MAX_H2O = 12
_LABEL_MIN_SEPARATION_A = 6.0

# Absorption features deeper than this (in the atlas' reference
# transmission) get text labels when labelling is enabled.
_LABEL_DEPTH_O2 = 0.4
_LABEL_DEPTH_H2O = 0.7

# Hover readout: telluric lines within this many plot pixels count as
# "near the cursor".
_HOVER_RADIUS_PX = 8.0

# Fraction of the spectrum trimmed from each end before display.  The
# edges of merged FIES spectra are dominated by noise blow-up, so the
# first and last 2% of the wavelength range are cut off.
_EDGE_TRIM_FRACTION = 0.02

# Telluric overlay line width (pixels).  User-adjustable; the spin box
# range below must stay in sync with these limits.
_DEFAULT_LINE_WIDTH = 0.75
_LINE_WIDTH_MIN = 0.25
_LINE_WIDTH_MAX = 5.0
_LINE_WIDTH_STEP = 0.25

# Spectrum curve line width (pixels), adjustable like the overlay
# width.  The default is the historical fixed width.
_SPECTRUM_COLOR = "#ffd54a"
_DEFAULT_SPECTRUM_WIDTH = 1.0
_SPECTRUM_WIDTH_MIN = 0.25
_SPECTRUM_WIDTH_MAX = 5.0
_SPECTRUM_WIDTH_STEP = 0.25


def _stored_float(key: str, default: float, lo: float, hi: float) -> float:
    """Persisted float setting, validated and clamped to ``[lo, hi]``.

    Falls back to *default* on a missing or invalid stored value.
    """
    try:
        value = float(settings_api.get_value(key, default))
    except (TypeError, ValueError):
        return default
    return min(max(value, lo), hi)


def stored_line_width() -> float:
    """Persisted telluric overlay line width (validated)."""
    return _stored_float(_KEY_LINE_WIDTH, _DEFAULT_LINE_WIDTH,
                         _LINE_WIDTH_MIN, _LINE_WIDTH_MAX)


def stored_spectrum_width() -> float:
    """Persisted spectrum curve line width (validated)."""
    return _stored_float(_KEY_SPECTRUM_WIDTH, _DEFAULT_SPECTRUM_WIDTH,
                         _SPECTRUM_WIDTH_MIN, _SPECTRUM_WIDTH_MAX)


# ======================================================================
# Pure helpers (module level so they can be tested without the GUI)
# ======================================================================

def vacuum_to_air(wl_vac):
    """Convert vacuum wavelengths (A) to air wavelengths (A).

    Uses the Morton (1991, ApJS 77, 119) dispersion relation, which is
    defined for the vacuum wavelength, so no iteration is needed.
    """
    wl_vac = np.asarray(wl_vac, dtype=float)
    n = 1.0 + 2.7357e-4 + 131.4182 / wl_vac ** 2 + 2.762e8 / wl_vac ** 4
    return wl_vac / n


def compute_wavelength_axis(header, n_pixels: int) -> np.ndarray:
    """Wavelength (air, A) of every pixel of a FIES merged spectrum.

    Implements the linear solution ``wl = CRVAL1 + (i + 1 - CRPIX1) *
    CDELT1`` with FITS 1-based pixel numbering, falling back to ``CD1_1``
    when ``CDELT1`` is absent.
    """
    crval = float(header["CRVAL1"])
    cdelt = header.get("CDELT1", header.get("CD1_1"))
    cdelt = float(cdelt)
    crpix = float(header.get("CRPIX1", 1.0))
    pixels = np.arange(1, n_pixels + 1, dtype=float)
    return crval + (pixels - crpix) * cdelt


def spectrum_error(header, data) -> Optional[str]:
    """Return a user-facing reason why *data* cannot be shown, or None.

    Checks the structure of a FIEStool merged 1-D product: 1-D flux
    array plus a linear wavelength solution.  Files from other
    instruments are accepted here (the viewer is wavelength-solution
    driven); :func:`is_fies_spectrum` gates the FIES-specific bits.
    """
    if data is None:
        return "the file contains no image data"
    if getattr(data, "ndim", 0) != 1:
        return (f"expected a 1-D merged spectrum, "
                f"got {data.ndim}-D data of shape {data.shape}")
    if "CRVAL1" not in header or ("CDELT1" not in header
                                   and "CD1_1" not in header):
        return ("no linear wavelength solution in the header "
                "(CRVAL1/CDELT1 missing)")
    dc_flag = header.get("DC-FLAG")
    if dc_flag not in (None, "", 0):
        return (f"non-linear wavelength solution (DC-FLAG = {dc_flag}) "
                "is not supported")
    return None


def is_fies_spectrum(header, data) -> bool:
    """True when *header*/*data* look like a FIES 1-D merged product."""
    instrume = str(header.get("INSTRUME", "") or "").strip().upper()
    return instrume == "FIES" and spectrum_error(header, data) is None


def read_spectrum(path: str):
    """Read a spectrum file; returns ``(data, wavelength, header)``.

    The first and last 2% of the wavelength range are trimmed (see
    ``_EDGE_TRIM_FRACTION``).  Raises ``ValueError`` with a user-facing
    message when the file cannot be displayed.
    """
    with fits.open(path, memmap=False) as hdul:
        data = hdul[0].data
        header = hdul[0].header
    reason = spectrum_error(header, data)
    if reason is None:
        reason = _sanity_check_axis(header, data)
    if reason is not None:
        raise ValueError(f"cannot display {os.path.basename(path)}: {reason}")
    wavelength = compute_wavelength_axis(header, int(data.shape[0]))
    flux = np.asarray(data, dtype=float)
    flux, wavelength = trim_edges(flux, wavelength)
    return flux, wavelength, header


def trim_edges(flux: np.ndarray, wavelength: np.ndarray):
    """Drop the first and last 2% of the pixels of a spectrum.

    On a merged FIES product the wavelength grid is linear, so this
    equals cutting 2% off each end of the wavelength range.
    """
    n = len(flux)
    cut = int(round(_EDGE_TRIM_FRACTION * n))
    if 2 * cut >= n:
        return flux, wavelength
    return flux[cut:n - cut], wavelength[cut:n - cut]


def _sanity_check_axis(header, data) -> Optional[str]:
    """Guard against a nonsense wavelength scale (all-NaN, zero span)."""
    wl = compute_wavelength_axis(header, int(data.shape[0]))
    if not np.all(np.isfinite(wl)):
        return "the wavelength solution contains non-finite values"
    if wl[-1] == wl[0]:
        return "the wavelength solution has zero dispersion"
    return None


def format_info_line(header) -> str:
    """One-line summary of the observation shown next to the file name."""
    parts = []
    object_ = str(header.get("OBJECT", "") or "").strip()
    if object_:
        parts.append(object_)
    fiber = str(header.get("FIFMSKNM", "") or "").strip()
    if fiber:
        parts.append(f"Fiber: {fiber}")
    exptime = header.get("EXPTIME")
    if exptime not in (None, ""):
        try:
            parts.append(f"Exp: {float(exptime):g} s")
        except (TypeError, ValueError):
            pass
    date_obs = str(header.get("DATE-OBS", "") or "").strip()
    if len(date_obs) >= 16:
        parts.append(date_obs[:10].replace("-", " ") + " "
                     + date_obs[11:16])
    vhelio = header.get("VHELIO")
    if vhelio not in (None, ""):
        try:
            parts.append(f"Vhelio: {float(vhelio):.2f} km/s")
        except (TypeError, ValueError):
            pass
    airmass = header.get("AIRMASS")
    if airmass not in (None, ""):
        try:
            parts.append(f"Airmass: {float(airmass):.2f}")
        except (TypeError, ValueError):
            pass
    return "   ".join(parts)


# ======================================================================
# Telluric line table (generated - see module docstring)
# ======================================================================

_TELLURIC_LINES_VACUUM = [
    (5899.793, "H2O", "H2O", 0.16),
    (5901.675, "H2O", "H2O", 0.15),
    (5903.103, "H2O", "H2O", 0.20),
    (5920.697, "H2O", "H2O", 0.17),
    (5921.284, "H2O", "H2O", 0.21),
    (5942.733, "H2O", "H2O", 0.17),
    (5944.213, "H2O", "H2O", 0.17),
    (5947.656, "H2O", "H2O", 0.15),
    (6278.334, "O2", "O2 gamma-band", 0.48),
    (6278.560, "O2", "O2 gamma-band", 0.61),
    (6279.056, "O2", "O2 gamma-band", 0.70),
    (6279.383, "O2", "O2 gamma-band", 0.56),
    (6279.816, "O2", "O2 gamma-band", 0.79),
    (6280.118, "O2", "O2 gamma-band", 0.17),
    (6280.620, "O2", "O2 gamma-band", 0.73),
    (6280.840, "O2", "O2 gamma-band", 0.80),
    (6281.637, "O2", "O2 gamma-band", 0.75),
    (6282.134, "O2", "O2 gamma-band", 0.78),
    (6282.919, "O2", "O2 gamma-band", 0.70),
    (6283.698, "O2", "O2 gamma-band", 0.70),
    (6284.465, "O2", "O2 gamma-band", 0.55),
    (6285.540, "O2", "O2 gamma-band", 0.49),
    (6286.281, "O2", "O2 gamma-band", 0.25),
    (6289.495, "O2", "O2 gamma-band", 0.43),
    (6291.142, "O2", "O2 gamma-band", 0.51),
    (6291.967, "O2", "O2 gamma-band", 0.66),
    (6293.905, "O2", "O2 gamma-band", 0.67),
    (6294.704, "O2", "O2 gamma-band", 0.76),
    (6296.927, "O2", "O2 gamma-band", 0.70),
    (6297.708, "O2", "O2 gamma-band", 0.78),
    (6300.202, "O2", "O2 gamma-band", 0.72),
    (6300.977, "O2", "O2 gamma-band", 0.76),
    (6303.750, "O2", "O2 gamma-band", 0.63),
    (6304.513, "O2", "O2 gamma-band", 0.68),
    (6307.559, "O2", "O2 gamma-band", 0.56),
    (6308.316, "O2", "O2 gamma-band", 0.59),
    (6311.634, "O2", "O2 gamma-band", 0.42),
    (6312.386, "O2", "O2 gamma-band", 0.46),
    (6315.985, "O2", "O2 gamma-band", 0.31),
    (6316.730, "O2", "O2 gamma-band", 0.32),
    (6320.604, "O2", "O2 gamma-band", 0.19),
    (6321.343, "O2", "O2 gamma-band", 0.20),
    (6477.609, "H2O", "H2O", 0.24),
    (6481.853, "H2O", "H2O", 0.17),
    (6485.037, "H2O", "H2O", 0.17),
    (6492.589, "H2O", "H2O", 0.19),
    (6497.669, "H2O", "H2O", 0.19),
    (6516.533, "H2O", "H2O", 0.21),
    (6518.351, "H2O", "H2O", 0.17),
    (6545.714, "H2O", "H2O", 0.20),
    (6554.438, "H2O", "H2O", 0.17),
    (6869.141, "O2", "O2 B-band", 1.00),
    (6869.450, "O2", "O2 B-band", 1.00),
    (6869.759, "O2", "O2 B-band", 0.64),
    (6870.007, "O2", "O2 B-band", 1.00),
    (6870.446, "O2", "O2 B-band", 1.00),
    (6870.818, "O2", "O2 B-band", 1.00),
    (6871.525, "O2", "O2 B-band", 0.38),
    (6871.875, "O2", "O2 B-band", 1.00),
    (6872.851, "O2", "O2 B-band", 1.00),
    (6873.188, "O2", "O2 B-band", 1.00),
    (6874.150, "O2", "O2 B-band", 1.00),
    (6874.748, "O2", "O2 B-band", 1.00),
    (6875.697, "O2", "O2 B-band", 1.00),
    (6876.557, "O2", "O2 B-band", 1.00),
    (6877.499, "O2", "O2 B-band", 1.00),
    (6878.620, "O2", "O2 B-band", 1.00),
    (6879.542, "O2", "O2 B-band", 1.00),
    (6880.945, "O2", "O2 B-band", 1.00),
    (6881.833, "O2", "O2 B-band", 1.00),
    (6885.736, "O2", "O2 B-band", 1.00),
    (6887.658, "O2", "O2 B-band", 1.00),
    (6888.650, "O2", "O2 B-band", 1.00),
    (6890.854, "O2", "O2 B-band", 1.00),
    (6891.812, "O2", "O2 B-band", 1.00),
    (6894.280, "O2", "O2 B-band", 1.00),
    (6895.217, "O2", "O2 B-band", 1.00),
    (6897.948, "O2", "O2 B-band", 1.00),
    (6898.873, "O2", "O2 B-band", 1.00),
    (6901.861, "O2", "O2 B-band", 1.00),
    (6902.779, "O2", "O2 B-band", 1.00),
    (6906.031, "O2", "O2 B-band", 1.00),
    (6906.935, "O2", "O2 B-band", 1.00),
    (6910.445, "O2", "O2 B-band", 1.00),
    (6911.343, "O2", "O2 B-band", 1.00),
    (6915.111, "O2", "O2 B-band", 1.00),
    (6916.003, "O2", "O2 B-band", 1.00),
    (6920.036, "O2", "O2 B-band", 1.00),
    (6920.915, "O2", "O2 B-band", 1.00),
    (6925.215, "O2", "O2 B-band", 0.97),
    (6926.087, "O2", "O2 B-band", 0.98),
    (6930.646, "O2", "O2 B-band", 0.82),
    (6931.228, "H2O", "H2O", 0.19),
    (6931.512, "O2", "O2 B-band", 0.83),
    (6935.735, "H2O", "H2O", 0.22),
    (6936.345, "O2", "O2 B-band", 0.55),
    (6937.205, "O2", "O2 B-band", 0.56),
    (6939.620, "H2O", "H2O", 0.31),
    (6941.535, "H2O", "H2O", 0.25),
    (6942.105, "H2O", "H2O", 0.35),
    (6943.146, "H2O", "H2O", 0.37),
    (6944.077, "H2O", "H2O", 0.29),
    (6944.292, "H2O", "H2O", 0.16),
    (6945.723, "H2O", "H2O", 0.37),
    (6949.460, "H2O", "H2O", 0.48),
    (6951.010, "H2O", "H2O", 0.19),
    (6952.685, "H2O", "H2O", 0.15),
    (6955.502, "H2O", "H2O", 0.20),
    (6958.333, "H2O", "H2O", 0.48),
    (6961.382, "H2O", "H2O", 0.30),
    (6963.185, "H2O", "H2O", 0.42),
    (6988.514, "H2O", "H2O", 0.33),
    (6990.918, "H2O", "H2O", 0.41),
    (6992.310, "H2O", "H2O", 0.15),
    (6995.450, "H2O", "H2O", 0.25),
    (6996.045, "H2O", "H2O", 0.18),
    (7000.902, "H2O", "H2O", 0.32),
    (7006.694, "H2O", "H2O", 0.26),
    (7018.398, "H2O", "H2O", 0.36),
    (7025.448, "H2O", "H2O", 0.29),
    (7028.884, "H2O", "H2O", 0.16),
    (7029.425, "H2O", "H2O", 0.29),
    (7039.484, "H2O", "H2O", 0.18),
    (7041.744, "H2O", "H2O", 0.24),
    (7052.809, "H2O", "H2O", 0.16),
    (7169.334, "H2O", "H2O", 0.28),
    (7169.886, "H2O", "H2O", 0.58),
    (7172.067, "H2O", "H2O", 0.29),
    (7172.303, "H2O", "H2O", 0.15),
    (7172.554, "H2O", "H2O", 0.38),
    (7174.706, "H2O", "H2O", 0.47),
    (7175.374, "H2O", "H2O", 0.60),
    (7175.754, "H2O", "H2O", 0.35),
    (7176.148, "H2O", "H2O", 0.34),
    (7178.137, "H2O", "H2O", 0.33),
    (7179.098, "H2O", "H2O", 0.41),
    (7179.350, "H2O", "H2O", 0.77),
    (7179.615, "H2O", "H2O", 0.39),
    (7180.405, "H2O", "H2O", 0.25),
    (7183.522, "H2O", "H2O", 0.88),
    (7183.744, "H2O", "H2O", 0.49),
    (7186.317, "H2O", "H2O", 0.25),
    (7186.518, "H2O", "H2O", 0.94),
    (7188.128, "H2O", "H2O", 0.64),
    (7188.365, "H2O", "H2O", 0.82),
    (7188.991, "H2O", "H2O", 0.52),
    (7189.386, "H2O", "H2O", 0.95),
    (7193.471, "H2O", "H2O", 0.96),
    (7193.852, "H2O", "H2O", 0.45),
    (7194.456, "H2O", "H2O", 0.30),
    (7195.550, "H2O", "H2O", 0.68),
    (7195.752, "H2O", "H2O", 0.69),
    (7197.032, "H2O", "H2O", 0.48),
    (7197.781, "H2O", "H2O", 0.22),
    (7199.221, "H2O", "H2O", 0.36),
    (7199.854, "H2O", "H2O", 0.41),
    (7200.423, "H2O", "H2O", 0.57),
    (7202.526, "H2O", "H2O", 0.78),
    (7203.189, "H2O", "H2O", 0.94),
    (7205.840, "H2O", "H2O", 0.25),
    (7206.294, "H2O", "H2O", 0.89),
    (7208.420, "H2O", "H2O", 0.96),
    (7211.513, "H2O", "H2O", 0.54),
    (7213.201, "H2O", "H2O", 0.21),
    (7225.625, "H2O", "H2O", 0.63),
    (7229.492, "H2O", "H2O", 0.54),
    (7234.236, "H2O", "H2O", 0.40),
    (7234.901, "H2O", "H2O", 0.91),
    (7236.399, "H2O", "H2O", 0.38),
    (7236.732, "H2O", "H2O", 0.96),
    (7238.129, "H2O", "H2O", 0.66),
    (7241.821, "H2O", "H2O", 0.28),
    (7242.632, "H2O", "H2O", 0.83),
    (7245.479, "H2O", "H2O", 0.57),
    (7245.711, "H2O", "H2O", 0.90),
    (7247.682, "H2O", "H2O", 0.61),
    (7249.241, "H2O", "H2O", 0.27),
    (7250.952, "H2O", "H2O", 0.19),
    (7252.228, "H2O", "H2O", 0.35),
    (7254.382, "H2O", "H2O", 0.87),
    (7254.861, "H2O", "H2O", 0.30),
    (7255.231, "H2O", "H2O", 0.57),
    (7255.724, "H2O", "H2O", 0.70),
    (7259.382, "H2O", "H2O", 0.34),
    (7259.948, "H2O", "H2O", 0.45),
    (7262.722, "H2O", "H2O", 0.21),
    (7263.434, "H2O", "H2O", 0.15),
    (7264.037, "H2O", "H2O", 0.27),
    (7264.981, "H2O", "H2O", 0.27),
    (7266.391, "H2O", "H2O", 0.22),
    (7266.609, "H2O", "H2O", 0.58),
    (7267.605, "H2O", "H2O", 0.91),
    (7271.756, "H2O", "H2O", 0.27),
    (7272.134, "H2O", "H2O", 0.22),
    (7274.963, "H2O", "H2O", 0.89),
    (7277.422, "H2O", "H2O", 0.67),
    (7278.325, "H2O", "H2O", 0.27),
    (7278.558, "H2O", "H2O", 0.25),
    (7278.856, "H2O", "H2O", 0.29),
    (7279.140, "H2O", "H2O", 0.17),
    (7279.402, "H2O", "H2O", 0.85),
    (7280.086, "H2O", "H2O", 0.36),
    (7281.703, "H2O", "H2O", 0.17),
    (7284.295, "H2O", "H2O", 0.33),
    (7289.396, "H2O", "H2O", 0.69),
    (7290.140, "H2O", "H2O", 0.43),
    (7292.407, "H2O", "H2O", 0.81),
    (7293.107, "H2O", "H2O", 0.30),
    (7294.187, "H2O", "H2O", 0.32),
    (7294.705, "H2O", "H2O", 0.37),
    (7295.390, "H2O", "H2O", 0.21),
    (7297.054, "H2O", "H2O", 0.37),
    (7301.652, "H2O", "H2O", 0.21),
    (7301.930, "H2O", "H2O", 0.33),
    (7304.143, "H2O", "H2O", 0.25),
    (7305.216, "H2O", "H2O", 0.62),
    (7306.218, "H2O", "H2O", 0.66),
    (7310.800, "H2O", "H2O", 0.28),
    (7311.538, "H2O", "H2O", 0.54),
    (7312.211, "H2O", "H2O", 0.18),
    (7314.631, "H2O", "H2O", 0.34),
    (7317.536, "H2O", "H2O", 0.35),
    (7319.314, "H2O", "H2O", 0.38),
    (7320.156, "H2O", "H2O", 0.15),
    (7320.398, "H2O", "H2O", 0.19),
    (7320.727, "H2O", "H2O", 0.48),
    (7322.894, "H2O", "H2O", 0.16),
    (7326.000, "H2O", "H2O", 0.22),
    (7329.393, "H2O", "H2O", 0.22),
    (7332.846, "H2O", "H2O", 0.17),
    (7335.721, "H2O", "H2O", 0.35),
    (7337.364, "H2O", "H2O", 0.22),
    (7345.983, "H2O", "H2O", 0.19),
    (7351.524, "H2O", "H2O", 0.21),
    (7362.383, "H2O", "H2O", 0.19),
    (7365.778, "H2O", "H2O", 0.22),
    (7370.501, "H2O", "H2O", 0.15),
    (7371.245, "H2O", "H2O", 0.22),
    (7385.759, "H2O", "H2O", 0.17),
    (7595.808, "O2", "O2 A-band", 1.00),
    (7596.089, "O2", "O2 A-band", 1.00),
    (7596.377, "O2", "O2 A-band", 0.39),
    (7596.605, "O2", "O2 A-band", 1.00),
    (7597.084, "O2", "O2 A-band", 1.00),
    (7597.350, "O2", "O2 A-band", 1.00),
    (7597.859, "O2", "O2 A-band", 1.00),
    (7598.315, "O2", "O2 A-band", 1.00),
    (7598.588, "O2", "O2 A-band", 1.00),
    (7598.862, "O2", "O2 A-band", 0.44),
    (7599.067, "O2", "O2 A-band", 0.36),
    (7599.553, "O2", "O2 A-band", 1.00),
    (7600.100, "O2", "O2 A-band", 0.59),
    (7601.324, "O2", "O2 A-band", 0.58),
    (7601.560, "O2", "O2 A-band", 0.65),
    (7602.366, "O2", "O2 A-band", 0.94),
    (7602.761, "O2", "O2 A-band", 0.92),
    (7603.217, "O2", "O2 A-band", 0.73),
    (7603.567, "O2", "O2 A-band", 0.86),
    (7604.130, "O2", "O2 A-band", 0.90),
    (7605.095, "O2", "O2 A-band", 0.82),
    (7605.324, "O2", "O2 A-band", 0.84),
    (7606.282, "O2", "O2 A-band", 0.92),
    (7607.172, "O2", "O2 A-band", 0.85),
    (7607.415, "O2", "O2 A-band", 0.73),
    (7608.290, "O2", "O2 A-band", 0.86),
    (7608.587, "O2", "O2 A-band", 0.84),
    (7609.462, "O2", "O2 A-band", 0.95),
    (7609.766, "O2", "O2 A-band", 0.83),
    (7610.687, "O2", "O2 A-band", 0.83),
    (7611.022, "O2", "O2 A-band", 0.68),
    (7611.966, "O2", "O2 A-band", 0.83),
    (7612.179, "O2", "O2 A-band", 0.63),
    (7613.108, "O2", "O2 A-band", 0.72),
    (7613.481, "O2", "O2 A-band", 0.41),
    (7613.702, "O2", "O2 A-band", 0.52),
    (7614.425, "O2", "O2 A-band", 0.85),
    (7614.684, "O2", "O2 A-band", 0.72),
    (7615.803, "O2", "O2 A-band", 0.54),
    (7616.123, "O2", "O2 A-band", 0.52),
    (7616.611, "O2", "O2 A-band", 0.30),
    (7617.707, "O2", "O2 A-band", 0.28),
    (7621.799, "O2", "O2 A-band", 0.41),
    (7622.180, "O2", "O2 A-band", 0.39),
    (7622.409, "O2", "O2 A-band", 0.20),
    (7623.430, "O2", "O2 A-band", 0.68),
    (7623.903, "O2", "O2 A-band", 0.51),
    (7624.605, "O2", "O2 A-band", 0.24),
    (7625.108, "O2", "O2 A-band", 0.81),
    (7625.657, "O2", "O2 A-band", 0.80),
    (7626.839, "O2", "O2 A-band", 0.91),
    (7627.456, "O2", "O2 A-band", 0.67),
    (7628.257, "O2", "O2 A-band", 0.31),
    (7628.624, "O2", "O2 A-band", 0.79),
    (7629.425, "O2", "O2 A-band", 0.76),
    (7631.203, "O2", "O2 A-band", 0.71),
    (7632.103, "O2", "O2 A-band", 0.32),
    (7632.347, "O2", "O2 A-band", 0.78),
    (7635.141, "O2", "O2 A-band", 0.73),
    (7636.279, "O2", "O2 A-band", 0.73),
    (7639.288, "O2", "O2 A-band", 0.68),
    (7640.411, "O2", "O2 A-band", 0.73),
    (7641.443, "O2", "O2 A-band", 0.87),
    (7642.559, "O2", "O2 A-band", 0.89),
    (7643.644, "O2", "O2 A-band", 0.60),
    (7644.752, "O2", "O2 A-band", 0.57),
    (7645.899, "O2", "O2 A-band", 0.64),
    (7647.008, "O2", "O2 A-band", 0.69),
    (7648.209, "O2", "O2 A-band", 0.48),
    (7649.310, "O2", "O2 A-band", 0.47),
    (7650.564, "O2", "O2 A-band", 0.44),
    (7651.666, "O2", "O2 A-band", 0.52),
    (7652.975, "O2", "O2 A-band", 0.35),
    (7654.069, "O2", "O2 A-band", 0.33),
    (7655.447, "O2", "O2 A-band", 0.28),
    (7656.534, "O2", "O2 A-band", 0.42),
    (7657.958, "O2", "O2 A-band", 0.24),
    (7659.046, "O2", "O2 A-band", 0.23),
    (7660.532, "O2", "O2 A-band", 0.18),
    (7661.474, "O2", "O2 A-band", 1.00),
    (7662.263, "O2", "O2 A-band", 0.19),
    (7662.562, "O2", "O2 A-band", 1.00),
    (7663.160, "O2", "O2 A-band", 0.16),
    (7666.985, "O2", "O2 A-band", 1.00),
    (7668.058, "O2", "O2 A-band", 1.00),
    (7672.722, "O2", "O2 A-band", 1.00),
    (7673.781, "O2", "O2 A-band", 1.00),
    (7678.678, "O2", "O2 A-band", 1.00),
    (7679.738, "O2", "O2 A-band", 1.00),
    (7684.877, "O2", "O2 A-band", 1.00),
    (7685.922, "O2", "O2 A-band", 1.00),
    (7691.304, "O2", "O2 A-band", 0.91),
    (7692.343, "O2", "O2 A-band", 0.91),
    (7697.968, "O2", "O2 A-band", 0.59),
    (7698.992, "O2", "O2 A-band", 0.62),
    (7704.861, "O2", "O2 A-band", 0.28),
    (7705.885, "O2", "O2 A-band", 0.30),
    (7895.684, "H2O", "H2O", 0.17),
    (7898.211, "H2O", "H2O", 0.20),
    (7903.947, "H2O", "H2O", 0.21),
    (7910.937, "H2O", "H2O", 0.25),
    (7922.852, "H2O", "H2O", 0.25),
    (7926.537, "H2O", "H2O", 0.17),
    (7930.803, "H2O", "H2O", 0.24),
    (7960.686, "H2O", "H2O", 0.22),
    (7962.932, "H2O", "H2O", 0.18),
    (7965.329, "H2O", "H2O", 0.18),
    (8002.501, "H2O", "H2O", 0.21),
    (8009.699, "H2O", "H2O", 0.19),
    (8015.148, "H2O", "H2O", 0.18),
    (8116.184, "H2O", "H2O", 0.17),
    (8127.677, "H2O", "H2O", 0.16),
    (8129.091, "H2O", "H2O", 0.20),
    (8132.254, "H2O", "H2O", 0.19),
    (8132.701, "H2O", "H2O", 0.33),
    (8133.457, "H2O", "H2O", 0.19),
    (8135.792, "H2O", "H2O", 0.18),
    (8136.020, "H2O", "H2O", 0.64),
    (8137.290, "H2O", "H2O", 0.50),
    (8138.763, "H2O", "H2O", 0.25),
    (8141.945, "H2O", "H2O", 0.26),
    (8142.898, "H2O", "H2O", 0.57),
    (8144.184, "H2O", "H2O", 0.65),
    (8146.034, "H2O", "H2O", 0.43),
    (8146.433, "H2O", "H2O", 0.18),
    (8148.453, "H2O", "H2O", 0.37),
    (8149.431, "H2O", "H2O", 0.43),
    (8150.311, "H2O", "H2O", 0.21),
    (8150.637, "H2O", "H2O", 0.65),
    (8151.518, "H2O", "H2O", 0.30),
    (8151.941, "H2O", "H2O", 0.50),
    (8154.746, "H2O", "H2O", 0.63),
    (8155.953, "H2O", "H2O", 0.46),
    (8156.655, "H2O", "H2O", 0.29),
    (8156.883, "H2O", "H2O", 0.86),
    (8160.269, "H2O", "H2O", 0.85),
    (8163.680, "H2O", "H2O", 0.94),
    (8164.227, "H2O", "H2O", 0.63),
    (8164.595, "H2O", "H2O", 0.90),
    (8166.416, "H2O", "H2O", 0.16),
    (8166.792, "H2O", "H2O", 0.97),
    (8167.584, "H2O", "H2O", 0.27),
    (8171.072, "H2O", "H2O", 0.71),
    (8171.636, "H2O", "H2O", 0.56),
    (8172.240, "H2O", "H2O", 0.98),
    (8179.231, "H2O", "H2O", 0.99),
    (8180.180, "H2O", "H2O", 0.72),
    (8180.744, "H2O", "H2O", 0.47),
    (8181.308, "H2O", "H2O", 0.77),
    (8184.099, "H2O", "H2O", 0.86),
    (8185.343, "H2O", "H2O", 0.19),
    (8188.626, "H2O", "H2O", 0.80),
    (8191.525, "H2O", "H2O", 0.99),
    (8195.368, "H2O", "H2O", 0.94),
    (8199.958, "H2O", "H2O", 0.97),
    (8201.270, "H2O", "H2O", 0.30),
    (8202.255, "H2O", "H2O", 0.39),
    (8202.952, "H2O", "H2O", 0.63),
    (8211.824, "H2O", "H2O", 0.20),
    (8220.377, "H2O", "H2O", 0.74),
    (8223.822, "H2O", "H2O", 0.52),
    (8226.257, "H2O", "H2O", 0.62),
    (8226.725, "H2O", "H2O", 0.17),
    (8227.951, "H2O", "H2O", 0.45),
    (8229.227, "H2O", "H2O", 0.98),
    (8230.239, "H2O", "H2O", 0.69),
    (8230.519, "H2O", "H2O", 0.97),
    (8231.005, "H2O", "H2O", 0.73),
    (8232.026, "H2O", "H2O", 0.72),
    (8232.750, "H2O", "H2O", 0.16),
    (8233.549, "H2O", "H2O", 0.94),
    (8233.993, "H2O", "H2O", 0.57),
    (8236.175, "H2O", "H2O", 0.80),
    (8236.900, "H2O", "H2O", 0.19),
    (8239.611, "H2O", "H2O", 0.41),
    (8242.198, "H2O", "H2O", 0.31),
    (8245.405, "H2O", "H2O", 0.21),
    (8245.760, "H2O", "H2O", 0.88),
    (8255.000, "H2O", "H2O", 0.42),
    (8258.790, "H2O", "H2O", 0.91),
    (8261.970, "H2O", "H2O", 0.62),
    (8265.721, "H2O", "H2O", 0.53),
    (8274.322, "H2O", "H2O", 0.64),
    (8276.632, "H2O", "H2O", 0.93),
    (8278.966, "H2O", "H2O", 0.57),
    (8281.880, "H2O", "H2O", 0.74),
    (8284.307, "H2O", "H2O", 0.92),
    (8290.233, "H2O", "H2O", 0.94),
    (8290.547, "H2O", "H2O", 0.17),
    (8291.808, "H2O", "H2O", 0.58),
    (8296.445, "H2O", "H2O", 0.79),
    (8296.826, "H2O", "H2O", 0.47),
    (8297.589, "H2O", "H2O", 0.18),
    (8302.652, "H2O", "H2O", 0.55),
    (8306.589, "H2O", "H2O", 0.52),
    (8307.378, "H2O", "H2O", 0.88),
    (8314.251, "H2O", "H2O", 0.43),
    (8316.164, "H2O", "H2O", 0.17),
    (8318.517, "H2O", "H2O", 0.40),
    (8320.431, "H2O", "H2O", 0.73),
    (8323.535, "H2O", "H2O", 0.78),
    (8323.868, "H2O", "H2O", 0.77),
    (8331.980, "H2O", "H2O", 0.51),
    (8335.879, "H2O", "H2O", 0.31),
    (8337.805, "H2O", "H2O", 0.19),
    (8341.349, "H2O", "H2O", 0.60),
    (8344.595, "H2O", "H2O", 0.16),
    (8351.465, "H2O", "H2O", 0.20),
    (8355.960, "H2O", "H2O", 0.18),
    (8359.336, "H2O", "H2O", 0.36),
    (8359.746, "H2O", "H2O", 0.16),
    (8364.604, "H2O", "H2O", 0.20),
    (8369.632, "H2O", "H2O", 0.22),
    (8378.685, "H2O", "H2O", 0.20),
    (8931.521, "H2O", "H2O", 0.21),
    (8932.727, "H2O", "H2O", 0.26),
    (8936.525, "H2O", "H2O", 0.28),
    (8942.657, "H2O", "H2O", 0.29),
    (8944.803, "H2O", "H2O", 0.23),
    (8948.811, "H2O", "H2O", 0.18),
    (8949.331, "H2O", "H2O", 0.41),
    (8951.067, "H2O", "H2O", 0.28),
    (8954.209, "H2O", "H2O", 0.25),
    (8954.639, "H2O", "H2O", 0.49),
    (8956.788, "H2O", "H2O", 0.62),
    (8957.424, "H2O", "H2O", 0.79),
    (8960.864, "H2O", "H2O", 0.39),
    (8964.799, "H2O", "H2O", 0.83),
    (8965.095, "H2O", "H2O", 0.52),
    (8965.955, "H2O", "H2O", 0.77),
    (8967.937, "H2O", "H2O", 0.88),
    (8968.879, "H2O", "H2O", 0.64),
    (8973.535, "H2O", "H2O", 0.62),
    (8974.091, "H2O", "H2O", 0.95),
    (8975.366, "H2O", "H2O", 0.29),
    (8976.730, "H2O", "H2O", 0.53),
    (8977.242, "H2O", "H2O", 0.69),
    (8978.894, "H2O", "H2O", 0.33),
    (8982.953, "H2O", "H2O", 0.98),
    (8983.870, "H2O", "H2O", 0.16),
    (8989.073, "H2O", "H2O", 0.96),
    (8989.801, "H2O", "H2O", 0.48),
    (8990.089, "H2O", "H2O", 0.99),
    (8990.619, "H2O", "H2O", 0.46),
    (8991.545, "H2O", "H2O", 0.80),
    (8992.012, "H2O", "H2O", 0.26),
    (8993.299, "H2O", "H2O", 0.97),
    (8994.171, "H2O", "H2O", 1.00),
    (8994.494, "H2O", "H2O", 0.83),
    (8995.520, "H2O", "H2O", 0.24),
]

# Precomputed (air wavelength, vacuum wavelength, species, band label,
# indicative depth) tuples sorted by air wavelength.  Built once at
# import time from the vacuum table above.
TELLURIC_LINES: List[Tuple[float, float, str, str, float]] = sorted(
    (float(vacuum_to_air(vac)), float(vac), species, label, float(depth))
    for vac, species, label, depth in _TELLURIC_LINES_VACUUM
)

# Air wavelengths as a flat array for nearest-line lookups on hover.
_TELLURIC_AIR = np.array([entry[0] for entry in TELLURIC_LINES])


# ======================================================================
# Plugin
# ======================================================================

class FiesSpectrumViewerPlugin(BaseModule):
    """Plugins ▸ NOT Toolkit ▸ FIES Spectrum Viewer."""

    def __init__(self):
        super().__init__(
            title="FIES Spectrum Viewer",
            category="Plugins",
            section="NOT Toolkit",
        )
        settings_api.define_setting(
            group=_SETTINGS_GROUP,
            key="last_folder",
            default_value="",
            value_type=str,
            description="Last folder a FIES spectrum was opened from.",
        )
        settings_api.define_setting(
            group=_SETTINGS_GROUP,
            key="line_width",
            default_value=_DEFAULT_LINE_WIDTH,
            value_type=float,
            description="Width of the telluric overlay lines (pixels).",
        )
        settings_api.define_setting(
            group=_SETTINGS_GROUP,
            key="spectrum_width",
            default_value=_DEFAULT_SPECTRUM_WIDTH,
            value_type=float,
            description="Width of the spectrum curve line (pixels).",
        )
        self._wavelength: Optional[np.ndarray] = None
        self._flux: Optional[np.ndarray] = None
        self._header = None
        self._path: Optional[str] = None
        self._telluric_items: List[pg.GraphicsObject] = []
        self._labelled_air: List[float] = []
        self._line_width: float = stored_line_width()
        self._spectrum_width: float = stored_spectrum_width()

    # -- lifecycle --------------------------------------------------------

    def show_gui(self):
        """Create the window lazily, then show/raise/focus it."""
        if self.gui_widget is None:
            self._build_window()
            self._try_autoload_current_frame()
        self.gui_widget.show()
        self.gui_widget.raise_()
        self.gui_widget.activateWindow()

    # -- GUI creation -----------------------------------------------------

    def _build_window(self):
        window = QtWidgets.QWidget()
        window.setWindowTitle("FIES Spectrum Viewer")
        window.resize(1000, 550)
        layout = QtWidgets.QVBoxLayout(window)

        # ---- File row ----
        file_row = QtWidgets.QHBoxLayout()
        self.open_button = QtWidgets.QPushButton("Open...")
        self.open_button.clicked.connect(self._choose_file)
        file_row.addWidget(self.open_button)
        self.file_label = QtWidgets.QLabel("No spectrum loaded")
        self.file_label.setStyleSheet("font-weight: bold;")
        file_row.addWidget(self.file_label, 1)
        layout.addLayout(file_row)

        # ---- Observation info row ----
        self.info_label = QtWidgets.QLabel("")
        self.info_label.setStyleSheet("color: gray;")
        layout.addWidget(self.info_label)

        # ---- Line width / telluric overlay controls ----
        overlay_row = QtWidgets.QHBoxLayout()
        spectrum_width_label = QtWidgets.QLabel("Spectrum width")
        spectrum_width_label.setStyleSheet("color: gray;")
        overlay_row.addWidget(spectrum_width_label)
        self.spectrum_width_spin = QtWidgets.QDoubleSpinBox()
        self.spectrum_width_spin.setRange(_SPECTRUM_WIDTH_MIN,
                                          _SPECTRUM_WIDTH_MAX)
        self.spectrum_width_spin.setSingleStep(_SPECTRUM_WIDTH_STEP)
        self.spectrum_width_spin.setDecimals(2)
        self.spectrum_width_spin.setValue(stored_spectrum_width())
        self.spectrum_width_spin.setToolTip(
            "Width of the spectrum curve line (pixels)")
        self.spectrum_width_spin.valueChanged.connect(
            self._on_spectrum_width_changed)
        overlay_row.addWidget(self.spectrum_width_spin)

        self.telluric_check = QtWidgets.QCheckBox("Telluric lines")
        self.telluric_check.setToolTip(
            "Overlay atmospheric absorption lines "
            "(vacuum wavelengths converted to air)")
        self.telluric_check.toggled.connect(self._on_overlay_toggled)
        overlay_row.addWidget(self.telluric_check)
        self.labels_check = QtWidgets.QCheckBox("Labels")
        self.labels_check.setToolTip(
            "Text labels on the strongest telluric features "
            "(hover any line for details)")
        self.labels_check.setEnabled(False)
        self.labels_check.toggled.connect(self._on_labels_toggled)
        overlay_row.addWidget(self.labels_check)

        width_label = QtWidgets.QLabel("Line width")
        width_label.setStyleSheet("color: gray;")
        overlay_row.addWidget(width_label)
        self.width_spin = QtWidgets.QDoubleSpinBox()
        self.width_spin.setRange(_LINE_WIDTH_MIN, _LINE_WIDTH_MAX)
        self.width_spin.setSingleStep(_LINE_WIDTH_STEP)
        self.width_spin.setDecimals(2)
        self.width_spin.setValue(stored_line_width())
        self.width_spin.setToolTip(
            "Width of the telluric overlay lines (pixels)")
        self.width_spin.setEnabled(False)
        self.width_spin.valueChanged.connect(self._on_line_width_changed)
        overlay_row.addWidget(self.width_spin)
        overlay_row.addStretch()
        self.legend_o2 = QtWidgets.QLabel("O2")
        self.legend_o2.setStyleSheet(f"color: {_SPECIES_COLORS['O2']};"
                                    " font-weight: bold;")
        self.legend_h2o = QtWidgets.QLabel("H2O")
        self.legend_h2o.setStyleSheet(f"color: {_SPECIES_COLORS['H2O']};"
                                       " font-weight: bold;")
        self.legend_o2.setVisible(False)
        self.legend_h2o.setVisible(False)
        overlay_row.addWidget(self.legend_o2)
        overlay_row.addWidget(self.legend_h2o)
        layout.addLayout(overlay_row)

        # ---- Plot ----
        self.plot_widget = pg.PlotWidget()
        self._plot_item = self.plot_widget.getPlotItem()
        self._plot_item.setLabel("bottom", "Wavelength (air)",
                                 units="Angstrom")
        self._plot_item.showGrid(x=True, y=True, alpha=0.25)
        self._curve = self._plot_item.plot(
            pen=pg.mkPen(_SPECTRUM_COLOR, width=self._spectrum_width))
        self._curve.setDownsampling(auto=True, method="peak")
        self._curve.setClipToView(True)
        self.plot_widget.scene().sigMouseMoved.connect(
            self._on_mouse_moved)
        layout.addWidget(self.plot_widget, stretch=1)

        # ---- Hover status line ----
        self.status_label = QtWidgets.QLabel(" ")
        self.status_label.setStyleSheet("color: gray;")
        layout.addWidget(self.status_label)

        self.gui_widget = window

    # -- File handling ----------------------------------------------------

    def _choose_file(self):
        last = settings_api.get_value(_KEY_LAST_FOLDER, "")
        path, _ = QtWidgets.QFileDialog.getOpenFileName(
            self.gui_widget, "Open FIES spectrum", last or "",
            "FITS spectra (*.fits *.fit *.fts *.gz);;All files (*)")
        if not path:
            return
        if self._load_file(path):
            settings_api.set_value(_KEY_LAST_FOLDER,
                                   os.path.dirname(path))

    def _try_autoload_current_frame(self):
        """Show the frame currently loaded in AMPA when it is a FIES
        1-D merged product."""
        path = ui_api.get_current_frame_path()
        if not path:
            return
        try:
            with fits.open(path, memmap=False) as hdul:
                header = hdul[0].header
                data = hdul[0].data
        except Exception:
            return
        if is_fies_spectrum(header, data):
            self._load_file(path)

    def _load_file(self, path: str) -> bool:
        try:
            flux, wavelength, header = read_spectrum(path)
        except Exception as exc:
            log(__name__, __LOGMODULE__, "warning",
                f"Could not load {path}: {exc}")
            ui_api.show_error_dialog(
                "FIES Spectrum Viewer", str(exc))
            return False

        self._path = path
        self._flux = flux
        self._wavelength = wavelength
        self._header = header

        name = os.path.basename(path)
        instrume = str(header.get("INSTRUME", "") or "").strip()
        if instrume.upper() != "FIES":
            shown = instrume or "unknown instrument"
            self.file_label.setText(f"{name}  ({shown})")
        else:
            self.file_label.setText(name)
        self.info_label.setText(format_info_line(header))
        self.gui_widget.setWindowTitle(f"FIES Spectrum Viewer - {name}")

        unit = str(header.get("BUNIT", "") or "").strip() or "Flux"
        self._plot_item.setLabel("left", unit)
        self._curve.setData(wavelength, flux)
        self._plot_item.autoRange()
        self._apply_telluric_overlay()
        log(__name__, __LOGMODULE__, "info",
            f"Loaded FIES spectrum {name} "
            f"({wavelength[0]:.1f}-{wavelength[-1]:.1f} A)")
        return True

    # -- Telluric overlay -------------------------------------------------

    def _on_overlay_toggled(self, checked: bool):
        self.labels_check.setEnabled(checked)
        self.width_spin.setEnabled(checked)
        self.legend_o2.setVisible(checked)
        self.legend_h2o.setVisible(checked)
        self._apply_telluric_overlay()

    def _on_labels_toggled(self, checked: bool):
        self._apply_telluric_overlay()

    def _on_line_width_changed(self, value: float):
        self._line_width = float(value)
        settings_api.set_value(_KEY_LINE_WIDTH, self._line_width)
        self._apply_telluric_overlay()

    def _on_spectrum_width_changed(self, value: float):
        self._spectrum_width = float(value)
        settings_api.set_value(_KEY_SPECTRUM_WIDTH,
                                self._spectrum_width)
        self._curve.setPen(
            pg.mkPen(_SPECTRUM_COLOR, width=self._spectrum_width))

    def _clear_telluric_overlay(self):
        for item in self._telluric_items:
            self._plot_item.removeItem(item)
        self._telluric_items = []
        self._labelled_air = []

    def _apply_telluric_overlay(self):
        """Rebuild the telluric lines for the current spectrum."""
        self._clear_telluric_overlay()
        if not self.telluric_check.isChecked() or self._wavelength is None:
            return
        lo, hi = sorted((float(self._wavelength[0]),
                         float(self._wavelength[-1])))
        show_labels = self.labels_check.isChecked()
        labelled: List[float] = []

        for air, vac, species, band, depth in TELLURIC_LINES:
            if not lo - 2.0 <= air <= hi + 2.0:
                continue
            colour = _SPECIES_COLORS.get(species, "#dddddd")
            label_text = None
            if show_labels and self._wants_label(air, species, depth,
                                                 labelled):
                label_text = f"{band} {air:.1f}"
                labelled.append(air)
            line = pg.InfiniteLine(
                pos=air, angle=90, movable=False,
                pen=pg.mkPen(colour, width=self._line_width),
                label=label_text,
                labelOpts={"position": 0.92, "color": colour})
            line.setToolTip(
                f"{band}\n{air:.3f} A (air)\n{vac:.3f} A (vacuum)\n"
                f"indicative depth {depth:.2f}")
            self._plot_item.addItem(line)
            self._telluric_items.append(line)

        self._labelled_air = labelled

    @staticmethod
    def _wants_label(air: float, species: str, depth: float,
                     labelled: List[float]) -> bool:
        """Label policy: every sufficiently deep O2 line keeps its band
        name (the bands are the landmark features), the strongest H2O
        lines get labelled with a minimum separation to avoid clutter."""
        if species == "O2":
            if depth < _LABEL_DEPTH_O2:
                return False
        else:
            if depth < _LABEL_DEPTH_H2O:
                return False
            if len(labelled) >= _LABEL_MAX_H2O:
                return False
        return all(abs(air - other) >= _LABEL_MIN_SEPARATION_A
                   for other in labelled)

    # -- Hover readout ------------------------------------------------------

    def _on_mouse_moved(self, scene_pos):
        if self._wavelength is None:
            return
        mouse_pt = self.plot_widget.mapFromScene(scene_pos)
        if not self.plot_widget.rect().contains(mouse_pt):
            return
        view_pos = self._plot_item.vb.mapSceneToView(scene_pos)
        x, y = float(view_pos.x()), float(view_pos.y())

        # nearest flux point for a value readout
        index = int(np.clip(np.searchsorted(self._wavelength, x), 0,
                            len(self._wavelength) - 1))
        text = (f"\u03bb = {self._wavelength[index]:.3f} A   "
                f"flux = {self._flux[index]:.4g}")

        if self.telluric_check.isChecked():
            near = self._nearest_telluric(x)
            if near is not None:
                air, vac, species, band, depth = near
                text += (f"   |   Telluric {band}: "
                         f"{air:.3f} A air / {vac:.3f} A vacuum, "
                         f"depth \u2248 {depth:.2f}")
        self.status_label.setText(text)

    def _nearest_telluric(self, x: float):
        """Closest overlay line within the hover radius, or None."""
        if len(_TELLURIC_AIR) == 0:
            return None
        index = int(np.clip(np.searchsorted(_TELLURIC_AIR, x), 0,
                            len(_TELLURIC_AIR) - 1))
        # pixels per Angstrom for the hover radius conversion
        view_range = self._plot_item.vb.viewRange()[0]
        view_width = self.plot_widget.width()
        if view_width <= 0 or view_range[1] <= view_range[0]:
            return None
        radius = (_HOVER_RADIUS_PX * (view_range[1] - view_range[0])
                  / view_width)

        best = None
        best_distance = radius
        for offset in (-1, 0, 1):
            i = index + offset
            if 0 <= i < len(TELLURIC_LINES):
                entry = TELLURIC_LINES[i]
                distance = abs(entry[0] - x)
                if distance <= best_distance:
                    best = entry
                    best_distance = distance
        return best
