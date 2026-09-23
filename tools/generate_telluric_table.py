#!/usr/bin/env python3
"""Generate the embedded telluric line table for the FIES Spectrum Viewer.

Downloads synthetic telluric transmission spectra from the DKIST Telluric
Atlas (https://github.com/tschad/dkist_telluric_atlas, Zenodo DOI
10.5281/zenodo.13922154).  The atlas files are LBLRTM/Py4CAtS calculations
built on the HITRAN2020 database (Gordon et al. 2022, JQSRT 277, 107949)
for the US standard atmosphere above a 3 km base height (DKIST site).

The molecular identity of each absorption feature is derived *from the
atlas itself*: two transmissions that differ only in precipitable water
vapor (PWV 1 mm vs PWV 10 mm, same airmass 1.5) are compared.  Features
that deepen strongly with PWV are H2O; features insensitive to PWV are
O2 (the only other significant absorber between 5800 and 9000 A in air).

Output wavelengths are stored in *vacuum* Angstroms (the native HITRAN
convention), converted from the atlas air grid with the inverted
Edlen/Morton dispersion relation.  The plugin converts them back to air
at runtime.

Run from the repository root:

    python3 tools/generate_telluric_table.py [--cache DIR] [--out FILE]

and paste the generated block between the ``_TELLURIC_BEGIN`` /
``_TELLURIC_END`` markers in ``fies_spectrum_viewer.py``.
"""

from __future__ import annotations

import argparse
import os
import sys
import urllib.request

import numpy as np

# -- Source files (DKIST Telluric Atlas, version v20260407) --------------

_BASE = ("https://raw.githubusercontent.com/tschad/dkist_telluric_atlas/"
        "main/atlases/")
_WV = "telluric_atlas_mainMol_USstd_wv_air_angstrom_v20260407.npy"
# Reference transmission: typical Nordic Optical Telescope conditions
# (La Palma, 2.4 km, median PWV ~ 3 mm, median airmass ~ 1.5).
_TR_REF = ("telluric_atlas_mainMol_USstd_CO2_416ppm-Base_3km-"
           "PWV_3__mm-Airmass_1.5_v20260407.npy")

# PWV pair used for the H2O/O2 discrimination.
_TR_DRY = ("telluric_atlas_mainMol_USstd_CO2_416ppm-Base_3km-"
           "PWV_1__mm-Airmass_1.5_v20260407.npy")
_TR_WET = ("telluric_atlas_mainMol_USstd_CO2_416ppm-Base_3km-"
           "PWV_10_mm-Airmass_1.5_v20260407.npy")

# -- Extraction parameters -------------------------------------------------

# Wavelength window searched for lines (air Angstroms, atlas grid domain).
_LO_AIR, _HI_AIR = 5800.0, 9000.0

# Minimum absorption depth (1 - T) of a line at the reference PWV.
_MIN_DEPTH = 0.15

# Minima closer than this (air Angstroms) are treated as one blended
# feature (FIES HiRes FWHM is ~0.11 A at 7600 A).
_BLEND_A = 0.20

# A feature must deepen by more than this (absolute transmission change)
# between PWV 10 mm and PWV 1 mm to count as H2O.  Only the water vapour
# column differs between the two files, so O2 features are bit-identical
# in both; any real deepening indicates water.
_H2O_DEEPENING = 0.05

# O2 band windows (air Angstroms) used to name the electronic transition.
# The windows cover the full R- and P-branch extent of each band.
_O2_BANDS = (
    ("O2 gamma-band", 6250.0, 6340.0),
    ("O2 B-band", 6840.0, 6960.0),
    ("O2 A-band", 7580.0, 7720.0),
)

# Number of decimal places for the emitted wavelengths.
_WL_DECIMALS = 3


# -- Air <-> vacuum conversion (Edlen/Morton) -----------------------------

def air_refractive_index(wl: "np.ndarray | float") -> "np.ndarray | float":
    """Refractive index of air at *wl* (Angstroms).

    Morton (1991, ApJS 77, 119) dispersion relation:
    n - 1 = 2.7357e-4 + 131.4182 / wl^2 + 2.762e8 / wl^4
    valid for the optical; accurate to ~1e-8 in n.
    """
    return 1.0 + 2.7357e-4 + 131.4182 / np.asarray(wl) ** 2 \
        + 2.762e8 / np.asarray(wl) ** 4


def air_to_vacuum(wl_air):
    """Convert air wavelengths (A) to vacuum wavelengths (A)."""
    wl_air = np.asarray(wl_air, dtype=float)
    wl_vac = wl_air * air_refractive_index(wl_air)
    for _ in range(3):  # converges to well below 1e-6 A
        wl_vac = wl_air * air_refractive_index(wl_vac)
    return wl_vac


# -- Helpers ---------------------------------------------------------------

def _download(url: str, dest: str) -> str:
    if not os.path.exists(dest):
        print(f"downloading {url}")
        urllib.request.urlretrieve(url, dest)
    return dest


def _load_cached(cache: str, name: str) -> np.ndarray:
    dest = os.path.join(cache, name)
    _download(_BASE + name, dest)
    return np.load(dest, allow_pickle=False)


def _local_minima(t: np.ndarray) -> np.ndarray:
    """Indices of strict local minima (window +-2 samples)."""
    n = len(t)
    idx = np.arange(2, n - 2)
    mask = (
        (t[idx] < t[idx - 1]) & (t[idx] < t[idx + 1])
        & (t[idx] <= t[idx - 2]) & (t[idx] <= t[idx + 2])
    )
    return idx[mask]


def _extract_features(wl: np.ndarray, t_ref: np.ndarray) -> list:
    """Return blended absorption features as (wl, depth) pairs."""
    m = (wl >= _LO_AIR) & (wl <= _HI_AIR)
    w, t = wl[m], t_ref[m]
    mins = _local_minima(t)
    feats = []
    for i in mins:
        depth = 1.0 - float(t[i])
        if depth < _MIN_DEPTH:
            continue
        pos = float(w[i])
        if feats and pos - feats[-1][0] < _BLEND_A:
            if depth > feats[-1][1]:
                feats[-1] = (pos, depth)
        else:
            feats.append((pos, depth))
    return feats


def _classify(feats: list, wl: np.ndarray, t_dry: np.ndarray,
              t_wet: np.ndarray) -> list:
    """Split features into H2O / O2 from their PWV response."""
    out = []
    for pos, depth in feats:
        j = int(np.searchsorted(wl, pos))
        j = min(max(j, 2), len(wl) - 3)
        # neighbourhood minimum: line cores can shift by a sample between
        # the two PWV files
        d_dry = 1.0 - float(t_dry[j - 2:j + 3].min())
        d_wet = 1.0 - float(t_wet[j - 2:j + 3].min())
        species = "H2O" if (d_wet - d_dry) > _H2O_DEEPENING else "O2"
        out.append((pos, depth, species))
    return out


def _band_label(pos_air: float, species: str) -> str:
    if species == "H2O":
        return "H2O"
    for name, lo, hi in _O2_BANDS:
        if lo <= pos_air <= hi:
            return name
    return "O2"


# -- Main -------------------------------------------------------------------

def generate(cache: str) -> str:
    wl = _load_cached(cache, _WV)
    t_ref = _load_cached(cache, _TR_REF)
    t_dry = _load_cached(cache, _TR_DRY)
    t_wet = _load_cached(cache, _TR_WET)

    feats = _extract_features(wl, t_ref)
    classified = _classify(feats, wl, t_dry, t_wet)

    lines = []
    for pos_air, depth, species in sorted(classified):
        pos_vac = float(air_to_vacuum(pos_air))
        lines.append((round(pos_vac, _WL_DECIMALS), species,
                      _band_label(pos_air, species), round(depth, 2)))

    n_h2o = sum(1 for _, s, _, _ in lines if s == "H2O")
    n_o2 = len(lines) - n_h2o
    print(f"extracted {len(lines)} lines "
          f"(O2: {n_o2}, H2O: {n_h2o}) "
          f"with depth >= {_MIN_DEPTH} at PWV 3 mm, airmass 1.5")

    rows = []
    for pos_vac, species, label, depth in lines:
        rows.append(f"    ({pos_vac:.3f}, \"{species}\", \"{label}\", {depth:.2f}),")
    body = "\n".join(rows)

    return (
        "# Telluric absorption lines, sorted by wavelength.  Tuples of\n"
        "# (vacuum wavelength [A], species, band label, indicative depth\n"
        "# at PWV 3 mm / airmass 1.5 / 3 km base height).  Derived from the\n"
        "# DKIST Telluric Atlas (HITRAN2020 / LBLRTM); regenerate with\n"
        "# tools/generate_telluric_table.py in the plugin repository.\n"
        "_TELLURIC_LINES_VACUUM = [\n" + body + "\n]\n"
    )


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        description="Generate the embedded telluric line table.")
    parser.add_argument("--cache", default="/tmp/dkist_telluric_atlas",
                        help="download cache directory")
    parser.add_argument("--out", default=None,
                        help="output file (default: stdout)")
    args = parser.parse_args(argv)

    os.makedirs(args.cache, exist_ok=True)
    text = generate(args.cache)
    if args.out:
        with open(args.out, "w", encoding="utf-8") as fh:
            fh.write(text)
        print(f"wrote {args.out}")
    else:
        print(text)
    return 0


if __name__ == "__main__":
    sys.exit(main())
