#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
openmx_banddat1_to_npy.py
=========================
Parse an OpenMX BANDDAT1 file after the 56-line header and write

  • <out>.npy        →  band energies, shape (Nb , Nk)
  • <out>.kdist.npy  →  k-point cumulative distance, shape (Nk,)

The routine is orientation-agnostic - it simply groups consecutive
numeric lines that share the same k-distance (to 1 x 10⁻⁶ Å⁻¹).
"""

from pathlib import Path
from collections import OrderedDict
import math, sys, numpy as np

HEADER_LINES = 56          # lines to skip verbatim
TOL          = 1e-6        # k‑distance merge tolerance


def parse_banddat1(fname: Path, skip: int = HEADER_LINES, tol: float = TOL):
    """
    Returns
    -------
    k_dist  : (Nk,)    float64
    bands   : (Nb,Nk)  float64     (shorter columns padded with nan)
    """
    buckets = OrderedDict()          # k_distance -> list[eivals]

    with fname.open("r", encoding="latin1", errors="ignore") as fh:
        # skip the fixed header
        for _ in range(skip):
            next(fh, None)

        # bucket numeric lines by (almost) equal k
        for raw in fh:
            line = raw.strip()
            if not line:
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            try:
                k_val = float(parts[0])
                e_val = float(parts[1])
            except ValueError:
                continue

            # find an existing key within tolerance
            for key in buckets:
                if abs(key - k_val) < tol:
                    buckets[key].append(e_val)
                    break
            else:
                buckets[k_val] = [e_val]

    if not buckets:
        raise RuntimeError("No numeric data recognised – wrong file?")

    # convert to arrays, pad with nan where a band is missing
    Nk      = len(buckets)
    Nb      = max(len(lst) for lst in buckets.values())
    k_dist  = np.fromiter(buckets.keys(), float, Nk)
    bands   = np.full((Nb, Nk), math.nan, float)

    for j, vals in enumerate(buckets.values()):
        bands[:len(vals), j] = vals

    return k_dist, bands


def main(argv):
    if len(argv) != 3:
        print("Usage:\n  python openmx_banddat1_to_npy.py  BANDDAT1  out.npy")
        sys.exit(1)

    infile  = Path(argv[1]).expanduser()
    outfile = Path(argv[2]).expanduser()

    kdist, bands = parse_banddat1(infile)

    np.save(outfile, bands)
    np.save(outfile.with_suffix(".kdist.npy"), kdist)

    nb, nk = bands.shape
    print(f"Parsed: Nk={nk}, Nb={nb} → saved '{outfile}' and '.kdist.npy'")
    if np.isnan(bands).any():
        print("  (columns padded with NaNs where a band was missing)")


if __name__ == "__main__":
    main(sys.argv)
