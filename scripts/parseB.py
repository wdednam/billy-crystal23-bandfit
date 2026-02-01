#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Enhanced CRYSTAL *.BAND* parser
--------------------------------
 • extracts the band-energy cube and the reduced-k abscissa
 • also captures the high-symmetry tick positions/labels (Γ,X,…) that
   CRYSTAL writes into metadata and stores them in two extra *.npy* files
   so that other tools (dong_cost) can verify that reference & candidate follow
   identical k-paths.

Usage
-----
    parse_crystal_band.py  cell.BAND  [out_stub]

If out_stub is omitted the default is cand.npy (matching billy's call).
The following files are produced:

    <stub>.npy          - bands        (n_b x n_k [x n_spin])   [eV]
    <stub>.kdist.npy    - k-path abscissa                (n_k,)
    <stub>.hs_k.npy     - subset of kdist for HS points  (n_hs,)
    <stub>.hs_lbl.npy   - matching labels ("Γ","X",…) (n_hs,)

"""
from __future__ import annotations
import re, sys
from pathlib import Path
import numpy as np

_HA2EV          = 27.2114
_NUM_RE         = re.compile(r"[+-]?(?:\d*\.)?\d+(?:[Ee][+-]?\d+)?")
_SPIN_RE        = re.compile(r"\bSPIN\s+\d", re.I)
_HEADER_RE      = re.compile(r"\bBAND", re.I)
_TICK_POS_RE    = re.compile(r"@\s+XAXIS\s+TICK\s+(\d+)\s*,\s*([+-]?\d*\.?\d+)")
_TICK_LBL_RE    = re.compile(r"@\s+XAXIS\s+TICKLABEL\s+(\d+)\s*,\s*\"([^\"]+)\"")

def _extract_ticks(lines: list[str]):
    """Return (hs_k, hs_lbl) from metadata or empties if absent."""
    pos, lab = {}, {}
    for L in lines:
        if m := _TICK_POS_RE.search(L):
            pos[int(m[1])] = float(m[2])
        elif m := _TICK_LBL_RE.search(L):
            lab[int(m[1])] = m[2].strip()

    if not pos:
        return np.empty(0), np.empty(0, dtype="U1")

    order   = sorted(pos)
    hs_k    = np.array([pos[i]             for i in order], dtype=float)
    hs_lbl  = np.array([lab.get(i, f"k{i}") for i in order], dtype="U32")
    return hs_k, hs_lbl

def read_band(path: str | Path):
    """Return kdist, bands, hs_k, hs_lbl (see doc-string)."""
    lines   = Path(path).read_text(encoding="utf-8", errors="ignore").splitlines()
    hs_k, hs_lbl = _extract_ticks(lines)

    kdist, blocks, cur = [], [], []
    for L in lines:
        if _HEADER_RE.search(L):
            continue
        if _SPIN_RE.search(L):
            if cur:
                blocks.append(cur);  cur = []
            continue
        if not L.strip():
            continue

        parts = L.split()
        if len(parts) < 3:
            continue

        try:
            k_val = next(float(tok) for tok in parts if _NUM_RE.fullmatch(tok))
        except StopIteration:
            continue

        tail   = L.split(")")[-1] if ")" in L else " ".join(parts[1:])
        e_vals = [float(tok) for tok in tail.split() if _NUM_RE.fullmatch(tok)]
        if not e_vals:
            continue

        if kdist and abs(k_val - kdist[-1]) < 1e-8:   # CRYSTAL duplicates ends
            continue

        kdist.append(k_val)
        cur.append(np.sort(e_vals))

    if cur:
        blocks.append(cur)

    n_spin = len(blocks)
    n_k    = len(blocks[0])
    n_b    = max(len(row) for blk in blocks for row in blk)
    bands  = np.full((n_b, n_k, n_spin), np.nan)

    for s, blk in enumerate(blocks):
        for i, row in enumerate(blk):
            bands[:len(row), i, s] = row

    bands *= _HA2EV                                   # Ha → eV
    return np.asarray(kdist), bands.squeeze(), hs_k, hs_lbl

def _cli():
    if len(sys.argv) < 2:
        sys.exit(__doc__)

    in_file = Path(sys.argv[1])
    stub    = Path(sys.argv[2]) if len(sys.argv) > 2 else Path("cand.npy")
    k, e, hs_k, hs_lbl = read_band(in_file)

    np.save(stub,                      e,       allow_pickle=False)
    np.save(stub.with_suffix(".kdist.npy"), k,  allow_pickle=False)
    np.save(stub.with_suffix(".hs_k.npy"),  hs_k,  allow_pickle=False)
    np.save(stub.with_suffix(".hs_lbl.npy"),hs_lbl,allow_pickle=False)

    print(f"[parseB] → {stub}, kdist, hs_k, hs_lbl (bands {e.shape})")


if __name__ == "__main__":
    _cli()
