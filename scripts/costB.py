#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
costB  -  RMS mis-fit metric (eV) with branch-wise k-grid alignment
          and robust edge handling inside a FIXED energy window.

USAGE
    costB  ref.npy  cand.npy  [options]

OPTIONS
    --nband N           use the N lowest bands that both sets share
                        (ignored if --erange is given)
    --erange Emin Emax  fixed energy window (eV); recommended e.g. -10 10
    --align Dgap        rigidly shift candidate by half the direct band gap (after window)
"""

from __future__ import annotations
import sys, argparse
from pathlib import Path
import numpy as np


# ----------------- IO helpers ------------------------------------------------

def load_bands(path: Path) -> tuple[np.ndarray, np.ndarray]:
    bands = np.load(path, allow_pickle=True)
    if bands.dtype == object:
        bands = np.vstack(bands).T
    kdist = path.with_suffix(".kdist.npy")
    if not kdist.exists():
        sys.exit(f"[costB] missing k-grid file {kdist!s}")
    kpts = np.load(kdist)

    # enforce (Nb, Nk) shape
    if bands.shape[1] != kpts.size and bands.shape[0] == kpts.size:
        bands = bands.T
    return bands.astype(float), kpts.astype(float)


# ----------------- k-path / branches ----------------------------------------

def branch_bounds(k_axis: np.ndarray, hs_k: np.ndarray) -> np.ndarray:
    """Indices that split k_axis at hs_k (last slice closed to the end)."""
    hs_k = np.asarray(hs_k).astype(float).ravel()
    if hs_k.size < 2:
        return np.array([0, k_axis.size], dtype=int)
    idx = np.searchsorted(k_axis, hs_k, side="left")
    idx[0] = 0
    idx[-1] = k_axis.size
    idx = np.clip(idx, 0, k_axis.size)
    return idx


def interp_to(target_k: np.ndarray, src_k: np.ndarray, src_bands: np.ndarray) -> np.ndarray:
    """
    NaN-safe 1D interpolation band-wise, NO EXTRAPOLATION.
    Returns (Nb, target_k.size) with NaNs outside the src_k span or when
    not enough finite points exist in a band.
    """
    out = np.full((src_bands.shape[0], target_k.size), np.nan, float)
    if src_k.size < 2:
        return out
    # Only interpolate on the overlap domain
    m = (target_k >= src_k[0] - 1e-12) & (target_k <= src_k[-1] + 1e-12)
    if not np.any(m):
        return out
    tk = target_k[m]
    for b in range(src_bands.shape[0]):
        y = np.asarray(src_bands[b], dtype=float)
        finite = np.isfinite(y)
        if finite.sum() >= 2:
            out[b, m] = np.interp(tk, src_k[finite], y[finite])
    return out


# ----------------- metrics ---------------------------------------------------

def rms(a: np.ndarray, b: np.ndarray) -> float:
    mask = np.isfinite(a) & np.isfinite(b)
    if not mask.any():
        return np.inf
    with np.errstate(all="ignore"):
        return float(np.sqrt(np.nanmean((a[mask] - b[mask])**2)))


def apply_window(a: np.ndarray, b: np.ndarray, lo: float, hi: float) -> tuple[np.ndarray, np.ndarray]:
    """
    Independently mask out-of-window entries in each array.
    Works even if a and b have different numbers of bands (rows).
    """
    a = a.copy(); b = b.copy()
    a[(a < lo) | (a > hi)] = np.nan
    b[(b < lo) | (b > hi)] = np.nan
    return a, b


def edge_weights(E: np.ndarray, Emin: float, Emax: float, delta: float = 1.0) -> np.ndarray:
    """
    Smooth taper near window edges to avoid de-emphasizing boundary bands.
    Linear taper of width 'delta' (eV) at both ends; weight=1 in the core.
    """
    w = np.ones_like(E, float)
    if delta <= 0:
        return w
    core_lo = Emin + delta
    core_hi = Emax - delta
    m_lo = E < core_lo
    m_hi = E > core_hi
    # linear ramps
    w[m_lo] = np.clip((E[m_lo] - Emin) / delta, 0.0, 1.0)
    w[m_hi] = np.clip((Emax - E[m_hi]) / delta, 0.0, 1.0)
    return w


def weighted_rms_sorted(ref_b: np.ndarray, can_b: np.ndarray,
                        Emin: float, Emax: float,
                        delta: float = 1.0) -> float:
    """
    Per-k comparison of energies inside the window:
    - sort both sets (ascending),
    - match by order statistic,
    - compute a weighted RMS using energy-based edge weights.

    This stays index-agnostic (robust to crossings) but stops ignoring edges.
    """
    Nb, Nk = ref_b.shape
    ssq = 0.0
    wsum = 0.0
    for j in range(Nk):
        r = ref_b[:, j]
        c = can_b[:, j]
        r = r[np.isfinite(r)]
        c = c[np.isfinite(c)]
        if r.size == 0 or c.size == 0:
            continue
        wr = edge_weights(r, Emin, Emax, delta)
        wc = edge_weights(c, Emin, Emax, delta)
        ir = np.argsort(r)
        ic = np.argsort(c)
        m = min(ir.size, ic.size)
        if m == 0:
            continue
        dr = r[ir[:m]] - c[ic[:m]]
        wj = 0.5 * (wr[ir[:m]] + wc[ic[:m]])  # symmetric weight
        ssq += float(np.dot(wj, dr * dr))
        wsum += float(np.sum(wj))
    if wsum == 0.0:
        return np.inf
    return float(np.sqrt(ssq / wsum))


# ----------------- core ------------------------------------------------------

def common_overlap(a: np.ndarray, b: np.ndarray) -> tuple[float, float]:
    """Return the common [lo, hi] interval covered by both monotone arrays a and b."""
    lo = max(a[0], b[0])
    hi = min(a[-1], b[-1])
    return lo, hi


def one_rms(ref_b0, can_b0, ref_k, can_k,
            erange: tuple[float, float] | None,
            align: float | None,
            taper_delta: float = 1.0) -> float:
    ref_b = ref_b0.copy()
    can_b = can_b0.copy()

    # high-symmetry breaks from the candidate (like the plotting script)
    hs_k_path = Path("cand.hs_k.npy")
    if hs_k_path.exists():
        hs_k = np.load(hs_k_path).astype(float).ravel()
    else:
        sys.exit(f"[costB] missing k path {hs_k_path!s}")

    ref_k = np.asarray(ref_k, float).ravel()
    can_k = np.asarray(can_k, float).ravel()

    r_seg = branch_bounds(ref_k, hs_k)
    c_seg = branch_bounds(can_k, hs_k)

    ref_parts, can_parts = [], []

    # walk branches ΓX, XW, WL, ... exactly like the plot
    for (rs, re), (cs, ce) in zip(zip(r_seg[:-1], r_seg[1:]),
                                  zip(c_seg[:-1], c_seg[1:])):
        if re <= rs or ce <= cs:
            continue
        r_k, c_k = ref_k[rs:re], can_k[cs:ce]
        r_b, c_b = ref_b[:, rs:re], can_b[:, cs:ce]

        # enforce symmetric branch overlap (avoid dangling edges)
        lo, hi = common_overlap(r_k, c_k)
        if hi <= lo:
            continue
        ri = (r_k >= lo - 1e-12) & (r_k <= hi + 1e-12)
        ci = (c_k >= lo - 1e-12) & (c_k <= hi + 1e-12)
        r_k2, c_k2 = r_k[ri], c_k[ci]
        r_b2, c_b2 = r_b[:, ri], c_b[:, ci]

        # choose the denser grid as target, NO extrapolation
        if r_k2.size >= c_k2.size:
            c_on_r = interp_to(r_k2, c_k2, c_b2)
            r_on_r = r_b2
            if not np.isfinite(c_on_r).any():
                continue
            ref_parts.append(r_on_r)
            can_parts.append(c_on_r)
        else:
            r_on_c = interp_to(c_k2, r_k2, r_b2)
            c_on_c = c_b2
            if not np.isfinite(r_on_c).any():
                continue
            ref_parts.append(r_on_c)
            can_parts.append(c_on_c)

    if not ref_parts:
        return np.nan

    ref_b = np.concatenate(ref_parts, axis=1)
    can_b = np.concatenate(can_parts, axis=1)

    # alignment (after symmetric trimming); require reasonable overlap
    if align is not None and abs(align) > 0.0:
        ref_has = np.isfinite(ref_b).any(axis=0)
        can_has = np.isfinite(can_b).any(axis=0)
        overlap_k = int(np.sum(ref_has & can_has))
        if overlap_k < 10:
            return np.nan
        can_b = can_b - 0.5 * align

    if erange is not None:
        Emin, Emax = erange
        ref_b, can_b = apply_window(ref_b, can_b, Emin, Emax)
        # weighted, index-agnostic RMS to keep edges relevant
        return weighted_rms_sorted(ref_b, can_b, Emin, Emax, delta=taper_delta)
    else:
        # no window: plain RMS on concatenated branches
        return rms(ref_b, can_b)


# ----------------- CLI -------------------------------------------------------

def main() -> None:
    p = argparse.ArgumentParser(add_help=False)
    p.add_argument("ref")
    p.add_argument("cand")
    p.add_argument("--nband", type=int)
    p.add_argument("--erange", nargs=2, type=float, metavar=("Emin", "Emax"))
    p.add_argument("--align", type=float, default=0.0)
    args = p.parse_args()

    ref_b, ref_k = load_bands(Path(args.ref))
    can_b, can_k = load_bands(Path(args.cand))

    nband = args.nband or min(ref_b.shape[0], can_b.shape[0])

    using_window = (args.erange is not None)

    # If a user-specified energy window is used (--erange),
    # keep *all* bands so selection is driven purely by the window.
    if not using_window:
        ref_b = ref_b[:nband]
        can_b = can_b[:nband]

    if args.erange:
        er = (args.erange[0], args.erange[1])
        r = one_rms(ref_b, can_b, ref_k, can_k, erange=er, align=args.align)
        print(f"{r:.8f}")
        return
    else:
        r = one_rms(ref_b, can_b, ref_k, can_k, erange=None, align=args.align)
        print(f"{r:.8f}")
        return


if __name__ == "__main__":
    main()
