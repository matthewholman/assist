#!/usr/bin/env python3
"""
Compare raw Chebyshev coefficient blocks between:
  - JPL ASCII-derived binary ephemeris (.440/.441)
  - NAIF SPK kernel (.bsp)

This script aligns coefficient *records* (Chebyshev blocks) and computes:
  - fraction of coefficients that are bit-identical
  - max abs / relative difference
  - max ULP difference + a coarse histogram of ULP differences
  - mantissa trailing-zero statistics for each format (proxy for quantization)

It is designed specifically for DE440/DE441 style files bundled with this repo:
  - data/linux_p1550p2650.440
  - data/de440.bsp
but should work for DE441 equivalents with the same layout.
"""

from __future__ import annotations

import argparse
import ctypes
import json
import os
import struct
import sysconfig
from dataclasses import asdict, dataclass
from typing import Dict, List, Optional, Tuple

import numpy as np
from numpy.lib.stride_tricks import as_strided


# -----------------------------
# Helpers: .440 header parsing
# -----------------------------


@dataclass(frozen=True)
class JplHeader:
    beg: float
    end: float
    inc_days: float
    ver: int
    num_constants: int
    cau_km: float
    cem: float
    # per-component tables (length 15)
    off: List[int]  # zero-based offset in double-words within a block
    ncf: List[int]
    niv: List[int]
    ncm: List[int]
    rec_doubles: int  # number of doubles per block record


def parse_jpl_header(path: str) -> JplHeader:
    """
    Parse enough of the ASCII-derived binary header to locate coefficient blocks.
    Mirrors the layout assumptions in `src/ascii_ephem.c` (`assist_ascii_init`).
    """

    JPL_N = 15
    JPL_NUT = 11
    JPL_TDB = 14

    with open(path, "rb") as f:
        f.seek(0x0A5C)
        beg, end, inc = struct.unpack("<ddd", f.read(24))
        num = struct.unpack("<i", f.read(4))[0]
        cau = struct.unpack("<d", f.read(8))[0]
        cem = struct.unpack("<d", f.read(8))[0]

        off: List[int] = []
        ncf: List[int] = []
        niv: List[int] = []

        # columns 1..12
        for _ in range(12):
            off.append(struct.unpack("<i", f.read(4))[0])
            ncf.append(struct.unpack("<i", f.read(4))[0])
            niv.append(struct.unpack("<i", f.read(4))[0])

        ver = struct.unpack("<i", f.read(4))[0]
        # column 13
        off.append(struct.unpack("<i", f.read(4))[0])
        ncf.append(struct.unpack("<i", f.read(4))[0])
        niv.append(struct.unpack("<i", f.read(4))[0])

        # Skip constant names, as in C:
        # first 400 names at 0x00FC, remaining names at 0x0B28.
        # After reading remaining names, the file position is where cols 14/15 live.
        f.seek(0x00FC)
        _ = f.read(400 * 6)
        f.seek(0x0B28)
        _ = f.read((num - 400) * 6)

        # columns 14..15
        for _ in range(2):
            off.append(struct.unpack("<i", f.read(4))[0])
            ncf.append(struct.unpack("<i", f.read(4))[0])
            niv.append(struct.unpack("<i", f.read(4))[0])

    ncm = [3] * JPL_N
    ncm[JPL_NUT] = 2
    ncm[JPL_TDB] = 1

    # convert to zero-based offsets
    off = [x - 1 for x in off]

    # record size in doubles: 2 doubles + sum over components (ncf*niv*ncm)
    rec = 2
    for p in range(JPL_N):
        rec += ncf[p] * niv[p] * ncm[p]

    return JplHeader(
        beg=beg,
        end=end,
        inc_days=inc,
        ver=ver,
        num_constants=num,
        cau_km=cau,
        cem=cem,
        off=off,
        ncf=ncf,
        niv=niv,
        ncm=ncm,
        rec_doubles=rec,
    )


# -----------------------------
# Helpers: SPK access via libassist
# -----------------------------


class SpkTarget(ctypes.Structure):
    _fields_ = [
        ("code", ctypes.c_int),
        ("cen", ctypes.c_int),
        ("mass", ctypes.c_double),
        ("beg", ctypes.c_double),
        ("end", ctypes.c_double),
        ("res", ctypes.c_double),
        ("one", ctypes.POINTER(ctypes.c_int)),
        ("two", ctypes.POINTER(ctypes.c_int)),
        ("ind", ctypes.c_int),
    ]


class SpkS(ctypes.Structure):
    _fields_ = [
        ("targets", ctypes.POINTER(SpkTarget)),
        ("num", ctypes.c_int),
        ("allocated_num", ctypes.c_int),
        ("map", ctypes.c_void_p),
        ("len", ctypes.c_size_t),
    ]


@dataclass
class SpkContext:
    lib: ctypes.CDLL
    spk: ctypes.POINTER(SpkS)
    doubles: np.ndarray  # float64 view of mmap'd file

    def close(self) -> None:
        self.lib.assist_spk_free(self.spk)


def load_spk_with_libassist(libassist_path: str, spk_path: str) -> SpkContext:
    lib = ctypes.CDLL(libassist_path)
    lib.assist_spk_init.argtypes = [ctypes.c_char_p]
    lib.assist_spk_init.restype = ctypes.POINTER(SpkS)
    lib.assist_spk_free.argtypes = [ctypes.POINTER(SpkS)]
    lib.assist_spk_free.restype = ctypes.c_int

    spk = lib.assist_spk_init(spk_path.encode("ascii"))
    if not spk:
        raise RuntimeError(f"assist_spk_init failed for {spk_path!r}")

    spk_v = spk.contents
    n_d = spk_v.len // 8
    arr = np.ctypeslib.as_array((ctypes.c_double * n_d).from_address(spk_v.map))
    return SpkContext(lib=lib, spk=spk, doubles=arr)


def find_target(spk: SpkContext, code: int, cen: Optional[int] = None) -> SpkTarget:
    v = spk.spk.contents
    for i in range(v.num):
        t = v.targets[i]
        if t.code == code and (cen is None or t.cen == cen):
            return t
    raise KeyError(f"SPK target code={code} cen={cen} not found")


@dataclass
class SpkDirectory:
    R: int
    P: int
    intlen_days: float
    N_records: int
    one: int
    two: int


def read_spk_directory(spk: SpkContext, target: SpkTarget) -> SpkDirectory:
    """
    Read directory tail values used by assist_spk_calc:
      val = map + two - 1
      N = val[0]
      intlen = val[-2]
      R = val[-1]
    """
    if target.ind != 0:
        raise RuntimeError(
            f"Expected single segment (ind==0) for code {target.code}, got {target.ind}"
        )
    one = int(target.one[0])
    two = int(target.two[0])
    dir_idx = two - 1
    R = int(spk.doubles[dir_idx - 1])
    intlen = float(spk.doubles[dir_idx - 2])  # seconds
    N = int(spk.doubles[dir_idx])  # number of records
    P = (R - 2) // 3
    return SpkDirectory(
        R=R, P=P, intlen_days=(intlen / 86400.0), N_records=N, one=one, two=two
    )


# -----------------------------
# Bit-level metrics
# -----------------------------


def mantissa_trailing_zeros(arr: np.ndarray) -> np.ndarray:
    """
    Return per-value trailing-zero count (0..52) of the mantissa bits.
    For mantissa==0 (exact power-of-two exponent with no fraction), returns 52.
    """
    bits = arr.view(np.uint64)
    mant = bits & np.uint64((1 << 52) - 1)
    # two's complement lowbit: mant & -mant
    lowbit = mant & (np.uint64(0) - mant)
    tz = np.where(mant == 0, 52, np.bitwise_count(lowbit - 1))
    return tz.astype(np.int64)


def ulp_distance(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """
    Compute ULP distance between float64 arrays a and b (same shape).
    Returns uint64 distances.
    """
    ua = a.view(np.uint64)
    ub = b.view(np.uint64)
    sign_a = ua >> np.uint64(63)
    sign_b = ub >> np.uint64(63)
    # ordered mapping: negatives invert, positives offset by 0x8000...
    oa = np.where(
        sign_a != 0,
        np.uint64(0xFFFFFFFFFFFFFFFF) - ua,
        ua + np.uint64(0x8000000000000000),
    )
    ob = np.where(
        sign_b != 0,
        np.uint64(0xFFFFFFFFFFFFFFFF) - ub,
        ub + np.uint64(0x8000000000000000),
    )
    return np.where(oa > ob, oa - ob, ob - oa)


# -----------------------------
# Comparison logic
# -----------------------------


@dataclass
class CompareStats:
    name: str
    spk_code: int
    spk_cen: int
    jpl_component: int
    scale: float
    coeffs_total: int = 0
    coeffs_equal: int = 0
    max_abs_diff: float = 0.0
    max_rel_diff: float = 0.0
    max_ulp_diff: int = 0
    ulp_hist: Dict[str, int] = None  # filled later
    tz_mean_spk: float = 0.0
    tz_mean_jpl: float = 0.0
    tz_ge_16_frac_spk: float = 0.0
    tz_ge_16_frac_jpl: float = 0.0


def update_max(a: float, b: float) -> float:
    return a if a >= b else b


def compare_component(
    spk: SpkContext,
    jpl_header: JplHeader,
    jpl_arr: np.ndarray,
    mapping_name: str,
    spk_code: int,
    spk_cen: int,
    jpl_component: int,
    scale: float = 1.0,
    max_blocks: Optional[int] = None,
) -> CompareStats:
    t = find_target(spk, spk_code, spk_cen)
    d = read_spk_directory(spk, t)

    ncf = jpl_header.ncf[jpl_component]
    niv = jpl_header.niv[jpl_component]
    off = jpl_header.off[jpl_component]
    ncm = jpl_header.ncm[jpl_component]
    inc = jpl_header.inc_days

    if d.P != ncf:
        raise RuntimeError(f"{mapping_name}: SPK P={d.P} != JPL ncf={ncf}")
    if abs((inc / d.intlen_days) - niv) > 1e-9:
        raise RuntimeError(
            f"{mapping_name}: SPK intlen_days={d.intlen_days} not consistent with "
            f"JPL inc/niv ({inc}/{niv})"
        )
    if ncm != 3:
        raise RuntimeError(
            f"{mapping_name}: unexpected ncm={ncm} (script currently expects 3-vector components)"
        )
    if d.N_records % niv != 0:
        raise RuntimeError(
            f"{mapping_name}: SPK N_records={d.N_records} not divisible by niv={niv}"
        )

    nblocks = d.N_records // niv
    if max_blocks is not None:
        nblocks = min(nblocks, max_blocks)

    stats = CompareStats(
        name=mapping_name,
        spk_code=spk_code,
        spk_cen=spk_cen,
        jpl_component=jpl_component,
        scale=scale,
        ulp_hist={},
    )

    # Histogram for ULP diffs (non-zero): exact counts for 1..8, plus a catch-all >8.
    # (In practice for DE440 we see max ULP <= ~5.)
    ulp_bins = np.zeros(10, dtype=np.int64)  # indices 0..9, where 9 is >8

    # Trailing-zero hist for mantissas (0..52)
    tz_hist_spk = np.zeros(53, dtype=np.int64)
    tz_hist_jpl = np.zeros(53, dtype=np.int64)

    for blk in range(nblocks):
        r0 = blk * niv
        # SPK block view: (niv, 3, P)
        spk_start = (d.one - 1) + r0 * d.R + 2
        spk_block = as_strided(
            spk.doubles[spk_start:],
            shape=(niv, 3, d.P),
            strides=(d.R * 8, d.P * 8, 8),
        )

        # JPL block view: (niv, 3, P)
        jpl_base = (blk + 2) * jpl_header.rec_doubles + off
        jpl_block = as_strided(
            jpl_arr[jpl_base:],
            shape=(niv, 3, d.P),
            strides=(3 * d.P * 8, d.P * 8, 8),
        )

        # Mantissa stats on raw blocks (pre-scale)
        tz_hist_spk += np.bincount(mantissa_trailing_zeros(spk_block).ravel(), minlength=53)
        tz_hist_jpl += np.bincount(mantissa_trailing_zeros(jpl_block).ravel(), minlength=53)

        jpl_cmp = jpl_block if scale == 1.0 else (jpl_block * np.float64(scale))

        eq = spk_block == jpl_cmp
        stats.coeffs_total += int(eq.size)
        stats.coeffs_equal += int(eq.sum())

        if (~eq).any():
            diff = spk_block - jpl_cmp
            absdiff = np.abs(diff)
            stats.max_abs_diff = update_max(stats.max_abs_diff, float(absdiff.max()))

            denom = np.maximum(np.abs(spk_block), np.float64(1e-300))
            reldiff = absdiff / denom
            stats.max_rel_diff = update_max(stats.max_rel_diff, float(reldiff.max()))

            # ULP diffs only on non-equal entries
            a = spk_block[~eq]
            b = jpl_cmp[~eq]
            u = ulp_distance(a, b)
            stats.max_ulp_diff = max(stats.max_ulp_diff, int(u.max()))

            u_non0 = u[u != 0].astype(np.uint64)
            if u_non0.size:
                capped = np.minimum(u_non0, np.uint64(9)).astype(np.int64)  # 9 means >8
                ulp_bins += np.bincount(capped, minlength=ulp_bins.size).astype(np.int64)

    # Persist histogram
    for i in range(1, 9):
        stats.ulp_hist[f"ulp=={i}"] = int(ulp_bins[i])
    stats.ulp_hist[">8"] = int(ulp_bins[9])

    # mantissa trailing-zero summaries
    stats.tz_mean_spk = float((tz_hist_spk * np.arange(53)).sum() / tz_hist_spk.sum())
    stats.tz_mean_jpl = float((tz_hist_jpl * np.arange(53)).sum() / tz_hist_jpl.sum())
    stats.tz_ge_16_frac_spk = float(tz_hist_spk[16:].sum() / tz_hist_spk.sum())
    stats.tz_ge_16_frac_jpl = float(tz_hist_jpl[16:].sum() / tz_hist_jpl.sum())

    return stats


def main(argv: Optional[List[str]] = None) -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--jpl", default="data/linux_p1550p2650.440", help="Path to ASCII-derived binary (.440/.441) file")
    p.add_argument("--spk", default="data/de440.bsp", help="Path to SPK .bsp file")
    p.add_argument("--libassist", default=None, help="Path to libassist shared library (defaults to ./libassist<EXT_SUFFIX>)")
    p.add_argument("--max-blocks", type=int, default=None, help="Limit number of 32-day blocks to analyze (debug)")
    p.add_argument("--json", default=None, help="Write JSON report to this path")
    args = p.parse_args(argv)

    suffix = sysconfig.get_config_var("EXT_SUFFIX") or ".so"
    libassist_path = args.libassist or os.path.join(os.getcwd(), "libassist" + suffix)
    if not os.path.exists(libassist_path):
        raise SystemExit(
            f"libassist not found at {libassist_path!r}. Build it first (python setup.py build_ext --inplace)."
        )

    jpl_header = parse_jpl_header(args.jpl)
    jpl_arr = np.memmap(args.jpl, dtype="<f8", mode="r")

    spk_ctx = load_spk_with_libassist(libassist_path, args.spk)
    try:
        # Mappings (SPK code, cen) -> JPL component index
        # JPL components: MER=0, VEN=1, EMB=2, MAR=3, JUP=4, SAT=5, URA=6, NEP=7, PLU=8, LUN=9, SUN=10
        mappings: List[Tuple[str, int, int, int, float]] = [
            ("Mercury barycentric", 1, 0, 0, 1.0),
            ("Venus barycentric", 2, 0, 1, 1.0),
            ("EMB barycentric", 3, 0, 2, 1.0),
            ("Mars barycentric", 4, 0, 3, 1.0),
            ("Jupiter barycentric", 5, 0, 4, 1.0),
            ("Saturn barycentric", 6, 0, 5, 1.0),
            ("Uranus barycentric", 7, 0, 6, 1.0),
            ("Neptune barycentric", 8, 0, 7, 1.0),
            ("Pluto barycentric", 9, 0, 8, 1.0),
            ("Sun barycentric", 10, 0, 10, 1.0),
        ]

        # Earth/Moon relative to EMB are stored in SPK (cen=3) but not directly in .440.
        # In .440, lunar geocentric vector LUN is used:
        #   Earth_rel_EMB = - LUN/(1+EMRAT)
        #   Moon_rel_EMB  =   LUN*EMRAT/(1+EMRAT)
        emrat = jpl_header.cem
        mappings += [
            ("Earth relative EMB (derived from .440 LUN)", 399, 3, 9, (-1.0 / (1.0 + emrat))),
            ("Moon relative EMB (derived from .440 LUN)", 301, 3, 9, (emrat / (1.0 + emrat))),
        ]

        results: List[CompareStats] = []
        for name, code, cen, comp, scale in mappings:
            results.append(
                compare_component(
                    spk=spk_ctx,
                    jpl_header=jpl_header,
                    jpl_arr=jpl_arr,
                    mapping_name=name,
                    spk_code=code,
                    spk_cen=cen,
                    jpl_component=comp,
                    scale=scale,
                    max_blocks=args.max_blocks,
                )
            )

        for s in results:
            frac = s.coeffs_equal / s.coeffs_total if s.coeffs_total else 0.0
            print(f"\n== {s.name} ==")
            print(f"coeffs: {s.coeffs_equal}/{s.coeffs_total} exact ({frac:.6f})")
            print(
                f"max abs diff: {s.max_abs_diff:.3e} | "
                f"max rel diff: {s.max_rel_diff:.3e} | "
                f"max ULP: {s.max_ulp_diff}"
            )
            print(
                f"mantissa tz mean: spk={s.tz_mean_spk:.2f} jpl={s.tz_mean_jpl:.2f} | "
                f"tz>=16 frac: spk={s.tz_ge_16_frac_spk:.3f} jpl={s.tz_ge_16_frac_jpl:.3f}"
            )
            nonzero = {k: v for k, v in s.ulp_hist.items() if v}
            if nonzero:
                print("ULP hist (nonzero bins):", nonzero)

        if args.json:
            payload = {
                "jpl": {"path": args.jpl, **asdict(jpl_header)},
                "spk": {"path": args.spk},
                "results": [asdict(r) for r in results],
            }
            with open(args.json, "w", encoding="utf-8") as f:
                json.dump(payload, f, indent=2, sort_keys=True)
            print(f"\nWrote JSON report to {args.json}")

        return 0
    finally:
        spk_ctx.close()


if __name__ == "__main__":
    raise SystemExit(main())


