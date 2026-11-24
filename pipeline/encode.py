"""Compact, fast-to-load binary files (gzip) for the browser.

File layout (little-endian), 64-byte header then raw arrays:

  offset  type   field
  0       4s     magic  b"GWX1"
  4       u8     kind   1 = cumulative totals, 2 = precipitation-type probabilities
  5       u8     narr    kind 1: number of cumulative u16 arrays that follow (4; files with 0 hold 2)
  6       u16    nmem   number of ensemble members (31)
  8       u16    ny     rows  (north -> south)
  10      u16    nx     cols  (west -> east)
  12      u16    fhr    forecast hour
  14      u16    step_h hours since the previous written step
  16      f32    lat_north (center of first row)
  20      f32    lon_west  (center of first column, degrees east, -180..180)
  24      f32    dlat
  28      f32    dlon
  32      f32    scale0  kind 1: mm of liquid per count (QPF)
  36      f32    scale1  kind 1: cm of snow per count
  40      f32    scale2  kind 1: mm of liquid per count, freezing rain
  44      f32    scale3  kind 1: mm of liquid per count, sleet / ice pellets
  48      u32    member bitmask, members 0-31: bit k set = member k (c00, p01 ...) has data in this file. 0 = legacy file, all present
  52      f32    scale4  kind 1: mm of liquid per count, snow liquid equivalent (any snow ratio can be applied to it)
  56      u32    member bitmask, members 32-63 (IFS has 50 members; 0 for the 31-member GEFS)
  60..63  zero padding

 kind 1 body:  u16 [nmem][ny][nx] cumulative QPF, snowfall, freezing rain (liquid eq.), sleet (liquid eq.), in that order.
               Any rolling window total is  cum(t) - cum(t - window)  -- two small files.
 kind 2 body:  u8  [4][nmem][ny][nx]   probabilities (0-100) of RA, SN, FZRA, PL.
               All zero where that member has no precipitation in the 3-h step.
"""
from __future__ import annotations

import gzip
import struct
from pathlib import Path

import numpy as np

from . import config as C

MAGIC = b"GWX1"
HEADER = struct.Struct("<4sBBHHHHHffffffffIfI")     # ... mask (members 0-31), the 5th scale (snow liquid equivalent), mask (members 32-63)
HEADER_SIZE = 64
KIND_CUM, KIND_PTYPE = 1, 2


def full_mask(nmem: int) -> int:
    return (1 << nmem) - 1


def mask_of(present_idx) -> int:
    m = 0
    for k in present_idx:
        m |= 1 << int(k)
    return m


def _header(kind, nmem, ny, nx, fhr, step_h, domain: C.Domain, scales=(0.0,) * 5, narr=0, mask=0) -> bytes:
    if nmem > 64:
        raise ValueError("the member bitmask holds at most 64 members")
    scales = tuple(scales) + (0.0,) * (5 - len(scales))
    raw = HEADER.pack(MAGIC, kind, narr, nmem, ny, nx, fhr, step_h,
                      domain.lat_north, domain.lon_west, domain.res, domain.res, *scales[:4], mask & 0xFFFFFFFF, scales[4],
                      (mask >> 32) & 0xFFFFFFFF)
    return raw + b"\x00" * (HEADER_SIZE - len(raw))


def _write_gz(path: Path, payload: bytes) -> int:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "wb") as raw:
        with gzip.GzipFile(fileobj=raw, mode="wb", compresslevel=9, mtime=0) as gz:
            gz.write(payload)
    return path.stat().st_size


def quantize(arr: np.ndarray, per_count: float) -> np.ndarray:
    return np.clip(np.rint(arr / per_count), 0, 65535).astype(np.uint16)


def write_cum(path: Path, fhr: int, step_h: int, qpf_cum_mm: np.ndarray, snow_cum_cm: np.ndarray,
              domain: C.Domain, zr_cum_mm: np.ndarray | None = None, pl_cum_mm: np.ndarray | None = None,
              mask: int | None = None, swe_cum_mm: np.ndarray | None = None) -> int:
    m, ny, nx = qpf_cum_mm.shape
    mask = full_mask(m) if mask is None else mask
    parts = [quantize(qpf_cum_mm, C.QPF_MM_PER_COUNT), quantize(snow_cum_cm, C.SNOW_CM_PER_COUNT)]
    scales = [C.QPF_MM_PER_COUNT, C.SNOW_CM_PER_COUNT, 0.0, 0.0, 0.0]
    if zr_cum_mm is not None and pl_cum_mm is not None:
        parts += [quantize(zr_cum_mm, C.ZR_MM_PER_COUNT), quantize(pl_cum_mm, C.PL_MM_PER_COUNT)]
        scales[2:4] = [C.ZR_MM_PER_COUNT, C.PL_MM_PER_COUNT]
        if swe_cum_mm is not None:
            parts.append(quantize(swe_cum_mm, C.SWE_MM_PER_COUNT))
            scales[4] = C.SWE_MM_PER_COUNT
    body = b"".join(a.astype("<u2").tobytes() for a in parts)
    return _write_gz(Path(path), _header(KIND_CUM, m, ny, nx, fhr, step_h, domain, tuple(scales), len(parts), mask) + body)


def write_ptype(path: Path, fhr: int, step_h: int, probs_u8: np.ndarray, domain: C.Domain, mask: int | None = None) -> int:
    """probs_u8: (4, nmem, ny, nx) Bourgouin percent [RA, SN, FZRA, PL]."""
    _, m, ny, nx = probs_u8.shape
    body = np.ascontiguousarray(probs_u8, dtype=np.uint8).tobytes()
    hdr = _header(KIND_PTYPE, m, ny, nx, fhr, step_h, domain, mask=full_mask(m) if mask is None else mask)
    return _write_gz(Path(path), hdr + body)


def read_bin(path: Path) -> dict:
    """Reference reader (used by tests; the browser has its own JS reader)."""
    raw = gzip.open(path, "rb").read()
    (magic, kind, narr, nmem, ny, nx, fhr, step_h, lat0, lon0, dlat, dlon, s0, s1, s2, s3, mask, s4, mask_hi) = HEADER.unpack_from(raw)
    mask |= mask_hi << 32
    if magic != MAGIC:
        raise ValueError("bad magic")
    n = nmem * ny * nx
    off = HEADER_SIZE
    out = {"kind": kind, "nmem": nmem, "ny": ny, "nx": nx, "fhr": fhr, "step_h": step_h,
           "lat0": lat0, "lon0": lon0, "dlat": dlat, "dlon": dlon, "mask": mask or full_mask(nmem)}
    if kind == KIND_CUM:
        narr = narr or 2
        arrs = [np.frombuffer(raw, "<u2", n, off + 2 * n * k).reshape(nmem, ny, nx) for k in range(narr)]
        out["qpf_mm"] = arrs[0].astype(np.float32) * s0
        out["snow_cm"] = arrs[1].astype(np.float32) * s1
        if narr >= 4:
            out["zr_mm"] = arrs[2].astype(np.float32) * s2
            out["pl_mm"] = arrs[3].astype(np.float32) * s3
        if narr >= 5:
            out["swe_mm"] = arrs[4].astype(np.float32) * s4
    else:
        out["probs"] = np.frombuffer(raw, np.uint8, 4 * n, off).reshape(4, nmem, ny, nx)
    return out
