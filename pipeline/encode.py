"""Compact, fast-to-load binary files (gzip) for the browser.

File layout (little-endian), 64-byte header then raw arrays:

  offset  type   field
  0       4s     magic  b"GWX1"
  4       u8     kind   1 = cumulative totals, 2 = precipitation-type probabilities
  5       u8     reserved
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
  40..63  zero padding

 kind 1 body:  u16 [nmem][ny][nx] cumulative QPF,  then u16 [nmem][ny][nx] cumulative snowfall.
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
HEADER = struct.Struct("<4sBBHHHHHffffff")
HEADER_SIZE = 64
KIND_CUM, KIND_PTYPE = 1, 2


def _header(kind, nmem, ny, nx, fhr, step_h, domain: C.Domain, s0=0.0, s1=0.0) -> bytes:
    raw = HEADER.pack(MAGIC, kind, 0, nmem, ny, nx, fhr, step_h,
                      domain.lat_north, domain.lon_west, domain.res, domain.res, s0, s1)
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
              domain: C.Domain) -> int:
    m, ny, nx = qpf_cum_mm.shape
    body = (quantize(qpf_cum_mm, C.QPF_MM_PER_COUNT).astype("<u2").tobytes()
            + quantize(snow_cum_cm, C.SNOW_CM_PER_COUNT).astype("<u2").tobytes())
    hdr = _header(KIND_CUM, m, ny, nx, fhr, step_h, domain, C.QPF_MM_PER_COUNT, C.SNOW_CM_PER_COUNT)
    return _write_gz(Path(path), hdr + body)


def write_ptype(path: Path, fhr: int, step_h: int, probs_u8: np.ndarray, domain: C.Domain) -> int:
    """probs_u8: (4, nmem, ny, nx)."""
    _, m, ny, nx = probs_u8.shape
    hdr = _header(KIND_PTYPE, m, ny, nx, fhr, step_h, domain)
    return _write_gz(Path(path), hdr + np.ascontiguousarray(probs_u8, dtype=np.uint8).tobytes())


def read_bin(path: Path) -> dict:
    """Reference reader (used by tests; the browser has its own JS reader)."""
    raw = gzip.open(path, "rb").read()
    (magic, kind, _r, nmem, ny, nx, fhr, step_h, lat0, lon0, dlat, dlon, s0, s1) = HEADER.unpack_from(raw)
    if magic != MAGIC:
        raise ValueError("bad magic")
    n = nmem * ny * nx
    off = HEADER_SIZE
    out = {"kind": kind, "nmem": nmem, "ny": ny, "nx": nx, "fhr": fhr, "step_h": step_h,
           "lat0": lat0, "lon0": lon0, "dlat": dlat, "dlon": dlon}
    if kind == KIND_CUM:
        q = np.frombuffer(raw, "<u2", n, off).reshape(nmem, ny, nx)
        s = np.frombuffer(raw, "<u2", n, off + 2 * n).reshape(nmem, ny, nx)
        out["qpf_mm"] = q.astype(np.float32) * s0
        out["snow_cm"] = s.astype(np.float32) * s1
    else:
        out["probs"] = np.frombuffer(raw, np.uint8, 4 * n, off).reshape(4, nmem, ny, nx)
    return out
