"""Decode GRIB2 messages and crop them to the regional domain."""
from __future__ import annotations

import numpy as np

from . import config as C

try:  # pygrib ships wheels with a bundled ecCodes
    import pygrib
except ImportError:  # pragma: no cover
    pygrib = None


def identify(msg) -> tuple[str, tuple[int, int] | None] | None:
    """Map a GRIB message to (internal key, accumulation window) or None if unwanted.

    Keys: T<hPa>, RH<hPa>, T2, RH2, PS, APCP.  Uses WMO discipline/category/number so it works
    the same for every source.
    """
    disc, cat, num = msg["discipline"], msg["parameterCategory"], msg["parameterNumber"]
    tol, lev = msg.typeOfLevel, msg.level
    if disc != 0:
        return None
    if (cat, num) == (0, 0):          # temperature
        if tol == "isobaricInhPa":
            return f"T{int(lev)}", None
        if tol == "heightAboveGround" and lev == 2:
            return "T2", None
    elif (cat, num) == (1, 1):        # relative humidity
        if tol == "isobaricInhPa":
            return f"RH{int(lev)}", None
        if tol == "heightAboveGround" and lev == 2:
            return "RH2", None
    elif (cat, num) == (3, 0) and tol == "surface":   # surface pressure
        return "PS", None
    elif (cat, num) == (1, 8) and tol == "surface":   # total precipitation (accumulated)
        return "APCP", (int(msg["startStep"]), int(msg["endStep"]))
    return None


class Cropper:
    """Index math to cut a regular lat/lon message down to the Domain (north->south, west->east).

    Works for both the full global 0.5 deg grid (S3/GCS/NOMADS-http) and the already-cropped
    grid returned by the NOMADS filter.
    """

    def __init__(self, domain: C.Domain):
        self.domain = domain

    def crop(self, msg) -> np.ndarray:
        d = self.domain
        di = float(msg["iDirectionIncrementInDegrees"])
        dj = float(msg["jDirectionIncrementInDegrees"])
        if abs(dj - d.res) > 1e-6 or abs(di - d.res) > 1e-6:
            raise ValueError(f"unexpected grid spacing {di}x{dj}, expected {d.res}")
        nj, ni = int(msg["Nj"]), int(msg["Ni"])
        lon1 = float(msg["longitudeOfFirstGridPointInDegrees"]) % 360.0
        vals = np.ma.filled(msg.values, np.nan).astype(np.float32)
        if bool(msg["jScansPositively"]):          # rows run south->north: flip to north->south
            vals = vals[::-1]
            lat_top = float(msg["latitudeOfLastGridPointInDegrees"])
        else:
            lat_top = float(msg["latitudeOfFirstGridPointInDegrees"])
        r0 = int(round((lat_top - d.lat_north) / dj))
        if r0 < 0 or r0 + d.ny > nj:
            raise ValueError("domain outside message latitude range")
        is_global = ni * di >= 359.9
        c0 = int(round(((d.lon_west % 360.0) - lon1) / di))
        if is_global:
            cols = (np.arange(d.nx) + c0) % ni
        else:
            cols = np.arange(d.nx) + c0
            if c0 < 0 or cols.max() >= ni:
                raise ValueError("domain outside message longitude range")
        return np.ascontiguousarray(vals[r0:r0 + d.ny][:, cols])


def decode_messages(messages: list[bytes], domain: C.Domain) -> tuple[dict[str, np.ndarray], tuple[int, int] | None]:
    """Decode raw GRIB messages -> {key: float32[ny, nx]} plus the APCP accumulation window."""
    if pygrib is None:  # pragma: no cover
        raise RuntimeError("pygrib is required: pip install pygrib")
    cropper = Cropper(domain)
    out: dict[str, np.ndarray] = {}
    window = None
    for b in messages:
        msg = pygrib.fromstring(b)
        ident = identify(msg)
        if ident is None:
            continue
        key, win = ident
        if key in out:
            continue
        out[key] = cropper.crop(msg)
        if win is not None:
            window = win
    return out, window
