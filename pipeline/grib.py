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

    Keys: T<hPa>, RH<hPa>, T2, RH2, PS, APCP, and the 2 m window extremes TMAX / TMIN (with their (start, end) window).
    Uses WMO discipline/category/number so it works the same for every source.
    """
    disc, cat, num = msg["discipline"], msg["parameterCategory"], msg["parameterNumber"]
    tol, lev = msg.typeOfLevel, msg.level
    if disc != 0:
        return None
    if (cat, num) in ((0, 4), (0, 5)) and tol == "heightAboveGround" and lev == 2:   # NCEP and ECCC: maximum / minimum temperature
        return ("TMAX" if num == 4 else "TMIN"), (int(msg["startStep"]), int(msg["endStep"]))
    if (cat, num) == (0, 0):          # temperature
        if tol == "isobaricInhPa":
            return f"T{int(lev)}", None
        if tol == "heightAboveGround" and lev == 2:
            # ECMWF writes its 3- and 6-hour maximum / minimum 2 m temperature (mx2t3, mn2t6, ...) with the same parameter
            # numbers as the instantaneous 2 m temperature, and tells them apart only by the statistical processing
            step_type = msg["stepType"] if msg.has_key("stepType") else "instant"
            if step_type in ("max", "min"):
                return ("TMAX" if step_type == "max" else "TMIN"), (int(msg["startStep"]), int(msg["endStep"]))
            return "T2", None
    elif (cat, num) == (0, 6) and tol == "heightAboveGround" and lev == 2:    # 2 m dew point (ECMWF sends this, not RH)
        return "TD2", None
    elif (cat, num) == (1, 1):        # relative humidity
        if tol == "isobaricInhPa":
            return f"RH{int(lev)}", None
        if tol == "heightAboveGround" and lev == 2:
            return "RH2", None
    elif (cat, num) == (3, 0) and tol == "surface":   # surface pressure
        return "PS", None
    elif (cat, num) in ((1, 8), (1, 193)) and tol == "surface":   # total precipitation (accumulated); ECMWF "tp" is (1, 193), in metres
        return "APCP", (int(msg["startStep"]), int(msg["endStep"]))
    elif (cat, num) == (1, 39) and tol == "surface":  # percent frozen precipitation (-50 = no precipitation)
        return "CPOFP", None
    return None


class Cropper:
    """Index math to cut a regular lat/lon message down to the Domain (north->south, west->east).

    Works for the global 0.5 deg grid (S3/GCS/NOMADS-http, GEPS), for the already-cropped grid returned by the NOMADS
    filter, and for finer grids whose spacing is a whole number of times finer than the domain's (the 0.25 deg IFS grid:
    every second point is a 0.5 deg point, so it is sub-sampled exactly, with no interpolation).
    """

    def __init__(self, domain: C.Domain):
        self.domain = domain

    def crop(self, msg) -> np.ndarray:
        d = self.domain
        di = float(msg["iDirectionIncrementInDegrees"])
        dj = float(msg["jDirectionIncrementInDegrees"])
        ki, kj = d.res / di, d.res / dj
        if abs(ki - round(ki)) > 1e-6 or abs(kj - round(kj)) > 1e-6 or round(ki) < 1 or round(kj) < 1:
            raise ValueError(f"unexpected grid spacing {di}x{dj}, expected {d.res} or a whole fraction of it")
        ki, kj = int(round(ki)), int(round(kj))
        nj, ni = int(msg["Nj"]), int(msg["Ni"])
        lon1 = float(msg["longitudeOfFirstGridPointInDegrees"]) % 360.0
        vals = np.ma.filled(msg.values, np.nan).astype(np.float32)
        if bool(msg["jScansPositively"]):          # rows run south->north: flip to north->south
            vals = vals[::-1]
            lat_top = float(msg["latitudeOfLastGridPointInDegrees"])
        else:
            lat_top = float(msg["latitudeOfFirstGridPointInDegrees"])
        r0 = int(round((lat_top - d.lat_north) / dj))
        rows = r0 + kj * np.arange(d.ny)
        if r0 < 0 or rows.max() >= nj:
            raise ValueError("domain outside message latitude range")
        is_global = ni * di >= 359.9
        c0 = int(round(((d.lon_west % 360.0) - lon1) / di))
        cols = c0 + ki * np.arange(d.nx)
        if is_global:
            cols = cols % ni
        elif c0 < 0 or cols.max() >= ni:
            raise ValueError("domain outside message longitude range")
        return np.ascontiguousarray(vals[rows][:, cols])


def decode_messages(messages: list[bytes], domain: C.Domain, windows: dict | None = None) -> tuple[dict[str, np.ndarray], tuple[int, int] | None]:
    """Decode raw GRIB messages -> {key: float32[ny, nx]} plus the APCP accumulation window.
    `windows`, if given, is filled with {key: (start, end)} for every field that has one (APCP, TMAX, TMIN)."""
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
            if windows is not None:
                windows[key] = win
            if key == "APCP":
                window = win
    return out, window


def decode_members(messages: list[bytes], domain: C.Domain, windows: dict | None = None) -> tuple[dict[str, dict[int, np.ndarray]], tuple[int, int] | None]:
    """Decode a file that holds every ensemble member (the GEPS layout) -> {key: {perturbation number: float32[ny, nx]}}.
    `windows` is filled as in decode_messages."""
    if pygrib is None:  # pragma: no cover
        raise RuntimeError("pygrib is required: pip install pygrib")
    cropper = Cropper(domain)
    out: dict[str, dict[int, np.ndarray]] = {}
    window = None
    for b in messages:
        msg = pygrib.fromstring(b)
        ident = identify(msg)
        if ident is None:
            continue
        key, win = ident
        n = int(msg["perturbationNumber"])
        if n in out.setdefault(key, {}):
            continue
        out[key][n] = cropper.crop(msg)
        if win is not None:
            if windows is not None:
                windows[key] = win
            if key == "APCP":
                window = win
    return out, window
