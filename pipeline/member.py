"""One GEFS member at one forecast hour: fetch -> decode -> Bourgouin p-type + snow ratio."""
from __future__ import annotations

import numpy as np

from . import config as C
from . import tempext
from .bourgouin import classify_profiles
from .grib import decode_messages
from .http import SourceError
from .sources import SourceChain, Target, required_fields


def fetch_fields(chain: SourceChain, date: str, cycle: str, member: str, fhr: int,
                 levels=C.PLEVELS, domain: C.Domain | None = None, allow_origin: bool = True, extremes: bool = False):
    """Download + decode all required fields. Returns (fields, apcp_window, {a|b: source}).
    With `extremes`, `fields` may also hold TMAX / TMIN (2 m window extremes) and TWIN, their (start, end) hours."""
    domain = domain or C.Domain()
    specs = required_fields(levels, extremes)
    fields: dict[str, np.ndarray] = {}
    wins: dict = {}
    window, used = None, {}
    for kind in ("a", "b"):
        sp = [s for s in specs if s.file == kind]
        if not sp:
            continue
        msgs, src = chain.fetch(Target(date, cycle, member, fhr, kind), sp, allow_origin=allow_origin)
        arrs, win = decode_messages(msgs, domain, wins)
        fields.update(arrs)
        used[kind] = src
        window = win or window
    tempext.attach(fields, fields, wins)
    for k in ("TMAX", "TMIN"):
        if k in fields and tempext.TWIN not in fields:
            del fields[k]                                  # a lone or mismatched extreme is of no use
    missing = [s.key for s in specs if s.key not in fields and not s.optional]
    if missing:
        raise SourceError(f"{member} f{fhr:03d}: decoded GRIB is missing {missing}")
    if window is None:
        raise SourceError(f"{member} f{fhr:03d}: no APCP accumulation window")
    return fields, window, used


def compute_member_step(fields: dict[str, np.ndarray], levels=C.PLEVELS, slr_mode: str = C.SLR_MODE,
                        slr_min: float = C.SLR_MIN, slr_max: float = C.SLR_MAX):
    """Pure compute step (no I/O) -> (probs uint8 (4,ny,nx), slr float32 (ny,nx))."""
    t_lev = np.stack([fields[f"T{l}"] for l in levels])
    rh_lev = np.stack([fields[f"RH{l}"] for l in levels])
    probs, slr, _ = classify_profiles(
        np.asarray(levels, float),
        np.nan_to_num(t_lev, nan=273.15), np.nan_to_num(rh_lev, nan=50.0),
        np.nan_to_num(fields["T2"], nan=273.15), np.nan_to_num(fields["RH2"], nan=50.0),
        np.nan_to_num(fields["PS"], nan=101325.0),
        slr_mode=slr_mode, slr_min=slr_min, slr_max=slr_max,
    )
    return probs, slr


def cpofp_from(fields: dict[str, np.ndarray]):
    """Model-side percent frozen precipitation (CPOFP) -> float32 (ny, nx), NaN where the model has no precipitation
    (it uses -50) or the field was not in the files (None then)."""
    cp = fields.get("CPOFP")
    return None if cp is None else np.where((cp >= 0) & (cp <= 100), cp, np.nan).astype(np.float32)
