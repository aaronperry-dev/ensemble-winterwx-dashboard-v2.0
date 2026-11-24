"""One GEFS member at one forecast hour: fetch -> decode -> Bourgouin p-type + snow ratio."""
from __future__ import annotations

import numpy as np

from . import config as C
from .bourgouin import classify_profiles
from .grib import decode_messages
from .http import SourceError
from .sources import SourceChain, Target, required_fields


def fetch_fields(chain: SourceChain, date: str, cycle: str, member: str, fhr: int,
                 levels=C.PLEVELS, domain: C.Domain | None = None, allow_origin: bool = True):
    """Download + decode all required fields. Returns (fields, apcp_window, {a|b: source})."""
    domain = domain or C.Domain()
    specs = required_fields(levels)
    fields: dict[str, np.ndarray] = {}
    window, used = None, {}
    for kind in ("a", "b"):
        sp = [s for s in specs if s.file == kind]
        if not sp:
            continue
        msgs, src = chain.fetch(Target(date, cycle, member, fhr, kind), sp, allow_origin=allow_origin)
        arrs, win = decode_messages(msgs, domain)
        fields.update(arrs)
        used[kind] = src
        window = win or window
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


NCEP_PRIORITY = (0, 3, 1, 2)      # applied in this order, so the last one wins: freezing rain > snow > sleet > rain


def type_extras(fields: dict[str, np.ndarray]):
    """Optional model-side precip-type information -> (cpofp, ncep_code); either is None when it was not in the files.

    cpofp:      float32 (ny, nx) percent frozen precipitation; NaN where the model has no precipitation (it uses -50).
    ncep_code:  uint8 (ny, nx): the type from NCEP's own multi-algorithm dominant-type flags (0 RA, 1 SN, 2 FZRA, 3 PL),
                C.NCEP_NONE where no flag is set (no precipitation at the valid time).
    """
    cp = fields.get("CPOFP")
    cpofp = None if cp is None else np.where((cp >= 0) & (cp <= 100), cp, np.nan).astype(np.float32)
    flags = [fields.get(k) for k in ("CRAIN", "CSNOW", "CFRZR", "CICEP")]       # order = RA, SN, FZRA, PL
    if any(f is None for f in flags):
        return cpofp, None
    on = np.stack([np.nan_to_num(f, nan=0.0) > 0.5 for f in flags])
    code = np.full(on.shape[1:], C.NCEP_NONE, np.uint8)
    for t in NCEP_PRIORITY:
        code[on[t]] = t
    return cpofp, code
