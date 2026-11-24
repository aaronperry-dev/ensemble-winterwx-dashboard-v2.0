"""Turn per-member model output into 3-hourly QPF, snowfall and cumulative totals."""
from __future__ import annotations

import numpy as np

from . import config as C


def type_shares(probs: np.ndarray) -> np.ndarray:
    """Paper Sect. 5: split QPF among types in proportion to their probabilities.

    probs: (4, ...) in percent [RA, SN, FZRA, PL] -> shares that sum to 1 (0 where all zero).
    """
    p = probs.astype(np.float32)
    tot = p.sum(axis=0)
    return np.where(tot > 0, p / np.maximum(tot, 1e-6), 0.0).astype(np.float32)


class ApcpAccumulator:
    """Rebuild cumulative-since-initialization APCP from GEFS' bucketed accumulation windows.

    GEFS APCP windows look like 0-3, 0-6, 6-9, 6-12, 12-15, 12-18 ... (the 6-h bucket resets).
    For a message covering [a, b]:   cum[b] = cum[a] + value   (cum[0] = 0).
    Cumulatives are kept per forecast hour so any window total is a simple difference.
    """

    def __init__(self, n_members: int, ny: int, nx: int):
        self.shape = (n_members, ny, nx)
        self.cum: dict[int, np.ndarray] = {0: np.zeros(self.shape, np.float32)}

    def add(self, window: tuple[int, int], values: np.ndarray) -> np.ndarray:
        a, b = window
        if a not in self.cum:
            raise KeyError(f"accumulation window {a}-{b} needs cumulative at f{a:03d}, which is missing")
        self.cum[b] = (self.cum[a] + np.nan_to_num(values, nan=0.0)).astype(np.float32)
        return self.cum[b]

    def step(self, fhr: int, prev_fhr: int) -> np.ndarray:
        """Liquid-equivalent precipitation (mm) over (prev_fhr, fhr]; clipped at 0."""
        return np.maximum(self.cum[fhr] - self.cum[prev_fhr], 0.0)


def snowfall_step_cm(qpf_mm: np.ndarray, probs: np.ndarray, slr: np.ndarray) -> np.ndarray:
    """Snowfall (cm) for one step: QPF x snow share x snow-to-liquid ratio / 10 (mm -> cm)."""
    share_sn = type_shares(probs)[1]
    return (qpf_mm * share_sn * slr / 10.0).astype(np.float32)


def type_qpf_step_mm(qpf_mm: np.ndarray, probs: np.ndarray) -> np.ndarray:
    """Liquid-equivalent precipitation (mm) of each type for one step: QPF x type share.

    probs (4, M, ny, nx) in percent [RA, SN, FZRA, PL]  ->  (4, M, ny, nx) mm.  The four parts sum to QPF
    wherever any type probability is non-zero.
    """
    return (type_shares(probs) * qpf_mm[None]).astype(np.float32)


def surface_temperature_gates(probs: np.ndarray, t2_c: np.ndarray | None, fz_max: float = C.FZRA_MAX_T2_C,
                              pl_max: float = C.PL_MAX_T2_C) -> np.ndarray:
    """2 m temperature checks on the profile-based type (the last step of the RAP/HRRR precip-type tree, applied to Bourgouin).

    probs: (4, ...) percent [RA, SN, FZRA, PL]; t2_c: 2 m temperature in deg C, same trailing shape (NaN = unknown, no gate).
    Freezing rain needs a surface near or below freezing, so where T2m is above `fz_max` its share becomes rain; ice pellets
    need cold air at the surface, so above `pl_max` theirs becomes rain too.  The wet-bulb surface test inside Bourgouin
    only catches Tw > 0 C; in dry air Tw can be below 0 C while the air is a degree or two above it.  Snow is not gated:
    Bourgouin's own melting-energy test already handles it.
    """
    if t2_c is None:
        return probs
    p = probs.astype(np.float32)
    ok = np.isfinite(t2_c)
    t = np.where(ok, t2_c, 0.0)
    out = p.copy()
    fz = ok & (t > fz_max)
    pl = ok & (t > pl_max)
    out[0] = p[0] + np.where(fz, p[2], 0.0) + np.where(pl, p[3], 0.0)
    out[2] = np.where(fz, 0.0, p[2])
    out[3] = np.where(pl, 0.0, p[3])
    return np.rint(np.minimum(out, 100.0)).astype(np.uint8)


def cap_freezing_rain(probs: np.ndarray, cpofp: np.ndarray | None, min_frozen: float = C.CPOFP_CAP_MIN_PCT) -> np.ndarray:
    """Do not let freezing rain exceed the liquid fraction that the model's own microphysics reports.

    probs: (4, ...) percent [RA, SN, FZRA, PL] (Revised Bourgouin); cpofp: percent frozen precipitation from the GEFS
    output (same trailing shape; NaN = unknown).  Where Bourgouin puts more of the precipitation into freezing rain than
    the model says is liquid (100 - CPOFP), the excess moves to the frozen types in proportion to their current share
    (all of it to snow if both are zero).  One-sided on purpose: rain is never turned into snow, so marginal-temperature
    events do not gain the snow that NCEP's own 10:1 maps are criticized for (Manikin 2021).

    The cap applies only where the model says at least `min_frozen` percent is frozen.  Checked against ASOS present weather
    (14 freezing-rain days, 2023-26), CPOFP separates freezing rain from snow only at its extreme: where Bourgouin says
    freezing rain, ZR was observed 30-33 % of the time for CPOFP below 90 but only 17 % at 90 and above.  A full cap
    (min_frozen=0) halved freezing-rain detection (POD 0.50 -> 0.25) with almost no drop in false alarms.
    """
    if cpofp is None:
        return probs
    p = probs.astype(np.float32)
    known = np.isfinite(cpofp) & (cpofp >= max(min_frozen, 0.0)) & (cpofp <= 100)
    liquid_cap = np.where(known, (100.0 - np.where(known, cpofp, 0.0)) / 100.0, 1.0) * p.sum(axis=0)
    excess = np.where(p[2] > liquid_cap, p[2] - liquid_cap, 0.0)
    frozen = p[1] + p[3]
    to_sn = np.where(frozen > 0, p[1] / np.maximum(frozen, 1e-6), 1.0)
    out = p.copy()
    out[2] = p[2] - excess
    out[1] = np.minimum(p[1] + excess * to_sn, 100.0)
    out[3] = np.minimum(p[3] + excess * (1.0 - to_sn), 100.0)
    return np.rint(out).astype(np.uint8)


def step_average_shares(cur: np.ndarray, prev: np.ndarray | None) -> np.ndarray:
    """Type shares (as percent, float) to apply to a step's QPF: the mean of the shares at the start and end of the step.

    A single end-of-step profile stands in for the whole 3 hours otherwise, so a changeover inside the step is assigned
    to one side.  Where only one end has a type (precipitation starting or ending) that end is used.
    cur / prev: (4, ...) percent, already masked to zero where the step is dry.
    """
    c = 100.0 * type_shares(cur)
    if prev is None:
        return c
    pr = 100.0 * type_shares(prev)
    has_c, has_p = (cur.sum(axis=0) > 0)[None], (prev.sum(axis=0) > 0)[None]
    return np.where(has_c & has_p, 0.5 * (c + pr), np.where(has_c, c, pr)).astype(np.float32)
