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
