"""Daily maximum and minimum 2 m temperature, per member, from the models' own window extremes.

Each model publishes the highest and lowest 2 m temperature over a window ending at the forecast hour, so nothing is
estimated from the 3-hourly instantaneous temperature (which misses afternoon peaks and pre-dawn lows):

    GEFS      TMAX / TMIN   6-hour windows (0-6, 6-12, ...) in the 0.5 degree pgrb2a files; the 3-hour ones are not needed
    IFS ENS   mx2t3 / mn2t3 3-hour windows to F144, then mx2t6 / mn2t6 6-hour windows
    GEPS      TMAX / TMIN   6-hour windows, published only at the 6-hourly hours

Everything is reduced to 6-hour **blocks** (block a covers hours a..a+6 after the run time; IFS builds one from two 3-hour
windows).  The daily values follow the NWS convention for Eastern time, in fixed UTC windows:

    MaxT of a date   the highest temperature from 12Z to 00Z   (the 12-18Z and 18-24Z blocks;  7 am to 7 pm EST)
    MinT of a date   the lowest temperature from 00Z to 12Z    (the 00-06Z and 06-12Z blocks;  7 pm to 7 am EST, the overnight low)

The windows are the same in every column of the domain.  Only complete windows are shown: a window that a run cannot
cover with both of its blocks (the first or last half-window) is kept in the manifest (so a resumed run can finish it) but its
`expected` is below `blocks_per_window`, and the viewer shows a window only once both blocks are in.  A member that lacks the extremes for one block (field missing, or dropped from the run) is removed from that
window's member mask, because a maximum over half the window would be wrong.
"""
from __future__ import annotations

import datetime as dt
import logging
import warnings
from collections import Counter
from pathlib import Path

import numpy as np

from . import config as C
from . import encode

log = logging.getLogger("gefs.temps")
TWIN = "TWIN"          # key under which a member's field dict carries the (start, end) hours of its TMAX / TMIN window
PRODUCTS = ("tmax", "tmin")
BLOCKS_PER_WINDOW = {p: (w[1] - w[0]) // C.TEMP_BLOCK_H for p, w in C.TEMP_WINDOWS_UTC.items()}


def wanted(family: str, fhr: int) -> bool:
    """Does this family publish window extremes at this forecast hour?  (Asking for a file that does not exist would stall the hour.)"""
    if family == "ifs":
        return True                       # 3-hour windows to F144 and 6-hour windows after: every hour has them
    return fhr % C.TEMP_BLOCK_H == 0      # GEFS and GEPS: the 6-hour windows; GEPS has no file at the 3-hourly hours


def attach(fields: dict, arrs: dict, wins: dict) -> None:
    """Put decoded TMAX / TMIN (and their window) into a member's field dict when both are there and agree."""
    tx, tn = arrs.get("TMAX"), arrs.get("TMIN")
    wx, wn = wins.get("TMAX"), wins.get("TMIN")
    if tx is None or tn is None or wx is None or wx != wn:
        return
    fields["TMAX"], fields["TMIN"] = tx, tn
    fields[TWIN] = np.asarray(wx, np.int32)


def locate(cyc: dt.datetime, block_start_h: int) -> tuple[str, str]:
    """(product, date YYYY-MM-DD) of the window that the block starting `block_start_h` hours after `cyc` belongs to."""
    t0 = cyc + dt.timedelta(hours=block_start_h)
    for prod, (lo, hi) in C.TEMP_WINDOWS_UTC.items():
        if lo <= t0.hour < hi:
            return prod, t0.strftime("%Y-%m-%d")
    raise ValueError(f"block start {t0} is in no temperature window")


def plan(cyc: dt.datetime, fhrs: list[int]) -> dict[str, dict[str, int]]:
    """{product: {date: number of 6-hour blocks the run supplies}}, in date order.  A window with fewer blocks than
    BLOCKS_PER_WINDOW is only half covered by this run and is never shown."""
    n: dict[str, dict[str, int]] = {p: {} for p in PRODUCTS}
    for b in sorted(set(fhrs)):
        if b % C.TEMP_BLOCK_H or b < C.TEMP_BLOCK_H:
            continue
        prod, lab = locate(cyc, b - C.TEMP_BLOCK_H)
        n[prod][lab] = n[prod].get(lab, 0) + 1
    return {p: dict(sorted(d.items())) for p, d in n.items()}


def _fill(a: np.ndarray, ok: np.ndarray) -> np.ndarray:
    """Replace NaN cells of the valid members by the mean of the other valid members at that cell (a stray missing value must not
    become 'no data' for the whole map); a cell with no valid member at all gets 273.15 K."""
    out = a.copy()
    out[~ok] = np.nan
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)              # "mean of empty slice" where no member has a value
        mean = np.nanmean(out[ok], axis=0) if ok.any() else np.full(a.shape[1:], np.nan, np.float32)
    mean = np.where(np.isfinite(mean), mean, np.float32(273.15))
    for i in np.flatnonzero(ok):
        out[i] = np.where(np.isfinite(out[i]), out[i], mean)
    out[~ok] = 0.0
    return out.astype(np.float32)


class DailyExtremes:
    """Running per-member maximum (MaxT) and minimum (MinT) 2 m temperature, in Kelvin, for every window a run covers."""

    def __init__(self, cyc: dt.datetime, fhrs: list[int], n_members: int, ny: int, nx: int):
        self.cyc, self.M, self.ny, self.nx = cyc, n_members, ny, nx
        self.items: dict[tuple[str, str], dict] = {}
        for prod, dates in plan(cyc, fhrs).items():
            for k, (lab, exp) in enumerate(dates.items()):
                self.items[(prod, lab)] = {"product": prod, "index": k, "date": lab, "expected": exp, "blocks": 0, "upto": 0,
                                           "mask": encode.full_mask(n_members), "t": None}
        self._half: dict[int, dict] = {}        # IFS: the first 3-hour window of a block, until the second one arrives

    # ------------------------------------------------------------------ feeding
    def file_name(self, key: tuple[str, str]) -> str:
        it = self.items[key]
        return f"tmp/{it['product']}{it['index']:02d}.bin.gz"

    def step(self, fhr: int, got: dict[int, tuple]) -> list[tuple[str, str]]:
        """Take one forecast hour.  `got` maps member index -> (tmax_K, tmin_K, (start, end)) for the members that supplied
        extremes at this hour.  Returns the (product, date) keys whose arrays changed (a block was finished)."""
        wins = Counter(w for _, _, w in got.values())
        win = wins.most_common(1)[0][0] if wins else None
        sel = {i: (tx, tn) for i, (tx, tn, w) in got.items() if w == win}
        block = None
        half = C.TEMP_BLOCK_H // 2
        if win is not None and win[1] == fhr:
            length = win[1] - win[0]
            if length == C.TEMP_BLOCK_H and fhr % C.TEMP_BLOCK_H == 0:
                block = sel
            elif length == half:
                if win[0] % C.TEMP_BLOCK_H == 0:
                    self._half[win[0]] = sel                                  # first half: wait for the second
                elif fhr % C.TEMP_BLOCK_H == 0:
                    first = self._half.pop(win[0] - half, {})
                    block = {i: (np.fmax(first[i][0], tx), np.fmin(first[i][1], tn)) for i, (tx, tn) in sel.items() if i in first}
        if fhr % C.TEMP_BLOCK_H:
            return []
        self._half = {a: v for a, v in self._half.items() if a >= fhr - C.TEMP_BLOCK_H}
        return self._apply(fhr - C.TEMP_BLOCK_H, block or {})

    def _apply(self, a: int, block: dict[int, tuple]) -> list[tuple[str, str]]:
        key = locate(self.cyc, a)
        it = self.items.get(key)
        if it is None:
            return []
        is_max = it["product"] == "tmax"
        if it["t"] is None:
            it["t"] = np.full((self.M, self.ny, self.nx), np.nan, np.float32)
        for i in range(self.M):
            if i in block:
                it["t"][i] = (np.fmax if is_max else np.fmin)(it["t"][i], block[i][0 if is_max else 1])
            else:
                it["mask"] &= ~(1 << i)                                       # no extremes for this member in this block
        it["blocks"] += 1
        it["upto"] = a + C.TEMP_BLOCK_H
        return [key]

    # ------------------------------------------------------------------ output
    def n_members(self, key: tuple[str, str]) -> int:
        return bin(self.items[key]["mask"]).count("1") if self.items[key]["t"] is not None else 0

    def write_day(self, run_dir: Path, key: tuple[str, str], domain: C.Domain) -> int:
        """Write (or rewrite) one window's file.  A window with no usable member gets no file."""
        it = self.items[key]
        if it["t"] is None or not it["mask"]:
            return 0
        ok = np.array([(it["mask"] >> i) & 1 for i in range(self.M)], bool)
        return encode.write_temps(run_dir / self.file_name(key), it["index"], it["blocks"], _fill(it["t"], ok), domain, it["mask"])

    def manifest(self) -> dict:
        out: dict = {"windows_utc": {p: list(w) for p, w in C.TEMP_WINDOWS_UTC.items()}, "block_h": C.TEMP_BLOCK_H,
                     "blocks_per_window": BLOCKS_PER_WINDOW, "c_per_count": C.TEMP_C_PER_COUNT, "offset_c": C.TEMP_C_OFFSET,
                     "tmax": [], "tmin": []}
        for key, it in self.items.items():
            if not it["blocks"]:
                continue
            n = self.n_members(key)
            out[it["product"]].append({"date": it["date"], "index": it["index"], "file": self.file_name(key) if n else None,
                                       "expected": it["expected"], "blocks": it["blocks"], "n": n, "mask": it["mask"] if n else 0,
                                       "upto": it["upto"]})
        return out

    def available(self) -> bool:
        """Is there at least one finished window with members?  (The viewer offers the products only then.)"""
        return any(it["blocks"] >= BLOCKS_PER_WINDOW[it["product"]] and it["t"] is not None and it["mask"] for it in self.items.values())

    # ------------------------------------------------------------------ resuming
    def restore(self, listing: dict, run_dir: Path, last_fhr: int) -> None:
        """Reload the running state from the files that an earlier run of the same cycle published.  Blocks that ended by
        `last_fhr` count as done; maxima and minima are idempotent, so a file that already holds a later block does no harm."""
        for prod in PRODUCTS:
            for ent in listing.get(prod, []):
                it = self.items.get((prod, ent["date"]))
                if it is None or not ent.get("file"):
                    continue
                x = encode.read_bin(run_dir / ent["file"])
                it["t"] = np.where(np.isfinite(x["t_c"]), x["t_c"] + 273.15, np.nan).astype(np.float32)
                it["mask"] = int(x["mask"])
        for key, it in self.items.items():
            done = [b for b in range(C.TEMP_BLOCK_H, last_fhr + 1, C.TEMP_BLOCK_H) if locate(self.cyc, b - C.TEMP_BLOCK_H) == key]
            it["blocks"] = len(done)
            it["upto"] = done[-1] if done else 0
            if done and it["t"] is None:                                      # blocks were done but no member had data: stays empty
                it["mask"] = 0
