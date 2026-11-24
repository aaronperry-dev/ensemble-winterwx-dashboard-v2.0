"""ECMWF IFS ENS open data: one member at one forecast hour.

Files:  <base>/<YYYYMMDD>/<HH>z/ifs/0p25/enfo/<YYYYMMDDHH>0000-<step>h-enfo-ef.grib2  (+ a JSON-lines .index)
Every step file holds all 50 perturbed members (the control is not part of the open data) and every field, in no
useful order, so each wanted message is one ranged request.  Precipitation (tp, metres) is accumulated from the start of
the run; 2 m humidity arrives as dew point; the 0.25 degree grid is sub-sampled onto the 0.5 degree domain.
"""
from __future__ import annotations

import json
import logging
import threading

import numpy as np

from . import config as C
from .bourgouin import es_water
from .grib import decode_messages
from .http import HttpClient, NotAvailable, SourceError
from .idx import split_grib_messages

log = logging.getLogger("gefs.ifs")

SURFACE_PARAMS = ("2t", "2d", "sp", "tp")


def _stem(date: str, cycle: str, step: int) -> str:
    return f"{date}/{cycle}z/ifs/0p25/enfo/{date}{cycle}0000-{step}h-enfo-ef"


def wanted(row: dict, levels) -> bool:
    p = row.get("param")
    if row.get("type") != "pf":
        return False
    if p in SURFACE_PARAMS:
        return True
    return p in ("t", "r") and row.get("levtype") == "pl" and row.get("levelist") in {str(int(l)) for l in levels}


class IfsSource:
    """One host that serves the IFS open-data layout (data.ecmwf.int, or its AWS copy)."""

    def __init__(self, name: str, base: str, gate, http: HttpClient, missing=(404,)):
        self.name, self.base, self.gate, self.http, self.missing = name, base.rstrip("/"), gate, http, missing
        self._idx: dict[tuple, dict[int, list]] = {}
        self._lock = threading.Lock()

    def index(self, date: str, cycle: str, step: int, levels) -> dict[int, list[tuple[int, int]]]:
        """{member number: [(offset, length), ...]} for just the fields we use (cached: the full index is 2 MB)."""
        key = (date, cycle, step, tuple(levels))
        with self._lock:
            hit = self._idx.get(key)
        if hit is not None:
            return hit
        text = self.http.get_text(f"{self.base}/{_stem(date, cycle, step)}.index", self.gate, missing_statuses=self.missing)
        out: dict[int, list[tuple[int, int]]] = {}
        for line in text.splitlines():
            if not line.strip():
                continue
            row = json.loads(line)
            if wanted(row, levels):
                out.setdefault(int(row["number"]), []).append((int(row["_offset"]), int(row["_length"])))
        if not out:
            raise SourceError(f"{self.name}: index for step {step} lists none of the wanted fields")
        with self._lock:
            while len(self._idx) >= 3:
                self._idx.pop(next(iter(self._idx)))
            self._idx[key] = out
        return out

    def exists(self, date: str, cycle: str, step: int) -> bool:
        try:
            self.http.get_range(f"{self.base}/{_stem(date, cycle, step)}.index", self.gate, 0, 31, missing_statuses=self.missing)
            return True
        except NotAvailable:
            return False
        except SourceError:
            return False

    def messages(self, date: str, cycle: str, step: int, number: int, levels) -> list[bytes]:
        idx = self.index(date, cycle, step, levels)
        if number not in idx:
            raise NotAvailable(f"{self.name}: member {number} not in the step {step} index")
        url = f"{self.base}/{_stem(date, cycle, step)}.grib2"
        msgs: list[bytes] = []
        try:
            for off, ln in sorted(idx[number]):
                blob = self.http.get_range(url, self.gate, off, off + ln - 1, missing_statuses=self.missing)
                msgs.extend(split_grib_messages(blob))
        except ValueError as e:
            raise SourceError(f"{self.name}: {e}") from e
        if len(msgs) != len(idx[number]):
            raise SourceError(f"{self.name}: expected {len(idx[number])} messages, got {len(msgs)}")
        return msgs


def build_sources(gates, http=None, bases: dict | None = None) -> list[IfsSource]:
    http = http or HttpClient()
    bases = bases or {}
    return [IfsSource("ecmwf", bases.get("ecmwf", C.ECMWF_BASE), gates["ecmwf"], http),
            IfsSource("ecmwf_aws", bases.get("ecmwf_aws", C.ECMWF_AWS_BASE), gates["ecmwf_aws"], http, missing=(404, 403))]


def to_fields(arrs: dict[str, np.ndarray], window, levels) -> dict[str, np.ndarray]:
    """Decoded IFS arrays -> the common field set (K, %, Pa, mm)."""
    need = [f"T{l}" for l in levels] + [f"RH{l}" for l in levels] + ["T2", "TD2", "PS", "APCP"]
    miss = [k for k in need if k not in arrs]
    if miss:
        raise SourceError(f"decoded IFS GRIB is missing {miss}")
    out = {k: arrs[k] for k in need if k not in ("TD2", "APCP")}
    for l in levels:
        out[f"RH{l}"] = np.clip(out[f"RH{l}"], 0.0, 100.0)
    t_c, td_c = arrs["T2"] - 273.15, np.minimum(arrs["TD2"] - 273.15, arrs["T2"] - 273.15)
    out["RH2"] = np.clip(100.0 * es_water(td_c) / es_water(t_c), 0.0, 100.0).astype(np.float32)
    if window is None:
        raise SourceError("IFS precipitation message has no accumulation window")
    out["APCP"] = (arrs["APCP"] * 1000.0).astype(np.float32)            # metres -> mm
    return out, tuple(window)


class IfsChain:
    """data.ecmwf.int first, the AWS copy if that fails.  Same fall-through rules as the GEFS SourceChain."""

    def __init__(self, sources: list[IfsSource], domain: C.Domain, levels):
        self.sources, self.domain, self.levels = sources, domain, tuple(levels)

    def exists(self, date: str, cycle: str, step: int) -> bool:
        return any(s.exists(date, cycle, step) for s in self.sources)

    def fetch_member_step(self, date: str, cycle: str, member: str, fhr: int, allow_origin: bool = True):
        number = int(member[1:])
        errors, saw_error = [], False
        for s in self.sources:
            try:
                msgs = s.messages(date, cycle, fhr, number, self.levels)
                arrs, window = decode_messages(msgs, self.domain)
                fields, window = to_fields(arrs, window, self.levels)
                return fields, window, {"ifs": s.name}
            except NotAvailable as e:
                errors.append(f"{s.name}: not available")
                log.debug("%s", e)
            except (SourceError, ValueError) as e:
                saw_error = True
                errors.append(f"{s.name}: {e}")
                log.warning("%s failed for %s f%03d: %s -- trying next source", s.name, member, fhr, e)
        if saw_error:
            raise SourceError("; ".join(errors))
        raise NotAvailable("; ".join(errors))
