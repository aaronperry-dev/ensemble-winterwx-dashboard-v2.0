"""ECCC GEPS (Canadian ensemble) from the MSC datamart: all 21 members at one forecast hour.

Files:  <base>/<YYYYMMDD>/WXO-DD/ensemble/geps/grib2/raw/<HH>/<hhh>/CMC_geps-raw_<VAR>_<LEVEL>_latlon0p5x0p5_<YYYYMMDDHH>_P<hhh>_allmbrs.grib2
One file per variable and level, each holding the 21 members (perturbation number 0 = control) on the global 0.5 degree
grid, so there is no byte-range planning: a forecast hour is a handful of whole-file downloads.  APCP_SFC_0 is
accumulated from the start of the run, in mm.
"""
from __future__ import annotations

import logging

import numpy as np

from . import config as C
from .grib import decode_members
from .http import HttpClient, NotAvailable, SourceError
from .idx import split_grib_messages

log = logging.getLogger("gefs.geps")


def file_specs(levels) -> list[tuple[str, str]]:
    """[(variable, level code)] for one forecast hour."""
    specs = [("TMP", f"ISBL_{int(l):04d}") for l in levels] + [("RH", f"ISBL_{int(l):04d}") for l in levels]
    return specs + [("TMP", "TGL_2m"), ("RH", "TGL_2m"), ("PRES", "SFC_0"), ("APCP", "SFC_0")]


class GepsChain:
    def __init__(self, gates, domain: C.Domain, levels, http: HttpClient | None = None, base: str | None = None):
        self.http = http or HttpClient()
        self.gate = gates["msc"]
        self.base = (base or C.MSC_BASE).rstrip("/")
        self.domain, self.levels = domain, tuple(levels)

    def url(self, date: str, cycle: str, fhr: int, var: str, level: str) -> str:
        return (f"{self.base}/{date}/WXO-DD/ensemble/geps/grib2/raw/{cycle}/{fhr:03d}/"
                f"CMC_geps-raw_{var}_{level}_latlon0p5x0p5_{date}{cycle}_P{fhr:03d}_allmbrs.grib2")

    def exists(self, date: str, cycle: str, fhr: int) -> bool:
        try:
            self.http.get_range(self.url(date, cycle, fhr, "APCP", "SFC_0"), self.gate, 0, 31, missing_statuses=(404,))
            return True
        except (NotAvailable, SourceError):
            return False

    def fetch_step(self, date: str, cycle: str, fhr: int, allow_origin: bool = True):
        """-> ({member index: fields}, (a, b) accumulation window, {"geps": "msc"})."""
        merged: dict[str, dict[int, np.ndarray]] = {}
        window = None
        for var, level in file_specs(self.levels):
            blob = self.http.get_bytes(self.url(date, cycle, fhr, var, level), self.gate, missing_statuses=(404,))
            try:
                arrs, win = decode_members(split_grib_messages(blob), self.domain)
            except ValueError as e:
                raise SourceError(f"msc: {e}") from e
            for k, v in arrs.items():
                merged.setdefault(k, {}).update(v)
            window = win or window
        if window is None:
            raise SourceError(f"GEPS f{fhr:03d}: no precipitation accumulation window")
        need = [f"T{l}" for l in self.levels] + [f"RH{l}" for l in self.levels] + ["T2", "RH2", "PS", "APCP"]
        miss = [k for k in need if k not in merged]
        if miss:
            raise SourceError(f"GEPS f{fhr:03d}: decoded GRIB is missing {miss}")
        members = sorted(set.intersection(*(set(merged[k]) for k in need)))
        if not members:
            raise SourceError(f"GEPS f{fhr:03d}: no member has every field")
        out = {}
        for n in members:
            f = {k: merged[k][n] for k in need}
            for l in self.levels:
                f[f"RH{l}"] = np.clip(f[f"RH{l}"], 0.0, 100.0)
            f["RH2"] = np.clip(f["RH2"], 0.0, 100.0)
            out[n] = f
        return out, tuple(window), {"geps": "msc"}
