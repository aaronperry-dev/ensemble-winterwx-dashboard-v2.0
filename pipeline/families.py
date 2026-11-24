"""The ensemble families the viewer can show and blend: NCEP GEFS, ECMWF IFS ENS and ECCC GEPS.

Every family is reduced to the same thing before it is published: per-member cumulative QPF / snowfall / freezing rain /
sleet and per-member precipitation-type probabilities on one 0.5 degree regional grid (see encode.py).  What differs is
how the inputs are fetched and how many pressure levels the Revised Bourgouin profile gets to see.
"""
from __future__ import annotations

from dataclasses import dataclass

from . import config as C


@dataclass(frozen=True)
class Family:
    key: str                    # folder name under data/ and the id the viewer uses
    label: str
    members: tuple[str, ...]
    levels: tuple[int, ...]     # isobaric levels (hPa) that feed the wet-bulb profile
    cycle_hours: tuple[int, ...]
    units: str                  # "member": one fetch task per member and hour;  "step": one fetch task per hour (all members)
    default_max_fhr: int
    min_members: int            # fewest members worth publishing an hour for
    credit: str
    uses_cpofp: bool = False    # only GEFS publishes the model's own percent-frozen-precipitation field

    def forecast_hours(self, max_fhr: int | None = None) -> list[int]:
        return _HOURS[self.key](self.default_max_fhr if max_fhr is None else max_fhr)

    def task_units(self, members: list[str] | None = None) -> list[str]:
        return list(members or self.members) if self.units == "member" else ["all"]


def _gefs_hours(max_fhr: int) -> list[int]:
    return C.forecast_hours(max_fhr)


def _ifs_hours(max_fhr: int) -> list[int]:
    """IFS ENS open data: every 3 h to 144, then every 6 h to 360."""
    return [h for h in list(range(3, 145, 3)) + list(range(150, 361, 6)) if h <= max_fhr]


def _geps_hours(max_fhr: int) -> list[int]:
    """GEPS: every 3 h to 192, then every 6 h to 384."""
    return [h for h in list(range(3, 193, 3)) + list(range(198, 385, 6)) if h <= max_fhr]


_HOURS = {"gefs": _gefs_hours, "ifs": _ifs_hours, "geps": _geps_hours}

IFS_LEVELS: tuple[int, ...] = (1000, 925, 850, 700, 500)
GEPS_LEVELS: tuple[int, ...] = (1000, 925, 850, 700, 500)

FAMILIES: dict[str, Family] = {
    "gefs": Family("gefs", "GEFS", C.MEMBERS, C.PLEVELS, (0, 6, 12, 18), "member", C.DEFAULT_MAX_FHR, 20,
                   "NOAA/NCEP GEFS v12 (public domain)", uses_cpofp=True),
    "ifs": Family("ifs", "ECMWF IFS ENS", tuple(f"e{i:02d}" for i in range(1, 51)), IFS_LEVELS, (0, 12), "member",
                  C.DEFAULT_MAX_FHR, 33, "ECMWF IFS ENS open data (CC-BY-4.0)"),
    "geps": Family("geps", "GEPS", ("c00",) + tuple(f"p{i:02d}" for i in range(1, 21)), GEPS_LEVELS, (0, 12), "step",
                   C.DEFAULT_MAX_FHR, 1, "Environment and Climate Change Canada GEPS (Open Government Licence - Canada)"),
}


def get(key: str) -> Family:
    try:
        return FAMILIES[key]
    except KeyError:
        raise SystemExit(f"unknown family {key!r}; choose from {', '.join(FAMILIES)}") from None
