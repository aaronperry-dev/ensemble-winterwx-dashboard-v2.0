"""Central configuration: members, levels, domain, per-host rate policies."""
from __future__ import annotations

import os
from dataclasses import dataclass

# --------------------------------------------------------------------------- members
# GEFS has 1 control + 30 perturbed members = 31 total.
MEMBERS: tuple[str, ...] = ("c00",) + tuple(f"p{i:02d}" for i in range(1, 31))


def member_file_prefix(member: str) -> str:
    """c00 -> gec00, p07 -> gep07 (GEFS file naming)."""
    return "gec00" if member == "c00" else f"gep{member[1:]}"


# --------------------------------------------------------------------------- vertical levels
# TMP and RH on these isobaric levels feed the wet-bulb profile. pgrb2a (0.5 deg) only carries
# the "mandatory-ish" levels, pgrb2b carries the rest. Verified against the live S3 inventory.
LEVELS_A: tuple[int, ...] = (1000, 925, 850, 700, 500)
LEVELS_B: tuple[int, ...] = (975, 950, 900, 800, 750, 650, 600, 550, 450, 400)
PLEVELS: tuple[int, ...] = tuple(sorted(LEVELS_A + LEVELS_B, reverse=True))

# --------------------------------------------------------------------------- forecast hours
# GEFS atmos output is every 3 h to f240, then every 6 h to f384.
DEFAULT_MAX_FHR = 240


def forecast_hours(max_fhr: int = DEFAULT_MAX_FHR) -> list[int]:
    hrs = list(range(3, min(max_fhr, 240) + 1, 3))
    hrs += list(range(246, max_fhr + 1, 6))
    return hrs


# --------------------------------------------------------------------------- domain
@dataclass(frozen=True)
class Domain:
    """Regional crop of the global 0.5 deg grid (north-to-south rows, west-to-east columns)."""

    lat_north: float = 60.0
    lat_south: float = 20.0
    lon_west: float = -135.0
    lon_east: float = -50.0
    res: float = 0.5

    @property
    def ny(self) -> int:
        return int(round((self.lat_north - self.lat_south) / self.res)) + 1

    @property
    def nx(self) -> int:
        return int(round((self.lon_east - self.lon_west) / self.res)) + 1

    def to_dict(self) -> dict:
        return {
            "lat_north": self.lat_north,
            "lat_south": self.lat_south,
            "lon_west": self.lon_west,
            "lon_east": self.lon_east,
            "dlat": self.res,
            "dlon": self.res,
            "ny": self.ny,
            "nx": self.nx,
        }


# --------------------------------------------------------------------------- precipitation / snow
# 3-h liquid-equivalent totals below this (mm) are treated as "no precipitation" for p-type.
MIN_PRECIP_MM = 0.05

# Snow-to-liquid ratio mode: "kuchera" (default) or "fixed:<ratio>" e.g. "fixed:12.5".
SLR_MODE = os.environ.get("SLR_MODE", "kuchera")
SLR_MIN, SLR_MAX = 5.0, 30.0

# --------------------------------------------------------------------------- binary quantization
QPF_MM_PER_COUNT = 0.05   # uint16 cumulative liquid-equivalent  (max 3276 mm)
SNOW_CM_PER_COUNT = 0.1   # uint16 cumulative snowfall            (max 6553 cm)
ZR_MM_PER_COUNT = 0.01    # uint16 cumulative freezing-rain liquid equivalent (max 655 mm)
PL_MM_PER_COUNT = 0.01    # uint16 cumulative sleet / ice-pellet liquid equivalent (max 655 mm)
SWE_MM_PER_COUNT = 0.02   # uint16 cumulative snow liquid equivalent, any snow ratio can be applied in the viewer (max 1310 mm)
NCEP_NONE, NCEP_UNKNOWN = 254, 255   # codes in the ptype file's NCEP-type plane (0-3 = RA, SN, FZRA, PL)

# --------------------------------------------------------------------------- data sources
S3_BASE = "https://noaa-gefs-pds.s3.amazonaws.com"
GCS_BASE = "https://storage.googleapis.com/gfs-ensemble-forecast-system"
NOMADS_HTTP_BASE = "https://nomads.ncep.noaa.gov/pub/data/nccf/com/gens/prod"
NOMADS_FILTER_BASE = "https://nomads.ncep.noaa.gov/cgi-bin"

# Order requested: S3 -> Google Cloud -> NOMADS GRIB filter -> NOMADS HTTP.
DEFAULT_SOURCE_ORDER: tuple[str, ...] = ("aws", "gcs", "nomads_filter", "nomads_http")

# Rate policies are *totals for the whole job*; each worker process gets rate / n_workers.
#  - aws/gcs: public object stores, no published per-IP cap; we still self-limit politely.
#  - nomads: NCEP's published policy is 120 hits/minute per IP summed across all NCEP sites
#    (NWS SCN21-32). We stay at 90/min so we never trip it, and the filter + http endpoints
#    share ONE gate because the limit is per IP, not per endpoint.
HOST_POLICIES: dict[str, dict] = {
    "aws": {"rate": 80.0, "burst": 24, "concurrency": 32},
    "gcs": {"rate": 60.0, "burst": 20, "concurrency": 32},
    "nomads": {"rate": 1.5, "burst": 3, "concurrency": 2, "max_404_per_min": 8},
}

USER_AGENT = "gefs-winter-wx/1.0 (+https://github.com/{repo})".format(
    repo=os.environ.get("GITHUB_REPOSITORY", "local")
)
