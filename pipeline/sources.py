"""GEFS data sources with ordered fallback.

Order (as requested): AWS S3 -> Google Cloud -> NOMADS GRIB filter -> NOMADS HTTP.

* aws / gcs / nomads_http are *byte-range* sources: fetch the small .idx inventory, then request
  only the GRIB messages we need, merging byte-adjacent messages into one request.
* nomads_filter is a server-side subsetting CGI (variable + level + lat/lon box). Its output is
  generated on the fly, so there is no file to byte-range; the server-side subset is the
  equivalent bandwidth saver.

Every source returns raw GRIB2 messages. Fields are identified from the GRIB headers later
(grib.py), not from the source, so all four sources are interchangeable.
"""
from __future__ import annotations

import logging
import time
from dataclasses import dataclass

from . import config as C
from .http import HttpClient, NotAvailable, SourceError
from .idx import find_entry, merge_ranges, parse_idx, split_grib_messages
from .ratelimit import HostGate

log = logging.getLogger(__name__)


# --------------------------------------------------------------------------- data model
@dataclass(frozen=True)
class Target:
    """One GRIB file: (cycle, member, forecast hour, a|b product)."""

    date: str      # YYYYMMDD
    cycle: str     # "00" | "06" | "12" | "18"
    member: str    # c00, p01 ... p30
    fhr: int
    kind: str      # "a" (pgrb2a) or "b" (pgrb2b)

    @property
    def file(self) -> str:
        return (f"{C.member_file_prefix(self.member)}.t{self.cycle}z."
                f"pgrb2{self.kind}.0p50.f{self.fhr:03d}")

    @property
    def subdir(self) -> str:
        return f"gefs.{self.date}/{self.cycle}/atmos/pgrb2{self.kind}p5"


@dataclass(frozen=True)
class FieldSpec:
    file: str   # "a" | "b"
    var: str    # GRIB short name as it appears in the .idx ("TMP", "RH", "PRES", "APCP")
    level: str  # as it appears in the .idx ("850 mb", "2 m above ground", "surface")
    key: str    # our internal name: T850, RH850, T2, RH2, PS, APCP


def required_fields(levels=C.PLEVELS) -> list[FieldSpec]:
    specs = [
        FieldSpec("a", "TMP", "2 m above ground", "T2"),
        FieldSpec("a", "RH", "2 m above ground", "RH2"),
        FieldSpec("a", "PRES", "surface", "PS"),
        FieldSpec("a", "APCP", "surface", "APCP"),
    ]
    for lev in levels:
        f = "a" if lev in C.LEVELS_A else "b"
        specs.append(FieldSpec(f, "TMP", f"{lev} mb", f"T{lev}"))
        specs.append(FieldSpec(f, "RH", f"{lev} mb", f"RH{lev}"))
    return specs


# --------------------------------------------------------------------------- sources
class Source:
    name = "base"
    gate_key = "base"

    def __init__(self, http: HttpClient, gates: dict[str, HostGate], base: str | None = None):
        self.http = http
        self.gate: HostGate = gates[self.gate_key]
        self.base = base
        self._fail_streak = 0
        self._open_until = 0.0

    # circuit breaker ---------------------------------------------------------
    def available(self) -> bool:
        return time.monotonic() >= self._open_until

    def _ok(self):
        self._fail_streak = 0

    def _fail(self):
        self._fail_streak += 1
        if self._fail_streak >= 5:
            self._open_until = time.monotonic() + 120.0
            self._fail_streak = 0
            log.warning("%s: circuit open for 120 s after repeated failures", self.name)

    # API ---------------------------------------------------------------------
    def fetch(self, target: Target, specs: list[FieldSpec]) -> list[bytes]:
        raise NotImplementedError

    def exists(self, target: Target) -> bool:
        raise NotImplementedError


class IdxRangeSource(Source):
    """Object-store style source: <url> + <url>.idx, with HTTP Range support."""

    def url(self, t: Target) -> str:
        return f"{self.base}/{t.subdir}/{t.file}"

    @property
    def missing_statuses(self) -> tuple[int, ...]:
        return (404,)

    def exists(self, t: Target) -> bool:
        try:
            self.http.get_range(self.url(t) + ".idx", self.gate, 0, 31, missing_statuses=self.missing_statuses)
            self._ok()
            return True
        except NotAvailable:
            return False
        except SourceError:
            self._fail()
            return False

    def fetch(self, t: Target, specs: list[FieldSpec]) -> list[bytes]:
        url = self.url(t)
        try:
            entries = parse_idx(self.http.get_text(url + ".idx", self.gate,
                                                   missing_statuses=self.missing_statuses))
            chosen = []
            for s in specs:
                e = find_entry(entries, s.var, s.level)
                if e is None:
                    raise SourceError(f"{self.name}: {s.var}:{s.level} not in inventory of {t.file}")
                chosen.append(e)
            blobs = [self.http.get_range(url, self.gate, a, b, missing_statuses=self.missing_statuses)
                     for a, b in merge_ranges(chosen)]
            msgs = [m for blob in blobs for m in split_grib_messages(blob)]
            if len(msgs) != len(chosen):
                raise SourceError(f"{self.name}: expected {len(chosen)} messages, got {len(msgs)} for {t.file}")
        except NotAvailable:
            raise
        except SourceError:
            self._fail()
            raise
        except ValueError as e:  # corrupt GRIB framing
            self._fail()
            raise SourceError(f"{self.name}: {e}") from e
        self._ok()
        return msgs


class AwsSource(IdxRangeSource):
    name, gate_key = "aws", "aws"

    def __init__(self, http, gates, base=None):
        super().__init__(http, gates, base or C.S3_BASE)

    @property
    def missing_statuses(self):
        return (404, 403)  # S3 answers 403 for absent keys when listing is not granted


class GcsSource(IdxRangeSource):
    name, gate_key = "gcs", "gcs"

    def __init__(self, http, gates, base=None):
        super().__init__(http, gates, base or C.GCS_BASE)

    @property
    def missing_statuses(self):
        return (404, 403)


class NomadsHttpSource(IdxRangeSource):
    name, gate_key = "nomads_http", "nomads"

    def __init__(self, http, gates, base=None):
        super().__init__(http, gates, base or C.NOMADS_HTTP_BASE)


class NomadsFilterSource(Source):
    """NOMADS grib_filter CGI: server-side variable/level/sub-region subsetting."""

    name, gate_key = "nomads_filter", "nomads"

    def __init__(self, http, gates, base=None, domain: C.Domain | None = None):
        super().__init__(http, gates, base or C.NOMADS_FILTER_BASE)
        self.domain = domain or C.Domain()

    def _query(self, t: Target, specs: list[FieldSpec]) -> str:
        vars_ = sorted({s.var for s in specs if s.file == t.kind})
        levs = sorted({s.level for s in specs if s.file == t.kind})
        q = [f"dir=%2F{t.subdir.replace('/', '%2F')}", f"file={t.file}"]
        q += [f"var_{v}=on" for v in vars_]
        q += ["lev_" + lv.replace(" ", "_") + "=on" for lv in levs]
        d = self.domain
        lw, le = d.lon_west % 360.0, d.lon_east % 360.0
        q += ["subregion=", f"toplat={d.lat_north}", f"leftlon={lw}", f"rightlon={le}",
              f"bottomlat={d.lat_south}"]
        return f"{self.base}/filter_gefs_atmos_0p50{t.kind}.pl?" + "&".join(q)

    def exists(self, t: Target) -> bool:  # avoid probing: the filter returns an HTML error for missing files
        return False

    def fetch(self, t: Target, specs: list[FieldSpec]) -> list[bytes]:
        wanted = [s for s in specs if s.file == t.kind]
        if not wanted:
            return []
        try:
            blob = self.http.get_bytes(self._query(t, specs), self.gate)
            if blob[:4] != b"GRIB":
                raise SourceError(f"{self.name}: non-GRIB response for {t.file}: {blob[:80]!r}")
            msgs = split_grib_messages(blob)
        except NotAvailable:
            raise
        except (SourceError, ValueError) as e:
            self._fail()
            raise SourceError(str(e)) from e
        self._ok()
        return msgs


SOURCE_CLASSES = {
    "aws": AwsSource,
    "gcs": GcsSource,
    "nomads_filter": NomadsFilterSource,
    "nomads_http": NomadsHttpSource,
}


class SourceChain:
    """Try each source in order; fall through on NotAvailable *or* SourceError."""

    def __init__(self, sources: list[Source]):
        self.sources = sources
        self.stats: dict[str, int] = {s.name: 0 for s in sources}

    @classmethod
    def build(cls, order=C.DEFAULT_SOURCE_ORDER, gates=None, http=None, bases: dict | None = None,
              domain: C.Domain | None = None):
        http = http or HttpClient()
        bases = bases or {}
        srcs = []
        for name in order:
            klass = SOURCE_CLASSES[name]
            if name == "nomads_filter":
                srcs.append(klass(http, gates, bases.get(name), domain))
            else:
                srcs.append(klass(http, gates, bases.get(name)))
        return cls(srcs)

    def fetch(self, t: Target, specs: list[FieldSpec]) -> tuple[list[bytes], str]:
        """Return (messages, source_name). Raises NotAvailable if no source has it yet,
        SourceError if sources exist but all failed."""
        saw_error, errors = False, []
        for s in self.sources:
            if not s.available():
                continue
            try:
                msgs = s.fetch(t, specs)
                self.stats[s.name] += 1
                return msgs, s.name
            except NotAvailable as e:
                errors.append(f"{s.name}: not available")
                log.debug("%s", e)
            except SourceError as e:
                saw_error = True
                errors.append(f"{s.name}: {e}")
                log.warning("%s failed for %s: %s -- trying next source", s.name, t.file, e)
        if saw_error:
            raise SourceError("; ".join(errors))
        raise NotAvailable("; ".join(errors) or "no source available")

    def exists(self, t: Target) -> bool:
        for s in self.sources:
            if s.available() and s.exists(t):
                return True
        return False
