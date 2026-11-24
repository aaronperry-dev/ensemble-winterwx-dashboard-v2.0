"""Source-chain behaviour against a local mock server (no network)."""
import threading
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from urllib.parse import parse_qs, urlparse

import pytest

from pipeline import config as C
from pipeline.http import HttpClient, NotAvailable, SourceError
from pipeline.ratelimit import build_gates
from pipeline.sources import FieldSpec, SourceChain, Target

FAST = {k: {"rate": 5000.0, "burst": 50, "concurrency": 8} for k in ("aws", "gcs", "nomads")}


def msg(tag: bytes, n: int = 24) -> bytes:
    payload = tag.ljust(n, b".")
    return b"GRIB" + b"\0\0\0\2" + (16 + len(payload) + 4).to_bytes(8, "big") + payload + b"7777"


MESSAGES = [("HGT", "500 mb", msg(b"hgt")), ("TMP", "2 m above ground", msg(b"t2")),
            ("RH", "2 m above ground", msg(b"rh2")), ("PRES", "surface", msg(b"ps", 40)),
            ("APCP", "surface", msg(b"apcp"))]
BLOB = b"".join(m for _, _, m in MESSAGES)


def make_idx() -> str:
    pos, rows = 0, []
    for i, (var, lev, m) in enumerate(MESSAGES, 1):
        ft = "0-3 hour acc fcst" if var == "APCP" else "3 hour fcst"
        rows.append(f"{i}:{pos}:d=2026112412:{var}:{lev}:{ft}:ENS=low-res ctl")
        pos += len(m)
    return "\n".join(rows) + "\n"


class Handler(BaseHTTPRequestHandler):
    log: list = []
    behaviour: dict = {}
    honor_range = True

    def log_message(self, *a):
        pass

    def do_GET(self):
        u = urlparse(self.path)
        root = u.path.strip("/").split("/")[0]
        Handler.log.append((root, u.path, self.headers.get("Range"), u.query))
        mode = Handler.behaviour.get(root, "ok")
        if mode in ("403", "404", "500"):
            self.send_response(int(mode)); self.end_headers(); return
        if root == "nomads":                       # GRIB filter CGI: server-side subset, no ranges
            body = MESSAGES[1][2] + MESSAGES[2][2] + MESSAGES[4][2]
        elif u.path.endswith(".idx"):
            body = make_idx().encode()
        else:
            body, rng = BLOB, self.headers.get("Range")
            if rng and Handler.honor_range:
                a, b = rng.split("=")[1].split("-")
                a = int(a); b = int(b) if b else len(BLOB) - 1
                part = BLOB[a:b + 1]
                self.send_response(206); self.send_header("Content-Length", str(len(part))); self.end_headers()
                self.wfile.write(part); return
        self.send_response(200); self.send_header("Content-Length", str(len(body))); self.end_headers()
        self.wfile.write(body)


@pytest.fixture()
def server():
    Handler.log, Handler.behaviour, Handler.honor_range = [], {}, True
    srv = ThreadingHTTPServer(("127.0.0.1", 0), Handler)
    threading.Thread(target=srv.serve_forever, daemon=True).start()
    yield f"http://127.0.0.1:{srv.server_address[1]}"
    srv.shutdown()


T = Target("20261124", "12", "p01", 3, "a")
SPECS = [FieldSpec("a", "TMP", "2 m above ground", "T2"), FieldSpec("a", "RH", "2 m above ground", "RH2"),
         FieldSpec("a", "APCP", "surface", "APCP")]
WANT = [MESSAGES[1][2], MESSAGES[2][2], MESSAGES[4][2]]


def chain(server, order=("aws", "gcs", "nomads_filter", "nomads_http")):
    bases = {"aws": f"{server}/aws", "gcs": f"{server}/gcs", "nomads_filter": f"{server}/nomads",
             "nomads_http": f"{server}/nomadshttp"}
    http = HttpClient(max_retries=1, backoff_base=0.01, backoff_cap=0.02)
    return SourceChain.build(order, build_gates(FAST), http, bases)


def test_aws_primary_uses_byte_ranges_and_merges_adjacent(server):
    msgs, src = chain(server).fetch(T, SPECS)
    assert src == "aws" and msgs == WANT
    ranges = [r for root, p, r, q in Handler.log if r]
    assert len(ranges) == 2                         # TMP+RH adjacent -> one range, APCP -> another
    assert all(r.startswith("bytes=") for r in ranges)
    whole = [p for root, p, r, q in Handler.log if not r and not p.endswith(".idx")]
    assert not whole                                # never downloaded the whole object


def test_falls_through_aws_gcs_to_nomads_filter(server):
    Handler.behaviour = {"aws": "403", "gcs": "500"}
    msgs, src = chain(server).fetch(T, SPECS)
    assert src == "nomads_filter" and msgs == WANT
    roots = [r for r, *_ in Handler.log]
    assert roots.index("aws") < roots.index("gcs") < roots.index("nomads")
    q = next(q for r, p, rng, q in Handler.log if r == "nomads")
    assert "var_TMP=on" in q and "lev_2_m_above_ground=on" in q and "subregion=" in q
    assert "leftlon=225.0" in q and "rightlon=310.0" in q      # lon -135..-50 expressed 0-360


def test_final_fallback_nomads_http_with_ranges(server):
    Handler.behaviour = {"aws": "404", "gcs": "404", "nomads": "500"}
    msgs, src = chain(server).fetch(T, SPECS)
    assert src == "nomads_http" and msgs == WANT
    assert any(r for root, p, r, q in Handler.log if root == "nomadshttp")


def test_not_available_everywhere_is_not_an_error(server):
    Handler.behaviour = {k: "404" for k in ("aws", "gcs", "nomads", "nomadshttp")}
    with pytest.raises(NotAvailable):
        chain(server, ("aws", "gcs", "nomads_http")).fetch(T, SPECS)


def test_all_failing_raises_source_error(server):
    Handler.behaviour = {k: "500" for k in ("aws", "gcs", "nomads", "nomadshttp")}
    with pytest.raises(SourceError):
        chain(server).fetch(T, SPECS)


def test_server_ignoring_range_still_returns_correct_messages(server):
    Handler.honor_range = False
    msgs, src = chain(server).fetch(T, SPECS)
    assert src == "aws" and msgs == WANT


def test_missing_field_in_inventory_falls_through(server):
    bad = SPECS + [FieldSpec("a", "TMP", "850 mb", "T850")]       # not in the mock inventory
    msgs, src = chain(server).fetch(T, bad)
    assert src == "nomads_filter"                                 # filter path does not check the inventory


def test_circuit_breaker_opens_after_repeated_failures(server):
    Handler.behaviour = {"aws": "500"}
    c = chain(server, ("aws", "nomads_filter"))
    for _ in range(5):
        c.fetch(T, SPECS)
    aws = c.sources[0]
    assert not aws.available()
    n = len([1 for r, *_ in Handler.log if r == "aws"])
    c.fetch(T, SPECS)
    assert len([1 for r, *_ in Handler.log if r == "aws"]) == n   # skipped while the circuit is open


def test_non_grib_filter_reply_is_rejected(server):
    class H2(Handler):
        pass
    Handler.behaviour = {"aws": "404", "gcs": "404"}
    c = chain(server, ("nomads_filter",))
    c.sources[0].http.get_bytes = lambda *a, **k: b"<html>Data file is not present</html>"
    with pytest.raises(SourceError):
        c.fetch(T, SPECS)


def hits(*roots):
    return [r for r, *_ in Handler.log if r in roots]


def test_not_posted_yet_does_not_touch_nomads_until_allowed(server):
    Handler.behaviour = {"aws": "404", "gcs": "404"}
    c = chain(server)
    with pytest.raises(NotAvailable):
        c.fetch(T, SPECS, allow_origin=False)
    assert not hits("nomads", "nomadshttp")                       # mirrors merely lack it: NCEP is left alone
    msgs, src = c.fetch(T, SPECS, allow_origin=True)              # after the grace period NOMADS is asked
    assert src == "nomads_filter" and msgs == WANT


def test_mirror_errors_reach_nomads_even_when_origin_is_held_back(server):
    Handler.behaviour = {"aws": "500", "gcs": "404"}
    msgs, src = chain(server).fetch(T, SPECS, allow_origin=False)
    assert src == "nomads_filter" and msgs == WANT                # a broken mirror is a real outage, not "not posted yet"


def test_exhausted_404_budget_is_not_a_failure(server):
    Handler.behaviour = {k: "404" for k in ("aws", "gcs", "nomads", "nomadshttp")}
    gates = build_gates({**FAST, "nomads": {"rate": 5000.0, "burst": 50, "concurrency": 8, "max_404_per_min": 2}})
    http = HttpClient(max_retries=1, backoff_base=0.01, backoff_cap=0.02)
    bases = {"aws": f"{server}/aws", "gcs": f"{server}/gcs", "nomads_http": f"{server}/nomadshttp"}
    c = SourceChain.build(("aws", "gcs", "nomads_http"), gates, http, bases)
    for _ in range(12):                                           # far more misses than the budget allows
        with pytest.raises(NotAvailable) as ei:                   # never a SourceError, so no strike against the member
            c.fetch(T, SPECS)
        assert not isinstance(ei.value, SourceError)
    assert all(s.available() for s in c.sources)                  # and no circuit breaker opened
