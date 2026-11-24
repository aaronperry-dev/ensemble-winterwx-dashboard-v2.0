"""Model-side type information (CPOFP, NCEP flags), step-averaged type shares, and the snow liquid-equivalent array."""
import numpy as np
import pytest

from pipeline import config as C
from pipeline import encode
from pipeline.encode import read_bin
from pipeline.member import type_extras
from pipeline.products import cap_freezing_rain, step_average_shares, type_shares
from tests.test_policy import Fake, NX, NY, go
from tests.test_sources import FieldSpec, Handler, SPECS, T, WANT, chain, server  # noqa: F401  (server is a fixture)


def P(ra=0, sn=0, fz=0, pl=0):
    return np.array([ra, sn, fz, pl], np.uint8).reshape(4, 1)


# ------------------------------------------------------------------ freezing-rain cap
def test_cap_moves_excess_freezing_rain_to_snow_when_the_model_says_mostly_frozen():
    out = cap_freezing_rain(P(sn=22, fz=78), np.array([71.0], np.float32))          # the real Rockies case
    assert out[:, 0].tolist() == [0, 71, 29, 0]                                      # liquid share limited to 100 - 71 = 29 %


def test_cap_splits_the_excess_between_snow_and_sleet_by_their_shares():
    out = cap_freezing_rain(P(sn=30, fz=60, pl=10), np.array([100.0], np.float32))
    assert out[:, 0].tolist() == [0, 75, 0, 25]                                      # 60 -> 0; 30:10 split of the excess


def test_cap_is_one_sided_and_ignores_unknown():
    rain = P(ra=100)
    assert (cap_freezing_rain(rain, np.array([60.0], np.float32)) == rain).all()     # rain is never turned into snow
    zr = P(fz=100)
    assert (cap_freezing_rain(zr, np.array([np.nan], np.float32)) == zr).all()       # unknown fraction: leave it alone
    assert (cap_freezing_rain(zr, np.array([0.0], np.float32)) == zr).all()          # model says all liquid: nothing to cap
    assert (cap_freezing_rain(zr, None) == zr).all()


def test_cap_all_excess_goes_to_snow_when_no_frozen_type_yet():
    assert cap_freezing_rain(P(fz=100), np.array([80.0], np.float32))[:, 0].tolist() == [0, 80, 20, 0]


# ------------------------------------------------------------------ step averaging
def test_step_average_is_the_mean_of_start_and_end_shares():
    s = step_average_shares(P(sn=100), P(ra=100))
    assert s[:, 0].tolist() == [50.0, 50.0, 0.0, 0.0]


def test_step_average_uses_the_end_that_has_a_type_and_the_first_step_alone():
    assert step_average_shares(P(sn=100), P())[:, 0].tolist() == [0, 100, 0, 0]       # precipitation just started
    assert step_average_shares(P(), P(fz=100))[:, 0].tolist() == [0, 0, 100, 0]       # just ended
    assert step_average_shares(P(pl=40, sn=40), None)[:, 0].tolist() == [0, 50, 0, 50]
    assert step_average_shares(P(), P())[:, 0].tolist() == [0, 0, 0, 0]


# ------------------------------------------------------------------ NCEP flags -> one code per cell
def test_type_extras_priority_and_missing():
    z = np.zeros((1, 3), np.float32)
    f = {k: z.copy() for k in ("CRAIN", "CSNOW", "CFRZR", "CICEP")}
    f["CRAIN"][0, 0] = 1; f["CSNOW"][0, 1] = 1; f["CFRZR"][0, 2] = 1
    f["CRAIN"][0, 2] = 1                                                             # both set: freezing rain wins
    f["CPOFP"] = np.array([[-50.0, 40.0, 100.0]], np.float32)
    cp, code = type_extras(f)
    assert code[0].tolist() == [0, 1, 2] and np.isnan(cp[0, 0]) and cp[0, 1] == 40.0
    f["CRAIN"][:] = 0; f["CSNOW"][:] = 0; f["CFRZR"][:] = 0
    assert type_extras(f)[1][0].tolist() == [C.NCEP_NONE] * 3
    cp2, code2 = type_extras({"T2": z})                                              # none of the optional fields present
    assert cp2 is None and code2 is None


# ------------------------------------------------------------------ file format
def test_cum_and_ptype_files_round_trip_the_new_arrays(tmp_path):
    M, dom = 3, C.Domain()
    sh = (M, dom.ny, dom.nx)
    q, s, z, p, w = (np.full(sh, v, np.float32) for v in (10.0, 5.0, 1.0, 2.0, 7.0))
    encode.write_cum(tmp_path / "c.bin.gz", 6, 3, q, s, dom, z, p, None, swe_cum_mm=w)
    d = read_bin(tmp_path / "c.bin.gz")
    assert np.allclose(d["swe_mm"], 7.0, atol=0.02) and np.allclose(d["zr_mm"], 1.0) and np.allclose(d["qpf_mm"], 10.0, atol=0.05)
    probs = np.zeros((4, M, dom.ny, dom.nx), np.uint8)
    code = np.full(sh, 1, np.uint8)
    encode.write_ptype(tmp_path / "p.bin.gz", 6, 3, probs, dom, None, ncep_code=code)
    e = read_bin(tmp_path / "p.bin.gz")
    assert e["probs"].shape == (4, M, dom.ny, dom.nx) and (e["ncep"] == 1).all()
    encode.write_ptype(tmp_path / "p0.bin.gz", 6, 3, probs, dom)                     # older layout stays readable
    assert "ncep" not in read_bin(tmp_path / "p0.bin.gz")


# ------------------------------------------------------------------ optional fields
def test_optional_fields_missing_from_the_inventory_are_skipped(server):  # noqa: F811
    spec = SPECS + [FieldSpec("a", "CSNOW", "surface", "CSNOW", True)]               # the mock inventory has no CSNOW
    msgs, src = chain(server, ("aws",)).fetch(T, spec)
    assert msgs == WANT
    with pytest.raises(Exception):
        chain(server, ("aws",)).fetch(T, SPECS + [FieldSpec("a", "CSNOW", "surface", "CSNOW")])   # a required field still fails


# ------------------------------------------------------------------ through the whole run
class Weather(Fake):
    """Rain at f003, then snow from f006; a cold-top cell where Bourgouin says freezing rain but the model says 80 % frozen."""

    def __call__(self, task):
        r = super().__call__(task)
        if r["status"] != "ok":
            return r
        fhr = task[3]
        probs = np.zeros((4, NY, NX), np.uint8)
        probs[0 if fhr == 3 else 1] = 100
        probs[:, 0, 0] = [0, 30, 70, 0]                                              # cell (0, 0): the cold-top case
        cp = np.full((NY, NX), 50.0, np.float32); cp[0, 0] = 80.0
        code = np.full((NY, NX), 0 if fhr == 3 else 1, np.uint8)
        return {**r, "probs": probs, "cpofp": cp, "ncep": code}


def test_run_applies_cap_step_average_and_stores_ncep_and_swe(tmp_path):
    rc, man = go(tmp_path, Weather(), max_fhr=6)
    assert rc == 0
    run = tmp_path / man["base"]
    pt3, pt6 = read_bin(run / "ptype/f003.bin.gz"), read_bin(run / "ptype/f006.bin.gz")
    assert pt3["probs"][:, 0, 0, 0].tolist() == [0, 80, 20, 0]                       # cap: 70 -> 20, excess 50 -> snow
    assert pt6["probs"][:, 0, 5, 5].tolist() == [0, 100, 0, 0]                       # the map shows the end-of-step type
    assert (pt3["ncep"][:, 3, 3] == 0).all() and (pt6["ncep"][:, 3, 3] == 1).all()
    c3, c6 = read_bin(run / "cum/f003.bin.gz"), read_bin(run / "cum/f006.bin.gz")
    # member 0: 3 mm per 3 h.  f003 is all rain -> no snow; f006 averages rain (start) and snow (end) -> half snow
    assert np.isclose(c3["swe_mm"][0, 5, 5], 0.0) and np.isclose(c6["swe_mm"][0, 5, 5], 1.5, atol=0.02)
    assert np.isclose(c6["snow_cm"][0, 5, 5], 1.5 * 12 / 10, atol=0.1)               # SLR 12 here
    assert np.isclose(c6["zr_mm"][0, 0, 0], 2 * 3 * 0.2, atol=0.05)                  # capped to 20 % in both steps, 3 mm each
    # accumulations add up: types sum to the QPF wherever there is a type
    tot = c6["swe_mm"] + c6["zr_mm"] + c6["pl_mm"]
    assert (tot <= c6["qpf_mm"] + 0.1).all()


class Trace(Weather):
    """Cell (7, 7) gets a trace of precipitation (0.02 mm per step) while NCEP's flag says rain."""

    def __call__(self, task):
        r = super().__call__(task)
        if r["status"] == "ok":
            r["apcp"] = r["apcp"].copy(); r["apcp"][7, 7] = 0.02
        return r


def test_ncep_type_is_dropped_where_the_member_is_dry_by_the_bourgouin_wet_mask(tmp_path):
    rc, man = go(tmp_path, Trace(), max_fhr=6)
    assert rc == 0
    pt = read_bin(tmp_path / man["base"] / "ptype/f006.bin.gz")
    assert (pt["probs"][:, :, 7, 7] == 0).all()                                      # Bourgouin: dry
    assert (pt["ncep"][:, 7, 7] == C.NCEP_NONE).all()                                # NCEP follows the same mask
    assert (pt["ncep"][:, 5, 5] == 1).all()                                          # wet cells keep NCEP's type


def test_resuming_matches_an_uninterrupted_run_with_the_new_state(tmp_path):
    import functools, http.server, threading
    live = tmp_path / "live"
    go(live, Weather(), max_fhr=6)                                                   # stops after f006
    handler = functools.partial(http.server.SimpleHTTPRequestHandler, directory=str(live)); handler.log_message = lambda *a, **k: None
    srv = http.server.ThreadingHTTPServer(("127.0.0.1", 0), handler)
    threading.Thread(target=srv.serve_forever, daemon=True).start()
    try:
        rc, resumed = go(tmp_path / "resumed", Weather(), max_fhr=12, reuse_from=f"http://127.0.0.1:{srv.server_address[1]}")
    finally:
        srv.shutdown()
    _, full = go(tmp_path / "full", Weather(), max_fhr=12)
    assert resumed["reused_hours"] == 2
    for rel in ("cum/f009.bin.gz", "cum/f012.bin.gz", "ptype/f012.bin.gz"):
        a, b = read_bin(tmp_path / "resumed" / resumed["base"] / rel), read_bin(tmp_path / "full" / full["base"] / rel)
        for k, v in a.items():
            if isinstance(v, np.ndarray):
                assert np.allclose(v, b[k], atol=0.11 if k == "snow_cm" else 0.06), (rel, k)   # one file-quantization step
