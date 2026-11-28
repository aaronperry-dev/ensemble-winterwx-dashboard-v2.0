"""Daily MaxT (12Z-00Z) and MinT (00Z-12Z) from the models' window extremes."""
import datetime as dt
import functools
import http.server
import threading
from pathlib import Path

import numpy as np
import pytest

from pipeline import config as C
from pipeline import encode, grib, tempext
from pipeline.sources import required_fields
from tests.test_policy import NX, NY, Fake, go

UTC = dt.timezone.utc


def cyc(day, hour):
    return dt.datetime(2026, 11, day, hour, tzinfo=UTC)


def gefs_hours(max_fhr=240):
    return C.forecast_hours(max_fhr)


# ------------------------------------------------------------------ which windows a run covers
def test_blocks_belong_to_the_12z_00z_or_00z_12z_window():
    assert tempext.locate(cyc(28, 0), 0) == ("tmin", "2026-11-28")      # 00-06Z
    assert tempext.locate(cyc(28, 0), 6) == ("tmin", "2026-11-28")      # 06-12Z
    assert tempext.locate(cyc(28, 0), 12) == ("tmax", "2026-11-28")     # 12-18Z
    assert tempext.locate(cyc(28, 0), 18) == ("tmax", "2026-11-28")     # 18-24Z
    assert tempext.locate(cyc(28, 0), 24) == ("tmin", "2026-11-29")


def test_a_00z_run_gives_ten_maxima_and_ten_minima_and_no_partial_windows():
    p = tempext.plan(cyc(28, 0), gefs_hours())
    assert list(p["tmax"]) == [f"2026-11-{d}" for d in range(28, 31)] + [f"2026-12-{d:02d}" for d in range(1, 8)]
    assert list(p["tmin"]) == [f"2026-11-{d}" for d in range(28, 31)] + [f"2026-12-{d:02d}" for d in range(1, 8)]
    assert set(p["tmax"].values()) == {2} and set(p["tmin"].values()) == {2}


def test_a_12z_run_has_no_minimum_for_its_first_morning_and_a_final_overnight_one():
    p = tempext.plan(cyc(28, 12), gefs_hours())
    assert list(p["tmax"])[0] == "2026-11-28" and list(p["tmin"])[0] == "2026-11-29"      # 00-12Z of the 28th is already past
    assert list(p["tmax"])[-1] == "2026-12-07" and list(p["tmin"])[-1] == "2026-12-08"
    assert len(p["tmax"]) == 10 and len(p["tmin"]) == 10 and set(p["tmax"].values()) | set(p["tmin"].values()) == {2}


def test_a_half_covered_window_is_kept_but_marked_incomplete():
    p = tempext.plan(cyc(28, 6), gefs_hours(48))                          # 06Z run: the 00-06Z half of the morning is already past
    assert p["tmin"]["2026-11-28"] == 1 and p["tmax"]["2026-11-28"] == 2
    p = tempext.plan(cyc(28, 0), gefs_hours(18))                          # ends at 18Z: the afternoon window is half there
    assert p["tmin"] == {"2026-11-28": 2} and p["tmax"] == {"2026-11-28": 1}


# ------------------------------------------------------------------ building the daily fields
def block(value_max, value_min, w, M=3):
    return {i: (np.full((2, 2), value_max + i, np.float32), np.full((2, 2), value_min - i, np.float32), w) for i in range(M)}


def feed_gefs(d, hours, f=lambda h: (0.0, 0.0)):
    """GEFS-style: 6-hour windows at the even-6 hours; value of block b is 280 + b (max) / 270 - b (min)."""
    out = []
    for h in hours:
        if h % 6 == 0:
            out += d.step(h, block(280.0 + h, 270.0 - h, (h - 6, h)))
        else:
            out += d.step(h, {})
    return out


def test_maximum_is_the_highest_of_the_two_blocks_and_minimum_the_lowest():
    d = tempext.DailyExtremes(cyc(28, 0), gefs_hours(24), 3, 2, 2)
    feed_gefs(d, gefs_hours(24))
    mx, mn = d.items[("tmax", "2026-11-28")], d.items[("tmin", "2026-11-28")]
    assert mx["blocks"] == 2 and mn["blocks"] == 2
    # 12-18Z is block b=18 (value 298), 18-24Z is b=24 (304): max 304 + member index; 00-06Z is b=6 (264), 06-12Z b=12 (258): min 258 - i
    assert mx["t"][0, 0, 0] == 304.0 and mx["t"][2, 0, 0] == 306.0
    assert mn["t"][0, 0, 0] == 258.0 and mn["t"][2, 0, 0] == 256.0
    assert d.available()


def test_ifs_builds_a_block_from_two_three_hour_windows():
    d = tempext.DailyExtremes(cyc(28, 0), [3 * k for k in range(1, 9)], 2, 2, 2)      # f003 .. f024
    for h in range(3, 25, 3):
        d.step(h, {i: (np.full((2, 2), 270.0 + h + (3 if h % 6 == 0 else 0) * i, np.float32), np.full((2, 2), 260.0 - h, np.float32),
                       (h - 3, h)) for i in range(2)})
    mn = d.items[("tmin", "2026-11-28")]
    assert mn["blocks"] == 2 and mn["t"][0, 0, 0] == 260.0 - 12                       # lowest of the 3-hour minima up to F012
    mx = d.items[("tmax", "2026-11-28")]
    assert mx["t"][0, 0, 0] == 270.0 + 24                                              # highest 3-hour maximum, from 21-24Z


def test_a_member_without_extremes_for_one_block_is_left_out_of_that_window_only():
    d = tempext.DailyExtremes(cyc(28, 0), gefs_hours(12), 3, 2, 2)
    got = block(280.0, 270.0, (0, 6))
    del got[1]
    d.step(6, got)
    d.step(12, block(281.0, 271.0, (6, 12)))
    it = d.items[("tmin", "2026-11-28")]
    assert it["mask"] == 0b101 and d.n_members(("tmin", "2026-11-28")) == 2
    assert d.manifest()["tmin"][0]["n"] == 2


def test_windows_that_disagree_are_not_mixed():
    d = tempext.DailyExtremes(cyc(28, 0), gefs_hours(12), 3, 2, 2)
    got = block(280.0, 270.0, (0, 6))
    got[2] = (got[2][0], got[2][1], (3, 6))                                           # a member whose window is not the usual one
    d.step(6, got)
    assert d.items[("tmin", "2026-11-28")]["mask"] == 0b011


def test_missing_cells_are_filled_from_the_other_members():
    a = np.array([[[280.0, np.nan]], [[282.0, 284.0]], [[0.0, 0.0]]], np.float32)
    out = tempext._fill(a, np.array([True, True, False]))
    assert out[0, 0, 1] == 284.0 and out[2].max() == 0.0 and out[1, 0, 0] == 282.0


# ------------------------------------------------------------------ file format
def test_temperature_file_round_trips(tmp_path):
    dom = C.Domain()
    t = np.full((3, dom.ny, dom.nx), 268.15, np.float32)
    t[1, 0, 0] = 310.37
    t[2] = np.nan
    encode.write_temps(tmp_path / "t.bin.gz", 4, 2, t, dom, mask=0b011)
    d = encode.read_bin(tmp_path / "t.bin.gz")
    assert d["kind"] == encode.KIND_TEMP and d["day"] == 4 and d["blocks"] == 2 and d["mask"] == 0b011
    assert d["t_c"][0, 5, 5] == pytest.approx(-5.0, abs=0.03) and d["t_c"][1, 0, 0] == pytest.approx(37.22, abs=0.03)
    assert np.isnan(d["t_c"][2]).all()                                                 # no data stays no data


# ------------------------------------------------------------------ reading the GRIB messages
class Msg:
    def __init__(self, **kw):
        self.k = {"discipline": 0, "typeOfLevel": "heightAboveGround", "level": 2, "stepType": "instant", "startStep": 0, "endStep": 0, **kw}
        self.typeOfLevel, self.level = self.k["typeOfLevel"], self.k["level"]

    def __getitem__(self, key):
        return self.k[key]

    def has_key(self, key):
        return key in self.k


def test_ecmwf_extremes_are_not_mistaken_for_the_instantaneous_2m_temperature():
    inst = Msg(parameterCategory=0, parameterNumber=0)
    mx = Msg(parameterCategory=0, parameterNumber=0, stepType="max", startStep=3, endStep=6)       # mx2t3: same numbers as 2t
    mn = Msg(parameterCategory=0, parameterNumber=0, stepType="min", startStep=144, endStep=150)   # mn2t6
    assert grib.identify(inst) == ("T2", None)
    assert grib.identify(mx) == ("TMAX", (3, 6)) and grib.identify(mn) == ("TMIN", (144, 150))
    assert grib.identify(Msg(parameterCategory=0, parameterNumber=4, stepType="max", startStep=6, endStep=12)) == ("TMAX", (6, 12))   # NCEP / ECCC
    assert grib.identify(Msg(parameterCategory=0, parameterNumber=5, stepType="min", startStep=6, endStep=12)) == ("TMIN", (6, 12))
    assert grib.identify(Msg(parameterCategory=0, parameterNumber=0, typeOfLevel="isobaricInhPa", level=850)) == ("T850", None)


def test_extremes_are_requested_only_when_asked_and_are_optional():
    assert not any(s.var in ("TMAX", "TMIN") for s in required_fields())
    ex = [s for s in required_fields(extremes=True) if s.var in ("TMAX", "TMIN")]
    assert len(ex) == 2 and all(s.optional and s.file == "a" for s in ex)


def test_which_hours_have_extremes():
    assert tempext.wanted("ifs", 3) and tempext.wanted("ifs", 150)
    assert tempext.wanted("gefs", 6) and not tempext.wanted("gefs", 9)
    assert tempext.wanted("geps", 12) and not tempext.wanted("geps", 9)           # GEPS has no file at the 3-hourly hours


# ------------------------------------------------------------------ whole run with a fake model
class Warm(Fake):
    """Window extremes with a recognisable value per block and member: max = 280 + hour/10 + k, min = 270 - hour/10 - k."""

    def __call__(self, task):
        r = super().__call__(task)
        fhr = task[3]
        if r["status"] != "ok" or fhr % 6:
            return r
        k = C.MEMBERS.index(task[2])
        return {**r, "tmax": np.full((NY, NX), 280.0 + fhr / 10 + k, np.float32), "tmin": np.full((NY, NX), 270.0 - fhr / 10 - k, np.float32),
                "twin": np.array([fhr - 6, fhr], np.int32)}


def test_run_writes_a_file_per_complete_window_and_lists_them(tmp_path):
    rc, man = go(tmp_path, Warm(), cycle="2026112800", max_fhr=24)
    assert rc == 0 and man["complete"]
    t = man["temps"]
    assert [e["date"] for e in t["tmax"]] == ["2026-11-28"] and [e["date"] for e in t["tmin"]] == ["2026-11-28"]
    assert "tmax" in man["products"] and "tmin" in man["products"]
    mx = encode.read_bin(tmp_path / man["base"] / t["tmax"][0]["file"])
    mn = encode.read_bin(tmp_path / man["base"] / t["tmin"][0]["file"])
    # 12-18Z = block 12..18 (value at fhr 18), 18-24Z = block at fhr 24: max is the later one
    assert mx["t_c"][0, 3, 3] == pytest.approx(280.0 + 2.4 - 273.15, abs=0.03)
    assert mx["t_c"][30, 3, 3] == pytest.approx(280.0 + 2.4 + 30 - 273.15, abs=0.03)
    assert mn["t_c"][0, 3, 3] == pytest.approx(270.0 - 1.2 - 273.15, abs=0.03)             # 00-12Z: lowest is the f012 block
    assert mx["mask"] == (1 << 31) - 1 and man["temps"]["tmax"][0]["n"] == 31 and man["temps"]["tmax"][0]["blocks"] == 2


def test_a_model_without_extremes_publishes_no_temperature_products(tmp_path):
    rc, man = go(tmp_path, Fake(), cycle="2026112800", max_fhr=24)
    assert rc == 0 and "tmax" not in man["products"]
    assert all(e["file"] is None and e["n"] == 0 for e in man["temps"]["tmax"] + man["temps"]["tmin"])
    assert not list((tmp_path / man["base"]).glob("tmp/*"))


def test_a_dropped_member_is_left_out_of_the_windows_it_missed(tmp_path):
    rc, man = go(tmp_path, Warm(unavailable={"p07": 18}), cycle="2026112800", max_fhr=24)
    k = C.MEMBERS.index("p07")
    mx = man["temps"]["tmax"][0]
    assert mx["n"] == 30 and not (mx["mask"] >> k & 1)                                   # gone from 18Z on: the afternoon window is incomplete
    assert man["temps"]["tmin"][0]["n"] == 31                                            # the overnight window was finished before


def test_resuming_matches_an_uninterrupted_run(tmp_path):
    live = tmp_path / "live"
    go(live, Warm(), cycle="2026112800", max_fhr=18)                                     # stopped after f018: the 12-18Z block is in
    handler = functools.partial(http.server.SimpleHTTPRequestHandler, directory=str(live))
    handler.log_message = lambda *a, **k: None
    srv = http.server.ThreadingHTTPServer(("127.0.0.1", 0), handler)
    threading.Thread(target=srv.serve_forever, daemon=True).start()
    try:
        rc, resumed = go(tmp_path / "resumed", Warm(), cycle="2026112800", max_fhr=24, reuse_from=f"http://127.0.0.1:{srv.server_address[1]}")
    finally:
        srv.shutdown()
    _, full = go(tmp_path / "full", Warm(), cycle="2026112800", max_fhr=24)
    assert resumed["reused_hours"] > 0
    assert resumed["temps"]["tmax"] == full["temps"]["tmax"] and resumed["temps"]["tmin"] == full["temps"]["tmin"]
    for rel in ("tmp/tmax00.bin.gz", "tmp/tmin00.bin.gz"):
        a, b = encode.read_bin(tmp_path / "resumed" / resumed["base"] / rel), encode.read_bin(tmp_path / "full" / full["base"] / rel)
        assert np.allclose(a["t_c"], b["t_c"], atol=0.001, equal_nan=True), rel


def test_resuming_in_the_middle_of_a_block_goes_back_to_its_start(tmp_path):
    from pipeline import reuse
    prior = {"steps": [{"fhr": 3}, {"fhr": 6}, {"fhr": 9}], "complete": False}
    reuse.align_to_blocks(prior)
    assert [s["fhr"] for s in prior["steps"]] == [3, 6]
