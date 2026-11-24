"""IFS ENS and GEPS: forecast hours, grid sub-sampling, field conversion, fetchers, whole runs, the synoptic GEFS copy."""
import argparse
import concurrent.futures as cf
import datetime as dt
import json
from pathlib import Path

import numpy as np
import pytest

from pipeline import config as C
from pipeline import families, geps, ifs
from pipeline import run as R
from pipeline.encode import read_bin
from pipeline.grib import Cropper
from pipeline.http import NotAvailable
from pipeline.products import ApcpAccumulator

NY, NX = C.Domain().ny, C.Domain().nx


# ------------------------------------------------------------------ registry
def test_forecast_hours_follow_each_model_output_schedule():
    ifs_h = families.get("ifs").forecast_hours(240)
    assert ifs_h[:3] == [3, 6, 9] and 144 in ifs_h and 147 not in ifs_h and ifs_h[ifs_h.index(144) + 1] == 150 and ifs_h[-1] == 240
    geps_h = families.get("geps").forecast_hours(240)
    assert 192 in geps_h and 195 not in geps_h and geps_h[geps_h.index(192) + 1] == 198 and geps_h[-1] == 240
    assert families.get("gefs").forecast_hours(12) == [3, 6, 9, 12]
    assert families.get("ifs").forecast_hours(9) == [3, 6, 9]


def test_family_members_and_units():
    g, i, p = families.get("gefs"), families.get("ifs"), families.get("geps")
    assert len(g.members) == 31 and len(i.members) == 50 and len(p.members) == 21 and p.members[0] == "c00"
    assert i.members[0] == "e01" and i.members[-1] == "e50"
    assert p.task_units() == ["all"] and len(i.task_units()) == 50
    assert g.uses_cpofp and not i.uses_cpofp and not p.uses_cpofp
    with pytest.raises(SystemExit):
        families.get("icon")


# ------------------------------------------------------------------ 0.25 deg grid -> 0.5 deg domain
class Msg(dict):
    def __init__(self, vals, **kw):
        super().__init__(kw)
        self.values = vals


def global_msg(res, lon_first, jpos=False):
    nj, ni = int(round(180 / res)) + 1, int(round(360 / res))
    lats = 90 - res * np.arange(nj)
    lons = (lon_first + res * np.arange(ni)) % 360
    vals = (lats[:, None] * 1000.0 + lons[None, :]).astype(np.float32)          # every cell encodes its own lat and lon
    kw = dict(iDirectionIncrementInDegrees=res, jDirectionIncrementInDegrees=res, Nj=nj, Ni=ni,
              longitudeOfFirstGridPointInDegrees=lon_first, latitudeOfFirstGridPointInDegrees=90.0,
              latitudeOfLastGridPointInDegrees=-90.0, jScansPositively=0)
    if jpos:
        vals = vals[::-1]
        kw.update(jScansPositively=1, latitudeOfFirstGridPointInDegrees=-90.0, latitudeOfLastGridPointInDegrees=90.0)
    return Msg(vals, **kw), lats, lons


@pytest.mark.parametrize("res,lon_first,jpos", [(0.5, 0.0, True), (0.5, 0.0, False), (0.25, 180.0, False)])
def test_cropper_returns_the_same_half_degree_points_from_any_global_grid(res, lon_first, jpos):
    msg, _, _ = global_msg(res, lon_first, jpos)
    out = Cropper(C.Domain()).crop(msg)
    d = C.Domain()
    assert out.shape == (NY, NX)
    lat = d.lat_north - 0.5 * np.arange(NY)
    lon = (d.lon_west % 360) + 0.5 * np.arange(NX)
    assert np.allclose(out, lat[:, None] * 1000.0 + lon[None, :])


def test_cropper_rejects_a_grid_that_is_not_a_whole_fraction_of_the_domain_spacing():
    msg, _, _ = global_msg(0.5, 0.0)
    msg["iDirectionIncrementInDegrees"] = 0.3
    with pytest.raises(ValueError):
        Cropper(C.Domain()).crop(msg)


# ------------------------------------------------------------------ IFS field conversion and index selection
def test_ifs_fields_convert_units_and_derive_2m_humidity():
    lv = (1000, 850)
    full = lambda v: np.full((NY, NX), v, np.float32)  # noqa: E731
    arrs = {"T1000": full(275.0), "T850": full(270.0), "RH1000": full(104.0), "RH850": full(-3.0),
            "T2": full(283.15), "TD2": full(283.15 - 10.0), "PS": full(101000.0), "APCP": full(0.0125)}
    f, win = ifs.to_fields(arrs, (0, 24), lv)
    assert win == (0, 24) and float(f["APCP"][0, 0]) == pytest.approx(12.5)                 # metres -> mm
    assert f["RH1000"].max() == 100.0 and f["RH850"].min() == 0.0                             # clipped to 0-100
    assert 40 < float(f["RH2"][0, 0]) < 60                                                    # T 10 C, Td 0 C -> about 50 %
    assert "TD2" not in f
    arrs.pop("T850")
    with pytest.raises(Exception):
        ifs.to_fields(arrs, (0, 24), lv)


def test_ifs_index_selection_keeps_only_perturbed_members_and_needed_fields():
    ok = lambda **k: ifs.wanted({"type": "pf", "levtype": "pl", **k}, (1000, 850))  # noqa: E731
    assert ok(param="t", levelist="850") and ok(param="r", levelist="1000")
    assert not ok(param="t", levelist="925") and not ok(param="gh", levelist="850")
    assert ifs.wanted({"type": "pf", "levtype": "sfc", "param": "tp"}, ())
    assert ifs.wanted({"type": "pf", "levtype": "sfc", "param": "2d"}, ())
    assert not ifs.wanted({"type": "cf", "levtype": "sfc", "param": "tp"}, ())
    assert not ifs.wanted({"type": "pf", "levtype": "sol", "param": "sot"}, ())


def fake_grib(n: int) -> bytes:
    body = b"x" * n
    total = 16 + len(body) + 4
    return b"GRIB" + b"\0\0\0\2" + total.to_bytes(8, "big") + body + b"7777"


class FakeHttp:
    def __init__(self, pages: dict):
        self.pages, self.calls = pages, []

    def _find(self, url):
        for k, v in self.pages.items():
            if url.endswith(k):
                return v
        raise NotAvailable(url)

    def get_text(self, url, gate, **kw):
        self.calls.append(url)
        return self._find(url)

    def get_range(self, url, gate, a, b, **kw):
        self.calls.append((url, a, b))
        data = self._find(url)
        return data[a:(None if b is None else b + 1)]


def test_ifs_source_requests_exactly_the_wanted_messages_of_one_member():
    blobs = [fake_grib(10 + i) for i in range(6)]
    offs, pos = [], 0
    for b in blobs:
        offs.append((pos, len(b)))
        pos += len(b)
    names = [("t", "850", "1"), ("r", "850", "2"), ("tp", None, "1"), ("2t", None, "1"), ("sp", None, "1"), ("2d", None, "1")]
    rows = [{"type": "pf", "levtype": "pl" if lev else "sfc", "param": p, "levelist": lev, "number": n,
             "_offset": o, "_length": ln} for (p, lev, n), (o, ln) in zip(names, offs)]
    http = FakeHttp({"-6h-enfo-ef.index": "\n".join(json.dumps(r) for r in rows), "-6h-enfo-ef.grib2": b"".join(blobs)})
    src = ifs.IfsSource("ecmwf", "http://x", object(), http)
    msgs = src.messages("20261124", "12", 6, 1, (850,))
    assert len(msgs) == 5                                                       # member 1 has 5 of the 6 rows; member 2's RH is not ours
    assert [len(m) for m in msgs] == [len(b) for b, n in zip(blobs, names) if n[2] == "1"]
    src.messages("20261124", "12", 6, 1, (850,))
    assert sum(1 for c in http.calls if isinstance(c, str)) == 1                # the index is fetched once and cached
    with pytest.raises(NotAvailable):
        src.messages("20261124", "12", 6, 9, (850,))                            # member absent from the index
    with pytest.raises(NotAvailable):
        src.messages("20261124", "12", 9, 1, (850,))                            # step not posted: index 404


# ------------------------------------------------------------------ GEPS fetcher
def test_geps_file_names_and_member_assembly(monkeypatch):
    levels = (1000, 850)
    specs = geps.file_specs(levels)
    assert ("TMP", "ISBL_1000") in specs and ("RH", "TGL_2m") in specs and ("PRES", "SFC_0") in specs and ("APCP", "SFC_0") in specs
    http = FakeHttp({})
    seen = []

    def get_bytes(url, gate, **kw):
        seen.append(url)
        return fake_grib(8)

    http.get_bytes = get_bytes
    ch = geps.GepsChain({"msc": object()}, C.Domain(), levels, http, base="http://m")
    assert ch.url("20261124", "12", 24, "TMP", "ISBL_0850") == (
        "http://m/20261124/WXO-DD/ensemble/geps/grib2/raw/12/024/CMC_geps-raw_TMP_ISBL_0850_latlon0p5x0p5_2026112412_P024_allmbrs.grib2")

    def fake_decode(msgs, domain):
        var = seen[-1].split("CMC_geps-raw_")[1].split("_latlon")[0]
        key = {"TMP_ISBL_1000": "T1000", "TMP_ISBL_0850": "T850", "RH_ISBL_1000": "RH1000", "RH_ISBL_0850": "RH850",
               "TMP_TGL_2m": "T2", "RH_TGL_2m": "RH2", "PRES_SFC_0": "PS", "APCP_SFC_0": "APCP"}[var]
        n_members = 2 if key == "PS" else 3                                    # one file is short a member
        return {key: {n: np.full((NY, NX), float(n + 1), np.float32) for n in range(n_members)}}, ((0, 24) if key == "APCP" else None)

    monkeypatch.setattr(geps, "decode_members", fake_decode)
    out, win, used = ch.fetch_step("20261124", "12", 24)
    assert win == (0, 24) and used == {"geps": "msc"} and len(seen) == len(specs)
    assert sorted(out) == [0, 1]                                               # only members that have every field
    assert float(out[1]["T850"][0, 0]) == 2.0


# ------------------------------------------------------------------ accumulation that starts at hour 0
def test_accumulator_accepts_windows_that_start_at_zero_even_after_pruning():
    acc = ApcpAccumulator(2, 3, 3)
    acc.add((0, 3), np.full((2, 3, 3), 1.0, np.float32))
    acc.add((0, 6), np.full((2, 3, 3), 2.5, np.float32))
    del acc.cum[0]                                                              # finalize_step prunes old hours
    acc.add((0, 150), np.full((2, 3, 3), 9.0, np.float32))
    assert float(acc.step(150, 6)[0, 0, 0]) == pytest.approx(6.5)


# ------------------------------------------------------------------ whole runs with synthetic results
def ns(tmp_path, family, **kw):
    a = argparse.Namespace(out=str(tmp_path), cycle="2026112412", members=None, max_fhr=12, workers=2, lookahead=2,
                           wait_minutes=0.004, poll_seconds=0.05, sources="aws", base=None, slr_mode="kuchera",
                           min_members=None, straggler_minutes=0.0025, family=family)
    for k, v in kw.items():
        setattr(a, k, v)
    return a


def member_result(m, fhr, scale=1.0, t2=273.15):
    probs = np.zeros((4, NY, NX), np.uint8)
    probs[1] = 100
    return {"member": m, "fhr": fhr, "status": "ok", "window": (0, fhr), "apcp": np.full((NY, NX), fhr * scale, np.float32),
            "probs": probs, "slr": np.full((NY, NX), 10.0, np.float32), "cpofp": None,
            "t2": np.full((NY, NX), t2, np.float32), "sources": {"x": "x"}}


def read_manifest(tmp_path, name="latest.json"):
    return json.loads((Path(tmp_path) / name).read_text())


def test_geps_run_writes_one_task_per_hour_with_all_members(tmp_path):
    fam = families.get("geps")
    calls = []

    def task(t):
        date, cc, unit, fhr, *_ = t
        calls.append((unit, fhr))
        return {"member": unit, "fhr": fhr, "status": "ok", "window": (0, fhr), "sources": {"geps": "msc"},
                "members": {m: member_result(m, fhr, 1 + i / 100) for i, m in enumerate(fam.members)}}

    with cf.ThreadPoolExecutor(4) as ex:
        rc = R.run(ns(tmp_path, "geps", max_fhr=9), task_fn=task, executor=ex)
    man = read_manifest(tmp_path)
    assert rc == 0 and man["complete"] and man["family"] == "geps" and man["model"] == "GEPS"
    assert man["members"] == list(fam.members) and [s["n"] for s in man["steps"]] == [21, 21, 21]
    assert {u for u, _ in calls} == {"all"} and len(calls) == 3                 # one task per hour, not per member
    d = read_bin(Path(tmp_path) / man["base"] / "cum/f009.bin.gz")
    assert d["qpf_mm"].shape[0] == 21 and d["qpf_mm"][5, 0, 0] == pytest.approx(9 * 1.05, abs=0.06)
    assert d["snow_cm"][0, 0, 0] > 0                                            # snow accumulates from the cumulative window (0, fhr)


def test_ifs_run_handles_six_hourly_steps_and_cumulative_windows(tmp_path):
    fam = families.get("ifs")

    def task(t):
        date, cc, m, fhr, *_ = t
        return member_result(m, fhr, 0.1 + 0.01 * fam.members.index(m))

    with cf.ThreadPoolExecutor(4) as ex:
        rc = R.run(ns(tmp_path, "ifs", max_fhr=156, members=4, min_members=2), task_fn=task, executor=ex)
    man = read_manifest(tmp_path)
    fhrs = [s["fhr"] for s in man["steps"]]
    assert rc == 0 and man["complete"] and fhrs[-3:] == [144, 150, 156] and len(fhrs) == 50
    assert [s["step_h"] for s in man["steps"][-3:]] == [3, 6, 6]
    d = read_bin(Path(tmp_path) / man["base"] / "cum/f156.bin.gz")
    assert d["qpf_mm"][0, 0, 0] == pytest.approx(156 * 0.1, abs=0.06)           # the run never lost its zero base
    assert man["members"] == list(fam.members[:4]) and man["levels"] == [1000, 925, 850, 700, 500]


def test_a_missing_geps_hour_is_waited_for_then_the_run_publishes_what_it_has(tmp_path):
    fam = families.get("geps")

    def task(t):
        date, cc, unit, fhr, *_ = t
        if fhr >= 6:
            return {"member": unit, "fhr": fhr, "status": "unavailable", "msg": "not posted"}
        return {"member": unit, "fhr": fhr, "status": "ok", "window": (0, fhr), "sources": {"geps": "msc"},
                "members": {m: member_result(m, fhr) for m in fam.members}}

    with cf.ThreadPoolExecutor(4) as ex:
        rc = R.run(ns(tmp_path, "geps", max_fhr=9), task_fn=task, executor=ex)
    man = read_manifest(tmp_path)
    assert rc == 0 and not man["complete"] and [s["fhr"] for s in man["steps"]] == [3]


def test_alternate_manifest_name(tmp_path):
    with cf.ThreadPoolExecutor(2) as ex:
        R.run(ns(tmp_path, "ifs", max_fhr=6, members=3, min_members=2, manifest="synoptic.json"),
              task_fn=lambda t: member_result(t[2], t[3]), executor=ex)
    assert (Path(tmp_path) / "synoptic.json").exists() and not (Path(tmp_path) / "latest.json").exists()


# ------------------------------------------------------------------ cycles
class Probe:
    def __init__(self, have):
        self.have, self.asked = have, []

    def exists(self, *a):
        self.asked.append(a)
        return a[:2] in self.have


def test_ifs_and_geps_cycles_are_00_and_12z_only():
    now = dt.datetime(2026, 11, 24, 20, 0, tzinfo=dt.timezone.utc)
    p = Probe({("20261124", "12")})
    assert R.choose_cycle(p, now, family=families.get("ifs")) == dt.datetime(2026, 11, 24, 12, tzinfo=dt.timezone.utc)
    p = Probe({("20261124", "00")})
    assert R.choose_cycle(p, now, family=families.get("geps")) == dt.datetime(2026, 11, 24, 0, tzinfo=dt.timezone.utc)
    assert p.asked == [("20261124", "12", 3), ("20261124", "00", 3)]
    assert R.floor_cycle(dt.datetime(2026, 11, 24, 17, 5, tzinfo=dt.timezone.utc)).hour == 12


# ------------------------------------------------------------------ GEFS synoptic copy
def test_synoptic_copy_is_just_the_manifest_when_the_target_is_the_main_cycle(tmp_path):
    (tmp_path / "latest.json").write_text('{"run":"x"}')
    args = ns(tmp_path, "gefs")
    args.resolved_cycle = dt.datetime(2026, 11, 24, 12, tzinfo=dt.timezone.utc)
    assert R.publish_synoptic(args, target_fn=lambda main: main) == 0
    assert (tmp_path / "synoptic.json").read_text() == '{"run":"x"}'


def test_synoptic_pass_builds_the_cycle_the_other_models_have(tmp_path, monkeypatch):
    seen = {}
    monkeypatch.setattr(R, "run", lambda a, **k: seen.update(cycle=a.cycle, manifest=a.manifest) or 0)
    args = ns(tmp_path, "gefs")
    args.resolved_cycle = dt.datetime(2026, 11, 24, 18, tzinfo=dt.timezone.utc)
    assert R.publish_synoptic(args, target_fn=lambda main: dt.datetime(2026, 11, 24, 12, tzinfo=dt.timezone.utc)) == 0
    assert seen == {"cycle": "2026112412", "manifest": "synoptic.json"}


def test_synoptic_target_picks_the_newest_cycle_that_ifs_and_geps_both_have(monkeypatch):
    have = {"ifs": {("20261124", "00"), ("20261124", "12")}, "geps": {("20261124", "12")}}
    monkeypatch.setattr(R, "build_chain", lambda fam, *a, **k: Probe(have[fam.key]))
    now = dt.datetime(2026, 11, 24, 22, 30, tzinfo=dt.timezone.utc)
    main = dt.datetime(2026, 11, 24, 18, tzinfo=dt.timezone.utc)
    assert R.synoptic_target(main, now) == dt.datetime(2026, 11, 24, 12, tzinfo=dt.timezone.utc)
    have["geps"] = {("20261124", "00")}                                          # GEPS lags: IFS 12Z is not matched, 00Z is
    have["ifs"] = {("20261124", "00")}
    assert R.synoptic_target(main, now) == dt.datetime(2026, 11, 24, 0, tzinfo=dt.timezone.utc)
    have["ifs"], have["geps"] = set(), set()                                     # nothing reachable: the 12Z before the main cycle
    assert R.synoptic_target(main, now) == dt.datetime(2026, 11, 24, 12, tzinfo=dt.timezone.utc)
    assert R.synoptic_target(dt.datetime(2026, 11, 24, 6, tzinfo=dt.timezone.utc), now) == dt.datetime(2026, 11, 24, 0, tzinfo=dt.timezone.utc)


def test_an_explicit_off_cycle_is_floored_for_ifs(tmp_path):
    with cf.ThreadPoolExecutor(2) as ex:
        R.run(ns(tmp_path, "ifs", cycle="2026112418", max_fhr=3, members=3, min_members=2),
              task_fn=lambda t: member_result(t[2], t[3]), executor=ex)
    assert read_manifest(tmp_path)["run"] == "20261124_12"


# ------------------------------------------------------------------ reuse works per family and per manifest name
from tests.test_reuse import site  # noqa: E402,F401  (fixture: serves a directory like the live Pages site)


def test_reuse_reads_the_manifest_of_the_same_family_only(tmp_path, site):
    url, live = site
    with cf.ThreadPoolExecutor(2) as ex:
        R.run(ns(live, "ifs", max_fhr=9, members=3, min_members=2, manifest="synoptic.json"),
              task_fn=lambda t: member_result(t[2], t[3]), executor=ex)
    asked = []

    def task(t):
        asked.append(t[3])
        return member_result(t[2], t[3])

    with cf.ThreadPoolExecutor(2) as ex:                      # same family, same manifest name: nothing is downloaded
        rc = R.run(ns(tmp_path / "a", "ifs", max_fhr=9, members=3, min_members=2, manifest="synoptic.json", reuse_from=url),
                   task_fn=task, executor=ex)
    man = read_manifest(tmp_path / "a", "synoptic.json")
    assert rc == 0 and asked == [] and man["reused_hours"] == 3 and man["complete"]
    with cf.ThreadPoolExecutor(2) as ex:                      # the published file is synoptic.json, not latest.json: no reuse
        R.run(ns(tmp_path / "b", "ifs", max_fhr=9, members=3, min_members=2, reuse_from=url), task_fn=task, executor=ex)
    assert sorted(set(asked)) == [3, 6, 9]
    with cf.ThreadPoolExecutor(2) as ex:                      # another family never reuses it
        R.run(ns(tmp_path / "c", "geps", max_fhr=6, manifest="synoptic.json", reuse_from=url),
              task_fn=lambda t: {"member": "all", "fhr": t[3], "status": "ok", "window": (0, t[3]), "sources": {},
                                 "members": {m: member_result(m, t[3]) for m in families.get("geps").members}}, executor=ex)
    assert read_manifest(tmp_path / "c", "synoptic.json")["reused_hours"] == 0


# ------------------------------------------------------------------ 50-member masks (IFS) do not fit in 32 bits
def test_member_mask_round_trips_for_fifty_members(tmp_path):
    from pipeline import encode
    z = np.zeros((50, 3, 3), np.float32)
    present = [k for k in range(50) if k not in (3, 40, 49)]
    mask = encode.mask_of(present)
    assert mask > 2 ** 32
    encode.write_cum(tmp_path / "c.bin.gz", 3, 3, z, z, C.Domain(), z, z, mask, swe_cum_mm=z)
    encode.write_ptype(tmp_path / "p.bin.gz", 3, 3, np.zeros((4, 50, 3, 3), np.uint8), C.Domain(), mask)
    assert encode.read_bin(tmp_path / "c.bin.gz")["mask"] == mask and encode.read_bin(tmp_path / "p.bin.gz")["mask"] == mask
    full = encode.full_mask(50)
    encode.write_cum(tmp_path / "d.bin.gz", 3, 3, z, z, C.Domain(), z, z, None, swe_cum_mm=z)
    assert encode.read_bin(tmp_path / "d.bin.gz")["mask"] == full
