import gzip
import struct
import time

import numpy as np
import pytest

from pipeline import config as C
from pipeline.encode import quantize, read_bin, write_cum, write_ptype
from pipeline.idx import IdxEntry, merge_ranges, parse_acc_window, parse_idx, split_grib_messages
from pipeline.products import ApcpAccumulator
from pipeline.ratelimit import HostGate, RateLimiter, build_gates

IDX = """1:0:d=2026112412:HGT:500 mb:3 hour fcst:ENS=low-res ctl
2:1000:d=2026112412:TMP:850 mb:3 hour fcst:ENS=low-res ctl
3:2500:d=2026112412:RH:850 mb:3 hour fcst:ENS=low-res ctl
4:4000:d=2026112412:APCP:surface:0-3 hour acc fcst:ENS=low-res ctl
5:9000:d=2026112412:TMP:700 mb:3 hour fcst:ENS=low-res ctl
"""


def test_parse_idx_and_ranges():
    es = parse_idx(IDX)
    assert [e.var for e in es] == ["HGT", "TMP", "RH", "APCP", "TMP"]
    assert es[1].start == 1000 and es[1].end == 2499 and es[-1].end is None
    # adjacent TMP+RH merge into one request; the non-adjacent APCP stays separate
    assert merge_ranges([es[1], es[2]]) == [(1000, 3999)]
    assert merge_ranges([es[1], es[3]]) == [(1000, 2499), (4000, 8999)]
    assert merge_ranges([es[4]]) == [(9000, None)]
    assert merge_ranges([]) == []


def test_acc_window():
    assert parse_acc_window("0-3 hour acc fcst") == (0, 3)
    assert parse_acc_window("6-12 hour acc fcst") == (6, 12)
    assert parse_acc_window("3 hour fcst") is None


def fake_grib(payload=b"x" * 20):
    total = 16 + len(payload) + 4
    return b"GRIB" + b"\0\0\0\2" + total.to_bytes(8, "big") + payload + b"7777"


def test_split_grib_messages():
    a, b = fake_grib(b"a" * 5), fake_grib(b"b" * 30)
    assert split_grib_messages(a + b) == [a, b]
    with pytest.raises(ValueError):
        split_grib_messages((a + b)[:-3])


def test_token_bucket_paces_requests():
    rl = RateLimiter(rate_per_s=50.0, burst=2)
    t0 = time.monotonic()
    for _ in range(12):
        rl.acquire()
    el = time.monotonic() - t0
    assert 0.15 <= el < 0.6                         # 2 free + 10 paced at 50/s ~= 0.2 s


def test_penalize_blocks_for_retry_after():
    rl = RateLimiter(1000.0, burst=5)
    rl.penalize(0.25)
    t0 = time.monotonic()
    rl.acquire()
    assert time.monotonic() - t0 >= 0.2


def test_404_budget_and_job_wide_split():
    g = HostGate("nomads", 100, 5, 2, max_404_per_min=3)
    assert not g.exhausted_404_budget()
    for _ in range(3):
        g.note_404()
    assert g.exhausted_404_budget()
    gates = build_gates(C.HOST_POLICIES, n_workers=4)
    assert gates["nomads"].limiter.rate == pytest.approx(C.HOST_POLICIES["nomads"]["rate"] / 4)
    assert gates["aws"].limiter.rate == pytest.approx(C.HOST_POLICIES["aws"]["rate"] / 4)


def test_apcp_bucket_accumulation():
    acc = ApcpAccumulator(1, 1, 2)
    one = np.ones((1, 1, 2), np.float32)
    acc.add((0, 3), one * 1)          # f003 = 1
    acc.add((0, 6), one * 3)          # f006 = 3  (0-6 bucket)
    acc.add((6, 9), one * 2)          # f009 = f006 + 2 = 5
    acc.add((6, 12), one * 4)         # f012 = f006 + 4 = 7
    assert acc.cum[9][0, 0, 0] == 5 and acc.cum[12][0, 0, 0] == 7
    assert acc.step(6, 3)[0, 0, 0] == 2 and acc.step(9, 6)[0, 0, 0] == 2 and acc.step(12, 9)[0, 0, 0] == 2
    with pytest.raises(KeyError):
        acc.add((18, 21), one)        # needs f018 first


def test_encode_roundtrip(tmp_path):
    rng = np.random.default_rng(1)
    d = C.Domain()
    q = rng.uniform(0, 80, (31, d.ny, d.nx)).astype(np.float32)
    s = rng.uniform(0, 60, (31, d.ny, d.nx)).astype(np.float32)
    p = tmp_path / "f012.bin.gz"
    write_cum(p, 12, 3, q, s, d)
    assert gzip.open(p).read(4) == b"GWX1"
    r = read_bin(p)
    assert (r["nmem"], r["ny"], r["nx"], r["fhr"]) == (31, d.ny, d.nx, 12)
    assert np.abs(r["qpf_mm"] - q).max() <= C.QPF_MM_PER_COUNT / 2 + 1e-4
    assert np.abs(r["snow_cm"] - s).max() <= C.SNOW_CM_PER_COUNT / 2 + 1e-4
    z = rng.uniform(0, 40, (31, d.ny, d.nx)).astype(np.float32)
    pl = rng.uniform(0, 30, (31, d.ny, d.nx)).astype(np.float32)
    p3 = tmp_path / "f012_x.bin.gz"
    write_cum(p3, 12, 3, q, s, d, z, pl)
    r3 = read_bin(p3)
    assert np.abs(r3["zr_mm"] - z).max() <= C.ZR_MM_PER_COUNT / 2 + 1e-4
    assert np.abs(r3["pl_mm"] - pl).max() <= C.PL_MM_PER_COUNT / 2 + 1e-4
    assert np.abs(r3["qpf_mm"] - q).max() <= C.QPF_MM_PER_COUNT / 2 + 1e-4       # earlier arrays unaffected
    assert "zr_mm" not in r                                                      # 2-array files stay readable
    pr = rng.integers(0, 101, (4, 31, d.ny, d.nx)).astype(np.uint8)
    p2 = tmp_path / "pt.bin.gz"
    write_ptype(p2, 12, 3, pr, d)
    assert np.array_equal(read_bin(p2)["probs"], pr)


def test_quantize_saturates():
    assert quantize(np.array([1e9]), 0.05)[0] == 65535
    assert quantize(np.array([-3.0]), 0.05)[0] == 0


def test_default_run_covers_ten_days_three_hourly():
    hrs = C.forecast_hours()
    assert C.DEFAULT_MAX_FHR == 240 and hrs[0] == 3 and hrs[-1] == 240
    assert len(hrs) == 80 and all(b - a == 3 for a, b in zip(hrs, hrs[1:]))
    assert C.forecast_hours(252)[-2:] == [246, 252]        # 6-hourly beyond f240 if ever requested


def test_type_qpf_split_conserves_liquid():
    from pipeline.products import type_qpf_step_mm
    rng = np.random.default_rng(4)
    q = rng.uniform(0, 12, (3, 4, 5)).astype(np.float32)
    probs = rng.integers(0, 101, (4, 3, 4, 5)).astype(np.uint8)
    probs[:, 0, 0, 0] = 0                                    # a cell with no type information gets nothing
    parts = type_qpf_step_mm(q, probs)
    tot = parts.sum(axis=0)
    has = probs.sum(axis=0) > 0
    assert np.allclose(tot[has], q[has], atol=1e-4)
    assert parts[:, 0, 0, 0].sum() == 0
    only_zr = np.zeros_like(probs); only_zr[2] = 100
    assert np.allclose(type_qpf_step_mm(q, only_zr)[2], q) and type_qpf_step_mm(q, only_zr)[3].max() == 0


def test_site_ships_every_script_the_page_loads():
    import re
    from pipeline.site import ASSETS, ROOT
    html = (ROOT / "index.html").read_text()
    srcs = re.findall(r'<script src="([^"]+)"', html)
    assert srcs and all(s in ASSETS for s in srcs), srcs
    assert "statworker.js" in ASSETS            # loaded by viewer.js via new Worker(...)
    assert "importScripts('ensstats.js')" in (ROOT / "statworker.js").read_text()


def test_site_ships_the_station_list(tmp_path):
    import json
    from pipeline import site
    assert "stations.json" in site.ASSETS and "stations.js" in site.ASSETS
    site.build(tmp_path)
    lst = json.loads((tmp_path / "stations.json").read_text())
    assert len(lst) > 60 and {"id", "lat", "lon", "type"} <= set(lst[0])
