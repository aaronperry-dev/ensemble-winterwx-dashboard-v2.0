"""Reusing published hours of the same cycle instead of downloading them from NOAA again."""
import functools
import http.server
import json
import threading
from pathlib import Path

import numpy as np
import pytest

from pipeline import reuse
from pipeline.encode import read_bin
from tests.test_policy import Fake, go


class Boom(Fake):
    """A model that records which (member, hour) tasks were asked for."""

    def __init__(self, *a, **k):
        super().__init__(*a, **k)
        self.asked = []

    def __call__(self, task):
        self.asked.append(task[3])
        return super().__call__(task)


@pytest.fixture()
def site(tmp_path):
    """Serve a directory over HTTP, like the live GitHub Pages data folder."""
    handler = functools.partial(http.server.SimpleHTTPRequestHandler, directory=str(tmp_path / "live"))
    handler.log_message = lambda *a, **k: None
    srv = http.server.ThreadingHTTPServer(("127.0.0.1", 0), handler)
    threading.Thread(target=srv.serve_forever, daemon=True).start()
    yield f"http://127.0.0.1:{srv.server_address[1]}", tmp_path / "live"
    srv.shutdown()


def publish(live: Path, **kw):
    live.mkdir(parents=True, exist_ok=True)
    rc, man = go(live, Fake(), **kw)
    assert rc == 0
    return man


def same_files(a: Path, b: Path):
    for rel in sorted(p.relative_to(a) for p in a.glob("2026*/**/*.bin.gz")):
        x, y = read_bin(a / rel), read_bin(b / rel)
        for k in x:
            if isinstance(x[k], np.ndarray):
                assert np.allclose(x[k], y[k], atol=0.06), (rel, k)       # quantization of the seeded totals only


def test_complete_cycle_is_not_downloaded_again(tmp_path, site):
    url, live = site
    publish(live, max_fhr=12)
    fake = Boom()
    rc, man = go(tmp_path / "new", fake, max_fhr=12, reuse_from=url)
    assert rc == 0 and fake.asked == []                                     # not one member/hour was requested
    assert man["complete"] and man["reused_hours"] == 4 and [s["fhr"] for s in man["steps"]] == [3, 6, 9, 12]
    same_files(live, tmp_path / "new")


def test_partly_published_cycle_continues_at_the_next_hour(tmp_path, site):
    url, live = site
    publish(live, max_fhr=6)                                                # run stopped after f006
    fake = Boom()
    rc, man = go(tmp_path / "new", fake, max_fhr=12, reuse_from=url)
    assert rc == 0 and man["complete"] and man["reused_hours"] == 2
    assert sorted(set(fake.asked)) == [9, 12]                               # only the missing hours were downloaded
    full = tmp_path / "full"
    go(full, Fake(), max_fhr=12)                                            # same result as a run that never stopped
    same_files(full, tmp_path / "new")


def test_dropped_members_stay_dropped_when_resuming(tmp_path, site):
    url, live = site
    live.mkdir(parents=True, exist_ok=True)
    rc, man = go(live, Fake(unavailable={"p30": 6}), max_fhr=6)             # p30 never arrives from f006
    assert "p30" in man["dropped"]
    fake = Boom()
    rc, man2 = go(tmp_path / "new", fake, max_fhr=12, reuse_from=url)
    assert "p30" in man2["dropped"] and man2["steps"][-1]["n"] == 30


def test_changed_calculation_code_forces_a_fresh_download(tmp_path, site, monkeypatch):
    url, live = site
    publish(live, max_fhr=12)
    monkeypatch.setattr(reuse, "code_rev", lambda: "different")
    fake = Boom()
    rc, man = go(tmp_path / "new", fake, max_fhr=12, reuse_from=url)
    assert man["reused_hours"] == 0 and sorted(set(fake.asked)) == [3, 6, 9, 12]


def test_other_cycle_snow_mode_or_force_do_not_reuse(tmp_path, site):
    url, live = site
    publish(live, max_fhr=12)
    for i, kw in enumerate([{"cycle": "2026112418"}, {"slr_mode": "fixed:12.5"}, {"force": True}]):
        fake = Boom()
        rc, man = go(tmp_path / f"n{i}", fake, max_fhr=12, reuse_from=url, **kw)
        assert man["reused_hours"] == 0 and sorted(set(fake.asked)) == [3, 6, 9, 12], kw


def test_a_damaged_published_file_is_refetched_from_that_hour(tmp_path, site):
    url, live = site
    man0 = publish(live, max_fhr=12)
    bad = live / man0["base"] / "cum" / "f009.bin.gz"
    bad.write_bytes(bad.read_bytes()[:500])                                 # truncated download
    fake = Boom()
    rc, man = go(tmp_path / "new", fake, max_fhr=12, reuse_from=url)
    assert man["reused_hours"] == 2 and sorted(set(fake.asked)) == [9, 12] and man["complete"]


def test_no_published_site_just_runs_normally(tmp_path):
    fake = Boom()
    rc, man = go(tmp_path / "new", fake, max_fhr=6, reuse_from="http://127.0.0.1:9")    # nothing listening
    assert rc == 0 and man["reused_hours"] == 0 and sorted(set(fake.asked)) == [3, 6]


def test_code_rev_ignores_transport_code_but_not_the_calculation():
    assert "sources.py" not in reuse.REV_FILES and "http.py" not in reuse.REV_FILES
    assert {"products.py", "bourgouin.py", "encode.py"} <= set(reuse.REV_FILES)
