"""Failure policy: missing members / data must not stop the run while a quorum remains."""
import argparse
import concurrent.futures as cf
import json
from pathlib import Path

import numpy as np
import pytest

from pipeline import config as C
from pipeline import run as R
from pipeline.encode import read_bin

POL = R.Policy(min_members=20, grace_s=100.0, deadline_s=1000.0)
ALL = {f"m{i}" for i in range(31)}


# ------------------------------------------------------------------ pure decision function
def test_ready_when_everyone_is_in():
    assert R.step_decision(ALL, set(ALL), set(), 10, 0, None, POL)[0] == "ready"


def test_waits_for_stragglers_until_grace_expires():
    have = set(sorted(ALL)[:28])
    assert R.step_decision(ALL, have, set(), 50, 0, None, POL)[0] == "wait"          # 50 s after the last arrival
    act, who, _ = R.step_decision(ALL, have, set(), 101, 0, None, POL)
    assert act == "drop" and set(who) == ALL - have


def test_no_grace_drop_without_a_quorum():
    have = set(sorted(ALL)[:10])
    assert R.step_decision(ALL, have, set(), 5000, 0, None, POL)[0] == "wait"        # nothing says "run is over" yet


def test_hard_deadline_drops_with_quorum_and_stops_without():
    have = set(sorted(ALL)[:25])
    assert R.step_decision(ALL, have, set(), 20, 19, 0, POL)[0] == "wait"
    assert R.step_decision(ALL, have, set(), 1001, 1000, 0, POL)[0] == "drop"
    assert R.step_decision(ALL, set(sorted(ALL)[:10]), set(), 1001, 1000, 0, POL)[0] == "stop"


def test_failed_member_is_dropped_immediately_but_not_below_quorum():
    act, who, _ = R.step_decision(ALL, ALL - {"m30"}, {"m30"}, 1, 0, None, POL)
    assert act == "drop" and who == ["m30"]
    many = {f"m{i}" for i in range(12)}
    assert R.step_decision(ALL, ALL - many, many, 1, 0, None, POL)[0] == "stop"      # 12 failures would leave only 19


# ------------------------------------------------------------------ whole run with a fake model
NY, NX = C.Domain().ny, C.Domain().nx


def window(fhr):
    return (fhr - 3, fhr) if fhr % 6 == 3 else (fhr - 6, fhr)


class Fake:
    """Synthetic member results.  `rules` maps (member, fhr) -> 'unavailable' | 'error'; fhr None = every hour from `from`."""

    def __init__(self, unavailable=(), error=(), ):
        self.unavailable, self.error = dict(unavailable), dict(error)    # member -> first fhr affected (or {(m,fhr)} below)

    def __call__(self, task):
        date, cc, m, fhr = task
        base = {"member": m, "fhr": fhr}
        if m in self.unavailable and fhr >= self.unavailable[m]:
            return {**base, "status": "unavailable", "msg": "not posted"}
        if m in self.error and fhr >= self.error[m]:
            return {**base, "status": "error", "msg": "ValueError: corrupt"}
        k = C.MEMBERS.index(m)
        a, b = window(fhr)
        probs = np.zeros((4, NY, NX), np.uint8)
        probs[0] = 100
        return {**base, "status": "ok", "window": (a, b), "apcp": np.full((NY, NX), (b - a) * (1 + k / 100), np.float32),
                "probs": probs, "slr": np.full((NY, NX), 12.0, np.float32), "sources": {"fake": "fake"}}


def go(tmp_path, fake, **kw):
    args = argparse.Namespace(out=str(tmp_path), cycle="2026112412", members=31, max_fhr=12, workers=2, lookahead=2,
                              wait_minutes=0.004, poll_seconds=0.05, sources="aws", base=None, slr_mode="kuchera",
                              min_members=20, straggler_minutes=0.0025)
    for k, v in kw.items():
        setattr(args, k, v)
    with cf.ThreadPoolExecutor(8) as ex:
        rc = R.run(args, task_fn=fake, executor=ex)
    man = json.loads((Path(tmp_path) / "latest.json").read_text())
    return rc, man


def test_clean_run_has_all_members(tmp_path):
    rc, man = go(tmp_path, Fake())
    assert rc == 0 and man["complete"] and [s["n"] for s in man["steps"]] == [31] * 4 and man["dropped"] == {}


def test_member_missing_from_f006_is_dropped_and_run_finishes(tmp_path):
    rc, man = go(tmp_path, Fake(unavailable={"p07": 6}))
    assert rc == 0 and man["complete"], "run must finish all 10-day hours despite a missing member"
    assert [s["n"] for s in man["steps"]] == [31, 30, 30, 30]
    assert man["dropped"]["p07"]["from_fhr"] == 6
    k = C.MEMBERS.index("p07")
    f3 = read_bin(Path(tmp_path) / man["base"] / "cum/f003.bin.gz")
    f12 = read_bin(Path(tmp_path) / man["base"] / "cum/f012.bin.gz")
    assert f3["mask"] >> k & 1 and not (f12["mask"] >> k & 1)
    assert f12["qpf_mm"][k].max() == 0                                       # dropped member's plane is empty
    j = C.MEMBERS.index("p08")
    assert f12["qpf_mm"][j, 0, 0] == pytest.approx(12 * 1.08, abs=0.06)       # others accumulate normally


def test_repeated_errors_drop_the_member_instead_of_ending_the_run(tmp_path):
    rc, man = go(tmp_path, Fake(error={"p12": 9}))
    assert rc == 0 and man["complete"] and [s["n"] for s in man["steps"]] == [31, 31, 30, 30]
    assert man["dropped"]["p12"]["from_fhr"] == 9


def test_result_buffered_for_a_later_hour_cannot_resurrect_a_dropped_member(tmp_path):
    class OnlyF6(Fake):
        def __call__(self, task):
            if task[2] == "p05" and task[3] == 6:
                return {"member": "p05", "fhr": 6, "status": "unavailable", "msg": "late"}
            return super().__call__(task)
    rc, man = go(tmp_path, OnlyF6())                     # p05 is fine at f009/f012 (already buffered by lookahead)
    assert rc == 0 and man["complete"]
    k = C.MEMBERS.index("p05")
    for s in man["steps"][1:]:
        assert not (s["mask"] >> k & 1), f"p05 must stay excluded at f{s['fhr']:03d}"


def test_too_few_members_stops_but_keeps_what_was_written(tmp_path):
    gone = {m: 6 for m in C.MEMBERS[:15]}
    rc, man = go(tmp_path, Fake(unavailable=gone))
    assert rc == 0 and not man["complete"] and [s["fhr"] for s in man["steps"]] == [3]
    assert (Path(tmp_path) / man["base"] / "cum/f003.bin.gz").exists()


def test_nothing_published_at_all_returns_error_code(tmp_path):
    rc, man = go(tmp_path, Fake(unavailable={m: 3 for m in C.MEMBERS}))
    assert rc == 2 and man["steps"] == [] and not man["complete"]
