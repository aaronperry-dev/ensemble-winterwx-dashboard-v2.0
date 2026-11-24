"""Command-line entry point.

    python -m pipeline.run --out site/data                 # newest cycle, 31 members, f003..f240
    python -m pipeline.run --out site/data --max-fhr 24    # quick test
    python -m pipeline.run --out site/data --cycle 2026112406 --wait-minutes 0

Failure policy: a member that stays missing (or keeps failing) is dropped from that forecast hour onward and the
run continues with the rest, as long as at least --min-members remain; the files and manifest record which members
each step contains.  See step_decision().

Writes   <out>/<YYYYMMDD_CC>/cum/fHHH.bin.gz     cumulative QPF, snowfall, freezing rain, sleet, every member
         <out>/<YYYYMMDD_CC>/ptype/fHHH.bin.gz   Revised-Bourgouin type probabilities, every member
         <out>/latest.json                        manifest read by the web viewer
"""
from __future__ import annotations

import argparse
import concurrent.futures as cf
import datetime as dt
import json
import logging
import os
import sys
import time
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from . import config as C
from . import encode
from .http import HttpClient, NotAvailable, SourceError
from .member import compute_member_step, fetch_fields
from .products import ApcpAccumulator, snowfall_step_cm, type_qpf_step_mm
from .ratelimit import build_gates
from .sources import SourceChain, Target

log = logging.getLogger("gefs")

# ----------------------------------------------------------------------------- worker side
_W: dict = {}


def init_worker(cfg: dict) -> None:
    gates = build_gates(C.HOST_POLICIES, cfg["n_workers"])
    _W["chain"] = SourceChain.build(order=cfg["order"], gates=gates, http=HttpClient(),
                                    bases=cfg.get("bases"), domain=C.Domain(**cfg["domain"]))
    _W["cfg"] = cfg


def run_task(task: tuple[str, str, str, int]) -> dict:
    date, cycle, member, fhr = task
    cfg, chain = _W["cfg"], _W["chain"]
    base = {"member": member, "fhr": fhr}
    try:
        fields, window, used = fetch_fields(chain, date, cycle, member, fhr, domain=C.Domain(**cfg["domain"]))
        probs, slr = compute_member_step(fields, slr_mode=cfg["slr_mode"])
        return {**base, "status": "ok", "window": window, "apcp": fields["APCP"],
                "probs": probs, "slr": slr, "sources": used}
    except NotAvailable as e:
        return {**base, "status": "unavailable", "msg": str(e)}
    except Exception as e:  # noqa: BLE001 - report everything back to the scheduler
        return {**base, "status": "error", "msg": f"{type(e).__name__}: {e}"}


# ----------------------------------------------------------------------------- cycle choice
def floor_cycle(now: dt.datetime) -> dt.datetime:
    return now.replace(hour=now.hour - now.hour % 6, minute=0, second=0, microsecond=0)


def choose_cycle(chain: SourceChain, now: dt.datetime, tries: int = 4) -> dt.datetime | None:
    """Newest cycle whose last-posted member/product has its first step available."""
    cyc = floor_cycle(now)
    for _ in range(tries):
        t = Target(cyc.strftime("%Y%m%d"), cyc.strftime("%H"), "p30", 3, "b")
        if chain.exists(t):
            return cyc
        log.info("cycle %s not available yet; trying the previous one", cyc.strftime("%Y%m%d%H"))
        cyc -= dt.timedelta(hours=6)
    return None


# ----------------------------------------------------------------------------- failure policy
@dataclass(frozen=True)
class Policy:
    min_members: int = 20       # fewest members worth publishing a step for
    grace_s: float = 900.0      # once a quorum is in, how long stragglers may lag the last arrival before being dropped
    deadline_s: float = 6000.0  # hard cap on waiting for the head step (from the first "not posted yet")


def step_decision(alive: set, have: set, failed: set, now: float, last_progress: float,
                  blocked_since: float | None, pol: Policy):
    """Decide what to do about the forecast hour at the head of the queue.

    alive    members still in the ensemble          have     members whose data for this hour has arrived
    failed   members that exhausted their retries   last_progress  time of the latest arrival (or the hour became head)

    Returns (action, members, reason):
      ("ready", [], "")      every live member is in: write the step
      ("drop", [...], why)   remove these members from the ensemble for this hour and every later one, then re-evaluate
      ("wait", [], "")       keep polling
      ("stop", [], why)      cannot continue (too few members): publish what exists
    A member is never dropped for a single hour only: GEFS accumulations are chained per member, so a gap breaks it.
    """
    missing = alive - have
    if not missing:
        return "ready", [], ""
    hard = missing & failed
    if hard:
        if len(alive) - len(hard) < pol.min_members:
            return "stop", [], f"{len(hard)} member(s) failed and fewer than {pol.min_members} would remain"
        return "drop", sorted(hard), "retries exhausted"
    if len(have) >= pol.min_members and now - last_progress > pol.grace_s:
        return "drop", sorted(missing), f"still missing {pol.grace_s / 60:.0f} min after the rest arrived"
    if blocked_since is not None and now - blocked_since > pol.deadline_s:
        if len(have) >= pol.min_members:
            return "drop", sorted(missing), f"not posted after {pol.deadline_s / 60:.0f} min"
        return "stop", [], f"only {len(have)} of {len(alive)} members posted after {pol.deadline_s / 60:.0f} min"
    return "wait", [], ""


# ----------------------------------------------------------------------------- main loop
def run(args, task_fn=None, executor=None) -> int:
    """task_fn / executor are injectable for tests; by default a process pool runs run_task."""
    domain = C.Domain()
    members = list(C.MEMBERS[: args.members])
    fhrs = C.forecast_hours(args.max_fhr)
    bases = dict(kv.split("=", 1) for kv in (args.base or []))
    order = tuple(args.sources.split(","))
    pol = Policy(min_members=min(getattr(args, "min_members", 20), len(members)),
                 grace_s=getattr(args, "straggler_minutes", 15.0) * 60.0, deadline_s=args.wait_minutes * 60.0)

    probe = SourceChain.build(order=order, gates=build_gates(C.HOST_POLICIES, 1), bases=bases, domain=domain)
    now = dt.datetime.now(dt.timezone.utc)
    if args.cycle == "auto":
        cyc = choose_cycle(probe, now)
        if cyc is None:
            log.error("no GEFS cycle is available from any source")
            return 2
    else:
        cyc = dt.datetime.strptime(args.cycle, "%Y%m%d%H").replace(tzinfo=dt.timezone.utc)
    date, cc = cyc.strftime("%Y%m%d"), cyc.strftime("%H")
    run_id = f"{date}_{cc}"
    out = Path(args.out)
    run_dir = out / run_id
    log.info("GEFS %s: %d members x %d steps, sources=%s, workers=%d, min members=%d", run_id, len(members),
             len(fhrs), ",".join(order), args.workers, pol.min_members)

    cfg = {"n_workers": args.workers, "order": order, "bases": bases, "domain": domain.__dict__,
           "slr_mode": args.slr_mode}
    own_executor = executor is None
    ex = executor or cf.ProcessPoolExecutor(max_workers=args.workers, initializer=init_worker, initargs=(cfg,))
    task = task_fn or run_task

    acc = ApcpAccumulator(len(members), domain.ny, domain.nx)
    cum_x = {k: np.zeros((len(members), domain.ny, domain.nx), np.float32) for k in ("snow", "zr", "pl")}
    steps: list[dict] = []
    src_counts: dict[str, int] = {}
    total_bytes = 0

    alive: set[str] = set(members)
    dropped: dict[str, dict] = {}            # member -> {"from_fhr": ..., "why": ...}
    state = {f: {"res": {}, "err": {}, "failed": set()} for f in fhrs}
    futures: dict[cf.Future, tuple[int, str]] = {}
    retry_q: list[tuple[float, int, str]] = []
    submitted = -1
    head = 0
    blocked_since: float | None = None
    t_start = time.time()
    head_since = t_start
    fatal = None

    def submit(fhr: int, m: str):
        futures[ex.submit(task, (date, cc, m, fhr))] = (fhr, m)

    def write_manifest(complete: bool):
        man = {
            "version": 1, "model": "GEFS",
            "generated": dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
            "init": cyc.strftime("%Y-%m-%dT%H:%M:%SZ"), "run": run_id, "base": f"{run_id}/",
            "complete": complete, "members": members, "domain": domain.to_dict(),
            "ptype_order": ["RA", "SN", "FZRA", "PL"],
            "products": ["qpf", "snow", "zr", "pl", "ptype"],
            "scales": {"qpf_mm_per_count": C.QPF_MM_PER_COUNT, "snow_cm_per_count": C.SNOW_CM_PER_COUNT,
                       "zr_mm_per_count": C.ZR_MM_PER_COUNT, "pl_mm_per_count": C.PL_MM_PER_COUNT},
            "slr_mode": args.slr_mode, "min_precip_mm": C.MIN_PRECIP_MM,
            "min_members": pol.min_members, "dropped": dropped,
            "sources_used": src_counts, "steps": steps,
        }
        out.mkdir(parents=True, exist_ok=True)
        tmp = out / "latest.json.tmp"
        tmp.write_text(json.dumps(man, separators=(",", ":")))
        os.replace(tmp, out / "latest.json")

    try:
        while head < len(fhrs):
            while submitted < min(head + args.lookahead, len(fhrs) - 1):
                submitted += 1
                for m in members:
                    if m in alive:
                        submit(fhrs[submitted], m)
            now_t = time.time()
            for item in [r for r in retry_q if r[0] <= now_t]:
                retry_q.remove(item)
                if item[2] in alive:
                    submit(item[1], item[2])

            if futures:
                done, _ = cf.wait(list(futures), timeout=2.0, return_when=cf.FIRST_COMPLETED)
            else:
                done = set()
                time.sleep(min(0.2, max(0.01, args.poll_seconds / 5)))     # only retry timers pending: do not spin
            for fut in done:
                fhr, m = futures.pop(fut)
                r = fut.result()
                if m not in alive:
                    continue                       # dropped while this task was in flight
                st = state[fhr]
                if r["status"] == "ok":
                    st["res"][m] = r
                    st["t_ok"] = time.time()
                    for k, v in r["sources"].items():
                        src_counts[v] = src_counts.get(v, 0) + 1
                elif r["status"] == "unavailable":
                    retry_q.append((time.time() + args.poll_seconds, fhr, m))
                    if fhr == fhrs[head] and blocked_since is None:
                        blocked_since = time.time()
                else:
                    st["err"][m] = st["err"].get(m, 0) + 1
                    log.warning("f%03d %s failed (%d/3): %s", fhr, m, st["err"][m], r["msg"])
                    if st["err"][m] >= 3:
                        st["failed"].add(m)
                    else:
                        retry_q.append((time.time() + 5.0, fhr, m))

            fhr = fhrs[head]
            st = state[fhr]
            act, who, why = step_decision(alive, set(st["res"]), st["failed"], time.time(),
                                          max(st.get("t_ok", 0.0), head_since), blocked_since, pol)
            if act == "stop":
                fatal = f"f{fhr:03d}: {why}"
                break
            if act == "drop":
                for m in who:
                    alive.discard(m)
                    dropped[m] = {"from_fhr": fhr, "why": why}
                retry_q[:] = [r for r in retry_q if r[2] in alive]
                for stt in state.values():             # results already buffered for later hours must go too:
                    for m in who:                      # a member with a gap can never be accumulated again
                        stt["res"].pop(m, None)
                        stt["failed"].discard(m)
                log.warning("f%03d: dropping %s (%s); continuing with %d of %d members", fhr, ",".join(who), why,
                            len(alive), len(members))
                continue
            if act == "ready":
                total_bytes += finalize_step(fhr, fhrs[head - 1] if head else 0, members, st["res"], acc,
                                             cum_x, run_dir, domain, steps, cyc)
                state[fhr] = {"res": {}, "err": {}, "failed": set()}
                write_manifest(False)
                log.info("f%03d done (%d/%d, %d members)  %.0fs elapsed", fhr, head + 1, len(fhrs), len(alive),
                         time.time() - t_start)
                head += 1
                blocked_since = None
                head_since = time.time()
    finally:
        if own_executor:
            ex.shutdown(wait=False, cancel_futures=True)

    write_manifest(head >= len(fhrs))
    if fatal:
        log.warning("stopping early: %s", fatal)
    if dropped:
        log.warning("members dropped: %s", {m: d["from_fhr"] for m, d in dropped.items()})
    log.info("finished: %d/%d steps, sources=%s, %.1f MB written, %.0fs", len(steps), len(fhrs), src_counts,
             total_bytes / 1e6, time.time() - t_start)
    return 0 if steps else 2


def finalize_step(fhr, prev_fhr, members, res, acc, cum_x, run_dir, domain, steps, cyc) -> int:
    """Write one forecast hour.  `res` holds results only for the members present at this hour; the others keep
    zeros in every array and are flagged absent in the file header's member mask."""
    present = [m for m in members if m in res]
    idx = [members.index(m) for m in present]
    windows = {tuple(res[m]["window"]) for m in present}
    if len(windows) != 1:
        raise SourceError(f"members disagree on APCP window at f{fhr:03d}: {windows}")
    win = windows.pop()
    M, ny, nx = len(members), domain.ny, domain.nx
    apcp = np.zeros((M, ny, nx), np.float32)
    probs = np.zeros((4, M, ny, nx), np.uint8)
    slr = np.full((M, ny, nx), 12.0, np.float32)
    for i, m in zip(idx, present):
        apcp[i], probs[:, i], slr[i] = res[m]["apcp"], res[m]["probs"], res[m]["slr"]
    acc.add(win, apcp)
    qpf = acc.step(fhr, prev_fhr)                                         # (M, ny, nx) mm
    wet = qpf >= C.MIN_PRECIP_MM
    probs = (probs * wet[None]).astype(np.uint8)
    snow = snowfall_step_cm(qpf, probs, slr)
    cum_x["snow"] += snow
    by_type = type_qpf_step_mm(qpf, probs)                                # (4, M, ny, nx) mm liquid-equivalent
    cum_x["zr"] += by_type[2]
    cum_x["pl"] += by_type[3]
    for k in [k for k in acc.cum if k < fhr - 24]:
        del acc.cum[k]
    mask = encode.mask_of(idx)
    absent = [i for i in range(M) if i not in idx]
    cum_q = acc.cum[fhr].copy()
    cz, cp, cs = cum_x["zr"].copy(), cum_x["pl"].copy(), cum_x["snow"].copy()
    for arr in (cum_q, cz, cp, cs):
        arr[absent] = 0.0                                                 # keep dropped members' planes empty
    n = encode.write_cum(run_dir / "cum" / f"f{fhr:03d}.bin.gz", fhr, fhr - prev_fhr, cum_q, cs, domain, cz, cp, mask)
    n += encode.write_ptype(run_dir / "ptype" / f"f{fhr:03d}.bin.gz", fhr, fhr - prev_fhr, probs, domain, mask)
    steps.append({
        "fhr": fhr, "step_h": fhr - prev_fhr, "n": len(present), "mask": mask,
        "valid": (cyc + dt.timedelta(hours=fhr)).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "cum": f"cum/f{fhr:03d}.bin.gz", "ptype": f"ptype/f{fhr:03d}.bin.gz",
    })
    return n


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default="site/data")
    ap.add_argument("--cycle", default="auto", help="auto or YYYYMMDDHH")
    ap.add_argument("--members", type=int, default=31)
    ap.add_argument("--max-fhr", type=int, default=C.DEFAULT_MAX_FHR)
    ap.add_argument("--workers", type=int, default=int(os.environ.get("GEFS_WORKERS", 6)))
    ap.add_argument("--lookahead", type=int, default=2, help="forecast steps downloaded ahead of the writer")
    ap.add_argument("--wait-minutes", type=float, default=90.0, help="hard cap on waiting for not-yet-posted files at one forecast hour")
    ap.add_argument("--min-members", type=int, default=20, help="keep publishing while at least this many members remain")
    ap.add_argument("--straggler-minutes", type=float, default=15.0,
                    help="once --min-members have arrived, drop members that are still missing this long after the last arrival")
    ap.add_argument("--poll-seconds", type=float, default=45.0)
    ap.add_argument("--sources", default=",".join(C.DEFAULT_SOURCE_ORDER))
    ap.add_argument("--base", action="append", help="override a source base URL: name=url (testing)")
    ap.add_argument("--slr-mode", default=C.SLR_MODE, help='"kuchera" or "fixed:12.5"')
    ap.add_argument("-v", "--verbose", action="store_true")
    args = ap.parse_args(argv)
    logging.basicConfig(level=logging.DEBUG if args.verbose else logging.INFO,
                        format="%(asctime)s %(levelname)s %(message)s", datefmt="%H:%M:%S")
    return run(args)


if __name__ == "__main__":
    sys.exit(main())
