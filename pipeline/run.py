"""Command-line entry point.

    python -m pipeline.run --out site/data/gefs                       # newest GEFS cycle, 31 members, f003..f240
    python -m pipeline.run --family ifs  --out site/data/ifs          # ECMWF IFS ENS (50 members, 00/12Z)
    python -m pipeline.run --family geps --out site/data/geps         # ECCC GEPS (21 members, 00/12Z)
    python -m pipeline.run --out site/data/gefs --max-fhr 24          # quick test
    python -m pipeline.run --out site/data/gefs --cycle 2026112406 --wait-minutes 0

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
import shutil
import sys
import time
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from . import config as C
from . import encode, families, geps, ifs, reuse, tempext
from .http import HttpClient, NotAvailable, SourceError
from .member import compute_member_step, cpofp_from, fetch_fields
from .products import (ApcpAccumulator, cap_freezing_rain, snowfall_step_cm, step_average_shares, surface_temperature_gates,
                       type_qpf_step_mm, type_shares)
from .ratelimit import build_gates
from .sources import SourceChain, Target

log = logging.getLogger("gefs")

# ----------------------------------------------------------------------------- worker side
_W: dict = {}


def build_chain(fam: families.Family, gates, order, bases, domain: C.Domain):
    """The fetcher for one family: a source chain (GEFS), IFS open data, or the GEPS datamart."""
    if fam.key == "ifs":
        return ifs.IfsChain(ifs.build_sources(gates, HttpClient(), bases), domain, fam.levels)
    if fam.key == "geps":
        return geps.GepsChain(gates, domain, fam.levels, HttpClient(), (bases or {}).get("msc"))
    return SourceChain.build(order=order, gates=gates, http=HttpClient(), bases=bases, domain=domain)


def init_worker(cfg: dict) -> None:
    fam = families.get(cfg.get("family", "gefs"))
    gates = build_gates(C.HOST_POLICIES, cfg["n_workers"])
    _W["chain"] = build_chain(fam, gates, cfg["order"], cfg.get("bases"), C.Domain(**cfg["domain"]))
    _W["cfg"] = cfg
    _W["fam"] = fam


def _member_result(base: dict, fam: families.Family, fields: dict, window, used, slr_mode: str) -> dict:
    probs, slr = compute_member_step(fields, levels=fam.levels, slr_mode=slr_mode)
    return {**base, "status": "ok", "window": window, "apcp": fields["APCP"], "probs": probs, "slr": slr,
            "cpofp": cpofp_from(fields) if fam.uses_cpofp else None, "t2": fields["T2"], "sources": used,
            "tmax": fields.get("TMAX"), "tmin": fields.get("TMIN"), "twin": fields.get(tempext.TWIN)}


def run_task(task: tuple[str, str, str, int]) -> dict:
    date, cycle, unit, fhr, *rest = task
    allow_origin = rest[0] if rest else True
    cfg, chain, fam = _W["cfg"], _W["chain"], _W.get("fam") or families.get(_W["cfg"].get("family", "gefs"))
    base = {"member": unit, "fhr": fhr}
    try:
        if fam.key == "geps":
            by_n, window, used = chain.fetch_step(date, cycle, fhr, allow_origin=allow_origin)
            mem = {fam.members[n]: _member_result({"member": fam.members[n], "fhr": fhr}, fam, f, window, used, cfg["slr_mode"])
                   for n, f in by_n.items() if n < len(fam.members)}
            return {**base, "status": "ok", "window": window, "sources": used, "members": mem}
        if fam.key == "ifs":
            fields, window, used = chain.fetch_member_step(date, cycle, unit, fhr, allow_origin=allow_origin)
        else:
            fields, window, used = fetch_fields(chain, date, cycle, unit, fhr, levels=fam.levels,
                                               domain=C.Domain(**cfg["domain"]), allow_origin=allow_origin,
                                               extremes=tempext.wanted(fam.key, fhr))
        return _member_result(base, fam, fields, window, used, cfg["slr_mode"])
    except NotAvailable as e:
        return {**base, "status": "unavailable", "msg": str(e)}
    except Exception as e:  # noqa: BLE001 - report everything back to the scheduler
        return {**base, "status": "error", "msg": f"{type(e).__name__}: {e}"}


def flatten_results(res: dict) -> dict:
    """Per-task results -> per-member results (a GEPS task carries all of its members)."""
    out: dict = {}
    for unit, r in res.items():
        if "members" in r:
            out.update(r["members"])
        else:
            out[unit] = r
    return out


# ----------------------------------------------------------------------------- cycle choice
def floor_cycle(now: dt.datetime, every: int = 6) -> dt.datetime:
    return now.replace(hour=now.hour - now.hour % every, minute=0, second=0, microsecond=0)


def choose_cycle(chain, now: dt.datetime, tries: int = 4, family: families.Family | None = None) -> dt.datetime | None:
    """Newest cycle whose first step is available (GEFS: probes the last-posted member, p30)."""
    fam = family or families.get("gefs")
    every = 12 if fam.cycle_hours == (0, 12) else 6
    cyc = floor_cycle(now, every)
    for _ in range(tries):
        d, c = cyc.strftime("%Y%m%d"), cyc.strftime("%H")
        ok = chain.exists(d, c, 3) if fam.key in ("ifs", "geps") else chain.exists(Target(d, c, "p30", 3, "b"))
        if ok:
            return cyc
        log.info("%s cycle %s not available yet; trying the previous one", fam.label, cyc.strftime("%Y%m%d%H"))
        cyc -= dt.timedelta(hours=every)
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
    fam = families.get(getattr(args, "family", None) or "gefs")
    members = list(fam.members[: args.members]) if getattr(args, "members", None) else list(fam.members)
    units = fam.task_units(members)                    # what one fetch task covers: a member, or (GEPS) the whole hour
    fhrs = fam.forecast_hours(args.max_fhr if args.max_fhr else None)
    manifest_name = getattr(args, "manifest", None) or "latest.json"
    bases = dict(kv.split("=", 1) for kv in (args.base or []))
    order = tuple(args.sources.split(","))
    pol = Policy(min_members=min(getattr(args, "min_members", None) or fam.min_members, len(units)),
                 grace_s=getattr(args, "straggler_minutes", 15.0) * 60.0, deadline_s=args.wait_minutes * 60.0)

    now = dt.datetime.now(dt.timezone.utc)
    if args.cycle == "auto":
        probe = build_chain(fam, build_gates(C.HOST_POLICIES, 1), order, bases, domain)
        cyc = choose_cycle(probe, now, family=fam)
        if cyc is None:
            log.error("no %s cycle is available from any source", fam.label)
            return 2
    else:
        cyc = dt.datetime.strptime(args.cycle, "%Y%m%d%H").replace(tzinfo=dt.timezone.utc)
        if cyc.hour not in fam.cycle_hours:
            cyc = floor_cycle(cyc, 12)
            log.warning("%s runs at %s Z only: using %s", fam.label, "/".join(f"{h:02d}" for h in fam.cycle_hours), cyc.strftime("%Y%m%d%H"))
    date, cc = cyc.strftime("%Y%m%d"), cyc.strftime("%H")
    run_id = f"{date}_{cc}"
    out = Path(args.out)
    run_dir = out / run_id
    args.resolved_cycle = cyc
    log.info("%s %s: %d members x %d steps, sources=%s, workers=%d, min members=%d", fam.label, run_id, len(members),
             len(fhrs), ",".join(order) if fam.key == "gefs" else fam.key, args.workers, pol.min_members)

    cfg = {"n_workers": args.workers, "order": order, "bases": bases, "domain": domain.__dict__,
           "slr_mode": args.slr_mode, "family": fam.key}
    prior = None                              # hours of this same cycle that are already published: do not download them again
    if getattr(args, "reuse_from", None) and not getattr(args, "force", False):
        prior = reuse.load_prior(args.reuse_from, run_dir, cyc, members, domain, args.slr_mode, fhrs,
                                 family=fam.key, manifest=manifest_name)
    own_executor = executor is None
    ex = executor or (None if prior and prior["complete"] else
                      cf.ProcessPoolExecutor(max_workers=args.workers, initializer=init_worker, initargs=(cfg,)))
    task = task_fn or run_task

    if prior and not prior["complete"]:
        reuse.align_to_blocks(prior, cyc, fhrs)                 # resume at the end of a 6-hour block so no half-built block is lost
        if not prior["steps"]:
            prior = None
    acc = ApcpAccumulator(len(members), domain.ny, domain.nx)
    temps = tempext.DailyExtremes(cyc, fhrs, len(members), domain.ny, domain.nx)       # daily max / min 2 m temperature
    cum_x = {k: np.zeros((len(members), domain.ny, domain.nx), np.float32) for k in ("snow", "zr", "pl", "swe")}
    carry: dict = {"probs": None}            # previous step's type probabilities (for the start-of-step end of the average)
    steps: list[dict] = []
    src_counts: dict[str, int] = {}
    total_bytes = 0

    alive: set[str] = set(units)
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
    if prior:
        reuse.seed(prior, run_dir, acc, cum_x, carry, temps=temps)
        steps.extend(prior["steps"])
        dropped.update(prior["dropped"])
        alive -= set(dropped)
        src_counts.update(prior["sources_used"])
        head, submitted = len(steps), len(steps) - 1
        log.info("reusing %d of %d hours already published for %s (no NOAA downloads for them)%s", head, len(fhrs), run_id,
                 "" if head < len(fhrs) else "; nothing left to fetch")

    lookahead = args.lookahead if fam.units == "member" else max(args.lookahead, args.workers)   # one task is a whole hour
    first_miss: dict[tuple[int, str], float] = {}     # (hour, member) -> when the mirrors first said "not posted yet"
    origin_after = getattr(args, "origin_after_seconds", 120.0)

    def submit(fhr: int, m: str):
        # NOMADS is asked only for the hour we are waiting on, and only once the mirrors have lacked the file a while
        allow = head < len(fhrs) and fhr == fhrs[head] and (time.time() - first_miss.get((fhr, m), time.time())) >= origin_after
        futures[ex.submit(task, (date, cc, m, fhr, allow))] = (fhr, m)

    def write_manifest(complete: bool):
        man = {
            "version": 1, "model": fam.label, "family": fam.key, "levels": list(fam.levels), "credit": fam.credit,
            "generated": dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
            "init": cyc.strftime("%Y-%m-%dT%H:%M:%SZ"), "run": run_id, "base": f"{run_id}/",
            "complete": complete, "members": members, "domain": domain.to_dict(),
            "ptype_order": ["RA", "SN", "FZRA", "PL"],
            "products": ["qpf", "snow", "zr", "pl", "ptype"] + (["tmax", "tmin"] if temps.available() else []),
            "scales": {"qpf_mm_per_count": C.QPF_MM_PER_COUNT, "snow_cm_per_count": C.SNOW_CM_PER_COUNT,
                       "zr_mm_per_count": C.ZR_MM_PER_COUNT, "pl_mm_per_count": C.PL_MM_PER_COUNT,
                       "swe_mm_per_count": C.SWE_MM_PER_COUNT},
            "slr_mode": args.slr_mode, "min_precip_mm": C.MIN_PRECIP_MM, "code_rev": reuse.code_rev(),
            "reused_hours": len(prior["steps"]) if prior else 0,
            "min_members": pol.min_members, "dropped": dropped,
            "sources_used": src_counts, "steps": steps, "temps": temps.manifest(),
        }
        out.mkdir(parents=True, exist_ok=True)
        tmp = out / (manifest_name + ".tmp")
        tmp.write_text(json.dumps(man, separators=(",", ":")))
        os.replace(tmp, out / manifest_name)

    try:
        while head < len(fhrs):
            while submitted < min(head + lookahead, len(fhrs) - 1):
                submitted += 1
                for m in units:
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
                    first_miss.pop((fhr, m), None)
                    for k, v in r["sources"].items():
                        src_counts[v] = src_counts.get(v, 0) + 1
                elif r["status"] == "unavailable":
                    first_miss.setdefault((fhr, m), time.time())
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
                total_bytes += finalize_step(fhr, fhrs[head - 1] if head else 0, members, flatten_results(st["res"]), acc,
                                             cum_x, run_dir, domain, steps, cyc, carry, temps)
                state[fhr] = {"res": {}, "err": {}, "failed": set()}
                write_manifest(False)
                log.info("f%03d done (%d/%d, %d members)  %.0fs elapsed", fhr, head + 1, len(fhrs), steps[-1]["n"],
                         time.time() - t_start)
                head += 1
                blocked_since = None
                head_since = time.time()
    finally:
        if own_executor and ex is not None:
            ex.shutdown(wait=False, cancel_futures=True)

    write_manifest(head >= len(fhrs))
    if fatal:
        log.warning("stopping early: %s", fatal)
    if dropped:
        log.warning("members dropped: %s", {m: d["from_fhr"] for m, d in dropped.items()})
    log.info("finished: %d/%d steps, sources=%s, %.1f MB written, %.0fs", len(steps), len(fhrs), src_counts,
             total_bytes / 1e6, time.time() - t_start)
    return 0 if steps else 2


def finalize_step(fhr, prev_fhr, members, res, acc, cum_x, run_dir, domain, steps, cyc, carry=None, temps=None) -> int:
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
    cpofp = np.full((M, ny, nx), np.nan, np.float32)                      # model % frozen precipitation (NaN = unknown)
    t2c = np.full((M, ny, nx), np.nan, np.float32)                        # 2 m temperature, deg C (NaN = unknown)
    for i, m in zip(idx, present):
        apcp[i], probs[:, i], slr[i] = res[m]["apcp"], res[m]["probs"], res[m]["slr"]
        if res[m].get("cpofp") is not None:
            cpofp[i] = res[m]["cpofp"]
        if res[m].get("t2") is not None:
            t2c[i] = res[m]["t2"] - 273.15
    acc.add(win, apcp)
    qpf = acc.step(fhr, prev_fhr)                                         # (M, ny, nx) mm
    wet = qpf >= C.MIN_PRECIP_MM
    probs = surface_temperature_gates(probs, t2c if np.isfinite(t2c).any() else None)       # 1. 2 m temperature checks
    probs = cap_freezing_rain(probs, cpofp if np.isfinite(cpofp).any() else None)            # 2. the model's own liquid fraction
    probs = (probs * wet[None]).astype(np.uint8)                          # what the ptype map shows: end-of-step type
    carry = carry if carry is not None else {"probs": None}
    shares = step_average_shares(probs, carry["probs"])                   # what the accumulations use: mean over the step
    carry["probs"] = probs.copy()
    snow = snowfall_step_cm(qpf, shares, slr)
    cum_x["snow"] += snow
    cum_x["swe"] += qpf * type_shares(shares)[1]                          # snow liquid equivalent: any ratio can be applied later
    by_type = type_qpf_step_mm(qpf, shares)                               # (4, M, ny, nx) mm liquid-equivalent
    cum_x["zr"] += by_type[2]
    cum_x["pl"] += by_type[3]
    for k in [k for k in acc.cum if k < fhr - 24]:
        del acc.cum[k]
    mask = encode.mask_of(idx)
    absent = [i for i in range(M) if i not in idx]
    cum_q = acc.cum[fhr].copy()
    cz, cp, cs, cw = cum_x["zr"].copy(), cum_x["pl"].copy(), cum_x["snow"].copy(), cum_x["swe"].copy()
    for arr in (cum_q, cz, cp, cs, cw):
        arr[absent] = 0.0                                                 # keep dropped members' planes empty
    n = encode.write_cum(run_dir / "cum" / f"f{fhr:03d}.bin.gz", fhr, fhr - prev_fhr, cum_q, cs, domain, cz, cp, mask, swe_cum_mm=cw)
    n += encode.write_ptype(run_dir / "ptype" / f"f{fhr:03d}.bin.gz", fhr, fhr - prev_fhr, probs, domain, mask)
    if temps is not None:                                                 # daily maximum / minimum 2 m temperature
        got = {i: (res[m]["tmax"], res[m]["tmin"], tuple(int(v) for v in res[m]["twin"]))
               for i, m in zip(idx, present) if res[m].get("tmax") is not None and res[m].get("twin") is not None}
        for label in temps.step(fhr, got):
            n += temps.write_day(run_dir, label, domain)
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
    ap.add_argument("--family", default="gefs", choices=sorted(families.FAMILIES), help="which ensemble to process")
    ap.add_argument("--manifest", default="latest.json", help="manifest file name under --out")
    ap.add_argument("--synoptic-too", action="store_true",
                    help="GEFS only: make sure a 00Z/12Z cycle is also published (synoptic.json), because IFS and GEPS exist only for those")
    ap.add_argument("--members", type=int, default=None, help="use only the first N members (testing)")
    ap.add_argument("--max-fhr", type=int, default=C.DEFAULT_MAX_FHR)
    ap.add_argument("--workers", type=int, default=int(os.environ.get("GEFS_WORKERS", 6)))
    ap.add_argument("--lookahead", type=int, default=2, help="forecast steps downloaded ahead of the writer")
    ap.add_argument("--wait-minutes", type=float, default=90.0, help="hard cap on waiting for not-yet-posted files at one forecast hour")
    ap.add_argument("--min-members", type=int, default=None, help="keep publishing while at least this many members remain (default: per family)")
    ap.add_argument("--straggler-minutes", type=float, default=15.0,
                    help="once --min-members have arrived, drop members that are still missing this long after the last arrival")
    ap.add_argument("--poll-seconds", type=float, default=45.0)
    ap.add_argument("--origin-after-seconds", type=float, default=120.0,
                    help="ask NOMADS for a file only after AWS/GCS have lacked it this long (it is simply not posted yet)")
    ap.add_argument("--sources", default=",".join(C.DEFAULT_SOURCE_ORDER))
    ap.add_argument("--base", action="append", help="override a source base URL: name=url (testing)")
    ap.add_argument("--reuse-from", default=os.environ.get("GEFS_REUSE_FROM"), metavar="URL",
                    help="base URL of the published data folder (e.g. https://user.github.io/repo/data). Hours of this same "
                         "cycle that it already holds are copied from there instead of being downloaded from NOAA again")
    ap.add_argument("--force", action="store_true", help="ignore published data and download everything again")
    ap.add_argument("--slr-mode", default=C.SLR_MODE, help='"kuchera" or "fixed:12.5"')
    ap.add_argument("-v", "--verbose", action="store_true")
    args = ap.parse_args(argv)
    logging.basicConfig(level=logging.DEBUG if args.verbose else logging.INFO,
                        format="%(asctime)s %(levelname)s %(message)s", datefmt="%H:%M:%S")
    rc = run(args)
    if args.synoptic_too and args.family == "gefs" and rc == 0 and getattr(args, "resolved_cycle", None) is not None:
        rc = publish_synoptic(args)
    return rc


def synoptic_target(main: dt.datetime, now: dt.datetime | None = None, bases: dict | None = None) -> dt.datetime:
    """The 00/12Z cycle GEFS should also publish for the blend: the newest one that IFS and GEPS both already have, so
    the three models line up; the 00/12Z cycle at or before the main GEFS cycle when they cannot be reached."""
    fallback = floor_cycle(main, 12)
    now = now or dt.datetime.now(dt.timezone.utc)
    try:
        domain = C.Domain()
        chains = [build_chain(families.get(k), build_gates(C.HOST_POLICIES, 1), (), bases, domain) for k in ("ifs", "geps")]
        cyc = floor_cycle(now, 12)
        for _ in range(4):
            if cyc <= fallback and all(ch.exists(cyc.strftime("%Y%m%d"), cyc.strftime("%H"), 3) for ch in chains):
                return cyc
            cyc -= dt.timedelta(hours=12)
    except Exception as e:  # noqa: BLE001 - the blend is optional; never let the probe stop the GEFS build
        log.warning("could not check which IFS/GEPS cycles exist (%s); using %s", e, fallback.strftime("%Y%m%d%H"))
    return fallback


def publish_synoptic(args, target_fn=None) -> int:
    """GEFS also runs at 06/18Z but IFS and GEPS do not: keep a 00/12Z cycle as synoptic.json for the blend."""
    main = args.resolved_cycle
    out = Path(args.out)
    target = (target_fn or synoptic_target)(main)
    if target == main:
        shutil.copyfile(out / "latest.json", out / "synoptic.json")
        return 0
    args2 = argparse.Namespace(**vars(args))
    args2.cycle = target.strftime("%Y%m%d%H")
    args2.manifest = "synoptic.json"
    log.info("GEFS %s is not the cycle IFS/GEPS have; also building %s for the multi-model blend", main.strftime("%Y%m%d%H"), args2.cycle)
    return run(args2)


if __name__ == "__main__":
    sys.exit(main())
