"""Command-line entry point.

    python -m pipeline.run --out site/data                 # newest cycle, 31 members, f003..f240
    python -m pipeline.run --out site/data --max-fhr 24    # quick test
    python -m pipeline.run --out site/data --cycle 2026112406 --wait-minutes 0

Writes   <out>/<YYYYMMDD_CC>/cum/fHHH.bin.gz     cumulative QPF + snowfall, every member
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


# ----------------------------------------------------------------------------- main loop
def run(args) -> int:
    domain = C.Domain()
    members = list(C.MEMBERS[: args.members])
    fhrs = C.forecast_hours(args.max_fhr)
    bases = dict(kv.split("=", 1) for kv in (args.base or []))
    order = tuple(args.sources.split(","))

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
    log.info("GEFS %s: %d members x %d steps, sources=%s, workers=%d", run_id, len(members), len(fhrs),
             ",".join(order), args.workers)

    cfg = {"n_workers": args.workers, "order": order, "bases": bases, "domain": domain.__dict__,
           "slr_mode": args.slr_mode}
    ex = cf.ProcessPoolExecutor(max_workers=args.workers, initializer=init_worker, initargs=(cfg,))

    acc = ApcpAccumulator(len(members), domain.ny, domain.nx)
    cum_x = {k: np.zeros((len(members), domain.ny, domain.nx), np.float32) for k in ("snow", "zr", "pl")}
    steps: list[dict] = []
    src_counts: dict[str, int] = {}
    total_bytes = 0

    state = {f: {"res": {}, "err": {}} for f in fhrs}
    futures: dict[cf.Future, tuple[int, str]] = {}
    retry_q: list[tuple[float, int, str]] = []
    submitted = -1
    head = 0
    blocked_since: float | None = None
    deadline_s = args.wait_minutes * 60.0
    t_start = time.time()
    fatal = None

    def submit(fhr: int, m: str):
        futures[ex.submit(run_task, (date, cc, m, fhr))] = (fhr, m)

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
                    submit(fhrs[submitted], m)
            now_t = time.time()
            for item in [r for r in retry_q if r[0] <= now_t]:
                retry_q.remove(item)
                submit(item[1], item[2])

            done, _ = cf.wait(list(futures), timeout=2.0, return_when=cf.FIRST_COMPLETED)
            for fut in done:
                fhr, m = futures.pop(fut)
                r = fut.result()
                st = state[fhr]
                if r["status"] == "ok":
                    st["res"][m] = r
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
                        fatal = f"f{fhr:03d} {m}: {r['msg']}"
                    else:
                        retry_q.append((time.time() + 5.0, fhr, m))
            if fatal:
                break

            fhr = fhrs[head]
            st = state[fhr]
            if len(st["res"]) == len(members):
                total_bytes += finalize_step(fhr, fhrs[head - 1] if head else 0, members, st["res"], acc,
                                             cum_x, run_dir, domain, steps, cyc)
                state[fhr] = {"res": {}, "err": {}}
                write_manifest(False)
                log.info("f%03d done (%d/%d)  %.0fs elapsed", fhr, head + 1, len(fhrs), time.time() - t_start)
                head += 1
                blocked_since = None
                continue
            if blocked_since is not None and time.time() - blocked_since > deadline_s:
                fatal = f"f{fhr:03d} still unavailable after {args.wait_minutes} min"
                break
    finally:
        ex.shutdown(wait=False, cancel_futures=True)

    write_manifest(head >= len(fhrs))
    if fatal:
        log.warning("stopping early: %s", fatal)
    log.info("finished: %d/%d steps, sources=%s, %.1f MB written, %.0fs", len(steps), len(fhrs), src_counts,
             total_bytes / 1e6, time.time() - t_start)
    return 0 if steps else 2


def finalize_step(fhr, prev_fhr, members, res, acc, cum_x, run_dir, domain, steps, cyc) -> int:
    windows = {tuple(res[m]["window"]) for m in members}
    if len(windows) != 1:
        raise SourceError(f"members disagree on APCP window at f{fhr:03d}: {windows}")
    win = windows.pop()
    acc.add(win, np.stack([res[m]["apcp"] for m in members]))
    qpf = acc.step(fhr, prev_fhr)                                         # (M, ny, nx) mm
    probs = np.stack([res[m]["probs"] for m in members], axis=1)          # (4, M, ny, nx)
    slr = np.stack([res[m]["slr"] for m in members])                      # (M, ny, nx)
    wet = qpf >= C.MIN_PRECIP_MM
    probs = (probs * wet[None]).astype(np.uint8)
    snow = snowfall_step_cm(qpf, probs, slr)
    cum_x["snow"] += snow
    by_type = type_qpf_step_mm(qpf, probs)                                # (4, M, ny, nx) mm liquid-equivalent
    cum_x["zr"] += by_type[2]
    cum_x["pl"] += by_type[3]
    for k in [k for k in acc.cum if k < fhr - 24]:
        del acc.cum[k]
    n = encode.write_cum(run_dir / "cum" / f"f{fhr:03d}.bin.gz", fhr, fhr - prev_fhr, acc.cum[fhr], cum_x["snow"], domain,
                         cum_x["zr"], cum_x["pl"])
    n += encode.write_ptype(run_dir / "ptype" / f"f{fhr:03d}.bin.gz", fhr, fhr - prev_fhr, probs, domain)
    steps.append({
        "fhr": fhr, "step_h": fhr - prev_fhr,
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
    ap.add_argument("--wait-minutes", type=float, default=90.0, help="how long to wait for not-yet-posted files")
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
