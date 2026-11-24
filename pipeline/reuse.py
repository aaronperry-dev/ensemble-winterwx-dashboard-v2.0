"""Reuse GEFS data that is already published for the same cycle instead of downloading it from NOAA again.

Each deploy builds a fresh site folder, so the only copy of a finished cycle is the one on the live GitHub Pages site.
`load_prior` reads that site's `latest.json`; if it is the same cycle, made by the same calculation code, for the same
members / domain / snow-ratio mode, it downloads the published hours (from our own site, not NOAA), checks them, and
returns the longest valid run of consecutive hours.  `seed` then restores the running accumulators from those files so
the pipeline carries on at the next hour exactly as if it had never stopped.
"""
from __future__ import annotations

import concurrent.futures as cf
import hashlib
import logging
import time
from pathlib import Path

import numpy as np
import requests

from . import config as C
from . import encode

log = logging.getLogger("gefs.reuse")
ROOT = Path(__file__).resolve().parent
# Files whose contents change the numbers we publish.  Transport code (sources, http, ratelimit) and site assembly do not.
REV_FILES = ("bourgouin.py", "config.py", "encode.py", "grib.py", "member.py", "products.py", "run.py")


def code_rev() -> str:
    """Short fingerprint of the calculation code: a published run is reused only if it matches."""
    h = hashlib.sha1()
    for name in REV_FILES:
        h.update(name.encode())
        h.update((ROOT / name).read_bytes().replace(b"\r\n", b"\n"))
    return h.hexdigest()[:12]


def _get(session: requests.Session, url: str, tries: int = 3, timeout=(10, 60)) -> bytes | None:
    for k in range(tries):
        try:
            r = session.get(url, timeout=timeout, headers={"Cache-Control": "no-cache"})
            if r.status_code == 200:
                return r.content
            if r.status_code == 404:
                return None
        except requests.RequestException:
            pass
        time.sleep(0.5 * (k + 1))
    return None


def _check_file(path: Path, kind: int, fhr: int, nmem: int) -> bool:
    try:
        d = encode.read_bin(path)
    except Exception:  # noqa: BLE001 - truncated gzip, bad magic, wrong size: all mean "do not trust it"
        return False
    return d["kind"] == kind and d["fhr"] == fhr and d["nmem"] == nmem


def load_prior(base_url: str, run_dir: Path, cyc, members: list[str], domain: C.Domain, slr_mode: str,
               fhrs: list[int], session: requests.Session | None = None) -> dict | None:
    """Download and validate the published hours of this cycle.  Returns None when nothing can be reused."""
    import json
    session = session or requests.Session()
    base = base_url.rstrip("/")
    raw = _get(session, f"{base}/latest.json?cb={int(time.time())}")
    if raw is None:
        log.info("reuse: no published data at %s", base)
        return None
    try:
        man = json.loads(raw)
    except ValueError:
        log.info("reuse: published latest.json is not valid JSON")
        return None
    want = {"init": cyc.strftime("%Y-%m-%dT%H:%M:%SZ"), "members": members, "domain": domain.to_dict(),
            "slr_mode": slr_mode, "code_rev": code_rev(), "version": 1}
    for k, v in want.items():
        if man.get(k) != v:
            why = "different cycle" if k == "init" else f"{k} differs" + (" (calculation code changed)" if k == "code_rev" else "")
            log.info("reuse: published data not reused: %s", why)
            return None
    scales = {"qpf_mm_per_count": C.QPF_MM_PER_COUNT, "snow_cm_per_count": C.SNOW_CM_PER_COUNT,
              "zr_mm_per_count": C.ZR_MM_PER_COUNT, "pl_mm_per_count": C.PL_MM_PER_COUNT}
    if man.get("scales") != scales:
        log.info("reuse: published data not reused: file scales differ")
        return None

    steps = []
    for want_f, st in zip(fhrs, man.get("steps", [])):          # longest prefix that matches the hours we want
        if st.get("fhr") != want_f:
            break
        steps.append(st)
    run_id = man["run"]

    def pull(rel: str) -> bool:
        dest = run_dir / rel
        if not dest.exists() or dest.stat().st_size < 64:
            blob = _get(session, f"{base}/{run_id}/{rel}")
            if blob is None:
                return False
            dest.parent.mkdir(parents=True, exist_ok=True)
            dest.write_bytes(blob)
        return True

    def fetch_step(st: dict) -> bool:
        nmem = len(members)
        return (pull(st["cum"]) and _check_file(run_dir / st["cum"], encode.KIND_CUM, st["fhr"], nmem)
                and pull(st["ptype"]) and _check_file(run_dir / st["ptype"], encode.KIND_PTYPE, st["fhr"], nmem))

    with cf.ThreadPoolExecutor(8) as ex:
        ok = list(ex.map(fetch_step, steps))
    good = 0
    while good < len(ok) and ok[good]:
        good += 1
    if good < len(steps):
        log.warning("reuse: hour f%03d could not be fetched or failed validation; reusing only the %d hours before it",
                    steps[good]["fhr"], good)
    steps = steps[:good]
    if not steps:
        return None
    return {"steps": steps, "dropped": man.get("dropped", {}), "sources_used": man.get("sources_used", {}),
            "complete": len(steps) == len(fhrs)}


def seed(prior: dict, run_dir: Path, acc, cum_x: dict, keep_hours: int = 24) -> None:
    """Restore the running accumulators from published files so the next hour can be added to them."""
    last = prior["steps"][-1]["fhr"]
    for st in prior["steps"]:
        if st["fhr"] >= last - keep_hours:                        # same horizon finalize_step keeps in memory
            acc.cum[st["fhr"]] = encode.read_bin(run_dir / st["cum"])["qpf_mm"].copy()
    d = encode.read_bin(run_dir / prior["steps"][-1]["cum"])
    cum_x["snow"][:] = d["snow_cm"]
    cum_x["zr"][:] = d["zr_mm"]
    cum_x["pl"][:] = d["pl_mm"]
