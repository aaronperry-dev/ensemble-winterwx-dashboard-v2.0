"""Parse wgrib2-style .idx inventories and plan byte-range requests."""
from __future__ import annotations

import re
from dataclasses import dataclass


@dataclass(frozen=True)
class IdxEntry:
    num: int
    start: int
    end: int | None  # inclusive; None = to end of file (last message)
    var: str
    level: str
    ftype: str       # e.g. "3 hour fcst", "0-3 hour acc fcst"


def parse_idx(text: str) -> list[IdxEntry]:
    rows = []
    for line in text.strip().splitlines():
        parts = line.split(":")
        if len(parts) < 6:
            continue
        rows.append((int(parts[0]), int(parts[1]), parts[3], parts[4], parts[5]))
    out = []
    for i, (num, start, var, level, ftype) in enumerate(rows):
        end = rows[i + 1][1] - 1 if i + 1 < len(rows) else None
        out.append(IdxEntry(num, start, end, var, level, ftype))
    return out


def find_entry(entries: list[IdxEntry], var: str, level: str) -> IdxEntry | None:
    for e in entries:
        if e.var == var and e.level == level:
            return e
    return None


def merge_ranges(entries: list[IdxEntry]) -> list[tuple[int, int | None]]:
    """Merge byte-adjacent entries into the fewest contiguous ranges.

    S3 and most CDNs do not support multi-range requests, so only exactly-adjacent messages are
    merged (zero wasted bytes).
    """
    if not entries:
        return []
    es = sorted(entries, key=lambda e: e.start)
    ranges: list[list] = [[es[0].start, es[0].end]]
    for e in es[1:]:
        last = ranges[-1]
        if last[1] is not None and e.start == last[1] + 1:
            last[1] = e.end
        else:
            ranges.append([e.start, e.end])
    return [(a, b) for a, b in ranges]


_ACC_RE = re.compile(r"(\d+)-(\d+) hour acc")


def parse_acc_window(ftype: str) -> tuple[int, int] | None:
    m = _ACC_RE.search(ftype)
    return (int(m.group(1)), int(m.group(2))) if m else None


def split_grib_messages(blob: bytes) -> list[bytes]:
    """Split concatenated GRIB2 messages using the 64-bit length in section 0."""
    msgs, pos, n = [], 0, len(blob)
    while pos < n:
        if blob[pos:pos + 4] != b"GRIB":
            nxt = blob.find(b"GRIB", pos + 1)
            if nxt < 0:
                break
            pos = nxt
            continue
        if pos + 16 > n:
            raise ValueError("truncated GRIB header")
        total = int.from_bytes(blob[pos + 8:pos + 16], "big")
        end = pos + total
        if end > n or blob[end - 4:end] != b"7777":
            raise ValueError("truncated or corrupt GRIB message")
        msgs.append(blob[pos:end])
        pos = end
    return msgs
