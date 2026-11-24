"""HTTP client with retries, exponential backoff, Retry-After handling and range support."""
from __future__ import annotations

import logging
import random
import time

import requests

from .config import USER_AGENT
from .ratelimit import HostGate

log = logging.getLogger(__name__)


class NotAvailable(Exception):
    """The object does not exist (yet) on this source (HTTP 404, or 403 on object stores)."""


class SourceError(Exception):
    """Transient or structural failure talking to a source (after retries)."""


RETRY_STATUS = {429, 500, 502, 503, 504}


class HttpClient:
    def __init__(self, session: requests.Session | None = None, timeout=(10, 90),
                 max_retries: int = 4, backoff_base: float = 0.6, backoff_cap: float = 30.0):
        self.session = session or requests.Session()
        self.session.headers.setdefault("User-Agent", USER_AGENT)
        self.timeout = timeout
        self.max_retries = max_retries
        self.backoff_base = backoff_base
        self.backoff_cap = backoff_cap

    # ------------------------------------------------------------------ core
    def request(self, url: str, gate: HostGate, headers: dict | None = None,
                missing_statuses: tuple[int, ...] = (404,)) -> requests.Response:
        last_err: Exception | None = None
        for attempt in range(self.max_retries + 1):
            if gate.exhausted_404_budget():
                raise SourceError(f"{gate.name}: 404 budget exhausted, backing off")
            try:
                with gate:
                    resp = self.session.get(url, headers=headers, timeout=self.timeout)
                    body = resp.content  # read inside the gate so concurrency cap covers transfer
            except (requests.ConnectionError, requests.Timeout, requests.exceptions.ChunkedEncodingError) as e:
                last_err = e
                self._sleep(attempt, None, gate)
                continue
            status = resp.status_code
            if status in (200, 206):
                gate.bytes += len(body)
                return resp
            if status in missing_statuses:
                gate.note_404()
                raise NotAvailable(f"{status} {url}")
            if status in RETRY_STATUS:
                last_err = SourceError(f"HTTP {status} {url}")
                self._sleep(attempt, resp.headers.get("Retry-After"), gate, penalize=True)
                continue
            raise SourceError(f"HTTP {status} {url}")
        raise SourceError(f"giving up after {self.max_retries + 1} attempts: {url}: {last_err}")

    def _sleep(self, attempt: int, retry_after: str | None, gate: HostGate, penalize: bool = False):
        delay = min(self.backoff_cap, self.backoff_base * (2 ** attempt)) * (0.75 + 0.5 * random.random())
        if retry_after:
            try:
                delay = max(delay, float(retry_after))
            except ValueError:
                pass
        if penalize:
            # A 429/503 means the server wants fewer hits: pause the whole host gate, not just this thread.
            gate.limiter.penalize(delay)
        else:
            time.sleep(delay)

    # ------------------------------------------------------------------ helpers
    def get_text(self, url: str, gate: HostGate, **kw) -> str:
        return self.request(url, gate, **kw).text

    def get_bytes(self, url: str, gate: HostGate, **kw) -> bytes:
        return self.request(url, gate, **kw).content

    def get_range(self, url: str, gate: HostGate, start: int, end: int | None, **kw) -> bytes:
        """Fetch bytes [start, end] inclusive (end=None -> to EOF). Verifies length."""
        rng = f"bytes={start}-{'' if end is None else end}"
        resp = self.request(url, gate, headers={"Range": rng}, **kw)
        data = resp.content
        if resp.status_code == 200:
            # Server ignored Range and sent the whole object. Slice locally (correct but wasteful).
            log.warning("%s ignored Range header; slicing %d bytes locally", gate.name, len(data))
            data = data[start:(None if end is None else end + 1)]
        elif end is not None and len(data) != end - start + 1:
            raise SourceError(f"short read from {url}: wanted {end - start + 1}, got {len(data)}")
        return data
