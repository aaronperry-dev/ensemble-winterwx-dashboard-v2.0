"""Thread-safe token-bucket rate limiting, per-host gates, and a 404 budget."""
from __future__ import annotations

import threading
import time
from collections import deque


class RateLimiter:
    """Token bucket. `acquire()` blocks until a token is available.

    `penalize(seconds)` pushes the whole bucket into a cool-down (used for 429/503 + Retry-After).
    """

    def __init__(self, rate_per_s: float, burst: int = 1):
        if rate_per_s <= 0:
            raise ValueError("rate must be positive")
        self.rate = float(rate_per_s)
        self.capacity = max(1.0, float(burst))
        self._tokens = self.capacity
        self._stamp = time.monotonic()
        self._blocked_until = 0.0
        self._lock = threading.Lock()

    def acquire(self) -> float:
        """Block until one token is available; return seconds waited."""
        waited = 0.0
        while True:
            with self._lock:
                now = time.monotonic()
                if now < self._blocked_until:
                    sleep = self._blocked_until - now
                else:
                    self._tokens = min(self.capacity, self._tokens + (now - self._stamp) * self.rate)
                    self._stamp = now
                    if self._tokens >= 1.0:
                        self._tokens -= 1.0
                        return waited
                    sleep = (1.0 - self._tokens) / self.rate
            time.sleep(sleep)
            waited += sleep

    def penalize(self, seconds: float) -> None:
        with self._lock:
            self._blocked_until = max(self._blocked_until, time.monotonic() + max(0.0, seconds))
            self._tokens = 0.0
            self._stamp = self._blocked_until


class HostGate:
    """Rate limit + concurrency cap + optional 404-per-minute budget for one logical host."""

    def __init__(self, name: str, rate_per_s: float, burst: int, concurrency: int,
                 max_404_per_min: int | None = None):
        self.name = name
        self.limiter = RateLimiter(rate_per_s, burst)
        self.sem = threading.BoundedSemaphore(max(1, concurrency))
        self.max_404_per_min = max_404_per_min
        self._404s: deque[float] = deque()
        self._lock = threading.Lock()
        self.requests = 0
        self.bytes = 0

    # 404 budget (NCEP firewalls IPs that generate too many "missing file" hits)
    def note_404(self) -> None:
        with self._lock:
            self._404s.append(time.monotonic())

    def exhausted_404_budget(self) -> bool:
        if not self.max_404_per_min:
            return False
        with self._lock:
            cutoff = time.monotonic() - 60.0
            while self._404s and self._404s[0] < cutoff:
                self._404s.popleft()
            return len(self._404s) >= self.max_404_per_min

    def __enter__(self):
        self.sem.acquire()
        try:
            self.limiter.acquire()
        except BaseException:
            self.sem.release()
            raise
        with self._lock:
            self.requests += 1
        return self

    def __exit__(self, *exc):
        self.sem.release()
        return False


def build_gates(policies: dict[str, dict], n_workers: int = 1) -> dict[str, HostGate]:
    """Create one gate per policy, dividing the job-wide budget across `n_workers` processes."""
    n = max(1, n_workers)
    gates = {}
    for name, p in policies.items():
        rate = p["rate"] / n
        gates[name] = HostGate(
            name,
            rate_per_s=rate,
            burst=max(1, int(round(p.get("burst", 1) / n))),
            concurrency=max(1, int(round(p.get("concurrency", 1) / n))),
            max_404_per_min=(max(1, int(round(p["max_404_per_min"] / n)))
                             if p.get("max_404_per_min") else None),
        )
    return gates
