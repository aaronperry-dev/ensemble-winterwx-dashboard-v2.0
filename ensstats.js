/* Ensemble statistics for the GEFS viewer (pure functions; runs in the browser and in Node tests).
 *
 * Layout of `data`: Float32Array [member][y][x], length M*ny*nx.
 *
 *   min, max, mean, p<N>   per-cell statistics (percentiles use linear interpolation between order statistics)
 *   pmm                    probability-matched mean (Ebert 2001): the ensemble-mean field supplies the spatial
 *                          pattern, the pooled member values (all members, whole domain) supply the amplitudes.
 *   lpmm                   local probability-matched mean (Clark 2017; Snook et al.): the same matching done
 *                          independently in a square neighborhood around every grid cell.
 */
(function (root, factory) {
  const api = factory();
  if (typeof module === 'object' && module.exports) module.exports = api; else root.EnsStats = api;
})(typeof self !== 'undefined' ? self : this, function () {
  'use strict';

  const PCTS = [5, 10, 25, 50, 75, 90, 95];
  const TIP_ORDER = ['min', 'p5', 'p10', 'p25', 'p50', 'mean', 'p75', 'p90', 'p95', 'max'];

  const isPct = (s) => /^p\d+(\.\d+)?$/.test(s);
  const pctOf = (s) => parseFloat(s.slice(1)) / 100;

  /** In-place Floyd-Rivest selection: afterwards a[k] holds the k-th smallest of a[left..right]. */
  function quickselect(a, k, left, right) {
    while (right > left) {
      if (right - left > 600) {
        const n = right - left + 1, m = k - left + 1, z = Math.log(n), s = 0.5 * Math.exp((2 * z) / 3);
        const sd = 0.5 * Math.sqrt((z * s * (n - s)) / n) * (m - n / 2 < 0 ? -1 : 1);
        quickselect(a, k, Math.max(left, Math.floor(k - (m * s) / n + sd)), Math.min(right, Math.floor(k + ((n - m) * s) / n + sd)));
      }
      const t = a[k];
      let i = left, j = right, tmp;
      tmp = a[left]; a[left] = a[k]; a[k] = tmp;
      if (a[right] > t) { tmp = a[left]; a[left] = a[right]; a[right] = tmp; }
      while (i < j) {
        tmp = a[i]; a[i] = a[j]; a[j] = tmp; i++; j--;
        while (a[i] < t) i++;
        while (a[j] > t) j--;
      }
      if (a[left] === t) { tmp = a[left]; a[left] = a[j]; a[j] = tmp; }
      else { j++; tmp = a[j]; a[j] = a[right]; a[right] = tmp; }
      if (j <= k) left = j + 1;
      if (k <= j) right = j - 1;
    }
  }

  function insertionSort(a, n) {
    for (let i = 1; i < n; i++) { const x = a[i]; let j = i - 1; while (j >= 0 && a[j] > x) { a[j + 1] = a[j]; j--; } a[j + 1] = x; }
  }

  /** p in [0,1] on an ascending array of n values. */
  function percentile(sorted, n, p) {
    const f = p * (n - 1), lo = Math.floor(f), hi = Math.ceil(f);
    return sorted[lo] + (sorted[hi] - sorted[lo]) * (f - lo);
  }

  function meanField(data, M, N) {
    const out = new Float32Array(N);
    for (let m = 0; m < M; m++) { const o = m * N; for (let i = 0; i < N; i++) out[i] += data[o + i]; }
    for (let i = 0; i < N; i++) out[i] /= M;
    return out;
  }

  /** min / max / mean / any percentile, one value per cell. */
  function perCell(data, M, N, stat) {
    if (stat === 'mean') return meanField(data, M, N);
    const out = new Float32Array(N), sc = new Float32Array(M), p = isPct(stat) ? pctOf(stat) : null;
    if (p === null && stat !== 'min' && stat !== 'max') throw new Error(`unknown statistic ${stat}`);
    for (let i = 0; i < N; i++) {
      let mx = -Infinity, mn = Infinity;
      for (let m = 0; m < M; m++) { const v = data[m * N + i]; sc[m] = v; if (v > mx) mx = v; if (v < mn) mn = v; }
      if (stat === 'max') out[i] = mx;
      else if (stat === 'min') out[i] = mn;
      else if (mx === mn) out[i] = mx;                  // all members agree (commonly 0): skip the sort
      else { insertionSort(sc, M); out[i] = percentile(sc, M, p); }
    }
    return out;
  }

  /** Probability-matched mean over the whole domain. */
  function pmm(data, M, N) {
    const mean = meanField(data, M, N), pooled = Float32Array.from(data).sort();
    const order = new Uint32Array(N);
    for (let i = 0; i < N; i++) order[i] = i;
    order.sort((a, b) => mean[a] - mean[b]);
    const out = new Float32Array(N), off = M >> 1;
    for (let j = 0; j < N; j++) out[order[j]] = pooled[j * M + off];
    return out;
  }

  /**
   * Local probability-matched mean. For each cell the pooled member values inside the
   * (2R+1)x(2R+1) window are ranked; the cell receives the pooled value whose rank matches the rank of
   * its ensemble mean among the window's ensemble means. Cells whose ensemble mean is zero stay zero.
   */
  function lpmm(data, M, ny, nx, R) {
    const N = ny * nx, mean = meanField(data, M, N), out = new Float32Array(N);
    const buf = new Float32Array(M * (2 * R + 1) * (2 * R + 1)), off = M >> 1;
    for (let y = 0; y < ny; y++) {
      const y0 = Math.max(0, y - R), y1 = Math.min(ny - 1, y + R);
      for (let x = 0; x < nx; x++) {
        const i = y * nx + x, m0 = mean[i];
        if (!(m0 > 0)) continue;
        const x0 = Math.max(0, x - R), x1 = Math.min(nx - 1, x + R);
        let less = 0, eq = 0;
        for (let yy = y0; yy <= y1; yy++) for (let xx = x0; xx <= x1; xx++) { const mj = mean[yy * nx + xx]; if (mj < m0) less++; else if (mj === m0) eq++; }
        let n = 0;
        for (let m = 0; m < M; m++) for (let yy = y0; yy <= y1; yy++) { const row = m * N + yy * nx; for (let xx = x0; xx <= x1; xx++) buf[n++] = data[row + xx]; }
        const idx = Math.min(n - 1, (less + ((eq - 1) >> 1)) * M + off);
        quickselect(buf, idx, 0, n - 1);
        out[i] = buf[idx];
      }
    }
    return out;
  }

  /** stat: min | max | mean | p<N> | pmm | lpmm. Returns Float32Array(ny*nx) in the units of `data`. */
  function compute(data, M, ny, nx, stat, opts) {
    const N = ny * nx;
    if (stat === 'pmm') return pmm(data, M, N);
    if (stat === 'lpmm') return lpmm(data, M, ny, nx, (opts && opts.radius) || 4);
    return perCell(data, M, N, stat);
  }

  return { PCTS, TIP_ORDER, isPct, pctOf, quickselect, percentile, insertionSort, meanField, perCell, pmm, lpmm, compute };
});
