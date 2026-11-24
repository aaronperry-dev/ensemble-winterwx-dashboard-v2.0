/* Ensemble statistics for the ensemble viewer (pure functions; runs in the browser and in Node tests).
 *
 * Layout of `data`: Float32Array [member][y][x], length M*ny*nx.
 *
 *   min, max, mean, p<N>   per-cell statistics (percentiles use linear interpolation between order statistics)
 *   pmm                    probability-matched mean (Ebert 2001): the ensemble-mean field supplies the spatial
 *                          pattern, the pooled member values (all members, whole domain) supply the amplitudes.
 *   lpmm                   local probability-matched mean (Clark 2017; Snook et al.): the same matching done
 *                          independently in a square neighborhood around every grid cell.
 *
 * Weights (opts.weights, one per member, any scale): used when the members come from several models and each model
 * should count the same however many members it has.  Mean and probabilities are weighted averages; percentiles are
 * the weighted analogue of linear interpolation between order statistics (below); PMM and LPMM match amplitudes by
 * cumulative weight instead of by count.  Equal weights take the unweighted code path, so single-model results do not change.
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

  // ------------------------------------------------------------------ weighted versions
  /** Weights scaled to sum to 1, or null when they are all equal (the unweighted code is then exact and faster). */
  function normWeights(w, M) {
    if (!w || w.length < M) return null;
    let sum = 0, lo = Infinity, hi = -Infinity;
    for (let m = 0; m < M; m++) { const v = w[m]; if (!(v >= 0)) return null; sum += v; if (v < lo) lo = v; if (v > hi) hi = v; }
    if (!(sum > 0) || hi - lo <= 1e-9 * hi) return null;
    const out = new Float64Array(M);
    for (let m = 0; m < M; m++) out[m] = w[m] / sum;
    return out;
  }

  /** Runs of consecutive members that share a weight: [{s: first member, c: count, w}]. */
  function weightRuns(w, M) {
    const runs = [];
    for (let m = 0; m < M; m++) {
      const r = runs[runs.length - 1];
      if (r && Math.abs(r.w - w[m]) <= 1e-12) r.c++; else runs.push({ s: m, c: 1, w: w[m] });
    }
    return runs;
  }

  /**
   * Weighted percentile, p in [0,1], of n values with weights summing to 1 (both ascending-sorted together, in place).
   * Each value sits at the cumulative weight of its centre, rescaled so the smallest is at 0 and the largest at 1 -- for
   * equal weights that is exactly the usual (i-1)/(n-1) interpolation.
   */
  function wPercentileSorted(v, w, n, p) {
    if (n === 1) return v[0];
    const first = 0.5 * w[0];
    let C = 0, c0 = first;
    // centre of the last value: 1 - w[n-1]/2
    const span = (1 - 0.5 * w[n - 1]) - first;
    if (!(span > 0)) return v[0];
    let prevQ = 0, prevV = v[0];
    for (let i = 0; i < n; i++) {
      const q = (C + 0.5 * w[i] - c0) / span;
      if (q >= p) return i === 0 || q === prevQ ? v[i] : prevV + (v[i] - prevV) * ((p - prevQ) / (q - prevQ));
      prevQ = q; prevV = v[i]; C += w[i];
    }
    return v[n - 1];
  }

  /** Weighted percentile of arbitrary (unsorted) values; p in [0,1].  `weights` need not be normalised. */
  function wPercentile(values, weights, p) {
    const n = values.length, idx = Array.from({ length: n }, (_, i) => i).sort((a, b) => values[a] - values[b]);
    let sum = 0; for (let i = 0; i < n; i++) sum += weights[i];
    const v = new Float64Array(n), w = new Float64Array(n);
    idx.forEach((k, j) => { v[j] = values[k]; w[j] = weights[k] / sum; });
    return wPercentileSorted(v, w, n, p);
  }

  function meanFieldW(data, M, N, w) {
    const out = new Float32Array(N);
    for (let m = 0; m < M; m++) { const o = m * N, wm = w[m]; for (let i = 0; i < N; i++) out[i] += wm * data[o + i]; }
    return out;
  }

  function perCellW(data, M, N, stat, w) {
    if (stat === 'mean') return meanFieldW(data, M, N, w);
    const p = isPct(stat) ? pctOf(stat) : null;
    if (p === null && stat !== 'min' && stat !== 'max') throw new Error(`unknown statistic ${stat}`);
    const out = new Float32Array(N), v = new Float64Array(M), ws = new Float64Array(M);
    for (let i = 0; i < N; i++) {
      let mx = -Infinity, mn = Infinity;
      for (let m = 0; m < M; m++) { const x = data[m * N + i]; if (x > mx) mx = x; if (x < mn) mn = x; }
      if (stat === 'max') out[i] = mx;
      else if (stat === 'min') out[i] = mn;
      else if (mx === mn) out[i] = mx;
      else {
        let n = 0;                                            // insertion sort of (value, weight) pairs
        for (let m = 0; m < M; m++) {
          const x = data[m * N + i], wm = w[m]; let j = n++;
          while (j > 0 && v[j - 1] > x) { v[j] = v[j - 1]; ws[j] = ws[j - 1]; j--; }
          v[j] = x; ws[j] = wm;
        }
        out[i] = wPercentileSorted(v, ws, M, p);
      }
    }
    return out;
  }

  /** Probability-matched mean where each pooled value carries its member's weight (whole-domain mass = number of cells). */
  function pmmW(data, M, N, w) {
    const mean = meanFieldW(data, M, N, w), runs = weightRuns(w, M);
    const sorted = runs.map((r) => Float32Array.from(data.subarray(r.s * N, (r.s + r.c) * N)).sort());
    const order = new Uint32Array(N);
    for (let i = 0; i < N; i++) order[i] = i;
    order.sort((a, b) => mean[a] - mean[b]);
    const out = new Float32Array(N), ptr = new Array(runs.length).fill(0);
    let cum = 0, j = 0, last = 0;
    while (j < N) {
      let bg = -1, bv = Infinity;
      for (let g = 0; g < runs.length; g++) if (ptr[g] < sorted[g].length && sorted[g][ptr[g]] < bv) { bv = sorted[g][ptr[g]]; bg = g; }
      if (bg < 0) break;
      ptr[bg]++; cum += runs[bg].w; last = bv;
      while (j < N && j + 0.5 < cum) out[order[j++]] = bv;
    }
    while (j < N) out[order[j++]] = last;
    return out;
  }

  function lpmmW(data, M, ny, nx, R, w) {
    const N = ny * nx, mean = meanFieldW(data, M, N, w), out = new Float32Array(N), runs = weightRuns(w, M);
    const cells = (2 * R + 1) * (2 * R + 1), buf = new Float32Array(M * cells), ptr = new Array(runs.length), seg = new Array(runs.length);
    for (let y = 0; y < ny; y++) {
      const y0 = Math.max(0, y - R), y1 = Math.min(ny - 1, y + R);
      for (let x = 0; x < nx; x++) {
        const i = y * nx + x, m0 = mean[i];
        if (!(m0 > 0)) continue;
        const x0 = Math.max(0, x - R), x1 = Math.min(nx - 1, x + R);
        let less = 0, eq = 0;
        for (let yy = y0; yy <= y1; yy++) for (let xx = x0; xx <= x1; xx++) { const mj = mean[yy * nx + xx]; if (mj < m0) less++; else if (mj === m0) eq++; }
        let n = 0;
        runs.forEach((r, g) => {
          const a = n;
          for (let m = r.s; m < r.s + r.c; m++) for (let yy = y0; yy <= y1; yy++) { const row = m * N + yy * nx; for (let xx = x0; xx <= x1; xx++) buf[n++] = data[row + xx]; }
          seg[g] = buf.subarray(a, n).sort(); ptr[g] = 0;
        });
        const target = less + ((eq - 1) >> 1) + 0.5;
        let cum = 0, val = 0;
        for (;;) {
          let bg = -1, bv = Infinity;
          for (let g = 0; g < runs.length; g++) if (ptr[g] < seg[g].length && seg[g][ptr[g]] < bv) { bv = seg[g][ptr[g]]; bg = g; }
          if (bg < 0) break;
          ptr[bg]++; cum += runs[bg].w; val = bv;
          if (target < cum) break;
        }
        out[i] = val;
      }
    }
    return out;
  }

  /** stat: min | max | mean | p<N> | pmm | lpmm. Returns Float32Array(ny*nx) in the units of `data`. */
  function compute(data, M, ny, nx, stat, opts) {
    const N = ny * nx, w = normWeights(opts && opts.weights, M);
    if (stat === 'pmm') return w ? pmmW(data, M, N, w) : pmm(data, M, N);
    if (stat === 'lpmm') return w ? lpmmW(data, M, ny, nx, (opts && opts.radius) || 4, w) : lpmm(data, M, ny, nx, (opts && opts.radius) || 4);
    return w ? perCellW(data, M, N, stat, w) : perCell(data, M, N, stat);
  }

  /** Indices of the members flagged present in a file's bitmask (0 = legacy file: everyone). */
  function membersFromMask(mask, nmem) {
    const out = [];
    for (let k = 0; k < nmem; k++) if (!mask || Math.floor(mask / 2 ** k) % 2 === 1) out.push(k);
    return out;
  }

  return { PCTS, TIP_ORDER, membersFromMask, isPct, pctOf, quickselect, percentile, insertionSort, meanField, perCell, pmm, lpmm, compute,
    normWeights, weightRuns, wPercentile, wPercentileSorted, meanFieldW, perCellW, pmmW, lpmmW };
});
