// Run with:  node --test tests/
const test = require('node:test');
const assert = require('node:assert/strict');
const E = require('../ensstats.js');

function rng(seed) { let s = seed >>> 0; return () => ((s = (Math.imul(s, 1664525) + 1013904223) >>> 0) / 4294967296); }

/** Synthetic precipitation-like ensemble: blobs of random position/amplitude per member, with dry areas. */
function synth(M, ny, nx, seed = 7) {
  const r = rng(seed), N = ny * nx, d = new Float32Array(M * N);
  for (let m = 0; m < M; m++) {
    for (let b = 0; b < 5; b++) {
      const cx = r() * nx, cy = r() * ny, a = Math.pow(r(), 2) * 40, s = 2 + r() * 5;
      for (let y = 0; y < ny; y++) for (let x = 0; x < nx; x++) d[m * N + y * nx + x] += a * Math.exp(-((x - cx) ** 2 + (y - cy) ** 2) / (2 * s * s));
    }
    for (let i = 0; i < N; i++) if (d[m * N + i] < 0.5) d[m * N + i] = 0;       // dry cells are exactly zero
  }
  return d;
}

test('quickselect matches a full sort, including heavy duplicates', () => {
  const r = rng(3);
  for (const dup of [false, true]) {
    for (let trial = 0; trial < 30; trial++) {
      const n = 1 + Math.floor(r() * 3000), a = Float32Array.from({ length: n }, () => (dup ? Math.floor(r() * 4) * 0 + (r() < 0.7 ? 0 : Math.floor(r() * 5)) : r()));
      const k = Math.floor(r() * n), ref = Float32Array.from(a).sort()[k], b = Float32Array.from(a);
      E.quickselect(b, k, 0, n - 1);
      assert.equal(b[k], ref);
    }
  }
});

test('percentiles use linear interpolation on 31 members', () => {
  const M = 31, d = Float32Array.from({ length: M }, (_, i) => i);        // 0..30, one cell
  const get = (s) => E.compute(d, M, 1, 1, s)[0];
  assert.equal(get('min'), 0); assert.equal(get('max'), 30); assert.equal(get('mean'), 15);
  assert.equal(get('p50'), 15); assert.equal(get('p10'), 3); assert.equal(get('p90'), 27);
  assert.ok(Math.abs(get('p5') - 1.5) < 1e-6 && Math.abs(get('p95') - 28.5) < 1e-5);
  assert.ok(Math.abs(get('p25') - 7.5) < 1e-6 && Math.abs(get('p75') - 22.5) < 1e-5);
});

test('statistics are ordered min <= p5 <= ... <= p95 <= max everywhere', () => {
  const M = 31, ny = 20, nx = 30, d = synth(M, ny, nx), f = {};
  for (const s of E.TIP_ORDER) f[s] = E.compute(d, M, ny, nx, s);
  const keys = ['min', 'p5', 'p10', 'p25', 'p50', 'p75', 'p90', 'p95', 'max'];
  for (let i = 0; i < ny * nx; i++) for (let k = 1; k < keys.length; k++) assert.ok(f[keys[k]][i] >= f[keys[k - 1]][i] - 1e-6);
});

test('PMM keeps the mean field\'s ordering and takes its amplitudes from the pooled members', () => {
  const M = 31, ny = 25, nx = 40, N = ny * nx, d = synth(M, ny, nx), mean = E.meanField(d, M, N), pm = E.pmm(d, M, N);
  const idx = Array.from({ length: N }, (_, i) => i).sort((a, b) => mean[a] - mean[b]);
  for (let j = 1; j < N; j++) assert.ok(pm[idx[j]] >= pm[idx[j - 1]], 'monotone in the mean');
  const pooled = Float32Array.from(d).sort();
  const sortedPm = Float32Array.from(pm).sort();
  for (let j = 0; j < N; j++) assert.equal(sortedPm[j], pooled[j * M + (M >> 1)]);   // distribution matched exactly
  for (let i = 0; i < N; i++) if (mean[i] === 0) assert.equal(pm[i], 0);
  const sum = (a) => a.reduce((s, v) => s + v, 0);
  assert.ok(Math.abs(sum(pm) - sum(mean)) / sum(mean) < 0.03, 'domain total preserved to ~3%');
  assert.ok(Math.max(...pm) > Math.max(...mean) * 1.3, 'amplitudes are not smoothed away like the mean');
});

test('LPMM with a window covering the whole domain equals PMM', () => {
  const M = 31, ny = 12, nx = 16, d = synth(M, ny, nx, 11);
  const a = E.lpmm(d, M, ny, nx, 40), b = E.pmm(d, M, ny * nx);
  for (let i = 0; i < ny * nx; i++) assert.ok(Math.abs(a[i] - b[i]) < 1e-5, `cell ${i}: ${a[i]} vs ${b[i]}`);
});

test('LPMM is local, non-negative, zero where the mean is zero, and bounded by the local member max', () => {
  const M = 31, ny = 30, nx = 45, N = ny * nx, R = 3, d = synth(M, ny, nx, 5), mean = E.meanField(d, M, N), lp = E.lpmm(d, M, ny, nx, R);
  for (let y = 0; y < ny; y++) for (let x = 0; x < nx; x++) {
    const i = y * nx + x; assert.ok(lp[i] >= 0);
    if (mean[i] === 0) { assert.equal(lp[i], 0); continue; }
    let mx = 0;
    for (let m = 0; m < M; m++) for (let yy = Math.max(0, y - R); yy <= Math.min(ny - 1, y + R); yy++) for (let xx = Math.max(0, x - R); xx <= Math.min(nx - 1, x + R); xx++) mx = Math.max(mx, d[m * N + yy * nx + xx]);
    assert.ok(lp[i] <= mx + 1e-5);
  }
});

test('LPMM matches a brute-force implementation on a small grid', () => {
  const M = 9, ny = 8, nx = 11, N = ny * nx, R = 2, d = synth(M, ny, nx, 21), mean = E.meanField(d, M, N), lp = E.lpmm(d, M, ny, nx, R);
  for (let y = 0; y < ny; y++) for (let x = 0; x < nx; x++) {
    const i = y * nx + x; if (!(mean[i] > 0)) continue;
    const vals = [], means = [];
    for (let yy = Math.max(0, y - R); yy <= Math.min(ny - 1, y + R); yy++) for (let xx = Math.max(0, x - R); xx <= Math.min(nx - 1, x + R); xx++) {
      means.push(mean[yy * nx + xx]); for (let m = 0; m < M; m++) vals.push(d[m * N + yy * nx + xx]);
    }
    vals.sort((a, b) => a - b);
    const less = means.filter((v) => v < mean[i]).length, eq = means.filter((v) => v === mean[i]).length;
    assert.equal(lp[i], vals[Math.min(vals.length - 1, (less + ((eq - 1) >> 1)) * M + (M >> 1))]);
  }
});

test('LPMM speed on a full GEFS-sized domain is interactive', () => {
  const M = 31, ny = 81, nx = 171, d = synth(M, ny, nx, 2);
  for (const R of [2, 4, 6]) { const t = performance.now(); E.lpmm(d, M, ny, nx, R); const ms = performance.now() - t; console.log(`LPMM R=${R}: ${ms.toFixed(0)} ms`); assert.ok(ms < 4000); }
  const t = performance.now(); E.pmm(d, M, ny * nx); console.log(`PMM: ${(performance.now() - t).toFixed(0)} ms`);
});

test('membersFromMask lists present members (bit k = member k); 0 means a legacy file with everyone', () => {
  assert.deepEqual(E.membersFromMask(0, 4), [0, 1, 2, 3]);
  assert.deepEqual(E.membersFromMask(0b1011, 4), [0, 1, 3]);
  const full = 2 ** 31 - 1, m7 = full - 2 ** 7;
  const got = E.membersFromMask(m7, 31);
  assert.equal(got.length, 30); assert.ok(!got.includes(7) && got.includes(30) && got.includes(0));
  assert.equal(E.membersFromMask(2 ** 31 + 1, 32).join(), '0,31');         // bit 31 must not go negative
});

test('statistics over a compacted member set ignore the dropped member entirely', () => {
  const M = 31, ny = 6, nx = 8, N = ny * nx, d = synth(M, ny, nx, 9), keep = E.membersFromMask(2 ** 31 - 1 - 2 ** 7, 31);
  const poisoned = Float32Array.from(d); for (let i = 0; i < N; i++) poisoned[7 * N + i] = 1e6;       // the dropped member's plane is garbage
  const compact = new Float32Array(keep.length * N); keep.forEach((m, j) => compact.set(poisoned.subarray(m * N, (m + 1) * N), j * N));
  const ref = new Float32Array(keep.length * N); keep.forEach((m, j) => ref.set(d.subarray(m * N, (m + 1) * N), j * N));
  for (const s of ['mean', 'max', 'p50', 'p95', 'pmm', 'lpmm']) {
    const a = E.compute(compact, keep.length, ny, nx, s, { radius: 2 }), b = E.compute(ref, keep.length, ny, nx, s, { radius: 2 });
    for (let i = 0; i < N; i++) assert.equal(a[i], b[i]);
    assert.ok(Math.max(...a) < 1e5);
  }
});

// ---------------------------------------------------------------------------------------------- weighted statistics
function weightsFor(sizes, byModel = true) {            // e.g. [31, 50, 21]: each model counts the same, or every member the same
  const w = []; sizes.forEach((n) => { for (let k = 0; k < n; k++) w.push(byModel ? 1 / (sizes.length * n) : 1 / sizes.reduce((a, b) => a + b, 0)); });
  return Float32Array.from(w);
}

test('equal weights give exactly the unweighted results (single-model views do not change)', () => {
  const M = 31, ny = 8, nx = 10, d = synth(M, ny, nx, 4), w = new Float32Array(M).fill(0.2);
  for (const s of ['mean', 'min', 'max', 'p10', 'p50', 'p95', 'pmm', 'lpmm']) {
    const a = E.compute(d, M, ny, nx, s, { radius: 2 }), b = E.compute(d, M, ny, nx, s, { radius: 2, weights: w });
    for (let i = 0; i < ny * nx; i++) assert.equal(a[i], b[i], s);
  }
  assert.equal(E.normWeights(w, M), null);
  assert.equal(E.normWeights(weightsFor([31, 50, 21], false), 102), null);          // equal-by-member is the unweighted path too
  assert.ok(E.normWeights(weightsFor([31, 50, 21], true), 102));
});

test('weighted percentile equals the unweighted one for equal weights and moves toward heavy members otherwise', () => {
  const v = Array.from({ length: 31 }, (_, i) => i), eq = new Array(31).fill(1);
  for (const p of [0, 0.05, 0.25, 0.5, 0.9, 1]) assert.ok(Math.abs(E.wPercentile(v, eq, p) - E.percentile(Float32Array.from(v), 31, p)) < 1e-9, `p=${p}`);
  assert.equal(E.wPercentile([0, 5, 10], [0.1, 0.1, 0.8], 0), 0);
  assert.equal(E.wPercentile([0, 5, 10], [0.1, 0.1, 0.8], 1), 10);
  const med = E.wPercentile([0, 5, 10], [0.1, 0.1, 0.8], 0.5);
  assert.ok(med > 5 && med < 10, `median ${med} should lean toward the heavy member`);
  assert.ok(Math.abs(E.wPercentile([7, 1, 3], [2, 5, 3], 0.5) - E.wPercentile([1, 3, 7], [5, 3, 2], 0.5)) < 1e-12, 'order of input does not matter');
});

test('weighted mean and probabilities-by-weight match hand calculations; a heavy model pulls the mean', () => {
  const M = 3, d = Float32Array.from([0, 6, 12]), w = Float32Array.from([1, 1, 2]);
  assert.ok(Math.abs(E.compute(d, M, 1, 1, 'mean', { weights: w })[0] - 7.5) < 1e-6);
  assert.equal(E.compute(d, M, 1, 1, 'max', { weights: w })[0], 12);
  assert.equal(E.compute(d, M, 1, 1, 'min', { weights: w })[0], 0);
});

test('weighted percentiles stay ordered and inside the member range', () => {
  const sizes = [31, 50, 21], M = 102, ny = 14, nx = 18, N = ny * nx, d = synth(M, ny, nx, 13), w = weightsFor(sizes);
  const keys = ['min', 'p5', 'p10', 'p25', 'p50', 'p75', 'p90', 'p95', 'max'], f = {};
  for (const k of keys) f[k] = E.compute(d, M, ny, nx, k, { weights: w });
  for (let i = 0; i < N; i++) for (let k = 1; k < keys.length; k++) assert.ok(f[keys[k]][i] >= f[keys[k - 1]][i] - 1e-5, `${keys[k]} at ${i}`);
});

test('weighted PMM: monotone in the weighted mean, amplitudes from the weighted pool, zero where the mean is zero', () => {
  const sizes = [6, 10, 4], M = 20, ny = 16, nx = 22, N = ny * nx, d = synth(M, ny, nx, 17), w = weightsFor(sizes);
  const mean = E.meanFieldW(d, M, N, E.normWeights(w, M)), pm = E.compute(d, M, ny, nx, 'pmm', { weights: w });
  const idx = Array.from({ length: N }, (_, i) => i).sort((a, b) => mean[a] - mean[b]);
  for (let j = 1; j < N; j++) assert.ok(pm[idx[j]] >= pm[idx[j - 1]]);
  for (let i = 0; i < N; i++) if (mean[i] === 0) assert.equal(pm[i], 0);
  // brute force: pooled (value, weight) pairs sorted; the cell of rank j gets the value covering mass j + 0.5
  const wn = E.normWeights(w, M), pairs = []; for (let m = 0; m < M; m++) for (let i = 0; i < N; i++) pairs.push([d[m * N + i], wn[m]]);
  pairs.sort((a, b) => a[0] - b[0]);
  const ref = new Float32Array(N); let cum = 0, j = 0;
  for (const [v, wt] of pairs) { cum += wt; while (j < N && j + 0.5 < cum) ref[idx[j++]] = v; }
  while (j < N) ref[idx[j++]] = pairs[pairs.length - 1][0];
  for (let i = 0; i < N; i++) assert.ok(Math.abs(pm[i] - ref[i]) < 1e-5, `cell ${i}: ${pm[i]} vs ${ref[i]}`);
  assert.ok(Math.max(...pm) > Math.max(...mean) * 1.2);
});

test('weighted LPMM with a window covering the whole domain equals weighted PMM, and matches brute force locally', () => {
  const sizes = [5, 8, 3], M = 16, ny = 9, nx = 12, N = ny * nx, d = synth(M, ny, nx, 23), w = weightsFor(sizes);
  const a = E.compute(d, M, ny, nx, 'lpmm', { radius: 40, weights: w }), b = E.compute(d, M, ny, nx, 'pmm', { weights: w });
  for (let i = 0; i < N; i++) assert.ok(Math.abs(a[i] - b[i]) < 1e-5, `cell ${i}`);
  const R = 2, wn = E.normWeights(w, M), lp = E.compute(d, M, ny, nx, 'lpmm', { radius: R, weights: w }), mean = E.meanFieldW(d, M, N, wn);
  for (let y = 0; y < ny; y++) for (let x = 0; x < nx; x++) {
    const i = y * nx + x; if (!(mean[i] > 0)) continue;
    const pairs = [], means = [];
    for (let yy = Math.max(0, y - R); yy <= Math.min(ny - 1, y + R); yy++) for (let xx = Math.max(0, x - R); xx <= Math.min(nx - 1, x + R); xx++) {
      means.push(mean[yy * nx + xx]); for (let m = 0; m < M; m++) pairs.push([d[m * N + yy * nx + xx], wn[m]]);
    }
    pairs.sort((p, q) => p[0] - q[0]);
    const less = means.filter((v) => v < mean[i]).length, eq = means.filter((v) => v === mean[i]).length, t = less + ((eq - 1) >> 1) + 0.5;
    let cum = 0, val = pairs[pairs.length - 1][0];
    for (const [v, wt] of pairs) { cum += wt; if (t < cum) { val = v; break; } }
    assert.ok(Math.abs(lp[i] - val) < 1e-5, `cell ${i}: ${lp[i]} vs ${val}`);
  }
});

test('weighted statistics on a full 102-member blend are fast enough for a background worker', () => {
  const sizes = [31, 50, 21], M = 102, ny = 81, nx = 171, d = synth(M, ny, nx, 2), w = weightsFor(sizes);
  for (const [s, R, lim] of [['p50', 4, 4000], ['mean', 4, 2000], ['pmm', 4, 12000], ['lpmm', 4, 40000]]) {
    const t = performance.now(); E.compute(d, M, ny, nx, s, { radius: R, weights: w }); const ms = performance.now() - t;
    console.log(`weighted ${s}${s === 'lpmm' ? ` R=${R}` : ''} (102 members): ${ms.toFixed(0)} ms`); assert.ok(ms < lim);
  }
});
