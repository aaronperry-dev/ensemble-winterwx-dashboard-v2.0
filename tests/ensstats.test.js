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
