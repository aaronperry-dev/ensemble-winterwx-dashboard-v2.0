/* Ensemble winter-weather viewer: NOAA GEFS, ECMWF IFS ENS, ECCC GEPS, and a blend of them.
 * Loads small gzipped binary files (cumulative QPF/snow per member, p-type probabilities per member),
 * derives any window total / ensemble statistic in the browser, and paints it with a WebGL layer.
 * The blend is made here, not in the pipeline: the members of every model that shares the run time are pooled and
 * (by default) weighted so that each model counts the same however many members it has.
 */
'use strict';

const $ = (s, el = document) => el.querySelector(s);
const isMem = (s) => /^m\d+$/.test(s);                     // m<global member index>
const MM_PER_IN = 25.4, CM_PER_IN = 2.54;
const TIE_ORDER = [2, 3, 1, 0];            // FZRA > PL > SN > RA  (indices into RA, SN, FZRA, PL)
const PERIODS = [['step', 'Step'], [6, '6 hr'], [12, '12 hr'], [24, '24 hr'], [48, '48 hr']];
const PROB_THRESH = { qpf: [0.01, 0.10, 0.25, 0.50, 1.00, 2.00], snow: [1, 3, 6, 12, 18, 24],
  zr: [0.01, 0.05, 0.10, 0.25, 0.50, 1.00], pl: [0.01, 0.05, 0.10, 0.25, 0.50, 1.00] };   // inches
/* Accumulation products: array key in the cumulative file, header scale slot, and SI unit. Snowfall is in cm; the rest are mm of liquid. */
const PRODUCTS = { qpf: { name: 'QPF', key: 'q', scale: 's0', si: 'mm' }, snow: { name: 'Snowfall', key: 's', scale: 's1', si: 'cm' },
  zr: { name: 'Freezing rain', key: 'z', scale: 's2', si: 'mm' }, pl: { name: 'Sleet', key: 'p', scale: 's3', si: 'mm' } };
/* Where a product's numbers come from and how they turn into inches.  Snowfall defaults to the Kuchera-ratio depth in the
 * file; any fixed ratio (10:1, 12:1, 13:1) is applied here to the snow liquid equivalent (array 'w'), so changing the
 * ratio needs no new data. */
function src(kind) {
  const P = PRODUCTS[kind];
  if (kind === 'snow' && S.slr !== 'kuchera' && S.hasSwe) return { key: 'w', scale: 's4', perIn: MM_PER_IN / +S.slr };
  return { key: P.key, scale: P.scale, perIn: P.si === 'cm' ? CM_PER_IN : MM_PER_IN };
}
const siPerIn = (kind) => src(kind).perIn;
const isProb = (s) => s.startsWith('prob:');       // map shows a percentage
const PROB_RAMP = {
  qpf: ['#e5f5e0', '#c7e9c0', '#a1d99b', '#74c476', '#41ab5d', '#31a353', '#238b45', '#006d2c', '#00441b'],
  snow: ['#deebf7', '#c6dbef', '#bdd7e7', '#9ecae1', '#6baed6', '#4292c6', '#3182bd', '#08519c', '#082694'],
};
const PROB_EDGES = [10, 20, 30, 40, 50, 60, 70, 80, 90];

/* The models.  `color` marks a model's members in the tooltip's dot strip. */
const FAMS = {
  gefs: { label: 'GEFS', long: 'NOAA GEFS', color: '#4b5563' },
  ifs: { label: 'IFS ENS', long: 'ECMWF IFS ENS', color: '#0f766e' },
  geps: { label: 'GEPS', long: 'ECCC GEPS', color: '#7e22ce' },
};
const FAM_KEYS = ['gefs', 'ifs', 'geps'];
const sig = () => `${S.fam}|${S.wmode}|${S.src.map((s) => s.key + s.man.init + s.man.generated).join(',')}`;     // identifies the pooled ensemble (cache key)
const S = {
  sets: {}, fam: 'blend', wmode: 'model', note: '',          // loaded manifests; selected model; blend weighting
  src: [], all: [], off: {}, init: '', domain: null, gen: '', complete: true, products: null,
  steps: [], byFhr: new Map(), idx: 0,
  stn: false, stations: null,                                 // station-value layer on/off; New England METAR/ASOS list (stations.json)
  tdays: {}, day: '',                                       // daily MaxT / MinT: per product, the dates every pooled model has; the selected date
  product: 'qpf', period: 24, stat: 'mean', lpmmR: 4, slr: 'kuchera', hasSwe: false, opacity: 0.82,
  cur: null, token: 0, pin: null, playing: false, timer: null,
};
let map, wx, labels, labelData, pinMarker;

/* ============================================================================= data loading */
const CACHE = new Map(), SIZES = new Map();
const BUDGET = (navigator.deviceMemory && navigator.deviceMemory <= 4 ? 90 : 220) * 1048576;   // decoded bytes kept in memory
function cached(url, make) {
  if (CACHE.has(url)) { const p = CACHE.get(url); CACHE.delete(url); CACHE.set(url, p); return p; }
  const p = make().then((o) => {
    SIZES.set(url, (o.q || o.p || o.t).buffer.byteLength);
    let tot = 0; SIZES.forEach((v) => (tot += v));
    for (const k of CACHE.keys()) { if (tot <= BUDGET || CACHE.size <= 8) break; if (k === url) continue; tot -= SIZES.get(k) || 0; CACHE.delete(k); SIZES.delete(k); }
    return o;
  }).catch((e) => { CACHE.delete(url); SIZES.delete(url); throw e; });
  CACHE.set(url, p);
  return p;
}

async function fetchBin(url) {
  return cached(url, async () => {
    const r = await fetch(url);
    if (!r.ok) throw new Error(`${r.status} ${url}`);
    let buf = await r.arrayBuffer();
    const u = new Uint8Array(buf, 0, 2);
    if (u[0] === 0x1f && u[1] === 0x8b) {            // still gzipped (server did not set Content-Encoding)
      buf = await new Response(new Blob([buf]).stream().pipeThrough(new DecompressionStream('gzip'))).arrayBuffer();
    }
    const dv = new DataView(buf);
    const h = { kind: dv.getUint8(4), narr: dv.getUint8(5) || 2, nmem: dv.getUint16(6, true), ny: dv.getUint16(8, true), nx: dv.getUint16(10, true),
      fhr: dv.getUint16(12, true), stepH: dv.getUint16(14, true), s0: dv.getFloat32(32, true), s1: dv.getFloat32(36, true),
      s2: dv.getFloat32(40, true), s3: dv.getFloat32(44, true), s4: dv.getFloat32(52, true),
      mask: dv.getUint32(48, true) + dv.getUint32(56, true) * 4294967296 };      // members 0-31 | members 32-63; 0 = legacy file, all members
    const n = h.nmem * h.ny * h.nx;
    if (h.kind === 1) {
      const o = { h }; ['q', 's', 'z', 'p', 'w'].slice(0, h.narr).forEach((k, i) => { o[k] = new Uint16Array(buf, 64 + 2 * n * i, n); });
      return o;
    }
    if (h.kind === 3) return { h, t: new Uint16Array(buf, 64, n) };          // daily max / min 2 m temperature, counts of s0 deg C from s1
    return { h, p: new Uint8Array(buf, 64, 4 * n) };
  });
}
const loadCum = (s, fhr) => fetchBin(s.base + s.byFhr.get(fhr).cum);
const loadPtype = (s, fhr) => fetchBin(s.base + s.byFhr.get(fhr).ptype);

/* ============================================================================= windows & statistics */
function windowFor(i, period) {
  const st = S.steps[i], end = st.fhr;
  let start = period === 'step' ? end - st.step_h : end - period, since = false;
  if (start < 0 || (start > 0 && !S.byFhr.has(start))) { since = true; start = 0; }
  return { end, start, since, hours: end - start };
}

/** Member totals over the window, raw units (mm for QPF, cm for snow), pooled over the selected model(s). Float32Array(M*N).
 *  `extra`: {source key: member bitmask} of members to leave out (e.g. those without p-type data).
 *  `members` are global indices (S.all); `parts` says which rows came from which source; `weights` is null for equal weighting. */
async function windowValues(kind, i, period, extra) {
  const w = windowFor(i, period), sc = src(kind), key = sc.key;
  const loaded = await Promise.all(S.src.map(async (s) => {
    const [a, b] = await Promise.all([loadCum(s, w.end), w.start > 0 ? loadCum(s, w.start) : null]);
    return { s, a, b };
  }));
  // Members missing from either file (or from `extra`) are left out; everything downstream works on the compacted set.
  const parts = [];
  for (const { s, a, b } of loaded) {
    if (!a[key]) throw new Error(`${PRODUCTS[kind].name}${key === 'w' ? ' at a fixed ratio' : ''} is not in ${FAMS[s.key].label}'s files yet (the next pipeline run adds it)`);
    const full = 2 ** a.h.nmem - 1, mk = (h) => h.mask || full;
    let mask = mk(a.h); if (b) mask = Number(BigInt(mask) & BigInt(mk(b.h))); if (extra && extra[s.key]) mask = Number(BigInt(mask) & BigInt(extra[s.key]));
    const ks = EnsStats.membersFromMask(mask, a.h.nmem);
    if (ks.length) parts.push({ s, a, b, ks });
  }
  if (!parts.length) throw new Error('no ensemble member has data for this window');
  const N = parts[0].a.h.ny * parts[0].a.h.nx, M = parts.reduce((t, p) => t + p.ks.length, 0), out = new Float32Array(M * N);
  const members = [], weights = parts.length > 1 && S.wmode === 'model' ? new Float32Array(M) : null;
  let j = 0;
  for (const p of parts) {
    const scale = p.a.h[sc.scale], A = p.a[key], B = p.b && p.b[key], j0 = j;
    p.j0 = j0;
    for (const k of p.ks) {
      const o = k * N, d = j * N;
      if (B) for (let q = 0; q < N; q++) out[d + q] = Math.max(0, A[o + q] - B[o + q]) * scale;
      else for (let q = 0; q < N; q++) out[d + q] = A[o + q] * scale;
      members.push(S.off[p.s.key] + k);
      if (weights) weights[j] = 1 / (parts.length * p.ks.length);       // each model counts the same
      j++;
    }
  }
  return { data: out, w, h: { ...parts[0].a.h, nmem: M }, members, weights, parts: parts.map((p) => ({ s: p.s, ks: p.ks, j0: p.j0 })) };
}

/* ---- daily maximum / minimum temperature (one file per product and date; values are deg F in memory) ---- */
const tempMode = () => TempDays.isTemp(S.product);
const currentDay = () => (S.tdays[S.product] || []).find((d) => d.date === S.day);
const dayPos = () => Math.max(0, (S.tdays[S.product] || []).findIndex((d) => d.date === S.day));
function fixDay() {                                 // keep the selected date valid for the product / model in view
  const d = S.tdays[S.product]; if (!d || !d.length || d.some((x) => x.date === S.day)) return;
  S.day = d[TempDays.nearestDay(d, S.day || d[0].date)].date;
}

async function tempValues(kind, date) {
  const day = (S.tdays[kind] || []).find((d) => d.date === date);
  if (!day) throw new Error('that day is not available');
  const loaded = await Promise.all(S.src.map(async (s, j) => ({ s, a: await fetchBin(s.base + day.entries[j].file) })));
  const parts = [];
  for (const { s, a } of loaded) {
    if (a.h.kind !== 3) throw new Error(`unexpected file for ${FAMS[s.key].label}`);
    const ks = EnsStats.membersFromMask(a.h.mask || 2 ** a.h.nmem - 1, a.h.nmem);
    if (ks.length) parts.push({ s, a, ks });
  }
  if (!parts.length) throw new Error('no ensemble member has data for this day');
  const N = parts[0].a.h.ny * parts[0].a.h.nx, M = parts.reduce((t, p) => t + p.ks.length, 0), out = new Float32Array(M * N);
  const members = [], weights = parts.length > 1 && S.wmode === 'model' ? new Float32Array(M) : null;
  let j = 0;
  for (const p of parts) {
    const h = p.a.h, f0 = h.s1 * 1.8 + 32, fs = h.s0 * 1.8, T = p.a.t;        // deg C = s1 + count * s0  ->  deg F
    p.j0 = j;
    for (const k of p.ks) {
      const o = k * N, d = j * N;
      for (let q = 0; q < N; q++) out[d + q] = f0 + T[o + q] * fs;
      members.push(S.off[p.s.key] + k);
      if (weights) weights[j] = 1 / (parts.length * p.ks.length);
      j++;
    }
  }
  return { data: out, w: { day: date, product: kind }, h: { ...parts[0].a.h, nmem: M }, members, weights, parts: parts.map((p) => ({ s: p.s, ks: p.ks, j0: p.j0 })) };
}

/** One deg F field from a day's member values: ensemble statistic, one member, or the chance of being at/below or at/above a threshold. */
function reduceTemp(win, stat) {
  const { nmem: M, ny, nx } = win.h, N = ny * nx, data = win.data, wt = win.weights;
  if (isMem(stat) && win.members.includes(+stat.slice(1))) return data.slice(win.members.indexOf(+stat.slice(1)) * N, (win.members.indexOf(+stat.slice(1)) + 1) * N);
  if (isMem(stat)) stat = 'mean';
  const pr = TempDays.parseProb(stat);
  if (pr) {
    const [dir, thr] = pr, le = dir === 'le', out = new Float32Array(N);
    for (let i = 0; i < N; i++) { let k = 0; for (let m = 0; m < M; m++) { const v = data[m * N + i]; if (le ? v <= thr : v >= thr) k += wt ? wt[m] : 1; } out[i] = wt ? 100 * k : (100 * k) / M; }
    return out;
  }
  return EnsStats.compute(data, M, ny, nx, stat, { weights: wt });
}

const toDisp = (kind, v) => v / siPerIn(kind);   // raw (mm or cm) -> inches

const pct = EnsStats.percentile;
const lpmmMi = () => Math.round((S.lpmmR * 34.5) / 10) * 10;   // radius in cells of 0.5° -> miles, rounded to match the menu

/** PMM / LPMM run in a Web Worker (falls back to the main thread). A newer request cancels an older one still running. */
let worker = null, wJob = 0, wPending = null;
function computeSlow(data, M, ny, nx, stat, radius, weights) {
  const here = () => EnsStats.compute(data, M, ny, nx, stat, { radius, weights });
  if (typeof Worker === 'undefined' || worker === false) return Promise.resolve(here());
  if (wPending) { wPending(null); wPending = null; worker.terminate(); worker = null; }      // drop the stale job
  try { worker = worker || new Worker('statworker.js'); } catch (e) { worker = false; return Promise.resolve(here()); }
  const id = ++wJob, w = worker;
  return new Promise((resolve) => {
    wPending = resolve;
    w.onmessage = (e) => { if (e.data.id !== id) return; wPending = null; resolve(e.data.error ? here() : e.data.out); };
    w.onerror = () => { wPending = null; worker = false; resolve(here()); };
    w.postMessage({ id, data, M, ny, nx, stat, radius, weights });
  });
}
const FIELD_MEMO = new Map();          // PMM / LPMM fields are the expensive ones: keep the last few

/** One display-unit (inches) field from the member totals of a window. */
async function reduceStat(kind, win, stat) {
  const { nmem: M, ny, nx } = win.h, N = ny * nx, data = win.data, c = 1 / siPerIn(kind);
  if (isMem(stat) && win.members.includes(+stat.slice(1))) { const out = new Float32Array(N), o = win.members.indexOf(+stat.slice(1)) * N; for (let i = 0; i < N; i++) out[i] = data[o + i] * c; return out; }
  if (isMem(stat)) stat = 'mean';                   // that member has no data for this window
  const wt = win.weights;
  if (stat.startsWith('prob:')) {
    const out = new Float32Array(N), thr = (+stat.slice(5)) / c;                 // inches -> raw
    for (let i = 0; i < N; i++) { let k = 0; for (let m = 0; m < M; m++) if (data[m * N + i] >= thr) k += wt ? wt[m] : 1; out[i] = wt ? 100 * k : (100 * k) / M; }
    return out;
  }
  const slow = stat === 'pmm' || stat === 'lpmm';
  const key = slow ? [sig(), kind, win.w.start, win.w.end, stat, S.lpmmR, kind === 'snow' ? S.slr : ''].join('|') : null;
  if (key && FIELD_MEMO.has(key)) return FIELD_MEMO.get(key);
  const out = slow ? await computeSlow(data, M, ny, nx, stat, S.lpmmR, wt) : EnsStats.compute(data, M, ny, nx, stat, { radius: S.lpmmR, weights: wt });
  if (!out) return null;                              // superseded while computing
  for (let i = 0; i < N; i++) out[i] *= c;
  if (key) { FIELD_MEMO.set(key, out); while (FIELD_MEMO.size > 10) FIELD_MEMO.delete(FIELD_MEMO.keys().next().value); }
  return out;
}

/* --- precipitation type ------------------------------------------------------------------ */
const domType = (p, N, m, i) => {                 // dominant type index for member m at cell i (-1 = dry)
  let best = -1, bv = 0;
  for (const t of TIE_ORDER) { const v = p[t * (p.length / 4) + m * N + i]; if (v > bv) { bv = v; best = t; } }
  return best;
};
const dbzOf = (rate, type) => (rate <= 0 ? 0 : type === 1 ? 32.5 + 22.1 * Math.log10(rate) : 23 + 16 * Math.log10(rate));

async function ptypeData() {
  const st = S.steps[S.idx];
  const pts = await Promise.all(S.src.map((s) => loadPtype(s, st.fhr)));
  const extra = {}; S.src.forEach((s, j) => { extra[s.key] = pts[j].h.mask; });
  const qw = await windowValues('qpf', S.idx, 'step', extra);
  const M = qw.h.nmem, N = qw.h.ny * qw.h.nx, p = new Uint8Array(4 * M * N);
  for (const part of qw.parts) {
    const pt = pts[S.src.indexOf(part.s)], full = pt.h.nmem * N;
    part.ks.forEach((k, jj) => { for (let t = 0; t < 4; t++) p.set(pt.p.subarray(t * full + k * N, t * full + (k + 1) * N), t * M * N + (part.j0 + jj) * N); });
  }
  return { p, q: qw.data, h: { ...qw.h }, w: qw.w, members: qw.members, weights: qw.weights, parts: qw.parts };
}

function buildPtype(d, stat) {
  if (isMem(stat) && !d.members.includes(+stat.slice(1))) stat = 'ens';
  const { p, q, h, w } = d, M = h.nmem, N = h.ny * h.nx, rateDiv = w.hours || 3, wt = d.weights, wm = (m) => (wt ? wt[m] : 1 / M);
  const field = new Float32Array(N), types = new Float32Array(N), tot = M * N;
  const share = (m, i, t) => { let s = 0; for (let k = 0; k < 4; k++) s += p[k * tot + m * N + i]; return s ? p[t * tot + m * N + i] / s : 0; };
  for (let i = 0; i < N; i++) {
    if (isMem(stat) && d.members.includes(+stat.slice(1))) {
      const m = d.members.indexOf(+stat.slice(1)), t = domType(p, N, m, i);
      if (t >= 0) { field[i] = dbzOf(q[m * N + i] / rateDiv, t); types[i] = t; }
    } else if (stat.startsWith('prob:')) {
      const want = +stat.slice(5); let k = 0;
      for (let m = 0; m < M; m++) if (domType(p, N, m, i) === want) k += wm(m);
      field[i] = 100 * k; types[i] = want;
    } else {                                              // 'ens': QPF-weighted dominant type, mean-QPF intensity
      const wsum = [0, 0, 0, 0]; let qs = 0;
      for (let m = 0; m < M; m++) { const qq = q[m * N + i] * wm(m); qs += qq; if (qq > 0) for (let t = 0; t < 4; t++) wsum[t] += qq * share(m, i, t); }
      let best = -1, bv = 0; for (const t of TIE_ORDER) if (wsum[t] > bv) { bv = wsum[t]; best = t; }
      if (best >= 0) { field[i] = dbzOf(qs / rateDiv, best); types[i] = best; }
    }
  }
  return { field, types };
}

/* ============================================================================= rendering */
const PROB_COLD = ['#eff3ff', '#c6dbef', '#9ecae1', '#6baed6', '#4292c6', '#2171b5', '#08519c', '#08306b', '#041f4a'];
const PROB_WARM = ['#fff5eb', '#fee6ce', '#fdd0a2', '#fdae6b', '#fd8d3c', '#f16913', '#d94801', '#a63603', '#7f2704'];
function classesFor() {
  if (tempMode()) {
    if (!isProb(S.stat)) return { ...TempDays.CLASSES, title: TempDays.NAMES[S.product], temp: true };
    const pr = TempDays.parseProb(S.stat) || ['ge', 0];
    return { edges: PROB_EDGES, colors: pr[0] === 'le' ? PROB_COLD : PROB_WARM, labels: PROB_EDGES.map((e) => `${e}%+`), prob: true };
  }
  if (isProb(S.stat)) {
    const k = PROB_RAMP[S.product] || ptypeProbColors(S.product === 'zr' ? 2 : S.product === 'pl' ? 3 : +S.stat.split(':')[1]);
    return { edges: PROB_EDGES, colors: k, labels: PROB_EDGES.map((e) => `${e}%+`), prob: true };
  }
  return RAMPS[S.product];
}
function ptypeProbColors(t) { const r = PTYPE.ramps[PTYPE.names[t]]; return PROB_EDGES.map((_, k) => r[Math.round(4 + (k * (r.length - 12)) / 8)]); }

async function refresh() {
  const tok = ++S.token; $('#busy').hidden = false;
  try {
    if (S.product === 'ptype') {
      const d = await ptypeData(); if (tok !== S.token) return;
      const { field, types } = buildPtype(d, S.stat);
      if (isProb(S.stat)) { const c = classesFor(); wx.setClasses(c.edges, c.colors); }
      else wx.setPtypeLut(PTYPE.ramps, PTYPE.names, PTYPE.dbzMin, PTYPE.dbzStep);
      wx.setField(field, types); S.cur = { kind: 'ptype', d, field, types };
    } else if (tempMode()) {
      const win = await tempValues(S.product, S.day); if (tok !== S.token) return;
      const field = reduceTemp(win, S.stat);
      const c = classesFor(); wx.setClasses(c.edges, c.colors); wx.setField(field);
      S.cur = { kind: S.product, win, field };
    } else {
      const win = await windowValues(S.product, S.idx, S.period); if (tok !== S.token) return;
      const field = await reduceStat(S.product, win, S.stat); if (tok !== S.token || !field) return;
      const c = classesFor(); wx.setClasses(c.edges, c.colors); wx.setField(field);
      S.cur = { kind: S.product, win, field };
    }
    updateLegend(); updateTitle(); updateTip(); updateStations(); prefetch();
  } catch (e) { console.error(e); toast(`Could not load data: ${e.message}`); }
  finally { if (tok === S.token) $('#busy').hidden = true; }
}

/* Idle-time prefetch: quietly load the files the next click is most likely to need, one at a time, lowest priority.
 * Order: neighbouring hours for the current view, the start files of the other accumulation periods, then precip type. */
let pfGen = 0, lastIdx = 0, dir = 1;
function prefetch() {
  const pos = tempMode() ? dayPos() : S.idx;
  if (pos !== lastIdx) dir = Math.sign(pos - lastIdx) || dir; lastIdx = pos;
  const gen = ++pfGen, nav = navigator.connection || {};
  if (nav.saveData || /(^|-)2g$/.test(nav.effectiveType || '')) return;
  if (tempMode()) {                                   // neighbouring days (travel direction first), then the same day of the other temperature product
    const turls = [], days = S.tdays[S.product], i0 = dayPos(), other = S.tdays[S.product === 'tmax' ? 'tmin' : 'tmax'];
    [1, 2, -1].map((d) => d * dir).forEach((d) => { const e = days[i0 + d]; if (e) e.entries.forEach((en, j) => turls.push(S.src[j].base + en.file)); });
    const o = other && other.find((x) => x.date === S.day); if (o) o.entries.forEach((en, j) => turls.push(S.src[j].base + en.file));
    return drain(gen, turls);
  }
  const urls = [], add = (fhr) => { if (fhr != null && S.byFhr.has(fhr)) S.src.forEach((s) => urls.push(s.base + s.byFhr.get(fhr).cum)); };
  const addPt = (fhr) => { if (S.byFhr.has(fhr)) S.src.forEach((s) => urls.push(s.base + s.byFhr.get(fhr).ptype)); };
  const per = S.product === 'ptype' ? 'step' : S.period;
  [1, 2, 3, -1].map((d) => d * dir).forEach((d) => { const j = S.idx + d; if (j < 0 || j >= S.steps.length) return; const w = windowFor(j, per); add(w.end); if (w.start > 0) add(w.start); if (S.product === 'ptype') addPt(w.end); });
  const here = S.steps[S.idx];
  [6, 12, 24, 48, 'step'].forEach((p) => { const w = windowFor(S.idx, p); if (w.start > 0) add(w.start); });
  addPt(here.fhr);
  [1, -1].forEach((d) => { const j = S.idx + d; if (j >= 0 && j < S.steps.length) addPt(S.steps[j].fhr); });
  drain(gen, urls);
}
function drain(gen, urls) {
  const queue = [...new Set(urls)].filter((u) => !CACHE.has(u));
  const idle = window.requestIdleCallback || ((f) => setTimeout(f, 200));
  const next = () => idle(() => {
    if (gen !== pfGen || !queue.length) return;
    fetchBin(queue.shift()).catch(() => {}).then(next);
  }, { timeout: 1500 });
  next();
}

/* ============================================================================= labels / titles */
const tz = Intl.DateTimeFormat().resolvedOptions().timeZone;
const fmtLocal = (d) => d.toLocaleString([], { weekday: 'short', month: 'short', day: 'numeric', hour: 'numeric', timeZoneName: 'short' });
const fmtUtc = (d) => `${d.toISOString().slice(5, 10).replace('-', '/')} ${d.toISOString().slice(11, 13)}Z`;
const stepN = (i) => (S.steps[i].n ?? S.all.length);
const memLabel = (g) => { const m = S.all[g]; return m ? (S.src.length > 1 ? `${FAMS[m.fam].label} ${m.id}` : m.id) : `member ${g}`; };
const initTxt = () => { const d = new Date(S.init); return `${d.toISOString().slice(0, 10)} ${d.toISOString().slice(11, 13)}Z`; };
const famTitle = () => (S.fam === 'blend' ? `Blend of ${S.src.map((s) => FAMS[s.key].label).join(' + ')}` : FAMS[S.fam].long);
const pad = (n, w = 3) => String(n).padStart(w, '0');

function productName() {
  if (S.product === 'ptype') return 'Precipitation type';
  const nm = PRODUCTS[S.product].name + (S.product === 'snow' && S.hasSwe ? (S.slr === 'kuchera' ? ' (Kuchera ratio)' : ` (${S.slr}:1)`) : '');
  const w = windowFor(S.idx, S.period);
  return `${w.since ? `${w.hours}-hr (since init)` : `${w.hours}-hr`} ${nm}`;
}
function statName() {
  const s = S.stat;
  if (s === 'mean') return 'Ensemble mean'; if (s === 'max') return 'Ensemble maximum'; if (s === 'min') return 'Ensemble minimum';
  if (s === 'pmm') return 'Probability-matched mean (PMM)';
  if (s === 'lpmm') return `Local prob.-matched mean (LPMM \u00b1${lpmmMi()} mi)`;
  if (EnsStats.isPct(s)) return s === 'p50' ? '50th percentile (median)' : `${s.slice(1)}th percentile`;
  if (s === 'ens') return S.src.length > 1 ? 'Most likely type (blend)' : 'Most likely type (ensemble)';
  const tp = TempDays.parseProb(s);
  if (tp) return `Probability ${tp[0] === 'le' ? '\u2264' : '\u2265'} ${tp[1]}\u00b0F`;
  if (s.startsWith('prob:')) {
    const v = +s.slice(5);
    return S.product === 'ptype' ? `Probability of ${PTYPE.long[v].toLowerCase()}` : `Probability ≥ ${v}″`;
  }
  const g = +s.slice(1), m = S.all[g];
  return m && m.id === 'c00' ? `Control (${memLabel(g)})` : `Member ${memLabel(g)}`;
}

function updateWarn() {
  const tm = tempMode(), st = tm ? null : S.steps[S.idx], el = $('#t-warn'), bits = [];
  S.src.forEach((s, j) => {
    const part = tm ? currentDay().entries[j] : st.parts[s.key], tot = s.man.members.length;
    if (!part || part.n === undefined || part.n >= tot) return;
    const gone = s.man.members.filter((m, k) => part.mask && !(Math.floor(part.mask / 2 ** k) % 2)).join(', ');
    bits.push(`${S.src.length > 1 ? FAMS[s.key].label + ': ' : ''}${part.n} of ${tot} members (missing: ${gone})`);
  });
  if (!bits.length) { el.hidden = true; return; }
  el.hidden = false; el.textContent = `\u26a0 ${bits.join('; ')} ${tm ? 'for this day' : 'at this hour'}. Statistics use the members that are available.`;
}

function updateTempTitle() {
  const day = currentDay(), [a, b] = TempDays.windowMs(S.product, S.day), A = new Date(a), B = new Date(b);
  const hh = (d) => `${String(d.getUTCHours()).padStart(2, '0')}Z`, short = (d) => d.toLocaleTimeString([], { hour: 'numeric', timeZoneName: 'short' });
  const n = day.entries.reduce((t, e) => t + e.n, 0);
  $('#t-main').textContent = `${statName()} \u00b7 ${TempDays.NAMES[S.product]}`;
  $('#t-valid').textContent = `${fmtLocal(A)} \u2192 ${fmtLocal(B)} (${fmtUtc(A)} \u2192 ${fmtUtc(B)})`;
  $('#t-init').textContent = `${famTitle()} \u00b7 ${initTxt()} run \u00b7 ${n} members${S.src.length > 1 ? (S.wmode === 'model' ? ' \u00b7 models weighted equally' : ' \u00b7 members weighted equally') : ''}`;
  updateWarn();
  $('#fh').textContent = new Date(`${S.day}T12:00:00Z`).toLocaleDateString([], { weekday: 'short', month: 'short', day: 'numeric', timeZone: 'UTC' });
  $('#fh-time').textContent = `${S.product === 'tmax' ? 'MaxT' : 'MinT'} ${hh(A)}\u2013${hh(B)}`;
  $('#fh-utc').textContent = `${short(A)} \u2192 ${short(B)}`;
  $('#slider').value = dayPos();
}

function updateTitle() {
  if (tempMode()) return updateTempTitle();
  const st = S.steps[S.idx], valid = new Date(st.valid), init = new Date(S.init);
  const w = windowFor(S.idx, S.product === 'ptype' ? 'step' : S.period);
  const startD = new Date(init.getTime() + w.start * 3600e3);
  $('#t-main').textContent = `${statName()} · ${productName()}`;
  $('#t-valid').textContent = S.product === 'ptype' || w.hours <= st.step_h
    ? `Valid ${fmtLocal(valid)} (${fmtUtc(valid)})`
    : `${fmtLocal(startD)} → ${fmtLocal(valid)}`;
  $('#t-init').textContent = `${famTitle()} · ${initTxt()} run · F${pad(st.fhr)} · ${stepN(S.idx)} members${S.src.length > 1 ? (S.wmode === 'model' ? ' · models weighted equally' : ' · members weighted equally') : ''}`;
  updateWarn();
  $('#fh').textContent = `F${pad(st.fhr)}`;
  $('#fh-time').textContent = `${fmtLocal(valid)}`;
  $('#fh-utc').textContent = fmtUtc(valid);
  $('#slider').value = S.idx;
}


function updateLegend() {
  const el = $('#legend-body'); el.innerHTML = '';
  const hdr = $('#legend-title');
  if (S.product === 'ptype' && !isProb(S.stat)) {
    hdr.textContent = 'Type color ramps — step QPF as dBZ';
    PTYPE.names.forEach((n, t) => {
      const row = document.createElement('div'); row.className = 'pt-row';
      row.innerHTML = `<span class="pt-name">${PTYPE.long[t]}</span>`;
      const bar = document.createElement('div'); bar.className = 'pt-bar';
      PTYPE.ramps[n].forEach((c, k) => { const s = document.createElement('i'); s.style.background = k === 0 ? 'repeating-linear-gradient(45deg,#fff,#fff 2px,#e7e7e7 2px,#e7e7e7 4px)' : c; bar.appendChild(s); });
      row.appendChild(bar); el.appendChild(row);
    });
    const ticks = document.createElement('div'); ticks.className = 'pt-ticks';
    [10, 20, 30, 40, 50, 60, 70].forEach((v) => { const s = document.createElement('span'); s.textContent = v; s.style.left = `${((v - 8) / 66) * 100}%`; ticks.appendChild(s); });
    el.appendChild(ticks);
    const note = document.createElement('div'); note.className = 'legend-note';
    note.textContent = 'dBZ ≈ rate from step QPF (rain/ice: Z=200R^1.6, snow: Z=1780R^2.21)'; el.appendChild(note); return;
  }
  const c = classesFor();
  if (c.temp) {                                    // 30 WPC classes as one bar, warmest on top, labelled every 20 deg F
    hdr.textContent = `${c.title} (\u00b0F)`;
    const wrap = document.createElement('div'), bar = document.createElement('div'); wrap.className = 'tl'; bar.className = 'tl-bar';
    for (let i = c.colors.length - 1; i >= 0; i--) { const q = document.createElement('i'); q.style.background = c.colors[i]; bar.appendChild(q); }
    const ticks = document.createElement('div'), n = c.colors.length; ticks.className = 'tl-ticks';
    for (let v = 100; v >= -40; v -= 20) { const t = document.createElement('span'); t.textContent = v; t.style.top = `${(((100 - v) / 5 + 1) / n) * 100}%`; ticks.appendChild(t); }
    wrap.append(bar, ticks); el.appendChild(wrap); return;
  }
  hdr.textContent = c.prob ? `Probability (% of ${stepN(S.idx)} members${S.src.length > 1 && S.wmode === 'model' ? ', models weighted equally' : ''})` :
    `${c.title} (in)`;
  for (let i = c.colors.length - 1; i >= 0; i--) {
    const row = document.createElement('div'); row.className = 'lg-row';
    row.innerHTML = `<i style="background:${c.colors[i]}"></i><span>${c.labels[i]}</span>`; el.appendChild(row);
  }
}

/* ============================================================================= station values */
/** Print the value of the map being shown at each New England METAR/ASOS site (same field, same grid interpolation as the colours). */
function updateStations() {
  if (!labels) return;
  const cur = S.cur;
  if (!S.stn || !S.stations || !cur || !cur.field || !S.domain) { labels.stations = null; labels.q(); return; }
  const N = S.domain.nx * S.domain.ny, prob = isProb(S.stat);
  labels.stations = S.stations.map((s) => {
    const c = cellValues({ lng: s.lon, lat: s.lat }, cur.field, 1, N);
    if (!c) return { id: s.id, lon: s.lon, lat: s.lat, txt: '' };
    const type = cur.kind === 'ptype' && cur.types ? cur.types[c.nearest] : -1;
    return { id: s.id, lon: s.lon, lat: s.lat, txt: Stations.format({ product: cur.kind, prob, v: c.v[0], type }) };
  });
  labels.q();
}
function loadStations() {
  return fetch('stations.json').then((r) => (r.ok ? r.json() : Promise.reject(new Error(r.status)))).then((l) => { S.stations = Stations.order(l); })
    .catch(() => { S.stations = null; const el = $('#tg-stations'); if (el) el.closest('label').hidden = true; });
}

/* ============================================================================= tooltip */
const tip = () => $('#tip');
let touchAt = 0;
addEventListener('touchstart', () => { touchAt = Date.now(); }, { capture: true, passive: true });
const recentTouch = () => Date.now() - touchAt < 1200;   // ignore emulated mouse events that follow a tap
let lastLL = null, raf = 0;

function cellValues(lngLat, data, M, N, h) {
  const d = S.domain, gx = (lngLat.lng - d.lon_west) / d.dlon, gy = (d.lat_north - lngLat.lat) / d.dlat;
  if (gx < 0 || gy < 0 || gx > d.nx - 1 || gy > d.ny - 1) return null;
  const x0 = Math.min(Math.floor(gx), d.nx - 2), y0 = Math.min(Math.floor(gy), d.ny - 2), fx = gx - x0, fy = gy - y0;
  const i00 = y0 * d.nx + x0, w = [(1 - fx) * (1 - fy), fx * (1 - fy), (1 - fx) * fy, fx * fy], off = [0, 1, d.nx, d.nx + 1];
  const out = new Float32Array(M);
  for (let m = 0; m < M; m++) { let s = 0; for (let k = 0; k < 4; k++) s += w[k] * data[m * N + i00 + off[k]]; out[m] = s; }
  return { v: out, nearest: i00 + (fx > 0.5 ? 1 : 0) + (fy > 0.5 ? d.nx : 0) };
}

function placeName(pt) {
  const f = map.queryRenderedFeatures(pt, { layers: ['counties-hit', 'admin1-hit'] });
  const county = f.find((x) => x.layer.id === 'counties-hit'), st = f.find((x) => x.layer.id === 'admin1-hit');
  const cName = county && county.properties.n, stName = st && st.properties.n;
  const cc = st && ({ USA: 'United States', CAN: 'Canada', MEX: 'Mexico' }[st.properties.c]);
  return { line1: cName ? `${cName}${county.properties.s ? ', ' + county.properties.s : ''}` : stName || '', line2: [cName ? stName : '', cc].filter(Boolean).join(' · ') };
}

/* Display precision: snowfall to tenths of an inch; QPF, freezing rain and sleet (liquid) to hundredths. */
const fmtVal = (kind, inches) => (TempDays.isTemp(kind) ? `${inches.toFixed(1)}\u00b0F` : `${inches.toFixed(kind === 'snow' ? 1 : 2)}″`);   // temperatures to tenths

function stripSVG(vals, mean, med, kind, selected, famOf, lo = 0, hi = null) {      // lo/hi: axis range (temperatures); precipitation runs 0..max
  const W = 268, H = S.src.length > 1 ? 52 : 40, pad = 6, n = vals.length, mx = hi == null ? Math.max(...vals, 1e-9) * 1.04 : hi;
  const sx = (v) => pad + ((v - lo) / (mx - lo || 1)) * (W - 2 * pad), bins = new Map(), dots = [], step = n > 60 ? 4.4 : 5.2;
  vals.forEach((v, m) => { const b = Math.round(sx(v) / 5); const k = bins.get(b) || 0; bins.set(b, k + 1); dots.push([sx(v), H - 8 - k * step, m]); });
  let s = `<svg viewBox="0 0 ${W} ${H}" width="100%" role="img" aria-label="Ensemble member values">`;
  s += `<line x1="${pad}" x2="${W - pad}" y1="${H - 3}" y2="${H - 3}" stroke="#9aa3ab" stroke-width="1"/>`;
  dots.forEach(([x, y, m]) => { s += `<circle cx="${x.toFixed(1)}" cy="${Math.max(4, y).toFixed(1)}" r="2.3" fill="${m === selected ? '#111' : famOf ? FAMS[famOf[m]].color : '#5b6670'}" ${m === selected ? 'stroke="#fff" stroke-width="1"' : ''}/>`; });
  s += `<line x1="${sx(mean)}" x2="${sx(mean)}" y1="2" y2="${H - 3}" stroke="#c2410c" stroke-width="1.5"/>`;
  s += `<line x1="${sx(med)}" x2="${sx(med)}" y1="2" y2="${H - 3}" stroke="#1d4ed8" stroke-width="1.5" stroke-dasharray="2 2"/>`;
  return s + '</svg>';
}

/** Model of every pooled member (`members` are global indices) and a one-line "GEFS 0.52″ · IFS ..." summary of per-model means. */
function famInfo(members, vals, kind) {
  const famOf = members.map((g) => S.all[g].fam);
  if (S.src.length < 2) return { famOf: null, line: '' };
  const acc = {}; famOf.forEach((f, j) => { const a = acc[f] || (acc[f] = [0, 0]); a[0] += vals[j]; a[1]++; });
  const line = S.src.map((s) => acc[s.key] ? `<span><i class="sw" style="background:${FAMS[s.key].color}"></i>${FAMS[s.key].label} <b>${fmtVal(kind, acc[s.key][0] / acc[s.key][1])}</b> <small>(${acc[s.key][1]})</small></span>` : '').join('');
  return { famOf, line: `<div class="tip-fam">${line}</div>` };
}

function tipHTML(ll, pt) {
  const cur = S.cur; if (!cur) return '';
  const nm = placeName(pt);
  const head = `<div class="tip-h"><b>${nm.line1 || 'Open water / outside boundaries'}</b>${nm.line2 ? `<span>${nm.line2}</span>` : ''}
    <small>${Math.abs(ll.lat).toFixed(2)}°${ll.lat >= 0 ? 'N' : 'S'}, ${Math.abs(ll.lng).toFixed(2)}°${ll.lng >= 0 ? 'E' : 'W'}</small></div>`;
  if (cur.kind === 'ptype') return head + ptypeTip(ll, cur.d);
  if (TempDays.isTemp(cur.kind)) return head + tempTip(ll, cur);
  const { win } = cur, M = win.h.nmem, N = win.h.ny * win.h.nx, cell = cellValues(ll, win.data, M, N, win.h);
  if (!cell) return head + '<div class="tip-b muted">Outside the domain.</div>';
  const wt = win.weights, wsum = wt ? Array.from(wt).reduce((a, b) => a + b, 0) : M;
  const vals = Array.from(cell.v, (x) => toDisp(S.product, x)), sorted = Float32Array.from(vals).sort();
  const mean = wt ? vals.reduce((a, b, m) => a + b * wt[m], 0) / wsum : vals.reduce((a, b) => a + b, 0) / M;
  const wp = (p) => (wt ? EnsStats.wPercentile(vals, wt, p) : pct(sorted, M, p)), med = wp(0.5);
  const kind = S.product, f = (x) => fmtVal(kind, x);
  const sel = isMem(S.stat) ? win.members.indexOf(+S.stat.slice(1)) : -1;
  const val = { min: sorted[0], max: sorted[M - 1], mean };
  EnsStats.PCTS.forEach((p) => { val[`p${p}`] = wp(p / 100); });
  const nm2 = (k) => (k === 'mean' ? 'Mean' : k === 'min' ? 'Min' : k === 'max' ? 'Max' : k.toUpperCase());
  const ord = EnsStats.TIP_ORDER, half = ord.length / 2, fi = famInfo(win.members, vals, kind);
  let h = head + `<div class="tip-t">${productName()} \u00b7 ${$('#t-valid').textContent}</div>`;
  if ((S.stat === 'pmm' || S.stat === 'lpmm') && cur.field) {
    const fv = cellValues(ll, cur.field, 1, N, win.h);
    if (fv) h += `<div class="tip-map"><span>${S.stat === 'pmm' ? 'PMM' : `LPMM \u00b1${lpmmMi()} mi`}</span><b>${f(fv.v[0])}</b></div>`;
  }
  h += '<div class="tip-grid">';
  for (let r = 0; r < half; r++) h += `<span>${nm2(ord[r])}</span><b>${f(val[ord[r]])}</b><span>${nm2(ord[r + half])}</span><b>${f(val[ord[r + half]])}</b>`;
  h += `</div>${fi.line}<div class="tip-s">Chance of at least</div><div class="tip-p">`;
  PROB_THRESH[kind].forEach((t) => {
    let k = 0; vals.forEach((v, m) => { if (v >= t) k += wt ? wt[m] : 1; });
    const p = Math.round((100 * k) / wsum);
    h += `<span>${t}″</span><i><u style="width:${p}%"></u></i><b>${p}%</b>`;
  });
  h += `</div><div class="tip-s">All ${M} members <span class="key"><i class="k1"></i>mean <i class="k2"></i>median</span></div>${stripSVG(vals, mean, med, kind, sel, fi.famOf)}`;
  h += `<div class="tip-ax"><span>0</span><span>${f(Math.max(...vals))}</span></div>`;
  return h;
}

function tempTip(ll, cur) {
  const { win } = cur, M = win.h.nmem, N = win.h.ny * win.h.nx, cell = cellValues(ll, win.data, M, N, win.h);
  if (!cell) return '<div class="tip-b muted">Outside the domain.</div>';
  const wt = win.weights, wsum = wt ? Array.from(wt).reduce((a, b) => a + b, 0) : M, kind = S.product, f = (x) => fmtVal(kind, x);
  const vals = Array.from(cell.v), sorted = Float32Array.from(vals).sort(), lo = sorted[0], hi = sorted[M - 1];
  const mean = wt ? vals.reduce((a, b, m) => a + b * wt[m], 0) / wsum : vals.reduce((a, b) => a + b, 0) / M;
  const wp = (p) => (wt ? EnsStats.wPercentile(vals, wt, p) : pct(sorted, M, p)), med = wp(0.5);
  const sel = isMem(S.stat) ? win.members.indexOf(+S.stat.slice(1)) : -1, val = { min: lo, max: hi, mean };
  EnsStats.PCTS.forEach((p) => { val[`p${p}`] = wp(p / 100); });
  const nm2 = (k) => (k === 'mean' ? 'Mean' : k === 'min' ? 'Min' : k === 'max' ? 'Max' : k.toUpperCase());
  const ord = EnsStats.TIP_ORDER, half = ord.length / 2, fi = famInfo(win.members, vals, kind);
  let h = `<div class="tip-t">${TempDays.NAMES[kind]} \u00b7 ${$('#fh').textContent} ${$('#fh-time').textContent.split(' ')[1]}</div><div class="tip-grid">`;
  for (let r = 0; r < half; r++) h += `<span>${nm2(ord[r])}</span><b>${f(val[ord[r]])}</b><span>${nm2(ord[r + half])}</span><b>${f(val[ord[r + half]])}</b>`;
  h += `</div>${fi.line}<div class="tip-s">Chance of</div><div class="tip-p">`;
  TempDays.PROB[kind].forEach((pr) => {
    let k = 0; vals.forEach((v, m) => { if (pr[0] === 'le' ? v <= pr[1] : v >= pr[1]) k += wt ? wt[m] : 1; });
    const p = Math.round((100 * k) / wsum);
    h += `<span>${TempDays.probLabel(pr)}</span><i><u style="width:${p}%"></u></i><b>${p}%</b>`;
  });
  h += `</div><div class="tip-s">All ${M} members <span class="key"><i class="k1"></i>mean <i class="k2"></i>median</span></div>${stripSVG(vals, mean, med, kind, sel, fi.famOf, lo, hi > lo ? hi : lo + 1)}`;
  return h + `<div class="tip-ax"><span>${f(lo)}</span><span>${f(hi)}</span></div>`;
}

function ptypeTip(ll, d) {
  const { p, q, h } = d, M = h.nmem, N = h.ny * h.nx, c = cellValues(ll, q, M, N, h), wt = d.weights, wm = (m) => (wt ? wt[m] : 1 / M);
  if (!c) return '<div class="tip-b muted">Outside the domain.</div>';
  const i = c.nearest, cnt = [0, 0, 0, 0], pr = [0, 0, 0, 0], avg = [0, 0, 0, 0]; let wet = 0, qs = 0, wetW = 0;
  for (let m = 0; m < M; m++) {
    qs += c.v[m] * wm(m);
    const t = domType(p, N, m, i); if (t < 0) continue; cnt[t]++; pr[t] += wm(m); wet++; wetW += wm(m);
    for (let k = 0; k < 4; k++) avg[k] += wm(m) * p[k * M * N + m * N + i];
  }
  const order = [1, 0, 2, 3], pc = (x) => Math.round(100 * x);
  let h2 = `<div class="tip-t">Precip type · ${$('#t-valid').textContent}</div>`;
  h2 += `<div class="tip-grid"><span>Members with precip</span><b>${wet} / ${M}</b><span>Mean ${d.w.hours}-hr QPF</span><b>${fmtVal('qpf', qs / MM_PER_IN)}</b></div>`;
  h2 += `<div class="tip-s">Dominant type (${wt ? 'members \u00b7 share of the weighted ensemble' : 'members'})</div><div class="tip-p">`;
  order.forEach((t) => { h2 += `<span><i class="sw" style="background:${PTYPE.ramps[PTYPE.names[t]][18]}"></i>${PTYPE.long[t]}</span><i><u style="width:${wt ? pc(pr[t]) : (100 * cnt[t]) / M}%"></u></i><b>${wt ? `${cnt[t]} \u00b7 ${pc(pr[t])}%` : cnt[t]}</b>`; });
  h2 += '</div>';
  if (S.src.length > 1) {
    const by = S.src.map((s) => { const part = d.parts.find((x) => x.s === s); if (!part) return ''; let n = 0; part.ks.forEach((k, jj) => { if (domType(p, N, part.j0 + jj, i) >= 0) n++; }); return `<span><i class="sw" style="background:${FAMS[s.key].color}"></i>${FAMS[s.key].label} <b>${n}</b>/${part.ks.length} wet</span>`; }).join('');
    h2 += `<div class="tip-fam">${by}</div>`;
  }
  h2 += '<div class="tip-s">Mean Bourgouin probability (wet members)</div><div class="tip-grid">';
  order.forEach((t) => { h2 += `<span>${PTYPE.long[t]}</span><b>${wetW ? Math.round(avg[t] / wetW) : 0}%</b>`; });
  return h2 + '</div>';
}

function updateTip() {
  const el = tip(); const ll = S.pin || lastLL;
  if (!ll || (!S.pin && !el.classList.contains('show'))) return;
  const pt = map.project(ll); el.innerHTML = tipHTML(ll, pt) + (S.pin ? '<button id="tip-x" aria-label="Close">×</button>' : '');
  if (S.pin) { const x = $('#tip-x'); if (x) x.onclick = unpin; }
}
function showTip(e) {
  if (S.pin || !S.cur) return;
  lastLL = e.lngLat; const el = tip();
  cancelAnimationFrame(raf);
  raf = requestAnimationFrame(() => {
    el.innerHTML = tipHTML(e.lngLat, e.point); el.classList.add('show');
    const r = el.getBoundingClientRect(), W = innerWidth, H = innerHeight;
    let x = e.originalEvent.clientX + 18, y = e.originalEvent.clientY + 18;
    if (x + r.width > W - 8) x = e.originalEvent.clientX - r.width - 18;
    if (y + r.height > H - 76) y = Math.max(8, H - r.height - 76);
    el.style.left = `${Math.max(8, x)}px`; el.style.top = `${y}px`;
  });
}
function hideTip() { if (!S.pin) { tip().classList.remove('show'); lastLL = null; } }
function pin(e) {
  if (S.pin && !recentTouch()) { unpin(); return; }
  S.pin = e.lngLat; tip().classList.add('show', 'pinned');
  if (!pinMarker) pinMarker = new maplibregl.Marker({ element: Object.assign(document.createElement('div'), { className: 'pin' }) });
  pinMarker.setLngLat(e.lngLat).addTo(map); updateTip();
  const el = tip(); el.style.left = `${Math.max(8, Math.min(e.originalEvent.clientX + 18, innerWidth - 330))}px`; el.style.top = `${Math.max(8, Math.min(e.originalEvent.clientY + 18, innerHeight - 470))}px`;
}
function unpin() { S.pin = null; tip().classList.remove('show', 'pinned'); if (pinMarker) pinMarker.remove(); }

/* ============================================================================= models, sources and steps */
const hhZ = (iso) => `${iso.slice(5, 10).replace('-', '/')} ${iso.slice(11, 13)}Z`;
function makeSource(key, man) {
  const byFhr = new Map(); man.steps.forEach((st) => byFhr.set(st.fhr, st));
  return { key, man, base: `${man.legacy ? 'data' : `data/${key}`}/${man.base}`, byFhr, steps: man.steps };
}
const loadedFams = () => FAM_KEYS.filter((k) => S.sets[k] && S.sets[k].latest && S.sets[k].latest.steps.length);

/** Which model(s) feed a selection.  The blend takes the run time shared by the most models (newest on a tie): IFS and GEPS
 *  exist only at 00/12Z, so GEFS contributes its 00/12Z cycle (synoptic.json) and a model still on an older run is left out. */
function pickSources(fam) {
  if (fam !== 'blend') { const m = S.sets[fam] && S.sets[fam].latest; return { src: m && m.steps.length ? [makeSource(fam, m)] : [], note: '' }; }
  const cand = FAM_KEYS.map((k) => { const x = S.sets[k]; const m = x && (k === 'gefs' ? (x.syn || x.latest) : x.latest); return m && m.steps.length ? { k, m } : null; }).filter(Boolean);
  const groups = new Map(); cand.forEach((c) => groups.set(c.m.init, [...(groups.get(c.m.init) || []), c]));
  let best = null; groups.forEach((g, init) => { if (!best || g.length > best.g.length || (g.length === best.g.length && init > best.init)) best = { init, g }; });
  if (!best || best.g.length < 2) return { src: [], note: 'The blend needs at least two models with the same run time.' };
  const left = cand.filter((c) => !best.g.includes(c)).map((c) => `${FAMS[c.k].label} is still on its ${hhZ(c.m.init)} run`);
  return { src: best.g.map((c) => makeSource(c.k, c.m)), note: left.length ? `Not in the blend yet: ${left.join('; ')}.` : '' };
}

function applySources(src) {
  S.src = src; S.all = []; S.off = {};
  src.forEach((s) => { S.off[s.key] = S.all.length; s.man.members.forEach((id) => S.all.push({ fam: s.key, id })); });
  S.init = src[0].man.init; S.domain = src[0].man.domain;
  S.complete = src.every((s) => s.man.complete);
  S.gen = src.map((s) => s.man.generated).sort().pop();
  const sets = src.map((s) => new Set(s.man.products || ['qpf', 'snow', 'ptype']));
  S.products = [...sets[0]].filter((p) => sets.every((x) => x.has(p)));
  S.tdays = {};                                   // dates every pooled model has complete, per temperature product
  TempDays.PRODUCTS.forEach((p) => { S.tdays[p] = TempDays.commonDays(p, src.map((s) => s.man.temps)); });
  S.products = S.products.filter((p) => !TempDays.isTemp(p) || S.tdays[p].length);
  S.hasSwe = src.every((s) => !!(s.man.scales && s.man.scales.swe_mm_per_count));
  if (src.length === 1) S.steps = src[0].steps.map((st) => ({ ...st, parts: { [src[0].key]: st } }));
  else {                                          // hours every selected model has; a step spans back to the previous shared hour
    const common = src[0].steps.map((st) => st.fhr).filter((f) => src.every((s) => s.byFhr.has(f)));
    S.steps = common.map((f, i) => {
      const parts = {}; src.forEach((s) => { parts[s.key] = s.byFhr.get(f); });
      return { fhr: f, step_h: f - (i ? common[i - 1] : 0), valid: parts[src[0].key].valid, parts,
        n: src.reduce((t, s) => t + (parts[s.key].n ?? s.man.members.length), 0) };
    });
  }
  S.byFhr = new Map(S.steps.map((st) => [st.fhr, st]));
}
const nearestIdx = (fhr) => { let best = S.steps.length - 1, d = Infinity; S.steps.forEach((s, i) => { const x = Math.abs(s.fhr - fhr); if (x < d) { d = x; best = i; } }); return best; };

/** Switch the selected model: returns false (and says why) when it has no data yet. */
function setFamily(fam, quiet) {
  const pick = pickSources(fam);
  if (!pick.src.length) { if (!quiet) toast(pick.note || 'No data for that model yet'); return false; }
  const fhr = S.steps.length ? S.steps[S.idx].fhr : 24;
  S.fam = fam; S.note = pick.note; applySources(pick.src);
  S.idx = nearestIdx(fhr);
  if (!S.products.includes(S.product)) S.product = 'qpf';
  fixDay();
  if (!S.hasSwe) S.slr = 'kuchera';
  FIELD_MEMO.clear(); S.cur = null;
  return true;
}

/* ============================================================================= UI wiring */
function toast(msg) { const t = $('#toast'); t.textContent = msg; t.hidden = false; setTimeout(() => (t.hidden = true), 6000); }

function statOptions() {
  const opts = [];
  const g = (label, items) => opts.push(`<optgroup label="${label}">${items.map(([v, t]) => `<option value="${v}">${t}</option>`).join('')}</optgroup>`);
  if (S.product === 'ptype') {
    g('Ensemble', [['ens', 'Most likely type (QPF-weighted)']]);
    g('Probability of type', PTYPE.names.map((n, t) => [`prob:${t}`, `${PTYPE.long[t]} (% of ${S.src.length > 1 && S.wmode === 'model' ? 'weighted ' : ''}members)`]));
  } else {
    g('Ensemble statistics', [['min', 'Minimum'], ['p5', '5th percentile'], ['p10', '10th percentile'], ['p25', '25th percentile'], ['p50', '50th percentile (median)'],
      ['mean', 'Mean'], ['p75', '75th percentile'], ['p90', '90th percentile'], ['p95', '95th percentile'], ['max', 'Maximum']]);
    if (tempMode()) g('Probability', TempDays.PROB[S.product].map((pr) => [`prob:${pr[0]}:${pr[1]}`, `Chance ${TempDays.probLabel(pr)}`]));
    else {
      g('Probability-matched mean', [['pmm', 'PMM (whole map)'], ['lpmm', 'LPMM (local neighborhood)']]);
      g('Probability of at least', PROB_THRESH[S.product].map((t) => [`prob:${t}`, `${t}″`]));
    }
  }
  S.src.forEach((s, j) => {
    const here = EnsStats.membersFromMask(tempMode() ? currentDay().entries[j].mask : S.steps[S.idx].parts[s.key].mask, s.man.members.length);
    g(S.src.length > 1 ? `Single member \u2014 ${FAMS[s.key].label}` : 'Single member', here.map((k) => [`m${S.off[s.key] + k}`, s.man.members[k] === 'c00' ? 'Control (c00)' : `Member ${s.man.members[k]}`]));
  });
  return opts.join('');
}

function syncControls() {
  $('#stat').innerHTML = statOptions();
  const valid = [...$('#stat').options].map((o) => o.value);
  if (!valid.includes(S.stat)) S.stat = S.product === 'ptype' ? 'ens' : 'mean';
  $('#stat').value = S.stat;
  $('#lpmm-wrap').hidden = S.stat !== 'lpmm'; $('#lpmm-r').value = String(S.lpmmR);
  $('#slr-wrap').hidden = !(S.product === 'snow' && S.hasSwe); $('#slr').value = S.slr;
  document.querySelectorAll('#fam button').forEach((b) => b.classList.toggle('on', b.dataset.v === S.fam));
  $('#wmode-wrap').hidden = S.src.length < 2; $('#wmode').value = S.wmode;
  $('#wmode-help').textContent = S.wmode === 'model' ? 'Each model carries the same total weight, however many members it has.' : 'Every member counts the same, so a model with more members counts more.';
  $('#fam-note').innerHTML = S.src.map((s) => `<span><i class="sw" style="background:${FAMS[s.key].color}"></i>${FAMS[s.key].long} ${hhZ(s.man.init)} \u00b7 ${s.man.members.length} members</span>`).join('') + (S.note ? `<em>${S.note}</em>` : '');
  renderProducts();
  const help = { pmm: 'Ensemble-mean pattern with amplitudes re-assigned from all members pooled over the whole map, so peaks are not smoothed away.',
    lpmm: `Same matching done separately inside a \u00b1${lpmmMi()}-mile neighborhood around each point: keeps local detail and avoids mixing different weather regimes.` }[S.stat];
  $('#stat-help').textContent = help || ''; $('#stat-help').hidden = !help;
  document.querySelectorAll('#product button').forEach((b) => b.classList.toggle('on', b.dataset.v === S.product));
  document.querySelectorAll('#period button').forEach((b) => { b.classList.toggle('on', !tempMode() && String(b.dataset.v) === String(S.period)); b.disabled = tempMode() || (S.product === 'ptype' && b.dataset.v !== 'step'); });
  $('#period-note').hidden = S.product !== 'ptype' && !tempMode(); $('#period-note').textContent = tempMode() ? '(not used: one value per day)' : '(type is per output step)';
  location.hash = `m=${S.fam}${S.src.length > 1 && S.wmode !== 'model' ? `&w=${S.wmode}` : ''}&p=${S.product}${tempMode() ? '' : `&per=${S.period}`}&s=${S.stat}${S.stat === 'lpmm' ? `&r=${S.lpmmR}` : ''}${S.product === 'snow' && S.hasSwe && S.slr !== 'kuchera' ? `&slr=${S.slr}` : ''}${tempMode() ? `&d=${S.day.replace(/-/g, '')}` : `&f=${S.steps[S.idx].fhr}`}${S.stn ? '&st=1' : ''}`;
}

function renderProducts() {
  const have = new Set(S.products), el = $('#product'), html = [['qpf', 'QPF'], ['snow', 'Snowfall'], ['zr', 'Freezing rain'], ['pl', 'Sleet'], ['ptype', 'Precip type'], ['tmax', 'Daily MaxT'], ['tmin', 'Daily MinT']].filter(([v]) => have.has(v)).map(([v, t]) => `<button data-v="${v}">${t}</button>`).join('');
  if (el.dataset.have !== html) { el.innerHTML = html; el.dataset.have = html; }
}

/** Day labels under the slider (they depend on which hours the selected model(s) have). */
function buildDays() {
  const days = $('#days'); days.innerHTML = ''; let prev = '';
  if (tempMode()) {                                // one tick per day (labelled by the UTC date the window starts on)
    const d = S.tdays[S.product], every = innerWidth < 760 && d.length > 6 ? 2 : 1;
    d.forEach((x, i) => {
      if (i % every) return;
      const sp = document.createElement('span'), fx = d.length > 1 ? i / (d.length - 1) : 0;
      sp.textContent = new Date(`${x.date}T12:00:00Z`).toLocaleDateString([], { weekday: 'short', day: 'numeric', timeZone: 'UTC' });
      sp.style.left = `${fx * 100}%`; if (fx > 0.93) sp.style.transform = 'translateX(-100%)'; days.appendChild(sp);
    });
    $('#slider').max = Math.max(0, d.length - 1); return;
  }
  const every = innerWidth < 760 && S.steps.length > 40 ? 2 : 1; let nd = 0;
  S.steps.forEach((s, i) => { const d = new Date(s.valid), key = d.toDateString(); if (key !== prev) { prev = key; if (nd++ % every) return; const sp = document.createElement('span'); sp.textContent = d.toLocaleDateString([], { weekday: 'short', day: 'numeric' }); const fx = i / (S.steps.length - 1); sp.style.left = `${fx * 100}%`; if (fx > 0.93) sp.style.transform = 'translateX(-100%)'; days.appendChild(sp); } });
  $('#slider').max = S.steps.length - 1;
}

function buildUI() {
  const famBtn = [...(pickSources('blend').src.length ? [['blend', 'Blend']] : []), ...loadedFams().map((k) => [k, FAMS[k].label])];
  $('#fam').innerHTML = famBtn.map(([v, t]) => `<button data-v="${v}">${t}</button>`).join('');
  $('#fam').onclick = (e) => { const b = e.target.closest('button'); if (!b || b.dataset.v === S.fam) return; if (setFamily(b.dataset.v)) { buildDays(); syncControls(); refresh(); } };
  $('#wmode').onchange = (e) => { S.wmode = e.target.value; FIELD_MEMO.clear(); syncControls(); refresh(); };
  $('#period').innerHTML = PERIODS.map(([v, t]) => `<button data-v="${v}">${t}</button>`).join('');
  $('#product').onclick = (e) => { const b = e.target.closest('button'); if (!b) return; S.product = b.dataset.v; if (S.product === 'ptype') S.period = 'step'; else if (S.period === 'step') S.period = 24; fixDay(); buildDays(); syncControls(); refresh(); };
  $('#period').onclick = (e) => { const b = e.target.closest('button'); if (!b || b.disabled) return; S.period = b.dataset.v === 'step' ? 'step' : +b.dataset.v; syncControls(); refresh(); };
  $('#stat').onchange = (e) => { S.stat = e.target.value; syncControls(); refresh(); };
  $('#lpmm-r').onchange = (e) => { S.lpmmR = +e.target.value; syncControls(); refresh(); };
  $('#slr').onchange = (e) => { S.slr = e.target.value; syncControls(); refresh(); };
  $('#opacity').oninput = (e) => wx.setOpacity(+e.target.value / 100);
  const vis = (id, layers) => ($(id).onchange = (e) => layers.forEach((l) => map.setLayoutProperty(l, 'visibility', e.target.checked ? 'visible' : 'none')));
  vis('#tg-counties', ['counties-line']); vis('#tg-states', ['admin1-line']); vis('#tg-rivers', ['rivers']); vis('#tg-grat', ['graticule']);
  $('#tg-labels').onchange = (e) => { labels.show.places = labels.show.admin1 = labels.show.water = labels.show.countries = e.target.checked; labels.q(); };
  $('#tg-stations').onchange = (e) => { S.stn = e.target.checked; labels.show.stations = S.stn; updateStations(); syncControls(); };
  $('#tg-grat').addEventListener('change', (e) => { labels.show.graticule = e.target.checked; labels.q(); });
  $('#views').onclick = (e) => {
    const b = e.target.closest('button'); if (!b) return;
    const v = { ne: [[-74.5, 41], [-66.5, 47.8]], conus: [[-125.5, 24], [-66, 50]], na: [[-135, 20], [-52, 58]], se: [[-92, 24.5], [-75, 37.5]], mw: [[-98, 36], [-80, 49]] }[b.dataset.v];
    map.fitBounds(v, { padding: 20, duration: 600 });
  };
  const sl = $('#slider');
  sl.oninput = () => { if (tempMode()) S.day = S.tdays[S.product][+sl.value].date; else S.idx = +sl.value; syncControls(); refresh(); };
  $('#prev').onclick = () => step(-1); $('#next').onclick = () => step(1);
  $('#play').onclick = togglePlay;
  addEventListener('keydown', (e) => { if (/INPUT|SELECT|TEXTAREA/.test(document.activeElement.tagName) && document.activeElement.type !== 'range') return; if (e.key === 'ArrowRight') step(1); else if (e.key === 'ArrowLeft') step(-1); else if (e.key === ' ') { e.preventDefault(); togglePlay(); } });
  $('#menu').onclick = () => $('#panel').classList.toggle('open');
  if (innerWidth > 760) document.body.classList.add('panel-open');
  $('#legend-toggle').onclick = () => $('#legend').classList.toggle('min');
  buildDays();
}
const navLen = () => (tempMode() ? S.tdays[S.product].length : S.steps.length), navPos = () => (tempMode() ? dayPos() : S.idx);
function navSet(j) { if (tempMode()) S.day = S.tdays[S.product][j].date; else S.idx = j; }
function step(d) { const j = Math.max(0, Math.min(navLen() - 1, navPos() + d)); if (j === navPos()) { if (S.playing) togglePlay(); return; } navSet(j); syncControls(); refresh(); }
function togglePlay() {
  S.playing = !S.playing; $('#play').textContent = S.playing ? '❚❚' : '▶'; $('#play').setAttribute('aria-label', S.playing ? 'Pause' : 'Play');
  clearInterval(S.timer);
  if (S.playing) { if (navPos() >= navLen() - 1) navSet(0); S.timer = setInterval(() => { if (!$('#busy').hidden) return; step(1); }, 650); }
}

/* ============================================================================= boot */
async function boot() {
  map = createMap('map');
  const getJson = (u) => fetch(`${u}?t=${Date.now()}`, { cache: 'no-store' }).then((r) => (r.ok ? r.json() : null)).catch(() => null);
  let [gl, gs, il, pl, lab] = await Promise.all([getJson('data/gefs/latest.json'), getJson('data/gefs/synoptic.json'), getJson('data/ifs/latest.json'),
    getJson('data/geps/latest.json'), fetch('basemap_labels.json').then((r) => r.json())]);
  if (!gl) { const old = await getJson('data/latest.json'); if (old) { old.legacy = true; gl = old; } }   // single-model layout written by the older workflow
  S.sets = { gefs: { latest: gl, syn: gs }, ifs: { latest: il }, geps: { latest: pl } };
  if (!loadedFams().length) throw new Error('No model data found yet \u2014 the pipeline has not published a run');
  lab.places.sort((a, b) => b[3] - a[3]); labelData = lab;
  const hs = Object.fromEntries(new URLSearchParams(location.hash.slice(1)));
  if (hs.w === 'member') S.wmode = 'member';
  const first = [hs.m, 'blend', ...FAM_KEYS].find((f) => f && setFamily(f, true));           // the requested view, else the best one that has data
  const want = +hs.f || 24; S.idx = nearestIdx(want);
  S.product = ['qpf', 'snow', 'zr', 'pl', 'ptype', 'tmax', 'tmin'].includes(hs.p) && S.products.includes(hs.p) ? hs.p : 'qpf';
  if (/^\d{8}$/.test(hs.d || '')) S.day = `${hs.d.slice(0, 4)}-${hs.d.slice(4, 6)}-${hs.d.slice(6)}`;
  fixDay();
  S.period = S.product === 'ptype' ? 'step' : hs.per === 'step' ? 'step' : [6, 12, 24, 48].includes(+hs.per) ? +hs.per : 24;
  if (hs.s) S.stat = hs.s === 'median' ? 'p50' : hs.s;
  if ([2, 4, 6, 8].includes(+hs.r)) S.lpmmR = +hs.r;
  if (S.hasSwe && ['10', '12', '13'].includes(hs.slr)) S.slr = hs.slr;
  S.stn = hs.st === '1';
  const stationsReady = loadStations();
  await new Promise((res) => (map.loaded() ? res() : map.once('load', res)));
  wx = new WxLayer({ opacity: S.opacity }); wx.setDomain(S.domain);
  map.addLayer(wx, 'counties-hit');
  labels = new LabelLayer(map, $('#labels'), lab);
  const grat = () => { labels.step = updateGraticule(map); };
  map.on('moveend', grat); grat();
  map.addControl(new maplibregl.NavigationControl({ showCompass: false }), 'top-right');
  map.addControl(new maplibregl.ScaleControl({ unit: 'imperial', maxWidth: 110 }), 'bottom-left');
  map.on('mousemove', (e) => { if (!recentTouch()) showTip(e); });
  map.on('mouseout', () => { if (!recentTouch()) hideTip(); });
  map.on('click', (e) => { if (recentTouch() && S.pin) unpin(); pin(e); });
  map.on('dragstart', () => hideTip());
  buildUI(); $('#tg-stations').checked = S.stn; labels.show.stations = S.stn; syncControls();
  stationsReady.then(() => updateStations());
  $('#splash').hidden = true;
  $('#t-gen').textContent = `Built ${new Date(S.gen).toLocaleString([], { month: 'short', day: 'numeric', hour: 'numeric', minute: '2-digit', timeZoneName: 'short' })}${S.complete ? '' : ' \u00b7 run still filling in'}`;
  refresh();
}
boot().catch((e) => { console.error(e); $('#splash').hidden = false; $('#splash-msg').textContent = e.message; });
