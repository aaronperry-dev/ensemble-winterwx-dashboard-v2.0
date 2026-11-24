/* GEFS ensemble winter-weather viewer.
 * Loads small gzipped binary files (cumulative QPF/snow per member, p-type probabilities per member),
 * derives any window total / ensemble statistic in the browser, and paints it with a WebGL layer.
 */
'use strict';

const $ = (s, el = document) => el.querySelector(s);
const isMem = (s) => /^m\d+$/.test(s);
const MM_PER_IN = 25.4, CM_PER_IN = 2.54;
const TIE_ORDER = [2, 3, 1, 0];            // FZRA > PL > SN > RA  (indices into RA, SN, FZRA, PL)
const PERIODS = [['step', '3 hr'], [6, '6 hr'], [12, '12 hr'], [24, '24 hr'], [48, '48 hr']];
const PROB_THRESH = { qpf: [0.01, 0.10, 0.25, 0.50, 1.00, 2.00], snow: [1, 3, 6, 12, 18, 24],
  zr: [0.01, 0.05, 0.10, 0.25, 0.50, 1.00], pl: [0.01, 0.05, 0.10, 0.25, 0.50, 1.00] };   // inches
/* Accumulation products: array key in the cumulative file, header scale slot, and SI unit. Snowfall is in cm; the rest are mm of liquid. */
const PRODUCTS = { qpf: { name: 'QPF', key: 'q', scale: 's0', si: 'mm' }, snow: { name: 'Snowfall', key: 's', scale: 's1', si: 'cm' },
  zr: { name: 'Freezing rain', key: 'z', scale: 's2', si: 'mm' }, pl: { name: 'Sleet', key: 'p', scale: 's3', si: 'mm' } };
const siPerIn = (kind) => (PRODUCTS[kind].si === 'cm' ? CM_PER_IN : MM_PER_IN);
const PROB_RAMP = {
  qpf: ['#e5f5e0', '#c7e9c0', '#a1d99b', '#74c476', '#41ab5d', '#31a353', '#238b45', '#006d2c', '#00441b'],
  snow: ['#deebf7', '#c6dbef', '#bdd7e7', '#9ecae1', '#6baed6', '#4292c6', '#3182bd', '#08519c', '#082694'],
};
const PROB_EDGES = [10, 20, 30, 40, 50, 60, 70, 80, 90];

const S = {
  man: null, base: '', steps: [], byFhr: new Map(), idx: 0,
  product: 'qpf', period: 24, stat: 'mean', lpmmR: 4, units: 'in', opacity: 0.82,
  cur: null, token: 0, pin: null, playing: false, timer: null,
};
let map, wx, labels, labelData, pinMarker;

/* ============================================================================= data loading */
const CACHE = new Map();
function cached(url, make) {
  if (CACHE.has(url)) { const p = CACHE.get(url); CACHE.delete(url); CACHE.set(url, p); return p; }
  const p = make().catch((e) => { CACHE.delete(url); throw e; });
  CACHE.set(url, p);
  while (CACHE.size > 36) CACHE.delete(CACHE.keys().next().value);
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
      s2: dv.getFloat32(40, true), s3: dv.getFloat32(44, true) };
    const n = h.nmem * h.ny * h.nx;
    if (h.kind === 1) {
      const o = { h }; ['q', 's', 'z', 'p'].slice(0, h.narr).forEach((k, i) => { o[k] = new Uint16Array(buf, 64 + 2 * n * i, n); });
      return o;
    }
    return { h, p: new Uint8Array(buf, 64, 4 * n) };
  });
}
const loadCum = (fhr) => fetchBin(S.base + S.byFhr.get(fhr).cum);
const loadPtype = (fhr) => fetchBin(S.base + S.byFhr.get(fhr).ptype);

/* ============================================================================= windows & statistics */
function windowFor(i, period) {
  const st = S.steps[i], end = st.fhr;
  let start = period === 'step' ? end - st.step_h : end - period, since = false;
  if (start < 0 || (start > 0 && !S.byFhr.has(start))) { since = true; start = 0; }
  return { end, start, since, hours: end - start };
}

/** Member totals over the window, raw units (mm for QPF, cm for snow). Float32Array(M*N). */
async function windowValues(kind, i, period) {
  const w = windowFor(i, period), key = PRODUCTS[kind].key;
  const [a, b] = await Promise.all([loadCum(w.end), w.start > 0 ? loadCum(w.start) : null]);
  if (!a[key]) throw new Error(`${PRODUCTS[kind].name} is not in this run's files yet (the next pipeline run adds it)`);
  const scale = a.h[PRODUCTS[kind].scale], n = a[key].length, out = new Float32Array(n), A = a[key];
  if (b) { const B = b[key]; for (let k = 0; k < n; k++) out[k] = Math.max(0, A[k] - B[k]) * scale; }
  else for (let k = 0; k < n; k++) out[k] = A[k] * scale;
  return { data: out, w, h: a.h };
}

const toDisp = (kind, v) => v / siPerIn(kind);   // raw (mm or cm) -> inches

const pct = EnsStats.percentile;
const FIELD_MEMO = new Map();          // PMM / LPMM fields are the expensive ones: keep the last few

/** One display-unit (inches) field from the member totals of a window. */
function reduceStat(kind, win, stat) {
  const { nmem: M, ny, nx } = win.h, N = ny * nx, data = win.data, c = 1 / siPerIn(kind);
  if (isMem(stat)) { const out = new Float32Array(N), o = +stat.slice(1) * N; for (let i = 0; i < N; i++) out[i] = data[o + i] * c; return out; }
  if (stat.startsWith('prob:')) {
    const out = new Float32Array(N), thr = (+stat.slice(5)) / c;                 // inches -> raw
    for (let i = 0; i < N; i++) { let k = 0; for (let m = 0; m < M; m++) if (data[m * N + i] >= thr) k++; out[i] = (100 * k) / M; }
    return out;
  }
  const slow = stat === 'pmm' || stat === 'lpmm';
  const key = slow ? [kind, S.man.init, win.w.start, win.w.end, stat, S.lpmmR].join('|') : null;
  if (key && FIELD_MEMO.has(key)) return FIELD_MEMO.get(key);
  const out = EnsStats.compute(data, M, ny, nx, stat, { radius: S.lpmmR });
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
  const [pt, qw] = await Promise.all([loadPtype(st.fhr), windowValues('qpf', S.idx, 'step')]);
  return { p: pt.p, q: qw.data, h: pt.h, w: qw.w };
}

function buildPtype(d, stat) {
  const { p, q, h, w } = d, M = h.nmem, N = h.ny * h.nx, rateDiv = w.hours || 3;
  const field = new Float32Array(N), types = new Float32Array(N), tot = M * N;
  const share = (m, i, t) => { let s = 0; for (let k = 0; k < 4; k++) s += p[k * tot + m * N + i]; return s ? p[t * tot + m * N + i] / s : 0; };
  for (let i = 0; i < N; i++) {
    if (isMem(stat)) {
      const m = +stat.slice(1), t = domType(p, N, m, i);
      if (t >= 0) { field[i] = dbzOf(q[m * N + i] / rateDiv, t); types[i] = t; }
    } else if (stat.startsWith('prob:')) {
      const want = +stat.slice(5); let k = 0;
      for (let m = 0; m < M; m++) if (domType(p, N, m, i) === want) k++;
      field[i] = (100 * k) / M; types[i] = want;
    } else {                                              // 'ens': QPF-weighted dominant type, mean-QPF intensity
      const wsum = [0, 0, 0, 0]; let qs = 0;
      for (let m = 0; m < M; m++) { const qq = q[m * N + i]; qs += qq; if (qq > 0) for (let t = 0; t < 4; t++) wsum[t] += qq * share(m, i, t); }
      let best = -1, bv = 0; for (const t of TIE_ORDER) if (wsum[t] > bv) { bv = wsum[t]; best = t; }
      if (best >= 0) { field[i] = dbzOf(qs / M / rateDiv, best); types[i] = best; }
    }
  }
  return { field, types };
}

/* ============================================================================= rendering */
function classesFor() {
  if (S.stat.startsWith('prob:')) {
    const k = PROB_RAMP[S.product] || ptypeProbColors(S.product === 'zr' ? 2 : S.product === 'pl' ? 3 : +S.stat.slice(5));
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
      if (S.stat.startsWith('prob:')) { const c = classesFor(); wx.setClasses(c.edges, c.colors); }
      else wx.setPtypeLut(PTYPE.ramps, PTYPE.names, PTYPE.dbzMin, PTYPE.dbzStep);
      wx.setField(field, types); S.cur = { kind: 'ptype', d, field };
    } else {
      const win = await windowValues(S.product, S.idx, S.period); if (tok !== S.token) return;
      const field = reduceStat(S.product, win, S.stat);
      const c = classesFor(); wx.setClasses(c.edges, c.colors); wx.setField(field);
      S.cur = { kind: S.product, win, field };
    }
    updateLegend(); updateTitle(); updateTip(); prefetch();
  } catch (e) { console.error(e); toast(`Could not load data: ${e.message}`); }
  finally { if (tok === S.token) $('#busy').hidden = true; }
}

function prefetch() {
  [S.idx + 1, S.idx - 1].forEach((j) => {
    if (j < 0 || j >= S.steps.length) return;
    const f = S.steps[j].fhr; loadCum(f).catch(() => {});
    if (S.product === 'ptype') loadPtype(f).catch(() => {});
  });
}

/* ============================================================================= labels / titles */
const tz = Intl.DateTimeFormat().resolvedOptions().timeZone;
const fmtLocal = (d) => d.toLocaleString([], { weekday: 'short', month: 'short', day: 'numeric', hour: 'numeric', timeZoneName: 'short' });
const fmtUtc = (d) => `${d.toISOString().slice(5, 10).replace('-', '/')} ${d.toISOString().slice(11, 13)}Z`;
const pad = (n, w = 3) => String(n).padStart(w, '0');

function productName() {
  if (S.product === 'ptype') return 'Precipitation type';
  const nm = PRODUCTS[S.product].name;
  const w = windowFor(S.idx, S.period);
  return `${w.since ? `${w.hours}-hr (since init)` : `${w.hours}-hr`} ${nm}`;
}
function statName() {
  const s = S.stat;
  if (s === 'mean') return 'Ensemble mean'; if (s === 'max') return 'Ensemble maximum'; if (s === 'min') return 'Ensemble minimum';
  if (s === 'pmm') return 'Probability-matched mean (PMM)';
  if (s === 'lpmm') return `Local prob.-matched mean (LPMM \u00b1${S.lpmmR / 2}\u00b0)`;
  if (EnsStats.isPct(s)) return s === 'p50' ? '50th percentile (median)' : `${s.slice(1)}th percentile`;
  if (s === 'ens') return 'Most likely type (ensemble)';
  if (s.startsWith('prob:')) {
    const v = +s.slice(5);
    return S.product === 'ptype' ? `Probability of ${PTYPE.long[v].toLowerCase()}` : `Probability ≥ ${v}″`;
  }
  const m = +s.slice(1); return m === 0 ? 'Control (c00)' : `Member p${pad(m, 2)}`;
}

function updateTitle() {
  const st = S.steps[S.idx], valid = new Date(st.valid), init = new Date(S.man.init);
  const w = windowFor(S.idx, S.product === 'ptype' ? 'step' : S.period);
  const startD = new Date(init.getTime() + w.start * 3600e3);
  $('#t-main').textContent = `${statName()} · ${productName()}`;
  $('#t-valid').textContent = S.product === 'ptype' || w.hours <= st.step_h
    ? `Valid ${fmtLocal(valid)} (${fmtUtc(valid)})`
    : `${fmtLocal(startD)} → ${fmtLocal(valid)}`;
  $('#t-init').textContent = `GEFS ${init.toISOString().slice(0, 10)} ${init.toISOString().slice(11, 13)}Z run · F${pad(st.fhr)} · ${S.man.members.length} members`;
  $('#fh').textContent = `F${pad(st.fhr)}`;
  $('#fh-time').textContent = `${fmtLocal(valid)}`;
  $('#fh-utc').textContent = fmtUtc(valid);
  $('#slider').value = S.idx;
}

function fmtEdge(v, kind) {
  if (S.units === 'in') return null;
  const x = v * siPerIn(kind);
  return x < 10 ? x.toFixed(x < 1 ? 2 : 1) : Math.round(x);
}
const fmtEdge2b = (t, kind) => `${fmtEdge(t, kind) ?? t} ${PRODUCTS[kind].si}`;

function updateLegend() {
  const el = $('#legend-body'); el.innerHTML = '';
  const hdr = $('#legend-title');
  if (S.product === 'ptype' && !S.stat.startsWith('prob:')) {
    hdr.textContent = 'Type color ramps — 3-hr QPF as dBZ';
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
    note.textContent = 'dBZ ≈ rate from 3-hr QPF (rain/ice: Z=200R^1.6, snow: Z=1780R^2.21)'; el.appendChild(note); return;
  }
  const c = classesFor();
  const kind = S.product === 'ptype' ? 'ptype' : S.product;
  hdr.textContent = c.prob ? `Probability (% of ${S.man.members.length} members)` :
    `${c.title} (${S.units === 'in' ? 'in' : PRODUCTS[S.product].si})`;
  for (let i = c.colors.length - 1; i >= 0; i--) {
    const row = document.createElement('div'); row.className = 'lg-row';
    let label = c.labels[i];
    if (!c.prob && S.units !== 'in') {
      const lo = fmtEdge(c.edges[i], kind), hi = i + 1 < c.edges.length ? fmtEdge(c.edges[i + 1], kind) : null;
      label = hi === null ? `${lo}+` : S.product === 'snow' && i === 0 ? `< ${hi}` : `${lo}–${hi}`;
    }
    row.innerHTML = `<i style="background:${c.colors[i]}"></i><span>${label}</span>`; el.appendChild(row);
  }
}

/* ============================================================================= tooltip */
const tip = () => $('#tip');
let touchAt = 0;
addEventListener('touchstart', () => { touchAt = Date.now(); }, { capture: true, passive: true });
const recentTouch = () => Date.now() - touchAt < 1200;   // ignore emulated mouse events that follow a tap
let lastLL = null, raf = 0;

function cellValues(lngLat, data, M, N, h) {
  const d = S.man.domain, gx = (lngLat.lng - d.lon_west) / d.dlon, gy = (d.lat_north - lngLat.lat) / d.dlat;
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

const fmtVal = (kind, inches) => {
  if (S.units === 'in') return `${inches < 0.1 && inches > 0 ? inches.toFixed(3) : inches.toFixed(inches < 10 ? 2 : 1)}″`;
  const x = inches * siPerIn(kind);
  return `${x < 10 ? x.toFixed(x < 1 ? 2 : 1) : Math.round(x)} ${PRODUCTS[kind].si}`;
};

function stripSVG(vals, mean, med, kind, selected) {
  const W = 268, H = 40, pad = 6, n = vals.length, mx = Math.max(...vals, 1e-9) * 1.04;
  const sx = (v) => pad + (v / mx) * (W - 2 * pad), bins = new Map(), dots = [];
  vals.forEach((v, m) => { const b = Math.round(sx(v) / 5); const k = bins.get(b) || 0; bins.set(b, k + 1); dots.push([sx(v), H - 8 - k * 5.2, m]); });
  let s = `<svg viewBox="0 0 ${W} ${H}" width="100%" role="img" aria-label="Ensemble member values">`;
  s += `<line x1="${pad}" x2="${W - pad}" y1="${H - 3}" y2="${H - 3}" stroke="#9aa3ab" stroke-width="1"/>`;
  dots.forEach(([x, y, m]) => { s += `<circle cx="${x.toFixed(1)}" cy="${Math.max(4, y).toFixed(1)}" r="2.3" fill="${m === selected ? '#111' : '#5b6670'}" ${m === selected ? 'stroke="#fff" stroke-width="1"' : ''}/>`; });
  s += `<line x1="${sx(mean)}" x2="${sx(mean)}" y1="2" y2="${H - 3}" stroke="#c2410c" stroke-width="1.5"/>`;
  s += `<line x1="${sx(med)}" x2="${sx(med)}" y1="2" y2="${H - 3}" stroke="#1d4ed8" stroke-width="1.5" stroke-dasharray="2 2"/>`;
  return s + '</svg>';
}

function tipHTML(ll, pt) {
  const cur = S.cur; if (!cur) return '';
  const nm = placeName(pt);
  const head = `<div class="tip-h"><b>${nm.line1 || 'Open water / outside boundaries'}</b>${nm.line2 ? `<span>${nm.line2}</span>` : ''}
    <small>${Math.abs(ll.lat).toFixed(2)}°${ll.lat >= 0 ? 'N' : 'S'}, ${Math.abs(ll.lng).toFixed(2)}°${ll.lng >= 0 ? 'E' : 'W'}</small></div>`;
  if (cur.kind === 'ptype') return head + ptypeTip(ll, cur.d);
  const { win } = cur, M = win.h.nmem, N = win.h.ny * win.h.nx, cell = cellValues(ll, win.data, M, N, win.h);
  if (!cell) return head + '<div class="tip-b muted">Outside the GEFS domain.</div>';
  const vals = Array.from(cell.v, (x) => toDisp(S.product, x)), sorted = Float32Array.from(vals).sort();
  const mean = vals.reduce((a, b) => a + b, 0) / M, med = pct(sorted, M, 0.5);
  const kind = S.product, f = (x) => fmtVal(kind, x);
  const sel = isMem(S.stat) ? +S.stat.slice(1) : -1;
  const val = { min: sorted[0], max: sorted[M - 1], mean };
  EnsStats.PCTS.forEach((p) => { val[`p${p}`] = pct(sorted, M, p / 100); });
  const nm2 = (k) => (k === 'mean' ? 'Mean' : k === 'min' ? 'Min' : k === 'max' ? 'Max' : k.toUpperCase());
  const ord = EnsStats.TIP_ORDER, half = ord.length / 2;
  let h = head + `<div class="tip-t">${productName()} \u00b7 ${$('#t-valid').textContent}</div>`;
  if ((S.stat === 'pmm' || S.stat === 'lpmm') && cur.field) {
    const fv = cellValues(ll, cur.field, 1, N, win.h);
    if (fv) h += `<div class="tip-map"><span>${S.stat === 'pmm' ? 'PMM' : `LPMM \u00b1${S.lpmmR / 2}\u00b0`}</span><b>${f(fv.v[0])}</b></div>`;
  }
  h += '<div class="tip-grid">';
  for (let r = 0; r < half; r++) h += `<span>${nm2(ord[r])}</span><b>${f(val[ord[r]])}</b><span>${nm2(ord[r + half])}</span><b>${f(val[ord[r + half]])}</b>`;
  h += '</div><div class="tip-s">Chance of at least</div><div class="tip-p">';
  PROB_THRESH[kind].forEach((t) => {
    const p = Math.round((100 * vals.filter((v) => v >= t).length) / M);
    const lab = S.units === 'in' ? `${t}″` : `${fmtEdge2b(t, kind)}`;
    h += `<span>${lab}</span><i><u style="width:${p}%"></u></i><b>${p}%</b>`;
  });
  h += `</div><div class="tip-s">All ${M} members <span class="key"><i class="k1"></i>mean <i class="k2"></i>median</span></div>${stripSVG(vals, mean, med, kind, sel)}`;
  h += `<div class="tip-ax"><span>0</span><span>${f(Math.max(...vals))}</span></div>`;
  return h;
}

function ptypeTip(ll, d) {
  const { p, q, h } = d, M = h.nmem, N = h.ny * h.nx, c = cellValues(ll, q, M, N, h);
  if (!c) return '<div class="tip-b muted">Outside the GEFS domain.</div>';
  const i = c.nearest, cnt = [0, 0, 0, 0], avg = [0, 0, 0, 0]; let wet = 0, qs = 0;
  for (let m = 0; m < M; m++) {
    const t = domType(p, N, m, i); if (t < 0) continue; cnt[t]++; wet++;
    for (let k = 0; k < 4; k++) avg[k] += p[k * M * N + m * N + i];
  }
  for (let m = 0; m < M; m++) qs += c.v[m];
  const qm = qs / M, order = [1, 0, 2, 3];
  let h2 = `<div class="tip-t">Precip type · ${$('#t-valid').textContent}</div>`;
  h2 += `<div class="tip-grid"><span>Members with precip</span><b>${wet} / ${M}</b><span>Mean 3-hr QPF</span><b>${fmtVal('qpf', qm / MM_PER_IN)}</b></div>`;
  h2 += '<div class="tip-s">Dominant type (members)</div><div class="tip-p">';
  order.forEach((t) => { h2 += `<span><i class="sw" style="background:${PTYPE.ramps[PTYPE.names[t]][18]}"></i>${PTYPE.long[t]}</span><i><u style="width:${(100 * cnt[t]) / M}%"></u></i><b>${cnt[t]}</b>`; });
  h2 += '</div><div class="tip-s">Mean Bourgouin probability (wet members)</div><div class="tip-grid">';
  order.forEach((t) => { h2 += `<span>${PTYPE.long[t]}</span><b>${wet ? Math.round(avg[t] / wet) : 0}%</b>`; });
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

/* ============================================================================= UI wiring */
function toast(msg) { const t = $('#toast'); t.textContent = msg; t.hidden = false; setTimeout(() => (t.hidden = true), 6000); }

function statOptions() {
  const M = S.man.members, opts = [];
  const g = (label, items) => opts.push(`<optgroup label="${label}">${items.map(([v, t]) => `<option value="${v}">${t}</option>`).join('')}</optgroup>`);
  if (S.product === 'ptype') {
    g('Ensemble', [['ens', 'Most likely type (QPF-weighted)']]);
    g('Probability of type', PTYPE.names.map((n, t) => [`prob:${t}`, `${PTYPE.long[t]} (% of members)`]));
  } else {
    g('Ensemble statistics', [['min', 'Minimum'], ['p5', '5th percentile'], ['p10', '10th percentile'], ['p25', '25th percentile'], ['p50', '50th percentile (median)'],
      ['mean', 'Mean'], ['p75', '75th percentile'], ['p90', '90th percentile'], ['p95', '95th percentile'], ['max', 'Maximum']]);
    g('Probability-matched mean', [['pmm', 'PMM (whole map)'], ['lpmm', 'LPMM (local neighborhood)']]);
    g('Probability of at least', PROB_THRESH[S.product].map((t) => [`prob:${t}`, `${t}″ (${fmtEdge2(t)})`]));
  }
  g('Single member', M.map((m, k) => [`m${k}`, k === 0 ? 'Control (c00)' : `Member ${m}`]));
  return opts.join('');
}
const fmtEdge2 = (t) => { const x = t * siPerIn(S.product); return `${x < 10 ? x.toFixed(x < 1 ? 2 : 1) : Math.round(x)} ${PRODUCTS[S.product].si}`; };

function syncControls() {
  $('#stat').innerHTML = statOptions();
  const valid = [...$('#stat').options].map((o) => o.value);
  if (!valid.includes(S.stat)) S.stat = S.product === 'ptype' ? 'ens' : 'mean';
  $('#stat').value = S.stat;
  $('#lpmm-wrap').hidden = S.stat !== 'lpmm'; $('#lpmm-r').value = String(S.lpmmR);
  const help = { pmm: 'Ensemble-mean pattern with amplitudes re-assigned from all members pooled over the whole map, so peaks are not smoothed away.',
    lpmm: `Same matching done separately inside a \u00b1${S.lpmmR / 2}\u00b0 neighborhood around each point: keeps local detail and avoids mixing different weather regimes.` }[S.stat];
  $('#stat-help').textContent = help || ''; $('#stat-help').hidden = !help;
  document.querySelectorAll('#product button').forEach((b) => b.classList.toggle('on', b.dataset.v === S.product));
  document.querySelectorAll('#period button').forEach((b) => { b.classList.toggle('on', String(b.dataset.v) === String(S.period)); b.disabled = S.product === 'ptype' && b.dataset.v !== 'step'; });
  $('#period-note').hidden = S.product !== 'ptype';
  location.hash = `p=${S.product}&per=${S.period}&s=${S.stat}${S.stat === 'lpmm' ? `&r=${S.lpmmR}` : ''}&f=${S.steps[S.idx].fhr}`;
}

function buildUI() {
  const have = new Set(S.man.products || ['qpf', 'snow', 'ptype']);
  $('#product').innerHTML = [['qpf', 'QPF'], ['snow', 'Snowfall'], ['zr', 'Freezing rain'], ['pl', 'Sleet'], ['ptype', 'Precip type']].filter(([v]) => have.has(v)).map(([v, t]) => `<button data-v="${v}">${t}</button>`).join('');
  $('#period').innerHTML = PERIODS.map(([v, t]) => `<button data-v="${v}">${t}</button>`).join('');
  $('#product').onclick = (e) => { const b = e.target.closest('button'); if (!b) return; S.product = b.dataset.v; if (S.product === 'ptype') S.period = 'step'; else if (S.period === 'step') S.period = 24; syncControls(); refresh(); };
  $('#period').onclick = (e) => { const b = e.target.closest('button'); if (!b || b.disabled) return; S.period = b.dataset.v === 'step' ? 'step' : +b.dataset.v; syncControls(); refresh(); };
  $('#stat').onchange = (e) => { S.stat = e.target.value; syncControls(); refresh(); };
  $('#lpmm-r').onchange = (e) => { S.lpmmR = +e.target.value; syncControls(); refresh(); };
  $('#units').onchange = (e) => { S.units = e.target.value; updateLegend(); updateTip(); statRefresh(); };
  $('#opacity').oninput = (e) => wx.setOpacity(+e.target.value / 100);
  const vis = (id, layers) => ($(id).onchange = (e) => layers.forEach((l) => map.setLayoutProperty(l, 'visibility', e.target.checked ? 'visible' : 'none')));
  vis('#tg-counties', ['counties-line']); vis('#tg-states', ['admin1-line']); vis('#tg-rivers', ['rivers']); vis('#tg-grat', ['graticule']);
  $('#tg-labels').onchange = (e) => { labels.show.places = labels.show.admin1 = labels.show.water = labels.show.countries = e.target.checked; labels.q(); };
  $('#tg-grat').addEventListener('change', (e) => { labels.show.graticule = e.target.checked; labels.q(); });
  $('#views').onclick = (e) => {
    const b = e.target.closest('button'); if (!b) return;
    const v = { ne: [[-74.5, 41], [-66.5, 47.8]], conus: [[-125.5, 24], [-66, 50]], na: [[-135, 20], [-52, 58]], se: [[-92, 24.5], [-75, 37.5]], mw: [[-98, 36], [-80, 49]] }[b.dataset.v];
    map.fitBounds(v, { padding: 20, duration: 600 });
  };
  const sl = $('#slider'); sl.max = S.steps.length - 1;
  sl.oninput = () => { S.idx = +sl.value; syncControls(); refresh(); };
  $('#prev').onclick = () => step(-1); $('#next').onclick = () => step(1);
  $('#play').onclick = togglePlay;
  addEventListener('keydown', (e) => { if (/INPUT|SELECT|TEXTAREA/.test(document.activeElement.tagName) && document.activeElement.type !== 'range') return; if (e.key === 'ArrowRight') step(1); else if (e.key === 'ArrowLeft') step(-1); else if (e.key === ' ') { e.preventDefault(); togglePlay(); } });
  $('#menu').onclick = () => $('#panel').classList.toggle('open');
  if (innerWidth > 760) document.body.classList.add('panel-open');
  $('#legend-toggle').onclick = () => $('#legend').classList.toggle('min');
  // day labels under slider
  const days = $('#days'); days.innerHTML = ''; let prev = '';
  const every = innerWidth < 760 && S.steps.length > 40 ? 2 : 1; let nd = 0;
  S.steps.forEach((s, i) => { const d = new Date(s.valid), key = d.toDateString(); if (key !== prev) { prev = key; if (nd++ % every) return; const sp = document.createElement('span'); sp.textContent = d.toLocaleDateString([], { weekday: 'short', day: 'numeric' }); const fx = i / (S.steps.length - 1); sp.style.left = `${fx * 100}%`; if (fx > 0.93) sp.style.transform = 'translateX(-100%)'; days.appendChild(sp); } });
}
const statRefresh = () => refresh();
function step(d) { const j = Math.max(0, Math.min(S.steps.length - 1, S.idx + d)); if (j === S.idx) { if (S.playing) togglePlay(); return; } S.idx = j; syncControls(); refresh(); }
function togglePlay() {
  S.playing = !S.playing; $('#play').textContent = S.playing ? '❚❚' : '▶'; $('#play').setAttribute('aria-label', S.playing ? 'Pause' : 'Play');
  clearInterval(S.timer);
  if (S.playing) { if (S.idx >= S.steps.length - 1) { S.idx = 0; } S.timer = setInterval(() => { if (!$('#busy').hidden) return; step(1); }, 650); }
}

/* ============================================================================= boot */
async function boot() {
  map = createMap('map');
  const [man, lab] = await Promise.all([
    fetch(`data/latest.json?t=${Date.now()}`, { cache: 'no-store' }).then((r) => { if (!r.ok) throw new Error('data/latest.json not found yet — the pipeline has not published a run'); return r.json(); }),
    fetch('basemap_labels.json').then((r) => r.json()),
  ]);
  S.man = man; S.base = `data/${man.base}`; S.steps = man.steps; man.steps.forEach((s) => S.byFhr.set(s.fhr, s));
  lab.places.sort((a, b) => b[3] - a[3]); labelData = lab;
  const hs = Object.fromEntries(new URLSearchParams(location.hash.slice(1)));
  S.product = ['qpf', 'snow', 'zr', 'pl', 'ptype'].includes(hs.p) && (!man.products || man.products.includes(hs.p)) ? hs.p : 'qpf';
  S.period = S.product === 'ptype' ? 'step' : hs.per === 'step' ? 'step' : [6, 12, 24, 48].includes(+hs.per) ? +hs.per : 24;
  const want = +hs.f || 24; let best = S.steps.length - 1; S.steps.forEach((s, i) => { if (s.fhr === want) best = i; }); S.idx = Math.min(best, S.steps.length - 1);
  if (hs.s) S.stat = hs.s === 'median' ? 'p50' : hs.s;
  if ([2, 4, 6, 8].includes(+hs.r)) S.lpmmR = +hs.r;
  await new Promise((res) => (map.loaded() ? res() : map.once('load', res)));
  wx = new WxLayer({ opacity: S.opacity }); wx.setDomain(man.domain);
  map.addLayer(wx, 'counties-hit');
  labels = new LabelLayer(map, $('#labels'), lab);
  const grat = () => { labels.step = updateGraticule(map); };
  map.on('moveend', grat); grat();
  map.addControl(new maplibregl.NavigationControl({ showCompass: false }), 'top-right');
  map.addControl(new maplibregl.ScaleControl({ unit: 'imperial', maxWidth: 110 }), 'bottom-left');
  map.addControl(new maplibregl.ScaleControl({ unit: 'metric', maxWidth: 110 }), 'bottom-left');
  map.on('mousemove', (e) => { if (!recentTouch()) showTip(e); });
  map.on('mouseout', () => { if (!recentTouch()) hideTip(); });
  map.on('click', (e) => { if (recentTouch() && S.pin) unpin(); pin(e); });
  map.on('dragstart', () => hideTip());
  buildUI(); syncControls();
  $('#splash').hidden = true;
  $('#t-gen').textContent = `Built ${new Date(man.generated).toLocaleString([], { month: 'short', day: 'numeric', hour: 'numeric', minute: '2-digit', timeZoneName: 'short' })}${man.complete ? '' : ' · run still filling in'}`;
  refresh();
}
boot().catch((e) => { console.error(e); $('#splash').hidden = false; $('#splash-msg').textContent = e.message; });
