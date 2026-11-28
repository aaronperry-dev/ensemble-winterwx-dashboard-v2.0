/* Daily maximum / minimum temperature helpers for the viewer (pure functions; runs in the browser and in Node tests).
 *
 *   MaxT of a date = highest 2 m temperature from 12Z to 00Z (7 am to 7 pm EST)
 *   MinT of a date = lowest 2 m temperature from 00Z to 12Z (7 pm to 7 am EST)
 * The pipeline publishes one file per product and date (manifest `temps.tmax` / `temps.tmin`); a window is shown only when
 * every pooled model has all of its 6-hour blocks.  Colours are the WPC maximum / minimum temperature legend (5 deg F classes).
 */
(function (root, factory) {
  const api = factory();
  if (typeof module === 'object' && module.exports) module.exports = api; else root.TempDays = api;
})(typeof self !== 'undefined' ? self : this, function () {
  'use strict';

  const PRODUCTS = ['tmax', 'tmin'];
  const WINDOWS_UTC = { tmax: [12, 24], tmin: [0, 12] };
  const NAMES = { tmax: 'Daily MaxT', tmin: 'Daily MinT' };
  const isTemp = (p) => p === 'tmax' || p === 'tmin';

  /* WPC legend, warmest first: 100+, 95-100, ... , -35 to -40, below -40 (30 classes, sampled from the supplied images). */
  const WPC_TOP_DOWN = ['#e71818', '#ff4d29', '#ff8e39', '#ffae29', '#ffcf73', '#ffeb9c', '#ffffb5', '#d6ffbd', '#a5d79c', '#73c373',
    '#31a252', '#086d29', '#10555a', '#009a94', '#31cfc6', '#adffff', '#73d7ff', '#4ac7ff', '#299eff', '#0065c6', '#103c9c', '#10007b',
    '#52288c', '#7369ad', '#9c9ac6', '#bdbede', '#d6dbe7', '#ffffff', '#ffbee7', '#ff71de'];
  /* Class edges in deg F, ascending; class 0 is "below -40" (its edge is far below any real temperature so it always colours). */
  const EDGES = [-200]; for (let t = -40; t <= 100; t += 5) EDGES.push(t);
  const COLORS = WPC_TOP_DOWN.slice().reverse();
  const CLASSES = { edges: EDGES, colors: COLORS };

  const cToF = (c) => c * 1.8 + 32;
  const fToC = (f) => (f - 32) / 1.8;

  /* "Chance of" thresholds (deg F): 'le' = at or below, 'ge' = at or above. */
  const PROB = {
    tmax: [['le', 20], ['le', 32], ['ge', 50], ['ge', 70], ['ge', 90], ['ge', 100]],
    tmin: [['le', -10], ['le', 0], ['le', 20], ['le', 32], ['ge', 60], ['ge', 70]],
  };
  const probLabel = ([dir, t]) => `${dir === 'le' ? '≤' : '≥'} ${t}°F`;
  const parseProb = (stat) => { const m = /^prob:(le|ge):(-?\d+(?:\.\d+)?)$/.exec(stat); return m ? [m[1], +m[2]] : null; };

  /** Window of a date, [start, end) in ms since the epoch (UTC). */
  function windowMs(product, date) {
    const [lo, hi] = WINDOWS_UTC[product], t0 = Date.parse(`${date}T00:00:00Z`);
    return [t0 + lo * 3600e3, t0 + hi * 3600e3];
  }

  /**
   * Dates that every pooled model can show for a product: the same date listed in each model's manifest, with all of its blocks
   * in (blocks === expected), a full window (expected >= blocks_per_window) and at least one member.
   * `temps` holds one manifest `temps` object per model (null/undefined when a model has none).
   * Returns [{date, entries: [one manifest entry per model, in the order given]}] by date.
   */
  function commonDays(product, temps) {
    if (!temps.length || temps.some((t) => !t || !Array.isArray(t[product]))) return [];
    const need = (t) => (t.blocks_per_window && t.blocks_per_window[product]) || 2;
    const ok = (t, e) => e.file && e.n > 0 && e.blocks >= e.expected && e.expected >= need(t);
    const maps = temps.map((t) => new Map(t[product].filter((e) => ok(t, e)).map((e) => [e.date, e])));
    return [...maps[0].keys()].filter((d) => maps.every((m) => m.has(d))).sort().map((date) => ({ date, entries: maps.map((m) => m.get(date)) }));
  }

  /** Index of the day nearest (in dates) to `date` (YYYY-MM-DD); 0 when there are none. */
  function nearestDay(days, date) {
    let best = 0, bd = Infinity; const t = Date.parse(`${date}T00:00:00Z`);
    days.forEach((d, i) => { const x = Math.abs(Date.parse(`${d.date}T00:00:00Z`) - t); if (x < bd) { bd = x; best = i; } });
    return best;
  }

  return { PRODUCTS, WINDOWS_UTC, NAMES, isTemp, WPC_TOP_DOWN, EDGES, COLORS, CLASSES, cToF, fToC, PROB, probLabel, parseProb, windowMs,
    commonDays, nearestDay };
});
