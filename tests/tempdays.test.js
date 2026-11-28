// Run with:  node --test tests/
const test = require('node:test');
const assert = require('node:assert/strict');
const T = require('../tempdays.js');

const ent = (date, o = {}) => ({ date, index: 0, file: `tmp/tmax_${date}.bin.gz`, expected: 2, blocks: 2, n: 31, mask: 1, ...o });
const man = (list, bpw = { tmax: 2, tmin: 2 }) => ({ blocks_per_window: bpw, tmax: list, tmin: [] });

test('the WPC ramp has one colour per class, warmest on top, edges every 5 deg F', () => {
  assert.equal(T.COLORS.length, T.EDGES.length);
  assert.equal(T.WPC_TOP_DOWN.length, 30);
  assert.deepEqual(T.EDGES.slice(0, 3), [-200, -40, -35]);
  assert.equal(T.EDGES.at(-1), 100);
  T.EDGES.slice(1).forEach((e, i, a) => { if (i) assert.equal(e - a[i - 1], 5); });
  assert.equal(T.COLORS.at(-1), '#e71818');          // 100 F and up
  assert.equal(T.COLORS[0], '#ff71de');              // below -40 F
  assert.equal(T.COLORS[T.EDGES.indexOf(32 - 2)], T.COLORS[T.EDGES.indexOf(30)]);
  assert.equal(T.COLORS[T.EDGES.indexOf(-35)], '#ffffff');
});

test('windows: MaxT is 12Z-00Z and MinT is 00Z-12Z of the labelled UTC date', () => {
  const [a, b] = T.windowMs('tmax', '2026-11-28'), [c, d] = T.windowMs('tmin', '2026-11-28');
  assert.equal(new Date(a).toISOString(), '2026-11-28T12:00:00.000Z');
  assert.equal(new Date(b).toISOString(), '2026-11-29T00:00:00.000Z');
  assert.equal(new Date(c).toISOString(), '2026-11-28T00:00:00.000Z');
  assert.equal(new Date(d).toISOString(), '2026-11-28T12:00:00.000Z');
});

test('unit conversion and probability parsing', () => {
  assert.equal(T.cToF(0), 32); assert.ok(Math.abs(T.fToC(212) - 100) < 1e-9);
  assert.deepEqual(T.parseProb('prob:le:-10'), ['le', -10]);
  assert.deepEqual(T.parseProb('prob:ge:32.5'), ['ge', 32.5]);
  assert.equal(T.parseProb('prob:0.5'), null); assert.equal(T.parseProb('mean'), null);
  assert.equal(T.probLabel(['le', 32]), '≤ 32°F');
  T.PRODUCTS.forEach((p) => T.PROB[p].forEach(([d, t]) => assert.equal(T.probLabel([d, t]).length > 0, true)));
});

test('commonDays keeps only complete windows that every pooled model has', () => {
  const g = man([ent('2026-11-28'), ent('2026-11-29'), ent('2026-11-30', { blocks: 1 }), ent('2026-12-01', { expected: 1, blocks: 1 })]);
  const i = man([ent('2026-11-29'), ent('2026-11-30'), ent('2026-12-01', { n: 0, file: null })]);
  assert.deepEqual(T.commonDays('tmax', [g]).map((d) => d.date), ['2026-11-28', '2026-11-29']);   // half window (expected 1) and unfinished blocks are out
  const both = T.commonDays('tmax', [g, i]);
  assert.deepEqual(both.map((d) => d.date), ['2026-11-29']);
  assert.equal(both[0].entries.length, 2);                                                          // one entry per model, in order
  assert.deepEqual(T.commonDays('tmin', [g, i]), []);
});

test('commonDays tolerates models without temperatures', () => {
  assert.deepEqual(T.commonDays('tmax', [man([ent('2026-11-28')]), null]), []);
  assert.deepEqual(T.commonDays('tmax', [man([ent('2026-11-28')]), {}]), []);
  assert.deepEqual(T.commonDays('tmax', []), []);
  assert.equal(T.commonDays('tmax', [man([ent('2026-11-28')], { tmax: 3 })]).length, 0);            // a window needing 3 blocks is not done at 2
});

test('nearestDay picks the closest date', () => {
  const d = [{ date: '2026-11-28' }, { date: '2026-11-29' }, { date: '2026-12-02' }];
  assert.equal(T.nearestDay(d, '2026-11-29'), 1); assert.equal(T.nearestDay(d, '2026-12-01'), 2);
  assert.equal(T.nearestDay(d, '2026-01-01'), 0); assert.equal(T.nearestDay([], '2026-01-01'), 0);
});
