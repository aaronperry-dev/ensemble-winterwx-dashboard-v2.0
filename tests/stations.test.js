// Run with:  node --test tests/
const test = require('node:test');
const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const S = require('../stations.js');

test('precision: snow tenths, liquid hundredths, temperatures whole degrees, probabilities whole percent', () => {
  assert.equal(S.format({ product: 'snow', v: 3.449 }), '3.4');
  assert.equal(S.format({ product: 'qpf', v: 0.256 }), '0.26');
  assert.equal(S.format({ product: 'zr', v: 0 }), '0.00');
  assert.equal(S.format({ product: 'pl', v: 0.1 }), '0.10');
  assert.equal(S.format({ product: 'tmax', v: 31.6 }), '32');
  assert.equal(S.format({ product: 'tmin', v: -0.4 }), '0');            // never "-0"
  assert.equal(S.format({ product: 'tmin', v: -12.5 }), '-12');
  assert.equal(S.format({ product: 'qpf', prob: true, v: 42.4 }), '42%');
  assert.equal(S.format({ product: 'tmax', prob: true, v: 99.6 }), '100%');
});

test('precip type prints the dominant type only where there is precipitation', () => {
  assert.equal(S.format({ product: 'ptype', v: 20, type: 1 }), 'SN');
  assert.equal(S.format({ product: 'ptype', v: 20, type: 2 }), 'ZR');
  assert.equal(S.format({ product: 'ptype', v: 0, type: 1 }), '');
  assert.equal(S.format({ product: 'ptype', v: 20, type: -1 }), '');
  assert.equal(S.format({ product: 'ptype', prob: true, v: 61, type: 1 }), '61%');
});

test('missing values print nothing', () => {
  for (const v of [null, undefined, NaN, Infinity]) assert.equal(S.format({ product: 'qpf', v }), '');
});

test('order puts ASOS first, then the longest records', () => {
  const l = [{ id: 'B', type: 'AWOS', seasons: 9 }, { id: 'C', type: 'ASOS', seasons: 2 }, { id: 'A', type: 'ASOS', seasons: 11 }, { id: 'D', type: 'ASOS', seasons: 11 }];
  assert.deepEqual(S.order(l).map((s) => s.id), ['A', 'D', 'C', 'B']);
  assert.equal(l[0].id, 'B');                                              // input is not modified
});

test('stations.json is a clean New England list', () => {
  const l = JSON.parse(fs.readFileSync(path.join(__dirname, '..', 'stations.json'), 'utf8'));
  assert.ok(l.length > 60);
  assert.equal(new Set(l.map((s) => s.id)).size, l.length);
  for (const s of l) {
    assert.ok(['CT', 'MA', 'ME', 'NH', 'RI', 'VT'].includes(s.state), s.id);
    assert.ok(s.lat > 40.9 && s.lat < 47.6 && s.lon > -73.8 && s.lon < -66.8, s.id);
    assert.ok(s.id && s.name && ['ASOS', 'AWOS'].includes(s.type), s.id);
  }
  for (const id of ['BOS', 'PWM', 'BGR', 'BTV', 'PVD', 'BDL', 'MHT', 'CAR']) assert.ok(l.some((s) => s.id === id), id);
});
