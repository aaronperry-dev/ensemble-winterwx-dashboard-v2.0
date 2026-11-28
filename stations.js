/* Station-value helpers for the viewer (pure functions; browser and Node tests).
 * The station layer prints the value of the map being shown (whatever product / statistic is selected) at each
 * New England METAR/ASOS site, so a station reads exactly what the map shows at that point. */
(function (root, factory) {
  const api = factory();
  if (typeof module === 'object' && module.exports) module.exports = api; else root.Stations = api;
})(typeof self !== 'undefined' ? self : this, function () {
  'use strict';

  const TYPE = ['RA', 'SN', 'ZR', 'PL'];                 // precip-type ids 0-3 (rain, snow, freezing rain, sleet)

  /** Drawing priority (the label layer skips a station whose label would overlap one drawn earlier): ASOS before AWOS, then longest record. */
  function order(list) {
    return list.slice().sort((a, b) => (a.type === 'ASOS' ? 0 : 1) - (b.type === 'ASOS' ? 0 : 1) || (b.seasons || 0) - (a.seasons || 0) || a.id.localeCompare(b.id));
  }

  /**
   * Label text for one station.
   *   product: 'qpf' | 'snow' | 'zr' | 'pl' | 'ptype' | 'tmax' | 'tmin'
   *   prob:    the map shows a percentage
   *   v:       value of the map's field at the station (inches, deg F or %), type: precip-type id for the 'ptype' map
   * Precision follows the tooltip: snowfall to tenths, liquid (QPF, freezing rain, sleet) to hundredths, temperatures to whole degrees
   * (a 0.5 deg grid does not support more), probabilities to whole percent.  Returns '' when there is nothing to print.
   */
  function format({ product, prob, v, type }) {
    if (v == null || !Number.isFinite(v)) return '';
    if (prob) return `${Math.round(v)}%`;
    if (product === 'ptype') return v > 0 && type >= 0 && type < TYPE.length ? TYPE[type] : '';
    if (product === 'tmax' || product === 'tmin') { const r = Math.round(v); return String(Object.is(r, -0) ? 0 : r); }
    return v.toFixed(product === 'snow' ? 1 : 2);
  }

  return { TYPE, order, format };
});
