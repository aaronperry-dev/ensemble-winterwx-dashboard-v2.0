/* Color ramps sampled from the supplied NWS / WeatherBELL legends. */
'use strict';

/* Classed legends. `min` is the lower edge of each class in INCHES; class 0 (below the first edge) is transparent. */
const RAMPS = {
  qpf: {
    title: 'Liquid-equivalent precipitation (QPF)',
    unit: 'in',
    edges: [0.01, 0.10, 0.25, 0.50, 1.00, 1.50, 2.00, 3.00, 4.00, 6.00, 8.00, 10.00, 15.00, 20.00, 30.00],
    labels: ['0.01\u20130.09', '0.10\u20130.24', '0.25\u20130.49', '0.50\u20130.99', '1.00\u20131.49', '1.50\u20131.99',
             '2.00\u20132.99', '3.00\u20133.99', '4.00\u20135.99', '6.00\u20137.99', '8.00\u20139.99', '10.00\u201314.99',
             '15.00\u201319.99', '20.00\u201329.99', '30.00+'],
    colors: ['#c7e9c0', '#a1d99b', '#74c476', '#31a353', '#006d2c', '#ffff8a', '#ffcc4f', '#fe8d3c', '#fc4e2a',
             '#d61a1c', '#ad0026', '#700026', '#3b0030', '#4c0073', '#ffdbff'],
  },
  snow: {
    title: 'Snowfall',
    unit: 'in',
    edges: [0.05, 1, 2, 3, 4, 6, 8, 12, 18, 24, 30, 36],
    labels: ['< 1', '1\u20132', '2\u20133', '3\u20134', '4\u20136', '6\u20138', '8\u201312', '12\u201318', '18\u201324',
             '24\u201330', '30\u201336', '36+'],
    colors: ['#bdd7e7', '#6baed6', '#3182bd', '#08519c', '#082694', '#ffff96', '#ffc400', '#ff8700', '#db1400',
             '#9e0000', '#690000', '#360000'],
  },
};

/* Precipitation-type ramps: 33 steps, 8 -> 72 dBZ in 2 dBZ increments (HRRR-style legend). */
const PTYPE = {
  names: ['RA', 'SN', 'FZRA', 'PL'],
  long: ['Rain', 'Snow', 'Freezing rain', 'Sleet'],
  dbzMin: 8, dbzStep: 2,
  ramps: {
    RA: ["#ffffff", "#c5c5c6", "#a3a3a4", "#767676", "#5a5a5a", "#448944", "#28bf28", "#0ee70e", "#0af50a", "#24db24", "#30c930", "#28a728", "#2e911e", "#74b314", "#dde502", "#efe736", "#efe162", "#eec168", "#ee8d37", "#ee6e22", "#f75626", "#fb442b", "#f13435", "#b42f5a", "#b13060", "#b13060", "#b13060", "#b43a84", "#c14dd7", "#dd24e8", "#ee16f5", "#fe3ffe", "#ffd5ff"],
    SN: ["#c9ffff", "#bbf5ff", "#adebff", "#a1e1ff", "#93d7ff", "#93d7ff", "#85cdff", "#78c1ff", "#6ab7ff", "#5cadff", "#50a3ff", "#4299ff", "#348fff", "#348fff", "#2885ff", "#1a7aff", "#0064ff", "#005ef5", "#0058eb", "#0050e1", "#024bd5", "#024bd5", "#0042cb", "#003cc0", "#0036b7", "#002eab", "#0028a1", "#002297", "#001a8d", "#001a8d", "#001481", "#000d76", "#00076c"],
    FZRA: ["#fc92ef", "#fc92ef", "#f386e9", "#f386e9", "#f386e9", "#eb78df", "#eb78df", "#eb78df", "#eb78df", "#e36cd7", "#e36cd7", "#dc5fcf", "#dc5fcf", "#d14dc3", "#d14dc3", "#d14dc3", "#d14dc3", "#c740ba", "#c740ba", "#c740ba", "#c740ba", "#bd2faf", "#bd2faf", "#bd2faf", "#b524a7", "#b524a7", "#b31ea3", "#af189f", "#a9119b", "#a9119b", "#a9119b", "#8f0c91", "#5e0881"],
    PL: ["#f9ddcc", "#f9ddcc", "#f9ddcc", "#f9ddcc", "#f9d0c1", "#f9d0c1", "#f9c5b7", "#f9c5b7", "#f9c5b7", "#f9bbae", "#f9bbae", "#f8b5a4", "#f9aca0", "#f9aca0", "#f9aca0", "#f99f94", "#f99f94", "#f9958b", "#f9958b", "#f9958b", "#f88d85", "#f88d85", "#f98078", "#f98078", "#f98078", "#f98078", "#f8766f", "#f7706e", "#f96963", "#f96963", "#f96963", "#f95957", "#f95957"],
  },
};

/* Freezing rain and sleet accumulations (liquid equivalent, inches). No legend was supplied for these, so the
 * classes follow the magenta (freezing rain) and salmon (sleet) families of the precipitation-type legend. */
RAMPS.zr = {
  title: 'Freezing rain (liquid equiv.)',
  unit: 'in',
  edges: [0.01, 0.05, 0.10, 0.15, 0.25, 0.35, 0.50, 0.75, 1.00],
  labels: ['0.01–0.04', '0.05–0.09', '0.10–0.14', '0.15–0.24', '0.25–0.34', '0.35–0.49', '0.50–0.74', '0.75–0.99', '1.00+'],
  colors: ['#f9c6f3', '#fc92ef', '#eb78df', '#dc5fcf', '#c740ba', '#b524a7', '#a9119b', '#8f0c91', '#5e0881'],
};
RAMPS.pl = {
  title: 'Sleet / ice pellets (liquid equiv.)',
  unit: 'in',
  edges: [0.01, 0.05, 0.10, 0.25, 0.50, 0.75, 1.00, 1.50],
  labels: ['0.01–0.04', '0.05–0.09', '0.10–0.24', '0.25–0.49', '0.50–0.74', '0.75–0.99', '1.00–1.49', '1.50+'],
  colors: ['#f9ddcc', '#f9c5b7', '#f9aca0', '#f98078', '#f96963', '#e8433f', '#c92a2a', '#8f1a1a'],
};
