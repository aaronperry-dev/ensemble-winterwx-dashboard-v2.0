/* Custom cartographic basemap: Natural Earth + US Census boundaries drawn by MapLibre GL, plus a
 * collision-aware canvas label layer and graticule. No tile server, no API key, no external fonts. */
'use strict';

const BM = {
  ocean: '#cde1ee', land: '#f4f2ec', lake: '#cde1ee', river: '#a9cbe0', coast: '#6f93a8',
  county: '#c8c3b8', state: '#8d8779', border: '#3d3d3d', graticule: 'rgba(60,80,100,0.28)',
  font: '"Inter","Segoe UI",Roboto,"Helvetica Neue",Arial,sans-serif',
};
const MAX_BOUNDS = [[-165, 8], [-35, 74]];

function createMap(container) {
  const gj = (f) => ({ type: 'geojson', data: f, buffer: 0, tolerance: 0.2 });
  const style = {
    version: 8,
    sources: {
      land: gj('basemap_land.json'), lakes: gj('basemap_lakes.json'), rivers: gj('basemap_rivers.json'),
      coast: gj('basemap_coast.json'), borders: gj('basemap_borders.json'), admin1: gj('basemap_admin1.json'),
      counties: { type: 'geojson', data: 'basemap_counties.json', tolerance: 0.1 },
      graticule: { type: 'geojson', data: { type: 'FeatureCollection', features: [] } },
    },
    layers: [
      { id: 'ocean', type: 'background', paint: { 'background-color': BM.ocean } },
      { id: 'land', type: 'fill', source: 'land', paint: { 'fill-color': BM.land } },
      { id: 'lakes', type: 'fill', source: 'lakes', paint: { 'fill-color': BM.lake } },
      { id: 'rivers', type: 'line', source: 'rivers',
        paint: { 'line-color': BM.river, 'line-width': ['interpolate', ['linear'], ['zoom'], 3, 0.4, 8, 1.4] } },
      /* weather layer is inserted here (before 'counties-hit') */
      { id: 'counties-hit', type: 'fill', source: 'counties', paint: { 'fill-color': '#000', 'fill-opacity': 0.001 } },
      { id: 'counties-line', type: 'line', source: 'counties', minzoom: 4.8,
        paint: { 'line-color': BM.county,
                 'line-width': ['interpolate', ['linear'], ['zoom'], 5, 0.3, 8, 0.7, 10, 1.1],
                 'line-opacity': ['interpolate', ['linear'], ['zoom'], 4.8, 0, 5.6, 1] } },
      { id: 'admin1-hit', type: 'fill', source: 'admin1', paint: { 'fill-color': '#000', 'fill-opacity': 0.001 } },
      { id: 'admin1-line', type: 'line', source: 'admin1', layout: { 'line-join': 'round' },
        paint: { 'line-color': BM.state, 'line-width': ['interpolate', ['linear'], ['zoom'], 3, 0.5, 6, 0.9, 10, 1.6] } },
      { id: 'graticule', type: 'line', source: 'graticule',
        paint: { 'line-color': BM.graticule, 'line-width': 0.6, 'line-dasharray': [2, 3] } },
      { id: 'coast', type: 'line', source: 'coast',
        paint: { 'line-color': BM.coast, 'line-width': ['interpolate', ['linear'], ['zoom'], 3, 0.5, 8, 1] } },
      { id: 'borders-casing', type: 'line', source: 'borders', layout: { 'line-join': 'round' },
        paint: { 'line-color': 'rgba(255,255,255,0.75)', 'line-width': ['interpolate', ['linear'], ['zoom'], 3, 2.4, 8, 4.2] } },
      { id: 'borders', type: 'line', source: 'borders', layout: { 'line-join': 'round' },
        paint: { 'line-color': BM.border, 'line-width': ['interpolate', ['linear'], ['zoom'], 3, 1.1, 8, 2.2] } },
    ],
  };
  const map = new maplibregl.Map({
    container, style, center: [-72.5, 43.5], zoom: 5.2, minZoom: 2.6, maxZoom: 10.5, maxBounds: MAX_BOUNDS,
    dragRotate: false, pitchWithRotate: false, touchPitch: false, attributionControl: false,
    renderWorldCopies: false, fadeDuration: 0,
  });
  map.touchZoomRotate.disableRotation();
  map.keyboard.disableRotation();
  return map;
}

/* --------------------------------------------------------------------------- graticule */
function graticuleStep(zoom) { return zoom < 3.6 ? 10 : zoom < 5.2 ? 5 : zoom < 7 ? 2 : zoom < 8.6 ? 1 : 0.5; }

function updateGraticule(map) {
  const b = map.getBounds(), step = graticuleStep(map.getZoom()), feats = [];
  const w = Math.max(-179, Math.floor(b.getWest() / step) * step), e = Math.min(179, Math.ceil(b.getEast() / step) * step);
  const s = Math.max(-80, Math.floor(b.getSouth() / step) * step), n = Math.min(80, Math.ceil(b.getNorth() / step) * step);
  for (let lat = s; lat <= n; lat += step) feats.push({ type: 'Feature', properties: {}, geometry: { type: 'LineString', coordinates: [[w, lat], [e, lat]] } });
  for (let lon = w; lon <= e; lon += step) feats.push({ type: 'Feature', properties: {}, geometry: { type: 'LineString', coordinates: [[lon, s], [lon, n]] } });
  map.getSource('graticule').setData({ type: 'FeatureCollection', features: feats });
  return step;
}

/* --------------------------------------------------------------------------- labels */
class LabelLayer {
  constructor(map, canvas, data) {
    this.map = map; this.cv = canvas; this.ctx = canvas.getContext('2d'); this.d = data;
    this.show = { places: true, admin1: true, water: true, countries: true, graticule: true };
    this.step = 10; this.queued = false;
    const q = () => { if (!this.queued) { this.queued = true; requestAnimationFrame(() => { this.queued = false; this.draw(); }); } };
    map.on('render', q); window.addEventListener('resize', q);
    this.q = q;
  }

  resize() {
    const r = this.map.getContainer().getBoundingClientRect(), dpr = window.devicePixelRatio || 1;
    this.cv.width = Math.round(r.width * dpr); this.cv.height = Math.round(r.height * dpr);
    this.cv.style.width = r.width + 'px'; this.cv.style.height = r.height + 'px';
    this.w = r.width; this.h = r.height; this.dpr = dpr;
  }

  draw() {
    const { map, ctx } = this, z = map.getZoom();
    this.resize();
    ctx.setTransform(this.dpr, 0, 0, this.dpr, 0, 0);
    ctx.clearRect(0, 0, this.w, this.h);
    ctx.lineJoin = 'round'; ctx.textBaseline = 'middle';
    const placed = [], pad = 3;
    const hit = (b) => placed.some((p) => !(b.x2 < p.x1 || b.x1 > p.x2 || b.y2 < p.y1 || b.y1 > p.y2));
    const inView = (p) => p.x > -40 && p.y > -20 && p.x < this.w + 40 && p.y < this.h + 20;
    const text = (str, x, y, o) => {
      ctx.font = o.font; ctx.textAlign = o.align || 'center';
      if (o.spacing && 'letterSpacing' in ctx) ctx.letterSpacing = o.spacing; else if ('letterSpacing' in ctx) ctx.letterSpacing = '0px';
      const wd = ctx.measureText(str).width, ht = o.size + 2;
      const x1 = o.align === 'left' ? x : o.align === 'right' ? x - wd : x - wd / 2;
      const box = { x1: x1 - pad, x2: x1 + wd + pad, y1: y - ht / 2 - pad, y2: y + ht / 2 + pad };
      if (o.extra) { box.x1 = Math.min(box.x1, o.extra.x1); box.x2 = Math.max(box.x2, o.extra.x2); }
      if (hit(box)) return false;
      ctx.lineWidth = o.halo ?? 3; ctx.strokeStyle = o.haloColor || 'rgba(255,255,255,0.9)'; ctx.strokeText(str, x, y);
      ctx.fillStyle = o.color; ctx.fillText(str, x, y); placed.push(box); return true;
    };
    const proj = (lon, lat) => map.project([lon, lat]);
    const F = (size, w = 400, style = '') => `${style} ${w} ${size}px ${BM.font}`.trim();

    /* priority 1: country + water names */
    if (this.show.countries && z < 7) this.d.countries.forEach(([n, lon, lat, rank]) => {
      if (rank > 3 + (z > 4 ? 0 : -1)) return; const p = proj(lon, lat); if (!inView(p)) return;
      text(n.toUpperCase(), p.x, p.y, { font: F(z < 4 ? 11 : 13, 600), size: 13, color: '#6b6b6b', spacing: '3px', halo: 3 });
    });
    if (this.show.water && z < 6.5) {
      this.d.marine.forEach(([n, lon, lat, minz, cls]) => {
        const p = proj(lon, lat); if (!inView(p)) return;
        text(n, p.x, p.y, { font: F(cls === 'ocean' ? 13 : 11.5, 400, 'italic'), size: 13, color: '#4b7893', spacing: cls === 'ocean' ? '4px' : '1px', halo: 2.5 });
      });
      this.d.lakes.forEach(([n, lon, lat, minz]) => {
        if (z < Math.max(3.8, minz - 1.5)) return; const p = proj(lon, lat); if (!inView(p)) return;
        text(n, p.x, p.y, { font: F(11, 400, 'italic'), size: 11, color: '#4b7893', spacing: '1px', halo: 2.5 });
      });
    }
    /* priority 2: states / provinces */
    if (this.show.admin1 && z >= 3.9 && z < 9.2) this.d.admin1.forEach(([n, lon, lat, minz, c]) => {
      if (z < Math.max(3.9, minz - 1.6)) return; const p = proj(lon, lat); if (!inView(p)) return;
      text(n.toUpperCase(), p.x, p.y, { font: F(z < 5 ? 10 : 11.5, 600), size: 12, color: '#867f70', spacing: '2px', halo: 2.5 });
    });
    /* priority 3: populated places by population */
    if (this.show.places) {
      let drawn = 0;
      for (const [n, lon, lat, pop, minz, cap] of this.d.places) {
        if (drawn > 220) break;
        if (z < minz - 0.4 + (cap ? -1.2 : 0)) continue;
        const p = proj(lon, lat); if (!inView(p)) continue;
        const big = pop > 1e6, r = cap === 2 ? 4.2 : big ? 3.6 : pop > 1e5 ? 3 : 2.3;
        const size = big ? 13 : pop > 1e5 ? 12 : 11;
        const dot = { x1: p.x - r - 1, x2: p.x + r + 1, y1: p.y - r - 1, y2: p.y + r + 1 };
        if (hit(dot)) continue;
        ctx.font = F(size, big ? 600 : 500);
        const wd = ctx.measureText(n).width;
        const right = p.x + r + 4 + wd < this.w;
        const ok = text(n, right ? p.x + r + 4 : p.x - r - 4, p.y, { font: F(size, big ? 600 : 500), size, color: '#242424', align: right ? 'left' : 'right', halo: 3, extra: dot });
        if (!ok) continue;
        ctx.beginPath();
        if (cap === 2) this.star(p.x, p.y, r + 1); else ctx.arc(p.x, p.y, r, 0, Math.PI * 2);
        ctx.fillStyle = cap ? '#111' : '#2b2b2b'; ctx.fill(); ctx.lineWidth = 1.2; ctx.strokeStyle = '#fff'; ctx.stroke();
        placed.push(dot); drawn++;
      }
    }
    /* graticule labels along top and left edges */
    if (this.show.graticule) this.gratLabels(z);
  }

  star(x, y, r) {
    const c = this.ctx; c.moveTo(x, y - r);
    for (let i = 1; i < 10; i++) { const a = -Math.PI / 2 + (i * Math.PI) / 5, rr = i % 2 ? r * 0.45 : r; c.lineTo(x + Math.cos(a) * rr, y + Math.sin(a) * rr); }
    c.closePath();
  }

  gratLabels(z) {
    const { map, ctx } = this, b = map.getBounds(), step = this.step;
    ctx.font = `500 10px ${BM.font}`; ctx.fillStyle = '#4c5a66'; ctx.strokeStyle = 'rgba(255,255,255,0.85)'; ctx.lineWidth = 2.5;
    ctx.textAlign = 'left'; ctx.textBaseline = 'middle';
    for (let lat = Math.ceil(b.getSouth() / step) * step; lat <= b.getNorth(); lat += step) {
      const p = map.project([b.getWest(), lat]); if (p.y < 70 || p.y > this.h - 80) continue;
      const s = `${Math.abs(lat)}°${lat >= 0 ? 'N' : 'S'}`; ctx.strokeText(s, 4, p.y - 7); ctx.fillText(s, 4, p.y - 7);
    }
    ctx.textAlign = 'center';
    for (let lon = Math.ceil(b.getWest() / step) * step; lon <= b.getEast(); lon += step) {
      const p = map.project([lon, b.getNorth()]); if (p.x < 40 || p.x > this.w - 60) continue;
      const s = `${Math.abs(lon)}°${lon >= 0 ? 'E' : 'W'}`; ctx.strokeText(s, p.x, 10); ctx.fillText(s, p.x, 10);
    }
  }
}
