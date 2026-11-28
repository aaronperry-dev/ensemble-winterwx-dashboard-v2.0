/* MapLibre custom WebGL2 layer: bilinear-interpolated, classed raster of a lat/lon grid.
 * The grid (0.5 deg) is uploaded as an R32F texture; the fragment shader inverts Web-Mercator per pixel,
 * interpolates the grid and classifies, so contours stay crisp at every zoom and updates only re-upload
 * one small texture.
 *   mode 0: value -> classed color (edges/colors uniforms)
 *   mode 1: precip type: value = dBZ-equivalent intensity, type texture selects the color ramp row
 */
'use strict';

const WX_VS = `#version 300 es
uniform mat4 u_matrix;
in vec2 a_merc;
out vec2 v_merc;
void main() { v_merc = a_merc; gl_Position = u_matrix * vec4(a_merc, 0.0, 1.0); }`;

const WX_FS = `#version 300 es
precision highp float; precision highp int;
uniform sampler2D u_val; uniform sampler2D u_typ; uniform sampler2D u_lut;
uniform vec4 u_dom;      // lon_west, lat_north, dlon, dlat
uniform vec2 u_dim;      // nx, ny
uniform float u_opacity; uniform int u_mode; uniform int u_n;
uniform float u_edges[32]; uniform vec4 u_cols[32];
uniform vec3 u_dbz;      // dbz min, step, count
in vec2 v_merc; out vec4 o;
const float PI = 3.14159265358979;
float F(ivec2 p) { p = clamp(p, ivec2(0), ivec2(u_dim) - 1); return texelFetch(u_val, p, 0).r; }
void main() {
  float lon = v_merc.x * 360.0 - 180.0;
  float lat = degrees(atan(sinh(PI * (1.0 - 2.0 * v_merc.y))));
  float gx = (lon - u_dom.x) / u_dom.z, gy = (u_dom.y - lat) / u_dom.w;
  if (gx < 0.0 || gy < 0.0 || gx > u_dim.x - 1.0 || gy > u_dim.y - 1.0) discard;
  ivec2 i0 = ivec2(floor(gx), floor(gy)); vec2 f = vec2(gx, gy) - vec2(i0);
  float v = mix(mix(F(i0), F(i0 + ivec2(1, 0)), f.x), mix(F(i0 + ivec2(0, 1)), F(i0 + ivec2(1, 1)), f.x), f.y);
  vec4 c = vec4(0.0);
  if (u_mode == 0) {
    for (int i = 0; i < 32; i++) { if (i < u_n && v >= u_edges[i]) c = u_cols[i]; }
  } else {
    if (v >= u_dbz.x) {
      ivec2 nn = clamp(ivec2(floor(vec2(gx, gy) + 0.5)), ivec2(0), ivec2(u_dim) - 1);
      int t = int(texelFetch(u_typ, nn, 0).r + 0.5);
      int k = clamp(int(floor((v - u_dbz.x) / u_dbz.y)), 0, int(u_dbz.z) - 1);
      c = texelFetch(u_lut, ivec2(k, t), 0);
    }
  }
  float a = c.a * u_opacity;
  if (a <= 0.0) discard;
  o = vec4(c.rgb * a, a);
}`;

class WxLayer {
  constructor(opts = {}) {
    this.id = 'wx'; this.type = 'custom'; this.renderingMode = '2d';
    this.opacity = opts.opacity ?? 0.82;
    this.mode = 0; this.n = 0;
    this.edges = new Float32Array(32); this.cols = new Float32Array(128);
    this.dom = null; this.ready = false; this.hasField = false;
    this.lutDims = [33, 4];
  }

  onAdd(map, gl) {
    this.map = map; this.gl = gl;
    if (!(gl instanceof WebGL2RenderingContext)) { console.error('WebGL2 is required'); return; }
    const sh = (type, src) => {
      const s = gl.createShader(type); gl.shaderSource(s, src); gl.compileShader(s);
      if (!gl.getShaderParameter(s, gl.COMPILE_STATUS)) throw new Error(gl.getShaderInfoLog(s));
      return s;
    };
    const p = gl.createProgram();
    gl.attachShader(p, sh(gl.VERTEX_SHADER, WX_VS)); gl.attachShader(p, sh(gl.FRAGMENT_SHADER, WX_FS));
    gl.linkProgram(p);
    if (!gl.getProgramParameter(p, gl.LINK_STATUS)) throw new Error(gl.getProgramInfoLog(p));
    this.prog = p;
    this.u = {};
    ['u_matrix', 'u_val', 'u_typ', 'u_lut', 'u_dom', 'u_dim', 'u_opacity', 'u_mode', 'u_n', 'u_edges', 'u_cols', 'u_dbz']
      .forEach((n) => { this.u[n] = gl.getUniformLocation(p, n); });
    this.aMerc = gl.getAttribLocation(p, 'a_merc');
    this.buf = gl.createBuffer();
    const mk = () => {
      const t = gl.createTexture(); gl.bindTexture(gl.TEXTURE_2D, t);
      gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.NEAREST);
      gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.NEAREST);
      gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE);
      gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);
      return t;
    };
    this.texVal = mk(); this.texTyp = mk(); this.texLut = mk();
    this.ready = true;
    if (this.dom) this._buildQuad();
  }

  static merc(lon, lat) {
    return [(lon + 180) / 360, 0.5 - Math.log(Math.tan(Math.PI / 4 + (lat * Math.PI) / 360)) / (2 * Math.PI)];
  }

  setDomain(d) {                    // {lat_north, lon_west, dlat, dlon, ny, nx}
    this.dom = d; if (this.ready) this._buildQuad();
  }

  _buildQuad() {
    const d = this.dom, gl = this.gl;
    const lonE = d.lon_west + (d.nx - 1) * d.dlon, latS = d.lat_north - (d.ny - 1) * d.dlat;
    const a = WxLayer.merc(d.lon_west, d.lat_north), b = WxLayer.merc(lonE, latS);
    const v = new Float32Array([a[0], a[1], b[0], a[1], a[0], b[1], a[0], b[1], b[0], a[1], b[0], b[1]]);
    gl.bindBuffer(gl.ARRAY_BUFFER, this.buf); gl.bufferData(gl.ARRAY_BUFFER, v, gl.STATIC_DRAW);
  }

  _upload(tex, data, w, h) {
    const gl = this.gl; gl.bindTexture(gl.TEXTURE_2D, tex);
    gl.pixelStorei(gl.UNPACK_ALIGNMENT, 1);
    gl.texImage2D(gl.TEXTURE_2D, 0, gl.R32F, w, h, 0, gl.RED, gl.FLOAT, data);
  }

  /** values/types: Float32Array(ny*nx) north->south rows. types only used in mode 1. */
  setField(values, types) {
    if (!this.ready) return;
    this._upload(this.texVal, values, this.dom.nx, this.dom.ny);
    this._upload(this.texTyp, types || new Float32Array(values.length), this.dom.nx, this.dom.ny);
    this.hasField = true; this.map && this.map.triggerRepaint();
  }

  /** Classed colors: edges ascending (display units), colors as CSS hex strings. */
  setClasses(edges, colors) {
    this.mode = 0; this.n = Math.min(32, edges.length);
    this.edges.fill(0); this.cols.fill(0);
    for (let i = 0; i < this.n; i++) {
      this.edges[i] = edges[i];
      const c = WxLayer.hex(colors[i]); this.cols.set([c[0], c[1], c[2], 1], i * 4);
    }
    this.map && this.map.triggerRepaint();
  }

  /** ramps: {RA:[hex x33], SN:..., FZRA:..., PL:...}; order of type ids: RA, SN, FZRA, PL. */
  setPtypeLut(ramps, order, dbzMin, dbzStep) {
    this.mode = 1; this.dbz = [dbzMin, dbzStep, ramps[order[0]].length];
    const w = ramps[order[0]].length, data = new Uint8Array(w * order.length * 4);
    order.forEach((name, t) => ramps[name].forEach((hex, k) => {
      const c = WxLayer.hex(hex), o = (t * w + k) * 4;
      const white = c[0] > 0.99 && c[1] > 0.99 && c[2] > 0.99;       // legend's blank first box -> transparent
      data.set([Math.round(c[0] * 255), Math.round(c[1] * 255), Math.round(c[2] * 255), white ? 0 : 255], o);
    }));
    if (this.ready) {
      const gl = this.gl; gl.bindTexture(gl.TEXTURE_2D, this.texLut);
      gl.pixelStorei(gl.UNPACK_ALIGNMENT, 1);
      gl.texImage2D(gl.TEXTURE_2D, 0, gl.RGBA8, w, order.length, 0, gl.RGBA, gl.UNSIGNED_BYTE, data);
    }
    this.map && this.map.triggerRepaint();
  }

  setOpacity(o) { this.opacity = o; this.map && this.map.triggerRepaint(); }

  static hex(h) {
    const n = parseInt(h.slice(1), 16); return [(n >> 16 & 255) / 255, (n >> 8 & 255) / 255, (n & 255) / 255];
  }

  render(gl, arg) {
    if (!this.ready || !this.hasField || !this.dom) return;
    const m = arg && arg.defaultProjectionData ? arg.defaultProjectionData.mainMatrix : arg;   // v4 + v5
    gl.useProgram(this.prog);
    gl.uniformMatrix4fv(this.u.u_matrix, false, m);
    const d = this.dom;
    gl.uniform4f(this.u.u_dom, d.lon_west, d.lat_north, d.dlon, d.dlat);
    gl.uniform2f(this.u.u_dim, d.nx, d.ny);
    gl.uniform1f(this.u.u_opacity, this.opacity);
    gl.uniform1i(this.u.u_mode, this.mode);
    gl.uniform1i(this.u.u_n, this.n);
    gl.uniform1fv(this.u.u_edges, this.edges);
    gl.uniform4fv(this.u.u_cols, this.cols);
    gl.uniform3f(this.u.u_dbz, ...(this.dbz || [8, 2, 33]));
    [[this.texVal, 0, 'u_val'], [this.texTyp, 1, 'u_typ'], [this.texLut, 2, 'u_lut']].forEach(([t, i, n]) => {
      gl.activeTexture(gl.TEXTURE0 + i); gl.bindTexture(gl.TEXTURE_2D, t); gl.uniform1i(this.u[n], i);
    });
    gl.activeTexture(gl.TEXTURE0);
    gl.enable(gl.BLEND); gl.blendFunc(gl.ONE, gl.ONE_MINUS_SRC_ALPHA);
    gl.disable(gl.DEPTH_TEST);
    gl.bindBuffer(gl.ARRAY_BUFFER, this.buf);
    gl.enableVertexAttribArray(this.aMerc);
    gl.vertexAttribPointer(this.aMerc, 2, gl.FLOAT, false, 0, 0);
    gl.drawArrays(gl.TRIANGLES, 0, 6);
  }

  onRemove() { /* GL objects are released with the context */ }
}
