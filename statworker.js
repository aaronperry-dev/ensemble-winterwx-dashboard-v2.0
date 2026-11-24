/* Background worker for the slow ensemble statistics (PMM, LPMM) so the map never freezes while they compute. */
'use strict';
importScripts('ensstats.js');
self.onmessage = (e) => {
  const { id, data, M, ny, nx, stat, radius, weights } = e.data;
  try { const out = EnsStats.compute(data, M, ny, nx, stat, { radius, weights }); self.postMessage({ id, out }, [out.buffer]); }
  catch (err) { self.postMessage({ id, error: String(err && err.message || err) }); }
};
