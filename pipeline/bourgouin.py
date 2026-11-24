"""Revised (Modified) Bourgouin precipitation-type algorithm.

Implements Birk, Lenning, Donofrio & Friedlein (2021), "A Revised Bourgouin Precipitation-Type
Algorithm", Wea. Forecasting 36, 425-438 (equation numbers below refer to that paper), fully
vectorised over grid columns.

Inputs are model vertical profiles; outputs are the probabilities (0-100) of the four basic
types RA, SN, FZRA, PL for every column, plus diagnostics.

Implementation notes / documented choices
-----------------------------------------
* Energies use the WET-BULB profile (Eq. 1 with Tw in place of T), integrated with trapezoids and
  exact linear zero-crossings, in J/kg.  Heights come from hypsometric integration of the virtual
  temperature profile from the surface (so no model HGT fields have to be downloaded).
* METw (melting energy) / RETw (refreezing energy) selection follows the paper's appendix:
    - ProbSN uses ALL melting energy in the column (below the precipitation-generation-layer top).
    - ProbPL / ProbFZRA use RETw of the lowest cold layer that has a warm layer somewhere above it
      ("near-surface cold layer") and METw = melting energy above that cold layer. If there is no
      such cold layer, RETw = 0 and METw = all melting energy (surface-based melting layer).
* ProbFZRA is relabelled ProbRA when surface Tw > 0 degC (paper, Appendix c).
* Ice nucleation (ProbIce, Eq. 2): precipitation-generation layer = RH_ice > 75 % layer deeper
  than 1 km; a further such layer separated by a > 1.5 km dry gap is dropped (sublimation).
  Where no layer is deeper than 1 km, the coldest saturated level of any depth is used; where
  nothing is saturated at all, ice is assumed present (ProbIce = 100).
* Model RH is relative to water; RH_ice is derived from it for T < 0 degC.
"""
from __future__ import annotations

import numpy as np

G = 9.80665          # m s-2
RD = 287.04          # J kg-1 K-1
T0 = 273.15          # K
S_MAX = 8            # max alternating warm/cold segments tracked per column

RA, SN, FZRA, PL = 0, 1, 2, 3
TYPE_NAMES = ("RA", "SN", "FZRA", "PL")


# --------------------------------------------------------------------------- thermodynamics
def es_water(tc):
    """Saturation vapour pressure over water (hPa), Bolton (1980)."""
    return 6.112 * np.exp(17.67 * tc / (tc + 243.5))


def es_ice(tc):
    """Saturation vapour pressure over ice (hPa)."""
    return 6.112 * np.exp(22.46 * tc / (272.62 + tc))


def wet_bulb_c(tc, rh, p_hpa, iters: int = 8):
    """Isobaric wet-bulb temperature (degC) by Newton iteration on the psychrometer equation."""
    tc = np.asarray(tc, dtype=np.float64)
    rh = np.clip(np.asarray(rh, dtype=np.float64), 1.0, 100.0)
    p = np.asarray(p_hpa, dtype=np.float64)
    e = rh / 100.0 * es_water(tc)
    tw = tc - (100.0 - rh) / 5.0
    for _ in range(iters):
        tw = np.clip(tw, tc - 40.0, tc)
        esw = es_water(tw)
        a = 6.6e-4 * (1.0 + 1.15e-3 * tw)
        f = esw - a * p * (tc - tw) - e
        fp = esw * 17.67 * 243.5 / (tw + 243.5) ** 2 + a * p
        tw = tw - f / fp
    return np.minimum(tw, tc)


def prob_ice_from_temp_c(t):
    """Eq. 2a-c: probability of ice in the precipitation-generation layer from its minimum T (degC)."""
    t = np.asarray(t, dtype=np.float64)
    poly = -0.065 * t**4 - 3.1544 * t**3 - 56.414 * t**2 - 449.6 * t - 1308.0
    return np.where(t <= -15.0, 100.0, np.where(t >= -7.0, 0.0, np.clip(poly, 0.0, 100.0)))


# --------------------------------------------------------------------------- profile prep
def build_profile(levels_hpa, t_lev_k, rh_lev, t2_k, rh2, ps_hpa):
    """Stack surface + pressure levels into (K, C) arrays ordered by decreasing pressure.

    Levels at/below the ground are clamped onto the surface (zero thickness) so they contribute
    nothing.  Returns p (hPa), tc (degC), rh (%).
    """
    levels = np.asarray(levels_hpa, dtype=np.float64)
    L, C = t_lev_k.shape
    ps = np.asarray(ps_hpa, dtype=np.float64)
    p = np.empty((L + 1, C))
    tc = np.empty((L + 1, C))
    rh = np.empty((L + 1, C))
    p[0], tc[0], rh[0] = ps, t2_k - T0, np.clip(rh2, 1.0, 100.0)
    for k in range(L):
        below = levels[k] >= ps - 0.25
        p[k + 1] = np.where(below, ps, levels[k])
        tc[k + 1] = np.where(below, tc[0], t_lev_k[k] - T0)
        rh[k + 1] = np.where(below, rh[0], np.clip(rh_lev[k], 1.0, 100.0))
    return p, tc, rh


def layer_heights(p, tc, rh):
    """Heights above ground (m) of each profile level via the hypsometric equation with Tv."""
    e = rh / 100.0 * es_water(tc)
    q = 0.622 * e / (p - 0.378 * e)
    tv = (tc + T0) * (1.0 + 0.608 * q)
    dz = (RD / G) * 0.5 * (tv[:-1] + tv[1:]) * np.log(p[:-1] / p[1:])
    z = np.zeros_like(p)
    z[1:] = np.cumsum(dz, axis=0)
    return z


# --------------------------------------------------------------------------- ice nucleation
def ice_nucleation(tc, rh, z):
    """Return (ProbIce 0-100, index of the top level of the precipitation-generation layer)."""
    K, C = tc.shape
    rhi = np.where(tc < 0.0, rh * es_water(tc) / es_ice(tc), rh)
    sat = rhi > 75.0

    in_run = np.zeros(C, bool)
    run_z0 = np.zeros(C)
    run_tmin = np.full(C, np.inf)
    run_end_z = np.zeros(C)
    run_end_idx = np.zeros(C, int)
    have = np.zeros(C, bool)
    blocked = np.zeros(C, bool)
    last_top_z = np.zeros(C)
    gen_tmin = np.full(C, np.inf)
    gen_top = np.zeros(C, int)
    any_tmin = np.full(C, np.inf)
    any_top = np.zeros(C, int)

    for k in range(K + 1):                       # one extra pass closes runs that reach the top
        sat_k = sat[k] if k < K else np.zeros(C, bool)
        tk = tc[k] if k < K else np.zeros(C)
        zk = z[k] if k < K else np.zeros(C)
        start = sat_k & ~in_run
        cont = sat_k & in_run
        run_z0 = np.where(start, zk, run_z0)
        run_tmin = np.where(start, tk, np.where(cont, np.minimum(run_tmin, tk), run_tmin))
        run_end_z = np.where(sat_k, zk, run_end_z)
        run_end_idx = np.where(sat_k, k, run_end_idx)

        closing = ~sat_k & in_run
        qual = closing & ((run_end_z - run_z0) > 1000.0)
        gap = run_z0 - last_top_z
        accept = qual & ~blocked & (~have | (gap <= 1500.0))
        blocked |= qual & ~blocked & have & (gap > 1500.0)      # dry gap > 1.5 km: upper layer sublimates
        gen_tmin = np.where(accept, np.minimum(gen_tmin, run_tmin), gen_tmin)
        gen_top = np.where(accept, run_end_idx, gen_top)
        last_top_z = np.where(accept, run_end_z, last_top_z)
        have |= accept
        any_tmin = np.where(closing, np.minimum(any_tmin, run_tmin), any_tmin)
        any_top = np.where(closing, np.maximum(any_top, run_end_idx), any_top)
        in_run = sat_k

    tmin = np.where(have, gen_tmin, any_tmin)
    top = np.where(have, gen_top, np.where(np.isfinite(any_tmin), any_top, K - 1))
    prob_ice = np.where(np.isfinite(tmin), prob_ice_from_temp_c(np.where(np.isfinite(tmin), tmin, 0.0)), 100.0)
    return prob_ice, np.maximum(top, 1)


# --------------------------------------------------------------------------- energies
def layer_energies(tw, z, top):
    """Alternating warm/cold segment energies (J/kg, positive magnitudes), bottom to top.

    Returns seg_E (S, C), seg_warm (S, C) bool, nseg (C,) = index of the last segment (-1 = none).
    """
    K, C = tw.shape
    seg_E = np.zeros((S_MAX, C))
    seg_warm = np.zeros((S_MAX, C), bool)
    nseg = np.full(C, -1, int)
    cur = np.full(C, -1, np.int8)
    gk = G / T0

    def feed(warm, area, valid):
        valid = valid & (area > 1e-9)
        sign = warm.astype(np.int8)
        new = np.nonzero(valid & (sign != cur))[0]
        nseg[new] = np.minimum(nseg[new] + 1, S_MAX - 1)
        seg_warm[nseg[new], new] = warm[new]
        cur[new] = sign[new]
        v = np.nonzero(valid)[0]
        seg_E[nseg[v], v] += area[v]

    for k in range(K - 1):
        inc = (k + 1) <= top
        w0, w1 = tw[k], tw[k + 1]
        dz = z[k + 1] - z[k]
        same = (w0 > 0) == (w1 > 0)
        # non-crossing layer: one trapezoid
        feed(w0 > 0, np.abs(0.5 * (w0 + w1)) * dz * gk, inc & same)
        # crossing layer: two triangles split at the zero crossing
        denom = np.where(w0 == w1, 1.0, w0 - w1)
        f = np.clip(w0 / denom, 0.0, 1.0)
        cross = inc & ~same
        feed(w0 > 0, 0.5 * np.abs(w0) * f * dz * gk, cross)
        feed(w1 > 0, 0.5 * np.abs(w1) * (1.0 - f) * dz * gk, cross)
    return seg_E, seg_warm, nseg


def summarize_energies(seg_E, seg_warm, nseg):
    """Collapse segments into (total ME, near-surface-cold RE, ME above that layer)."""
    S, C = seg_E.shape
    sidx = np.arange(S)[:, None]
    valid = sidx <= nseg[None, :]
    warm = valid & seg_warm
    cold = valid & ~seg_warm
    total_me = np.where(warm, seg_E, 0.0).sum(0)

    warm_above = np.zeros((S, C), bool)           # is there any warm segment above segment s?
    run = np.zeros(C, bool)
    for s in range(S - 1, -1, -1):
        warm_above[s] = run
        run = run | warm[s]
    key = cold & warm_above
    has_key = key.any(0)
    key_idx = key.argmax(0)
    cols = np.arange(C)
    re = np.where(has_key, seg_E[key_idx, cols], 0.0)
    me_above = np.where(has_key, np.where(warm & (sidx > key_idx[None, :]), seg_E, 0.0).sum(0), total_me)
    return total_me, re, me_above


# --------------------------------------------------------------------------- probabilities
def type_probabilities(total_me, re, me_above, prob_ice, surface_warm):
    """Appendix equations -> float arrays RA, SN, FZRA, PL in [0, 100]."""
    pi = prob_ice / 100.0
    sn = pi * np.clip(1540.0 * np.exp(-0.29 * total_me), 0.0, 100.0)                    # Eq. 9

    pl_i = np.where((me_above > 0) & (re > 0),
                    2.3 * re - 42.0 * np.log(me_above + 1.0) + 3.0, 0.0)               # Eq. 8
    pl = pi * np.clip(pl_i, 0.0, 100.0)

    fz_i = np.clip(-2.1 * re + 0.2 * me_above + 458.0, 0.0, 100.0)                      # Eq. 7
    fz_i = np.where(me_above < 5.0, fz_i * 0.2 * me_above, fz_i)
    fz = np.clip((100.0 - prob_ice) + pi * fz_i, 0.0, 100.0)

    ra = np.where(surface_warm, fz, 0.0)
    fzra = np.where(surface_warm, 0.0, fz)
    return ra, np.clip(sn, 0, 100), fzra, np.clip(pl, 0, 100)


def birk_from_profile(p, tc, rh):
    """Run the full algorithm on a prepared (K, C) profile. Returns dict of (C,) arrays."""
    z = layer_heights(p, tc, rh)
    tw = wet_bulb_c(tc, rh, p)
    prob_ice, top = ice_nucleation(tc, rh, z)
    seg_E, seg_warm, nseg = layer_energies(tw, z, top)
    total_me, re, me_above = summarize_energies(seg_E, seg_warm, nseg)
    ra, sn, fzra, pl = type_probabilities(total_me, re, me_above, prob_ice, tw[0] > 0.0)
    return {
        "probs": np.stack([ra, sn, fzra, pl]),
        "prob_ice": prob_ice, "me_total": total_me, "me_above": me_above, "re": re, "tw_sfc": tw[0],
    }


def kuchera_slr(p, tc, mode: str = "kuchera", lo: float = 5.0, hi: float = 30.0):
    """Snow-to-liquid ratio per column.

    "kuchera": warmest temperature between the surface and 500 hPa (K):
        SLR = 12 + 2 (271.16 - Tmax)  if Tmax > 271.16 K   else   12 + (271.16 - Tmax)
    "fixed:<r>": constant ratio r.   Result clamped to [lo, hi].
    """
    if mode.startswith("fixed:"):
        return np.full(tc.shape[1], float(mode.split(":", 1)[1]), dtype=np.float64)
    tmax = np.where(p >= 499.9, tc + T0, -np.inf).max(axis=0)
    slr = np.where(tmax > 271.16, 12.0 + 2.0 * (271.16 - tmax), 12.0 + (271.16 - tmax))
    return np.clip(slr, lo, hi)


def classify_profiles(levels_hpa, t_lev_k, rh_lev, t2_k, rh2, ps_pa, slr_mode="kuchera",
                      slr_min=5.0, slr_max=30.0):
    """High-level entry point for 2-D fields.

    Parameters are arrays shaped (L, ny, nx) for levels and (ny, nx) for surface fields
    (T in K, RH in %, surface pressure in Pa).  Returns
        probs  uint8 (4, ny, nx)   [RA, SN, FZRA, PL] in percent
        slr    float32 (ny, nx)
        diag   dict of float32 (ny, nx) arrays
    """
    L, ny, nx = t_lev_k.shape
    C = ny * nx
    p, tc, rh = build_profile(levels_hpa, t_lev_k.reshape(L, C), rh_lev.reshape(L, C),
                              t2_k.reshape(C), rh2.reshape(C), ps_pa.reshape(C) / 100.0)
    out = birk_from_profile(p, tc, rh)
    slr = kuchera_slr(p, tc, slr_mode, slr_min, slr_max)
    probs = np.rint(out["probs"]).astype(np.uint8).reshape(4, ny, nx)
    diag = {k: v.astype(np.float32).reshape(ny, nx) for k, v in out.items() if k != "probs"}
    return probs, slr.astype(np.float32).reshape(ny, nx), diag


def dominant_type(probs):
    """Index of most likely type; ties broken by impact: FZRA > PL > SN > RA."""
    order = np.array([FZRA, PL, SN, RA])
    pr = np.asarray(probs)[order].astype(np.float64)
    return order[np.argmax(pr, axis=0)]
