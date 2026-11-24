import numpy as np
import pytest

from pipeline import bourgouin as B
from pipeline import config as C
from pipeline.products import snowfall_step_cm, type_shares

LEV = np.array(C.PLEVELS, float)
RA, SN, FZRA, PL = 0, 1, 2, 3


def classify(tfun, rh=95.0, ps=1013.0, slr_mode="kuchera"):
    t = (np.array([tfun(p) for p in LEV]) + 273.15).reshape(-1, 1, 1)
    probs, slr, diag = B.classify_profiles(
        LEV, t, np.full(t.shape, rh), np.array([[tfun(ps) + 273.15]]), np.array([[rh]]),
        np.array([[ps * 100.0]]), slr_mode=slr_mode)
    return probs[:, 0, 0].astype(int), float(slr[0, 0]), diag


def test_wet_bulb_values_and_bounds():
    assert B.wet_bulb_c(np.array([5.0]), np.array([100.0]), np.array([1000.0]))[0] == pytest.approx(5.0, abs=1e-6)
    assert B.wet_bulb_c(np.array([20.0]), np.array([50.0]), np.array([1000.0]))[0] == pytest.approx(13.8, abs=0.3)
    t = np.linspace(-20, 30, 11)
    tw = B.wet_bulb_c(t, np.full_like(t, 40.0), np.full_like(t, 900.0))
    assert np.all(tw <= t + 1e-9)


def test_prob_ice_eq2():
    assert B.prob_ice_from_temp_c(-6.9) == 0.0
    assert B.prob_ice_from_temp_c(-15.0) == 100.0
    assert B.prob_ice_from_temp_c(-30.0) == 100.0
    mid = [float(B.prob_ice_from_temp_c(t)) for t in (-8, -10, -12, -14)]
    assert mid == sorted(mid) and 0.0 < mid[0] < mid[-1] <= 100.0


def test_energy_integral_isothermal_cold_layer():
    # 1 km isothermal layer at Tw = -5 C  ->  g/T0 * 5 * 1000 = 179.5 J/kg of refreezing energy
    z = np.array([[0.0], [1000.0]])
    tw = np.array([[-5.0], [-5.0]])
    seg_E, seg_warm, nseg = B.layer_energies(tw, z, np.array([1]))
    assert nseg[0] == 0 and not seg_warm[0, 0]
    assert seg_E[0, 0] == pytest.approx(179.5, abs=0.1)


def test_energy_splits_at_zero_crossing():
    z = np.array([[0.0], [1000.0]])
    tw = np.array([[-2.0], [2.0]])
    seg_E, seg_warm, nseg = B.layer_energies(tw, z, np.array([1]))
    assert nseg[0] == 1 and not seg_warm[0, 0] and seg_warm[1, 0]
    assert seg_E[0, 0] == pytest.approx(seg_E[1, 0], rel=1e-6)


def test_cold_column_is_snow():
    p, slr, _ = classify(lambda pr: -3 - (1013 - pr) * 0.045)
    assert p.argmax() == SN and p[SN] >= 95
    assert 12.0 <= slr <= 30.0


def test_warm_column_is_rain():
    p, _, _ = classify(lambda pr: 12 - (1013 - pr) * 0.045)
    assert p.argmax() == RA and p[RA] >= 95 and p[SN] == 0


def test_shallow_cold_layer_is_freezing_rain():
    """Regression for the Eq. 7 coefficient (-2.1, not -22.1): shallow sub-freezing layer under a warm nose."""
    p, _, _ = classify(lambda pr: -2 if pr >= 950 else 4 - (950 - pr) * 0.06)
    assert p.argmax() == FZRA and p[FZRA] >= 90 and p[PL] == 0


def test_deep_cold_layer_under_warm_nose_is_sleet():
    def t(pr):
        if pr >= 800: return -9.0
        if pr >= 700: return 2.0
        return 2.0 - (700 - pr) * 0.12
    p, _, _ = classify(t)
    assert p.argmax() == PL and p[PL] >= 50


def test_dry_column_gives_ice_probability_from_coldest_level():
    # no saturated level at all -> falls back to the coldest level for ProbIce
    p, _, diag = classify(lambda pr: -3 - (1013 - pr) * 0.045, rh=30.0)
    assert p.sum() > 0 and p.argmax() == SN


def test_below_ground_levels_do_not_contribute():
    # surface at 850 hPa (high terrain): 925/1000 levels are underground and must be clamped
    p, _, _ = classify(lambda pr: -5 - (850 - pr) * 0.05 if pr <= 850 else 8.0, ps=850.0)
    assert p.argmax() == SN


def test_probabilities_bounded_and_types_split_qpf():
    probs = np.array([[40], [30], [20], [10]], np.uint8)
    sh = type_shares(probs)
    assert sh.sum(axis=0)[0] == pytest.approx(1.0)
    assert sh[1, 0] == pytest.approx(0.3)
    assert type_shares(np.zeros((4, 1), np.uint8)).sum() == 0


def test_dominant_type_tie_breaks_by_impact():
    assert B.dominant_type(np.array([[50], [50], [0], [0]]))[0] == SN          # SN beats RA
    assert B.dominant_type(np.array([[0], [50], [50], [50]]))[0] == FZRA       # FZRA > PL > SN


def test_snowfall_and_slr_modes():
    q = np.array([10.0], np.float32)                                   # 10 mm liquid
    probs = np.array([[0], [100], [0], [0]], np.uint8)
    assert snowfall_step_cm(q, probs, np.array([12.5], np.float32))[0] == pytest.approx(12.5)
    _, slr, _ = classify(lambda pr: -3 - (1013 - pr) * 0.045, slr_mode="fixed:12.5")
    assert slr == pytest.approx(12.5)
    # Kuchera: colder column -> higher ratio, warmer -> lower, clamped to [5, 30]
    _, cold, _ = classify(lambda pr: -20 - (1013 - pr) * 0.03)
    _, warm, _ = classify(lambda pr: -0.5 - (1013 - pr) * 0.045)
    assert cold > warm and 5.0 <= warm and cold <= 30.0
