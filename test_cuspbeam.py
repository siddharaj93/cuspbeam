"""Tests for cuspbeam.py.  Run:  python -m pytest test_cuspbeam.py -v     (about 20 seconds)
These check internal consistency and reproduce the numbers in README.md; they do NOT prove the physics (that rests on the comparison with
exact numerical integrals described in the README)."""
import numpy as np
from math import pi
import cuspbeam as cb

GEOMS = [(1, 1, pi / 2), (2, 1, 1.0), (4, 0.5, 0.5), (0.5, 3, 2.2)]


def test_K_constant():
    assert abs(cb.K_ONAXIS - 0.7161) < 2e-4


def test_profile_on_axis_equals_published_law():
    for g in GEOMS:
        v = cb.relative_profile(*g, 512, np.array([0.0]), np.array([0.0]))[0]
        assert abs(v - 1) < 1e-9


def test_reference_constants():
    n = 1e8
    assert abs(cb.effective_solid_angle(1, 1, pi / 2, n, 1.5) * n ** (2 / 3) / cb.OMEGA32_REF - 1) < 2e-3
    assert abs(cb.effective_solid_angle(1, 1, pi / 2, n, 1.0) * n ** (2 / 3) / cb.OMEGA1_REF - 1) < 2e-3


def test_closed_form_matches_numerical_integral():
    for g in GEOMS:
        num, cf = cb.effective_solid_angle(*g, 1e8), cb.effective_solid_angle_scaling(*g, 1e8)
        assert abs(num / cf - 1) < 2e-3, g


def test_solid_angle_scales_as_n_to_minus_two_thirds():
    r = cb.effective_solid_angle(2, 1, 1.0, 1e6) / cb.effective_solid_angle(2, 1, 1.0, 1e8)
    assert abs(r / (1e8 / 1e6) ** (2 / 3) - 1) < 5e-3          # ratio of n's is 100, so the solid angle ratio is 100^(2/3) = 21.5


def test_symmetric_in_the_two_curves():
    assert abs(cb.effective_solid_angle(3, 0.7, 0.9, 1e7) / cb.effective_solid_angle(0.7, 3, 0.9, 1e7) - 1) < 2e-3


def test_detection_weight_scaling():
    assert abs(cb.detection_weight(2, 2, 1.0) / cb.detection_weight(1, 1, 1.0) - 4 ** (-2 / 3)) < 1e-12
    assert abs(cb.detection_weight(1, 1, 0.5) / cb.detection_weight(1, 1, pi / 2) - 1 / np.sin(0.5)) < 1e-12


def test_noise_curves():
    f = np.geomspace(10, 3000, 600)
    assert 150 < f[np.argmin(cb.S_ligo(f))] < 300
    g = np.geomspace(1e-4, 0.1, 600)
    assert 3e-3 < g[np.argmin(cb.S_lisa(g))] < 1.5e-2
    assert np.all(cb.S_lisa_full(g) >= cb.S_lisa(g))


def test_detector_weighted_scaling_matches_geometry_factor():
    for det in ("LIGO", "LISA+confusion"):
        a = cb.detector_weighted_solid_angle(2, 1, 1.0, det, 1e8)
        b = cb.detector_weighted_solid_angle(1, 1, pi / 2, det, 1e8)
        assert abs(a / b / (2 ** (1 / 3) / np.sin(1.0)) - 1) < 5e-3


def test_correction_scales_as_n_to_minus_two_thirds_and_vanishes():
    a, b = (2.0, 1.5, 0.5, 2.0), (1.2, -0.8, 0.4, 1.0)
    c1, c2 = cb.onaxis_correction(1e6, a, b), cb.onaxis_correction(8e6, a, b)
    assert abs(c1 / c2 / 4 - 1) < 0.05 and abs(c2) < abs(c1) < 0.01


def test_is_isolated():
    c = np.array([0, 0, 1.0])
    others = np.array([[0, 0, 1.0], [0.2, 0, 1.0], [0, 1, 0.0]])
    assert not cb.is_isolated(c, others, 0.5)
    assert cb.is_isolated(c, others[[0, 2]], 0.5)


def test_validity_warnings():
    assert cb.check_validity(2000, 2.0, 1.5, 1.0) == []
    assert len(cb.check_validity(100, 0.2, 20.0, 0.1, correction=0.5)) == 5
