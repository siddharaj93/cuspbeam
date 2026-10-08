"""cuspbeam -- a cheap, validated local model of the gravitational-wave beam of a cosmic-string cusp.

Single file, needs only numpy and scipy.  Status: research code by a non-specialist, derivations NOT checked by others.
See README.md for the validation record and, importantly, the limits.

GEOMETRY.  A cusp is where the two tangent-vector curves a(u) and b(v) of a loop meet on the unit sphere.  At the cusp:
    s_a = |da/du|, s_b = |db/dv|     slopes ("sharpness") of the two curves
    psi                              crossing angle between the two tangent directions (sin psi = |a' x b'| / (s_a s_b))
    n                                harmonic number (frequency f = 2 n / l for a loop of length l)
The beam is described in the plane perpendicular to the cusp direction: theta = angle from the cusp direction, alpha = azimuth
measured from the tangent direction of curve a, tangent of curve b at angle psi.

NORMALISATION.  B = n^(2/3) dP_n/dOmega  (units of G mu^2, loop length 2 pi, as in Blanco-Pillado & Olum arXiv:1709.02693).
On axis B = K (s_a s_b)^(-2/3), K = Gamma(2/3)^4 6^(8/3) / (18 pi^3) = 0.7161 (their appendix A, theta -> 0).

WHAT IT IS:  leading-order (cubic phase) Airy model of each curve near the cusp, accurate for ISOLATED cusps of smooth loops at n >~ 500.
WHAT IT IS NOT:  it knows nothing about other beams (neighbouring cusps, pseudocusps), kinks, or wiggly loops.  Use is_isolated() and
check_validity() before trusting a result."""
from __future__ import annotations
import numpy as np
from math import pi, gamma
from scipy.special import airy
from scipy.integrate import quad

__all__ = ["K_ONAXIS", "on_axis_B", "beam_profile", "relative_profile", "lobe_radius", "effective_solid_angle",
           "effective_solid_angle_scaling", "detection_weight", "onaxis_correction", "is_isolated", "check_validity",
           "S_ligo", "S_lisa", "S_lisa_full", "detector_weighted_solid_angle", "angle_factor_ratio"]

K_ONAXIS = gamma(2 / 3) ** 4 * 6 ** (8 / 3) / (18 * pi ** 3)          # = 0.7161
# n^(2/3) * (effective solid angle) for the reference cusp s_a = s_b = 1, psi = pi/2, single harmonic (filled in below, verified by tests)
OMEGA32_REF = 10.868044   # weight power 3/2:  integral (P/P_axis)^(3/2) dOmega  (detection volume); re-derived by test_reference_constants
OMEGA1_REF = 9.404667    # weight power 1:    integral (P/P_axis) dOmega         (energy)


# ------------------------------------------------------------------ local beam model
def on_axis_B(sa, sb):
    """B = n^(2/3) dP_n/dOmega on the cusp axis."""
    return K_ONAXIS * (np.asarray(sa) * np.asarray(sb)) ** (-2 / 3)


def _side_vector(s, t, n, th, al):
    nh = np.array([0, 0, 1.0])
    e = np.stack([np.cos(al), np.sin(al), np.zeros_like(al)], -1)
    c = e @ t
    mu = (2 * s / n) ** (1 / 3)
    x = n * th ** 2 * (1 - c ** 2) / (2 * s) * mu
    Ai, Aip, _, _ = airy(x)
    V0 = nh * (1 - th ** 2 * c ** 2 / 2)[..., None] + (th * c)[..., None] * t
    V1 = t - (th * c)[..., None] * nh
    V2 = -nh / 2
    return (mu / s) * (V0 * Ai[..., None] - 1j * mu * V1 * Aip[..., None] - (mu ** 2) * V2 * (x * Ai)[..., None])


def beam_profile(sa, sb, psi, n, theta, alpha):
    """B(theta, alpha) = n^(2/3) dP_n/dOmega near the cusp axis (arrays allowed for theta, alpha; theta in radians, small)."""
    th, al = np.broadcast_arrays(np.asarray(theta, float), np.asarray(alpha, float))
    ta, tb = np.array([1.0, 0, 0]), np.array([np.cos(psi), np.sin(psi), 0])
    f1 = np.stack([np.cos(th) * np.cos(al), np.cos(th) * np.sin(al), -np.sin(th)], -1)
    f2 = np.stack([-np.sin(al), np.cos(al), np.zeros_like(al)], -1)
    comps = lambda I: (np.sum(I * f1, -1), np.sum(I * f2, -1))
    Ax, Ay = comps(_side_vector(sa, ta, n, th, al))
    Bx, By = comps(_side_vector(sb, tb, n, th, al))
    S = (abs(Ax) ** 2 + abs(Ay) ** 2) * (abs(Bx) ** 2 + abs(By) ** 2) + 4 * (Ax * np.conj(Ay)).imag * (Bx * np.conj(By)).imag
    return n ** (2 / 3) * 8 * pi * n ** 2 * S


def relative_profile(sa, sb, psi, n, theta, alpha):
    """B(theta, alpha) / B(axis).  Note it can exceed 1: the beam is not always brightest on the axis."""
    return beam_profile(sa, sb, psi, n, theta, alpha) / on_axis_B(sa, sb)


def _grid(sa, sb, psi, n, nth=400, nal=360, theta_max=None):
    w = (max(sa, sb) / n) ** (1 / 3)
    thmax = min(10.0 * w / max(abs(np.sin(psi)), 0.05), 1.0) if theta_max is None else float(theta_max)
    th = np.concatenate([[0.0], np.geomspace(1e-3 * thmax, thmax, nth)])
    al = np.linspace(0, 2 * pi, nal, endpoint=False)
    return th, al, np.meshgrid(th, al, indexing="ij")


def lobe_radius(sa, sb, psi, n, fraction=0.99):
    """Radius (rad) containing `fraction` of the model's detection weight; use it to decide how far away another beam must be."""
    th, al, (TH, AL) = _grid(sa, sb, psi, n)
    rp = relative_profile(sa, sb, psi, n, TH, AL) ** 1.5
    ring = rp.mean(axis=1) * 2 * pi * np.sin(th)
    cum = np.concatenate([[0.0], np.cumsum((ring[1:] + ring[:-1]) / 2 * np.diff(th))])
    return float(th[np.searchsorted(cum, fraction * cum[-1])])


def effective_solid_angle(sa, sb, psi, n, power=1.5, theta_max=None):
    """Numerical integral of (B/B_axis)^power over the beam (steradians).  power = 1.5: detection weight (Euclidean volume ~ amplitude^3);
    power = 1: energy.  Scales as n^(-2/3).  theta_max (rad) restricts the integral to a cone around the cusp direction."""
    th, al, (TH, AL) = _grid(sa, sb, psi, n, theta_max=theta_max)
    ring = (relative_profile(sa, sb, psi, n, TH, AL) ** power).mean(axis=1) * 2 * pi * np.sin(th)
    return float(np.trapezoid(ring, th))


def effective_solid_angle_scaling(sa, sb, psi, n, power=1.5):
    """Closed form  Omega_ref * (s_a s_b)^(1/3) / |sin psi| * n^(-2/3)   (asymptotic scaling; verified numerically, see tests)."""
    ref = _ref(power)
    return ref * (np.asarray(sa) * np.asarray(sb)) ** (1 / 3) / np.abs(np.sin(psi)) * np.asarray(n, float) ** (-2 / 3)


def detection_weight(sa, sb, psi):
    """Toy detection weight of one cusp, n^(2/3) * integral B^(3/2) dOmega  = K^(3/2) Omega32_ref (s_a s_b)^(-2/3) / |sin psi|.
    Euclidean volume, fixed threshold, one harmonic.  Only ratios between cusps are meaningful."""
    return K_ONAXIS ** 1.5 * _ref(1.5) * (np.asarray(sa) * np.asarray(sb)) ** (-2 / 3) / np.abs(np.sin(psi))


def _ref(power):
    global OMEGA32_REF, OMEGA1_REF
    if OMEGA32_REF is None or OMEGA1_REF is None:                       # (only if the constants were cleared)
        n = 1e8
        OMEGA32_REF = effective_solid_angle(1.0, 1.0, pi / 2, n, 1.5) * n ** (2 / 3)
        OMEGA1_REF = effective_solid_angle(1.0, 1.0, pi / 2, n, 1.0) * n ** (2 / 3)
    return OMEGA32_REF if power == 1.5 else OMEGA1_REF


# ------------------------------------------------------------------ next-order correction (isolated cusp, on axis)
def _J(m):
    return 3 ** ((m - 2) / 3) * gamma((m + 1) / 3) * (np.exp(1j * pi * (m + 1) / 6) + (-1) ** m * np.exp(-1j * pi * (m + 1) / 6))


def _side_log(n, s, tau, k, w):
    lam = (2 / (n * s ** 2)) ** (1 / 3)
    a1, a2, e1 = tau * lam / (2 * s), w * lam ** 2 / (6 * s), lam * tau / (4 * s)
    e2 = lam ** 2 * (4 * s * w + 3 * (tau ** 2 + k ** 2 + s ** 4)) / (60 * s ** 2)
    V = _J(1) + a1 * _J(2) + a2 * _J(3) + 1j * e1 * _J(5) + 1j * e2 * _J(6) + 1j * a1 * e1 * _J(6) - 0.5 * e1 ** 2 * _J(9)
    return float(np.log(abs(V) ** 2 / abs(_J(1)) ** 2))


def onaxis_correction(n, side_a, side_b):
    """Predicted log(B_n / [K (s_a s_b)^(-2/3)]) for an isolated cusp.  Each side = (s, tau, k, w) with s = |a'|, tau = a''.t,
    k = a''.m, w = a'''.t (t = a'/s, m = nhat x t).  Parameter-free; relative order n^(-2/3).  Trust it only if |result| < 0.3."""
    return _side_log(n, *side_a) + _side_log(n, *side_b)


# ------------------------------------------------------------------ validity helpers
def is_isolated(cusp_direction, other_cusp_directions, min_sep=0.5):
    """True if no other exact cusp lies within min_sep radians of this cusp's direction.  NOT sufficient on its own: a near-coincidence
    of the two curves a couple of beam widths away (pseudocusp-like) is not detected; see README."""
    d = np.asarray(other_cusp_directions, float).reshape(-1, 3)
    c = np.asarray(cusp_direction, float)
    ang = np.arccos(np.clip(d @ c / (np.linalg.norm(d, axis=1) * np.linalg.norm(c)), -1, 1))
    ang = ang[ang > 1e-6]
    return bool(len(ang) == 0 or ang.min() >= min_sep)


def check_validity(n, sa, sb, psi, correction=None):
    """List of human-readable warnings; an empty list means 'no warning from these checks', not 'verified'."""
    w = []
    if n < 300:
        w.append("n < 300: model corrections are tens of percent; use exact numerics")
    if abs(np.sin(psi)) < 0.1:
        w.append("|sin psi| < 0.1: beam arms longer than ~1 rad, outside the small-angle model (exact/model weight 0.45 -> 0.38 for n = 1024 -> 4096 at sin psi = 0.048)")
    elif abs(np.sin(psi)) < 0.4:
        w.append("|sin psi| < 0.4: weakly tested (Turok loop at sin psi = 0.26: exact/model weight 0.70 at n = 1024, 1.21 at n = 4096)")
    if min(sa, sb) < 0.3:
        w.append("min slope < 0.3: tested only down to s = 0.22, where the model sat about 12% low at n = 8192")
    if max(sa, sb) > 12:
        w.append("max slope > 12: outside the tested range")
    if correction is not None and abs(correction) > 0.3:
        w.append("predicted next-order correction > 0.3: the truncated expansion is not reliable")
    return w


# ------------------------------------------------------------------ detector noise (analytic approximations; CHECK against official curves)
def S_ligo(f):
    """aLIGO design PSD, the common analytic fit (215 Hz)."""
    x = np.asarray(f, float) / 215.0
    return 1e-49 * (x ** -4.14 - 5 * x ** -2 + 111 * (1 - x ** 2 + 0.5 * x ** 4) / (1 + 0.5 * x ** 2))


def S_lisa(f):
    """LISA sky-averaged instrument noise, Robson, Cornish & Liu 2019 (analytic)."""
    f = np.asarray(f, float)
    L, fs = 2.5e9, 19.09e-3
    Poms = (1.5e-11) ** 2 * (1 + (2e-3 / f) ** 4)
    Pacc = (3e-15) ** 2 * (1 + (0.4e-3 / f) ** 2) * (1 + (f / 8e-3) ** 4)
    return (10 / (3 * L ** 2)) * (Poms + 2 * (1 + np.cos(f / fs) ** 2) * Pacc / (2 * pi * f) ** 4) * (1 + 0.6 * (f / fs) ** 2)


def S_lisa_full(f):
    """S_lisa + Galactic confusion noise for a 4-year mission (Robson et al. 2019: A = 9e-45, alpha = 0.138, beta = -221,
    kappa = 521, gamma = 1680, f_k = 0.00113 Hz)."""
    f = np.asarray(f, float)
    A, al, be, ka, ga, fk = 9e-45, 0.138, -221.0, 521.0, 1680.0, 0.00113
    return S_lisa(f) + A * f ** (-7 / 3) * np.exp(-f ** al + be * f * np.sin(ka * f)) * (1 + np.tanh(ga * (fk - f)))


_DETECTORS = {"LIGO": (S_ligo, 15.0, 1500.0, 100.0), "LISA": (S_lisa, 1e-4, 5e-2, 3e-3), "LISA+confusion": (S_lisa_full, 1e-4, 5e-2, 3e-3)}


def detector_weighted_solid_angle(sa, sb, psi, detector="LIGO", n_ref=1e8, nf=12):
    """Effective solid angle of the cusp burst weighted by a detector noise curve (cusp spectrum f^(-4/3), Euclidean, SNR^3 weighting).
    n = n_ref * f / f_ref at the detector's reference frequency (LIGO 100 Hz, LISA 3 mHz).  Scales as n_ref^(-2/3) (s_a s_b)^(1/3) / |sin psi|."""
    S, fmin, fmax, fref = _DETECTORS[detector]
    f = np.geomspace(fmin, fmax, nf)
    n = n_ref * f / fref
    wgt = f ** (-8 / 3) / S(f) * f
    th, al, (TH, AL) = _grid(sa, sb, psi, float(n.min()))
    X = sum(wi * relative_profile(sa, sb, psi, ni, TH, AL) for ni, wi in zip(n, wgt))
    ring = ((X / wgt.sum()) ** 1.5).mean(axis=1) * 2 * pi * np.sin(th)
    return float(np.trapezoid(ring, th))


def angle_factor_ratio(sa, sb, psi, floor=0.1):
    """Mean toy weight WITH the 1/sin(psi) factor divided by the mean WITHOUT it, over a list of cusps (floor caps small sin psi)."""
    sa, sb, psi = map(np.asarray, (sa, sb, psi))
    base = (sa * sb) ** (-2 / 3)
    return float((base / np.maximum(np.abs(np.sin(psi)), floor)).mean() / base.mean())
