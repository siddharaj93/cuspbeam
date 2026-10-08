"""turok_onaxis_check.py -- reproduce the on-axis table for Turok loops (Suresh & Chernoff arXiv:2310.00825, Eq. 64, l = 2 pi).
Self-contained (numpy + scipy + cuspbeam.py).  Run:  python examples/turok_onaxis_check.py     (about 1 minute)
For every cusp it prints the exact log-error of the leading on-axis law, res_n = log(B_n / [K (s_a s_b)^(-2/3)]), from a direct numerical
integral, next to the parameter-free next-order prediction cb.onaxis_correction.  Cusps contaminated by pseudocusps (near-coincidences of the
curves close to the cusp direction) disagree at low n, as the paper describes."""
import os
import sys
import numpy as np
from math import pi
from scipy.optimize import least_squares

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import cuspbeam as cb

TWO_PI = 2 * pi


def turok_curves(alpha, Phi):
    r = 2 * np.sqrt(alpha * (1 - alpha))
    a = lambda u: np.stack([-(1 - alpha) * np.cos(u) - alpha * np.cos(3 * u), (1 - alpha) * np.sin(u) + alpha * np.sin(3 * u), r * np.sin(np.asarray(u, float))], -1)
    b = lambda u: np.stack([np.cos(u), np.sin(u) * np.cos(Phi), np.sin(u) * np.sin(Phi)], -1)
    return a, b


def find_cusps(a, b, grid=720):
    g = np.linspace(0, TWO_PI, grid, endpoint=False)
    A, B = a(g), b(g)
    d2 = 2 - 2 * A @ B.T
    mins = np.ones_like(d2, dtype=bool)
    for di in (-1, 0, 1):
        for dj in (-1, 0, 1):
            if di or dj:
                mins &= d2 <= np.roll(np.roll(d2, di, 0), dj, 1)
    cusps = []
    for i, j in np.argwhere(mins & (d2 < 0.05 ** 2)):
        sol = least_squares(lambda x: a(x[0]) - b(x[1]), [g[i], g[j]], xtol=1e-15, ftol=1e-15, gtol=1e-15)
        if np.linalg.norm(sol.fun) < 1e-9:
            u, v = sol.x % TWO_PI
            if all(abs((u - u2 + pi) % TWO_PI - pi) > 1e-6 or abs((v - v2 + pi) % TWO_PI - pi) > 1e-6 for u2, v2 in cusps):
                cusps.append((u, v))
    return cusps


def derivs(fun, x, h=0.02, npts=9):
    from math import factorial
    k = np.arange(npts) - npts // 2
    f = fun(x + k * h)
    V = np.vander(k, increasing=True).T.astype(float)
    return tuple((np.linalg.solve(V, np.eye(npts)[m] * factorial(m)) / h ** m) @ f for m in (1, 2, 3))


def side(fun, s0, nh):
    d1, d2, d3 = derivs(fun, s0)
    s = np.linalg.norm(d1); t = d1 / s; m = np.cross(nh, t)
    return (s, float(d2 @ t), float(d2 @ m), float(d3 @ t))


def exact_B(fa, fb, nh, ns, N=1 << 18):
    """B_n = n^(2/3) dP_n/dOmega on the cusp axis from direct numerical integrals (Blanco-Pillado & Olum Eq. 45-50, l = 2 pi)."""
    s = np.arange(N) * TWO_PI / N
    e1 = np.cross(nh, [0.3, 0.5, 0.8]); e1 /= np.linalg.norm(e1); e2 = np.cross(nh, e1)
    I = []
    for fun in (fa, fb):
        a = fun(s)
        Pz = np.concatenate([[0], np.cumsum((a[:-1] + a[1:]) / 2 @ nh)]) * (TWO_PI / N)
        I.append({n: (a * np.exp(1j * n * (s - Pz))[:, None]).sum(0) / N for n in ns})
    out = {}
    for n in ns:
        A = np.array([I[0][n] @ e1, I[0][n] @ e2]); B = np.array([I[1][n] @ e1, I[1][n] @ e2])
        im = lambda V: (V[0] * np.conj(V[1])).imag
        S = (abs(A) ** 2).sum() * (abs(B) ** 2).sum() + 4 * im(A) * im(B)
        out[n] = n ** (2 / 3) * 8 * pi * n ** 2 * S
    return out


def main(ns=(100, 300, 1000, 3000)):
    loops = {"(alpha, Phi) = (1/5, -pi/2): two well-separated cusps": (0.2, -pi / 2),
             "(1/5, 3pi/20): acute crossing": (0.2, 3 * pi / 20),
             "(3/10, pi/4): six cusps": (0.3, pi / 4),
             "(4/5, 2pi/5): six well-separated cusps": (0.8, 2 * pi / 5)}
    for name, (al, ph) in loops.items():
        print(f"\nTurok loop {name}")
        fa, fb = turok_curves(al, ph)
        for i, (u, v) in enumerate(find_cusps(fa, fb)):
            nh = fa(u)
            A, B = side(fa, u, nh), side(fb, v, nh)
            da, db = derivs(fa, u)[0], derivs(fb, v)[0]
            sinpsi = abs(np.cross(da / A[0], db / B[0]) @ nh)
            Bn = exact_B(fa, fb, nh, list(ns))
            print(f"  cusp {i}: s_a = {A[0]:.3f}  s_b = {B[0]:.3f}  sin psi = {sinpsi:.3f}")
            print("      n:         " + "".join(f"{n:>10d}" for n in ns))
            print("      exact res: " + "".join(f"{np.log(Bn[n] / cb.on_axis_B(A[0], B[0])):>+10.4f}" for n in ns))
            print("      predicted: " + "".join(f"{cb.onaxis_correction(n, A, B):>+10.4f}" for n in ns))


if __name__ == "__main__":
    main()
