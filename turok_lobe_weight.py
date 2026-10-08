"""turok_lobe_weight.py -- reproduce the lobe-weight comparison on Turok loops (Suresh & Chernoff arXiv:2310.00825, Eq. 64, l = 2 pi).
Self-contained (numpy + scipy + cuspbeam.py + turok_onaxis_check.py).  Run:
    python examples/turok_lobe_weight.py            (n = 1024, about 1 minute)
    python examples/turok_lobe_weight.py 4096       (n = 4096, a few minutes)
For one cusp of each loop: exact weight W = n^(2/3) * integral (B_exact)^(3/2) dOmega, integrated on a grid of directions inside the model's own
lobe (out to 1.15 x the radius holding 99% of the model weight), divided by the same quantity from the local model in cuspbeam.py.
A ratio near 1 means the model describes the lobe; a nearby second beam or a very small crossing angle spoils it."""
import os
import sys
import time
import numpy as np
from math import pi

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, ".."))
import cuspbeam as cb
from turok_onaxis_check import turok_curves, find_cusps, derivs, side

TWO_PI = 2 * pi


def exact_B_on_directions(fa, fb, n, K, N, chunk=16):
    """B = n^(2/3) dP_n/dOmega in each direction of K (D x 3), direct numerical integrals."""
    s = np.arange(N) * TWO_PI / N
    a, b = fa(s), fb(s)
    cum = lambda x: np.vstack([np.zeros((1, 3)), np.cumsum((x[:-1] + x[1:]) / 2, axis=0)]) * (TWO_PI / N)
    Pa, Pb = cum(a), cum(b)
    ref = np.array([0.3, 0.5, 0.8])
    out = np.empty(len(K))
    for i in range(0, len(K), chunk):
        k = K[i:i + chunk]
        f1 = np.cross(k, ref); f1 /= np.linalg.norm(f1, axis=1, keepdims=True)
        f2 = np.cross(k, f1)
        Ia = a.T @ np.exp(1j * n * (s[:, None] - Pa @ k.T)) / N
        Ib = b.T @ np.exp(1j * n * (s[:, None] - Pb @ k.T)) / N
        comp = lambda I: (np.einsum('id,di->d', I, f1), np.einsum('id,di->d', I, f2))
        Ax, Ay = comp(Ia); Bx, By = comp(Ib)
        S = (abs(Ax) ** 2 + abs(Ay) ** 2) * (abs(Bx) ** 2 + abs(By) ** 2) + 4 * (Ax * np.conj(Ay)).imag * (Bx * np.conj(By)).imag
        out[i:i + chunk] = n ** (2 / 3) * 8 * pi * n ** 2 * S
    return out


def main(n=1024):
    cases = [("(1/5, -pi/2)", 0.2, -pi / 2), ("(4/5, 2pi/5)", 0.8, 2 * pi / 5), ("(3/10, pi/4)", 0.3, pi / 4), ("(1/5, 3pi/20)", 0.2, 3 * pi / 20)]
    t0 = time.time()
    for name, al, ph in cases:
        fa, fb = turok_curves(al, ph)
        u, v = find_cusps(fa, fb)[0]
        nh = fa(u)
        A, B = side(fa, u, nh), side(fb, v, nh)
        sa, sb = A[0], B[0]
        da, db = derivs(fa, u)[0], derivs(fb, v)[0]
        ta = da / sa; ta -= (ta @ nh) * nh; ta /= np.linalg.norm(ta); yax = np.cross(nh, ta)
        tb = db / sb; tb -= (tb @ nh) * nh; tb /= np.linalg.norm(tb)
        psi = np.arctan2(tb @ yax, tb @ ta)
        r99 = cb.lobe_radius(sa, sb, psi, n)
        rmax = min(1.15 * r99, 0.8)
        nth, nal = 40, 72
        th = np.linspace(0, rmax, nth); al_ = np.arange(nal) * TWO_PI / nal
        TH, AL = np.meshgrid(th[1:], al_, indexing="ij")
        Kdir = (np.sin(TH) * np.cos(AL))[..., None] * ta + (np.sin(TH) * np.sin(AL))[..., None] * yax + np.cos(TH)[..., None] * nh
        N = 1 << int(np.ceil(np.log2(max(65536, 32 * n))))
        B_ex = exact_B_on_directions(fa, fb, n, Kdir.reshape(-1, 3), N).reshape(nth - 1, nal)
        ring = (B_ex ** 1.5).mean(axis=1) * TWO_PI * np.sin(th[1:])
        W_exact = np.trapezoid(np.r_[0.0, ring], th) * n ** (2 / 3)
        W_model = cb.on_axis_B(sa, sb) ** 1.5 * cb.effective_solid_angle(sa, sb, psi, n, 1.5, theta_max=rmax) * n ** (2 / 3)
        print(f"Turok {name:14s} n={n}: s_a={sa:.2f} s_b={sb:.2f} sin psi={abs(np.sin(psi)):.3f}  lobe radius {r99:.3f} rad  "
              f"exact/model weight = {W_exact / W_model:.3f}   [{time.time() - t0:.0f}s]", flush=True)


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 1024)
