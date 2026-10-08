# Comparison with Suresh & Chernoff (arXiv:2310.00825)

**What this is.** A comparison between my cusp-beam work (`cuspbeam.py`) and the published formalism of Suresh & Chernoff, "Modeling the beam of gravitational radiation from a cosmic string loop". I used their *published formulas and test loops* only; I did not have their code. The derivations are mine and unchecked by others. Everything below was computed in a sandbox with direct numerical integrals as the reference.

## 1. Where the two approaches coincide (so it is not new)

| Item | Suresh & Chernoff | My work |
|---|---|---|
| Local expansion at a cusp | Cubic phase, linear prefactor, Airy-type integrals (their Sec. III C, Eq. 41, App. B) | Cubic phase with a prefactor to second order, Airy integrals |
| On-axis scaling | I1 ~ m^(-1/3), I2 ~ m^(-2/3), dP/dOmega ~ m^(-2/3) | Same; constant K = 0.7161 reproduced |
| Beam shape | More complicated than a filled cone; two prominent sub-peaks; asymptotically self-similar | Beam not brightest on axis; profile depends on n^(1/3) * theta |
| Acute crossing angle | "Extended emission region near the cusp" for an acute crossing (their Sec. VI B) | Beam area ~ 1/sin(psi) (quantitative version of the same statement) |
| Validity | Bound m_low from 4th/5th-order phase terms (Eq. 57); work exclusively at cubic order and quantify errors empirically | Same terms, used to derive a closed-form correction |

## 2. What I could test on their own loops (Turok loops, their Eq. 64)

On-axis, isolated cusps. Exact on-axis log-error of the leading law versus my parameter-free prediction:

| Loop (alpha, Phi) | m = 100 | m = 300 | m = 1000 | m = 3000 |
|---|---|---|---|---|
| (1/5, -pi/2), exact / predicted | -0.0375 / -0.0377 | -0.0180 / -0.0180 | -0.0081 / -0.0081 | -0.0039 / -0.0039 |

- The correction does not depend on the crossing angle: the acute-crossing loop (sin psi = 0.048) gave the same numbers, as the theory requires.
- A well-separated six-cusp loop (4/5, 2pi/5) agreed to four digits from m = 300; at m = 100 the exact value was -0.159 against -0.081 predicted (expansion parameter too large).
- Four cusps of the six-cusp loop (3/10, pi/4) had exact values of +0.49 and +0.25 at m = 100 and 300 against about -0.03 predicted: pseudocusp contamination, as their paper describes; the agreement was restored by m = 1000 to 3000.

Beam weight inside the model's lobe (exact / model, n = 1024 -> 4096):

| Cusp | sin psi | n = 1024 | n = 4096 |
|---|---|---|---|
| (1/5, -pi/2) | 0.868 | 0.990 | not run |
| (4/5, 2pi/5) | 0.818 | 0.984 | not run |
| (3/10, pi/4) | 0.262 | 0.699 | 1.207 |
| (1/5, 3pi/20) | 0.048 | 0.450 | 0.381 |

At sin psi of 0.8 and above the local model matches the exact weight to about 1.5% on an independent loop family. At sin psi = 0.26 it moves toward 1 as n grows (overshooting). At sin psi = 0.048 the beam arms are longer than about 1 rad, outside the small-angle model, and it fails.

## 3. Their validity bound versus actual accuracy

For 245 random-loop cusps, Eq. 57 gives m_low with a median of about 30 (16-84%: 7 to 260). The harmonic needed for the cubic law to reach 10% accuracy has a median of about 220 (7 times m_low), for 1% about 7,000 (240 times m_low), and for 0.1% about 2e5 (7,500 times m_low). On the Turok loops m_low is 0.24 to 2.4 while errors are still 4 to 16% at m = 100. The bound marks where the expansion is not diverging, not where it is accurate; the error then falls as m^(-2/3), so the harmonic for accuracy epsilon scales as epsilon^(-3/2).

## 4. What they have that I do not
Multipoint emission over the full sphere including pseudocusps and kinks (a general loop, not only the cusp direction), a complex-asymptotic treatment, an FFT-based exact calculation, and a validated treatment of Turok loops. My work covers only the cusp direction and its lobe.

## 5. What remains possibly new (not confirmed)
1. The closed-form n^(-2/3) correction to the on-axis amplitude in terms of four local quantities per curve (their paper works at cubic order and does not give it).
2. The explicit scaling of the beam's effective solid angle, (s_a s_b)^(1/3) / |sin psi|, and its independence of the detector noise curve.
3. The quantified validity limits above (the tested crossing-angle range, the gap between Eq. 57 and actual accuracy).
Items 1 and 2 may be routine to specialists; I have not been able to establish that.

## 6. Questions for the authors
Is the on-axis n^(-2/3) correction known to you? Does your multipoint code reproduce the effective-solid-angle scaling for acute crossing angles? Would a comparison on your loop set be useful?
