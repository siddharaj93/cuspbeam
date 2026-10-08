# Cusp beams of cosmic-string loops: a verified next-order correction and a local beam model

**Status: working draft. Numerical results on smooth, centrally symmetric loops only. The derivations are the author's own and have not been checked by anyone else. A limited literature search found neither result, which is weak evidence of novelty, not proof. Both results are small refinements, not discoveries; relevance to real searches is modest (section 5).**

## 1. Summary

1. **Next-order correction (isolated cusps).** The on-axis cusp beam amplitude B_n = n^(2/3) dP_n/dΩ = K (s_a s_b)^(-2/3), K = 0.7161 (the on-axis limit of Blanco-Pillado & Olum, arXiv:1709.02693, App. A) has a closed-form, parameter-free correction of relative order n^(-2/3), verified on fresh seeds.
2. **Local beam model and detection weight.** A cubic-phase (Airy) model of each curve near the cusp reproduces the full off-axis beam. A toy detection weight W = ∫ (dP/dΩ)^(3/2) dΩ scales as (s_a s_b)^(-2/3)/sin ψ in the model (ψ = crossing angle of the two tangent directions); the model weight matches exact numerics to a few percent at n = 8192 for isolated cusps (36 cusps over three runs, four loop families; the 24 selected with the neighbour rule all lie within [0.88, 1.15]).
3. **Where both fail:** overlapping beams (neighbouring cusps, pseudocusp-like near-coincidences) and wigglier loops.

## 2. Next-order correction

Expanding the phase to fifth order and the amplitude to third order about the cusp, the Airy-type moments J_2 and J_5 vanish exactly, the cross-polarization term is fourth order, and the amplitude stays purely imaginary, so the first correction to the power is relative order n^(-2/3). To leading order
res_n ≈ 1.59 n^(-2/3) Σ_{a,b} [ −0.27 w/s^(7/3) + 0.48 τ²/s^(10/3) − 0.55 k²/s^(10/3) − 0.55 s^(2/3) ],
with s = |a'|, τ = a''·t, k = a''·m, w = a'''·t (t = a'/s, m = n̂×t). Coefficients are rounded; the code uses exact constants.

Tests (random odd-harmonic loops, direct integration, "clean" cusps: no other near-approach of either curve within 0.8 rad and, from the later runs, no other exact cusp within 0.5 rad):

| Seed set | Cusps | observed/predicted | median error at n=2048: leading law → with correction |
|---|---|---|---|
| 1000+ | 76 | 1.00 | reduced 165× |
| 6000+ | 64 | 1.00 | reduced 128× |
| 7000+ | 67 | 1.00 | reduced 157× |

The correction is stable to sampling resolution once the local derivatives are converged (an early 5-point stencil inflated the leftover). The leftover error falls with a log–log slope between −1.37 and −1.7 across samples; the next-order prediction (−1.33) is neither confirmed nor excluded. On wigglier loops (8–16 harmonics, 57–92 cusps) the single-cusp formula removes little of the error (reduction 1–2×); a first multipoint extension (adding an Airy term for each other near-stationary point) lowered the median error from 0.54 to 0.46 only.

## 3. Local beam model and detection weight

Each curve confines the beam to a strip along its tangent, of half-width about (s/n)^(1/3); the beam is the parallelogram where the strips cross, with area ∝ (s_a s_b)^(1/3)/sin ψ n^(-2/3). The model reproduces the published on-axis value exactly (1.0000 in all tests). Fitting the model weight over 96 combinations of (s_a, s_b, ψ): exponent −0.659 on s_a s_b (predicted −0.667), −1.044 on sin ψ (predicted −1), no s_a/s_b dependence, rms residual 1.7%.

Against exact numerics (weights integrated inside the model's own lobe; pass criteria fixed beforehand):

| Run | Cusps | exact/model at n = 8192 |
|---|---|---|
| Seeds 2000+ (first run) | 12 | median 1.01; 11 within 10%; one outlier (2.27) had another exact cusp 0.099 rad away |
| Seeds 4000+ (neighbour rule, fresh) | 12 | median 1.01; all 12 within [0.98, 1.15]; 11 within 4% |
| Seeds 5000+ (neighbour rule, fresh, three cusps from each of four loop families) | 12 | median 1.01; all 12 within [0.88, 1.04]; 11 within 4%; the smallest-slope cusp (s_a = 0.22) stays at 0.88 for n = 2048 and 8192 instead of converging |

The second run's cusps all came from the simplest loop family (two harmonics); the third run drew three cusps from each of four families (two to four harmonics, spectral decay 0.5 to 1.5). A cusp with a pseudocusp-like near-coincidence 0.26 rad from its axis went 1.41 → 1.74 → 1.15 for n = 512, 2048, 8192, as that feature left the shrinking lobe. The beam is not brightest on the axis: for one cusp, directions about 0.6 widths away reach 2.5× the on-axis value, which the model reproduces.

**Detector-weighted version (model only).** Weighting the burst by an approximate noise curve (SNR(Omega)^2 = integral f^(-8/3) Pr(Omega; n(f)) / S_n(f) df; detection volume ~ integral SNR^3 dOmega; Euclidean, no antenna pattern or redshift) leaves the geometric dependence of the effective beam solid angle unchanged: Omega_eff ~ (s_a s_b)^(+0.333) (sin psi)^(-1.000) for aLIGO-design and LISA analytic noise curves at reference harmonic n = 10^8, with exponents +0.335 and -1.014 at n = 10^5 (rms of the fit at most 0.004 in log). So the crossing-angle factor is universal, not detector-specific; the detector enters only through the overall constant. Omega_eff scales as n^(-2/3) (ratio 99.8 between n = 10^5 and 10^8), and for unit slopes and perpendicular tangents it equals about 3.5 times the solid angle pi n^(-2/3) of a cone of half-angle n^(-1/3); that ratio depends on how the cone angle is defined (a factor of order 2 in the angle changes it by 4), so it is not a prediction for existing forecasts. The noise curves are analytic approximations written from memory and should be checked against the official curves. This is a property of the model, which was validated against exact numerics only up to n = 8192; the detector bands correspond to far higher harmonics.

**Ensemble (toy).** For 666 cusps of random odd-harmonic loops (highest harmonic 3, 7, 11): the 1/sin ψ factor raises the mean weight by 1.44–1.53×; relative to unit slopes the mean slope-weight is 0.59, 0.38, 0.28; mean amplitude factor ⟨(s_a s_b)^(-1/3)⟩ = 0.73, 0.57, 0.50, against published means 0.668, 0.513, 0.446 (Pazouli et al., arXiv:2008.13693) for a different loop generator. The top 5% of cusps carry 18–23% of the weight.

## 4. Literature context

- On-axis amplitude ∝ (|a''||b''|)^(-1/3) and a beaming angle set by the smaller second derivative are standard (Damour & Vilenkin; Pazouli et al.).
- Burst-rate formulas typically carry O(1) "ignorance constants" for beam geometry.
- Multipoint (nearly stationary points) treatment: Suresh & Chernoff, arXiv:2310.00825; pseudocusps: Stott, Elghozi & Sakellariadou, arXiv:1612.07599.
- I did not find the closed-form n^(-2/3) correction, nor a crossing-angle dependence of the beam area, in the papers I could read. I did not read the full text of the cusp/pseudocusp thesis (arXiv:2108.08242), which derives the waveform in arbitrary observation directions.

## 4b. Comparison with Suresh & Chernoff (added after reading the full paper)

Their local cubic-phase/Airy expansion is the same method as my local model; they also state the m^(-2/3) on-axis scaling, a self-similar beam with two prominent sub-peaks (so "the beam is not brightest on the axis" is known), and an extended emission region for acute crossing angles. On their Turok loops my parameter-free on-axis correction agreed with exact numerics to four digits for isolated cusps from m = 300 (for example -0.0377 predicted against -0.0375 exact at m = 100 for loop (1/5, -pi/2)), and the beam weight agreed to 1.6% and 1.0% at sin psi = 0.82 and 0.87; at sin psi = 0.26 the exact/model weight went 0.70 (n = 1024) to 1.21 (n = 4096), and the model fails at sin psi = 0.048. Their validity bound (Eq. 57) is a very loose criterion: for random-loop cusps its median is about 30, against about 220, 7,000 and 2e5 for 10%, 1% and 0.1% accuracy. What may remain new is only the closed-form next-order correction and the explicit scaling of the beam area with slopes and crossing angle; both may be routine to specialists. Details: comparison_with_suresh_chernoff.md.

## 5. Relevance and limitations

- The correction is 1–7% at n = 256–2048 and smaller at higher n; the angle factor changes the toy ensemble average by about 1.5×. Loop population, string tension, and cusp counts matter far more for detection.
- Useful mainly as a verified benchmark for cusp-waveform code, a cheap beam model for burst Monte Carlos (milliseconds against minutes per cusp), and a map of where single-cusp formulas fail.
- Detector-relevant harmonic numbers are much higher than those tested (an estimate: n of order 10^6 or more). The n^(-2/3) correction is then below about 10^-4 times a geometry factor, so the beam model is the only part of this work that does not shrink with n.
- Limits: smooth centrally symmetric loops; isolated cusps; small samples (36 cusps over three runs for the weight check, loops with at most four harmonics); the clean-cusp rules are not complete (a near-coincidence beside the cusp is missed); detection weight is a toy (Euclidean volume, fixed threshold, no redshift or detector response).

## 6. Next steps

Extend the model to overlapping beams; test richer loop families and loops with kinks; compute rate multipliers with a real loop population; read the Pazouli thesis chapter on rates; ask a specialist whether either result is already known.
