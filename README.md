# cuspbeam

A small, tested local model of the gravitational-wave beam of a cosmic-string cusp, with a closed-form next-order correction to the on-axis amplitude and a detector-weighted effective beam solid angle.

**Status.** Research code by a non-specialist. The derivations have not been checked by anyone else, and I could not confirm whether the results are already known. Treat it as a candidate benchmark and cheap approximation, not as a validated waveform model. It is single-file (`cuspbeam.py`), needs only `numpy` and `scipy`, and runs in milliseconds per cusp.

## Repository contents

|File|What it is|
|-|-|
|`cuspbeam.py`|the module (numpy + scipy only)|
|`test\\\_cuspbeam.py`|12 self-consistency tests: `python -m pytest test\\\_cuspbeam.py -v` (about 15 seconds; `pip install pytest` first)|
|`examples/turok\\\_onaxis\\\_check.py`|reproduces the on-axis table for Turok loops (Suresh \& Chernoff, arXiv:2310.00825, Eq. 64); about 1 minute|
|`examples/turok\\\_lobe\\\_weight.py`|reproduces the exact/model lobe-weight comparison on those loops; about 1.5 minutes for n = 1024|
|`docs/comparison\\\_with\\\_suresh\\\_chernoff.md`|what coincides with their paper, what was tested on their loops, what may remain new|
|`docs/writeup.md`|the longer write-up with all results and limits|

The random-loop validation behind the numbers below used scripts that depend on a separate loop-generation package that is not part of this repository, so those runs cannot be reproduced from here. The Turok-loop examples are self-contained.

### Expected output of the examples (for the loop (1/5, -pi/2), n = 100, 300, 1000, 3000)

```
exact res:    -0.0375   -0.0180   -0.0081   -0.0039
predicted:    -0.0377   -0.0181   -0.0081   -0.0039
```

and, for `turok\\\_lobe\\\_weight.py` at n = 1024, exact/model weights of 0.990, 0.984, 0.700 and 0.450 for the loops (1/5, -pi/2), (4/5, 2pi/5), (3/10, pi/4) and (1/5, 3pi/20).

## Licence

MIT (see `LICENSE`). Check and edit the name and year.

## What it computes

For one cusp with tangent-curve slopes `sa`, `sb` and crossing angle `psi` (sin psi = |a' x b'| / (sa sb)), at harmonic number `n`:

|Function|Returns|
|-|-|
|`on\\\_axis\\\_B(sa, sb)`|B = n^(2/3) dP/dOmega on axis = K (sa sb)^(-2/3), K = 0.7161 (the theta -> 0 limit of Blanco-Pillado \& Olum, arXiv:1709.02693, App. A)|
|`beam\\\_profile(sa, sb, psi, n, theta, alpha)`|B off axis, from a cubic-phase Airy model of each curve|
|`relative\\\_profile(...)`|B / B\_axis (it can exceed 1: the beam is not always brightest on axis)|
|`effective\\\_solid\\\_angle(sa, sb, psi, n, power=1.5)`|integral of (B/B\_axis)^power over the beam; power 1.5 is the Euclidean detection-volume weight|
|`effective\\\_solid\\\_angle\\\_scaling(...)`|closed form: 10.868 (sa sb)^(1/3) / (\|sin psi\| n^(2/3)) for power 1.5|
|`detection\\\_weight(sa, sb, psi)`|toy weight, proportional to (sa sb)^(-2/3) / \|sin psi\||
|`detector\\\_weighted\\\_solid\\\_angle(sa, sb, psi, "LIGO" or "LISA" or "LISA+confusion", n\\\_ref)`|solid angle weighted by a noise curve (cusp spectrum f^(-4/3), SNR^3 weighting)|
|`onaxis\\\_correction(n, side\\\_a, side\\\_b)`|predicted log(B\_n / leading law), relative order n^(-2/3), parameter-free; each side = (s, tau, k, w)|
|`lobe\\\_radius`, `is\\\_isolated`, `check\\\_validity`|helpers for deciding whether the model should be trusted|

```python
import numpy as np
import cuspbeam as cb
sa, sb, psi, n = 2.0, 1.5, 1.0, 5e7
cb.on\\\_axis\\\_B(sa, sb)                                   # on-axis amplitude factor
cb.effective\\\_solid\\\_angle(sa, sb, psi, n)               # steradians, detection-volume weight
cb.detector\\\_weighted\\\_solid\\\_angle(sa, sb, psi, "LIGO", n\\\_ref=n)
cb.check\\\_validity(n, sa, sb, psi)                      # \\\[] means "no warning", not "verified"
```

## Validation record (what was actually checked)

All against **exact numerical integrals** of the beam for random odd-harmonic loops (loop length 2 pi), for **isolated cusps of smooth loops**:

* **On-axis law.** The model reproduces K (sa sb)^(-2/3) exactly (tests).
* **Next-order correction.** On three fresh seed sets (64 to 76 cusps each, clean cusps), the median error at n = 2048 fell by factors of 128 to 165 relative to the leading law, observed/predicted = 1.00. The leftover error falls with a log-log slope between -1.37 and -1.7 across samples; the predicted n^(-4/3) is neither confirmed nor excluded.
* **Beam weight.** The exact weight inside the model's own lobe was compared with the model at n up to 8192 on 36 cusps in three runs (fresh seeds, pass criteria fixed in advance). The 24 cusps chosen with the neighbour rule all lie within \[0.88, 1.15] at n = 8192 (median 1.01, 22 within 5%). The one outlier of the first run had another exact cusp 0.099 rad away, which is why `is\\\_isolated` exists.
* **Scaling.** Fits over 96 combinations of (sa, sb, psi): weight \~ (sa sb)^(-0.659) (sin psi)^(-1.044) at n = 10^4, and exactly +1/3, -1 for the effective solid angle in the n -> infinity limit; the same exponents hold with LIGO and LISA noise curves.
* **Self-consistency tests:** `python -m pytest test\\\_cuspbeam.py -v` (12 tests, about 15 seconds).

## Limits (read before using)

1. **Isolated cusps only.** It knows nothing about other beams. Neighbouring exact cusps (within about 0.5 rad), and also near-coincidences of the two curves a couple of beam widths from the cusp (pseudocusp-like), add their own radiation. `is\\\_isolated` catches the first kind only; the second kind was seen once in testing (it made the exact weight 1.4 to 1.7 times the model at n = 512 to 2048) and is not detected by any helper.
2. **Smooth, centrally symmetric loops.** Random odd-harmonic loops with up to 4 harmonics were used for the weight check. On wigglier loops (8 to 16 harmonics) the single-cusp formulas removed little of the error, and a first multipoint extension helped only modestly.
3. **Tested range.** n up to 8192 against exact numerics; slopes 0.22 to 15.7; |sin psi| from 0.26 (on Turok loops from Suresh \& Chernoff: exact/model weight 0.99 and 0.98 at sin psi = 0.87 and 0.82; 0.70 at n = 1024 and 1.21 at n = 4096 for sin psi = 0.26; outside the model for sin psi = 0.048). The smallest-slope cusp (0.22) sat about 12% low and did not converge to the model. The detector bands correspond to much higher n, so use at those harmonics is an extrapolation of an asymptotic model.
4. **The correction is small at detector harmonics.** It is of order n^(-2/3), so tiny for n of order 10^6 and above. It is mainly useful as a benchmark for numerical codes.
5. **Toy weights.** `detection\\\_weight` and `detector\\\_weighted\\\_solid\\\_angle` assume Euclidean volume, a fixed threshold, one burst, no redshift, no antenna pattern, and a f^(-4/3) spectrum. Only ratios between cusps are meaningful; they are not event rates.
6. **Noise curves** are analytic approximations (aLIGO design fit; LISA analytic sensitivity of Robson, Cornish \& Liu 2019, with their Galactic confusion-noise fit for 4 years). The LIGO and LISA instrument curves were written from memory; check them against official curves before quoting numbers.
7. **No kinks, no backreaction, no cosmology.**

## Not yet done (would make it more useful)

A comparison against Suresh \& Chernoff's code (only their published formulas and test loops were used; see docs/comparison\_with\_suresh\_chernoff.md), a rate calculation with a published loop population, tests on loops from network simulations, and an expert's judgement on novelty.

## References

* Blanco-Pillado \& Olum, arXiv:1709.02693 (on-axis cusp power, App. A).
* Suresh \& Chernoff, arXiv:2310.00825 (beam models, multipoint method).
* Pazouli, Avgoustidis \& Copeland, arXiv:2008.13693 (cusp sharpness statistics).
* Robson, Cornish \& Liu 2019 (LISA sensitivity curve and confusion noise).

## How to cite / contact

\[pandesamir01@gmail.com]

