# Corrections

Published numbers that do not survive re-analysis, and what replaces them.

`tests/test_reproduces_published.py` is frozen: it verifies the **published path
still yields the published number**, so the manuscript stays reproducible.
`tests/test_current.py` pins the **corrected** values. Nothing here is deleted.

| # | published | corrected | where |
|---|---|---|---|
| 1 | Fig. 4C one-sided **P = 0.019** | not significant; **rebuilt** as a sign reversal | `figures/figure4c_rebuilt.py` |
| 2 | Fig. 3B **rho = -0.43, n = 1545** | **rho = -0.690 (w5) / -0.813 (w10), n = 113 / 73** | `jacobian.dominant_eigenvalue` |
| 3 | "ninefold weakening of interactions" | amplitude, not interaction structure | `jacobian.corr` |
| 3b | Fig. 2C at a single window | holds; amplitude-free counterpart needs win >= 9 | `jacobian.summarise_corr` |
| 4 | "inhibitory-to-facilitative shift" | shift toward **symmetry**, never past it | `summarise_corr.frac_positive` |

---

## 1. Figure 4C — P = 0.019 does not survive

The comparison window was defined relative to the control span while the axis was
converted with the wrong clock, which silently narrowed the test to 1-7 d. Both clocks
give P = 0.019 on that window; over the full series both give P ~ 0.27.

Controls are indexed by DAY (days 1-10), colonised by SAMPLING SLOT (1 = 3 h, 2 = 6 h,
3 = 12 h, 4 = day 1, ... 19 = day 16). `timeaxis.to_days` requires `group` and has no
default. See `tests/test_timeaxis.py`.

**This is a submitted statistic. Co-authors need to be told.**

### Rebuilt — `figures/figure4c_rebuilt.py`

Re-running the test was not enough; three faults had to be fixed together.

**Clock.** Controls are indexed by DAY, colonised by SAMPLING SLOT. Fixed.

**Dimension.** Control states carry 16S families only (4-7 variables); colonised states
are dominated by barcode lineages (12). Mean J over matrix entries depends on dimension.
Rebuilding from the unfiltered `data/16s/family/` tables, the shared subspace across all
twelve mice is exactly **two families** - Enterobacteriaceae and Paenibacillaceae.
Controls carry 17-28 families above 0.1%, colonised only 3-6: the colonizer collapses
resident diversity, so no wider matched subspace exists. The window is also fixed at
**5 days** in both arms rather than 5 slots, which are different real durations.

**Amplitude.** Raw |J| says colonised are 0.07x controls (P = 2.2e-12); amplitude-free
|R| says 1.40x (P = 0.041). The magnitude difference reverses under normalisation, so
magnitude is the amplitude term - exactly what the published panel was testing.

**What survives is the SIGN, and it reverses in both directions:**

| | controls | colonised | P |
|---|---|---|---|
| R[Paeni <- Entero] | **-0.142** | **+0.265** | 2.8e-3 |
| R[Entero <- Paeni] | **+0.303** | **-0.259** | 6.5e-6 |

Significant under the overlap window (days 1-10) and the full series alike. Under
antibiotic alone Enterobacteriaceae suppresses Paenibacillaceae; with the colonizer
present that suppression is gone.

Caveat that belongs on the figure: "Enterobacteriaceae" is a native population in the
controls and the gavaged K12 in the colonised mice, so this compares community contexts,
not the same organism.

## 2. Figure 3B — confounded and pseudoreplicated

`jacobian.eigenvalues` reads the stored matrices in `data/jacobian`, estimated with an
**expanding** window. A covariance estimator has rank <= samples in window, so
`rank = min(window, S)` exactly: rank is a deterministic function of time, and the
null-space zeros sit at Re(lambda) = 0, above the mostly-negative real spectrum,
padding early timepoints upward. It also returns one row per eigenvalue, so a mouse
contributes windows x species to the correlation rather than windows.

`jacobian.dominant_eigenvalue` uses sliding-window matrices from `jacobian.full` and
returns one row per window. Under it, **rank no longer tracks time**: rho = +0.092
(p = 0.33) at w5 and -0.110 (p = 0.36) at w10. Only 1.4-1.8% of windows pick a
null-space zero as the maximum.

Re(lambda_max) remains **> 0 at 99.1% (w5) and 100% (w10) of timepoints**. It declines
without crossing zero. The defensible statement is "the community moves toward, but does
not reach, a stable fixed point within 16 days" - not that stability increases.

Independent check: the same fix in a separate Python implementation gives -0.691 / -0.815
against -0.690 / -0.813 here.

## 3. The ninefold weakening is amplitude

`J[i<-j] = cov(dz_i/dt, z_j)` estimates A*C, not A. As succession proceeds the
trajectories flatten, C shrinks, and |J| shrinks with it even if A never changes.
`jacobian.corr` divides the amplitude out.

Against time, all 8 mice:

| window | mean \|J\| | mean \|R\| |
|---|---|---|
| 5 (n=113) | **-0.892** (p=3.9e-40) | +0.150 (p=0.11) |
| 10 (n=73) | **-0.842** (p=1.1e-20) | -0.086 (p=0.47) |

Per-mouse at w5 in an independent implementation: mean |J| significant in 8/8,
mean |R| in 2/8 with **opposite signs**. The decline does not survive amplitude removal.

`J = A*C` is still the Jacobian, and community dynamics depend on J rather than A, so
"effective coupling weakens" remains literally true. What is not supported is that
interaction *structure* changes.

## 4. What the amplitude-free data does support - and what it does not

**Against TIME, `frac_positive` rises**: rho = +0.292 (w5, p = 1.7e-3) and +0.426
(w10, p = 1.7e-4), consistent with the independent sign test (rho = +0.293). Inhibitory
links become less PREVALENT.

**It stops at symmetry.** Low- vs high-diversity means: 0.412 -> 0.470, 0.327 -> 0.471,
0.427 -> 0.511, 0.345 -> 0.487. Never durably above 0.5. Facilitation is not reached.

**Against DIVERSITY the relationship is real but needs a wide enough window.** The
scale-free statistic is a Pearson correlation and needs more samples to estimate than a
covariance does. Sweeping the window over all 8 mice:

| win | n | raw J | scale-free R | median \|R\| | IQR \|R\| |
|---|---|---|---|---|---|
| 3 | 129 | +0.632 | -0.083 n.s. | 0.525 | 0.230 |
| 5 | 113 | **+0.732** | -0.074 n.s. | 0.390 | 0.156 |
| 7 | 97 | +0.652 | +0.073 n.s. | 0.318 | 0.124 |
| 8 | 89 | +0.609 | +0.204 n.s. | 0.309 | 0.101 |
| 9 | 81 | +0.581 | **+0.383** | 0.280 | 0.079 |
| 10 | 73 | +0.574 | **+0.502** | 0.287 | 0.083 |
| 11 | 65 | +0.500 | **+0.412** | 0.287 | 0.086 |
| 12 | 57 | +0.355 | **+0.412** | 0.261 | 0.084 |

Median |R| falls from 0.525 to 0.261 and its IQR from 0.230 to 0.084 as the window
widens - the signature of small-sample correlation noise, which inflates |R| and
attenuates its correlation with anything else. The scale-free relationship appears from
window 9, exactly where R stabilises.

Not a range effect: restricting both windows to the same diversity range (1D > 1.5)
leaves win 5 at -0.018 n.s. and win 10 at +0.500. Same data span, different window.

**Reading.** Figure 2C is supported, with a stated caveat: the raw statistic is
significant at every window (+0.63 to +0.73), and its amplitude-free counterpart is
significant wherever the window is wide enough to estimate a correlation (win >= 9,
rho = +0.38 to +0.50). Report the window sweep rather than a single window. Note also
that diversity here is close to a relabelled time axis - it is flat at 1.0 until the
Paenibacillaceae onset then plateaus near 2.1 - so the diversity and time analyses are
largely the same result seen twice.

## 5. Control Paenibacillaceae was filtered, not absent

Derived tables elsewhere dropped control families with `mean relative abundance < 0.005`.
Control Paenibacillaceae means are 0.0050, 0.0028, 0.0010, 0.0046 - **the threshold sat at
the level of the signal**, so the taxon vanished from the control arm and downstream code
substituted zeros.

`data/16s/family/c_m*_family.csv` in this repository is unfiltered and shows it present
in 4/4 controls, max 0.0036-0.0137, never blooming, against 0.53-0.67 in colonised mice.
"Present at ~1% and never taking over" is the measurement; "absent" was an artefact.
