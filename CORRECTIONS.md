# Corrections

Published numbers that do not survive re-analysis, and what replaces them.

`tests/test_reproduces_published.py` is frozen: it verifies the **published path
still yields the published number**, so the manuscript stays reproducible.
`tests/test_current.py` pins the **corrected** values. Nothing here is deleted.

| # | published | corrected | where |
|---|---|---|---|
| 1 | Fig. 4C one-sided **P = 0.019** | **P ~ 0.21-0.27**, not significant | `figures/figure4.py` |
| 2 | Fig. 3B **rho = -0.43, n = 1545** | **rho = -0.690 (w5) / -0.813 (w10), n = 113 / 73** | `jacobian.dominant_eigenvalue` |
| 3 | "ninefold weakening of interactions" | amplitude, not interaction structure | `jacobian.corr` |
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

**Against DIVERSITY, nothing scale-free is robust.** Every statistic flips sign between
windows or fails significance in most configurations:

| vs diversity | c1 w5 | c1 w10 | all8 w5 | all8 w10 |
|---|---|---|---|---|
| frac_positive | +0.250 n.s. | +0.286 n.s. | +0.228 | +0.168 n.s. |
| mean_absolute | +0.316 | **-0.455** | +0.167 n.s. | **-0.348** |
| mean_negative | -0.241 n.s. | **+0.572** | -0.074 n.s. | **+0.502** |

**Published Figure 2C (mean negative J vs diversity, rho = +0.737) therefore has no
robust amplitude-free counterpart.** Its scale-free equivalent runs -0.241 at window 5
and +0.572 at window 10. Report the temporal result; treat the diversity relationship as
window-dependent, or state the amplitude caveat with it.

## 5. Control Paenibacillaceae was filtered, not absent

Derived tables elsewhere dropped control families with `mean relative abundance < 0.005`.
Control Paenibacillaceae means are 0.0050, 0.0028, 0.0010, 0.0046 - **the threshold sat at
the level of the signal**, so the taxon vanished from the control arm and downstream code
substituted zeros.

`data/16s/family/c_m*_family.csv` in this repository is unfiltered and shows it present
in 4/4 controls, max 0.0036-0.0137, never blooming, against 0.53-0.67 in colonised mice.
"Present at ~1% and never taking over" is the measurement; "absent" was an artefact.
