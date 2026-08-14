"""Supplementary Figure 5 — eigenvalue spectra in the complex plane, all 8 mice.

Figure 3A shows one mouse. This is the same construction for every mouse on
shared axes, so the leftward drift of the cloud can be compared across animals.
Stability is Re(lambda) < 0, i.e. the shaded half-plane.

Run:  python figures/supplementary/figS5_eigenvalue_spectra.py
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from _common import COLONISED, jacobian, style

style.apply()


def main():
    eig = pd.concat([jacobian.eigenvalues(m) for m in COLONISED], ignore_index=True)
    xr = np.percentile(eig["re"], [1, 99]) * 1.15
    yr = np.abs(np.percentile(eig["im"], [1, 99])).max() * 1.1

    fig, axes = plt.subplots(2, 4, figsize=(15, 6.8), sharex=True, sharey=True)
    for ax, m in zip(axes.flat, COLONISED):
        e = eig[eig["mouse"] == m]
        ax.axvspan(xr[0], 0, color="#e8f4e8", alpha=0.45, zorder=0)
        ax.axhline(0, color="#aaaaaa", lw=0.5, zorder=1)
        ax.axvline(0, color="#cc2222", lw=0.8, ls="--", zorder=3)
        sc = ax.scatter(e["re"], e["im"], c=e["relative_time"], cmap="plasma",
                        s=13, alpha=0.8, lw=0, zorder=2)
        ax.set_title(m, fontweight="bold")
    axes[0, 0].set_xlim(*xr)
    axes[0, 0].set_ylim(-yr, yr)
    for ax in axes[1, :]:
        ax.set_xlabel(r"Re($\lambda$)")
    for ax in axes[:, 0]:
        ax.set_ylabel(r"Im($\lambda$)")
    cb = fig.colorbar(sc, ax=axes.ravel().tolist(), shrink=0.6, pad=0.015)
    cb.set_label("Time"); cb.set_ticks([0, 0.5, 1])
    cb.set_ticklabels(["Early", "Mid", "Late"])
    fig.suptitle("Supplementary Figure 5 — Jacobian eigenvalue spectra, all mice", y=0.98)
    return style.save(fig, "figS5_eigenvalue_spectra")


if __name__ == "__main__":
    print("Supplementary Figure 5")
    main()
