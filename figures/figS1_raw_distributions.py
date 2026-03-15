"""
Figure S1: null-distribution histograms for three conditions (D metric).

Panels:
  a) 10 amino acids, n=2 (doublet) — SGC inside null tail
  b) 10 amino acids, n=3 (triplet) — SGC far below null
  c) 20 amino acids, n=3 (triplet) — SGC far below null

Reads:
  results/null_A_doublet_10aa.npy   (col 0 = D, col 1 = E)
  results/null_B_triplet_20aa.npy
  results/null_C_triplet_10aa.npy
  results/sgc_baselines.json

Saves:
  figures/figS1_raw_distributions.pdf
  figures/figS1_raw_distributions.png  (300 dpi)
"""

import numpy as np
import json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

# ── paths ──────────────────────────────────────────────────────────────
root = Path(__file__).resolve().parent.parent
res  = root / "results"
out  = root / "figures"

# ── load data ──────────────────────────────────────────────────────────
null_A = np.load(res / "null_A_doublet_10aa.npy")[:, 0]   # D values
null_B = np.load(res / "null_B_triplet_20aa.npy")[:, 0]
null_C = np.load(res / "null_C_triplet_10aa.npy")[:, 0]

with open(res / "sgc_baselines.json") as f:
    bl = json.load(f)

sgc_A = bl["A"]["D"]   # 18.94
sgc_B = bl["B"]["D"]   # 23.92
sgc_C = bl["C"]["D"]   # 21.27

# ── Paul Tol palette ──────────────────────────────────────────────────
BLUE  = "#4477AA"
RED   = "#CC6677"
GREEN = "#228833"

# ── global font settings ──────────────────────────────────────────────
plt.rcParams.update({
    "font.family":    "sans-serif",
    "font.size":      8,
    "axes.labelsize": 8,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
})

# ── figure ─────────────────────────────────────────────────────────────
fig, axes = plt.subplots(3, 1, figsize=(3.5, 5.5))

# ── panel configuration ───────────────────────────────────────────────
# (data, sgc_val, color, title, label_letter, stats_text)
panels = [
    (null_A, sgc_A, BLUE,
     "10 amino acids, $n$ = 2",
     "a",
     "z = \u22123.8\np = 8.3 \u00d7 10\u207b\u2075\nCV = 3.4%"),
    (null_C, sgc_C, GREEN,
     "10 amino acids, $n$ = 3",
     "b",
     "z = \u221214.8\np < 10\u207b\u2076\nCV = 2.1%"),
    (null_B, sgc_B, RED,
     "20 amino acids, $n$ = 3",
     "c",
     "z = \u221216.6\np < 10\u207b\u2076\nCV = 1.6%"),
]

for ax, (data, sgc, color, title, letter, stats) in zip(axes, panels):

    # -- determine x range --
    # show full null distribution + SGC line with breathing room
    null_lo, null_hi = data.min(), data.max()
    span = null_hi - null_lo

    if sgc < null_lo:
        # sgc is below the null — pad left of sgc and right of null
        x_lo = sgc - 0.15 * span
        x_hi = null_hi + 0.05 * span
    else:
        # sgc falls within the null (panel A)
        x_lo = min(sgc, null_lo) - 0.05 * span
        x_hi = null_hi + 0.05 * span

    # -- histogram --
    ax.hist(data, bins=80, density=True, range=(null_lo, null_hi),
            color=color, alpha=0.6, edgecolor="black", linewidth=0.3)

    # -- SGC baseline (solid red) --
    ax.axvline(sgc, color="#CC6677", linewidth=1.2, zorder=5)

    # -- null minimum (dashed grey) --
    ax.axvline(null_lo, color="grey", linewidth=0.8, linestyle="--", zorder=4)

    # -- x limits --
    ax.set_xlim(x_lo, x_hi)

    # -- stats text box (top-right) --
    ax.text(0.97, 0.95, stats, transform=ax.transAxes,
            fontsize=6.5, verticalalignment="top", horizontalalignment="right",
            bbox=dict(boxstyle="round,pad=0.3", facecolor="white",
                      edgecolor="0.7", alpha=0.9))

    # -- panel label --
    ax.text(-0.12, 1.05, letter, transform=ax.transAxes,
            fontsize=10, fontweight="bold", va="top")

    # -- title --
    ax.set_title(title, fontsize=8, loc="center", pad=4)

    # -- axis labels --
    ax.set_xlabel("Noise distortion $D$", fontsize=8)

    # -- spines --
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # -- no gridlines --
    ax.grid(False)

# shared y-axis label on the left side of the figure
fig.text(0.01, 0.5, "Density", va="center", rotation="vertical", fontsize=8)

plt.tight_layout(rect=[0.05, 0, 1, 1], h_pad=1.5)

# ── save ───────────────────────────────────────────────────────────────
fig.savefig(out / "figS1_raw_distributions.pdf", bbox_inches="tight")
fig.savefig(out / "figS1_raw_distributions.png", dpi=300, bbox_inches="tight")
print("saved  figures/figS1_raw_distributions.pdf")
print("saved  figures/figS1_raw_distributions.png")
