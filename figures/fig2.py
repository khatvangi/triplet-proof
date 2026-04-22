#!/usr/bin/env python3
"""
Figure 2 — SGC is beyond the null.

three panels (double-column, 7.2 x 3.0 in):
  a) null histogram of Dirichlet energy E across 1M random codes, SGC at z = -14.6
  b) null histogram of noise distortion D, SGC at z = -16.6
  c) joint scatter of (D, E), 10k-point subsample, SGC as red star

reads:
  results/null_B_triplet_20aa.npy  (1M x 2, columns [D, E])
  results/sgc_baselines.json
  results/publication_controls.json

outputs:
  figures/fig2.pdf, fig2.png (300 dpi)
"""

import json
import pathlib
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
import seaborn as sns

ROOT = pathlib.Path(__file__).resolve().parent.parent
RES = ROOT / "results"
OUT = ROOT / "figures"

# ── style ──
sns.set_style("ticks")
mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Liberation Sans", "DejaVu Sans"],
    "font.size": 9,
    "axes.labelsize": 10,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "legend.fontsize": 8,
    "axes.linewidth": 0.6,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})

RED = "#CC6677"
GREY = "#BBBBBB"
GREY_LIGHT = "#DDDDDD"

# ── load data ──
null_B = np.load(RES / "null_B_triplet_20aa.npy")
D_null, E_null = null_B[:, 0], null_B[:, 1]

with open(RES / "sgc_baselines.json") as f:
    sgc = json.load(f)
D_sgc, E_sgc = sgc["B"]["D"], sgc["B"]["E"]

with open(RES / "publication_controls.json") as f:
    pub = json.load(f)
z_D = pub["B_20AA_n3"]["D"]["z_score"]
z_E = pub["B_20AA_n3"]["E"]["z_score"]
corr = pub["B_20AA_n3"]["DE_correlation"]

# ── figure: 3 panels in a row ──
fig = plt.figure(figsize=(7.5, 3.0))
gs = GridSpec(1, 3, figure=fig, width_ratios=[1, 1, 1.2],
              left=0.085, right=0.97, bottom=0.17, top=0.92, wspace=0.45)


def panel_hist(ax, data, sgc_val, z_score, xlabel, letter, p_note="$p < 10^{-6}$"):
    """draw a null histogram with SGC line + label."""
    ax.hist(data, bins=100, color=GREY, edgecolor="black", linewidth=0.2)
    ax.axvline(sgc_val, color=RED, linewidth=1.2, zorder=5)

    ymax = ax.get_ylim()[1]
    # SGC label + arrow
    ax.annotate(
        f"SGC\n$z = {z_score:.1f}$\n{p_note}",
        xy=(sgc_val, ymax * 0.72),
        xytext=(sgc_val + (data.max() - data.min()) * 0.05, ymax * 0.72),
        color=RED, fontsize=8, fontweight="bold",
        ha="left", va="center",
        arrowprops=dict(arrowstyle="->", color=RED, lw=0.7),
    )

    ax.set_xlabel(xlabel)
    ax.set_ylabel("count (of 10$^6$ random codes)")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.text(-0.20, 1.02, letter, transform=ax.transAxes,
            fontsize=11, fontweight="bold", va="bottom", ha="left")


# ── panel a: Dirichlet energy histogram ──
ax_a = fig.add_subplot(gs[0, 0])
panel_hist(ax_a, E_null, E_sgc, z_E, "Dirichlet energy $E$", "a")

# ── panel b: noise distortion histogram ──
ax_b = fig.add_subplot(gs[0, 1])
panel_hist(ax_b, D_null, D_sgc, z_D, "noise distortion $D$", "b")

# ── panel c: joint scatter with marginals ──
inner = GridSpecFromSubplotSpec(
    2, 2, subplot_spec=gs[0, 2],
    width_ratios=[4, 1], height_ratios=[1, 4],
    hspace=0.05, wspace=0.05,
)
ax_scat = fig.add_subplot(inner[1, 0])
ax_hx = fig.add_subplot(inner[0, 0], sharex=ax_scat)
ax_hy = fig.add_subplot(inner[1, 1], sharey=ax_scat)

# subsample 10k points from 1M
rng = np.random.default_rng(42)
idx = rng.choice(len(D_null), size=10000, replace=False)
ax_scat.scatter(
    D_null[idx], E_null[idx], s=4, c=GREY, alpha=0.3,
    edgecolors="none", rasterized=True,
)

# SGC star — make distinctly larger and offset so it's visually separate
ax_scat.plot(
    D_sgc, E_sgc, marker="*", markersize=14, color=RED,
    markeredgecolor="black", markeredgewidth=0.6, zorder=10,
)
ax_scat.annotate(
    "SGC",
    xy=(D_sgc, E_sgc),
    xytext=(D_sgc + 1.8, E_sgc + 250),
    fontsize=9, color=RED, fontweight="bold",
    arrowprops=dict(arrowstyle="->", color=RED, lw=0.7),
    ha="left", va="bottom",
)

# correlation
ax_scat.text(
    0.97, 0.04, f"$r = {corr:.3f}$",
    transform=ax_scat.transAxes,
    fontsize=8, ha="right", va="bottom", color="0.3",
)

ax_scat.set_xlabel("noise distortion $D$")
ax_scat.set_ylabel("Dirichlet energy $E$")
ax_scat.spines["top"].set_visible(False)
ax_scat.spines["right"].set_visible(False)

# marginals
ax_hx.hist(D_null, bins=120, color=GREY, edgecolor="none", alpha=0.8)
ax_hx.axvline(D_sgc, color=RED, linewidth=1.0)
ax_hx.tick_params(labelbottom=False, labelleft=False, left=False, bottom=False)
for s in ["top", "right", "left", "bottom"]:
    ax_hx.spines[s].set_visible(False)

ax_hy.hist(E_null, bins=120, orientation="horizontal",
           color=GREY, edgecolor="none", alpha=0.8)
ax_hy.axhline(E_sgc, color=RED, linewidth=1.0)
ax_hy.tick_params(labelleft=False, labelbottom=False, left=False, bottom=False)
for s in ["top", "right", "left", "bottom"]:
    ax_hy.spines[s].set_visible(False)

# panel c label on the top-marginal axis
ax_hx.text(-0.15, 1.02, "c", transform=ax_hx.transAxes,
           fontsize=11, fontweight="bold", va="bottom", ha="left")

# ── save ──
OUT.mkdir(exist_ok=True)
fig.savefig(OUT / "fig2.pdf", dpi=300)
fig.savefig(OUT / "fig2.png", dpi=300)
print(f"saved: {OUT / 'fig2.pdf'}")
print(f"saved: {OUT / 'fig2.png'}")
plt.close(fig)
