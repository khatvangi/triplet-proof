#!/usr/bin/env python3
"""
Fig 1 — SGC optimization is a factorial product of architecture and alphabet.

Three panels for a JME Letter (double-column, 7.2 × 2.8 in):
  a) Null distribution for 20-AA triplet SGC (Condition B) with SGC baseline
  b) Three stacked mini-histograms comparing Conditions A, C, B
  c) Position-specific synonymous fraction bar chart

Reads:
  results/null_A_doublet_10aa.npy   (1M × 2, columns [D, E])
  results/null_B_triplet_20aa.npy
  results/null_C_triplet_10aa.npy
  results/sgc_baselines.json

Outputs:
  figures/fig1_optimization_factorial.pdf
  figures/fig1_optimization_factorial.png (300 dpi)
"""

import json
import pathlib

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
from matplotlib.lines import Line2D

# ---------------------------------------------------------------------------
# paths
# ---------------------------------------------------------------------------
ROOT = pathlib.Path(__file__).resolve().parent.parent
RES = ROOT / "results"
OUT = ROOT / "figures"

# ---------------------------------------------------------------------------
# global style
# ---------------------------------------------------------------------------
mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Liberation Sans", "DejaVu Sans"],
    "font.size": 8,
    "axes.labelsize": 8,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "axes.linewidth": 0.6,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "xtick.major.size": 3,
    "ytick.major.size": 3,
    "pdf.fonttype": 42,       # editable text in PDF
    "ps.fonttype": 42,
})

# ---------------------------------------------------------------------------
# Paul Tol colorblind-safe palette
# ---------------------------------------------------------------------------
BLUE = "#4477AA"
RED = "#CC6677"
GREEN = "#228833"
GREY = "#DDDDDD"

# ---------------------------------------------------------------------------
# load data
# ---------------------------------------------------------------------------
null_A = np.load(RES / "null_A_doublet_10aa.npy")[:, 0]   # D column
null_B = np.load(RES / "null_B_triplet_20aa.npy")[:, 0]
null_C = np.load(RES / "null_C_triplet_10aa.npy")[:, 0]

with open(RES / "sgc_baselines.json") as fh:
    sgc = json.load(fh)

sgc_A = sgc["A"]["D"]
sgc_B = sgc["B"]["D"]
sgc_C = sgc["C"]["D"]

# ---------------------------------------------------------------------------
# helper: clean axes (left + bottom spines only, no grid)
# ---------------------------------------------------------------------------
def clean_ax(ax, left=True, bottom=True):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(left)
    ax.spines["bottom"].set_visible(bottom)
    ax.tick_params(axis="both", which="both", direction="out")
    ax.grid(False)

# ---------------------------------------------------------------------------
# figure layout — use manual subplot positioning for precise control
# ---------------------------------------------------------------------------
fig = plt.figure(figsize=(7.2, 2.8))

# column boundaries (in figure fraction)
# panel A: x = [0.07, 0.30]
# panel B: x = [0.38, 0.64]  (wider to accommodate brackets)
# panel C: x = [0.74, 0.96]

# ---------------------------------------------------------------------------
# Panel A: SGC optimization histogram (Condition B)
# ---------------------------------------------------------------------------
ax_a = fig.add_axes([0.07, 0.19, 0.23, 0.72])

ax_a.hist(null_B, bins=100, color=GREY, edgecolor="black", linewidth=0.3)
ax_a.axvline(sgc_B, color=RED, linewidth=1.2, zorder=5)

# text annotation — placed to the right of the SGC line so it doesn't clip
ymax_a = ax_a.get_ylim()[1]
ax_a.text(
    sgc_B + 0.3, ymax_a * 0.82,
    "SGC\n$z = -16.6$",
    color=RED, fontsize=7, fontweight="bold",
    ha="left", va="top",
)

# add a small arrow from text to line
ax_a.annotate(
    "", xy=(sgc_B, ymax_a * 0.72), xytext=(sgc_B + 0.25, ymax_a * 0.72),
    arrowprops=dict(arrowstyle="-", color=RED, lw=0.6),
)

ax_a.set_xlabel("Noise distortion $D$")
ax_a.set_ylabel("Count")
clean_ax(ax_a)
ax_a.text(
    -0.18, 1.05, "a", transform=ax_a.transAxes,
    fontsize=10, fontweight="bold", va="bottom",
)

# ---------------------------------------------------------------------------
# Panel B: Factorial comparison — three stacked mini-histograms
# each gets its own x-range since the distributions differ widely
# ---------------------------------------------------------------------------
panel_b_left = 0.40
panel_b_width = 0.22
panel_b_bottom = 0.19
panel_b_total_h = 0.72
mini_h = panel_b_total_h * 0.26       # height of each mini-hist
mini_gap = panel_b_total_h * 0.11     # gap between mini-hists

conditions = [
    # (data, sgc_val, color, label)
    (null_A, sgc_A, BLUE,  "10 AA, $n$=2,  $z=-3.8$"),
    (null_C, sgc_C, GREEN, "10 AA, $n$=3,  $z=-14.8$"),
    (null_B, sgc_B, RED,   "20 AA, $n$=3,  $z=-16.6$"),
]

ax_b = []
for i, (data, sgc_val, color, label) in enumerate(conditions):
    # stack from top to bottom: i=0 is top
    y0 = panel_b_bottom + panel_b_total_h - (i + 1) * mini_h - i * mini_gap
    ax = fig.add_axes([panel_b_left, y0, panel_b_width, mini_h])
    ax_b.append(ax)

    ax.hist(data, bins=80, color=color, alpha=0.7, edgecolor="none")
    ax.axvline(sgc_val, color=color, linewidth=1.0, linestyle="-")

    # set per-panel x-limits so each distribution fills its panel
    lo = min(sgc_val, data.min()) - 0.5
    hi = data.max() + 0.5
    ax.set_xlim(lo, hi)

    # label inside the panel, top-right
    ax.text(
        0.97, 0.85, label,
        transform=ax.transAxes,
        fontsize=5.5, ha="right", va="top", color="black",
    )

    clean_ax(ax)
    ax.tick_params(axis="x", labelsize=5.5)
    ax.tick_params(axis="y", labelsize=5.5)
    # reduce number of y ticks
    ax.locator_params(axis="y", nbins=3)

    # all mini-hists get x-tick labels since ranges differ
    if i < 2:
        ax.tick_params(axis="x", labelbottom=True)

ax_b[-1].set_xlabel("Noise distortion $D$", fontsize=7)
# add shared x-label text for upper panels too (just the bottom one has xlabel)

# panel label "b" above the top mini-hist
ax_b[0].text(
    -0.22, 1.12, "b", transform=ax_b[0].transAxes,
    fontsize=10, fontweight="bold", va="bottom",
)

# --- bracket annotations between mini-histograms ---
# force a draw so positions are finalized
fig.canvas.draw()

bracket_x = panel_b_left + panel_b_width + 0.015   # right of panel B

# architecture bracket: between A (top) and C (middle)
pos0 = ax_b[0].get_position()
pos1 = ax_b[1].get_position()
brace_top_1 = pos0.y0
brace_bot_1 = pos1.y1
brace_mid_1 = (brace_top_1 + brace_bot_1) / 2
bx = bracket_x

# vertical line
fig.add_artist(Line2D(
    [bx, bx], [brace_bot_1, brace_top_1],
    transform=fig.transFigure, color="black", linewidth=0.6, clip_on=False,
))
# top tick
fig.add_artist(Line2D(
    [bx, bx + 0.008], [brace_top_1, brace_top_1],
    transform=fig.transFigure, color="black", linewidth=0.6, clip_on=False,
))
# bottom tick
fig.add_artist(Line2D(
    [bx, bx + 0.008], [brace_bot_1, brace_bot_1],
    transform=fig.transFigure, color="black", linewidth=0.6, clip_on=False,
))
# text
fig.text(
    bx + 0.012, brace_mid_1,
    "architecture\n$\\Delta z = -11.0$",
    fontsize=5.5, ha="left", va="center", color="black",
    fontstyle="italic",
)

# alphabet bracket: between C (middle) and B (bottom)
pos2 = ax_b[2].get_position()
brace_top_2 = pos1.y0
brace_bot_2 = pos2.y1
brace_mid_2 = (brace_top_2 + brace_bot_2) / 2

fig.add_artist(Line2D(
    [bx, bx], [brace_bot_2, brace_top_2],
    transform=fig.transFigure, color="black", linewidth=0.6, clip_on=False,
))
fig.add_artist(Line2D(
    [bx, bx + 0.008], [brace_top_2, brace_top_2],
    transform=fig.transFigure, color="black", linewidth=0.6, clip_on=False,
))
fig.add_artist(Line2D(
    [bx, bx + 0.008], [brace_bot_2, brace_bot_2],
    transform=fig.transFigure, color="black", linewidth=0.6, clip_on=False,
))
fig.text(
    bx + 0.012, brace_mid_2,
    "alphabet\n$\\Delta z = -1.8$",
    fontsize=5.5, ha="left", va="center", color="black",
    fontstyle="italic",
)

# ---------------------------------------------------------------------------
# Panel C: Position-specific synonymy bar chart
# ---------------------------------------------------------------------------
ax_c = fig.add_axes([0.76, 0.19, 0.21, 0.72])

# verified synonymous fraction (%) values
syn = {
    "A": [33.3, 8.3, None],     # doublet: no Position 3
    "B": [4.4, 0.0, 68.9],      # 20-AA SGC
    "C": [26.2, 10.9, 76.5],    # 10-AA triplet
}

n_groups = 3
bar_w = 0.22
x = np.arange(n_groups)

# condition A (blue), B (red), C (green) — order within each group
bars_A = [syn["A"][0], syn["A"][1], 0]     # 0 placeholder for absent bar
bars_B = [syn["B"][0], syn["B"][1], syn["B"][2]]
bars_C = [syn["C"][0], syn["C"][1], syn["C"][2]]

rects_A = ax_c.bar(x - bar_w, bars_A, bar_w, color=BLUE, edgecolor="none",
                    label="10 AA, $n$=2")
rects_B = ax_c.bar(x, bars_B, bar_w, color=RED, edgecolor="none",
                    label="20 AA, $n$=3")
rects_C = ax_c.bar(x + bar_w, bars_C, bar_w, color=GREEN, edgecolor="none",
                    label="10 AA, $n$=3")

# mark Condition A Position 3 as N/A
ax_c.text(
    2 - bar_w, 1.5, "N/A", ha="center", va="bottom",
    fontsize=5, color=BLUE, fontstyle="italic",
)

ax_c.set_xticks(x)
ax_c.set_xticklabels(["Pos 1", "Pos 2", "Pos 3"], fontsize=6.5)
ax_c.set_ylabel("Synonymous fraction (%)")
ax_c.set_ylim(0, 95)

ax_c.legend(fontsize=5, loc="upper left", frameon=False, handlelength=1.0,
            handletextpad=0.4, labelspacing=0.3)
clean_ax(ax_c)

ax_c.text(
    -0.18, 1.05, "c", transform=ax_c.transAxes,
    fontsize=10, fontweight="bold", va="bottom",
)

# ---------------------------------------------------------------------------
# save
# ---------------------------------------------------------------------------
OUT.mkdir(exist_ok=True)

fig.savefig(OUT / "fig1_optimization_factorial.pdf", dpi=300)
fig.savefig(OUT / "fig1_optimization_factorial.png", dpi=300)
print(f"saved: {OUT / 'fig1_optimization_factorial.pdf'}")
print(f"saved: {OUT / 'fig1_optimization_factorial.png'}")
plt.close(fig)
