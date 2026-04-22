#!/usr/bin/env python3
"""
Figure 3 — triplet architecture enables deep optimization.

two panels (double-column, 7.2 x 3.0 in):
  a) three null distributions of D (noise distortion) on a shared x-axis:
     blue = condition A (doublet, 10 AA), green = condition C (triplet, 10 AA),
     red = condition B (SGC, 20 AA). each with its observed SGC line.
  b) position-specific synonymy fraction, grouped bar chart.

NOTE: no Δz "architecture" / "alphabet" brackets — manuscript moved away from
subtracting incommensurable z-scores. raw z-scores only.

reads:
  results/null_{A,B,C}_*.npy (1M x 2, columns [D, E])
  results/sgc_baselines.json
  results/publication_controls.json

outputs:
  figures/fig3.pdf, fig3.png (300 dpi)
"""

import json
import pathlib
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
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

BLUE = "#4477AA"    # condition A (doublet)
GREEN = "#228833"   # condition C (10-AA triplet)
RED = "#CC6677"     # condition B (SGC, 20-AA triplet)

# ── load data ──
null_A = np.load(RES / "null_A_doublet_10aa.npy")[:, 0]  # D values
null_B = np.load(RES / "null_B_triplet_20aa.npy")[:, 0]
null_C = np.load(RES / "null_C_triplet_10aa.npy")[:, 0]

with open(RES / "sgc_baselines.json") as f:
    sgc = json.load(f)
sgc_A, sgc_B, sgc_C = sgc["A"]["D"], sgc["B"]["D"], sgc["C"]["D"]

with open(RES / "publication_controls.json") as f:
    pub = json.load(f)
zA = pub["A_10AA_n2"]["D"]["z_score"]
zB = pub["B_20AA_n3"]["D"]["z_score"]
zC = pub["C_10AA_n3"]["D"]["z_score"]

# ── figure: 2 panels ──
fig = plt.figure(figsize=(7.2, 3.0))
gs = GridSpec(1, 2, figure=fig, width_ratios=[1.4, 1.0],
              left=0.08, right=0.97, bottom=0.17, top=0.92, wspace=0.30)

ax_a = fig.add_subplot(gs[0, 0])
ax_b = fig.add_subplot(gs[0, 1])

# ═════════════════════════════════════════════════════════════
# panel a — overlaid null distributions (shared x-axis)
# ═════════════════════════════════════════════════════════════
bins = np.linspace(17.5, 35.5, 120)

for data, color, label, sgc_val, z in [
    (null_A, BLUE,  f"A: 10 AA, $n$=2  ($z = {zA:.1f}$)", sgc_A, zA),
    (null_C, GREEN, f"C: 10 AA, $n$=3  ($z = {zC:.1f}$)", sgc_C, zC),
    (null_B, RED,   f"B: 20 AA, $n$=3 (SGC)  ($z = {zB:.1f}$)", sgc_B, zB),
]:
    ax_a.hist(data, bins=bins, color=color, alpha=0.5, edgecolor="none",
              label=label, zorder=2)
    ax_a.axvline(sgc_val, color=color, linewidth=1.4, zorder=5)

# annotate threshold crossing
ax_a.text(
    0.98, 0.98,
    "A lies inside its null tail\nB and C lie beyond",
    transform=ax_a.transAxes, fontsize=8, ha="right", va="top",
    fontstyle="italic", color="0.3",
)

ax_a.set_xlabel("noise distortion $D$")
ax_a.set_ylabel("count (of 10$^6$ random codes)")
ax_a.set_xlim(17.5, 36)
ax_a.spines["top"].set_visible(False)
ax_a.spines["right"].set_visible(False)

# legend with null hist + SGC line combined per condition
legend_handles = [
    (Patch(facecolor=BLUE, alpha=0.5),
     Line2D([0], [0], color=BLUE, lw=1.4)),
    (Patch(facecolor=GREEN, alpha=0.5),
     Line2D([0], [0], color=GREEN, lw=1.4)),
    (Patch(facecolor=RED, alpha=0.5),
     Line2D([0], [0], color=RED, lw=1.4)),
]
legend_labels = [
    f"A: 10 AA, $n$=2   ($z_D = {zA:.1f}$)",
    f"C: 10 AA, $n$=3   ($z_D = {zC:.1f}$)",
    f"B: 20 AA, $n$=3 (SGC)  ($z_D = {zB:.1f}$)",
]
ax_a.legend(
    legend_handles, legend_labels,
    loc="upper left", frameon=False, fontsize=7.5,
    handler_map={tuple: mpl.legend_handler.HandlerTuple(ndivide=None)},
    handlelength=2.0, handletextpad=0.5,
)

ax_a.text(-0.12, 1.02, "a", transform=ax_a.transAxes,
          fontsize=11, fontweight="bold", va="bottom", ha="left")

# ═════════════════════════════════════════════════════════════
# panel b — position-specific synonymy (grouped bars)
# ═════════════════════════════════════════════════════════════
syn = {
    "A": [33.3, 8.3, None],     # doublet: no pos 3
    "B": [4.4, 0.0, 68.9],      # 20-AA SGC
    "C": [26.2, 10.9, 76.5],    # 10-AA triplet
}

n_groups = 3
bar_w = 0.25
x = np.arange(n_groups)

bars_A = [syn["A"][0], syn["A"][1], 0]      # 0 placeholder, N/A
bars_B = [syn["B"][0], syn["B"][1], syn["B"][2]]
bars_C = [syn["C"][0], syn["C"][1], syn["C"][2]]

ax_b.bar(x - bar_w, bars_A, bar_w, color=BLUE, edgecolor="none",
         label="A: 10 AA, $n$=2")
ax_b.bar(x,         bars_B, bar_w, color=RED, edgecolor="none",
         label="B: 20 AA, $n$=3 (SGC)")
ax_b.bar(x + bar_w, bars_C, bar_w, color=GREEN, edgecolor="none",
         label="C: 10 AA, $n$=3")

# N/A marker for condition A at pos 3 — dark text, readable
ax_b.text(
    2 - bar_w, 2.0, "N/A", ha="center", va="bottom",
    fontsize=7, color="black",
)

ax_b.set_xticks(x)
ax_b.set_xticklabels(["pos 1", "pos 2", "pos 3"])
ax_b.set_ylabel("synonymous fraction (%)")
ax_b.set_ylim(0, 95)
ax_b.spines["top"].set_visible(False)
ax_b.spines["right"].set_visible(False)

ax_b.legend(loc="upper left", frameon=False, fontsize=7.5,
            handlelength=1.5, handletextpad=0.5)

ax_b.text(-0.18, 1.02, "b", transform=ax_b.transAxes,
          fontsize=11, fontweight="bold", va="bottom", ha="left")

# ── save ──
OUT.mkdir(exist_ok=True)
fig.savefig(OUT / "fig3.pdf", dpi=300)
fig.savefig(OUT / "fig3.png", dpi=300)
print(f"saved: {OUT / 'fig3.pdf'}")
print(f"saved: {OUT / 'fig3.png'}")
plt.close(fig)
