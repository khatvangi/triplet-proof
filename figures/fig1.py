#!/usr/bin/env python3
"""
Figure 1 — framework.

three panels (double-column, 7.2 x 2.8 in):
  a) codon Hamming graph (n=3), 61 sense codons, nodes colored by AA class
  b) 20 AAs in PC1-PC2 space, colored by same class scheme
  c) schematic subgraph centered on GCU (Ala) with real ||Δf|| edge distances

reads:
  codon_table.csv, src/io/aa_props_lib.py, data/processed/aa_props.parquet

outputs:
  figures/fig1.pdf, fig1.png (300 dpi)
"""

import os, sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import networkx as nx

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
OUT = ROOT / "figures"

from src.io.codon_io import load_codon_table_csv
from src.io.aa_props_io import load_aa_props, aa_prop_vector_map
from src.sims.codon_graph import build_graph_n

# ── style: match CLAUDE.md typography spec ──
mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Liberation Sans", "DejaVu Sans"],
    "font.size": 9,
    "axes.labelsize": 10,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "legend.fontsize": 8,
    "axes.linewidth": 0.6,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})

# ── RAA10 class mapping (Murphy 2000-style) ──
RAA10 = {
    "A": "Small nonpolar", "G": "Small nonpolar",
    "V": "Aliphatic", "L": "Aliphatic", "I": "Aliphatic", "M": "Aliphatic",
    "F": "Aromatic", "W": "Aromatic", "Y": "Aromatic",
    "S": "Hydroxylic", "T": "Hydroxylic",
    "N": "Amidic", "Q": "Amidic",
    "D": "Acidic", "E": "Acidic",
    "K": "Basic", "R": "Basic",
    "H": "Imidazole",
    "P": "Cyclic",
    "C": "Sulfur",
}
GROUP_ORDER = ["Small nonpolar", "Aliphatic", "Aromatic", "Hydroxylic",
               "Amidic", "Acidic", "Basic", "Imidazole", "Cyclic", "Sulfur"]
GROUP_COLORS = {
    "Small nonpolar": "#332288", "Aliphatic": "#88CCEE",
    "Aromatic": "#44AA99", "Hydroxylic": "#117733",
    "Amidic": "#999933", "Acidic": "#CC6677",
    "Basic": "#882255", "Imidazole": "#AA4499",
    "Cyclic": "#DDCC77", "Sulfur": "#EE8866",
}

# ── load data ──
codon_df = load_codon_table_csv()
aa_df, _ = load_aa_props("data/processed/aa_props.parquet")
aa_vec = aa_prop_vector_map(aa_df)
c2a = dict(zip(codon_df["codon"].str.upper(), codon_df["aa"].str.upper()))

# per-component variance for PC labels — recompute quickly
from src.io.aa_props_lib import AA
from src.io.aa_props_io import _numeric_cols, _drop_highly_correlated, _zscore
rows = [{"aa": aa, **props} for aa, props in AA.items()]
lib_df = pd.DataFrame(rows)
raw_cols = _numeric_cols(lib_df)
used_raw = _drop_highly_correlated(lib_df, raw_cols, 0.95)
sub = lib_df[["aa"] + used_raw].drop_duplicates("aa").set_index("aa").sort_index()
X = sub.to_numpy(dtype=float)
Xz, _, _ = _zscore(X)
_, S, _ = np.linalg.svd(Xz, full_matrices=False)
var_frac = (S**2) / (S**2).sum()

# ── figure canvas ──
fig = plt.figure(figsize=(7.2, 2.8))
gs = GridSpec(1, 3, figure=fig, width_ratios=[1.05, 1, 1.05],
              left=0.02, right=0.98, bottom=0.27, top=0.95, wspace=0.12)

ax_a = fig.add_subplot(gs[0, 0])
ax_b = fig.add_subplot(gs[0, 1])
ax_c = fig.add_subplot(gs[0, 2])

# ═════════════════════════════════════════════════════════════
# panel a — Hamming graph of 61 sense codons
# ═════════════════════════════════════════════════════════════
G = build_graph_n(3)
stops = set(codon_df.loc[codon_df["is_stop"], "codon"].str.upper())
G_sense = G.copy()
G_sense.remove_nodes_from(stops)

# spring layout with strong repulsion
pos = nx.spring_layout(G_sense, seed=7, k=0.35, iterations=300)

# draw edges (light grey, behind nodes)
nx.draw_networkx_edges(
    G_sense, pos, ax=ax_a, edge_color="0.85", width=0.3, alpha=0.7,
)

# draw nodes colored by AA class
node_colors = []
for codon in G_sense.nodes():
    aa = c2a.get(codon, "*")
    grp = RAA10.get(aa, None)
    if grp is None:
        node_colors.append("#CCCCCC")
    else:
        node_colors.append(GROUP_COLORS[grp])

nx.draw_networkx_nodes(
    G_sense, pos, ax=ax_a, node_size=35, node_color=node_colors,
    edgecolors="black", linewidths=0.3,
)

# anchor labels — peripheral codons only (AUG lands in dense center, skip it)
# positions chosen from spring layout at seed=7 for balanced coverage
ANCHORS = [
    ("UUU", "Phe",  +0.25, +0.15),   # top-right
    ("GGG", "Gly",  -0.10, +0.20),   # top-left
    ("CAG", "Gln",  -0.25, +0.02),   # left
    ("AAA", "Lys",  -0.10, -0.20),   # bottom-left
    ("UCC", "Ser",  +0.20, -0.02),   # right
]
for codon, aa_name, dx, dy in ANCHORS:
    if codon in pos:
        x, y = pos[codon]
        ax_a.annotate(
            f"{codon}\n{aa_name}",
            (x, y), xytext=(x + dx, y + dy),
            fontsize=6, ha="center", va="center",
            fontweight="bold",
            arrowprops=dict(arrowstyle="-", color="0.4", lw=0.4),
        )

ax_a.set_xlim(-1.25, 1.25)
ax_a.set_ylim(-1.25, 1.25)
ax_a.set_aspect("equal")
ax_a.axis("off")
ax_a.text(-0.05, 1.02, "a", transform=ax_a.transAxes,
          fontsize=11, fontweight="bold", va="bottom", ha="left")
ax_a.text(0.5, -0.02, "codon Hamming graph (61 sense codons)",
          transform=ax_a.transAxes, ha="center", va="top", fontsize=8)

# ═════════════════════════════════════════════════════════════
# panel b — PCA space of 20 AAs
# ═════════════════════════════════════════════════════════════
aa_list = sorted(aa_vec.keys())
pca_cols = [c for c in aa_df.columns if c.startswith("pca")]
aa_sorted_df = aa_df.drop_duplicates("aa").set_index("aa").sort_index()
pc1 = aa_sorted_df[pca_cols[0]].reindex(aa_list).values
pc2 = aa_sorted_df[pca_cols[1]].reindex(aa_list).values

for group in GROUP_ORDER:
    members = [aa for aa in aa_list if RAA10.get(aa) == group]
    idxs = [aa_list.index(m) for m in members]
    ax_b.scatter(pc1[idxs], pc2[idxs], c=GROUP_COLORS[group],
                 s=40, edgecolors="black", linewidths=0.3, zorder=5)

# hand-tuned label offsets (reused from fig1b_pca_space.py)
LABEL_OFFSETS = {
    "I": (-0.1, -0.55), "L": (+0.55, -0.05), "V": (-0.1, -0.55), "M": (-0.55, +0.10),
    "F": (-0.1, -0.50), "W": (-0.1, +0.50), "Y": (-0.55, +0.10),
    "N": (+0.50, +0.20), "D": (+0.50, -0.20), "S": (+0.10, -0.55), "T": (+0.50, +0.15),
    "E": (+0.50, -0.20), "Q": (-0.55, +0.10),
    "K": (-0.50, +0.20), "R": (+0.50, -0.10),
    "H": (-0.50, +0.15), "G": (+0.10, -0.50),
    "A": (-0.50, -0.10), "P": (-0.50, -0.10), "C": (-0.50, +0.10),
}
for i, aa in enumerate(aa_list):
    dx, dy = LABEL_OFFSETS.get(aa, (0.20, 0.20))
    ax_b.annotate(
        aa, (pc1[i], pc2[i]),
        xytext=(pc1[i] + dx, pc2[i] + dy),
        fontsize=6, fontweight="bold", ha="center", va="center",
        arrowprops=dict(arrowstyle="-", color="0.6", lw=0.3),
    )

ax_b.set_xlabel(f"PC1 ({var_frac[0]*100:.1f}% var)", fontsize=9)
ax_b.set_ylabel(f"PC2 ({var_frac[1]*100:.1f}% var)", fontsize=9)
ax_b.spines["top"].set_visible(False)
ax_b.spines["right"].set_visible(False)
ax_b.grid(False)
ax_b.text(-0.20, 1.02, "b", transform=ax_b.transAxes,
          fontsize=11, fontweight="bold", va="bottom", ha="left")

# ═════════════════════════════════════════════════════════════
# panel c — GCU neighborhood with real PC distances
# ═════════════════════════════════════════════════════════════
BASES = ["A", "U", "G", "C"]
center_codon = "GCU"
center_aa = c2a[center_codon]
center_vec = aa_vec[center_aa]

neighbors = []
for pos_i in range(3):
    for b in BASES:
        if b == center_codon[pos_i]:
            continue
        nc = center_codon[:pos_i] + b + center_codon[pos_i+1:]
        naa = c2a.get(nc, "*")
        if naa in aa_vec:
            d = np.linalg.norm(center_vec - aa_vec[naa])
            syn = (naa == center_aa)
        else:
            d = np.nan
            syn = False
        neighbors.append({
            "codon": nc, "aa": naa, "pos": pos_i + 1,
            "dist": d, "syn": syn,
        })

# place neighbors at polar positions — 3 groups of 3, by position
cx, cy = 0.0, 0.0
R = 1.0
R_node = 0.18
# group sectors: pos1 top, pos2 lower-left, pos3 lower-right
sector_centers = {1: 90, 2: 215, 3: 325}
spread = 32  # degrees within sector

def polar(cx, cy, angle_deg, r):
    a = np.radians(angle_deg)
    return cx + r * np.cos(a), cy + r * np.sin(a)

# draw edges with distance labels
for pos_i in [1, 2, 3]:
    sector = sector_centers[pos_i]
    group = [n for n in neighbors if n["pos"] == pos_i]
    angles = [sector + spread, sector, sector - spread]
    for n, ang in zip(group, angles):
        xy = polar(cx, cy, ang, R)
        # edge color: blue for synonymous, red for non-syn
        color = "#4477AA" if n["syn"] else "#CC6677"
        # draw edge from center to node edge
        dx, dy = xy[0] - cx, xy[1] - cy
        dist_r = np.hypot(dx, dy)
        ux, uy = dx / dist_r, dy / dist_r
        x0, y0 = cx + ux * 0.26, cy + uy * 0.26
        x1, y1 = xy[0] - ux * R_node, xy[1] - uy * R_node
        ax_c.plot([x0, x1], [y0, y1], color=color, linewidth=1.0,
                  zorder=1, solid_capstyle="round")
        # midpoint label for ||Δf||
        mx, my = (x0 + x1) / 2, (y0 + y1) / 2
        # offset label perpendicular to edge to avoid overlap
        perp_x, perp_y = -uy, ux
        lab_x = mx + perp_x * 0.10
        lab_y = my + perp_y * 0.10
        ax_c.text(
            lab_x, lab_y, f"{n['dist']:.1f}",
            fontsize=6, ha="center", va="center",
            color=color, fontweight="bold",
            bbox=dict(facecolor="white", edgecolor="none", pad=0.5, alpha=0.85),
            zorder=2,
        )
        # node
        ax_c.add_patch(plt.Circle(xy, R_node, facecolor="white",
                                   edgecolor=color, linewidth=0.8, zorder=3))
        ax_c.text(xy[0], xy[1] - 0.035, n["codon"], fontsize=5.5,
                  ha="center", va="center", fontweight="bold", zorder=4)
        ax_c.text(xy[0], xy[1] + 0.045, n["aa"], fontsize=5,
                  ha="center", va="center", color="0.3", zorder=4)

# center node
ax_c.add_patch(plt.Circle((cx, cy), 0.26, facecolor="white",
                          edgecolor="black", linewidth=1.2, zorder=3))
ax_c.text(cx, cy - 0.05, center_codon, fontsize=7,
          ha="center", va="center", fontweight="bold", zorder=4)
ax_c.text(cx, cy + 0.06, center_aa, fontsize=6,
          ha="center", va="center", color="0.3", zorder=4)

# sector labels
for pos_i in [1, 2, 3]:
    sector = sector_centers[pos_i]
    lx, ly = polar(cx, cy, sector, R + R_node + 0.22)
    label = f"pos {pos_i}" if pos_i < 3 else "pos 3\n(wobble)"
    color = "#4477AA" if pos_i == 3 else "#CC6677"
    ha = "center"
    if pos_i == 2:
        ha = "right"; lx -= 0.08
    elif pos_i == 3:
        ha = "left"; lx += 0.08
    ax_c.text(lx, ly, label, fontsize=7, color=color,
              ha=ha, va="center", fontstyle="italic")

# formula boxes beneath
ax_c.text(
    0, -1.55,
    r"$E = \sum_{(u,v)} \|f(u)-f(v)\|^2$" + "    " +
    r"$D = \frac{1}{|C|}\sum_c \sum_{c' \in N(c)} \|f(c)-f(c')\|$",
    fontsize=6, ha="center", va="center",
)

ax_c.set_xlim(-1.55, 1.55)
ax_c.set_ylim(-1.75, 1.55)
ax_c.set_aspect("equal")
ax_c.axis("off")
ax_c.text(-0.05, 1.02, "c", transform=ax_c.transAxes,
          fontsize=11, fontweight="bold", va="bottom", ha="left")

# ═════════════════════════════════════════════════════════════
# shared legend across bottom
# ═════════════════════════════════════════════════════════════
legend_handles = []
for group in GROUP_ORDER:
    h = plt.Line2D([0], [0], marker="o", color="w",
                   markerfacecolor=GROUP_COLORS[group], markeredgecolor="black",
                   markeredgewidth=0.3, markersize=6, label=group)
    legend_handles.append(h)

# synonymous/nonsynonymous edge handles for panel c
syn_handle = plt.Line2D([0], [0], color="#4477AA", linewidth=1.5,
                         label="synonymous (pos 3)")
nonsyn_handle = plt.Line2D([0], [0], color="#CC6677", linewidth=1.5,
                            label="non-syn edge")

fig.legend(
    handles=legend_handles,
    loc="lower center", bbox_to_anchor=(0.5, 0.005),
    ncol=5, fontsize=7, frameon=False,
    handlelength=1.2, handletextpad=0.4,
    columnspacing=1.0, labelspacing=0.3,
)

# ── save ──
OUT.mkdir(exist_ok=True)
fig.savefig(OUT / "fig1.pdf", dpi=300)
fig.savefig(OUT / "fig1.png", dpi=300)
print(f"saved: {OUT / 'fig1.pdf'}")
print(f"saved: {OUT / 'fig1.png'}")
plt.close(fig)
