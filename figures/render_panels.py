#!/usr/bin/env python3
"""
render_panels.py — emit each manuscript-figure panel as a separate SVG.

Composite figures (`fig1.py` … `figS1_raw_distributions.py`) remain the
canonical submission deliverables. This script produces the same panels
larger and standalone for hand-editing in Inkscape / Illustrator.

Output: figures/panels/fig{N}_{a,b,c,…}.svg

Conventions:
  - SVG with editable text (svg.fonttype = 'none').
  - When multiple panels share a metric (figS1 a/b/c are all D
    distributions for different conditions), they share x-axis range so
    the across-panel comparison is visually faithful.
  - Hamming graph (fig1_a) uses thicker, darker edges than the composite
    so the network structure is visible.

Usage: python figures/render_panels.py
"""

import json, os, sys
from pathlib import Path
from collections import defaultdict

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import networkx as nx

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

OUT = ROOT / "figures" / "panels"
OUT.mkdir(parents=True, exist_ok=True)
RES = ROOT / "results"

from src.io.codon_io import load_codon_table_csv
from src.io.aa_props_io import load_aa_props, aa_prop_vector_map
from src.io.aa_props_lib import AA
from src.sims.codon_graph import build_graph_n
from src.io.aa_props_io import _numeric_cols, _drop_highly_correlated, _zscore

# ── style ──
mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Liberation Sans", "DejaVu Sans"],
    "font.size": 11,
    "axes.labelsize": 12,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "legend.fontsize": 10,
    "axes.linewidth": 0.8,
    "xtick.major.width": 0.8,
    "ytick.major.width": 0.8,
    "xtick.major.size": 4,
    "ytick.major.size": 4,
    # SVG: keep text as text, not paths — editable in Inkscape/Illustrator
    "svg.fonttype": "none",
})

# Paul Tol palette — keep consistent with composite figures
BLUE  = "#4477AA"   # Condition A (doublet)
GREEN = "#228833"   # Condition C (10-AA triplet)
RED   = "#CC6677"   # Condition B (SGC, 20-AA triplet)
GREY  = "#BBBBBB"   # null / scatter cloud
GREY_LIGHT = "#DDDDDD"

# RAA10 grouping for AA chemical class color encoding
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
    "Aromatic": "#44AA99",       "Hydroxylic": "#117733",
    "Amidic": "#999933",         "Acidic": "#CC6677",
    "Basic": "#882255",          "Imidazole": "#AA4499",
    "Cyclic": "#DDCC77",         "Sulfur": "#EE8866",
}

# ── shared data load (reuses pipelines) ──
codon_df = load_codon_table_csv()
aa_df, _ = load_aa_props("data/processed/aa_props.parquet")
aa_vec   = aa_prop_vector_map(aa_df)
c2a      = dict(zip(codon_df["codon"].str.upper(), codon_df["aa"].str.upper()))

# PCA variance explained per component (for axis labels)
rows = [{"aa": a, **p} for a, p in AA.items()]
lib_df = pd.DataFrame(rows)
used_raw = _drop_highly_correlated(lib_df, _numeric_cols(lib_df), 0.95)
sub = lib_df[["aa"] + used_raw].drop_duplicates("aa").set_index("aa").sort_index()
Xz, _, _ = _zscore(sub.to_numpy(dtype=float))
_, S, _  = np.linalg.svd(Xz, full_matrices=False)
var_frac = (S**2) / (S**2).sum()

# null distributions (1M samples each, [D, E] columns)
null_A = np.load(RES / "null_A_doublet_10aa.npy")
null_B = np.load(RES / "null_B_triplet_20aa.npy")
null_C = np.load(RES / "null_C_triplet_10aa.npy")

with open(RES / "sgc_baselines.json") as f:
    sgc_base = json.load(f)
with open(RES / "publication_controls.json") as f:
    pub = json.load(f)

D_A, E_A, D_B, E_B, D_C, E_C = (null_A[:,0], null_A[:,1],
                                 null_B[:,0], null_B[:,1],
                                 null_C[:,0], null_C[:,1])
sgc_DA, sgc_EA = sgc_base["A"]["D"], sgc_base["A"]["E"]
sgc_DB, sgc_EB = sgc_base["B"]["D"], sgc_base["B"]["E"]
sgc_DC, sgc_EC = sgc_base["C"]["D"], sgc_base["C"]["E"]
zDA, zEA = pub["A_10AA_n2"]["D"]["z_score"], pub["A_10AA_n2"]["E"]["z_score"]
zDB, zEB = pub["B_20AA_n3"]["D"]["z_score"], pub["B_20AA_n3"]["E"]["z_score"]
zDC, zEC = pub["C_10AA_n3"]["D"]["z_score"], pub["C_10AA_n3"]["E"]["z_score"]
DE_corr_B = pub["B_20AA_n3"]["DE_correlation"]


def clean(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(False)


def save_svg(fig, name):
    p = OUT / f"{name}.svg"
    fig.savefig(p, format="svg", bbox_inches="tight")
    plt.close(fig)
    print(f"  saved: {p.relative_to(ROOT)}")


# ════════════════════════════════════════════════════════════════════
# Figure 1 — framework
# ════════════════════════════════════════════════════════════════════

def panel_fig1_a():
    """Hamming graph of 61 sense codons — large, edges visible."""
    fig, ax = plt.subplots(figsize=(8, 8))

    G = build_graph_n(3)
    stops = set(codon_df.loc[codon_df["is_stop"], "codon"].str.upper())
    G_sense = G.copy(); G_sense.remove_nodes_from(stops)
    pos = nx.spring_layout(G_sense, seed=7, k=0.35, iterations=300)

    # EDGES MORE VISIBLE: darker, thicker, less transparent
    nx.draw_networkx_edges(
        G_sense, pos, ax=ax,
        edge_color="0.55", width=0.6, alpha=0.85,
    )

    node_colors = []
    for codon in G_sense.nodes():
        aa = c2a.get(codon, "*")
        grp = RAA10.get(aa, None)
        node_colors.append(GROUP_COLORS.get(grp, "#CCCCCC"))

    nx.draw_networkx_nodes(
        G_sense, pos, ax=ax, node_size=130, node_color=node_colors,
        edgecolors="black", linewidths=0.6,
    )

    # peripheral anchor labels
    ANCHORS = [
        ("UUU", "Phe",  +0.25, +0.15),
        ("GGG", "Gly",  -0.10, +0.20),
        ("CAG", "Gln",  -0.25, +0.02),
        ("AAA", "Lys",  -0.10, -0.20),
        ("UCC", "Ser",  +0.20, -0.02),
    ]
    for codon, aa_name, dx, dy in ANCHORS:
        if codon in pos:
            x, y = pos[codon]
            ax.annotate(
                f"{codon}\n{aa_name}", (x, y),
                xytext=(x + dx, y + dy),
                fontsize=9, ha="center", va="center", fontweight="bold",
                arrowprops=dict(arrowstyle="-", color="0.3", lw=0.6),
            )

    ax.set_xlim(-1.30, 1.30); ax.set_ylim(-1.30, 1.30)
    ax.set_aspect("equal"); ax.axis("off")

    # legend below
    handles = [Line2D([0], [0], marker="o", color="w",
                      markerfacecolor=GROUP_COLORS[g], markeredgecolor="black",
                      markeredgewidth=0.4, markersize=9, label=g)
               for g in GROUP_ORDER]
    ax.legend(handles=handles, loc="upper center",
              bbox_to_anchor=(0.5, -0.02), ncol=5, frameon=False,
              handlelength=1.2, handletextpad=0.4, columnspacing=1.2)

    fig.suptitle("Codon Hamming graph (61 sense codons, 263 edges)",
                 fontsize=12, y=0.98)
    save_svg(fig, "fig1_a_hamming_graph")


def panel_fig1_b():
    """20 amino acids in PC1-PC2 space."""
    fig, ax = plt.subplots(figsize=(7, 6))

    aa_list = sorted(aa_vec.keys())
    pca_cols = [c for c in aa_df.columns if c.startswith("pca")]
    aa_sorted_df = aa_df.drop_duplicates("aa").set_index("aa").sort_index()
    pc1 = aa_sorted_df[pca_cols[0]].reindex(aa_list).values
    pc2 = aa_sorted_df[pca_cols[1]].reindex(aa_list).values

    for grp in GROUP_ORDER:
        members = [a for a in aa_list if RAA10.get(a) == grp]
        idxs = [aa_list.index(m) for m in members]
        ax.scatter(pc1[idxs], pc2[idxs], c=GROUP_COLORS[grp],
                   s=120, edgecolors="black", linewidths=0.5,
                   zorder=5, label=grp)

    OFF = {  # tuned offsets (from fig1.py)
        "I": (-0.1,-0.55), "L": (+0.55,-0.05), "V": (-0.1,-0.55), "M": (-0.55,+0.10),
        "F": (-0.1,-0.50), "W": (-0.1,+0.50), "Y": (-0.55,+0.10),
        "N": (+0.50,+0.20), "D": (+0.50,-0.20), "S": (+0.10,-0.55), "T": (+0.50,+0.15),
        "E": (+0.50,-0.20), "Q": (-0.55,+0.10),
        "K": (-0.50,+0.20), "R": (+0.50,-0.10),
        "H": (-0.50,+0.15), "G": (+0.10,-0.50),
        "A": (-0.50,-0.10), "P": (-0.50,-0.10), "C": (-0.50,+0.10),
    }
    for i, aa in enumerate(aa_list):
        dx, dy = OFF.get(aa, (0.20, 0.20))
        ax.annotate(aa, (pc1[i], pc2[i]),
                    xytext=(pc1[i]+dx, pc2[i]+dy),
                    fontsize=10, fontweight="bold",
                    ha="center", va="center",
                    arrowprops=dict(arrowstyle="-", color="0.55", lw=0.4))

    ax.set_xlabel(f"PC1 ({var_frac[0]*100:.1f} % variance)")
    ax.set_ylabel(f"PC2 ({var_frac[1]*100:.1f} % variance)")
    clean(ax)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.13),
              ncol=5, frameon=False, fontsize=8,
              handlelength=1.0, handletextpad=0.3, columnspacing=1.0,
              markerscale=0.7)
    fig.suptitle("Amino-acid property space (20 AAs, 8 PCs from 22 descriptors)",
                 fontsize=12, y=0.98)
    save_svg(fig, "fig1_b_pca_space")


def panel_fig1_c():
    """GCU Hamming-1 neighborhood with real ||Δf|| edge labels."""
    fig, ax = plt.subplots(figsize=(7, 7))

    BASES = ["A", "U", "G", "C"]
    center = "GCU"; center_aa = c2a[center]; cv = aa_vec[center_aa]
    neighbors = []
    for pos_i in range(3):
        for b in BASES:
            if b == center[pos_i]: continue
            nc  = center[:pos_i] + b + center[pos_i+1:]
            naa = c2a.get(nc, "*")
            if naa in aa_vec:
                d = float(np.linalg.norm(cv - aa_vec[naa]))
                syn = (naa == center_aa)
            else:
                d, syn = float("nan"), False
            neighbors.append({"codon": nc, "aa": naa, "pos": pos_i+1,
                              "dist": d, "syn": syn})

    cx, cy = 0.0, 0.0; R = 1.0; R_node = 0.20
    sec = {1: 90, 2: 215, 3: 325}; spread = 32

    def polar(cx, cy, ang_deg, r):
        a = np.radians(ang_deg)
        return cx + r*np.cos(a), cy + r*np.sin(a)

    for pos_i in [1, 2, 3]:
        a_c = sec[pos_i]
        group = [n for n in neighbors if n["pos"] == pos_i]
        angles = [a_c + spread, a_c, a_c - spread]
        for n, ang in zip(group, angles):
            xy = polar(cx, cy, ang, R)
            color = BLUE if n["syn"] else RED
            ux = (xy[0]-cx) / R; uy = (xy[1]-cy) / R
            x0, y0 = cx + ux*0.28, cy + uy*0.28
            x1, y1 = xy[0] - ux*R_node, xy[1] - uy*R_node
            ax.plot([x0, x1], [y0, y1], color=color, linewidth=1.5,
                    zorder=1, solid_capstyle="round")
            mx, my = (x0+x1)/2, (y0+y1)/2
            perp_x, perp_y = -uy, ux
            ax.text(mx + perp_x*0.10, my + perp_y*0.10,
                    f"{n['dist']:.1f}",
                    fontsize=10, ha="center", va="center",
                    color=color, fontweight="bold",
                    bbox=dict(facecolor="white", edgecolor="none",
                              pad=1.0, alpha=0.9), zorder=2)
            ax.add_patch(plt.Circle(xy, R_node, facecolor="white",
                                    edgecolor=color, linewidth=1.0, zorder=3))
            ax.text(xy[0], xy[1]-0.045, n["codon"], fontsize=9,
                    ha="center", va="center", fontweight="bold", zorder=4)
            ax.text(xy[0], xy[1]+0.06, n["aa"], fontsize=8,
                    ha="center", va="center", color="0.3", zorder=4)

    ax.add_patch(plt.Circle((cx, cy), 0.28, facecolor="white",
                            edgecolor="black", linewidth=1.5, zorder=3))
    ax.text(cx, cy-0.06, center, fontsize=11,
            ha="center", va="center", fontweight="bold", zorder=4)
    ax.text(cx, cy+0.07, center_aa, fontsize=10,
            ha="center", va="center", color="0.3", zorder=4)

    for pos_i in [1, 2, 3]:
        a = sec[pos_i]
        lx, ly = polar(cx, cy, a, R + R_node + 0.28)
        label = f"pos {pos_i}" if pos_i < 3 else "pos 3 (wobble)"
        color = BLUE if pos_i == 3 else RED
        ha = "center"
        if pos_i == 2: ha = "right"; lx -= 0.06
        elif pos_i == 3: ha = "left"; lx += 0.06
        ax.text(lx, ly, label, fontsize=11, color=color,
                ha=ha, va="center", fontstyle="italic")

    ax.text(0, -1.65,
            r"$E = \sum_{(u,v)} \|f(u)-f(v)\|^2$" +
            "        " +
            r"$D = \frac{1}{|C|}\sum_c \sum_{c' \in N(c)} \|f(c)-f(c')\|$",
            fontsize=10, ha="center", va="center")
    ax.set_xlim(-1.6, 1.6); ax.set_ylim(-1.85, 1.6)
    ax.set_aspect("equal"); ax.axis("off")
    fig.suptitle("GCU (Ala) Hamming-1 neighborhood with ‖Δf‖ in PC space",
                 fontsize=12, y=0.97)
    save_svg(fig, "fig1_c_gcu_schematic")


# ════════════════════════════════════════════════════════════════════
# Figure 2 — SGC null
# ════════════════════════════════════════════════════════════════════

def panel_fig2_a():
    """Dirichlet energy null histogram (Condition B)."""
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.hist(E_B, bins=120, color=GREY, edgecolor="black", linewidth=0.2)
    ax.axvline(sgc_EB, color=RED, linewidth=1.6, zorder=5)
    ymax = ax.get_ylim()[1]
    ax.annotate(
        f"SGC\n$z = {zEB:.1f}$\n$p < 10^{{-6}}$",
        xy=(sgc_EB, ymax*0.72),
        xytext=(sgc_EB + (E_B.max()-E_B.min())*0.05, ymax*0.72),
        color=RED, fontsize=12, fontweight="bold",
        ha="left", va="center",
        arrowprops=dict(arrowstyle="->", color=RED, lw=0.8),
    )
    ax.set_xlabel("Dirichlet energy $E$")
    ax.set_ylabel("Count (of 10$^6$ random codes)")
    clean(ax)
    fig.suptitle("Condition B: 20-AA SGC, Dirichlet energy null",
                 fontsize=12, y=0.97)
    save_svg(fig, "fig2_a_E_histogram")


def panel_fig2_b():
    """Noise distortion null histogram (Condition B)."""
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.hist(D_B, bins=120, color=GREY, edgecolor="black", linewidth=0.2)
    ax.axvline(sgc_DB, color=RED, linewidth=1.6, zorder=5)
    ymax = ax.get_ylim()[1]
    ax.annotate(
        f"SGC\n$z = {zDB:.1f}$\n$p < 10^{{-6}}$",
        xy=(sgc_DB, ymax*0.72),
        xytext=(sgc_DB + (D_B.max()-D_B.min())*0.05, ymax*0.72),
        color=RED, fontsize=12, fontweight="bold",
        ha="left", va="center",
        arrowprops=dict(arrowstyle="->", color=RED, lw=0.8),
    )
    ax.set_xlabel("Noise distortion $D$")
    ax.set_ylabel("Count (of 10$^6$ random codes)")
    clean(ax)
    fig.suptitle("Condition B: 20-AA SGC, noise distortion null",
                 fontsize=12, y=0.97)
    save_svg(fig, "fig2_b_D_histogram")


def panel_fig2_c():
    """Joint D-E scatter for Condition B with marginal histograms."""
    fig = plt.figure(figsize=(7, 7))
    gs = GridSpec(2, 2, figure=fig,
                  width_ratios=[4, 1], height_ratios=[1, 4],
                  hspace=0.04, wspace=0.04,
                  left=0.13, right=0.97, bottom=0.10, top=0.93)
    ax_s = fig.add_subplot(gs[1, 0])
    ax_x = fig.add_subplot(gs[0, 0], sharex=ax_s)
    ax_y = fig.add_subplot(gs[1, 1], sharey=ax_s)

    rng = np.random.default_rng(42)
    idx = rng.choice(len(D_B), size=10000, replace=False)
    ax_s.scatter(D_B[idx], E_B[idx], s=8, c=GREY, alpha=0.30,
                 edgecolors="none", rasterized=False)
    ax_s.plot(sgc_DB, sgc_EB, marker="*", markersize=22, color=RED,
              markeredgecolor="black", markeredgewidth=0.8, zorder=10)
    ax_s.annotate(
        f"SGC\n$z_D = {zDB:.1f}$\n$z_E = {zEB:.1f}$",
        xy=(sgc_DB, sgc_EB),
        xytext=(sgc_DB + 1.8, sgc_EB + 250),
        fontsize=11, color=RED, fontweight="bold",
        arrowprops=dict(arrowstyle="->", color=RED, lw=0.8),
        ha="left", va="bottom",
    )
    ax_s.text(0.97, 0.04, f"$r = {DE_corr_B:.3f}$",
              transform=ax_s.transAxes, fontsize=11,
              ha="right", va="bottom", color="0.3")
    ax_s.set_xlabel("Noise distortion $D$")
    ax_s.set_ylabel("Dirichlet energy $E$")
    clean(ax_s)

    ax_x.hist(D_B, bins=140, color=GREY, edgecolor="none", alpha=0.85)
    ax_x.axvline(sgc_DB, color=RED, linewidth=1.2)
    ax_x.tick_params(labelbottom=False, left=False, bottom=False, labelleft=False)
    for s in ["top","right","left","bottom"]: ax_x.spines[s].set_visible(False)

    ax_y.hist(E_B, bins=140, orientation="horizontal",
              color=GREY, edgecolor="none", alpha=0.85)
    ax_y.axhline(sgc_EB, color=RED, linewidth=1.2)
    ax_y.tick_params(labelleft=False, left=False, bottom=False, labelbottom=False)
    for s in ["top","right","left","bottom"]: ax_y.spines[s].set_visible(False)

    fig.suptitle("Condition B: joint (D, E) distribution; SGC marked",
                 fontsize=12, y=0.97)
    save_svg(fig, "fig2_c_joint_scatter")


# ════════════════════════════════════════════════════════════════════
# Figure 3 — factorial + synonymy
# ════════════════════════════════════════════════════════════════════

def panel_fig3_a():
    """Three null distributions of D on a shared x-axis (overlay)."""
    fig, ax = plt.subplots(figsize=(8, 5))
    bins = np.linspace(17.5, 35.5, 140)
    for data, color, label, sgc_v, z in [
        (D_A, BLUE,  f"A: 10 AA, $n$=2  ($z_D = {zDA:.1f}$)", sgc_DA, zDA),
        (D_C, GREEN, f"C: 10 AA, $n$=3  ($z_D = {zDC:.1f}$)", sgc_DC, zDC),
        (D_B, RED,   f"B: 20 AA, $n$=3 (SGC)  ($z_D = {zDB:.1f}$)", sgc_DB, zDB),
    ]:
        ax.hist(data, bins=bins, color=color, alpha=0.55,
                edgecolor="none", label=label, zorder=2)
        ax.axvline(sgc_v, color=color, linewidth=1.8, zorder=5)

    ax.text(0.98, 0.96,
            "A lies inside its null tail\nB and C lie beyond",
            transform=ax.transAxes, fontsize=10, ha="right", va="top",
            fontstyle="italic", color="0.3")

    handles = [
        (Patch(facecolor=BLUE,  alpha=0.55), Line2D([0],[0], color=BLUE,  lw=1.8)),
        (Patch(facecolor=GREEN, alpha=0.55), Line2D([0],[0], color=GREEN, lw=1.8)),
        (Patch(facecolor=RED,   alpha=0.55), Line2D([0],[0], color=RED,   lw=1.8)),
    ]
    labels = [
        f"A: 10 AA, $n$=2   ($z_D = {zDA:.1f}$)",
        f"C: 10 AA, $n$=3   ($z_D = {zDC:.1f}$)",
        f"B: 20 AA, $n$=3 (SGC)  ($z_D = {zDB:.1f}$)",
    ]
    ax.legend(handles, labels, loc="upper left", frameon=False, fontsize=9,
              handler_map={tuple: mpl.legend_handler.HandlerTuple(ndivide=None)},
              handlelength=2.4, handletextpad=0.5)
    ax.set_xlabel("Noise distortion $D$")
    ax.set_ylabel("Count (of 10$^6$ random codes)")
    ax.set_xlim(17.5, 36)
    clean(ax)
    fig.suptitle("Three null distributions on a shared D-axis",
                 fontsize=12, y=0.97)
    save_svg(fig, "fig3_a_factorial_overlay")


def panel_fig3_b():
    """Position-specific synonymy bars."""
    fig, ax = plt.subplots(figsize=(7, 5))
    syn = {"A": [33.3, 8.3, None], "B": [4.4, 0.0, 68.9], "C": [26.2, 10.9, 76.5]}
    x = np.arange(3); bw = 0.27
    bars_A = [syn["A"][0], syn["A"][1], 0]
    bars_B = [syn["B"][0], syn["B"][1], syn["B"][2]]
    bars_C = [syn["C"][0], syn["C"][1], syn["C"][2]]
    ax.bar(x - bw, bars_A, bw, color=BLUE,  edgecolor="none", label="A: 10 AA, $n$=2")
    ax.bar(x,      bars_B, bw, color=RED,   edgecolor="none", label="B: 20 AA, $n$=3 (SGC)")
    ax.bar(x + bw, bars_C, bw, color=GREEN, edgecolor="none", label="C: 10 AA, $n$=3")
    ax.text(2 - bw, 2.0, "N/A", ha="center", va="bottom",
            fontsize=10, color="black")
    ax.set_xticks(x); ax.set_xticklabels(["pos 1", "pos 2", "pos 3"])
    ax.set_ylabel("Synonymous fraction (%)")
    ax.set_ylim(0, 95)
    clean(ax)
    ax.legend(loc="upper left", frameon=False, fontsize=10,
              handlelength=1.5, handletextpad=0.5)
    fig.suptitle("Position-specific synonymous-substitution fraction",
                 fontsize=12, y=0.97)
    save_svg(fig, "fig3_b_synonymy_bars")


# ════════════════════════════════════════════════════════════════════
# Figure 4 — wobble mechanism
# ════════════════════════════════════════════════════════════════════

def _draw_wobble_panel(ax, center_codon, center_aa, sectors, title):
    """Draw a Hamming-1 neighborhood schematic. sectors: list of
    (pos_index_int_or_str, sector_angle_deg, list_of_(codon, aa, sub_pos, syn))."""
    R_center = 0.42; R_outer = 0.36; orbit = 1.85; spread = 30

    def polar(cx, cy, ang, r):
        a = np.radians(ang)
        return cx + r*np.cos(a), cy + r*np.sin(a)

    def draw_node(xy, codon, sub_pos, aa, fill, edge, sub_color):
        ax.add_patch(plt.Circle(xy, R_outer, facecolor=fill,
                                edgecolor=edge, linewidth=1.0, zorder=3,
                                clip_on=False))
        char_w = 0.13; n = len(codon)
        sx = xy[0] - char_w*n/2 + char_w/2
        for i, ch in enumerate(codon):
            ax.text(sx + i*char_w, xy[1] + 0.06, ch,
                    ha="center", va="center",
                    fontsize=10 if i != sub_pos else 12,
                    fontweight="bold",
                    color="black" if i != sub_pos else sub_color,
                    zorder=4)
        ax.text(xy[0], xy[1] - 0.10, aa, ha="center", va="center",
                fontsize=8, color="0.3", zorder=4)

    def draw_edge(xy_from, xy_to, color, lw=1.6):
        dx = xy_to[0] - xy_from[0]; dy = xy_to[1] - xy_from[1]
        d = np.hypot(dx, dy)
        if d == 0: return
        ux, uy = dx/d, dy/d
        x0, y0 = xy_from[0] + ux*R_center, xy_from[1] + uy*R_center
        x1, y1 = xy_to[0] - ux*R_outer, xy_to[1] - uy*R_outer
        ax.plot([x0, x1], [y0, y1], color=color, linewidth=lw,
                zorder=1, solid_capstyle="round")

    BLUE_FILL = (*mpl.colors.to_rgb(BLUE), 0.18)
    RED_FILL  = (*mpl.colors.to_rgb(RED),  0.18)

    for sec_label, sec_angle, group in sectors:
        n = len(group)
        if n == 1:    angles = [sec_angle]
        elif n == 3:  angles = [sec_angle + spread, sec_angle, sec_angle - spread]
        else:
            angles = [sec_angle - spread*(n-1)/2 + i*spread for i in range(n)]
        for (codon, aa, sp, syn), ang in zip(group, angles):
            xy = polar(0, 0, ang, orbit)
            color = BLUE if syn else RED
            fill  = BLUE_FILL if syn else RED_FILL
            draw_edge((0,0), xy, color)
            draw_node(xy, codon, sp, aa, fill, "black", color)
        # label sector
        lx, ly = polar(0, 0, sec_angle, orbit + R_outer + 0.30)
        ax.text(lx, ly, sec_label, fontsize=11, color="0.3",
                ha="center", va="center", fontstyle="italic")

    # center
    ax.add_patch(plt.Circle((0,0), R_center, facecolor="white",
                            edgecolor="black", linewidth=1.6, zorder=3))
    ax.text(0, -0.06, center_codon, ha="center", va="center",
            fontsize=12, fontweight="bold", zorder=4)
    ax.text(0, +0.08, center_aa, ha="center", va="center",
            fontsize=10, color="0.3", zorder=4)

    lim = orbit + R_outer + 1.05
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
    ax.set_aspect("equal"); ax.axis("off")
    return lim


def panel_fig4_a():
    """Triplet GCU neighborhood — 9 neighbors, pos 3 wobble."""
    fig, ax = plt.subplots(figsize=(7, 7))
    pos1 = [("ACU","Thr",0,False),("UCU","Ser",0,False),("CCU","Pro",0,False)]
    pos2 = [("GGU","Gly",1,False),("GUU","Val",1,False),("GAU","Asp",1,False)]
    pos3 = [("GCG","Ala",2,True), ("GCA","Ala",2,True), ("GCC","Ala",2,True)]
    sectors = [("pos 1",  90, pos1),
               ("pos 2", 215, pos2),
               ("pos 3 (wobble)", 325, pos3)]
    lim = _draw_wobble_panel(ax, "GCU", "Ala", sectors,
                             "Triplet: GCU (Ala) neighborhood")
    ax.text(0,  lim - 0.25, "Positions 1-2: all 6 non-synonymous",
            fontsize=11, ha="center", va="top", color=RED, fontweight="bold")
    ax.text(0, -lim + 0.40, "Position 3 (wobble): all 3 synonymous",
            fontsize=11, ha="center", va="bottom", color=BLUE, fontweight="bold")
    ax.text(0, -lim + 0.10, "total synonymous: 3/9 (33 %)",
            fontsize=10, ha="center", va="bottom", color="black")
    fig.suptitle("Triplet code: GCU Hamming-1 neighborhood",
                 fontsize=12, y=0.98)
    save_svg(fig, "fig4_a_triplet_GCU")


def panel_fig4_b():
    """Doublet GC neighborhood — 6 neighbors, no wobble."""
    fig, ax = plt.subplots(figsize=(7, 7))
    pos1 = [("AC","Ser/Thr",0,False),("UC","Ser/Thr",0,False),("CC","Pro",0,False)]
    pos2 = [("GA","Asp/Glu",1,False),("GU","VLIM",1,False),("GG","Ala/Gly",1,True)]
    sectors = [("pos 1",  90, pos1),
               ("pos 2 (mixed)", 270, pos2)]
    lim = _draw_wobble_panel(ax, "GC", "Ala/Gly", sectors,
                             "Doublet: GC (Ala/Gly) neighborhood")
    ax.text(0, -lim + 0.10,
            "1/6 synonymous (no dedicated position)",
            fontsize=11, ha="center", va="bottom", color=BLUE, fontweight="bold")
    fig.suptitle("Doublet code: GC Hamming-1 neighborhood",
                 fontsize=12, y=0.98)
    save_svg(fig, "fig4_b_doublet_GC")


# ════════════════════════════════════════════════════════════════════
# Figure S1 — raw distributions (SHARED X-AXIS for honest comparison)
# ════════════════════════════════════════════════════════════════════

# global x-range covering all three D nulls + their SGC marks
SHARED_D_LO = 17.0
SHARED_D_HI = 36.0


def _figS1_panel(data, sgc_v, z, color, condition_label, panel_letter, fname):
    fig, ax = plt.subplots(figsize=(8, 4.5))
    # density (not count) so panels with different sample counts compare fairly
    ax.hist(data, bins=120, color=color, alpha=0.70,
            edgecolor="black", linewidth=0.25, density=True)
    ax.axvline(sgc_v, color=RED, linewidth=1.6, zorder=5,
               label=f"SGC observed (z = {z:.1f})")
    ax.axvline(data.min(), color="grey", linewidth=1.0,
               linestyle="--", zorder=4,
               label=f"null minimum = {data.min():.2f}")

    ax.set_xlim(SHARED_D_LO, SHARED_D_HI)   # ← same x-range across panels
    ax.set_xlabel("Noise distortion $D$")
    ax.set_ylabel("Density")
    clean(ax)
    ax.legend(loc="upper right", frameon=False, fontsize=9)
    fig.suptitle(f"{panel_letter}: {condition_label}",
                 fontsize=12, y=0.97)
    save_svg(fig, fname)


def panel_figS1_a():
    _figS1_panel(D_A, sgc_DA, zDA, BLUE,
                 "Condition A — 10 AA, $n$ = 2 (doublet)", "A",
                 "figS1_a_condA_doublet")


def panel_figS1_b():
    _figS1_panel(D_C, sgc_DC, zDC, GREEN,
                 "Condition C — 10 AA, $n$ = 3 (RAA10 triplet)", "B",
                 "figS1_b_condC_10AA_triplet")


def panel_figS1_c():
    _figS1_panel(D_B, sgc_DB, zDB, RED,
                 "Condition B — 20 AA, $n$ = 3 (SGC)", "C",
                 "figS1_c_condB_SGC")


# ════════════════════════════════════════════════════════════════════
def main():
    print(f"writing SVG panels to: {OUT.relative_to(ROOT)}")
    print(f"\nFigure 1 (framework):")
    panel_fig1_a()
    panel_fig1_b()
    panel_fig1_c()
    print(f"\nFigure 2 (SGC null):")
    panel_fig2_a()
    panel_fig2_b()
    panel_fig2_c()
    print(f"\nFigure 3 (factorial + synonymy):")
    panel_fig3_a()
    panel_fig3_b()
    print(f"\nFigure 4 (wobble mechanism):")
    panel_fig4_a()
    panel_fig4_b()
    print(f"\nFigure S1 (raw distributions, shared x-axis):")
    panel_figS1_a()
    panel_figS1_b()
    panel_figS1_c()
    print(f"\ndone. {len(list(OUT.glob('*.svg')))} SVG files in "
          f"{OUT.relative_to(ROOT)}/")


if __name__ == "__main__":
    main()
