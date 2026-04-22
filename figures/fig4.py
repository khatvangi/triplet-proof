#!/usr/bin/env python3
"""
Figure 4 — wobble mechanism.

two panels (single-column, 3.5 x 5.0 in):
  a) triplet GCU neighborhood — 9 Hamming-1 neighbors, pos 1-2 in red
     (all nonsynonymous), pos 3 in blue (all synonymous).
     "Positions 1-2: all nonsynonymous", "Position 3 (wobble): all synonymous"
  b) doublet GC neighborhood — 6 neighbors, only 1/6 synonymous, no
     dedicated synonymy position.

substituted nucleotide in each neighbor codon is underlined.

outputs:
  figures/fig4.pdf, fig4.png (300 dpi)
"""

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "figures"

# ── style ──
mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Liberation Sans", "DejaVu Sans"],
    "font.size": 9,
    "axes.labelsize": 10,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})

BLUE = "#4477AA"   # synonymous / pos 3
RED = "#CC6677"    # non-synonymous / pos 1, 2
BLUE_FILL = (*mpl.colors.to_rgb(BLUE), 0.15)
RED_FILL = (*mpl.colors.to_rgb(RED), 0.15)

NODE_FONT = 6
LABEL_SIZE = 8
ANNOT_SIZE = 7.5


def polar(cx, cy, angle_deg, r):
    a = np.radians(angle_deg)
    return cx + r * np.cos(a), cy + r * np.sin(a)


def draw_node(ax, xy, r, codon, sub_pos, aa, fill, edge, sub_color=None, lw=1.0):
    """draw circle node with codon (substituted base in sub_color, bold) and AA below.

    sub_pos = -1 means no substitution (use for center node).
    sub_color defaults to the edge color so substituted base matches its edge.
    """
    if sub_color is None:
        sub_color = edge
    circ = plt.Circle(xy, r, facecolor=fill, edgecolor=edge,
                      linewidth=lw, zorder=3, clip_on=False)
    ax.add_patch(circ)

    # draw each codon character separately so the substituted one gets its color
    # char width in data coords, approximate
    char_w = 0.12
    n = len(codon)
    total_w = char_w * n
    start_x = xy[0] - total_w / 2 + char_w / 2
    y = xy[1] + 0.05
    for i, ch in enumerate(codon):
        if i == sub_pos:
            color = sub_color
            weight = "bold"
            fontsize = NODE_FONT + 1.5
        else:
            color = "black"
            weight = "bold"
            fontsize = NODE_FONT + 0.5
        ax.text(start_x + i * char_w, y, ch,
                ha="center", va="center",
                fontsize=fontsize, fontweight=weight, color=color, zorder=4)

    ax.text(xy[0], xy[1] - 0.07, aa,
            ha="center", va="center", fontsize=NODE_FONT - 0.5,
            color="0.3", zorder=4)


def draw_edge(ax, xy_from, xy_to, r_from, r_to, color, lw=1.0):
    dx = xy_to[0] - xy_from[0]
    dy = xy_to[1] - xy_from[1]
    d = np.hypot(dx, dy)
    if d == 0:
        return
    ux, uy = dx / d, dy / d
    x0 = xy_from[0] + ux * r_from
    y0 = xy_from[1] + uy * r_from
    x1 = xy_to[0] - ux * r_to
    y1 = xy_to[1] - uy * r_to
    ax.plot([x0, x1], [y0, y1], color=color, linewidth=lw, zorder=1,
            solid_capstyle="round")


# ═════════════════════════════════════════════════════════════
# panel a — triplet GCU neighborhood
# ═════════════════════════════════════════════════════════════

def draw_panel_a(ax):
    cx, cy = 0.0, 0.0
    R_center = 0.38
    R_outer = 0.32
    orbit = 1.75
    spread = 30

    # 3 sector groups: pos 1 (top), pos 2 (lower-left), pos 3 (lower-right)
    # each (codon, aa, sub_position_index)
    pos1 = [("ACU", "Thr", 0), ("UCU", "Ser", 0), ("CCU", "Pro", 0)]
    pos2 = [("GGU", "Gly", 1), ("GUU", "Val", 1), ("GAU", "Asp", 1)]
    pos3 = [("GCG", "Ala", 2), ("GCA", "Ala", 2), ("GCC", "Ala", 2)]

    sec1, sec2, sec3 = 90, 215, 325
    ang1 = [sec1 + spread, sec1, sec1 - spread]
    ang2 = [sec2 + spread, sec2, sec2 - spread]
    ang3 = [sec3 + spread, sec3, sec3 - spread]

    for (codon, aa, sp), ang in zip(pos1, ang1):
        xy = polar(cx, cy, ang, orbit)
        draw_edge(ax, (cx, cy), xy, R_center, R_outer, RED)
        draw_node(ax, xy, R_outer, codon, sp, aa, RED_FILL, "black",
                  sub_color=RED)

    for (codon, aa, sp), ang in zip(pos2, ang2):
        xy = polar(cx, cy, ang, orbit)
        draw_edge(ax, (cx, cy), xy, R_center, R_outer, RED)
        draw_node(ax, xy, R_outer, codon, sp, aa, RED_FILL, "black",
                  sub_color=RED)

    for (codon, aa, sp), ang in zip(pos3, ang3):
        xy = polar(cx, cy, ang, orbit)
        draw_edge(ax, (cx, cy), xy, R_center, R_outer, BLUE)
        draw_node(ax, xy, R_outer, codon, sp, aa, BLUE_FILL, "black",
                  sub_color=BLUE)

    # center node (no substitution)
    draw_node(ax, (cx, cy), R_center, "GCU", -1, "Ala",
              "white", "black", lw=1.5)

    # sector labels (outside outer ring)
    lx, ly = polar(cx, cy, sec1, orbit + R_outer + 0.25)
    ax.text(lx, ly, "pos 1", fontsize=ANNOT_SIZE, color=RED,
            ha="center", va="bottom", fontstyle="italic")
    lx, ly = polar(cx, cy, sec2, orbit + R_outer + 0.25)
    ax.text(lx - 0.10, ly, "pos 2", fontsize=ANNOT_SIZE, color=RED,
            ha="right", va="center", fontstyle="italic")
    lx, ly = polar(cx, cy, sec3, orbit + R_outer + 0.25)
    ax.text(lx + 0.10, ly, "pos 3\n(wobble)",
            fontsize=ANNOT_SIZE, color=BLUE,
            ha="left", va="center", fontstyle="italic", linespacing=1.0)

    # summary brackets — raise the top one above "pos 1" label
    ax.text(0.0, orbit + R_outer + 0.85,
            "Positions 1–2: all nonsynonymous (6/6)",
            fontsize=ANNOT_SIZE, ha="center", va="bottom", color=RED,
            fontweight="bold")
    ax.text(0.0, -(orbit + R_outer + 0.45),
            "Position 3 (wobble): all synonymous (3/3)",
            fontsize=ANNOT_SIZE, ha="center", va="top", color=BLUE,
            fontweight="bold")
    ax.text(0.0, -(orbit + R_outer + 0.85),
            "total synonymous: 3/9",
            fontsize=ANNOT_SIZE - 0.5, ha="center", va="top", color="black")

    lim = orbit + R_outer + 1.15
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-(orbit + R_outer + 1.05), lim)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title("Triplet: GCU (Ala) neighborhood", fontsize=10,
                 fontweight="bold", pad=8)


# ═════════════════════════════════════════════════════════════
# panel b — doublet GC neighborhood
# ═════════════════════════════════════════════════════════════

def draw_panel_b(ax):
    cx, cy = 0.0, 0.0
    R_center = 0.38
    R_outer = 0.34
    orbit = 1.60
    spread = 45

    # pos 1 (top), pos 2 (bottom)
    # GG maps to Ala/Gly under majority-vote projection → synonymous with GC
    pos1 = [("AC", "Ser/Thr", 0, False),
            ("UC", "Ser/Thr", 0, False),
            ("CC", "Pro",     0, False)]
    pos2 = [("GA", "Asp/Glu", 1, False),
            ("GU", "VLIM",    1, False),
            ("GG", "Ala/Gly", 1, True)]

    sec1, sec2 = 90, 270
    ang1 = [sec1 + spread, sec1, sec1 - spread]
    ang2 = [sec2 + spread, sec2, sec2 - spread]

    for (codon, aa, sp, syn), ang in zip(pos1, ang1):
        xy = polar(cx, cy, ang, orbit)
        c = BLUE if syn else RED
        fill = BLUE_FILL if syn else RED_FILL
        draw_edge(ax, (cx, cy), xy, R_center, R_outer, c)
        draw_node(ax, xy, R_outer, codon, sp, aa, fill, "black",
                  sub_color=c)

    for (codon, aa, sp, syn), ang in zip(pos2, ang2):
        xy = polar(cx, cy, ang, orbit)
        c = BLUE if syn else RED
        fill = BLUE_FILL if syn else RED_FILL
        draw_edge(ax, (cx, cy), xy, R_center, R_outer, c)
        draw_node(ax, xy, R_outer, codon, sp, aa, fill, "black",
                  sub_color=c)

    # center node (no substitution)
    draw_node(ax, (cx, cy), R_center, "GC", -1, "Ala/Gly",
              "white", "black", lw=1.5)

    # sector labels — move pos 2 slightly below to avoid node overlap
    lx, ly = polar(cx, cy, sec1, orbit + R_outer + 0.25)
    ax.text(lx, ly, "pos 1", fontsize=ANNOT_SIZE, color=RED,
            ha="center", va="bottom", fontstyle="italic")
    lx, ly = polar(cx, cy, sec2, orbit + R_outer + 0.30)
    # pos 2 is mixed: GG is synonymous (blue edge), GA/GU are not (red edges)
    ax.text(lx, ly, "pos 2 (mixed)", fontsize=ANNOT_SIZE,
            color="0.3",
            ha="center", va="top", fontstyle="italic")

    # summary
    ax.text(0.0, -(orbit + R_outer + 0.70),
            "1/6 synonymous (no dedicated position)",
            fontsize=ANNOT_SIZE, ha="center", va="top", color=BLUE,
            fontweight="bold")

    lim = orbit + R_outer + 1.00
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-(orbit + R_outer + 1.00), lim)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title("Doublet: GC (Ala/Gly) neighborhood", fontsize=10,
                 fontweight="bold", pad=8)


# ═════════════════════════════════════════════════════════════

def main():
    fig, (ax_a, ax_b) = plt.subplots(2, 1, figsize=(3.5, 5.5))
    fig.subplots_adjust(hspace=0.28, top=0.94, bottom=0.08,
                        left=0.02, right=0.98)

    draw_panel_a(ax_a)
    draw_panel_b(ax_b)

    # panel letters (outside axes, top-left)
    ax_a.text(0.02, 0.98, "a", transform=ax_a.transAxes,
              fontsize=12, fontweight="bold", va="top", ha="left")
    ax_b.text(0.02, 0.98, "b", transform=ax_b.transAxes,
              fontsize=12, fontweight="bold", va="top", ha="left")

    # shared edge-color legend at bottom
    from matplotlib.lines import Line2D
    legend_handles = [
        Line2D([0], [0], color=RED, linewidth=2.0,
               label="non-synonymous edge"),
        Line2D([0], [0], color=BLUE, linewidth=2.0,
               label="synonymous edge"),
    ]
    fig.legend(
        handles=legend_handles,
        loc="lower center", bbox_to_anchor=(0.5, 0.025),
        ncol=2, fontsize=7, frameon=False,
        handlelength=1.5, handletextpad=0.4,
        columnspacing=1.5,
    )
    # note about substituted base coloring
    fig.text(
        0.5, 0.002,
        "codon labels: substituted base shown in bold edge color",
        ha="center", va="bottom", fontsize=6.5, color="0.4", fontstyle="italic",
    )

    fig.savefig(OUT / "fig4.pdf", dpi=300, bbox_inches="tight")
    fig.savefig(OUT / "fig4.png", dpi=300, bbox_inches="tight")
    print(f"saved: {OUT / 'fig4.pdf'}")
    print(f"saved: {OUT / 'fig4.png'}")
    plt.close(fig)


if __name__ == "__main__":
    main()
