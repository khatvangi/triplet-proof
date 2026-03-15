"""
Fig 2 — Wobble mechanism schematic for JME Letter.

Panel A: Triplet code — Hamming-1 neighborhood of GCU (Ala).
  3/9 neighbors are synonymous, all via position-3 (wobble) changes.

Panel B: Doublet code — Hamming-1 neighborhood of GC (Ala/Gly).
  Only 1/6 neighbors is synonymous. No dedicated synonymy position.

Visual message: the third codon position creates an error-absorbing
corridor that doublet codes lack.
"""

import numpy as np
import matplotlib.pyplot as plt

# --- palette (Paul Tol colorblind-safe) ---
BLUE = "#4477AA"
RED = "#CC6677"
BLUE_FILL = (*plt.matplotlib.colors.to_rgb(BLUE), 0.15)
RED_FILL = (*plt.matplotlib.colors.to_rgb(RED), 0.15)

# --- global font settings ---
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.size": 8,
    "axes.titlesize": 9,
    "axes.titleweight": "bold",
})

LABEL_SIZE = 7.5
ANNOT_SIZE = 7
NODE_FONT = 5.5


def draw_node(ax, xy, r, label, fill="white", edgecolor="black", lw=1.0,
              fontsize=NODE_FONT, bold=False):
    """draw a circle node with a centered label."""
    circle = plt.Circle(xy, r, facecolor=fill, edgecolor=edgecolor,
                        linewidth=lw, zorder=3, clip_on=False)
    ax.add_patch(circle)
    weight = "bold" if bold else "normal"
    ax.text(xy[0], xy[1], label, ha="center", va="center",
            fontsize=fontsize, fontweight=weight, zorder=4)


def draw_edge(ax, xy_from, xy_to, r_from, r_to, color, lw=1.0):
    """draw a line between two nodes, clipped to circle boundaries."""
    dx = xy_to[0] - xy_from[0]
    dy = xy_to[1] - xy_from[1]
    dist = np.hypot(dx, dy)
    if dist == 0:
        return
    ux, uy = dx / dist, dy / dist
    x0 = xy_from[0] + ux * r_from
    y0 = xy_from[1] + uy * r_from
    x1 = xy_to[0] - ux * r_to
    y1 = xy_to[1] - uy * r_to
    ax.plot([x0, x1], [y0, y1], color=color, linewidth=lw, zorder=1,
            solid_capstyle="round")


def polar_xy(cx, cy, angle_deg, radius):
    """convert polar to cartesian."""
    a = np.radians(angle_deg)
    return (cx + radius * np.cos(a), cy + radius * np.sin(a))


# =========================================================================
# panel A — triplet: GCU neighborhood (9 neighbors)
# =========================================================================

def draw_panel_a(ax):
    cx, cy = 0.0, 0.0
    R_center = 0.45
    R_outer = 0.35
    orbit = 1.85

    # 3 groups of 3, centered at 90° (top), 215° (lower-left), 325° (lower-right).
    # shift lower-right group slightly right so "pos 3 (wobble)" label fits.
    spread = 30

    # pos 1 — top sector
    pos1 = [("ACU", "Thr"), ("UCU", "Ser"), ("CCU", "Pro")]
    a1_center = 90
    angles1 = [a1_center + spread, a1_center, a1_center - spread]

    # pos 2 — lower-left sector
    pos2 = [("GGU", "Gly"), ("GUU", "Val"), ("GAU", "Asp")]
    a2_center = 215
    angles2 = [a2_center + spread, a2_center, a2_center - spread]

    # pos 3 — lower-right sector (synonymous)
    pos3 = [("GCG", "Ala"), ("GCA", "Ala"), ("GCC", "Ala")]
    a3_center = 325
    angles3 = [a3_center + spread, a3_center, a3_center - spread]

    # draw edges + nodes: pos 1
    for (codon, aa), ang in zip(pos1, angles1):
        xy = polar_xy(cx, cy, ang, orbit)
        draw_edge(ax, (cx, cy), xy, R_center, R_outer, RED)
        draw_node(ax, xy, R_outer, f"{codon}\n{aa}", fill=RED_FILL)

    # pos 2
    for (codon, aa), ang in zip(pos2, angles2):
        xy = polar_xy(cx, cy, ang, orbit)
        draw_edge(ax, (cx, cy), xy, R_center, R_outer, RED)
        draw_node(ax, xy, R_outer, f"{codon}\n{aa}", fill=RED_FILL)

    # pos 3 (synonymous)
    for (codon, aa), ang in zip(pos3, angles3):
        xy = polar_xy(cx, cy, ang, orbit)
        draw_edge(ax, (cx, cy), xy, R_center, R_outer, BLUE)
        draw_node(ax, xy, R_outer, f"{codon}\n{aa}", fill=BLUE_FILL)

    # center node on top
    draw_node(ax, (cx, cy), R_center, "GCU\nAla", fill="white",
              edgecolor="black", lw=1.5, fontsize=LABEL_SIZE, bold=True)

    # group labels — outside the outer ring, positioned carefully
    # pos 1 label: above the top group
    lx1, ly1 = polar_xy(cx, cy, a1_center, orbit + R_outer + 0.15)
    ax.text(lx1, ly1, "pos 1",
            fontsize=ANNOT_SIZE, color=RED, ha="center", va="bottom",
            style="italic")

    # pos 2 label: to the left of lower-left group
    lx2, ly2 = polar_xy(cx, cy, a2_center, orbit + R_outer + 0.15)
    ax.text(lx2 - 0.15, ly2, "pos 2",
            fontsize=ANNOT_SIZE, color=RED, ha="right", va="center",
            style="italic")

    # pos 3 label: to the right of lower-right group
    lx3, ly3 = polar_xy(cx, cy, a3_center, orbit + R_outer + 0.15)
    ax.text(lx3 + 0.15, ly3, "pos 3\n(wobble)",
            fontsize=ANNOT_SIZE, color=BLUE, ha="left", va="center",
            style="italic", linespacing=1.0)

    # synonymy summary
    ax.text(0.0, -(orbit + R_outer + 0.55),
            "synonymous: 3/9 (all pos 3)",
            fontsize=ANNOT_SIZE, ha="center", va="top", color=BLUE,
            fontweight="bold")

    lim = orbit + R_outer + 0.90
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-(orbit + R_outer + 0.75), lim - 0.20)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title("Triplet: GCU neighborhood", pad=6)

    # panel label
    ax.text(-lim + 0.05, lim - 0.25, "a", fontsize=11, fontweight="bold",
            ha="left", va="top")


# =========================================================================
# panel B — doublet: GC neighborhood (6 neighbors)
# =========================================================================

def draw_panel_b(ax):
    cx, cy = 0.0, 0.0
    R_center = 0.45
    R_outer = 0.38
    orbit = 1.75

    # 2 groups of 3 at 90° (top) and 270° (bottom), spread ±45°
    spread = 45

    # pos 1 — top
    pos1 = [("AC", "Ser/Thr"), ("UC", "Ser/Thr"), ("CC", "Pro")]
    a1_center = 90
    angles1 = [a1_center + spread, a1_center, a1_center - spread]

    # pos 2 — bottom (GG is synonymous)
    pos2_data = [
        ("GA", "Asp/Glu", False),
        ("GU", "VLIM", False),
        ("GG", "Ala/Gly", True),
    ]
    a2_center = 270
    angles2 = [a2_center + spread, a2_center, a2_center - spread]

    # pos 1 nodes
    for (codon, aa), ang in zip(pos1, angles1):
        xy = polar_xy(cx, cy, ang, orbit)
        draw_edge(ax, (cx, cy), xy, R_center, R_outer, RED)
        draw_node(ax, xy, R_outer, f"{codon}\n{aa}", fill=RED_FILL)

    # pos 2 nodes
    for (codon, aa, syn), ang in zip(pos2_data, angles2):
        xy = polar_xy(cx, cy, ang, orbit)
        ec = BLUE if syn else RED
        fc = BLUE_FILL if syn else RED_FILL
        draw_edge(ax, (cx, cy), xy, R_center, R_outer, ec)
        draw_node(ax, xy, R_outer, f"{codon}\n{aa}", fill=fc)

    # center node
    draw_node(ax, (cx, cy), R_center, "GC\nAla/Gly", fill="white",
              edgecolor="black", lw=1.5, fontsize=LABEL_SIZE, bold=True)

    # group labels
    lx1, ly1 = polar_xy(cx, cy, a1_center, orbit + R_outer + 0.15)
    ax.text(lx1, ly1, "pos 1",
            fontsize=ANNOT_SIZE, color=RED, ha="center", va="bottom",
            style="italic")

    lx2, ly2 = polar_xy(cx, cy, a2_center, orbit + R_outer + 0.15)
    ax.text(lx2, ly2, "pos 2",
            fontsize=ANNOT_SIZE, color=RED, ha="center", va="top",
            style="italic")

    # synonymy summary
    ax.text(0.0, -(orbit + R_outer + 0.55),
            "synonymous: 1/6 (no dedicated position)",
            fontsize=ANNOT_SIZE, ha="center", va="top", color=BLUE,
            fontweight="bold")

    lim = orbit + R_outer + 0.80
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-(orbit + R_outer + 0.75), lim - 0.20)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title("Doublet: GC neighborhood", pad=6)

    # panel label
    ax.text(-lim + 0.05, lim - 0.25, "b", fontsize=11, fontweight="bold",
            ha="left", va="top")


# =========================================================================
# compose figure
# =========================================================================

def main():
    fig, (ax_a, ax_b) = plt.subplots(2, 1, figsize=(3.5, 4.5))
    fig.subplots_adjust(hspace=0.30, top=0.95, bottom=0.02, left=0.02, right=0.98)

    draw_panel_a(ax_a)
    draw_panel_b(ax_b)

    outbase = "figures/fig2_wobble_mechanism"
    fig.savefig(f"{outbase}.pdf", dpi=300, bbox_inches="tight")
    fig.savefig(f"{outbase}.png", dpi=300, bbox_inches="tight")
    print(f"saved {outbase}.pdf and {outbase}.png")
    plt.close(fig)


if __name__ == "__main__":
    main()
