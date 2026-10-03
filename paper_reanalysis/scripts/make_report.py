#!/usr/bin/env python3
"""figures and the results table for the reframed paper.  reads results/grid.csv."""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

HERE = Path(__file__).resolve().parent
PROJ = HERE.parent
OUT = PROJ / "results"
sns.set_theme(style="whitegrid", context="paper")

ORDER = ["pc8_whitened", "pc8_raw", "desc20_z", "prs1d", "blind"]
NICE = {"pc8_whitened": "8 PCs, whitened\n(what the paper computed)",
        "pc8_raw": "8 PCs, unwhitened",
        "desc20_z": "20 descriptors\n(pre-PCA)",
        "prs1d": "polar requirement\n(1-D)",
        "blind": "CHEMISTRY-BLIND\n(all AAs equidistant)"}


def main() -> None:
    t = pd.read_csv(OUT / "grid.csv")
    t["properties"] = pd.Categorical(t.properties, ORDER, ordered=True)
    t = t.sort_values(["objective", "null", "properties"])

    fig, ax = plt.subplots(1, 2, figsize=(13, 5), sharey=True)
    for k, objective in enumerate(("E", "D")):
        s = t[t.objective == objective].copy()
        s["label"] = s.properties.map(NICE)
        sns.barplot(data=s, x="label", y="z", hue="null", ax=ax[k],
                    palette={"codon": "#c0392b", "block": "#2980b9"})
        ax[k].axhline(0, color="black", lw=0.8)
        ax[k].axhline(-1.96, color="grey", ls="--", lw=0.8)
        ax[k].set_title(f"objective {objective}", fontsize=11)
        ax[k].set_xlabel("")
        ax[k].set_ylabel("z (negative = SGC better than null)" if k == 0 else "")
        ax[k].tick_params(axis="x", rotation=30, labelsize=8)
        for c in ax[k].containers:
            ax[k].bar_label(c, fmt="%.1f", fontsize=7, padding=1)
        # an invariant cell has no bar at all (z is NaN).  say so, or a reader
        # sees a gap and assumes the run failed.
        if (s.zero_variance & (s.null == "block")).any():
            ax[k].annotate("exactly invariant\n(p = 1)", xy=(4.2, 0), xytext=(4.2, -0.9),
                           ha="center", fontsize=7, color="#2980b9")
    ax[0].legend(title="null model", loc="lower left")
    if ax[1].get_legend():
        ax[1].get_legend().remove()
    fig.suptitle("Reported optimisation of the standard genetic code depends on "
                 "the null model, not on chemistry", fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    fig.savefig(OUT / "fig_null_vs_chemistry.png", dpi=220)
    fig.savefig(OUT / "fig_null_vs_chemistry.pdf")
    plt.close(fig)

    # the headline contrast, as a compact markdown table
    piv = t.pivot_table(index="properties", columns=["objective", "null"],
                        values="z", observed=True).round(2)
    piv.index = [NICE[i].replace("\n", " ") for i in piv.index]
    (OUT / "TABLE_z_grid.md").write_text(piv.to_markdown() + "\n")
    print(piv.to_string())
    print(f"\nwrote {OUT/'fig_null_vs_chemistry.png'}")


if __name__ == "__main__":
    main()
