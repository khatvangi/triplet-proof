#!/usr/bin/env python3
"""
build Table S1 — sensitivity of condition C (10-AA triplet) z-score to the
choice of representative amino acid for each RAA10 class.

three strategies:
  - closest to centroid (main text)      — from publication_controls.json
  - farthest from centroid               — from publication_controls.json
  - random (mean ± SD over 3 seeds)      — from publication_controls.json (seed 42)
                                           + sensitivity_extra_seeds.json (seeds 1, 2)

outputs:
  supplementary/table_S1.csv
  supplementary/table_S1.tex
"""

import json, os
import numpy as np
import pandas as pd

ROOT = os.path.abspath(os.path.dirname(__file__))
PUB = os.path.join(ROOT, "results", "publication_controls.json")
EXTRA = os.path.join(ROOT, "results", "sensitivity_extra_seeds.json")
SUPP = os.path.join(ROOT, "supplementary")
os.makedirs(SUPP, exist_ok=True)

with open(PUB) as f:
    pub = json.load(f)
with open(EXTRA) as f:
    extra = json.load(f)

sens = pub["sensitivity"]

# closest (main text)
closest_D = sens["closest"]["D"]["z_score"]
closest_E = sens["closest"]["E"]["z_score"]

# farthest
farthest_D = sens["farthest"]["D"]["z_score"]
farthest_E = sens["farthest"]["E"]["z_score"]

# random — 3 seeds: seed 42 from publication_controls, seeds 1 and 2 from extra
random_D_vals = [
    sens["random"]["D"]["z_score"],          # seed 42
    extra["random_seed_1"]["D_z"],           # seed 1
    extra["random_seed_2"]["D_z"],           # seed 2
]
random_E_vals = [
    sens["random"]["E"]["z_score"],
    extra["random_seed_1"]["E_z"],
    extra["random_seed_2"]["E_z"],
]

random_D_mean = float(np.mean(random_D_vals))
random_D_std = float(np.std(random_D_vals, ddof=1))   # sample SD
random_E_mean = float(np.mean(random_E_vals))
random_E_std = float(np.std(random_E_vals, ddof=1))

# z-score ranges across strategies (using random mean for the random row)
all_D = [closest_D, farthest_D, random_D_mean]
all_E = [closest_E, farthest_E, random_E_mean]
D_range = max(all_D) - min(all_D)
E_range = max(all_E) - min(all_E)

# ── assemble table ──
rows = [
    {
        "Strategy": "Closest to centroid (main text)",
        "Noise z": f"{closest_D:.2f}",
        "Dirichlet z": f"{closest_E:.2f}",
        "Notes": "single deterministic choice per class",
    },
    {
        "Strategy": "Farthest from centroid",
        "Noise z": f"{farthest_D:.2f}",
        "Dirichlet z": f"{farthest_E:.2f}",
        "Notes": "single deterministic choice per class",
    },
    {
        "Strategy": "Random member (mean ± SD, 3 seeds)",
        "Noise z": f"{random_D_mean:.2f} ± {random_D_std:.2f}",
        "Dirichlet z": f"{random_E_mean:.2f} ± {random_E_std:.2f}",
        "Notes": f"seeds 1, 2, 42; individual values: "
                 f"D=[{', '.join(f'{v:.2f}' for v in random_D_vals)}], "
                 f"E=[{', '.join(f'{v:.2f}' for v in random_E_vals)}]",
    },
]

df = pd.DataFrame(rows)

# ── CSV output ──
csv_path = os.path.join(SUPP, "table_S1.csv")
df.to_csv(csv_path, index=False)

# also write summary stats as a separate small CSV for downstream use
summary = {
    "D_range_across_strategies": round(D_range, 3),
    "E_range_across_strategies": round(E_range, 3),
    "random_D_individual_seeds": random_D_vals,
    "random_E_individual_seeds": random_E_vals,
    "random_D_mean": round(random_D_mean, 3),
    "random_D_sd": round(random_D_std, 3),
    "random_E_mean": round(random_E_mean, 3),
    "random_E_sd": round(random_E_std, 3),
    "n_sensitivity_per_run": 100000,
}
with open(os.path.join(SUPP, "table_S1_summary.json"), "w") as f:
    json.dump(summary, f, indent=2)

# ── LaTeX output (booktabs-style) ──
caption = (
    r"Sensitivity of the 10-AA triplet z-score (Condition C) to the choice "
    r"of representative amino acid for each reduced-alphabet (RAA10) group. "
    rf"Noise distortion z-score range across strategies: {D_range:.2f} units. "
    rf"Dirichlet energy z-score range: {E_range:.2f} units. "
    r"Each entry is computed against its own 100{,}000-shuffle null "
    r"distribution. SD is the sample standard deviation across the three "
    r"random seeds."
)

tex = r"""% Table S1 — representative-strategy sensitivity (Condition C, 10-AA triplet)
\begin{table}[htbp]
\centering
\caption{""" + caption + r"""}
\label{tab:S1}
\begin{tabular}{lccp{5cm}}
\toprule
Strategy & Noise $z_D$ & Dirichlet $z_E$ & Notes \\
\midrule
"""

for r in rows:
    strat = r["Strategy"]
    # escape ± for LaTeX
    D = r["Noise z"].replace("±", r"$\pm$")
    E = r["Dirichlet z"].replace("±", r"$\pm$")
    notes = r["Notes"]
    # LaTeX-escape notes
    notes = notes.replace("=[", r"= [").replace(",", ",\\,")
    tex += f"{strat} & {D} & {E} & \\small {notes} \\\\\n"

tex += r"""\bottomrule
\end{tabular}
\end{table}
"""

tex_path = os.path.join(SUPP, "table_S1.tex")
with open(tex_path, "w") as f:
    f.write(tex)

# ── console summary ──
print(f"\nTable S1 (Condition C sensitivity):")
print(df.to_string(index=False))
print()
print(f"D-range across strategies: {D_range:.3f} units  "
      f"({'<' if D_range < 0.5 else '>='} 0.5)")
print(f"E-range across strategies: {E_range:.3f} units  "
      f"({'<' if E_range < 0.5 else '>='} 0.5)")
print()
print(f"saved: {csv_path}")
print(f"saved: {tex_path}")
print(f"saved: {os.path.join(SUPP, 'table_S1_summary.json')}")
