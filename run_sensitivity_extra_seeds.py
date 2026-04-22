#!/usr/bin/env python3
"""
run 2 additional random-seed sensitivity runs for Table S1.

existing publication_controls.json has 1 'random' seed (42).
this script runs seeds 1 and 2 with identical setup and saves results.

writes: results/sensitivity_extra_seeds.json
"""
import os, sys, json, random
from datetime import datetime

sys.path.insert(0, os.path.abspath(os.path.dirname(__file__)))

import numpy as np
from src.io.codon_io import load_codon_table_csv
from src.io.aa_props_io import load_aa_props, aa_prop_vector_map
from src.sims.codon_graph import build_graph_n
from run_publication_controls import (
    RAA10, select_representatives, build_triplet_10aa, analyze_condition,
)

N_SENS = 100_000
WORKERS = 32

print(f"[{datetime.now()}] extra-seeds sensitivity run")

triplet_df = load_codon_table_csv()
aa_df, _ = load_aa_props("data/processed/aa_props.parquet")
vec_map_20 = aa_prop_vector_map(aa_df)

G3 = build_graph_n(3)
edges_n3 = list(G3.edges())

triplet_10aa_df = build_triplet_10aa(triplet_df, RAA10)

results = {"n_sensitivity": N_SENS, "seeds_run": []}

for seed in [1, 2]:
    print(f"\n=== random seed {seed} ===")
    rep_map, rep_info = select_representatives(vec_map_20, RAA10, "random", seed=seed)
    chosen = {cls: rep_info[cls]["chosen"] for cls in sorted(rep_info)}
    print(f"  chosen: {chosen}")

    result, _, _ = analyze_condition(
        f"C random seed={seed}", triplet_10aa_df, rep_map, edges_n3, N_SENS, WORKERS
    )
    results[f"random_seed_{seed}"] = {
        "chosen": chosen,
        "D_z": result["D"]["z_score"],
        "E_z": result["E"]["z_score"],
        "D_sgc": result["D"]["sgc"],
        "E_sgc": result["E"]["sgc"],
        "D_null_mean": result["D"]["null_mean"],
        "E_null_mean": result["E"]["null_mean"],
        "DE_corr": result["DE_correlation"],
    }
    results["seeds_run"].append(seed)

out = "results/sensitivity_extra_seeds.json"
with open(out, "w") as f:
    json.dump(results, f, indent=2)
print(f"\nsaved: {out}")
print(f"[{datetime.now()}] done")
