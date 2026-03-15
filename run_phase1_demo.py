# run_phase1_demo.py (root)
# phase-1 sanity check: small Monte Carlo on triplet SGC
import yaml, json
import numpy as np
import pandas as pd

from src.io.codon_io import load_codon_table_csv
from src.io.aa_props_io import load_aa_props, aa_prop_vector_map, prop_names
from src.sims.codon_graph import build_graph_n, label_nodes_with_aa
from src.metrics.dirichlet import dirichlet_energy
from src.metrics.mi import mutual_info_by_position
from src.sims.random_codes import fully_random_code, degeneracy_preserving_shuffle
from src.receiver.thermo_noise import expected_property_distortion_noise_multi as noise_multi

def main():
    cfg = yaml.safe_load(open("configs/phase1.yaml"))
    codon_df = load_codon_table_csv(cfg["data"]["codon_table_csv"])
    # load_aa_props returns (df, meta) tuple
    aa_df, _ = load_aa_props(cfg["data"]["aa_props_csv"])
    feat_cols = prop_names(aa_df)
    aa_map = aa_prop_vector_map(aa_df, feature_cols=feat_cols)

    # --- triplet graph + Dirichlet (phase-1)
    G = build_graph_n(3)
    sgc_map = {row["codon"]: (None if row["is_stop"] else row.get("aa")) for _, row in codon_df.iterrows()}
    label_nodes_with_aa(G, sgc_map, exclude_stops=cfg["options"]["exclude_stops_in_metrics"])

    prop_map = {node: (aa_map.get(G.nodes[node]["aa"]) if G.nodes[node]["aa"] else None) for node in G.nodes()}
    E_sgc = dirichlet_energy(G, prop_map)
    print(f"Dirichlet energy (SGC): {E_sgc:.4f}")

    mi = mutual_info_by_position(codon_df, aa_map, feat_cols)
    print("Mutual information by position (mean per property):")
    for pname, vals in mi.items():
        print(f"  {pname}: pos1={vals[0]:.4f}, pos2={vals[1]:.4f}, pos3={vals[2]:.4f}")

    # --- receiver/noise distortion (multi-property)
    N_sgc = noise_multi(codon_df, aa_map, RT=1.0, wobble_rules='strict', asymmetry=True)
    print(f"\nReceiver/noise distortion N (SGC, multi-prop): {N_sgc:.4f}")

    # --- nulls (demo sizes; scale later)
    nrand = int(cfg["options"]["n_random_codes_demo"])
    seed0 = int(cfg["options"]["random_seed"])

    energies_deg  = []
    noises_deg    = []

    for i in range(nrand):
        rnd_map2 = degeneracy_preserving_shuffle(codon_df, seed=seed0 + i)

        # label and Dirichlet
        label_nodes_with_aa(G, rnd_map2, exclude_stops=True)
        prop_map_r2 = {node: (aa_map.get(G.nodes[node]["aa"]) if G.nodes[node]["aa"] else None) for node in G.nodes()}
        energies_deg.append(dirichlet_energy(G, prop_map_r2))

        # noise: mutate codon_df's 'aa' temporarily
        tmp = codon_df.copy()
        tmp["aa"] = tmp["codon"].map(lambda c: rnd_map2.get(c, None))
        noises_deg.append(noise_multi(tmp, aa_map, RT=1.0, wobble_rules='strict', asymmetry=True))

    # restore SGC labels
    label_nodes_with_aa(G, sgc_map, exclude_stops=True)

    # percentile with Laplace smoothing (consistent with phase2 scripts)
    p_deg_E   = (np.sum(np.array(energies_deg) <= E_sgc) + 1) / (nrand + 1)
    p_deg_No  = (np.sum(np.array(noises_deg)  <= N_sgc) + 1) / (nrand + 1)

    print(f"\nSGC percentile (Dirichlet) vs degeneracy-preserving: {p_deg_E*100:.2f}%")
    print(f"SGC percentile (Noise) vs degeneracy-preserving: {p_deg_No*100:.2f}%")

    out = {
        "E_sgc": float(E_sgc),
        "N_sgc": float(N_sgc),
        "percentile_deg_E": float(p_deg_E),
        "percentile_deg_Noise": float(p_deg_No),
        "n_random": nrand
    }
    with open("results/phase1_demo.json","w") as f:
        json.dump(out, f, indent=2)
    print("\nSaved results → results/phase1_demo.json")

if __name__ == "__main__":
    main()

