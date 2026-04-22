#!/usr/bin/env python3
"""
run_publication_controls.py — all pre-manuscript verification in one run

three tasks:
1. 2×2 factorial (10/20 AA × n=2/n=3) at 1M samples, BOTH metrics (D + E)
2. Dirichlet energy concordance check (part of task 1)
3. sensitivity: 10AA@n=3 with alternative representative sets at 100K

computes both raw distortion (D) and Dirichlet energy (E) per shuffle
to avoid duplicating work and to check metric concordance.

Factorial structure — why n=4 is NOT a cell:
    The factorial is {10 AA, 20 AA} × {n=2, n=3}. The 20-AA × n=2 cell is
    structurally impossible (4^2 = 16 codons, minus 3 stops = 13 sense
    codons, which cannot encode 20 amino acids), so the factorial is
    triangular: conditions A (10-AA doublet), C (10-AA triplet), B (20-AA
    SGC). n=4 is excluded for three deliberate reasons:

      1. No biological SGC at n=4 — all extant life uses triplet codons.
         The doublet comparator (condition A) is a principled
         majority-vote projection of the real SGC onto its first two
         positions; quadruplet would require inventing a 4-position code
         from scratch, with no principled reduction from biology.
      2. The specific n=4 construction implemented in run_phase2_n2n4.py
         (`AA(b1 b2 b3 b4) := AA_SGC(b1 b2 b3)`, duplicated along b4)
         makes position 4 100 % synonymous by design: every single-base
         substitution at b4 is silent. This gifts the constructed SGC a
         fully synonymous position that random codes on the 4-mer graph
         do not automatically have. The resulting z-score (0/10^6 beyond
         null, see results/phase2_quadruplet.auto.json) is dominated by
         this construction artefact, not by any biological property.
      3. The factorial's question is whether triplet architecture is
         optimized relative to a SIMPLER (more constrained) architecture.
         n=2 is the right comparator; n=4 would test whether a longer
         code could be even better if it existed, which is a different
         question and not one the SGC-optimality claim needs to answer.

    n=4 machinery is kept in run_phase2_n2n4.py for exploratory work and
    as a sanity check that the SGC remains extreme when embedded in a
    4-mer graph — but its numbers should not be read as a factorial cell.
"""

import os, sys, json, random
from datetime import datetime
from collections import Counter, defaultdict
from itertools import product

sys.path.insert(0, os.path.abspath(os.path.dirname(__file__)))

import numpy as np
import pandas as pd
import multiprocessing as mp

from src.io.codon_io import load_codon_table_csv
from src.io.aa_props_io import load_aa_props, aa_prop_vector_map, prop_names
from src.sims.codon_graph import build_graph_n

RNA_BASES = ("A", "U", "G", "C")

RAA10 = {
    "A": "X1", "G": "X1",
    "V": "X2", "L": "X2", "I": "X2", "M": "X2",
    "F": "X3", "W": "X3", "Y": "X3",
    "S": "X4", "T": "X4",
    "N": "X5", "Q": "X5",
    "D": "X6", "E": "X6",
    "K": "X7", "R": "X7",
    "H": "X8",
    "P": "X9",
    "C": "X10",
}


# ── representative selection methods ──

def select_representatives(aa_vec_map, raa_map, method="closest", seed=42):
    """select one AA per RAA class using the given method.
    returns (vec_map, selection_info) where vec_map[cls] = vector."""
    class_members = defaultdict(dict)
    for aa, cls in raa_map.items():
        if aa.upper() in aa_vec_map:
            class_members[cls][aa.upper()] = aa_vec_map[aa.upper()]

    vec_map = {}
    info = {}
    rng = random.Random(seed)

    for cls, members in class_members.items():
        centroid = np.mean(list(members.values()), axis=0)
        dists = {aa: float(np.linalg.norm(v - centroid)) for aa, v in members.items()}

        if method == "closest":
            chosen = min(dists, key=dists.get)
        elif method == "farthest":
            chosen = max(dists, key=dists.get)
        elif method == "random":
            chosen = rng.choice(sorted(members.keys()))
        else:
            raise ValueError(f"unknown method: {method}")

        vec_map[cls] = members[chosen]
        info[cls] = {"members": sorted(members.keys()), "chosen": chosen,
                     "dist_to_centroid": dists[chosen]}

    return vec_map, info


# ── codon table builders ──

def majority_vote(labels):
    cnt = Counter(labels)
    maxc = max(cnt.values())
    return sorted(lab for lab, c in cnt.items() if c == maxc)[0]


def build_doublet_table(triplet_df, raa_map):
    """project triplet SGC → doublet via majority vote on pos 3."""
    t2c = {}
    for _, row in triplet_df.iterrows():
        aa = str(row["aa"]).strip().upper()
        codon = str(row["codon"]).strip().upper()
        is_stop = str(row.get("is_stop", "")).strip().upper() in ("TRUE", "1", "YES", "*")
        t2c[codon] = "*" if (is_stop or aa in ("*", "STOP")) else raa_map.get(aa, aa)
    rows = []
    for b1 in RNA_BASES:
        for b2 in RNA_BASES:
            pf = b1 + b2
            classes = [t2c.get(pf + b3, "*") for b3 in RNA_BASES]
            non_stop = [c for c in classes if c != "*"]
            cls = majority_vote(non_stop) if non_stop else "*"
            rows.append({"codon": pf, "aa": cls, "is_stop": cls == "*"})
    return pd.DataFrame(rows)


def build_triplet_10aa(triplet_df, raa_map):
    """project SGC 20 AAs → 10 RAA classes on same 64 codons."""
    rows = []
    for _, row in triplet_df.iterrows():
        aa = str(row["aa"]).strip().upper()
        codon = str(row["codon"]).strip().upper()
        is_stop = str(row.get("is_stop", "")).strip().upper() in ("TRUE", "1", "YES", "*")
        if is_stop or aa in ("*", "STOP"):
            rows.append({"codon": codon, "aa": "*", "is_stop": True})
        else:
            rows.append({"codon": codon, "aa": raa_map.get(aa, aa), "is_stop": False})
    return pd.DataFrame(rows)


# ── dual-metric computation ──

def compute_sgc_metrics(codon_df, aa_vec_map, edges):
    """compute both D (raw distortion) and E (Dirichlet energy) for a codon table."""
    sense = codon_df[~codon_df["is_stop"].astype(bool)]
    c2a = dict(zip(sense["codon"].str.upper(), sense["aa"].str.upper()))

    # D: raw distortion per codon
    total_d = 0.0
    n_sense = 0
    for codon, aa in c2a.items():
        if aa not in aa_vec_map:
            continue
        vec_c = aa_vec_map[aa]
        n_sense += 1
        for pos in range(len(codon)):
            for base in RNA_BASES:
                if base == codon[pos]:
                    continue
                mut = codon[:pos] + base + codon[pos + 1:]
                if mut in c2a:
                    mut_aa = c2a[mut]
                    if mut_aa != aa and mut_aa in aa_vec_map:
                        total_d += float(np.linalg.norm(vec_c - aa_vec_map[mut_aa]))
    D = total_d / n_sense if n_sense > 0 else 0

    # E: Dirichlet energy on the Hamming graph
    prop_map = {c: aa_vec_map[aa] for c, aa in c2a.items() if aa in aa_vec_map}
    E = 0.0
    for u, v in edges:
        pu = prop_map.get(u)
        pv = prop_map.get(v)
        if pu is not None and pv is not None:
            duv = pu - pv
            E += float(np.dot(duv, duv))

    return D, E


def dual_null_worker(args):
    """compute both D and E for n_iter shuffles. returns paired (D, E) lists."""
    seed, n_iter, codons, labels_template, aa_vec_map, edges, codon_len = args
    rng = random.Random(seed)

    D_vals = []
    E_vals = []

    for _ in range(n_iter):
        # shuffle AA labels (preserves degeneracy)
        labels = labels_template[:]
        rng.shuffle(labels)
        c2a = dict(zip(codons, labels))

        # D: raw distortion
        total_d = 0.0
        n_s = 0
        for codon, aa in c2a.items():
            if aa not in aa_vec_map:
                continue
            vec_c = aa_vec_map[aa]
            n_s += 1
            for pos in range(codon_len):
                for base in RNA_BASES:
                    if base == codon[pos]:
                        continue
                    mut = codon[:pos] + base + codon[pos + 1:]
                    if mut in c2a:
                        mut_aa = c2a[mut]
                        if mut_aa != aa and mut_aa in aa_vec_map:
                            total_d += float(np.linalg.norm(vec_c - aa_vec_map[mut_aa]))
        D_vals.append(total_d / n_s if n_s > 0 else 0)

        # E: Dirichlet energy
        prop_map = {c: aa_vec_map[aa] for c, aa in c2a.items() if aa in aa_vec_map}
        E = 0.0
        for u, v in edges:
            pu = prop_map.get(u)
            pv = prop_map.get(v)
            if pu is not None and pv is not None:
                duv = pu - pv
                E += float(np.dot(duv, duv))
        E_vals.append(E)

    return D_vals, E_vals


def run_dual_null(codons, labels, aa_vec_map, edges, codon_len, n_null, workers):
    """run null shuffles computing both D and E. returns paired lists."""
    chunk = n_null // workers
    seeds = [random.randint(0, 2**31) for _ in range(workers)]
    args = [(s, chunk, codons, labels, aa_vec_map, edges, codon_len) for s in seeds]

    all_D = []
    all_E = []
    with mp.Pool(workers) as pool:
        for D_vals, E_vals in pool.imap_unordered(dual_null_worker, args):
            all_D.extend(D_vals)
            all_E.extend(E_vals)
    return all_D, all_E


def stats_from_null(sgc_val, null_vals):
    """compute z-score, count, percentile from null distribution."""
    null_mean = np.mean(null_vals)
    null_std = np.std(null_vals)
    z = (sgc_val - null_mean) / null_std if null_std > 0 else 0
    count_le = sum(1 for v in null_vals if v <= sgc_val)
    pct = (count_le + 1) / (len(null_vals) + 1)
    return {
        "sgc": sgc_val, "null_mean": null_mean, "null_std": null_std,
        "null_min": min(null_vals), "null_max": max(null_vals),
        "null_cv": null_std / null_mean if null_mean > 0 else 0,
        "z_score": z, "count_le": count_le,
        "n_null": len(null_vals), "percentile": pct,
    }


def analyze_condition(label, codon_df, aa_vec_map, edges, n_null, workers):
    """run full dual-metric analysis for one condition.
    returns (summary_dict, raw_D_array, raw_E_array)."""
    sense = codon_df[~codon_df["is_stop"].astype(bool)]
    codons = sense["codon"].str.upper().tolist()
    labels = sense["aa"].str.upper().tolist()
    n_bases = len(codons[0])
    n_classes = len(set(labels))

    print(f"\n{'=' * 65}")
    print(f"  {label}")
    print(f"  {len(codons)} sense codons, {n_classes} classes, n={n_bases}")
    print(f"{'=' * 65}")

    # SGC metrics
    D_sgc, E_sgc = compute_sgc_metrics(codon_df, aa_vec_map, edges)
    print(f"  SGC: D={D_sgc:.4f}, E={E_sgc:.4f}")

    # null distribution (both metrics per shuffle)
    print(f"  running {n_null:,} shuffles (dual metric) on {workers} workers...")
    t0 = datetime.now()
    all_D, all_E = run_dual_null(codons, labels, aa_vec_map, edges, n_bases, n_null, workers)
    elapsed = (datetime.now() - t0).total_seconds()
    print(f"  completed in {elapsed:.1f}s")

    d_stats = stats_from_null(D_sgc, all_D)
    e_stats = stats_from_null(E_sgc, all_E)

    # D-E correlation within the null
    corr = float(np.corrcoef(all_D, all_E)[0, 1])

    for metric, s in [("D (distortion)", d_stats), ("E (Dirichlet)", e_stats)]:
        print(f"  {metric}: z={s['z_score']:.2f}, null_mean={s['null_mean']:.4f}, "
              f"null_std={s['null_std']:.4f}, count≤={s['count_le']}/{s['n_null']}")
    print(f"  D-E correlation in null: r={corr:.4f}")

    summary = {
        "n_bases": n_bases, "n_sense": len(codons), "n_classes": n_classes,
        "D": d_stats, "E": e_stats, "DE_correlation": corr,
    }
    return summary, np.array(all_D, dtype=np.float64), np.array(all_E, dtype=np.float64)


def main():
    import argparse
    ap = argparse.ArgumentParser(description="Publication controls: factorial + Dirichlet + sensitivity")
    ap.add_argument("--n-null", type=int, default=1000000, help="samples for main factorial")
    ap.add_argument("--n-sensitivity", type=int, default=100000, help="samples for sensitivity")
    ap.add_argument("--workers", type=int, default=32)
    args = ap.parse_args()

    print(f"[{datetime.now()}] publication controls")
    print(f"  factorial: {args.n_null:,} samples | sensitivity: {args.n_sensitivity:,} samples")

    # ── load data ──
    triplet_df = load_codon_table_csv()
    parquet = os.path.join(os.path.dirname(__file__), "data", "processed", "aa_props.parquet")
    aa_df, _ = load_aa_props(parquet)
    vec_map_20 = aa_prop_vector_map(aa_df)
    print(f"  features: {prop_names(aa_df)}")

    # build Hamming graphs and extract edge lists
    G2 = build_graph_n(2)
    G3 = build_graph_n(3)
    edges_n2 = list(G2.edges())
    edges_n3 = list(G3.edges())
    print(f"  n=2 graph: {G2.number_of_nodes()} nodes, {G2.number_of_edges()} edges")
    print(f"  n=3 graph: {G3.number_of_nodes()} nodes, {G3.number_of_edges()} edges")

    # build 10-AA representative vectors (closest to centroid)
    rep_closest, info_closest = select_representatives(vec_map_20, RAA10, "closest")
    print(f"\n  closest-to-centroid representatives:")
    for cls in sorted(info_closest.keys()):
        i = info_closest[cls]
        print(f"    {cls}: {i['members']} → {i['chosen']}")

    # build codon tables
    doublet_df = build_doublet_table(triplet_df, RAA10)
    triplet_10aa_df = build_triplet_10aa(triplet_df, RAA10)

    results = {"timestamp": str(datetime.now()), "n_null": args.n_null,
               "n_sensitivity": args.n_sensitivity}

    # ════════════════════════════════════════════════════
    # PART 1+2: 2×2 factorial with dual metrics
    # ════════════════════════════════════════════════════
    print(f"\n{'#' * 70}")
    print(f"  PART 1+2: 2×2 FACTORIAL (D + E), {args.n_null:,} samples")
    print(f"{'#' * 70}")

    # raw null distributions stored here for figure generation
    raw_dists = {}

    results["A_10AA_n2"], raw_dists["A_D"], raw_dists["A_E"] = analyze_condition(
        "A: 10 AAs × n=2 (doublet)", doublet_df, rep_closest, edges_n2,
        args.n_null, args.workers)

    results["B_20AA_n3"], raw_dists["B_D"], raw_dists["B_E"] = analyze_condition(
        "B: 20 AAs × n=3 (SGC)", triplet_df, vec_map_20, edges_n3,
        args.n_null, args.workers)

    results["C_10AA_n3"], raw_dists["C_D"], raw_dists["C_E"] = analyze_condition(
        "C: 10 AAs × n=3 (CONTROL)", triplet_10aa_df, rep_closest, edges_n3,
        args.n_null, args.workers)

    # ── factorial summary ──
    print(f"\n{'=' * 70}")
    print("  FACTORIAL SUMMARY")
    print(f"{'=' * 70}")

    zA_D = results["A_10AA_n2"]["D"]["z_score"]
    zB_D = results["B_20AA_n3"]["D"]["z_score"]
    zC_D = results["C_10AA_n3"]["D"]["z_score"]
    zA_E = results["A_10AA_n2"]["E"]["z_score"]
    zB_E = results["B_20AA_n3"]["E"]["z_score"]
    zC_E = results["C_10AA_n3"]["E"]["z_score"]

    print(f"\n  raw distortion (D):")
    print(f"  {'':>15} {'n=2':>10} {'n=3':>10} {'Δ arch':>10}")
    print(f"  {'10 AAs':<15} {zA_D:>10.2f} {zC_D:>10.2f} {zC_D - zA_D:>+10.2f}")
    print(f"  {'20 AAs':<15} {'—':>10} {zB_D:>10.2f}")
    print(f"  {'Δ alphabet':<15} {'':>10} {zB_D - zC_D:>+10.2f}")
    arch_D = abs(zC_D - zA_D)
    alpha_D = abs(zB_D - zC_D)
    ratio_D = arch_D / alpha_D if alpha_D > 0 else float('inf')

    print(f"\n  Dirichlet energy (E):")
    print(f"  {'':>15} {'n=2':>10} {'n=3':>10} {'Δ arch':>10}")
    print(f"  {'10 AAs':<15} {zA_E:>10.2f} {zC_E:>10.2f} {zC_E - zA_E:>+10.2f}")
    print(f"  {'20 AAs':<15} {'—':>10} {zB_E:>10.2f}")
    print(f"  {'Δ alphabet':<15} {'':>10} {zB_E - zC_E:>+10.2f}")
    arch_E = abs(zC_E - zA_E)
    alpha_E = abs(zB_E - zC_E)
    ratio_E = arch_E / alpha_E if alpha_E > 0 else float('inf')

    print(f"\n  architecture:alphabet ratio — D: {ratio_D:.1f}:1, E: {ratio_E:.1f}:1")

    # concordance check
    print(f"\n  D-E concordance:")
    for key, label in [("A_10AA_n2", "A"), ("B_20AA_n3", "B"), ("C_10AA_n3", "C")]:
        r = results[key]
        print(f"    {label}: D z={r['D']['z_score']:.2f}, E z={r['E']['z_score']:.2f}, "
              f"corr={r['DE_correlation']:.3f}")

    if abs(ratio_D - ratio_E) < 2:
        print(f"\n  ✓ both metrics agree: architecture dominates")
    else:
        print(f"\n  ⚠ metrics diverge — investigate before publishing")

    results["factorial_summary"] = {
        "arch_effect_D": zC_D - zA_D, "alpha_effect_D": zB_D - zC_D, "ratio_D": ratio_D,
        "arch_effect_E": zC_E - zA_E, "alpha_effect_E": zB_E - zC_E, "ratio_E": ratio_E,
    }

    # ════════════════════════════════════════════════════
    # PART 3: SENSITIVITY (alternative representative sets)
    # ════════════════════════════════════════════════════
    print(f"\n{'#' * 70}")
    print(f"  PART 3: SENSITIVITY (condition C with alternative reps), {args.n_sensitivity:,}")
    print(f"{'#' * 70}")

    methods = ["closest", "farthest", "random"]
    sensitivity = {}

    for method in methods:
        rep_map, rep_info = select_representatives(vec_map_20, RAA10, method)
        chosen = {cls: rep_info[cls]["chosen"] for cls in sorted(rep_info)}
        print(f"\n  representative set '{method}': {chosen}")

        result, sens_D, sens_E = analyze_condition(
            f"C ({method} reps): 10 AAs × n=3",
            triplet_10aa_df, rep_map, edges_n3, args.n_sensitivity, args.workers)
        sensitivity[method] = result
        # save sensitivity distributions too (smaller, 100K)
        raw_dists[f"sens_{method}_D"] = sens_D
        raw_dists[f"sens_{method}_E"] = sens_E

    results["sensitivity"] = sensitivity

    # ── sensitivity summary ──
    print(f"\n{'=' * 70}")
    print("  SENSITIVITY SUMMARY: condition C (10AA@n=3) across rep sets")
    print(f"{'=' * 70}")
    print(f"\n  {'method':<12} {'D z':>10} {'E z':>10} {'D-E corr':>10}")
    print(f"  {'─' * 42}")
    z_D_vals = []
    z_E_vals = []
    for method in methods:
        r = sensitivity[method]
        zd = r["D"]["z_score"]
        ze = r["E"]["z_score"]
        z_D_vals.append(zd)
        z_E_vals.append(ze)
        print(f"  {method:<12} {zd:>10.2f} {ze:>10.2f} {r['DE_correlation']:>10.3f}")

    d_range = max(z_D_vals) - min(z_D_vals)
    e_range = max(z_E_vals) - min(z_E_vals)
    print(f"\n  z-score range across rep sets:  D: {d_range:.2f},  E: {e_range:.2f}")

    if d_range < 3.0 and e_range < 3.0:
        print(f"  ✓ ROBUST: z-scores stable across representative choices")
    elif d_range < 5.0:
        print(f"  ⚠ MODERATE sensitivity to representative choice")
    else:
        print(f"  ✗ SENSITIVE: representative choice significantly affects results")

    results["sensitivity_summary"] = {
        "D_z_range": d_range, "E_z_range": e_range,
        "D_z_values": {m: sensitivity[m]["D"]["z_score"] for m in methods},
        "E_z_values": {m: sensitivity[m]["E"]["z_score"] for m in methods},
    }

    # ── save summary JSON ──
    outpath = "results/publication_controls.json"
    os.makedirs("results", exist_ok=True)
    with open(outpath, "w") as f:
        json.dump(results, f, indent=2, default=str)
    print(f"\n  saved to {outpath}")

    # ── save raw null distributions for figure generation ──
    dist_path = "results/null_distributions.npz"
    np.savez_compressed(dist_path, **raw_dists)
    size_mb = os.path.getsize(dist_path) / (1024 * 1024)
    print(f"  saved raw distributions to {dist_path} ({size_mb:.1f} MB)")
    print(f"  arrays: {sorted(raw_dists.keys())}")

    print(f"\n[{datetime.now()}] done")


if __name__ == "__main__":
    main()
