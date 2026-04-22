#!/usr/bin/env python3
"""
Preemptive defense — compute two additional metrics alongside the manuscript's
canonical E (Dirichlet energy, f^T L f) and D (uniform-weighted noise distortion)
so we have ready answers if a reviewer asks:

  E_squared = f^T L^2 f
      Uses the SQUARED graph Laplacian. Tlusty's phase-transition analysis
      (2007 Eq. 3; 2010 Sec. 3) works with r^2 = (Δ - Δ_C)^2, so a reviewer
      grounded in his framework may ask for the squared-Laplacian z-score.
      Formula: Σ_{v ∈ sense} ||Σ_{u ∈ N_sense(v)} (f(v) - f(u))||^2

  D_boltz = Boltzmann-weighted noise distortion (from src/receiver/thermo_noise.py)
      Uses positional penalty (K=3.0), wobble softening (strength=0.4),
      softmax over per-codon single-position mismatches (RT=1.0). This
      is the "mechanistic" variant; our manuscript D uses uniform weighting
      as a minimal-assumption baseline. Key check: does weighting change
      the qualitative z-score picture?

Runs 100 000 degeneracy-preserving shuffles per condition (A, B, C).
Reports z-scores, count-le, percentiles, and cross-metric correlations.

Outputs:
  results/preempt_metrics.json
"""

import os, sys, json, random
from datetime import datetime
from collections import defaultdict
from itertools import product

sys.path.insert(0, os.path.abspath(os.path.dirname(__file__)))

import numpy as np
import pandas as pd
import multiprocessing as mp

from src.io.codon_io import load_codon_table_csv
from src.io.aa_props_io import load_aa_props, aa_prop_vector_map
from src.sims.codon_graph import build_graph_n
from src.receiver.thermo_noise import (
    expected_property_distortion_noise_multi_general,
)
from run_publication_controls import (
    RAA10, select_representatives, build_doublet_table, build_triplet_10aa,
)

RNA_BASES = ("A", "U", "G", "C")


# ────────────────────────────────────────────────────────────
# Metric computations
# ────────────────────────────────────────────────────────────

def build_adjacency(edges):
    """Return dict node -> set of Hamming-1 neighbors from an edge list."""
    adj = defaultdict(set)
    for u, v in edges:
        adj[u].add(v)
        adj[v].add(u)
    return dict(adj)


def compute_metrics(codon_df, aa_vec_map, edges, adjacency, codon_len, wobble_pos):
    """Compute all 4 metrics for a single code.

    Returns (E_uniform, E_squared, D_uniform, D_boltz).

    All metrics restrict to sense codons with known amino-acid mapping.
    """
    sense = codon_df[~codon_df["is_stop"].astype(bool)]
    c2a = dict(zip(sense["codon"].str.upper(), sense["aa"].str.upper()))

    # codon -> property vector (sense codons with known AA only)
    codon_vec = {}
    for c, aa in c2a.items():
        if aa in aa_vec_map:
            codon_vec[c] = aa_vec_map[aa]

    # ── E_uniform = Σ_edges ||f(u) − f(v)||^2 ──
    E_uni = 0.0
    for u, v in edges:
        pu = codon_vec.get(u)
        pv = codon_vec.get(v)
        if pu is not None and pv is not None:
            duv = pu - pv
            E_uni += float(np.dot(duv, duv))

    # ── E_squared = Σ_v ||(Lf)(v)||^2 ──
    # where (Lf)(v) = Σ_{u ∈ N_sense(v)} (f(v) − f(u))
    E_sq = 0.0
    for v, fv in codon_vec.items():
        lf_v = np.zeros_like(fv)
        for u in adjacency.get(v, ()):
            fu = codon_vec.get(u)
            if fu is not None:
                lf_v = lf_v + (fv - fu)
        E_sq += float(np.dot(lf_v, lf_v))

    # ── D_uniform = (1/|C|) Σ_c Σ_{c' ∈ N(c), nonsyn sense} ||f(c) − f(c')|| ──
    total_d = 0.0
    n_sense = 0
    for codon, aa in c2a.items():
        if aa not in aa_vec_map:
            continue
        vec_c = aa_vec_map[aa]
        n_sense += 1
        for pos in range(codon_len):
            for base in RNA_BASES:
                if base == codon[pos]:
                    continue
                mut = codon[:pos] + base + codon[pos + 1:]
                if mut in c2a:
                    mut_aa = c2a[mut]
                    if mut_aa != aa and mut_aa in aa_vec_map:
                        total_d += float(np.linalg.norm(vec_c - aa_vec_map[mut_aa]))
    D_uni = total_d / n_sense if n_sense > 0 else 0.0

    # ── D_boltz = Boltzmann-weighted via thermo_noise (general-codon-length) ──
    D_b = expected_property_distortion_noise_multi_general(
        codon_df, aa_vec_map,
        RT=1.0, wobble_positions=wobble_pos,
        asymmetry=True, K=3.0, wobble_strength=0.4,
    )

    return E_uni, E_sq, D_uni, D_b


# ────────────────────────────────────────────────────────────
# Null generation
# ────────────────────────────────────────────────────────────

def null_worker(args):
    """Compute 4 metrics for n_iter shuffles."""
    (seed, n_iter, codons, labels_template, aa_vec_map, edges,
     adjacency, codon_len, wobble_pos, codon_df_template) = args
    rng = random.Random(seed)

    E_uni_list, E_sq_list, D_uni_list, D_b_list = [], [], [], []

    for _ in range(n_iter):
        labels = labels_template[:]
        rng.shuffle(labels)
        # build a DataFrame compatible with thermo_noise (needs codon, aa, is_stop)
        # we replace the sense-codon aa column; stops stay as they are
        shuffled_df = codon_df_template.copy()
        shuffled_df.loc[~shuffled_df["is_stop"].astype(bool), "aa"] = labels

        E_uni, E_sq, D_uni, D_b = compute_metrics(
            shuffled_df, aa_vec_map, edges, adjacency, codon_len, wobble_pos
        )
        E_uni_list.append(E_uni)
        E_sq_list.append(E_sq)
        D_uni_list.append(D_uni)
        D_b_list.append(D_b)

    return E_uni_list, E_sq_list, D_uni_list, D_b_list


def run_null(codon_df, aa_vec_map, edges, adjacency, codon_len, wobble_pos,
             n_null, workers):
    """Parallel null: returns dict of metric -> list of null values."""
    sense = codon_df[~codon_df["is_stop"].astype(bool)]
    codons = sense["codon"].str.upper().tolist()
    labels = sense["aa"].str.upper().tolist()

    chunk = max(1, n_null // workers)
    seeds = [random.randint(0, 2**31 - 1) for _ in range(workers)]
    # last chunk absorbs remainder
    remainders = [chunk] * workers
    remainders[-1] = n_null - chunk * (workers - 1)

    args_list = [
        (seeds[i], remainders[i], codons, labels, aa_vec_map,
         edges, adjacency, codon_len, wobble_pos, codon_df)
        for i in range(workers)
    ]

    all_E_uni, all_E_sq, all_D_uni, all_D_b = [], [], [], []
    with mp.Pool(workers) as pool:
        for Eu, Es, Du, Db in pool.imap_unordered(null_worker, args_list):
            all_E_uni.extend(Eu)
            all_E_sq.extend(Es)
            all_D_uni.extend(Du)
            all_D_b.extend(Db)
    return {
        "E_uniform": np.asarray(all_E_uni),
        "E_squared": np.asarray(all_E_sq),
        "D_uniform": np.asarray(all_D_uni),
        "D_boltz":   np.asarray(all_D_b),
    }


def z_stats(sgc_val, null_vals):
    mu = float(np.mean(null_vals))
    sd = float(np.std(null_vals))
    z = (sgc_val - mu) / sd if sd > 0 else 0.0
    cnt_le = int(np.sum(null_vals <= sgc_val))
    pct = (cnt_le + 1) / (len(null_vals) + 1)
    return {
        "sgc": float(sgc_val),
        "null_mean": mu,
        "null_std": sd,
        "z_score": z,
        "count_le": cnt_le,
        "n_null": int(len(null_vals)),
        "percentile": pct,
    }


def analyze(label, codon_df, aa_vec_map, edges, adjacency, codon_len,
            wobble_pos, n_null, workers):
    sense = codon_df[~codon_df["is_stop"].astype(bool)]
    n_classes = len(set(sense["aa"].str.upper()))

    print(f"\n{'=' * 65}")
    print(f"  {label}  (n_bases={codon_len}, {n_classes} classes, "
          f"wobble={sorted(wobble_pos) if wobble_pos else 'none'})")
    print(f"{'=' * 65}")

    # SGC baselines
    E_uni, E_sq, D_uni, D_b = compute_metrics(
        codon_df, aa_vec_map, edges, adjacency, codon_len, wobble_pos
    )
    print(f"  SGC: E_uni={E_uni:.3f}, E_sq={E_sq:.3f}, "
          f"D_uni={D_uni:.4f}, D_boltz={D_b:.4f}")

    t0 = datetime.now()
    print(f"  running {n_null:,} shuffles on {workers} workers...")
    nulls = run_null(codon_df, aa_vec_map, edges, adjacency, codon_len,
                    wobble_pos, n_null, workers)
    elapsed = (datetime.now() - t0).total_seconds()
    print(f"  done in {elapsed:.1f} s")

    summary = {
        "n_bases": codon_len,
        "n_classes": n_classes,
        "wobble_positions": sorted(wobble_pos) if wobble_pos else [],
        "E_uniform": z_stats(E_uni, nulls["E_uniform"]),
        "E_squared": z_stats(E_sq, nulls["E_squared"]),
        "D_uniform": z_stats(D_uni, nulls["D_uniform"]),
        "D_boltz":   z_stats(D_b,  nulls["D_boltz"]),
    }

    # cross-metric correlations in the null
    correlations = {}
    for m1 in ["E_uniform", "E_squared", "D_uniform", "D_boltz"]:
        for m2 in ["E_uniform", "E_squared", "D_uniform", "D_boltz"]:
            if m1 < m2:
                correlations[f"{m1}__{m2}"] = float(
                    np.corrcoef(nulls[m1], nulls[m2])[0, 1]
                )
    summary["null_correlations"] = correlations

    for metric in ["E_uniform", "E_squared", "D_uniform", "D_boltz"]:
        s = summary[metric]
        print(f"  {metric:>10s}: z={s['z_score']:+.2f}  "
              f"μ={s['null_mean']:.3f}  σ={s['null_std']:.3f}  "
              f"count≤={s['count_le']}/{s['n_null']}  "
              f"p={s['percentile']:.2e}")
    print(f"  correlations in null:")
    for k, v in correlations.items():
        print(f"    r({k.replace('__', ', ')}) = {v:+.3f}")

    return summary


def main():
    import argparse
    ap = argparse.ArgumentParser(description="Preemptive-defense metrics")
    ap.add_argument("--n-null", type=int, default=100_000)
    ap.add_argument("--workers", type=int, default=32)
    ap.add_argument("--out", type=str, default="results/preempt_metrics.json")
    args = ap.parse_args()

    print(f"[{datetime.now()}] preemptive metrics: "
          f"n_null={args.n_null:,}, workers={args.workers}")

    triplet_df = load_codon_table_csv()
    aa_df, _ = load_aa_props("data/processed/aa_props.parquet")
    vec_map_20 = aa_prop_vector_map(aa_df)

    G2 = build_graph_n(2)
    G3 = build_graph_n(3)
    edges_n2 = list(G2.edges())
    edges_n3 = list(G3.edges())
    adj_n2 = build_adjacency(edges_n2)
    adj_n3 = build_adjacency(edges_n3)

    # closest-to-centroid representatives (main text)
    rep_closest, _ = select_representatives(vec_map_20, RAA10, "closest")

    doublet_df = build_doublet_table(triplet_df, RAA10)
    triplet_10aa_df = build_triplet_10aa(triplet_df, RAA10)

    results = {
        "timestamp": str(datetime.now()),
        "n_null": args.n_null,
    }

    # condition A: doublet, 10 AA. no wobble (n=2).
    results["A_10AA_n2"] = analyze(
        "A: 10 AA × n=2 (doublet)", doublet_df, rep_closest, edges_n2,
        adj_n2, 2, set(), args.n_null, args.workers,
    )
    # condition B: SGC, 20 AA. wobble at pos 3 (index 2).
    results["B_20AA_n3"] = analyze(
        "B: 20 AA × n=3 (SGC)", triplet_df, vec_map_20, edges_n3,
        adj_n3, 3, {2}, args.n_null, args.workers,
    )
    # condition C: 10 AA triplet. wobble at pos 3.
    results["C_10AA_n3"] = analyze(
        "C: 10 AA × n=3 (RAA10 triplet)", triplet_10aa_df, rep_closest,
        edges_n3, adj_n3, 3, {2}, args.n_null, args.workers,
    )

    # final summary table
    print(f"\n{'#' * 70}")
    print("  Z-SCORE COMPARISON ACROSS METRICS")
    print(f"{'#' * 70}\n")
    header = f"  {'condition':<20s} {'E_uniform':>11s} {'E_squared':>11s} {'D_uniform':>11s} {'D_boltz':>11s}"
    print(header)
    print("  " + "─" * len(header.strip()))
    for k, lab in [("A_10AA_n2", "A: 10AA n=2"),
                   ("C_10AA_n3", "C: 10AA n=3"),
                   ("B_20AA_n3", "B: 20AA n=3 (SGC)")]:
        r = results[k]
        print(f"  {lab:<20s} "
              f"{r['E_uniform']['z_score']:>+11.2f} "
              f"{r['E_squared']['z_score']:>+11.2f} "
              f"{r['D_uniform']['z_score']:>+11.2f} "
              f"{r['D_boltz']['z_score']:>+11.2f}")

    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    with open(args.out, "w") as f:
        json.dump(results, f, indent=2, default=str)
    print(f"\n[{datetime.now()}] saved: {args.out}")


if __name__ == "__main__":
    main()
