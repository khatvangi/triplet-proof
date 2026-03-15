#!/usr/bin/env python3
"""
run_phase2_n2n4.py — Monte Carlo for doublet (n=2) and quadruplet (n=4)

Design:
  n=2 (doublet):
    • Build a reduced-AA mapping (RAA10 or RAA13).
    • Collapse SGC triplets → reduced classes via AA→RAA.
    • Induce a doublet baseline by majority vote over the 3rd base:
        class(doublet XY) = mode( class(XY{A,U,G,C}) ).
      Ties broken by a fixed deterministic order.
    • Random codes: preserve class counts across 16 doublets.

  n=4 (quadruplet):
    • Baseline mapping via projection to triplet:
        AA(b1 b2 b3 b4) := AA_SGC(b1 b2 b3)  (duplicate along b4)
    • Random codes: preserve AA counts across 256 4-mers
        (i.e., each AA count is 4× its triplet count).
    • NOTE: This evaluates "longer code with extra dimension" under the triplet
      AA usage, and randomized deg-preserving alternatives on the 4-mer graph.

Metrics:
  • Dirichlet energy (E) on the Hamming graph of size 4^n
  • Receiver/noise distortion (N) using generalized proxy with wobble positions:
      - n=2: wobble_positions = set()   (no wobble)
      - n=4: wobble_positions = {2}     (treat only position 3 as wobble)

Output JSON mirrors the triplet runner.
"""

import os
import sys
import json
import math
import argparse
import random
from datetime import datetime
from collections import Counter, defaultdict

sys.path.insert(0, os.path.abspath(os.path.dirname(__file__)))

import numpy as np
import pandas as pd
import multiprocessing as mp
from itertools import product

# repo imports
from src.io.codon_io import load_codon_table_csv
from src.io.aa_props_io import load_aa_props, aa_prop_vector_map, prop_names
from src.sims.codon_graph import build_graph_n, label_nodes_with_aa
from src.metrics.dirichlet import dirichlet_energy
from src.receiver.thermo_noise import expected_property_distortion_noise_multi_general as noise_general

RNA_BASES = ("A","U","G","C")

def ts(): return datetime.now().strftime("%Y-%m-%d %H:%M:%S")

def wilson_ci(k, n, z=1.96):
    if n == 0: return (0.0, 0.0)
    p = k/n
    denom = 1 + (z*z)/n
    center = (p + (z*z)/(2*n)) / denom
    margin = (z/denom) * math.sqrt((p*(1-p)/n) + ((z*z)/(4*n*n)))
    lo = max(0.0, center-margin); hi = min(1.0, center+margin)
    return (lo, hi)

# -----------------------------
# Reduced amino-acid mappings
# -----------------------------
# Simple, literature-aligned RAA10/RAA13 buckets (coarse but serviceable).
# You can swap with your preferred scheme later. Keys are one-letter AAs.
RAA10 = {
    "A":"X1","G":"X1",
    "V":"X2","L":"X2","I":"X2","M":"X2",
    "F":"X3","W":"X3","Y":"X3",
    "S":"X4","T":"X4",
    "N":"X5","Q":"X5",
    "D":"X6","E":"X6",
    "K":"X7","R":"X7",
    "H":"X8",
    "P":"X9",
    "C":"X10",
}
RAA13 = {
    "A":"Y1","G":"Y1",
    "V":"Y2","I":"Y2","L":"Y2",
    "M":"Y3",
    "F":"Y4","Y":"Y4",
    "W":"Y5",
    "S":"Y6","T":"Y6",
    "N":"Y7","Q":"Y7",
    "D":"Y8","E":"Y8",
    "K":"Y9","R":"Y9",
    "H":"Y10",
    "P":"Y11",
    "C":"Y12",
    # keep one slot open for flexibility / mapping stability if needed:
    # "U" or "O" could be routed to Y13 in extended sets
}

def get_raa_map(name: str):
    name = name.lower()
    if name == "raa10": return RAA10, "X"
    if name == "raa13": return RAA13, "Y"
    raise ValueError("RAA scheme not supported. Use raa10 or raa13.")

# -----------------------------
# Utilities
# -----------------------------
def build_kmer_df(k: int) -> pd.DataFrame:
    kmers = [''.join(p) for p in product(RNA_BASES, repeat=k)]
    return pd.DataFrame({"codon": kmers})

def project4_to_triplet(c4: str) -> str:
    # drop the 4th base
    return c4[:3]

def majority_vote(labels, tie_key=None):
    cnt = Counter(labels)
    maxc = max(cnt.values())
    winners = [lab for lab,c in cnt.items() if c == maxc]
    if len(winners) == 1:
        return winners[0]
    # deterministic tie-break
    winners.sort(key=(tie_key or (lambda x: x)))
    return winners[0]

def degeneracy_preserving_assign(codons, aa_counts, rng):
    """
    Assign each codon a label while preserving the total count per label in aa_counts.
    Returns dict codon->label.
    """
    labels = []
    for aa, count in aa_counts.items():
        labels.extend([aa] * count)
    if len(labels) != len(codons):
        raise ValueError("aa_counts total does not match number of codons.")
    labels = labels[:]
    rng.shuffle(labels)
    return dict(zip(codons, labels))

def _random_eval_pack(args_tuple):
    return _random_eval(*args_tuple)

# -----------------------------
# Workers
# -----------------------------
def _random_eval(seed, n_iter, codon_df, G, aa_map, E_sgc, N_sgc, wobble_positions):
    rng = random.Random(seed)
    codons = list(codon_df["codon"])
    # counts by AA/class for degeneracy preservation
    target_counts = Counter(codon_df["aa"])
    target_counts.pop("*", None)

    count_E = 0
    count_N = 0
    for _ in range(n_iter):
        sd = rng.randrange(1 << 63)

        # build a random assignment preserving AA/class totals
        rnd_map = degeneracy_preserving_assign(codons, target_counts, rng)

        # label graph
        label_nodes_with_aa(G, rnd_map, exclude_stops=True)

        # E
        prop_map_r = {
            node: (aa_map.get(G.nodes[node]["aa"]) if G.nodes[node]["aa"] else None)
            for node in G.nodes()
        }
        E = dirichlet_energy(G, prop_map_r)
        if E <= E_sgc:
            count_E += 1

        # N (noise)
        tmp = codon_df.copy()
        tmp["aa"] = tmp["codon"].map(lambda c: rnd_map.get(c, None))
        N = noise_general(tmp, aa_map, RT=1.0, wobble_positions=wobble_positions)
        if N <= N_sgc:
            count_N += 1

    return (count_E, count_N, n_iter)

# -----------------------------
# Main
# -----------------------------
def main():
    ap = argparse.ArgumentParser(description="Phase-2 n=2/n=4 Monte Carlo.")
    ap.add_argument("--n-bases", type=int, required=True, choices=[2,4], help="codon length: 2 or 4")
    ap.add_argument("--n", type=int, default=1_000_000, help="number of random codes")
    ap.add_argument("--workers", type=int, default=32, help="processes")
    ap.add_argument("--seed", type=int, default=2025, help="master RNG seed")
    ap.add_argument("--codon-csv", type=str, default="codon_table.csv", help="SGC triplet CSV (for baseline)")
    ap.add_argument("--raa", type=str, default="raa10", choices=["raa10","raa13"], help="doublet scheme")
    ap.add_argument("--aa-parquet", type=str, default="data/processed/aa_props.parquet", help="AA features")
    ap.add_argument("--out", type=str, default="results/phase2_n2n4.auto.json", help="output JSON")
    args = ap.parse_args()

    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)

    # Features
    if os.path.exists(args.aa_parquet):
        aa_df = pd.read_parquet(args.aa_parquet)
    else:
        from src.io.aa_props_io import load_aa_props
        aa_df, _ = load_aa_props(auto_select=True, corr_thresh=0.95, pca_var=0.97)
    # The parquet defines the canonical feature scale; reuse it directly.
    aa_map_full = aa_prop_vector_map(aa_df, standardize=False)

    # Load SGC triplets
    sgc = load_codon_table_csv(args.codon_csv)  # 64 rows, 'codon','aa'
    triplet_to_aa = dict(zip(sgc["codon"], sgc["aa"]))

    n = args.n_bases
    G = build_graph_n(n=n)
    rng = random.Random(args.seed)

    if n == 2:
        # ---- DOUBLETS via reduced AA classes ----
        raa_map, _prefix = get_raa_map(args.raa)

        # Project triplet AA to RAA class
        def aa_to_cls(a):
            if not a or a == "*": return "*"
            return raa_map.get(a, None)

        # For each doublet XY, derive class by majority over XY{A,U,G,C}
        d_kmers = [''.join(p) for p in product(RNA_BASES, repeat=2)]
        rows = []
        for xy in d_kmers:
            labels = []
            for z in RNA_BASES:
                aa = triplet_to_aa.get(xy + z, "*")
                cls = aa_to_cls(aa)
                if cls and cls != "*":
                    labels.append(cls)
            if not labels:
                # if all were stops (unlikely), assign a fallback (pick most common class overall)
                labels.append("X1" if args.raa == "raa10" else "Y1")
            cls_xy = majority_vote(labels)
            rows.append((xy, cls_xy))

        df = pd.DataFrame(rows, columns=["codon","aa"])
        # count per class for degeneracy-preserving randomization
        # (this creates the null with same class usage across 16 doublets)
        # Build AA vectors for classes by averaging member amino acids in aa_map_full
        # Make class vectors as centroids
        class_members = defaultdict(list)
        for aa,cls in raa_map.items():
            if aa in aa_map_full:
                class_members[cls].append(aa_map_full[aa])
        class_vec = {cls: np.mean(np.stack(v, axis=0), axis=0) for cls,v in class_members.items()}
        aa_map = class_vec  # swap in class vectors

        # Label graph with baseline
        label_nodes_with_aa(G, dict(zip(df["codon"], df["aa"])), exclude_stops=True)

        # Metrics for baseline
        prop_map = {node: (aa_map.get(G.nodes[node]["aa"]) if G.nodes[node]["aa"] else None) for node in G.nodes()}
        E_sgc = dirichlet_energy(G, prop_map)

        wobble_positions = set()   # no wobble in doublet model
        N_sgc = noise_general(df, aa_map, RT=1.0, wobble_positions=wobble_positions)

        codon_df = df  # for nulls

    elif n == 4:
        # ---- QUADRUPLETS by projection to SGC triplet along base4 ----
        q_kmers = [''.join(p) for p in product(RNA_BASES, repeat=4)]
        aa_assign = []
        for kmer in q_kmers:
            aa = triplet_to_aa.get(project4_to_triplet(kmer), "*")
            aa_assign.append(aa)
        df = pd.DataFrame({"codon": q_kmers, "aa": aa_assign})

        # feature map is full amino-acid map
        aa_map = aa_map_full

        # Label graph baseline
        label_nodes_with_aa(G, dict(zip(df["codon"], df["aa"])), exclude_stops=True)

        # Metrics baseline
        prop_map = {node: (aa_map.get(G.nodes[node]["aa"]) if G.nodes[node]["aa"] else None) for node in G.nodes()}
        E_sgc = dirichlet_energy(G, prop_map)

        wobble_positions = {2}  # treat only pos3 as wobble; pos4 strong
        N_sgc = noise_general(df, aa_map, RT=1.0, wobble_positions=wobble_positions)

        codon_df = df

    else:
        raise ValueError("n-bases must be 2 or 4.")

    print(f"[{ts()}] n={n} baseline computed. E_sgc={E_sgc:.6f}, N_sgc={N_sgc:.6f}")

    # ---- Monte Carlo (degeneracy-preserving by AA/class totals) ----
    sizes = [args.n // args.workers] * args.workers
    for i in range(args.n % args.workers): sizes[i] += 1
    seeds = [args.seed + i * 10007 for i in range(args.workers)]

    count_E = 0; count_N = 0; seen = 0
    codon_df = codon_df.copy()

    tasks = [(seeds[i], sizes[i], codon_df, G, aa_map, E_sgc, N_sgc, wobble_positions) for i in range(args.workers)]

    def _consume(results_iter):
        nonlocal count_E, count_N, seen
        for i, (cE, cN, ni) in enumerate(results_iter):
            count_E += cE
            count_N += cN
            seen += ni
            if i % max(1, args.workers // 8) == 0:
                print(f"[{ts()}] progress: {seen}/{args.n}")

    if args.workers <= 1:
        _consume(_random_eval_pack(task) for task in tasks)
    else:
        try:
            with mp.Pool(processes=args.workers) as pool:
                _consume(pool.imap_unordered(_random_eval_pack, tasks, chunksize=1))
        except PermissionError:
            print(f"[{ts()}] [warn] multiprocessing unavailable; falling back to serial execution")
            _consume(_random_eval_pack(task) for task in tasks)

    pct_E = (count_E + 1) / (args.n + 1)
    pct_N = (count_N + 1) / (args.n + 1)
    ciE = wilson_ci(count_E, args.n)
    ciN = wilson_ci(count_N, args.n)

    out = {
        "mode": f"n={n}",
        "raa": (args.raa if n==2 else None),
        "n_random": args.n,
        "workers": args.workers,
        "seed": args.seed,
        "features": prop_names(aa_df) if n==4 else "class-centroids",
        "E_sgc": float(E_sgc),
        "N_sgc": float(N_sgc),
        "count_le_E": int(count_E),
        "count_le_N": int(count_N),
        "percentile_E": float(pct_E),
        "percentile_N": float(pct_N),
        "ci95_E": [float(ciE[0]), float(ciE[1])],
        "ci95_N": [float(ciN[0]), float(ciN[1])],
    }
    with open(args.out, "w") as f:
        json.dump(out, f, indent=2)
    print(f"[{ts()}] Wrote {args.out}")


if __name__ == "__main__":
    main()
