#!/usr/bin/env python
# run_phase2_batch.py
# Triplet (n=3) batch: Dirichlet energy (E) + receiver/noise (N) with degeneracy-preserving nulls.

from __future__ import annotations
import os, sys, json, time, argparse, math, random
from typing import Dict, List, Tuple
import numpy as np
import pandas as pd
from multiprocessing import Pool

# Ensure we can import from src/ when run from repo root
sys.path.insert(0, os.path.abspath('src'))

# ---- our modules ----
from src.io.codon_io import build_sgc_map  # returns pd.DataFrame columns: codon, aa, is_stop
from src.io.aa_props_io import load_aa_props, prop_names, aa_prop_vector_map
from src.sims.codon_graph import build_graph_n
from src.metrics.dirichlet import dirichlet_energy
from src.receiver.thermo_noise import expected_property_distortion_noise_multi as noise_multi


# ----------------- helpers -----------------


def _codon_vectors(codon_df: pd.DataFrame, aa_vec_map: dict[str, list[float]]) -> dict[str, list[float]]:
    """
    Build codon→feature-vector map for *sense* codons only.
    Skips STOPs and rows with missing/unknown AA (NaN, '*', 'STOP').
    """
    out: dict[str, list[float]] = {}
    skipped = []

    # Normalize expected columns
    cols_needed = {"codon", "aa", "is_stop"}
    for c in cols_needed - set(codon_df.columns):
        if c == "is_stop":
            codon_df["is_stop"] = False
        else:
            raise ValueError(f"codon_df missing required column: {c}")

    for _, r in codon_df.iterrows():
        cod = str(r["codon"]).upper()
        aa_raw = r["aa"]

        # Drop NaN AA
        if pd.isna(aa_raw):
            skipped.append((cod, "NaN"))
            continue

        aa = str(aa_raw).strip().upper()

        # Skip STOPs by flag or token
        if bool(r["is_stop"]) or aa in {"*", "STOP"}:
            continue

        vec = aa_vec_map.get(aa)
        if vec is None:
            skipped.append((cod, aa))
            continue

        out[cod] = vec

    if skipped:
        print(f"[warn] skipped {len(skipped)} codons lacking a usable AA mapping (e.g., {skipped[:5]})")

    return out


def _degeneracy_counts(codon_df: pd.DataFrame) -> Dict[str, int]:
    """Counts of codons per amino acid (sense only)."""
    sense = codon_df[~codon_df['is_stop']]
    return sense.groupby('aa').size().to_dict()


def _degeneracy_preserving_shuffle(codon_df: pd.DataFrame, rng: random.Random) -> pd.DataFrame:
    """Return a new DataFrame with same codon list; STOP fixed; sense codons reassigned
    to AAs with the same per-AA counts as SGC (degeneracy-preserving)."""
    df = codon_df.copy()
    mask_stop = df['is_stop'].astype(bool)
    sense = df.loc[~mask_stop].copy()

    target_counts = sense.groupby('aa').size().to_dict()  # e.g., {'LEU':6,'SER':6,...}
    aa_list = []
    for aa, k in target_counts.items():
        aa_list.extend([aa] * k)
    rng.shuffle(aa_list)
    sense = sense.assign(aa=aa_list)
    df.loc[~mask_stop, 'aa'] = sense['aa'].values
    return df


def _wilson_ci(k: int, n: int, alpha: float = 0.05) -> Tuple[float, float]:
    """Wilson score interval for binomial proportion."""
    if n <= 0:
        return (0.0, 0.0)
    z = 1.959963984540054  # ~95%
    p = k / n
    denom = 1 + z*z/n
    center = (p + z*z/(2*n)) / denom
    half = (z/denom) * math.sqrt((p*(1-p)/n) + (z*z/(4*n*n)))
    lo = max(0.0, center - half)
    hi = min(1.0, center + half)
    return (lo, hi)


# ------------- worker payload (picklable) -------------

class WorkerConfig:
    __slots__ = ('codons', 'stops_idx', 'aa_vec_map', 'G', 'seed')
    def __init__(self, codon_df: pd.DataFrame, aa_vec_map: Dict[str, np.ndarray], G, seed: int):
        # split codons into sense vs stop for fast reassignment
        self.codons: List[str] = list(codon_df['codon'].astype(str))
        self.stops_idx: List[int] = [i for i, aa in enumerate(codon_df['aa']) if str(aa).upper() in ('*','STOP')]
        self.aa_vec_map = aa_vec_map
        self.G = G
        self.seed = seed


def _pack_for_eval(codon_df: pd.DataFrame) -> Tuple[List[str], List[int]]:
    codons = list(codon_df['codon'].astype(str))
    stops_idx = [i for i, aa in enumerate(codon_df['aa']) if str(aa).upper() in ('*','STOP')]
    return codons, stops_idx


def _random_trial(args) -> Tuple[float, float, int]:
    """
    One random degeneracy-preserving reassignment → returns (E, N, 1)
    """
    codon_df, aa_vec_map, G, seed = args  # pass by value; args prepared by main
    rng = random.Random(seed)
    r_df = _degeneracy_preserving_shuffle(codon_df, rng)

    # E: Dirichlet on graph
    codon_to_vec = _codon_vectors(r_df, aa_vec_map)
    E = dirichlet_energy(G, codon_to_vec)

    # N: receiver/noise with multi-property vectors
    N = noise_multi(r_df, aa_vec_map, RT=1.0, wobble_rules='strict', asymmetry=True)
    return (E, N, 1)


def _batched_tasks(codon_df: pd.DataFrame, aa_vec_map: Dict[str, np.ndarray], G, n_jobs: int, base_seed: int):
    for i in range(n_jobs):
        yield (codon_df, aa_vec_map, G, base_seed + i + 12345)


# ----------------- main -----------------

def main():
    ap = argparse.ArgumentParser(description="Phase-2 triplet batch (E & N) with degeneracy-preserving nulls.")
    ap.add_argument('--n', type=int, required=True, help='Number of random codes')
    ap.add_argument('--workers', type=int, default=max(1, os.cpu_count() or 1), help='Parallel workers')
    ap.add_argument('--out', type=str, required=True, help='Output JSON path')
    ap.add_argument('--seed', type=int, default=1234, help='Base RNG seed')
    ap.add_argument('--aa-parquet', type=str, default='data/processed/aa_props.parquet',
                    help='Pre-built AA features parquet (from build_aa_props.py); falls back to aa_props.csv')
    args = ap.parse_args()

    n_random = int(args.n)
    workers = int(args.workers)
    seed = int(args.seed)
    out_path = args.out

    # 0) Load AA properties — prefer parquet (22 props, 8 PCs) over CSV (6 props, 4 PCs)
    if os.path.exists(args.aa_parquet):
        aa_df = pd.read_parquet(args.aa_parquet)
        meta = {}
    else:
        aa_df, meta = load_aa_props(auto_select=True, corr_thresh=0.95, pca_var=0.97)
    feat_cols = prop_names(aa_df)
    # The parquet stores the canonical PCA feature space; do not z-score it again.
    aa_map = aa_prop_vector_map(aa_df, feature_cols=feat_cols, standardize=False)
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] AA feature cols: {feat_cols}")

    # 1) Load SGC codon table (codon, aa, is_stop) and build graph
    codon_df = build_sgc_map()  # uses your codon_table.csv inside src/io/codon_io.py
    # assure required columns exist
    for col in ('codon', 'aa', 'is_stop'):
        if col not in codon_df.columns:
            raise ValueError(f"build_sgc_map() missing column '{col}'")
    # codon graph (triplet, 4-letter alphabet)
    G = build_graph_n(3)

    # 2) Baselines for SGC
    codon_to_vec_sgc = _codon_vectors(codon_df, aa_map)
    E_sgc = dirichlet_energy(G, codon_to_vec_sgc)
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] E_sgc (Dirichlet): {E_sgc:.6f}")
    N_sgc = noise_multi(codon_df, aa_map, RT=1.0, wobble_rules='strict', asymmetry=True)
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] N_sgc (receiver multi-prop): {N_sgc:.6f}")

    # 3) Random nulls (degeneracy-preserving), parallel
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] Starting nulls: n={n_random:,}, workers={workers}, n_bases=3")
    count_le_E = 0
    count_le_N = 0
    processed = 0

    # Chunk by workers for fewer dispatches
    chunk = max(1, n_random // (workers * 4))
    seeds = (seed + 1000000*np.arange(n_random)).astype(int)

    def _consume(results_iter):
        nonlocal count_le_E, count_le_N, processed
        for i, (E, N, _) in enumerate(results_iter, 1):
            if E <= E_sgc:
                count_le_E += 1
            if N <= N_sgc:
                count_le_N += 1
            processed = i
            if (i % max(1, n_random // 8)) == 0 or i == n_random:
                print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] progress: {processed}/{n_random}")

    args_iter = ((codon_df, aa_map, G, int(seeds[i])) for i in range(n_random))
    if workers <= 1:
        _consume(_random_trial(args) for args in args_iter)
    else:
        try:
            with Pool(processes=workers) as pool:
                _consume(pool.imap_unordered(_random_trial, args_iter, chunksize=chunk))
        except PermissionError:
            print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] [warn] multiprocessing unavailable; falling back to serial execution")
            _consume(_random_trial(args) for args in args_iter)

    # Laplace-smoothed percentile (consistent with phase1 and n2n4 scripts)
    p_E = (count_le_E + 1) / (n_random + 1)
    p_N = (count_le_N + 1) / (n_random + 1)
    ciE = _wilson_ci(count_le_E, n_random)
    ciN = _wilson_ci(count_le_N, n_random)

    # 4) Write results
    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)
    result = {
        "n_random": n_random,
        "workers": workers,
        "seed": seed,
        "n_bases": 3,
        "features": feat_cols,
        "E_sgc": float(E_sgc),
        "N_sgc": float(N_sgc),
        "count_le_E": int(count_le_E),
        "count_le_N": int(count_le_N),
        "percentile_E": float(p_E),
        "percentile_N": float(p_N),
        "ci95_E": [float(ciE[0]), float(ciE[1])],
        "ci95_N": [float(ciN[0]), float(ciN[1])],
    }
    with open(out_path, "w") as f:
        json.dump(result, f, indent=2)
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] Wrote {out_path}")


if __name__ == "__main__":
    main()
