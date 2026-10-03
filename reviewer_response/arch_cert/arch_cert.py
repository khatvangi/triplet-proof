"""certify per-residue optima of D on H(n,4) for doublet (n=2) vs triplet (n=3).

D = (1/N_sense) * sum_{sense c} sum_{sense Hamming-1 nbr c'} dist(label(c), label(c'))
  = 2*H1/N_sense.

scenarios (each at n=2 and n=3):
  S1 chemistry-blind, dist = 1 for distinct classes, k in {4,8,10,16}, balanced sizes, no stops
  S2 chemistry, k=10, balanced sizes, no stops
  S3 chemistry, class sizes fixed to the projection code (doublet: 9 classes; triplet: RAA10, stops excluded)
  S4 chemistry, k=10, free sizes (each class >= 1); stops excluded for n=3 only

method: simulated annealing (swap moves for fixed sizes; single-codon recolor + swap for S4),
geometric cooling, then steepest descent over all swaps/recolors to a local optimum.

run from repo root:  python reviewer_response/arch_cert/arch_cert.py
"""
import itertools
import json
import math
import os
import sys
import time
from multiprocessing import get_context

import numpy as np
import pandas as pd
from numba import njit

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, ROOT)
from run_publication_controls import RAA10, build_doublet_table  # noqa: E402

ALPHABET = "ACGU"
STOPS = ("UAA", "UAG", "UGA")
REPS = {"X1": "A", "X2": "L", "X3": "W", "X4": "S", "X5": "N",
        "X6": "D", "X7": "R", "X8": "H", "X9": "P", "X10": "C"}
CLASS_ORDER = [f"X{i}" for i in range(1, 11)]
GATE_TRIPLET = 21.268639884706705
GATE_DOUBLET = 18.937298834660997
N_RESTARTS = 256
N_WORKERS = 48
OUT_JSON = os.path.join(HERE, "RESULTS.json")


# ── inputs ──

def loadClassVectors():
    """z-score pca1..pca8 across the 20 AAs (ddof=0), return {class: vector} for the 10 reps."""
    df = pd.read_csv(os.path.join(ROOT, "reviewer_response", "aa_props_export.csv"))
    cols = [f"pca{i}" for i in range(1, 9)]
    x = df[cols].to_numpy(float)
    x = (x - x.mean(axis=0)) / x.std(axis=0, ddof=0)
    aa_vec = {aa: x[i] for i, aa in enumerate(df["aa"].str.upper())}
    return {cls: aa_vec[aa] for cls, aa in REPS.items()}


def buildSenseGraph(n, excluded=()):
    """codons of H(n,4) minus excluded, with a padded (-1) neighbour array restricted to them."""
    codons = ["".join(t) for t in itertools.product(ALPHABET, repeat=n)]
    codons = [c for c in codons if c not in excluded]
    idx = {c: i for i, c in enumerate(codons)}
    deg = 3 * n
    nbr = -np.ones((len(codons), deg), dtype=np.int64)
    for c, i in idx.items():
        j = 0
        for p in range(n):
            for b in ALPHABET:
                if b != c[p]:
                    m = c[:p] + b + c[p + 1:]
                    if m in idx:
                        nbr[i, j] = idx[m]
                        j += 1
    return codons, nbr


def distMatrix(classes, vecs):
    k = len(classes)
    d = np.zeros((k, k))
    for a in range(k):
        for b in range(k):
            d[a, b] = np.linalg.norm(vecs[classes[a]] - vecs[classes[b]])
    return d


def balancedSizes(N, k):
    base, extra = divmod(N, k)
    return [base + 1 if i < extra else base for i in range(k)]


# ── numba kernels ──

@njit(cache=True)
def totalCost(lab, nbr, dist):
    """sum over nodes of sum over neighbours = 2*H1."""
    s = 0.0
    for i in range(nbr.shape[0]):
        for j in range(nbr.shape[1]):
            v = nbr[i, j]
            if v >= 0:
                s += dist[lab[i], lab[v]]
    return s


@njit(cache=True)
def nodeSum(lab, nbr, dist, i, c, skip):
    """sum_{nbr v of i, v != skip} dist(c, lab[v])."""
    s = 0.0
    for j in range(nbr.shape[1]):
        v = nbr[i, j]
        if v >= 0 and v != skip:
            s += dist[c, lab[v]]
    return s


@njit(cache=True)
def swapDelta(lab, nbr, dist, a, b):
    # change in 2*H1; the a-b edge (if any) is unchanged by the swap so it is skipped
    la, lb = lab[a], lab[b]
    d = nodeSum(lab, nbr, dist, a, lb, b) - nodeSum(lab, nbr, dist, a, la, b)
    d += nodeSum(lab, nbr, dist, b, la, a) - nodeSum(lab, nbr, dist, b, lb, a)
    return 2.0 * d


@njit(cache=True)
def recolorDelta(lab, nbr, dist, a, c):
    return 2.0 * (nodeSum(lab, nbr, dist, a, c, -1) - nodeSum(lab, nbr, dist, a, lab[a], -1))


@njit(cache=True)
def annealSwap(lab, nbr, dist, n_iter, seed):
    np.random.seed(seed)
    N = lab.shape[0]
    # initial temperature: mean |delta| over random swaps
    acc, cnt = 0.0, 0
    for _ in range(2000):
        a, b = np.random.randint(N), np.random.randint(N)
        if lab[a] != lab[b]:
            acc += abs(swapDelta(lab, nbr, dist, a, b))
            cnt += 1
    T = acc / max(cnt, 1) + 1e-12
    alpha = math.exp(math.log(1e-4) / n_iter)  # cool by 1e4 over the run
    for _ in range(n_iter):
        a, b = np.random.randint(N), np.random.randint(N)
        if lab[a] != lab[b]:
            dl = swapDelta(lab, nbr, dist, a, b)
            if dl <= 0.0 or np.random.random() < math.exp(-dl / T):
                t = lab[a]
                lab[a] = lab[b]
                lab[b] = t
        T *= alpha
    # steepest descent over all pairwise swaps
    while True:
        best, ba, bb = -1e-12, -1, -1
        for a in range(N):
            for b in range(a + 1, N):
                if lab[a] != lab[b]:
                    dl = swapDelta(lab, nbr, dist, a, b)
                    if dl < best:
                        best, ba, bb = dl, a, b
        if ba < 0:
            break
        t = lab[ba]
        lab[ba] = lab[bb]
        lab[bb] = t
    return lab


@njit(cache=True)
def annealRecolor(lab, nbr, dist, k, n_iter, seed):
    np.random.seed(seed)
    N = lab.shape[0]
    counts = np.zeros(k, dtype=np.int64)
    for i in range(N):
        counts[lab[i]] += 1
    acc, cnt = 0.0, 0
    for _ in range(2000):
        a, c = np.random.randint(N), np.random.randint(k)
        if c != lab[a]:
            acc += abs(recolorDelta(lab, nbr, dist, a, c))
            cnt += 1
    T = acc / max(cnt, 1) + 1e-12
    alpha = math.exp(math.log(1e-4) / n_iter)
    for _ in range(n_iter):
        # half the moves are size-preserving swaps: with recolor alone a singleton class
        # can never relocate (recoloring it would empty the class), so it freezes in place
        if np.random.random() < 0.5:
            a, b = np.random.randint(N), np.random.randint(N)
            if lab[a] != lab[b]:
                dl = swapDelta(lab, nbr, dist, a, b)
                if dl <= 0.0 or np.random.random() < math.exp(-dl / T):
                    t = lab[a]
                    lab[a] = lab[b]
                    lab[b] = t
            T *= alpha
            continue
        a, c = np.random.randint(N), np.random.randint(k)
        # recolor only if it keeps every class nonempty
        if c != lab[a] and counts[lab[a]] > 1:
            dl = recolorDelta(lab, nbr, dist, a, c)
            if dl <= 0.0 or np.random.random() < math.exp(-dl / T):
                counts[lab[a]] -= 1
                counts[c] += 1
                lab[a] = c
        T *= alpha
    # steepest descent over all legal recolors and all pairwise swaps
    while True:
        best, ba, bc, bb = -1e-12, -1, -1, -1
        for a in range(N):
            if counts[lab[a]] > 1:
                for c in range(k):
                    if c != lab[a]:
                        dl = recolorDelta(lab, nbr, dist, a, c)
                        if dl < best:
                            best, ba, bc, bb = dl, a, c, -1
            for b in range(a + 1, N):
                if lab[a] != lab[b]:
                    dl = swapDelta(lab, nbr, dist, a, b)
                    if dl < best:
                        best, ba, bc, bb = dl, a, -1, b
        if ba < 0:
            break
        if bb >= 0:
            t = lab[ba]
            lab[ba] = lab[bb]
            lab[bb] = t
        else:
            counts[lab[ba]] -= 1
            counts[bc] += 1
            lab[ba] = bc
    return lab


# ── worker ──

G = {}  # per-process job state, set before fork


def runRestart(seed):
    job = G["job"]
    rng = np.random.default_rng(seed)
    N, k = job["N"], job["k"]
    if job["free"]:
        lab = np.empty(N, dtype=np.int64)
        perm = rng.permutation(N)
        lab[perm[:k]] = np.arange(k)
        lab[perm[k:]] = rng.integers(0, k, N - k)
        lab = annealRecolor(lab, job["nbr"], job["dist"], k, job["n_iter"], seed)
    else:
        lab = rng.permutation(job["base_lab"]).astype(np.int64)
        lab = annealSwap(lab, job["nbr"], job["dist"], job["n_iter"], seed)
    D = totalCost(lab, job["nbr"], job["dist"]) / N
    return D, lab.tolist()


def runScenario(job, n_restarts):
    G["job"] = job
    t0 = time.time()
    with get_context("fork").Pool(N_WORKERS) as pool:
        out = pool.map(runRestart, [job["seed0"] + r for r in range(n_restarts)], chunksize=1)
    Ds = np.array([d for d, _ in out])
    best_i = int(np.argmin(Ds))
    best = float(Ds[best_i])
    lab = out[best_i][1]
    # sanity: the saved labeling respects the size constraints
    counts = np.bincount(lab, minlength=job["k"])
    assert counts.min() >= 1
    if not job["free"]:
        assert sorted(counts.tolist()) == sorted(np.bincount(job["base_lab"], minlength=job["k"]).tolist())
    return {
        "best_D": best,
        "hits": int(np.sum(np.abs(Ds - best) <= 1e-9)),
        "restarts": n_restarts,
        "p5": float(np.percentile(Ds, 5)),
        "p50": float(np.percentile(Ds, 50)),
        "max": float(Ds.max()),
        "N_sense": job["N"],
        "n_iter": job["n_iter"],
        "class_sizes_best": {job["classes"][c]: int(counts[c]) for c in range(job["k"])},
        "best_labeling": {job["codons"][i]: job["classes"][lab[i]] for i in range(job["N"])},
        "seconds": round(time.time() - t0, 1),
    }


# ── gates ──

def labelsFromMap(codons, cmap, classes):
    ci = {c: i for i, c in enumerate(classes)}
    return np.array([ci[cmap[c]] for c in codons], dtype=np.int64)


def main():
    vecs = loadClassVectors()
    tri = pd.read_csv(os.path.join(ROOT, "reviewer_response", "codon_table_export.csv"))
    tri["codon"] = tri["codon"].str.upper()
    stops = set(tri.loc[tri["is_stop"].astype(bool), "codon"])
    assert stops == set(STOPS), stops

    # triplet projection: RAA10 on the 61 sense codons
    tri_map = {r.codon: RAA10[r.aa.upper()] for r in tri.itertuples() if not r.is_stop}
    tri_codons, tri_nbr = buildSenseGraph(3, STOPS)
    tri_lab = labelsFromMap(tri_codons, tri_map, CLASS_ORDER)
    D_tri = totalCost(tri_lab, tri_nbr, distMatrix(CLASS_ORDER, vecs)) / len(tri_codons)

    # doublet projection via the pipeline's own majority-vote builder
    dbl = build_doublet_table(tri, RAA10)
    assert not dbl["is_stop"].any()
    dbl_map = dict(zip(dbl["codon"], dbl["aa"]))
    dbl_classes = [c for c in CLASS_ORDER if c in set(dbl_map.values())]
    dbl_codons, dbl_nbr = buildSenseGraph(2)
    dbl_lab = labelsFromMap(dbl_codons, dbl_map, dbl_classes)
    D_dbl = totalCost(dbl_lab, dbl_nbr, distMatrix(dbl_classes, vecs)) / len(dbl_codons)

    gates = {"triplet_projection_D": D_tri, "triplet_expected": GATE_TRIPLET,
             "triplet_pass": abs(D_tri - GATE_TRIPLET) <= 1e-9,
             "doublet_projection_D": D_dbl, "doublet_expected": GATE_DOUBLET,
             "doublet_pass": abs(D_dbl - GATE_DOUBLET) <= 1e-9,
             "doublet_classes": dbl_classes,
             "doublet_missing_class": [c for c in CLASS_ORDER if c not in dbl_classes]}
    print("gates:", json.dumps(gates, indent=1))
    results = {"gates": gates, "representatives": REPS, "scenarios": []}
    if not (gates["triplet_pass"] and gates["doublet_pass"]):
        results["aborted"] = "gate failure"
        json.dump(results, open(OUT_JSON, "w"), indent=1)
        sys.exit("ABORT: gate failure")

    # ── build scenario jobs ──
    jobs = []
    seed0 = 1000
    for n in (2, 3):
        codons, nbr = buildSenseGraph(n)
        N = len(codons)
        n_iter = 1000 * N * N
        for k in (4, 8, 10, 16):
            classes = [f"c{i}" for i in range(k)]
            dist = 1.0 - np.eye(k)
            base = np.repeat(np.arange(k), balancedSizes(N, k))
            jobs.append(dict(scenario="S1", n=n, k=k, codons=codons, nbr=nbr, N=N, classes=classes,
                             dist=dist, base_lab=base, free=False, n_iter=n_iter))
        dist = distMatrix(CLASS_ORDER, vecs)
        base = np.repeat(np.arange(10), balancedSizes(N, 10))
        jobs.append(dict(scenario="S2", n=n, k=10, codons=codons, nbr=nbr, N=N, classes=CLASS_ORDER,
                         dist=dist, base_lab=base, free=False, n_iter=5 * n_iter))
    # S3: sizes of the projection codes
    jobs.append(dict(scenario="S3", n=2, k=len(dbl_classes), codons=dbl_codons, nbr=dbl_nbr,
                     N=len(dbl_codons), classes=dbl_classes, dist=distMatrix(dbl_classes, vecs),
                     base_lab=np.sort(dbl_lab), free=False, n_iter=5000 * 16 * 16,
                     reference_D=D_dbl))
    jobs.append(dict(scenario="S3", n=3, k=10, codons=tri_codons, nbr=tri_nbr, N=len(tri_codons),
                     classes=CLASS_ORDER, dist=distMatrix(CLASS_ORDER, vecs),
                     base_lab=np.sort(tri_lab), free=False, n_iter=5000 * 61 * 61, reference_D=D_tri))
    # S4: free sizes
    for n, (codons, nbr) in ((2, (dbl_codons, dbl_nbr)), (3, (tri_codons, tri_nbr))):
        N = len(codons)
        jobs.append(dict(scenario="S4", n=n, k=10, codons=codons, nbr=nbr, N=N, classes=CLASS_ORDER,
                         dist=distMatrix(CLASS_ORDER, vecs), base_lab=None, free=True,
                         n_iter=10000 * N * N))

    # warm the jit in the parent so forked workers inherit compiled code
    j = jobs[0]
    annealSwap(j["base_lab"].astype(np.int64).copy(), j["nbr"], j["dist"], 10, 0)
    annealRecolor(np.arange(16) % 10, dbl_nbr, distMatrix(CLASS_ORDER, vecs), 10, 10, 0)

    t_start = time.time()
    for i, job in enumerate(jobs):
        job["seed0"] = seed0 + 100000 * i
        # chemistry landscapes are glassy (a probe showed hits saturate with schedule length),
        # so they get more restarts: 2048 at n=2 (cheap), 1024 at n=3 (512 for the 2x-long S4)
        n_r = N_RESTARTS if job["scenario"] == "S1" else (2048 if job["n"] == 2 else
                                                       512 if job["scenario"] == "S4" else 1024)
        # budget guard: drop to 64 restarts if we are running long
        if time.time() - t_start > 15 * 60:
            n_r = 64
        res = runScenario(job, n_r)
        res.update(scenario=job["scenario"], n=job["n"], k=job["k"])
        if "reference_D" in job:
            res["projection_code_D"] = job["reference_D"]
        if job["scenario"] == "S1":
            lk = math.log(job["k"], 4)
            if abs(lk - round(lk)) < 1e-12:
                res["analytic_3log4k"] = 3.0 * round(lk)
                res["equals_analytic"] = abs(res["best_D"] - res["analytic_3log4k"]) <= 1e-9
        results["scenarios"].append(res)
        results["elapsed_seconds"] = round(time.time() - t_start, 1)
        json.dump(results, open(OUT_JSON, "w"), indent=1)
        print(f"{job['scenario']} k={job['k']} n={job['n']}: best={res['best_D']:.9f} "
              f"hits={res['hits']}/{n_r} p5={res['p5']:.6f} p50={res['p50']:.6f} "
              f"({res['seconds']}s)", flush=True)


if __name__ == "__main__":
    main()
