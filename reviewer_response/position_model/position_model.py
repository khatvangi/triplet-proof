"""position-weighted, usage-weighted doublet-vs-triplet optima.

Dpos = sum_a f_a (1/n_a) sum_{c in block(a)} sum_p w_p sum_{c' = c with position p changed, c' sense}
       dist(a, label(c'))

this is usage_cert's Dw with every arc multiplied by the weight of the position it changes. an
arc and its reverse change the same position, so a single per-arc weight matrix W[i, j] (weight
of the arc i -> nbr[i, j]) serves both directions in the incremental updates.

decoder hypotheses (the doublet weights are an ASSUMPTION: no doublet decoder exists):
  H-inv  : doublet (w1, w2),  triplet (w1, w2, w3)
  H-term : doublet (w1, w3),  triplet (w1, w2, w3)   (last-read base inherits the wobble rate)
with w2 = 1, r1 = w1 in {1, 3}, r3 = w3 in {1, 2, 5, 10, 30, 100}.

the triplet weights are the same under both hypotheses, so one n=3 search per (r1, r3, usage)
serves both. the n=2 searches are indexed by the doublet weight pair (r1, x): H-inv uses x = 1,
H-term uses x = r3, so 12 distinct n=2 searches per usage cover both hypotheses.

kernels are copied from usage_cert/usage_cert.py with the weight factor added; usage vectors,
class vectors and graphs come from the same sources (arch_cert, the saved Swiss-Prot page).

run from repo root:  python reviewer_response/position_model/position_model.py
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
# keep numba caches and .pyc files out of usage_cert/ and arch_cert/
os.environ["NUMBA_CACHE_DIR"] = os.path.join(HERE, "numba_cache")
sys.dont_write_bytecode = True

import itertools  # noqa: E402
import json  # noqa: E402
import math  # noqa: E402
import time  # noqa: E402
from multiprocessing import get_context  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from numba import njit  # noqa: E402

ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "reviewer_response", "arch_cert"))
sys.path.insert(0, os.path.join(ROOT, "reviewer_response", "usage_cert"))
from run_publication_controls import RAA10, build_doublet_table  # noqa: E402
from arch_cert import STOPS, CLASS_ORDER, REPS, loadClassVectors, buildSenseGraph, distMatrix, labelsFromMap  # noqa: E402
import usage_cert as UC  # noqa: E402  (reference kernels + usage builders, imported read-only)

N_WORKERS = 60
N_RANDOM, N_WARM = 128, 32
N_RANDOM_FLOOR, N_WARM_FLOOR = 64, 16
BUDGET_S = 23 * 60
R1S = (1.0, 3.0)
R3S = (1.0, 2.0, 5.0, 10.0, 30.0, 100.0)
OUT_JSON = os.path.join(HERE, "RESULTS.json")
OUT_LAB = os.path.join(HERE, "best_labelings.json")
USAGE_CERT_JSON = os.path.join(ROOT, "reviewer_response", "usage_cert", "RESULTS.json")


# ── weighted numba kernels (usage_cert kernels with the per-arc weight W) ──

@njit(cache=True)
def arcSums(lab, nbr, W, dist, k):
    T = np.zeros(k)
    for i in range(nbr.shape[0]):
        a = lab[i]
        for j in range(nbr.shape[1]):
            v = nbr[i, j]
            if v >= 0:
                T[a] += W[i, j] * dist[a, lab[v]]
    return T


@njit(cache=True)
def evalD(T, cnt, f):
    s = 0.0
    for a in range(T.shape[0]):
        if cnt[a] > 0:
            s += f[a] * T[a] / cnt[a]
    return s


@njit(cache=True)
def fullD(lab, nbr, W, dist, f, k):
    cnt = np.zeros(k, dtype=np.int64)
    for i in range(lab.shape[0]):
        cnt[lab[i]] += 1
    return evalD(arcSums(lab, nbr, W, dist, k), cnt, f)


@njit(cache=True)
def applyRecolor(lab, nbr, W, dist, T, cnt, i, b):
    a = lab[i]
    for j in range(nbr.shape[1]):
        v = nbr[i, j]
        if v >= 0:
            lv = lab[v]
            w = W[i, j]  # the reverse arc v -> i changes the same position
            T[a] -= w * dist[a, lv]
            T[b] += w * dist[b, lv]
            T[lv] += w * (dist[lv, b] - dist[lv, a])
    cnt[a] -= 1
    cnt[b] += 1
    lab[i] = b


@njit(cache=True)
def recolorValue(lab, nbr, W, dist, f, T, cnt, i, b, sT, sc):
    sT[:] = T
    sc[:] = cnt
    a = lab[i]
    for j in range(nbr.shape[1]):
        v = nbr[i, j]
        if v >= 0:
            lv = lab[v]
            w = W[i, j]
            sT[a] -= w * dist[a, lv]
            sT[b] += w * dist[b, lv]
            sT[lv] += w * (dist[lv, b] - dist[lv, a])
    sc[a] -= 1
    sc[b] += 1
    return evalD(sT, sc, f)


@njit(cache=True)
def swapValue(lab, nbr, W, dist, f, T, cnt, i, j, sT, sc):
    sT[:] = T
    sc[:] = cnt
    li, lj = lab[i], lab[j]
    applyRecolor(lab, nbr, W, dist, sT, sc, i, lj)
    applyRecolor(lab, nbr, W, dist, sT, sc, j, li)
    val = evalD(sT, sc, f)
    lab[i] = li
    lab[j] = lj
    return val


@njit(cache=True)
def descend(lab, nbr, W, dist, f, k):
    """steepest descent over all single recolors (every class kept nonempty) and all swaps."""
    N = lab.shape[0]
    cnt = np.zeros(k, dtype=np.int64)
    for i in range(N):
        cnt[lab[i]] += 1
    sT = np.zeros(k)
    sc = np.zeros(k, dtype=np.int64)
    n_moves = 0
    while True:
        T = arcSums(lab, nbr, W, dist, k)
        cur = evalD(T, cnt, f)
        # relative threshold: at r3 = 100 the objective is ~100x larger than at w = 1
        best, kind, bi, bj = -1e-12 * max(1.0, abs(cur)), -1, -1, -1
        for i in range(N):
            if cnt[lab[i]] > 1:
                for b in range(k):
                    if b != lab[i]:
                        d = recolorValue(lab, nbr, W, dist, f, T, cnt, i, b, sT, sc) - cur
                        if d < best:
                            best, kind, bi, bj = d, 0, i, b
        for i in range(N):
            for j in range(i + 1, N):
                if lab[i] != lab[j]:
                    d = swapValue(lab, nbr, W, dist, f, T, cnt, i, j, sT, sc) - cur
                    if d < best:
                        best, kind, bi, bj = d, 1, i, j
        if kind < 0:
            break
        if kind == 0:
            cnt[lab[bi]] -= 1
            cnt[bj] += 1
            lab[bi] = bj
        else:
            t = lab[bi]
            lab[bi] = lab[bj]
            lab[bj] = t
        n_moves += 1
    return n_moves


@njit(cache=True)
def anneal(lab, nbr, W, dist, f, k, n_iter, t_scale, seed):
    """single-codon recolor SA (t_scale = 0: descent only), then steepest descent."""
    np.random.seed(seed)
    N = lab.shape[0]
    cnt = np.zeros(k, dtype=np.int64)
    for i in range(N):
        cnt[lab[i]] += 1
    sT = np.zeros(k)
    sc = np.zeros(k, dtype=np.int64)
    if t_scale > 0.0:
        T = arcSums(lab, nbr, W, dist, k)
        cur = evalD(T, cnt, f)
        acc, m = 0.0, 0
        for _ in range(2000):
            i, b = np.random.randint(N), np.random.randint(k)
            if b != lab[i] and cnt[lab[i]] > 1:
                acc += abs(recolorValue(lab, nbr, W, dist, f, T, cnt, i, b, sT, sc) - cur)
                m += 1
        temp = t_scale * acc / max(m, 1) + 1e-12
        alpha = math.exp(math.log(1e-4) / n_iter)
        for it in range(n_iter):
            i, b = np.random.randint(N), np.random.randint(k)
            if b != lab[i] and cnt[lab[i]] > 1:
                d = recolorValue(lab, nbr, W, dist, f, T, cnt, i, b, sT, sc) - cur
                if d <= 0.0 or np.random.random() < math.exp(-d / temp):
                    applyRecolor(lab, nbr, W, dist, T, cnt, i, b)
                    cur += d
            temp *= alpha
            if (it & 0xFFFFF) == 0:
                T = arcSums(lab, nbr, W, dist, k)
                cur = evalD(T, cnt, f)
    return descend(lab, nbr, W, dist, f, k)


# ── graph helpers ──

def arcPositions(codons, nbr):
    """pos[i, j] = 0-based position changed by arc i -> nbr[i, j] (-1 for padding)."""
    pos = -np.ones(nbr.shape, dtype=np.int64)
    for i in range(nbr.shape[0]):
        for j in range(nbr.shape[1]):
            v = nbr[i, j]
            if v >= 0:
                diff = [p for p in range(len(codons[i])) if codons[i][p] != codons[v][p]]
                assert len(diff) == 1
                pos[i, j] = diff[0]
    return pos


def weightMatrix(pos, w):
    w = np.asarray(w, dtype=float)
    return np.where(pos >= 0, w[np.maximum(pos, 0)], 0.0)


def bruteD(lab, codons, w, dist, f, k):
    """independent pure-python Dpos straight from the definition (gates only)."""
    idx = {c: i for i, c in enumerate(codons)}
    cnt = np.bincount(lab, minlength=k)
    s = 0.0
    for i, c in enumerate(codons):
        a = lab[i]
        for p in range(len(c)):
            for b in "ACGU":
                if b != c[p]:
                    m = c[:p] + b + c[p + 1:]
                    if m in idx:
                        s += f[a] / cnt[a] * w[p] * dist[a, lab[idx[m]]]
    return s


def pos3Synonymy(lab, nbr, pos, f, k):
    """fraction of position-3 arcs that are synonymous: raw, and usage-weighted (f_a/n_a)."""
    cnt = np.bincount(lab, minlength=k)
    n_all = n_syn = wa = ws = 0.0
    for i in range(nbr.shape[0]):
        for j in range(nbr.shape[1]):
            if pos[i, j] == 2:
                a = lab[i]
                syn = lab[nbr[i, j]] == a
                n_all += 1
                n_syn += syn
                wa += f[a] / cnt[a]
                ws += syn * f[a] / cnt[a]
    return {"raw": n_syn / n_all, "usage_weighted": ws / wa, "n_arcs": int(n_all)}


# ── worker ──

G = {}


def runRestart(task):
    seed, src, rank = task
    job = G["job"]
    rng = np.random.default_rng(seed)
    if src < 0:
        lab, t_scale = UC.randomStart(job["N"], job["k"], rng), 1.0
    else:
        lab = job["warm_labs"][src].copy()
        # rank 0 of each warm source is pure descent; the rest span 1e-3 .. 1e-1
        nr = job["warm_per_src"]
        t_scale = 0.0 if rank == 0 else 10.0 ** (-3.0 + 2.0 * (rank - 1) / max(nr - 2, 1))
    anneal(lab, job["nbr"], job["W"], job["dist"], job["f"], job["k"], job["n_iter"], t_scale, seed)
    return fullD(lab, job["nbr"], job["W"], job["dist"], job["f"], job["k"]), lab.tolist()


def summarise(vals, best):
    v = np.array(vals)
    tol = 1e-9 * max(1.0, abs(best))
    return {"restarts": len(v), "best": float(v.min()), "hits": int(np.sum(np.abs(v - best) <= tol)),
            "p50": float(np.percentile(v, 50))}


def runSearch(job, n_random, n_warm):
    """n_warm restarts split evenly across job['warm_labs'] sources."""
    G["job"] = job
    t0 = time.time()
    n_src = len(job["warm_labs"])
    per = -(-n_warm // n_src) if n_src else 0  # ceil: never fewer than n_warm in total
    job["warm_per_src"] = per
    tasks = [(job["seed0"] + r, -1, 0) for r in range(n_random)]
    tasks += [(job["seed0"] + 50000 + 1000 * s + r, s, r) for s in range(n_src) for r in range(per)]
    with get_context("fork").Pool(N_WORKERS) as pool:
        out = pool.map(runRestart, tasks, chunksize=1)
    vals = [d for d, _ in out]
    bi = int(np.argmin(vals))
    best, lab = float(vals[bi]), np.array(out[bi][1], dtype=np.int64)
    assert np.bincount(lab, minlength=job["k"]).min() >= 1
    res = {"best": best, "random": summarise(vals[:n_random], best),
           "best_from": "random" if bi < n_random else f"warm:{job['warm_names'][(bi - n_random) // per]}",
           "n_iter": job["n_iter"], "seconds": round(time.time() - t0, 1),
           "class_sizes_best": dict(zip(CLASS_ORDER, np.bincount(lab, minlength=job["k"]).tolist()))}
    if per:
        res["warm"] = {}
        for s in range(n_src):
            res["warm"][job["warm_names"][s]] = summarise(vals[n_random + s * per:n_random + (s + 1) * per], best)
    return res, lab


# ── main ──

def main():
    t_start = time.time()
    vecs = loadClassVectors()
    chem = distMatrix(CLASS_ORDER, vecs)
    comp, prot_meta = UC.loadProteome()
    assert not prot_meta["fallback_used"], prot_meta
    uc = json.load(open(USAGE_CERT_JSON))
    assert prot_meta["release"] == uc["proteome"]["release"]

    tri = pd.read_csv(os.path.join(ROOT, "reviewer_response", "codon_table_export.csv"))
    tri["codon"] = tri["codon"].str.upper()
    tri_map = {r.codon: RAA10[r.aa.upper()] for r in tri.itertuples() if not r.is_stop}
    c61, nbr61 = buildSenseGraph(3, STOPS)
    c64, nbr64 = buildSenseGraph(3)
    c16, nbr16 = buildSenseGraph(2)
    pos61, pos64, pos16 = arcPositions(c61, nbr61), arcPositions(c64, nbr64), arcPositions(c16, nbr16)
    tri_lab = labelsFromMap(c61, tri_map, CLASS_ORDER)
    dbl = build_doublet_table(tri, RAA10)
    dbl_lab = labelsFromMap(c16, dict(zip(dbl["codon"], dbl["aa"])), CLASS_ORDER)
    dbl_present = sorted(set(dbl_lab.tolist()))

    fU = UC.uniformUsage(10, list(range(10)))
    fP = UC.classUsage(comp, set(CLASS_ORDER))
    fU9 = UC.uniformUsage(10, dbl_present)
    fP9 = UC.classUsage(comp, {CLASS_ORDER[i] for i in dbl_present})
    assert np.allclose(fP, [uc["usage_P_classes"][c] for c in CLASS_ORDER], atol=1e-15, rtol=0)

    # ── gate A: w = 1 reproduces usage_cert Dw (kernel call and stored RESULTS values) ──
    ref = uc["reference_codes"]
    W61_1, W16_1 = weightMatrix(pos61, (1, 1, 1)), weightMatrix(pos16, (1, 1))
    cases = [("sgc_triplet_U", fullD(tri_lab, nbr61, W61_1, chem, fU, 10),
              UC.fullDw(tri_lab, nbr61, chem, fU, 10), ref["sgc_triplet_projection_stops_excluded"]["Dw_U"]),
             ("sgc_triplet_P", fullD(tri_lab, nbr61, W61_1, chem, fP, 10),
              UC.fullDw(tri_lab, nbr61, chem, fP, 10), ref["sgc_triplet_projection_stops_excluded"]["Dw_P"]),
             ("doublet_proj_U", fullD(dbl_lab, nbr16, W16_1, chem, fU9, 10),
              UC.fullDw(dbl_lab, nbr16, chem, fU9, 10), ref["doublet_projection"]["Dw_U"]),
             ("doublet_proj_P", fullD(dbl_lab, nbr16, W16_1, chem, fP9, 10),
              UC.fullDw(dbl_lab, nbr16, chem, fP9, 10), ref["doublet_projection"]["Dw_P"])]
    gateA = {name: {"Dpos_w1": a, "usage_cert_kernel": b, "usage_cert_RESULTS": c,
                    "absdiff_max": max(abs(a - b), abs(a - c))} for name, a, b, c in cases}
    gateA["pass"] = all(v["absdiff_max"] <= 1e-9 for v in gateA.values())

    # ── gate B: weighted kernel vs brute-force definition, and incremental recolor/swap ──
    rng = np.random.default_rng(11)
    diffs, inc = [], []
    for _ in range(40):
        for codons, nbr, pos, n in ((c61, nbr61, pos61, 3), (c64, nbr64, pos64, 3), (c16, nbr16, pos16, 2)):
            w = rng.uniform(0.1, 100.0, n)
            W = weightMatrix(pos, w)
            lab = UC.randomStart(len(codons), 10, rng)
            f = rng.dirichlet(np.ones(10))
            d = fullD(lab, nbr, W, chem, f, 10)
            diffs.append(abs(d - bruteD(lab, codons, w, chem, f, 10)) / max(1.0, d))
            T = arcSums(lab, nbr, W, chem, 10)
            cnt = np.bincount(lab, minlength=10)
            sT, sc = np.zeros(10), np.zeros(10, dtype=np.int64)
            i, b = int(rng.integers(len(codons))), int(rng.integers(10))
            if b != lab[i] and cnt[lab[i]] > 1:
                v = recolorValue(lab, nbr, W, chem, f, T, cnt, i, b, sT, sc)
                l2 = lab.copy()
                l2[i] = b
                inc.append(abs(v - fullD(l2, nbr, W, chem, f, 10)) / max(1.0, v))
            i, j = rng.choice(len(codons), 2, replace=False)
            v = swapValue(lab, nbr, W, chem, f, T, cnt, i, j, sT, sc)
            l2 = lab.copy()
            l2[i], l2[j] = lab[j], lab[i]
            inc.append(abs(v - fullD(l2, nbr, W, chem, f, 10)) / max(1.0, v))
    gateB = {"n_cases": len(diffs), "max_reldiff_brute": float(max(diffs)),
             "n_incremental": len(inc), "max_reldiff_incremental": float(max(inc)), "tolerance": 1e-12}
    gateB["pass"] = gateB["max_reldiff_brute"] <= 1e-12 and gateB["max_reldiff_incremental"] <= 1e-12

    # ── gate C: under H-inv weights the silent lift is exact, Dpos3(L3) = Dpos2(L2) ──
    cdiffs = []
    for _ in range(100):
        L2 = rng.integers(0, 10, 16)
        present = sorted(set(L2.tolist()))
        f = np.zeros(10)
        f[present] = rng.dirichlet(np.ones(len(present)))
        w = (rng.uniform(0.1, 10), 1.0, rng.uniform(0.1, 100))
        L3 = UC.liftLabels(L2, c16, c64)
        d2 = fullD(L2, nbr16, weightMatrix(pos16, w[:2]), chem, f, 10)
        d3 = fullD(L3, nbr64, weightMatrix(pos64, w), chem, f, 10)
        cdiffs.append(abs(d2 - d3) / max(1.0, d2))
    gateC = {"n_random_labelings": 100, "max_reldiff": float(max(cdiffs)), "tolerance": 1e-12}
    gateC["pass"] = gateC["max_reldiff"] <= 1e-12

    results = {"gateA": gateA, "gateB": gateB, "gateC": gateC, "proteome": prot_meta,
               "representatives": REPS, "usage_U": fU.tolist(), "usage_P": fP.tolist(),
               "config": {"R1S": R1S, "R3S": R3S, "N_RANDOM": N_RANDOM, "N_WARM": N_WARM,
                          "floors": [N_RANDOM_FLOOR, N_WARM_FLOOR], "workers": N_WORKERS,
                          "budget_s": BUDGET_S, "n_iter_n2": "10000*N^2", "n_iter_n3": "2000*N^2"},
               "n2": [], "n3": [], "sgc": [], "secondary": [], "reductions": []}
    print("gateA:", gateA["pass"], {k: v["absdiff_max"] for k, v in gateA.items() if k != "pass"}, flush=True)
    print("gateB:", gateB, "\ngateC:", gateC, flush=True)
    best_labs = {}

    def dump():
        results["elapsed_seconds"] = round(time.time() - t_start, 1)
        json.dump(results, open(OUT_JSON + ".part", "w"), indent=1)
        os.replace(OUT_JSON + ".part", OUT_JSON)
        json.dump(best_labs, open(OUT_LAB + ".part", "w"), indent=1)
        os.replace(OUT_LAB + ".part", OUT_LAB)

    if not (gateA["pass"] and gateB["pass"] and gateC["pass"]):
        results["aborted"] = "gate failure"
        dump()
        sys.exit("ABORT: gate failure")

    # warm the jit before forking
    anneal(UC.randomStart(16, 10, np.random.default_rng(0)), nbr16, W16_1, chem, fU, 10, 100, 1.0, 0)

    # ── SGC under every triplet weight setting (61 sense codons) ──
    for usage, f in (("U", fU), ("P", fP)):
        for r1 in R1S:
            for r3 in R3S:
                w = (r1, 1.0, r3)
                results["sgc"].append({"usage": usage, "r1": r1, "r3": r3,
                                       "Dpos": fullD(tri_lab, nbr61, weightMatrix(pos61, w), chem, f, 10),
                                       "pos3_synonymy": pos3Synonymy(tri_lab, nbr61, pos61, f, 10)})
    dump()

    # cell order: U then P; within a usage, all n=2 first (cheap), then n=3 by (r1, r3)
    n_cells_total = 2 * (len(R1S) * len(R3S)) * 2
    seed0 = 7000

    def plan(n):
        """restart counts for the next cell; cut to the floors if the clock is running down.
        reserve the last measured full-size n=3 cell time for every n=3 cell still to run."""
        left = BUDGET_S - (time.time() - t_start)
        n3_left = n_cells_total // 2 - sum(1 for _ in results["n3"])
        if n == 3 and left < n3_left * results.get("_n3_sec", 0.0):
            return N_RANDOM_FLOOR, N_WARM_FLOOR, True
        return N_RANDOM, N_WARM, False

    for usage, f in (("U", fU), ("P", fP)):
        best2 = {}
        # n = 2, doublet weights (r1, x)
        for r1 in R1S:
            prev = None
            for x in R3S:
                w = (r1, x)
                W = weightMatrix(pos16, w)
                warm, names = ([prev], ["n2_neighbour"]) if prev is not None else ([], [])
                job = dict(nbr=nbr16, W=W, dist=chem, f=f, k=10, N=16, n_iter=10000 * 256,
                           seed0=seed0, warm_labs=warm, warm_names=names)
                seed0 += 100000
                nr, nw, cut = plan(2)
                res, lab = runSearch(job, nr, nw if warm else 0)
                res.update(usage=usage, n=2, w=list(w), restarts_cut=cut)
                if cut:
                    results["reductions"].append(f"n2 {usage} w={w}: {nr}+{nw}")
                best2[w] = {"D": res["best"], "lab": lab}
                prev = lab
                results["n2"].append(res)
                best_labs[f"n2_{usage}_w{w}"] = {c16[i]: CLASS_ORDER[lab[i]] for i in range(16)}
                dump()
                print(f"n=2 {usage} w={w}: best={res['best']:.9f} random {res['random']['hits']}/"
                      f"{res['random']['restarts']} ({res['seconds']}s, t={time.time() - t_start:.0f}s)", flush=True)

        # n = 3, triplet weights (r1, 1, r3); warm = lift(H-inv n2), lift(H-term n2), neighbour n3
        prev3 = None
        for r1 in R1S:
            for r3 in R3S:
                w = (r1, 1.0, r3)
                W = weightMatrix(pos64, w)
                lift_inv = UC.liftLabels(best2[(r1, 1.0)]["lab"], c16, c64)
                lift_term = UC.liftLabels(best2[(r1, r3)]["lab"], c16, c64)
                warm, names = [lift_inv, lift_term], ["lift_Hinv", "lift_Hterm"]
                if prev3 is not None:
                    warm.append(prev3[1])
                    names.append(f"n3_neighbour{prev3[0]}")
                job = dict(nbr=nbr64, W=W, dist=chem, f=f, k=10, N=64, n_iter=2000 * 64 * 64,
                           seed0=seed0, warm_labs=warm, warm_names=names)
                seed0 += 100000
                nr, nw, cut = plan(3)
                # at least 1 descent-only + some annealed restart per warm source
                nw = max(nw, 2 * len(warm))
                res, lab = runSearch(job, nr, nw)
                if not cut:
                    results["_n3_sec"] = res["seconds"]
                if cut:
                    results["reductions"].append(f"n3 {usage} w={w}: {nr}+{nw}")
                D3 = res["best"]
                d2_inv, d2_term = best2[(r1, 1.0)]["D"], best2[(r1, r3)]["D"]
                lift_inv_val = fullD(lift_inv, nbr64, W, chem, f, 10)
                lift_term_val = fullD(lift_term, nbr64, W, chem, f, 10)
                sgc = next(s for s in results["sgc"] if s["usage"] == usage and s["r1"] == r1 and s["r3"] == r3)
                res.update(usage=usage, n=3, w=list(w), r1=r1, r3=r3, restarts_cut=cut,
                           pos3_synonymy=pos3Synonymy(lab, nbr64, pos64, f, 10),
                           H_inv={"D2": d2_inv, "D3": D3, "A": 1.0 - D3 / d2_inv,
                                  "lift_D3": lift_inv_val, "lift_exact_absdiff": abs(lift_inv_val - d2_inv),
                                  "certified": D3 <= lift_inv_val + 1e-9 * max(1.0, lift_inv_val),
                                  "random_only_certified": res["random"]["best"] <= lift_inv_val + 1e-9 * max(1.0, lift_inv_val)},
                           H_term={"D2": d2_term, "D3": D3, "A": 1.0 - D3 / d2_term,
                                   "lift_D3_under_triplet_weights": lift_term_val},
                           sgc={"Dpos": sgc["Dpos"], "excess": sgc["Dpos"] - D3,
                                "excess_percent": 100.0 * (sgc["Dpos"] - D3) / D3})
                if not res["H_inv"]["certified"]:
                    res["SEARCH_FAILURE"] = True
                results["n3"].append(res)
                best_labs[f"n3_{usage}_w{w}"] = {c64[i]: CLASS_ORDER[lab[i]] for i in range(64)}
                prev3 = (w, lab)
                dump()
                wh = " ".join(f"{k}:{v['hits']}/{v['restarts']}" for k, v in res.get("warm", {}).items())
                print(f"n=3 {usage} w={w}: best={D3:.9f} random {res['random']['hits']}/{res['random']['restarts']} "
                      f"{wh} A_inv={res['H_inv']['A']:+.5f} A_term={res['H_term']['A']:+.5f} "
                      f"cert={res['H_inv']['certified']} syn3={res['pos3_synonymy']['raw']:.3f} "
                      f"({res['seconds']}s, t={time.time() - t_start:.0f}s)", flush=True)

    # ── secondary: uniform weights, 61 sense codons, U and P (only if time remains) ──
    W61 = weightMatrix(pos61, (1, 1, 1))
    for usage, f, d2 in (("U", fU, 20.017516292), ("P", fP, 18.835846661)):
        left = BUDGET_S - (time.time() - t_start)
        if left < 1.3 * results.get("_n3_sec", 60.0):
            results["secondary"].append({"usage": usage, "skipped": f"time: {left:.0f}s left"})
            dump()
            continue
        n2lab = best_labs[f"n2_{usage}_w{(1.0, 1.0)}"]
        lift = labelsFromMap(c61, {c: n2lab[c[:2]] for c in c61}, CLASS_ORDER)
        # the restricted lift may leave a class empty; warm starts must keep every class nonempty
        warm, names = [], []
        if np.bincount(lift, minlength=10).min() >= 1:
            warm, names = [lift], ["restricted_lift_n2"]
        job = dict(nbr=nbr61, W=W61, dist=chem, f=f, k=10, N=61, n_iter=2000 * 61 * 61,
                   seed0=seed0, warm_labs=warm, warm_names=names)
        seed0 += 100000
        res, lab = runSearch(job, N_RANDOM, N_WARM if warm else 0)
        res.update(usage=usage, n=3, stops_excluded=True, w=[1, 1, 1], D2_best_usage_cert=d2,
                   gap=d2 - res["best"], gap_percent=100.0 * (d2 - res["best"]) / d2,
                   usage_cert_R4=next(r["best_Dw"] for r in uc["rows"] if r["row"] == f"R4({'R2' if usage == 'U' else 'R3'})"))
        results["secondary"].append(res)
        best_labs[f"secondary_{usage}_61"] = {c61[i]: CLASS_ORDER[lab[i]] for i in range(61)}
        dump()
        print(f"secondary {usage}: best={res['best']:.9f} random {res['random']['hits']}/{res['random']['restarts']} "
              f"gap={res['gap']:.6f} ({res['seconds']}s, t={time.time() - t_start:.0f}s)", flush=True)
    results["finished"] = True
    dump()


if __name__ == "__main__":
    main()
