"""usage-weighted doublet-vs-triplet optima with lift certification.

Dw = sum_a f_a * (1/n_a) * sum_{c in block(a)} sum_{c' in N(c)} dist(a, label(c'))
   = sum_a f_a * T_a / n_a,   T_a = arc sum out of block a

with f_a = n_a/N_sense this is the old per-residue D (gate 1). the silent lift
L3(xyz) = L2(xy) multiplies every n_a and every T_a by 4, so Dw(L3) = Dw(L2) for any f
(gate 2), and the best n=2 labeling lifted to n=3 is an upper bound for the n=3 optimum.

search: single-codon recolor simulated annealing (every class kept nonempty), geometric
cooling, then steepest descent over all single recolors AND all pairwise label swaps.

run from repo root:  python reviewer_response/usage_cert/usage_cert.py
"""
import html
import itertools
import json
import math
import os
import re
import sys
import time
from multiprocessing import get_context

import numpy as np
import pandas as pd
from numba import njit

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.dont_write_bytecode = True  # do not write .pyc into arch_cert/ or repo root
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "reviewer_response", "arch_cert"))
from run_publication_controls import RAA10, build_doublet_table  # noqa: E402
from arch_cert import (ALPHABET, STOPS, CLASS_ORDER, GATE_TRIPLET, GATE_DOUBLET,  # noqa: E402
                       loadClassVectors, buildSenseGraph, distMatrix, labelsFromMap)

N_WORKERS = 60
N_RANDOM = 512
N_RANDOM_FLOOR = 128
N_WARM = 64
BUDGET_S = 25 * 60
OUT_JSON = os.path.join(HERE, "RESULTS.json")
OUT_LAB = os.path.join(HERE, "best_labelings.json")
RELSTAT = os.path.join(HERE, "relstat.html")
RELSTAT_URL = "https://web.expasy.org/docs/relnotes/relstat.html"
# fallback composition (percent), used only if the release page cannot be parsed
FALLBACK_P = dict(A=8.25, R=5.53, N=4.06, D=5.45, C=1.37, Q=3.93, E=6.75, G=7.07, H=2.27, I=5.96,
                  L=9.66, K=5.84, M=2.42, F=3.86, P=4.70, S=6.56, T=5.34, W=1.08, Y=2.92, V=6.87)


# ── numba kernels ──

@njit(cache=True)
def arcSums(lab, nbr, dist, k):
    """T_a = sum over codons c of class a of sum over sense nbrs of dist(a, lab[nbr])."""
    T = np.zeros(k)
    for i in range(nbr.shape[0]):
        a = lab[i]
        for j in range(nbr.shape[1]):
            v = nbr[i, j]
            if v >= 0:
                T[a] += dist[a, lab[v]]
    return T


@njit(cache=True)
def evalDw(T, cnt, f):
    # empty classes contribute nothing (callers guarantee f_a = 0 there, or cnt_a >= 1)
    s = 0.0
    for a in range(T.shape[0]):
        if cnt[a] > 0:
            s += f[a] * T[a] / cnt[a]
    return s


@njit(cache=True)
def fullDw(lab, nbr, dist, f, k):
    cnt = np.zeros(k, dtype=np.int64)
    for i in range(lab.shape[0]):
        cnt[lab[i]] += 1
    return evalDw(arcSums(lab, nbr, dist, k), cnt, f)


@njit(cache=True)
def applyRecolor(lab, nbr, dist, T, cnt, i, b):
    """recolor codon i to class b, updating T and cnt in place (lab[i] too)."""
    a = lab[i]
    for j in range(nbr.shape[1]):
        v = nbr[i, j]
        if v >= 0:
            lv = lab[v]
            T[a] -= dist[a, lv]          # arc i->v leaves block a
            T[b] += dist[b, lv]          # and enters block b
            T[lv] += dist[lv, b] - dist[lv, a]  # arc v->i now points at class b
    cnt[a] -= 1
    cnt[b] += 1
    lab[i] = b


@njit(cache=True)
def recolorValue(lab, nbr, dist, f, T, cnt, i, b, sT, sc):
    """Dw after recoloring i to b, computed on scratch copies sT, sc (lab untouched)."""
    sT[:] = T
    sc[:] = cnt
    a = lab[i]
    for j in range(nbr.shape[1]):
        v = nbr[i, j]
        if v >= 0:
            lv = lab[v]
            sT[a] -= dist[a, lv]
            sT[b] += dist[b, lv]
            sT[lv] += dist[lv, b] - dist[lv, a]
    sc[a] -= 1
    sc[b] += 1
    return evalDw(sT, sc, f)


@njit(cache=True)
def swapValue(lab, nbr, dist, f, T, cnt, i, j, sT, sc):
    """Dw after swapping the labels of i and j (lab restored before return)."""
    sT[:] = T
    sc[:] = cnt
    li, lj = lab[i], lab[j]
    applyRecolor(lab, nbr, dist, sT, sc, i, lj)
    applyRecolor(lab, nbr, dist, sT, sc, j, li)
    val = evalDw(sT, sc, f)
    lab[i] = li
    lab[j] = lj
    return val


@njit(cache=True)
def descend(lab, nbr, dist, f, k):
    """steepest descent to a strict local optimum under all single recolors (keeping every
    class nonempty) and all pairwise label swaps. T is recomputed exactly after each move."""
    N = lab.shape[0]
    cnt = np.zeros(k, dtype=np.int64)
    for i in range(N):
        cnt[lab[i]] += 1
    sT = np.zeros(k)
    sc = np.zeros(k, dtype=np.int64)
    n_moves = 0
    while True:
        T = arcSums(lab, nbr, dist, k)
        cur = evalDw(T, cnt, f)
        best, kind, bi, bj = -1e-12, -1, -1, -1
        for i in range(N):
            if cnt[lab[i]] > 1:
                for b in range(k):
                    if b != lab[i]:
                        d = recolorValue(lab, nbr, dist, f, T, cnt, i, b, sT, sc) - cur
                        if d < best:
                            best, kind, bi, bj = d, 0, i, b
        for i in range(N):
            for j in range(i + 1, N):
                if lab[i] != lab[j]:
                    d = swapValue(lab, nbr, dist, f, T, cnt, i, j, sT, sc) - cur
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
def anneal(lab, nbr, dist, f, k, n_iter, t_scale, seed):
    """single-codon recolor SA; t_scale multiplies the automatic initial temperature
    (1 for random starts, < 1 for warm starts, 0 = skip SA and only descend)."""
    np.random.seed(seed)
    N = lab.shape[0]
    cnt = np.zeros(k, dtype=np.int64)
    for i in range(N):
        cnt[lab[i]] += 1
    sT = np.zeros(k)
    sc = np.zeros(k, dtype=np.int64)
    if t_scale > 0.0:
        T = arcSums(lab, nbr, dist, k)
        cur = evalDw(T, cnt, f)
        # initial temperature: mean |delta| over random legal recolors
        acc, m = 0.0, 0
        for _ in range(2000):
            i, b = np.random.randint(N), np.random.randint(k)
            if b != lab[i] and cnt[lab[i]] > 1:
                acc += abs(recolorValue(lab, nbr, dist, f, T, cnt, i, b, sT, sc) - cur)
                m += 1
        temp = t_scale * acc / max(m, 1) + 1e-12
        alpha = math.exp(math.log(1e-4) / n_iter)  # cool by 1e4 over the run
        for it in range(n_iter):
            i, b = np.random.randint(N), np.random.randint(k)
            if b != lab[i] and cnt[lab[i]] > 1:
                d = recolorValue(lab, nbr, dist, f, T, cnt, i, b, sT, sc) - cur
                if d <= 0.0 or np.random.random() < math.exp(-d / temp):
                    applyRecolor(lab, nbr, dist, T, cnt, i, b)
                    cur += d
            temp *= alpha
            # re-anchor the incremental sums periodically against float drift
            if (it & 0xFFFFF) == 0:
                T = arcSums(lab, nbr, dist, k)
                cur = evalDw(T, cnt, f)
    n_moves = descend(lab, nbr, dist, f, k)
    return n_moves


# ── worker ──

G = {}  # per-process job state, set before fork


def randomStart(N, k, rng):
    lab = np.empty(N, dtype=np.int64)
    perm = rng.permutation(N)
    lab[perm[:k]] = np.arange(k)
    lab[perm[k:]] = rng.integers(0, k, N - k)
    return lab


def runRestart(task):
    seed, warm_idx = task
    job = G["job"]
    rng = np.random.default_rng(seed)
    if warm_idx < 0:
        lab, t_scale = randomStart(job["N"], job["k"], rng), 1.0
    else:
        lab = job["warm_lab"].copy()
        # warm restart 0 is pure descent; the rest span initial temperatures 1e-3 .. 1e-1
        t_scale = 0.0 if warm_idx == 0 else 10.0 ** (-3.0 + 2.0 * (warm_idx - 1) / max(N_WARM - 2, 1))
    anneal(lab, job["nbr"], job["dist"], job["f"], job["k"], job["n_iter"], t_scale, seed)
    return fullDw(lab, job["nbr"], job["dist"], job["f"], job["k"]), lab.tolist()


def summarise(vals, best):
    v = np.array(vals)
    return {"restarts": len(v), "best": float(v.min()), "hits": int(np.sum(np.abs(v - best) <= 1e-9)),
            "p5": float(np.percentile(v, 5)), "p50": float(np.percentile(v, 50)), "max": float(v.max())}


def runSearch(job, n_random, n_warm):
    G["job"] = job
    t0 = time.time()
    tasks = [(job["seed0"] + r, -1) for r in range(n_random)]
    tasks += [(job["seed0"] + 50000 + r, r) for r in range(n_warm)]
    with get_context("fork").Pool(N_WORKERS) as pool:
        out = pool.map(runRestart, tasks, chunksize=1)
    vals = [d for d, _ in out]
    bi = int(np.argmin(vals))
    best, lab = float(vals[bi]), np.array(out[bi][1], dtype=np.int64)
    cnt = np.bincount(lab, minlength=job["k"])
    assert cnt.min() >= 1
    res = {"best_Dw": best, "random": summarise(vals[:n_random], best),
           "best_from": "random" if bi < n_random else "warm",
           "N_sense": job["N"], "n_iter": job["n_iter"],
           "class_sizes_best": {job["classes"][c]: int(cnt[c]) for c in range(job["k"])},
           "seconds": round(time.time() - t0, 1)}
    if n_warm:
        res["warm"] = summarise(vals[n_random:], best)
    return res, lab


# ── inputs ──

def loadProteome():
    """Swiss-Prot composition (percent) from the saved release page; fallback if unparsable."""
    try:
        if not os.path.exists(RELSTAT):
            import urllib.request
            urllib.request.urlretrieve(RELSTAT_URL, RELSTAT)
        s = html.unescape(re.sub(r"<[^>]+>", " ", open(RELSTAT, errors="replace").read()))
        s = re.sub(r"\s+", " ", s.replace("\xa0", " "))  # headers use irregular spacing
        rel = re.search(r"Release (\d{4}_\d{2}) of (\d{2}-\w{3}-\d{4})", s)
        sec = s[s.find("6.1 Composition"):s.find("6.2 Classification")]
        comp = {aa: float(x) for aa, x in re.findall(r"\(([A-Z])\)\s+([\d.]+)", sec)}
        comp = {aa: v for aa, v in comp.items() if aa in RAA10}
        assert len(comp) == 20 and rel
        return comp, {"source": RELSTAT_URL, "release": rel.group(1), "date": rel.group(2),
                      "fallback_used": False, "sum_percent": round(sum(comp.values()), 4)}
    except Exception as e:  # noqa: BLE001
        return dict(FALLBACK_P), {"source": "hardcoded fallback", "fallback_used": True,
                                  "error": repr(e), "sum_percent": round(sum(FALLBACK_P.values()), 4)}


def classUsage(comp, present):
    """class frequencies over CLASS_ORDER, zero outside `present`, renormalised."""
    f = np.zeros(len(CLASS_ORDER))
    for aa, cls in RAA10.items():
        if cls in present:
            f[CLASS_ORDER.index(cls)] += comp[aa]
    return f / f.sum()


def uniformUsage(k, present_idx):
    f = np.zeros(k)
    f[present_idx] = 1.0 / len(present_idx)
    return f


def liftLabels(lab2, codons2, codons3):
    """silent lift L3(xyz) = L2(xy) onto the given triplet codon list."""
    m = dict(zip(codons2, lab2))
    return np.array([m[c[:2]] for c in codons3], dtype=np.int64)


def main():
    t_start = time.time()
    vecs = loadClassVectors()
    chem = distMatrix(CLASS_ORDER, vecs)
    comp, prot_meta = loadProteome()
    print("proteome:", prot_meta, flush=True)

    tri = pd.read_csv(os.path.join(ROOT, "reviewer_response", "codon_table_export.csv"))
    tri["codon"] = tri["codon"].str.upper()
    assert set(tri.loc[tri["is_stop"].astype(bool), "codon"]) == set(STOPS)
    tri_map = {r.codon: RAA10[r.aa.upper()] for r in tri.itertuples() if not r.is_stop}
    c61, nbr61 = buildSenseGraph(3, STOPS)
    c64, nbr64 = buildSenseGraph(3)
    c16, nbr16 = buildSenseGraph(2)
    tri_lab = labelsFromMap(c61, tri_map, CLASS_ORDER)
    dbl = build_doublet_table(tri, RAA10)
    assert not dbl["is_stop"].any()
    dbl_lab = labelsFromMap(c16, dict(zip(dbl["codon"], dbl["aa"])), CLASS_ORDER)
    dbl_present = sorted(set(dbl_lab.tolist()))

    # ── gate 1: f = n_a / N reproduces the old D ──
    def countUsage(lab):
        return np.bincount(lab, minlength=10) / len(lab)
    D_tri = fullDw(tri_lab, nbr61, chem, countUsage(tri_lab), 10)
    D_dbl = fullDw(dbl_lab, nbr16, chem, countUsage(dbl_lab), 10)
    gate1 = {"triplet_Dw": D_tri, "triplet_expected": GATE_TRIPLET, "triplet_absdiff": abs(D_tri - GATE_TRIPLET),
             "doublet_Dw": D_dbl, "doublet_expected": GATE_DOUBLET, "doublet_absdiff": abs(D_dbl - GATE_DOUBLET)}
    gate1["pass"] = gate1["triplet_absdiff"] <= 1e-9 and gate1["doublet_absdiff"] <= 1e-9

    # ── gate 2: silent lift is exact for any f ──
    diffs = []
    lab3 = liftLabels(dbl_lab, c16, c64)
    f9 = uniformUsage(10, dbl_present)
    d2, d3 = fullDw(dbl_lab, nbr16, chem, f9, 10), fullDw(lab3, nbr64, chem, f9, 10)
    diffs.append(abs(d2 - d3))
    g2_proj = {"doublet_uniform9_n2": d2, "lift_n3": d3}
    rng = np.random.default_rng(7)
    for _ in range(100):
        L2 = rng.integers(0, 10, 16)
        present = sorted(set(L2.tolist()))
        f = np.zeros(10)
        f[present] = rng.dirichlet(np.ones(len(present)))
        L3 = liftLabels(L2, c16, c64)
        diffs.append(abs(fullDw(L2, nbr16, chem, f, 10) - fullDw(L3, nbr64, chem, f, 10)))
    gate2 = {"projection": g2_proj, "n_random_labelings": 100, "random_f": "Dirichlet(1) over present classes",
             "max_absdiff": float(max(diffs)), "tolerance": 1e-12}
    gate2["pass"] = gate2["max_absdiff"] <= 1e-12

    results = {"gate1": gate1, "gate2": gate2, "proteome": prot_meta, "proteome_percent": comp,
               "representatives": dict(zip(CLASS_ORDER, "ALWSNDRHPC")), "rows": []}
    print("gate1:", gate1, "\ngate2:", gate2, flush=True)
    if not (gate1["pass"] and gate2["pass"]):
        results["aborted"] = "gate failure"
        json.dump(results, open(OUT_JSON, "w"), indent=1)
        sys.exit("ABORT: gate failure")

    fP = classUsage(comp, set(CLASS_ORDER))
    fP9 = classUsage(comp, {CLASS_ORDER[i] for i in dbl_present})
    fU = uniformUsage(10, list(range(10)))
    results["usage_P_classes"] = dict(zip(CLASS_ORDER, fP.tolist()))
    results["usage_P_doublet9"] = {CLASS_ORDER[i]: fP9[i] for i in dbl_present}

    # ── reference codes ──
    results["reference_codes"] = {
        "doublet_projection": {"classes": [CLASS_ORDER[i] for i in dbl_present], "missing": ["X8"],
                               "note": "f renormalised over the 9 present classes",
                               "Dw_U": fullDw(dbl_lab, nbr16, chem, uniformUsage(10, dbl_present), 10),
                               "Dw_P": fullDw(dbl_lab, nbr16, chem, fP9, 10)},
        "sgc_triplet_projection_stops_excluded": {
            "Dw_U": fullDw(tri_lab, nbr61, chem, fU, 10), "Dw_P": fullDw(tri_lab, nbr61, chem, fP, 10)},
    }
    print("reference:", results["reference_codes"], flush=True)

    def dump():
        results["elapsed_seconds"] = round(time.time() - t_start, 1)
        json.dump(results, open(OUT_JSON, "w"), indent=1)
        json.dump(best_labs, open(OUT_LAB, "w"), indent=1)

    best_labs = {}
    graphs = {2: (c16, nbr16), 3: (c64, nbr64), "3s": (c61, nbr61)}
    scen = [("R1", k, 1.0 - np.eye(k), uniformUsage(k, list(range(k))), "U") for k in (4, 8, 10, 16)]
    scen += [("R2", 10, chem, fU, "U"), ("R3", 10, chem, fP, "P")]

    # warm the jit in the parent so forked workers inherit compiled code
    anneal(randomStart(16, 10, np.random.default_rng(0)), nbr16, chem, fU, 10, 100, 1.0, 0)

    best2 = {}
    seed0 = 1000
    order = [(s, 2) for s in scen] + [(s, 3) for s in scen] + [(s, "3s") for s in scen if s[0] != "R1"]
    for idx, ((name, k, dist, f, usage), gkey) in enumerate(order):
        codons, nbr = graphs[gkey]
        N = len(codons)
        n = 2 if gkey == 2 else 3
        classes = CLASS_ORDER if name != "R1" else [f"c{i}" for i in range(k)]
        chemistry = name != "R1"
        # schedule lengths set from a timing probe on the loaded machine (probe.py):
        # n=3 chemistry at 2000*N^2 is ~12 s/restart; hit rate roughly doubles 1000->4000*N^2
        n_iter = (1000 if not chemistry else (10000 if n == 2 else 2000)) * N * N
        job = dict(nbr=nbr, dist=dist, f=f, k=k, N=N, n_iter=n_iter, classes=classes,
                   seed0=seed0 + 100000 * idx, warm_lab=None)
        row_name = name if gkey != "3s" else f"R4({name})"
        n_random, n_warm = N_RANDOM, 0
        if n == 3 and chemistry:
            job["warm_lab"] = liftLabels(best2[(name, k)]["lab"], c16, codons)
            n_warm = N_WARM
        # budget guard: cut random restarts (never warm starts) if the clock is running down
        if n == 3 and time.time() - t_start > 0.55 * BUDGET_S:
            n_random = N_RANDOM_FLOOR
        res, lab = runSearch(job, n_random, n_warm)
        res.update(row=row_name, scenario=name, n=n, k=k, usage=usage,
                   stops_excluded=(gkey == "3s"), chemistry=chemistry,
                   random_restarts_cut=(n_random < N_RANDOM))
        if n == 2:
            # keyed by (scenario, k): the four R1 rows share a name
            best2[(name, k)] = {"Dw": res["best_Dw"], "lab": lab}
        else:
            lift_lab = liftLabels(best2[(name, k)]["lab"], c16, codons)
            # numba does not bounds-check: a label >= k silently corrupts memory
            assert 0 <= lift_lab.min() and lift_lab.max() < k
            lift_val = fullDw(lift_lab, nbr, dist, f, k)
            d2 = best2[(name, k)]["Dw"]
            G2 = d2 - res["best_Dw"]
            res.update(Dw_n2_best=d2, lift_Dw=lift_val, gap_G=G2, gap_G_percent=100.0 * G2 / d2)
            if gkey == 3:
                res["lift_exact_check"] = abs(lift_val - d2)
                res["certified"] = res["best_Dw"] <= lift_val + 1e-9
                res["random_only_certified"] = res["random"]["best"] <= lift_val + 1e-9
                if not res["certified"]:
                    res["SEARCH_FAILURE"] = True
                    res["n3_upper_bound"] = lift_val
            else:
                res["lift_note"] = "lift restricted to 61 sense codons is not exact; report only"
                res["below_restricted_lift"] = res["best_Dw"] <= lift_val + 1e-9
        if name == "R1" and k in (4, 16):
            res["analytic_3log4k"] = 3.0 * round(math.log(k, 4))
            res["reached_analytic"] = abs(res["best_Dw"] - res["analytic_3log4k"]) <= 1e-9
        best_labs[f"{row_name}_k{k}_n{n}"] = {codons[i]: classes[lab[i]] for i in range(N)}
        results["rows"].append(res)
        dump()
        extra = f" warm {res['warm']['hits']}/{res['warm']['restarts']}" if "warm" in res else ""
        print(f"{row_name} n={n} k={k} {usage}: best={res['best_Dw']:.9f} "
              f"random {res['random']['hits']}/{res['random']['restarts']}{extra} "
              f"p5={res['random']['p5']:.6f} p50={res['random']['p50']:.6f} "
              f"G={res.get('gap_G', float('nan')):.6f} ({res['seconds']}s, t={time.time() - t_start:.0f}s)",
              flush=True)


if __name__ == "__main__":
    main()
