"""triplet (n=3) vs quadruplet (n=4) Dpos optima at the full 20-residue alphabet.

Dpos = sum_a f_a (1/n_a) sum_{c in block(a)} sum_p w_p sum_{c' = c with position p changed} dist(a, label(c'))

k = 20 residues, each its own vector (pca1..pca8 z-scored across the 20, ddof=0); H(n,4) with all
4^n codons sense (no stops). usage U (1/20) or P (Swiss-Prot 2026_03, usage_cert's saved page).

quadruplet decoder hypotheses built from triplet weights (w1, w2, w3):
  Q-mon : (w1, w2, w2, w3)     Q-inv : (w1, w2, w3, w3)

the kernels (fullD, anneal, descend, bruteD, weightMatrix, arcPositions, runRestart) are imported
read-only from position_model; the graph builder from arch_cert; usage builders from usage_cert.

n=3 optima depend on the weight vector only up to a permutation of positions (a position permutation
is an automorphism of H(3,4)), so each n=3 search is run once on the sorted weights and the code for
any ordering is recovered by permuting codon positions (checked numerically every time).

run from repo root:  python reviewer_response/n4/n4.py
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
# numba caches and .pyc files stay inside n4/ (numba must see this before position_model resets it)
os.environ["NUMBA_CACHE_DIR"] = os.path.join(HERE, "numba_cache")
sys.dont_write_bytecode = True
import numba  # noqa: E402,F401

import itertools  # noqa: E402
import json  # noqa: E402
import time  # noqa: E402
from multiprocessing import get_context  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "reviewer_response", "position_model"))
import position_model as PM  # noqa: E402  (also puts arch_cert/usage_cert on sys.path)
os.environ["NUMBA_CACHE_DIR"] = os.path.join(HERE, "numba_cache")
from arch_cert import STOPS, buildSenseGraph  # noqa: E402
import usage_cert as UC  # noqa: E402

N_WORKERS = int(os.environ.get("N4_WORKERS", 56))
N3_RANDOM = 256
N4_RANDOM, N4_WARM = 128, 40
N4_RANDOM_FLOOR, N4_WARM_FLOOR = 96, 16
ITER3 = int(os.environ.get("N4_ITER3", 150))      # SA iterations per restart = ITER * N^2
ITER4 = int(os.environ.get("N4_ITER4", 15))
BUDGET_S = float(os.environ.get("N4_BUDGET", 25 * 60))
GATE_D = 23.922359743253192
OUT_JSON = os.path.join(HERE, "RESULTS.json")
OUT_LAB = os.path.join(HERE, "best_labelings.json")

WSETS = [("UNIF", (1.0, 1.0, 1.0)), ("K-str", (0.98, 1.0, 1.39)), ("K-int", (0.99, 1.0, 1.79)),
         ("L-Ec", (1.37, 1.0, 0.71)), ("w3=2", (1.0, 1.0, 2.0)), ("w3=5", (1.0, 1.0, 5.0)),
         ("w3=10", (1.0, 1.0, 10.0))]


def quadWeights(w, hyp):
    w1, w2, w3 = w
    return (w1, w2, w2, w3) if hyp == "Q-mon" else (w1, w2, w3, w3)


def loadAAVectors():
    """z-score pca1..pca8 across the 20 AAs (ddof=0), exactly as arch_cert.loadClassVectors."""
    df = pd.read_csv(os.path.join(ROOT, "reviewer_response", "aa_props_export.csv"))
    x = df[[f"pca{i}" for i in range(1, 9)]].to_numpy(float)
    x = (x - x.mean(axis=0)) / x.std(axis=0, ddof=0)
    aas = df["aa"].str.upper().tolist()
    assert len(aas) == 20 and len(set(aas)) == 20
    return aas, x


def insertLift(lab3, c3, c4, p):
    """lift a triplet code by inserting a silent base at 0-based position p of the quadruplet."""
    m = dict(zip(c3, lab3))
    return np.array([m[c[:p] + c[p + 1:]] for c in c4], dtype=np.int64)


def permuteCode(lab_sorted, c3, order):
    """code for weights w from the code for sorted weights w[order]: L(c) = Ls(c[order])."""
    m = dict(zip(c3, lab_sorted))
    return np.array([m["".join(c[order[q]] for q in range(3))] for c in c3], dtype=np.int64)


def posSynonymy(lab, nbr, pos, f, n):
    """per-position fraction of synonymous arcs: raw and usage-weighted (arc from class a weighted f_a/n_a)."""
    cnt = np.bincount(lab, minlength=len(f))
    out = []
    for p in range(n):
        m = pos == p
        src = np.repeat(np.arange(nbr.shape[0])[:, None], nbr.shape[1], axis=1)[m]
        dst = nbr[m]
        syn = lab[src] == lab[dst]
        wt = f[lab[src]] / cnt[lab[src]]
        out.append({"raw": float(syn.mean()), "usage_weighted": float((wt * syn).sum() / wt.sum())})
    return out


def runSearch(job, n_random, n_warm):
    """n_random random starts (full SA) + n_warm warm starts split evenly over job['warm_labs']."""
    PM.G["job"] = job
    t0 = time.time()
    n_src = len(job["warm_labs"])
    per = -(-n_warm // n_src) if n_src else 0
    job["warm_per_src"] = per
    tasks = [(job["seed0"] + r, -1, 0) for r in range(n_random)]
    tasks += [(job["seed0"] + 50000 + 1000 * s + r, s, r) for s in range(n_src) for r in range(per)]
    with get_context("fork").Pool(N_WORKERS) as pool:
        out = pool.map(PM.runRestart, tasks, chunksize=1)
    vals = [d for d, _ in out]
    bi = int(np.argmin(vals))
    best, lab = float(vals[bi]), np.array(out[bi][1], dtype=np.int64)
    assert np.bincount(lab, minlength=job["k"]).min() >= 1
    res = {"best": best, "random": PM.summarise(vals[:n_random], best),
           "best_from": "random" if bi < n_random else f"warm:{job['warm_names'][(bi - n_random) // per]}",
           "n_iter": job["n_iter"], "seconds": round(time.time() - t0, 1),
           "class_sizes_best": np.bincount(lab, minlength=job["k"]).tolist()}
    if per:
        res["warm"] = {job["warm_names"][s]: PM.summarise(vals[n_random + s * per:n_random + (s + 1) * per], best)
                       for s in range(n_src)}
    return res, lab


def main():
    t_start = time.time()
    aas, X = loadAAVectors()
    k = 20
    chem = np.array([[np.linalg.norm(X[a] - X[b]) for b in range(k)] for a in range(k)])
    blind = 1.0 - np.eye(k)
    comp, prot_meta = UC.loadProteome()
    assert not prot_meta["fallback_used"] and prot_meta["release"] == "2026_03", prot_meta
    fU = np.full(k, 1.0 / k)
    fP = np.array([comp[a] for a in aas])
    fP = fP / fP.sum()

    c3, nbr3 = buildSenseGraph(3)
    c4, nbr4 = buildSenseGraph(4)
    c61, nbr61 = buildSenseGraph(3, STOPS)
    pos3, pos4, pos61 = PM.arcPositions(c3, nbr3), PM.arcPositions(c4, nbr4), PM.arcPositions(c61, nbr61)
    rng = np.random.default_rng(20260923)

    # ── gate 1: SGC, stops removed, f_a = n_a/N, weights 1 -> published D ──
    tri = pd.read_csv(os.path.join(ROOT, "reviewer_response", "codon_table_export.csv"))
    tri["codon"] = tri["codon"].str.upper()
    sgc = {r.codon: r.aa.upper() for r in tri.itertuples() if not r.is_stop}
    sgc_lab = np.array([aas.index(sgc[c]) for c in c61], dtype=np.int64)
    fN = np.bincount(sgc_lab, minlength=k) / 61.0
    g1 = PM.fullD(sgc_lab, nbr61, PM.weightMatrix(pos61, (1, 1, 1)), chem, fN, k)
    g1b = PM.bruteD(sgc_lab, c61, (1, 1, 1), chem, fN, k)
    gate1 = {"Dpos": g1, "brute": g1b, "published": GATE_D, "absdiff": abs(g1 - GATE_D),
             "pass": abs(g1 - GATE_D) <= 1e-9 and abs(g1b - GATE_D) <= 1e-9}

    # ── gate 2: silent lift at position 4 is exact for any w4; also insertion at every position p ──
    d2 = []
    for _ in range(100):
        L3 = UC.randomStart(64, k, rng)
        f = rng.dirichlet(np.ones(k))
        w3 = rng.uniform(0.1, 10.0, 3)
        d3 = PM.fullD(L3, nbr3, PM.weightMatrix(pos3, w3), chem, f, k)
        for p in range(4):
            w4 = np.insert(w3, p, rng.uniform(0.1, 100.0))
            d4 = PM.fullD(insertLift(L3, c3, c4, p), nbr4, PM.weightMatrix(pos4, w4), chem, f, k)
            d2.append((p, abs(d4 - d3) / max(1.0, d3)))
    gate2 = {"n_labelings": 100, "max_reldiff_pos4": max(x for p, x in d2 if p == 3),
             "max_reldiff_any_insert_position": max(x for _, x in d2), "tolerance": 1e-12}
    gate2["pass"] = gate2["max_reldiff_any_insert_position"] <= 1e-12

    # ── gate 3: weighted kernel vs brute force at n = 3, 4 with k = 20; permuteCode invariance ──
    d3s, dperm = [], []
    for _ in range(10):
        for cod, nbr, pos, n in ((c3, nbr3, pos3, 3), (c4, nbr4, pos4, 4)):
            w = rng.uniform(0.1, 10.0, n)
            lab = UC.randomStart(len(cod), k, rng)
            f = rng.dirichlet(np.ones(k))
            d = PM.fullD(lab, nbr, PM.weightMatrix(pos, w), chem, f, k)
            d3s.append(abs(d - PM.bruteD(lab, cod, w, chem, f, k)) / max(1.0, d))
        w = rng.uniform(0.1, 10.0, 3)
        order = np.argsort(w, kind="stable")
        Ls = UC.randomStart(64, k, rng)
        f = rng.dirichlet(np.ones(k))
        a = PM.fullD(Ls, nbr3, PM.weightMatrix(pos3, w[order]), chem, f, k)
        b = PM.fullD(permuteCode(Ls, c3, order), nbr3, PM.weightMatrix(pos3, w), chem, f, k)
        dperm.append(abs(a - b) / max(1.0, a))
    gate3 = {"max_reldiff_brute": max(d3s), "max_reldiff_permutation": max(dperm), "tolerance": 1e-12}
    gate3["pass"] = max(d3s) <= 1e-12 and max(dperm) <= 1e-12

    results = {"gate1": gate1, "gate2": gate2, "gate3": gate3, "proteome": prot_meta, "aa_order": aas,
               "usage_U": fU.tolist(), "usage_P": fP.tolist(), "weight_sets": dict(WSETS),
               "config": {"N3_RANDOM": N3_RANDOM, "N4_RANDOM": N4_RANDOM, "N4_WARM": N4_WARM,
                          "floors_n4": [N4_RANDOM_FLOOR, N4_WARM_FLOOR], "workers": N_WORKERS,
                          "n_iter_n3": f"{ITER3}*N^2", "n_iter_n4": f"{ITER4}*N^2", "budget_s": BUDGET_S},
               "n3": [], "cells": [], "blind": [], "reductions": []}
    best_labs = {}

    def dump():
        results["elapsed_seconds"] = round(time.time() - t_start, 1)
        json.dump(results, open(OUT_JSON + ".part", "w"), indent=1, default=lambda o: o.item())  # numpy scalars
        os.replace(OUT_JSON + ".part", OUT_JSON)
        json.dump(best_labs, open(OUT_LAB + ".part", "w"))
        os.replace(OUT_LAB + ".part", OUT_LAB)

    print("gate1:", gate1, "\ngate2:", gate2, "\ngate3:", gate3, flush=True)
    if not (gate1["pass"] and gate2["pass"] and gate3["pass"]):
        results["aborted"] = "gate failure"
        dump()
        sys.exit("ABORT: gate failure")
    dump()

    # warm the jit before forking
    PM.anneal(UC.randomStart(64, k, np.random.default_rng(0)), nbr3, PM.weightMatrix(pos3, (1, 1, 1)),
              chem, fU, k, 100, 1.0, 0)
    seed = [100000]

    def nextSeed():
        seed[0] += 100000
        return seed[0]

    n3_cache = {}  # (usage tag, sorted weights) -> (D, sorted-weight code)

    def bestTriplet(tag, f, dist, w):
        """best n=3 code for ordered triplet weights w (search on sorted weights, cached)."""
        w = np.asarray(w, float)
        order = np.argsort(w, kind="stable")
        key = (tag, tuple(np.round(w[order], 12)))
        if key not in n3_cache:
            ws = w[order]
            job = dict(nbr=nbr3, W=PM.weightMatrix(pos3, ws), dist=dist, f=f, k=k, N=64,
                       n_iter=ITER3 * 64 * 64, seed0=nextSeed(), warm_labs=[], warm_names=[])
            res, lab = runSearch(job, N3_RANDOM, 0)
            res.update(usage=tag, n=3, w_sorted=ws.tolist(),
                       pos_synonymy=posSynonymy(lab, nbr3, pos3, f, 3))
            results["n3"].append(res)
            best_labs[f"n3_{tag}_{ws.tolist()}"] = {c3[i]: aas[lab[i]] for i in range(64)}
            n3_cache[key] = (res["best"], lab)
            dump()
            print(f"  n=3 {tag} w={ws.tolist()}: D3={res['best']:.9f} hits {res['random']['hits']}/"
                  f"{res['random']['restarts']} ({res['seconds']}s, t={time.time() - t_start:.0f}s)", flush=True)
        D, Ls = n3_cache[key]
        lab = permuteCode(Ls, c3, order)
        chk = PM.fullD(lab, nbr3, PM.weightMatrix(pos3, w), dist, f, k)
        assert abs(chk - D) <= 1e-9 * max(1.0, D), (chk, D)
        return D, lab

    def plan4(n4_left):
        """n=4 restart counts; cut to the floors if the last full cell time won't fit."""
        left = BUDGET_S - (time.time() - t_start)
        if left < n4_left * results.get("_n4_sec", 0.0):
            return N4_RANDOM_FLOOR, N4_WARM_FLOOR, True
        return N4_RANDOM, N4_WARM, False

    def quadCell(tag, f, dist, w3, hyp, wname, n4_left):
        w4 = quadWeights(w3, hyp)
        W4 = PM.weightMatrix(pos4, w4)
        D3, L3 = bestTriplet(tag, f, dist, w3)
        # warm sources: silent lift of the best triplet code at position 4, and for each quadruplet
        # position p the best triplet code for w4 with p deleted, lifted by inserting p back
        srcs = [("lift4_of_triplet_best", insertLift(L3, c3, c4, 3), None)]
        dele = []
        for p in range(4):
            wd = tuple(w4[:p] + w4[p + 1:])
            Dd, Ld = bestTriplet(tag, f, dist, wd)
            lift = insertLift(Ld, c3, c4, p)
            dele.append({"p": p + 1, "w_deleted": list(wd), "D3": Dd})
            srcs.append((f"del{p + 1}", lift, Dd))
        # every lift is evaluated under the quadruplet weights; deletion lifts must reproduce D3 exactly
        lifts = {}
        for name, lab, Dd in srcs:
            v = PM.fullD(lab, nbr4, W4, dist, f, k)
            lifts[name] = v
            if Dd is not None:
                assert abs(v - Dd) <= 1e-9 * max(1.0, Dd), (name, v, Dd)
        # dedupe identical warm labelings (the triplet lift can coincide with a deletion lift)
        warm, names, seen = [], [], set()
        for name, lab, _ in srcs:
            key = lab.tobytes()
            if key not in seen:
                seen.add(key)
                warm.append(lab)
                names.append(name)
        nr, nw, cut = plan4(n4_left)
        job = dict(nbr=nbr4, W=W4, dist=dist, f=f, k=k, N=256, n_iter=ITER4 * 256 * 256,
                   seed0=nextSeed(), warm_labs=warm, warm_names=names)
        res, lab = runSearch(job, nr, max(nw, len(warm)))
        if not cut:
            results["_n4_sec"] = res["seconds"]
        else:
            results["reductions"].append(f"n4 {tag} {wname} {hyp}: {nr}+{nw}")
        D4 = res["best"]
        min_lift = min(lifts.values())
        min_del = min(d["D3"] for d in dele)
        tol = 1e-9 * max(1.0, min_lift)
        cert = D4 <= min_lift + tol
        cell = dict(res, usage=tag, weight_set=wname, hypothesis=hyp, w3=list(w3), w4=list(w4),
                    restarts_cut=cut, D3=D3, D4_search=D4, lift_values=lifts, deletions=dele,
                    certified=cert, D4=D4 if cert else min_lift,
                    D4_equals_min_deletion=abs(min(D4, min_lift) - min_del) <= tol,
                    D4_below_min_deletion_by=min_del - min(D4, min_lift),
                    pos_synonymy_n4=posSynonymy(lab, nbr4, pos4, f, 4),
                    pos_synonymy_n3=posSynonymy(L3, nbr3, pos3, f, 3))
        if not cert:
            cell["SEARCH_FAILURE"] = True
            best_lift = min(lifts, key=lifts.get)
            cell["D4_upper_bound_from"] = best_lift
            lab = dict((n, l) for n, l, _ in srcs)[best_lift]
            cell["pos_synonymy_n4"] = posSynonymy(lab, nbr4, pos4, f, 4)
        cell["G"] = 1.0 - cell["D4"] / D3
        results["cells"].append(cell)
        best_labs[f"n4_{tag}_{wname}_{hyp}"] = "".join(aas[x] for x in lab)  # one letter per codon, codon order ACGU^4
        dump()
        print(f"n=4 {tag} {wname} {hyp} w4={w4}: D3={D3:.6f} D4={cell['D4']:.6f} G={cell['G']:+.6f} "
              f"cert={cert} =mindel={cell['D4_equals_min_deletion']} rand {res['random']['hits']}/"
              f"{res['random']['restarts']} best_from={res['best_from']} ({res['seconds']}s, "
              f"t={time.time() - t_start:.0f}s)", flush=True)

    n4_total = 2 * len(WSETS) * 2 + 1
    done4 = [0]
    for tag, f in (("U", fU), ("P", fP)):
        for wname, w3 in WSETS:
            for hyp in ("Q-mon", "Q-inv"):
                quadCell(tag, f, chem, w3, hyp, wname, n4_total - done4[0])
                done4[0] += 1
        if tag == "U":
            # chemistry-blind, U usage, uniform weights: n=3 and n=4 (one n=4 cell; Q-mon = Q-inv here)
            quadCell("U-blind", fU, blind, (1.0, 1.0, 1.0), "Q-mon", "UNIF", n4_total - done4[0])
            done4[0] += 1
            c = results["cells"][-1]
            results["blind"] = {"k": 20, "D3": c["D3"], "D4": c["D4"], "certified": c["certified"],
                                "lower_bound_3log4k": 3.0 * np.log(20) / np.log(4)}
            dump()
    results["finished"] = True
    dump()


if __name__ == "__main__":
    main()
