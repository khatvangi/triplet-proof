"""COPY of rebuild_opt20/run.py for usage P only (rebuild_grid2). edits, all marked "rebuild_grid2":
  plan = P cells only; seed counter resumes where opt20's U block ended (29 seeds used -> 12_900_000);
  n=4 warm restarts >= 40 (floor 96 random + 40 warm); gate B (22 s, old vectors) skipped because
  rebuild_opt20 already passed it with identical code; gate A kept. outputs go to rebuild_grid2/.

n=3 vs n=4 Dpos optima at k = 20 with the rebuilt descriptor matrix (NEW_VECTORS.csv).

reuses reviewer_response/n4/n4.py (imported read-only): quadWeights, insertLift, permuteCode,
posSynonymy, runSearch, and through it the position_model kernels / arch_cert graph / usage_cert usage.

differences from n4.main():
  - vectors: descriptor_rebuild/NEW_VECTORS.csv (8 PCs, already z-scored; asserted)
  - weight sets: UNIF, K-str, K-int, w3=2, w3=5, w3=10 (no L-Ec); chemistry-blind skipped
    (it is matrix-independent)
  - every n=3 search is two independent batches of 256 random restarts (seed blocks A and B);
    agreement to 1e-6 is recorded, and D3* = min of the two
  - gate: old vectors must reproduce n4's gate 1 and n4's first n=3 search (U, w = 1,1,1, same
    seed block 200000) before anything runs with the new vectors

run from anywhere:  python reviewer_response/rebuild_opt20/run.py
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
# numba cache and bytecode stay inside rebuild_opt20/: numba reads NUMBA_CACHE_DIR once, at its own
# import, so importing numba here pins the cache dir before n4 / position_model reassign the variable
os.environ["NUMBA_CACHE_DIR"] = os.path.join(HERE, "numba_cache")
sys.dont_write_bytecode = True
import numba  # noqa: E402

import json  # noqa: E402
import time  # noqa: E402
from multiprocessing import get_context  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
os.environ.setdefault("N4_WORKERS", "60")
sys.path.insert(0, os.path.join(ROOT, "reviewer_response", "n4"))
import n4  # noqa: E402  (read-only reuse)
PM, UC = n4.PM, n4.UC
os.environ["NUMBA_CACHE_DIR"] = os.path.join(HERE, "numba_cache")
assert numba.core.config.CACHE_DIR == os.path.join(HERE, "numba_cache"), numba.core.config.CACHE_DIR

N3_BATCH = 256
ITER3, ITER4 = 500, 150                         # as in the n4 run (config n_iter_n3/n4 = 500/150 * N^2)
N4_RANDOM, N4_WARM = n4.N4_RANDOM, n4.N4_WARM   # 128, 40
N4_RANDOM_FLOOR, N4_WARM_FLOOR = 96, 40  # rebuild_grid2: warm floor raised to 40 per task
BUDGET_S = float(os.environ.get("OPT20_BUDGET", 20 * 60))
SEED_B_OFFSET = 50_000_000                      # batch B seeds never overlap batch A or n4's blocks
AGREE_TOL = 1e-6
WSETS = [("UNIF", (1.0, 1.0, 1.0)), ("K-str", (0.98, 1.0, 1.39)), ("K-int", (0.99, 1.0, 1.79)),
         ("w3=2", (1.0, 1.0, 2.0)), ("w3=5", (1.0, 1.0, 5.0)), ("w3=10", (1.0, 1.0, 10.0))]
OUT_JSON = os.path.join(HERE, "opt20P_RESULTS.json")  # rebuild_grid2
OUT_LAB = os.path.join(HERE, "best_labelings_opt20P.json")  # rebuild_grid2


def loadNewVectors(aas):
    """NEW_VECTORS.csv rows reordered to the old aa order; already z-scored (ddof=0), checked."""
    df = pd.read_csv(os.path.join(ROOT, "reviewer_response", "descriptor_rebuild", "NEW_VECTORS.csv"))
    df = df.set_index(df["aa"].str.upper()).loc[aas]
    x = df[[f"pc{i}" for i in range(1, 9)]].to_numpy(float)
    assert np.abs(x.mean(axis=0)).max() < 1e-9 and np.abs(x.std(axis=0, ddof=0) - 1).max() < 1e-9
    return x


def distOf(X):
    return np.array([[np.linalg.norm(X[a] - X[b]) for b in range(len(X))] for a in range(len(X))])


def searchTwoBatches(job):
    """two independent batches of N3_BATCH random restarts in one pool map; returns per-batch stats."""
    PM.G["job"] = job
    job["warm_per_src"] = 0
    t0 = time.time()
    tasks = [(job["seed0"] + r, -1, 0) for r in range(N3_BATCH)]
    tasks += [(job["seed0"] + SEED_B_OFFSET + r, -1, 0) for r in range(N3_BATCH)]
    with get_context("fork").Pool(n4.N_WORKERS) as pool:
        out = pool.map(PM.runRestart, tasks, chunksize=1)
    vals = [d for d, _ in out]
    va, vb = vals[:N3_BATCH], vals[N3_BATCH:]
    ia, ib = int(np.argmin(va)), N3_BATCH + int(np.argmin(vb))
    bi = ia if vals[ia] <= vals[ib] else ib
    lab = np.array(out[bi][1], dtype=np.int64)
    best = float(vals[bi])
    res = {"best": best, "batchA": PM.summarise(va, best), "batchB": PM.summarise(vb, best),
           "batchA_best": float(vals[ia]), "batchB_best": float(vals[ib]),
           "batch_absdiff": abs(vals[ia] - vals[ib]), "batches_agree": abs(vals[ia] - vals[ib]) <= AGREE_TOL,
           "best_from": "A" if bi == ia else "B", "seed0_A": job["seed0"],
           "seed0_B": job["seed0"] + SEED_B_OFFSET, "n_iter": job["n_iter"],
           "seconds": round(time.time() - t0, 1)}
    assert np.bincount(lab, minlength=job["k"]).min() >= 1
    return res, lab


def main():
    t_start = time.time()
    aas, Xold = n4.loadAAVectors()
    k = 20
    c3, nbr3 = n4.buildSenseGraph(3)
    c4, nbr4 = n4.buildSenseGraph(4)
    c61, nbr61 = n4.buildSenseGraph(3, n4.STOPS)
    pos3, pos4, pos61 = PM.arcPositions(c3, nbr3), PM.arcPositions(c4, nbr4), PM.arcPositions(c61, nbr61)
    comp, prot_meta = UC.loadProteome()
    assert not prot_meta["fallback_used"] and prot_meta["release"] == "2026_03", prot_meta
    fU = np.full(k, 1.0 / k)
    fP = np.array([comp[a] for a in aas])
    fP = fP / fP.sum()

    results = {"proteome": prot_meta, "aa_order": aas, "weight_sets": dict(WSETS),
               "config": {"N3_BATCHES": [N3_BATCH, N3_BATCH], "N4_RANDOM": N4_RANDOM, "N4_WARM": N4_WARM,
                          "floors_n4": [N4_RANDOM_FLOOR, N4_WARM_FLOOR], "workers": n4.N_WORKERS,
                          "n_iter_n3": f"{ITER3}*N^2", "n_iter_n4": f"{ITER4}*N^2", "budget_s": BUDGET_S,
                          "seed_B_offset": SEED_B_OFFSET, "agree_tol": AGREE_TOL},
               "n3": [], "cells": [], "reductions": []}
    best_labs = {}

    def dump():
        results["elapsed_seconds"] = round(time.time() - t_start, 1)
        json.dump(results, open(OUT_JSON + ".part", "w"), indent=1, default=lambda o: o.item())
        os.replace(OUT_JSON + ".part", OUT_JSON)
        json.dump(best_labs, open(OUT_LAB + ".part", "w"))
        os.replace(OUT_LAB + ".part", OUT_LAB)

    # ── gate A: old vectors, SGC (61 sense codons, w = 1, f = n_a/N) reproduces n4 gate 1 ──
    tri = pd.read_csv(os.path.join(ROOT, "reviewer_response", "codon_table_export.csv"))
    tri["codon"] = tri["codon"].str.upper()
    sgc = {r.codon: r.aa.upper() for r in tri.itertuples() if not r.is_stop}
    sgc_lab = np.array([aas.index(sgc[c]) for c in c61], dtype=np.int64)
    fN = np.bincount(sgc_lab, minlength=k) / 61.0
    chem_old = distOf(Xold)
    g1 = PM.fullD(sgc_lab, nbr61, PM.weightMatrix(pos61, (1, 1, 1)), chem_old, fN, k)
    gateA = {"Dpos": g1, "target": n4.GATE_D, "absdiff": abs(g1 - n4.GATE_D), "pass": abs(g1 - n4.GATE_D) <= 1e-9}
    results["gateA_sgc_old"] = gateA
    print("gateA:", gateA, flush=True)

    # warm the jit before forking
    PM.anneal(UC.randomStart(64, k, np.random.default_rng(0)), nbr3, PM.weightMatrix(pos3, (1, 1, 1)),
              chem_old, fU, k, 100, 1.0, 0)

    # ── gate B: old vectors, n=3 U w=(1,1,1), n4's seed block 200000, 256 restarts ──
    n4R = json.load(open(os.path.join(ROOT, "reviewer_response", "n4", "RESULTS.json")))
    ref = next(x for x in n4R["n3"] if x["usage"] == "U" and x["w_sorted"] == [1.0, 1.0, 1.0])
    PM.G["job"] = job = dict(nbr=nbr3, W=PM.weightMatrix(pos3, (1.0, 1.0, 1.0)), dist=chem_old, f=fU, k=k, N=64,
                             n_iter=ITER3 * 64 * 64, seed0=200000, warm_labs=[], warm_names=[])
    t0 = time.time()
    op = json.load(open(os.path.join(ROOT, "reviewer_response", "rebuild_opt20", "RESULTS.json")))  # rebuild_grid2
    resB = {"best": op["gateB_n3_old"]["best"], "random": {"hits": op["gateB_n3_old"]["hits"], "p50": op["gateB_n3_old"]["p50"]}}  # rebuild_grid2: taken from opt20, not rerun
    gB = {"best": resB["best"], "n4_best": ref["best"], "n4_p50": ref["random"]["p50"],
          "absdiff": abs(resB["best"] - ref["best"]), "hits": resB["random"]["hits"],
          "n4_hits": ref["random"]["hits"], "p50": resB["random"]["p50"], "seconds": round(time.time() - t0, 1)}
    # n4's reported spread for this search is [best, p50]; identical seeds should reproduce best exactly
    gB["within_spread"] = ref["best"] - 1e-9 <= resB["best"] <= ref["random"]["p50"]
    gB["pass"] = gB["within_spread"]
    results["gateB_n3_old"] = gB
    print("gateB:", gB, flush=True)
    if not (gateA["pass"] and gB["pass"]):
        results["aborted"] = "gate failure"
        dump()
        sys.exit("ABORT: gate failure")

    # ── new vectors ──
    Xnew = loadNewVectors(aas)
    chem = distOf(Xnew)
    results["sgc_D_new_info"] = PM.fullD(sgc_lab, nbr61, PM.weightMatrix(pos61, (1, 1, 1)), chem, fN, k)
    dump()

    seed = [12_900_000]  # rebuild_grid2: opt20 U block used 29 seeds from 10_000_000

    def nextSeed():
        seed[0] += 100_000
        return seed[0]

    n3_cache = {}
    n3_sec = [33.0]  # rebuild_grid2: opt20 measured ~30-39 s          # estimate for a two-batch search until measured
    n4_sec = [35.0]

    def bestTriplet(tag, f, w):
        w = np.asarray(w, float)
        order = np.argsort(w, kind="stable")
        key = (tag, tuple(np.round(w[order], 12)))
        if key not in n3_cache:
            ws = w[order]
            job = dict(nbr=nbr3, W=PM.weightMatrix(pos3, ws), dist=chem, f=f, k=k, N=64,
                       n_iter=ITER3 * 64 * 64, seed0=nextSeed(), warm_labs=[], warm_names=[])
            res, lab = searchTwoBatches(job)
            n3_sec[0] = res["seconds"]
            res.update(usage=tag, n=3, w_sorted=ws.tolist(), pos_synonymy=n4.posSynonymy(lab, nbr3, pos3, f, 3))
            results["n3"].append(res)
            best_labs[f"n3_{tag}_{ws.tolist()}"] = {c3[i]: aas[lab[i]] for i in range(64)}
            n3_cache[key] = (res["best"], lab)
            dump()
            print(f"  n=3 {tag} w={ws.tolist()}: A={res['batchA_best']:.9f} B={res['batchB_best']:.9f} "
                  f"agree={res['batches_agree']} hits A {res['batchA']['hits']} B {res['batchB']['hits']} "
                  f"({res['seconds']}s, t={time.time() - t_start:.0f}s)", flush=True)
        D, Ls = n3_cache[key]
        lab = n4.permuteCode(Ls, c3, order)
        chk = PM.fullD(lab, nbr3, PM.weightMatrix(pos3, w), chem, f, k)
        assert abs(chk - D) <= 1e-9 * max(1.0, D), (chk, D)
        return D, lab

    # n=3 searches still to run, per remaining cell, to budget the n=4 restarts
    def n3Pending(cells_left):
        keys = set()
        for tag, w3, hyp in cells_left:
            w4 = n4.quadWeights(w3, hyp)
            for wd in [w3] + [tuple(w4[:p] + w4[p + 1:]) for p in range(4)]:
                keys.add((tag, tuple(np.round(np.sort(np.asarray(wd, float)), 12))))
        return len(keys - set(n3_cache))

    def plan4(cells_left):
        left = BUDGET_S - (time.time() - t_start)
        need = len(cells_left) * n4_sec[0] + n3Pending(cells_left) * n3_sec[0]
        if left < need:
            return N4_RANDOM_FLOOR, N4_WARM_FLOOR, True
        return N4_RANDOM, N4_WARM, False

    def quadCell(tag, f, w3, hyp, wname, cells_left):
        w4 = n4.quadWeights(w3, hyp)
        W4 = PM.weightMatrix(pos4, w4)
        D3, L3 = bestTriplet(tag, f, w3)
        srcs = [("lift4_of_triplet_best", n4.insertLift(L3, c3, c4, 3), None)]
        dele = []
        for p in range(4):
            wd = tuple(w4[:p] + w4[p + 1:])
            Dd, Ld = bestTriplet(tag, f, wd)
            dele.append({"p": p + 1, "w_deleted": list(wd), "D3": Dd})
            srcs.append((f"del{p + 1}", n4.insertLift(Ld, c3, c4, p), Dd))
        lifts = {}
        for name, lab, Dd in srcs:
            v = PM.fullD(lab, nbr4, W4, chem, f, k)
            lifts[name] = v
            if Dd is not None:
                assert abs(v - Dd) <= 1e-9 * max(1.0, Dd), (name, v, Dd)
        warm, names, seen = [], [], set()
        for name, lab, _ in srcs:
            if lab.tobytes() not in seen:
                seen.add(lab.tobytes())
                warm.append(lab)
                names.append(name)
        nr, nw, cut = plan4(cells_left)
        job = dict(nbr=nbr4, W=W4, dist=chem, f=f, k=k, N=256, n_iter=ITER4 * 256 * 256,
                   seed0=nextSeed(), warm_labs=warm, warm_names=names)
        res, lab = n4.runSearch(job, nr, max(nw, len(warm)))
        if not cut:
            n4_sec[0] = res["seconds"]
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
                    pos_synonymy_n4=n4.posSynonymy(lab, nbr4, pos4, f, 4),
                    pos_synonymy_n3=n4.posSynonymy(L3, nbr3, pos3, f, 3))
        if not cert:
            cell["SEARCH_FAILURE"] = True
            best_lift = min(lifts, key=lifts.get)
            cell["D4_upper_bound_from"] = best_lift
            lab = dict((n, l) for n, l, _ in srcs)[best_lift]
            cell["pos_synonymy_n4"] = n4.posSynonymy(lab, nbr4, pos4, f, 4)
        cell["G"] = 1.0 - cell["D4"] / D3
        results["cells"].append(cell)
        best_labs[f"n4_{tag}_{wname}_{hyp}"] = "".join(aas[x] for x in lab)
        dump()
        print(f"n=4 {tag} {wname} {hyp} w4={w4}: D3={D3:.6f} D4={cell['D4']:.6f} G={cell['G']:+.6f} "
              f"cert={cert} =mindel={cell['D4_equals_min_deletion']} rand {res['random']['hits']}/"
              f"{res['random']['restarts']} best_from={res['best_from']} ({res['seconds']}s, "
              f"t={time.time() - t_start:.0f}s)", flush=True)

    # U, all cells, then P
    plan = [(tag, f, w3, hyp, wname) for tag, f in (("P", fP),)  # rebuild_grid2: P only
            for wname, w3 in WSETS for hyp in ("Q-mon", "Q-inv")]
    for i, (tag, f, w3, hyp, wname) in enumerate(plan):
        quadCell(tag, f, w3, hyp, wname, [(t, w, h) for t, _, w, h, _ in plan[i:]])
    results["finished"] = True
    dump()


if __name__ == "__main__":
    main()
