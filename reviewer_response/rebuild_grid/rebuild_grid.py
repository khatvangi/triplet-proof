"""ten-class position-weight grid (position_model) rerun with the rebuilt descriptors.

everything is imported read-only from reviewer_response/position_model/ (PM): kernels, runSearch,
restart counts (128 random + max(32, 2*sources) warm), warm sources (lift H-inv, lift H-term,
previous n=3 cell; previous n=2 cell for n=2), seed schedule (7000 + 100000*cell_index, same
cell order), certification (D3 <= H-inv lift value). the only change is the 10x10 distance
matrix, built from descriptor_rebuild/NEW_VECTORS.csv (pc1..pc8, z-scored ddof=0, reps as PM).

added (evaluation only, the searches are untouched):
  D3min2 = min_p D2*(w minus p)  (the best silent lift over all three dropped positions)
  silenced position = codon position(s) whose arcs are all synonymous in the best n=3 code
  (the doublet search is symmetric under swapping its two positions, so D2*(1, r3) needed for
   dropping p1 at r1 = 3 is the r1 = 1 H-term search at the same r3.)

gate: OLD vectors, U n=2 w=(1, 5) with its original seed0 and neighbour warm start must reproduce
position_model's best to 1e-6.

run from repo root:  python reviewer_response/rebuild_grid/rebuild_grid.py
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
os.environ["NUMBA_CACHE_DIR"] = os.path.join(HERE, "numba_cache")
sys.dont_write_bytecode = True
import numba.core.config  # noqa: E402

import json  # noqa: E402
import time  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

sys.path.insert(0, os.path.join(HERE, "..", "position_model"))
import position_model as PM  # noqa: E402

# position_model sets NUMBA_CACHE_DIR on import; point it back here
os.environ["NUMBA_CACHE_DIR"] = os.path.join(HERE, "numba_cache")
numba.core.config.reload_config()

UC, CLASS_ORDER, REPS = PM.UC, PM.CLASS_ORDER, PM.REPS
RR = os.path.join(PM.ROOT, "reviewer_response")
NEW_VECTORS = os.path.join(RR, "descriptor_rebuild", "NEW_VECTORS.csv")
PM.N_WORKERS = 62
T0 = time.time()
BUDGET_S = 19 * 60           # searches only; the rest of the 28 min is setup + report
FLOOR_R, FLOOR_W = 96, 16    # never below this
TOL = 1e-9
OUT_JSON = os.path.join(HERE, "RESULTS.json")
OUT_LAB = os.path.join(HERE, "best_labelings.json")


def loadNewVectors():
    """z-score pc1..pc8 across the 20 AAs (ddof=0, as arch_cert does for the old pca columns)."""
    df = pd.read_csv(NEW_VECTORS)
    x = df[[f"pc{i}" for i in range(1, 9)]].to_numpy(float)
    x = (x - x.mean(axis=0)) / x.std(axis=0, ddof=0)
    aa_vec = {aa: x[i] for i, aa in enumerate(df["aa"].str.upper())}
    return {cls: aa_vec[aa] for cls, aa in REPS.items()}


def silentLift(lab2, c16, c64, p, flip=False):
    """n=3 code whose label ignores position p (0-based): exact lift of the doublet code.
    flip: the doublet code was optimised for the reversed weight pair, so read its doublet backwards."""
    m = dict(zip(c16, lab2))
    key = (lambda d: d[::-1]) if flip else (lambda d: d)
    return np.array([m[key(c[:p] + c[p + 1:])] for c in c64], dtype=np.int64)


def dropPos(w, p):
    return tuple(float(x) for i, x in enumerate(w) if i != p)


def positionSynonymy(lab, nbr, pos, n):
    out = []
    for p in range(n):
        rows, cols = np.nonzero(pos == p)
        out.append(float((lab[nbr[rows, cols]] == lab[rows]).mean()))
    return out


def silencedPosition(syn):
    s = [p + 1 for p, v in enumerate(syn) if v == 1.0]
    return s[0] if len(s) == 1 else ("none" if not s else s)


def main():
    comp, prot_meta = UC.loadProteome()
    assert not prot_meta["fallback_used"], prot_meta
    tri = pd.read_csv(os.path.join(RR, "codon_table_export.csv"))
    tri["codon"] = tri["codon"].str.upper()
    tri_map = {r.codon: PM.RAA10[r.aa.upper()] for r in tri.itertuples() if not r.is_stop}
    c61, nbr61 = PM.buildSenseGraph(3, PM.STOPS)
    c64, nbr64 = PM.buildSenseGraph(3)
    c16, nbr16 = PM.buildSenseGraph(2)
    pos61, pos64, pos16 = PM.arcPositions(c61, nbr61), PM.arcPositions(c64, nbr64), PM.arcPositions(c16, nbr16)
    tri_lab = PM.labelsFromMap(c61, tri_map, CLASS_ORDER)
    fU = UC.uniformUsage(10, list(range(10)))
    fP = UC.classUsage(comp, set(CLASS_ORDER))
    pm = json.load(open(os.path.join(RR, "position_model", "RESULTS.json")))
    pm_lab = json.load(open(os.path.join(RR, "position_model", "best_labelings.json")))
    assert np.allclose(fP, pm["usage_P"], atol=1e-15, rtol=0)

    chem_old = PM.distMatrix(CLASS_ORDER, PM.loadClassVectors())
    chem = PM.distMatrix(CLASS_ORDER, loadNewVectors())
    results = {"vectors": os.path.relpath(NEW_VECTORS, PM.ROOT), "zscore": "pc1..pc8, ddof=0 across 20 AAs",
               "representatives": REPS, "proteome": prot_meta, "dist_new": chem.tolist(), "dist_old": chem_old.tolist(),
               "config": {"N_RANDOM": PM.N_RANDOM, "N_WARM": PM.N_WARM, "floor": [FLOOR_R, FLOOR_W],
                          "workers": PM.N_WORKERS, "budget_s": BUDGET_S, "seed_schedule": "7000 + 100000*cell_index"},
               "reductions": [], "sgc": [], "n2": [], "n3": []}
    best_labs = {}

    def dump():
        results["elapsed_seconds"] = round(time.time() - T0, 1)
        json.dump(results, open(OUT_JSON + ".part", "w"), indent=1, default=lambda o: o.item())
        os.replace(OUT_JSON + ".part", OUT_JSON)
        json.dump(best_labs, open(OUT_LAB + ".part", "w"), indent=1)
        os.replace(OUT_LAB + ".part", OUT_LAB)

    def toLab(m, codons):
        return np.array([CLASS_ORDER.index(m[c]) for c in codons], dtype=np.int64)

    # warm the jit in the parent so forked workers inherit compiled code
    PM.anneal(UC.randomStart(16, 10, np.random.default_rng(0)), nbr16, PM.weightMatrix(pos16, (1, 1)),
              chem, fU, 10, 100, 1.0, 0)

    # ── GATE: old vectors, U n=2 w=(1,5), cell index 2 of the original run ──
    ref = next(r for r in pm["n2"] if r["usage"] == "U" and r["w"] == [1.0, 5.0])
    job = dict(nbr=nbr16, W=PM.weightMatrix(pos16, (1.0, 5.0)), dist=chem_old, f=fU, k=10, N=16,
               n_iter=10000 * 256, seed0=7000 + 2 * 100000,
               warm_labs=[toLab(pm_lab["n2_U_w(1.0, 2.0)"], c16)], warm_names=["n2_neighbour"])
    g, _ = PM.runSearch(job, PM.N_RANDOM, PM.N_WARM)
    gate = {"cell": "U n=2 w=(1,5) (H-term of r1=1, r3=5)", "best_new": g["best"], "best_orig": ref["best"],
            "absdiff": abs(g["best"] - ref["best"]), "random_new": g["random"], "random_orig": ref["random"],
            "warm_new": g.get("warm"), "warm_orig": ref.get("warm"), "seconds": g["seconds"]}
    gate["pass"] = bool(gate["absdiff"] <= 1e-6)
    results["gate"] = gate
    dump()
    print("gate:", gate["pass"], g["best"], ref["best"], g["random"], ref["random"], f"({g['seconds']}s)", flush=True)
    if not gate["pass"]:
        results["aborted"] = "gate failure"
        dump()
        sys.exit("ABORT: gate failure")

    # ── SGC projection (61 sense codons) under every triplet weight, new vectors ──
    for usage, f in (("U", fU), ("P", fP)):
        for r1 in PM.R1S:
            for r3 in PM.R3S:
                lab = tri_lab
                results["sgc"].append({"usage": usage, "r1": r1, "r3": r3,
                                       "Dpos": PM.fullD(lab, nbr61, PM.weightMatrix(pos61, (r1, 1.0, r3)), chem, f, 10),
                                       "Dpos_old_vectors": PM.fullD(lab, nbr61, PM.weightMatrix(pos61, (r1, 1.0, r3)),
                                                                    chem_old, f, 10),
                                       "synonymy_by_position": positionSynonymy(lab, nbr61, pos61, 3)})
    dump()

    # ── restart planning: per-restart cost measured on the fly; cut to the floor if needed ──
    rate = {2: None, 3: None}  # seconds per restart
    n_left = {2: 24, 3: 24}

    def plan(n):
        left = BUDGET_S - (time.time() - T0)
        est = lambda m, nr, nw: (rate[m] or 0.0) * (nr + nw) * n_left[m]  # noqa: E731
        need = est(2, PM.N_RANDOM, PM.N_WARM) + est(3, PM.N_RANDOM, PM.N_WARM)
        if rate[3] is None and n == 2:
            return PM.N_RANDOM, PM.N_WARM, False
        return (PM.N_RANDOM, PM.N_WARM, False) if need <= left else (FLOOR_R, FLOOR_W, True)

    seed0 = 7000
    for usage, f in (("U", fU), ("P", fP)):
        best2 = {}
        for r1 in PM.R1S:
            prev = None
            for x in PM.R3S:
                w = (r1, x)
                warm, names = ([prev], ["n2_neighbour"]) if prev is not None else ([], [])
                job = dict(nbr=nbr16, W=PM.weightMatrix(pos16, w), dist=chem, f=f, k=10, N=16, n_iter=10000 * 256,
                           seed0=seed0, warm_labs=warm, warm_names=names)
                seed0 += 100000
                nr, nw, cut = plan(2)
                res, lab = PM.runSearch(job, nr, nw if warm else 0)
                rate[2] = res["seconds"] / (nr + (nw if warm else 0))
                n_left[2] -= 1
                res.update(usage=usage, n=2, w=list(w), restarts_cut=cut)
                if cut:
                    results["reductions"].append(f"n2 {usage} w={w}: {nr}+{nw}")
                best2[w] = {"D": res["best"], "lab": lab}
                prev = lab
                results["n2"].append(res)
                best_labs[f"n2_{usage}_w{w}"] = {c16[i]: CLASS_ORDER[lab[i]] for i in range(16)}
                dump()
                print(f"n=2 {usage} w={w}: best={res['best']:.9f} random {res['random']['hits']}/"
                      f"{res['random']['restarts']} ({res['seconds']}s, t={time.time() - T0:.0f}s)", flush=True)

        prev3 = None
        for r1 in PM.R1S:
            for r3 in PM.R3S:
                w = (r1, 1.0, r3)
                W = PM.weightMatrix(pos64, w)
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
                nw = max(nw, 2 * len(warm))
                res, lab = PM.runSearch(job, nr, nw)
                per = -(-nw // len(warm))
                rate[3] = res["seconds"] / (nr + per * len(warm))
                n_left[3] -= 1
                if cut:
                    results["reductions"].append(f"n3 {usage} w={w}: {nr}+{nw}")
                D3 = res["best"]
                d2_inv, d2_term = best2[(r1, 1.0)]["D"], best2[(r1, r3)]["D"]
                li = PM.fullD(lift_inv, nbr64, W, chem, f, 10)
                # min over dropped position p of D2*(w minus p); D2*(a, b) = D2*(b, a) by symmetry
                d2_drop, lift_drop = [], []
                for p in range(3):
                    wd = dropPos(w, p)
                    key = wd if wd in best2 else wd[::-1]
                    d2 = best2[key]["D"]
                    l2 = best2[key]["lab"]
                    lp = silentLift(l2, c16, c64, p, flip=key != wd)
                    d2_drop.append(d2)
                    lift_drop.append(PM.fullD(lp, nbr64, W, chem, f, 10))
                d3min2 = min(d2_drop)
                argmin = [p + 1 for p in range(3) if abs(d2_drop[p] - d3min2) <= TOL * max(1.0, d3min2)]
                syn = positionSynonymy(lab, nbr64, pos64, 3)
                sgc = next(s for s in results["sgc"] if s["usage"] == usage and s["r1"] == r1 and s["r3"] == r3)
                tol = TOL * max(1.0, li)
                res.update(usage=usage, n=3, w=list(w), r1=r1, r3=r3, restarts_cut=cut, restarts=[nr, nw],
                           synonymy_by_position=syn, silenced_position=silencedPosition(syn),
                           H_inv={"D2": d2_inv, "D3": D3, "A": 1.0 - D3 / d2_inv, "lift_D3": li,
                                  "lift_exact_absdiff": abs(li - d2_inv), "certified": bool(D3 <= li + tol),
                                  "random_only_certified": bool(res["random"]["best"] <= li + tol)},
                           H_term={"D2": d2_term, "D3": D3, "A": 1.0 - D3 / d2_term},
                           drop={"D2_drop_p1_p2_p3": d2_drop, "lift_value_under_w": lift_drop,
                                 "min_D2": d3min2, "argmin_positions": argmin,
                                 "D3_equals_min_D2": bool(abs(D3 - d3min2) <= TOL * max(1.0, d3min2)),
                                 "D3_minus_min_D2": D3 - d3min2,
                                 "certified_vs_min_D2": bool(D3 <= d3min2 + TOL * max(1.0, d3min2))},
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
                      f"minD2={d3min2:.9f} eq={res['drop']['D3_equals_min_D2']} sil={res['silenced_position']} "
                      f"cert={res['H_inv']['certified']} ({res['seconds']}s, t={time.time() - T0:.0f}s)", flush=True)
    results["finished"] = True
    dump()


if __name__ == "__main__":
    main()
