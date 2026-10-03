"""finish two pieces, importing (never modifying) usage_cert and position_model.

(1) R4 stop-codon variant: uniform position weights, k=10, RAA10 reps, free class sizes (>=1),
    usage U and P, n=3 on the 61-sense graph (UAA/UAG/UGA neither labelled nor neighbours).
    512 random restarts + 64 warm starts from the silent lift of the usage_cert best n=2
    labelings restricted to the 61 sense codons. kernels, schedule (2000*N^2) and the warm
    temperature ladder are usage_cert's own (UC.runSearch); only the seeds are new, so this is
    an independent replicate of usage_cert's R4 rows.

(2) position_model proteome (P) n=3 cells missing from its RESULTS.json, run with
    position_model's own runSearch, restart counts (128 random + max(32, 2*sources) warm),
    warm sources (lift H-inv, lift H-term, previous n=3 cell) and the seeds the original run
    would have used (7000 + cell_index*100000). the last completed cell is re-run first as a
    reproduction gate.

run from repo root:  python reviewer_response/finish/finish.py
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
os.environ["NUMBA_CACHE_DIR"] = os.path.join(HERE, "numba_cache")
sys.dont_write_bytecode = True
import numba  # noqa: E402,F401  (config read now, with our cache dir)

import json  # noqa: E402
import time  # noqa: E402

import numpy as np  # noqa: E402

ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
RR = os.path.join(ROOT, "reviewer_response")
sys.path.insert(0, ROOT)
for d in ("arch_cert", "usage_cert", "position_model"):
    sys.path.insert(0, os.path.join(RR, d))
import usage_cert as UC  # noqa: E402
import position_model as PM  # noqa: E402

# position_model overwrites NUMBA_CACHE_DIR on import; point caches back at finish/
os.environ["NUMBA_CACHE_DIR"] = os.path.join(HERE, "numba_cache")
numba.core.config.reload_config()

from arch_cert import STOPS, CLASS_ORDER, loadClassVectors, buildSenseGraph, distMatrix  # noqa: E402

OUT_JSON = os.path.join(HERE, "RESULTS.json")
OUT_LAB = os.path.join(HERE, "best_labelings.json")
D2_OPT = {"U": 20.017516292, "P": 18.835846661}
SEED_R4 = 9_000_000  # disjoint from usage_cert (1000 + 100000*idx, idx <= 13)


def positionSynonymy(lab, codons, nbr, f, k):
    """per codon position: fraction of sense->sense arcs whose endpoints share a class.
    raw counts every arc once; usage_weighted weights arc i->v by f_a / n_a (a = lab[i])."""
    cnt = np.bincount(lab, minlength=k)
    out = []
    for p in range(len(codons[0])):
        n_all = n_syn = wa = ws = 0.0
        for i in range(nbr.shape[0]):
            for v in nbr[i]:
                if v >= 0 and codons[i][p] != codons[v][p]:
                    a = lab[i]
                    syn = float(lab[v] == a)
                    n_all += 1
                    n_syn += syn
                    wa += f[a] / cnt[a]
                    ws += syn * f[a] / cnt[a]
        out.append({"position": p + 1, "raw": n_syn / n_all, "usage_weighted": ws / wa, "n_arcs": int(n_all)})
    return out


def main():
    t0 = time.time()
    chem = distMatrix(CLASS_ORDER, loadClassVectors())
    comp, prot = UC.loadProteome()
    assert not prot["fallback_used"], prot
    uc = json.load(open(os.path.join(RR, "usage_cert", "RESULTS.json")))
    uc_lab = json.load(open(os.path.join(RR, "usage_cert", "best_labelings.json")))
    pm = json.load(open(os.path.join(RR, "position_model", "RESULTS.json")))
    pm_lab = json.load(open(os.path.join(RR, "position_model", "best_labelings.json")))
    fU = UC.uniformUsage(10, list(range(10)))
    fP = UC.classUsage(comp, set(CLASS_ORDER))
    assert np.allclose(fP, [uc["usage_P_classes"][c] for c in CLASS_ORDER], atol=1e-15, rtol=0)
    assert np.allclose(fP, pm["usage_P"], atol=1e-15, rtol=0)

    c61, nbr61 = buildSenseGraph(3, STOPS)
    c64, nbr64 = buildSenseGraph(3)
    c16, nbr16 = buildSenseGraph(2)
    assert len(c61) == 61 and not set(STOPS) & set(c61)
    assert nbr61.max() < 61  # no arc points at a stop codon
    pos64 = PM.arcPositions(c64, nbr64)

    results = {"proteome": prot, "part1_R4": [], "part2_gate": None, "part2_cells": []}
    labs = {}

    def dump():
        results["elapsed_seconds"] = round(time.time() - t0, 1)
        for path, obj in ((OUT_JSON, results), (OUT_LAB, labs)):
            json.dump(obj, open(path + ".part", "w"), indent=1)
            os.replace(path + ".part", path)

    def toLab(m, codons):
        return np.array([CLASS_ORDER.index(m[c]) for c in codons], dtype=np.int64)

    # warm the jit before forking
    UC.anneal(UC.randomStart(16, 10, np.random.default_rng(0)), nbr16, chem, fU, 10, 100, 1.0, 0)
    PM.anneal(UC.randomStart(16, 10, np.random.default_rng(0)), nbr16, PM.weightMatrix(PM.arcPositions(c16, nbr16), (1, 1)),
              chem, fU, 10, 100, 1.0, 0)

    # ── part 1: R4 ──
    for idx, (usage, f, key, ucrow) in enumerate((("U", fU, "R2_k10_n2", "R4(R2)"), ("P", fP, "R3_k10_n2", "R4(R3)"))):
        lab2 = toLab(uc_lab[key], c16)
        d2 = UC.fullDw(lab2, nbr16, chem, f, 10)
        assert abs(d2 - D2_OPT[usage]) < 1e-8, (usage, d2)
        warm = UC.liftLabels(lab2, c16, c61)  # silent lift, stop codons simply dropped
        lift_val = UC.fullDw(warm, nbr61, chem, f, 10)
        job = dict(nbr=nbr61, dist=chem, f=f, k=10, N=61, n_iter=2000 * 61 * 61, classes=CLASS_ORDER,
                   seed0=SEED_R4 + 100000 * idx, warm_lab=warm)
        res, lab = UC.runSearch(job, UC.N_RANDOM, UC.N_WARM)
        prev = next(r for r in uc["rows"] if r["row"] == ucrow)
        best_all = min(res["best_Dw"], prev["best_Dw"])
        res.update(usage=usage, n=3, k=10, weights=[1, 1, 1], stops_excluded=True,
                   warm_lift_classes_present=int((np.bincount(warm, minlength=10) > 0).sum()),
                   warm_lift_Dw=lift_val, Dw_n2_opt=D2_OPT[usage],
                   gap=D2_OPT[usage] - res["best_Dw"], gap_percent=100 * (D2_OPT[usage] - res["best_Dw"]) / D2_OPT[usage],
                   usage_cert_R4_best=prev["best_Dw"], usage_cert_R4_hits=[prev["random"]["hits"], prev["warm"]["hits"]],
                   delta_vs_usage_cert=res["best_Dw"] - prev["best_Dw"], best_over_both_runs=best_all,
                   synonymy_by_position=positionSynonymy(lab, c61, nbr61, f, 10))
        results["part1_R4"].append(res)
        labs[f"R4_{usage}_61"] = {c61[i]: CLASS_ORDER[lab[i]] for i in range(61)}
        dump()
        syn = " ".join(f"p{s['position']}={s['raw']:.3f}" for s in res["synonymy_by_position"])
        print(f"R4 {usage}: best={res['best_Dw']:.9f} random {res['random']['hits']}/{res['random']['restarts']} "
              f"warm {res['warm']['hits']}/{res['warm']['restarts']} (usage_cert {prev['best_Dw']:.9f}) "
              f"gap={res['gap']:.6f} {syn} ({res['seconds']}s, t={time.time() - t0:.0f}s)", flush=True)

    # ── part 2: missing position_model P cells ──
    done = {(r["usage"], r["r1"], r["r3"]) for r in pm["n3"]}
    grid = [("P", r1, r3) for r1 in PM.R1S for r3 in PM.R3S]  # original order within usage P
    missing = [c for c in grid if c not in done]
    results["part2_missing"] = missing
    print("missing P n=3 cells:", missing, flush=True)
    best2 = {tuple(r["w"]): r["best"] for r in pm["n2"] if r["usage"] == "P"}

    def runCell(r1, r3, neighbour):
        w = (r1, 1.0, r3)
        W = PM.weightMatrix(pos64, w)
        j = grid.index(("P", r1, r3))
        seed0 = 7000 + (36 + j) * 100000  # U n2 (12), U n3 (12), P n2 (12) precede the P n3 cells
        lift_inv = UC.liftLabels(toLab(pm_lab[f"n2_P_w{(r1, 1.0)}"], c16), c16, c64)
        lift_term = UC.liftLabels(toLab(pm_lab[f"n2_P_w{(r1, r3)}"], c16), c16, c64)
        warm, names = [lift_inv, lift_term], ["lift_Hinv", "lift_Hterm"]
        if neighbour is not None:
            warm.append(neighbour[1])
            names.append(f"n3_neighbour{neighbour[0]}")
        job = dict(nbr=nbr64, W=W, dist=chem, f=fP, k=10, N=64, n_iter=2000 * 64 * 64,
                   seed0=seed0, warm_labs=warm, warm_names=names)
        res, lab = PM.runSearch(job, PM.N_RANDOM, max(PM.N_WARM, 2 * len(warm)))
        D3 = res["best"]
        d2_inv, d2_term = best2[(r1, 1.0)], best2[(r1, r3)]
        li, lt = PM.fullD(lift_inv, nbr64, W, chem, fP, 10), PM.fullD(lift_term, nbr64, W, chem, fP, 10)
        sgc = next(s for s in pm["sgc"] if s["usage"] == "P" and s["r1"] == r1 and s["r3"] == r3)
        tol = 1e-9 * max(1.0, li)
        res.update(usage="P", n=3, w=list(w), r1=r1, r3=r3, restarts_cut=False,
                   pos3_synonymy=PM.pos3Synonymy(lab, nbr64, pos64, fP, 10),
                   H_inv={"D2": d2_inv, "D3": D3, "A": 1.0 - D3 / d2_inv, "lift_D3": li,
                          "lift_exact_absdiff": abs(li - d2_inv), "certified": D3 <= li + tol,
                          "random_only_certified": res["random"]["best"] <= li + tol},
                   H_term={"D2": d2_term, "D3": D3, "A": 1.0 - D3 / d2_term, "lift_D3_under_triplet_weights": lt},
                   sgc={"Dpos": sgc["Dpos"], "excess": sgc["Dpos"] - D3,
                        "excess_percent": 100.0 * (sgc["Dpos"] - D3) / D3})
        if not res["H_inv"]["certified"]:
            res["SEARCH_FAILURE"] = True
        return res, lab

    # reproduction gate: re-run the last completed P cell (3, 1, 2) with its original neighbour (3, 1, 1)
    last = pm["n3"][-1]
    assert (last["usage"], last["r1"], last["r3"]) == ("P", 3.0, 2.0)
    nb = ((3.0, 1.0, 1.0), toLab(pm_lab["n3_P_w(3.0, 1.0, 1.0)"], c64))
    rg, lab_g = runCell(3.0, 2.0, nb)
    gate = {"cell": [3.0, 1.0, 2.0], "best_new": rg["best"], "best_orig": last["best"],
            "random_new": rg["random"], "random_orig": last["random"],
            "warm_new": rg["warm"], "warm_orig": last["warm"]}
    gate["pass"] = (abs(rg["best"] - last["best"]) <= 1e-9 * last["best"] and rg["random"] == last["random"]
                    and rg["warm"] == last["warm"])
    results["part2_gate"] = gate
    dump()
    print("reproduction gate:", gate["pass"], rg["best"], last["best"], f"t={time.time() - t0:.0f}s", flush=True)
    # continue the neighbour chain from the ORIGINAL (3,1,2) labeling, as the original run would have
    prev3 = ((3.0, 1.0, 2.0), toLab(pm_lab["n3_P_w(3.0, 1.0, 2.0)"], c64))
    for (_, r1, r3) in missing:
        res, lab = runCell(r1, r3, prev3)
        results["part2_cells"].append(res)
        labs[f"n3_P_w{(r1, 1.0, r3)}"] = {c64[i]: CLASS_ORDER[lab[i]] for i in range(64)}
        prev3 = ((r1, 1.0, r3), lab)
        dump()
        wh = " ".join(f"{k}:{v['hits']}/{v['restarts']}" for k, v in res["warm"].items())
        print(f"n=3 P w={res['w']}: best={res['best']:.9f} random {res['random']['hits']}/{res['random']['restarts']} "
              f"{wh} A_inv={res['H_inv']['A']:+.5f} A_term={res['H_term']['A']:+.5f} cert={res['H_inv']['certified']} "
              f"syn3={res['pos3_synonymy']['raw']:.3f} ({res['seconds']}s, t={time.time() - t0:.0f}s)", flush=True)
    results["finished"] = True
    dump()


if __name__ == "__main__":
    main()
