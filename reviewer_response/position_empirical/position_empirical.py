"""position-weighted optima at empirically measured misreading-rate vectors.

same objective, kernels, annealing, gates and certification as reviewer_response/position_model/
(imported read-only); only the weight vectors and the warm-start set differ.

triplet weights (w1, w2, w3) come from reviewer_response/position_rates/POSITION_RATES.md.
doublet hypotheses: H-inv = (w1, w2), H-term = (w1, w3).

silent-lift conjecture: D3*(w) = min_p D2*(w without p). the lift that silences position p
(label(c) = L2(c with p deleted)) reproduces D2 exactly, so D3* <= min_p D2* always holds; the
question is whether a non-silent triplet code can go strictly below it.

run from repo root:  python reviewer_response/position_empirical/position_empirical.py
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
# numba reads NUMBA_CACHE_DIR at import; set it before position_model overwrites os.environ,
# so every compiled kernel is cached here and nothing is written into position_model/
os.environ["NUMBA_CACHE_DIR"] = os.path.join(HERE, "numba_cache")
sys.dont_write_bytecode = True
import numba.core.config  # noqa: E402,F401

import json  # noqa: E402
import time  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

sys.path.insert(0, os.path.join(HERE, "..", "position_model"))
import position_model as PM  # noqa: E402  (read-only import)

os.environ["NUMBA_CACHE_DIR"] = os.path.join(HERE, "numba_cache")
assert numba.core.config.CACHE_DIR.startswith(HERE), numba.core.config.CACHE_DIR
UC, CLASS_ORDER, REPS = PM.UC, PM.CLASS_ORDER, PM.REPS

PM.N_WORKERS = 62
N2_RANDOM, N2_WARM = 256, 16
N3_RANDOM, N3_WARM = 256, 32
BUDGET_S = 24 * 60
TOL = 1e-9
WEIGHTS = [("L-Ec", (1.37, 1.0, 0.71), "Landerer 2024 E. coli"),
           ("L-Ec-LOO", (1.91, 1.0, 1.38), "Landerer E. coli, PXD014341 left out"),
           ("L-Sc", (1.40, 1.0, 0.71), "Landerer 2024 S. cerevisiae"),
           ("L-Sc-LOO", (1.16, 1.0, 0.57), "Landerer yeast, PXD026077 left out"),
           ("K-str", (0.98, 1.0, 1.39), "Katoh & Suga 2023 strict"),
           ("K-int", (0.99, 1.0, 1.79), "Katoh & Suga 2023 intensity-weighted"),
           ("UNIF", (1.0, 1.0, 1.0), "uniform")]
OUT_JSON = os.path.join(HERE, "RESULTS.json")
OUT_LAB = os.path.join(HERE, "best_labelings.json")
OUT_MD = os.path.join(HERE, "REPORT.md")


def silentLift(lab2, c16, c64, p):
    """triplet labeling with position p silent: label(c) = L2(c with position p deleted)."""
    m = dict(zip(c16, lab2))
    return np.array([m[c[:p] + c[p + 1:]] for c in c64], dtype=np.int64)


def dropPos(w, p):
    return tuple(float(x) for i, x in enumerate(w) if i != p)


def positionSynonymy(lab, nbr, pos, f, k, n):
    """per position: fraction of arcs changing that position that are synonymous (raw, usage-weighted)."""
    cnt = np.bincount(lab, minlength=k)
    out = []
    for p in range(n):
        m = pos == p
        rows, cols = np.nonzero(m)
        syn = lab[nbr[rows, cols]] == lab[rows]
        wt = f[lab[rows]] / cnt[lab[rows]]
        out.append({"raw": float(syn.mean()), "usage_weighted": float((syn * wt).sum() / wt.sum())})
    return out


def silencedPosition(syn):
    s = [p + 1 for p, d in enumerate(syn) if d["raw"] == 1.0]
    return s[0] if len(s) == 1 else ("none" if not s else s)


def writeReport(results):
    L = ["# Position-weighted optima at empirical misreading-rate vectors", "",
         "Objective, kernels, SA + steepest descent, and gates are imported unchanged from "
         "`reviewer_response/position_model/position_model.py`. k = 10 (RAA10 reps), class sizes free, "
         "no stop codons, 64 triplet / 16 doublet codons. Weights from `position_rates/POSITION_RATES.md`.", "",
         "H-inv doublet = (w1, w2); H-term doublet = (w1, w3). A = 1 - D3*/D2*. "
         "Conjecture: D3* = min_p D2*(w without p), tolerance 1e-9. "
         "Syn(p1/p2/p3) = raw fraction of arcs at each position that are synonymous in the best triplet code. "
         "SGC = 10-class projection of the standard code, 61 sense codons, same weights (not an optimum).", "",
         "| set | w | usage | D2*(H-inv) | D2*(H-term) | D2* silence p1 / p2 / p3 | D3* | silenced | syn p1/p2/p3 | A_inv | A_term | min_p D2* - D3* | conj. | cert. | SGC Dpos |",
         "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for r in results["n3"]:
        s = r["silent_D2"]
        syn = "/".join(f"{d['raw']:.3f}" for d in r["position_synonymy"])
        L.append(f"| {r['set']} | {tuple(r['w'])} | {r['usage']} | {r['D2_Hinv']:.6f} | {r['D2_Hterm']:.6f} | "
                 f"{s[0]:.6f} / {s[1]:.6f} / {s[2]:.6f} | {r['D3']:.6f} | {r['silenced_position']} | {syn} | "
                 f"{r['A_inv']:+.5f} | {r['A_term']:+.5f} | {r['conjecture']['gap']:.2e} | "
                 f"{'EQUAL' if r['conjecture']['equal'] else 'FALSE'} | {'yes' if r['certified'] else 'NO'} | "
                 f"{r['sgc_Dpos']:.6f} |")
    L += ["", f"Gates: A {results['gateA']['pass']}, B {results['gateB']['pass']}, "
          f"C (all three silent positions) {results['gateC']['pass']}. "
          f"Restarts: n=2 {N2_RANDOM} random (+{N2_WARM} warm from the previous n=2 optimum when one exists); "
          f"n=3 {N3_RANDOM} random + 33 warm (11 per silent lift). "
          f"Reductions: {results['reductions'] or 'none'}. Elapsed {results.get('elapsed_seconds')} s"
          f"{'' if results.get('finished') else ' (RUN NOT FINISHED)'}.", ""]
    open(OUT_MD + ".part", "w").write("\n".join(L))
    os.replace(OUT_MD + ".part", OUT_MD)


def main():
    t_start = time.time()
    root = PM.ROOT
    chem = PM.distMatrix(CLASS_ORDER, PM.loadClassVectors())
    comp, prot_meta = UC.loadProteome()
    assert not prot_meta["fallback_used"], prot_meta
    uc = json.load(open(PM.USAGE_CERT_JSON))

    tri = pd.read_csv(os.path.join(root, "reviewer_response", "codon_table_export.csv"))
    tri["codon"] = tri["codon"].str.upper()
    tri_map = {r.codon: PM.RAA10[r.aa.upper()] for r in tri.itertuples() if not r.is_stop}
    c61, nbr61 = PM.buildSenseGraph(3, PM.STOPS)
    c64, nbr64 = PM.buildSenseGraph(3)
    c16, nbr16 = PM.buildSenseGraph(2)
    pos61, pos64, pos16 = (PM.arcPositions(c, n) for c, n in ((c61, nbr61), (c64, nbr64), (c16, nbr16)))
    tri_lab = PM.labelsFromMap(c61, tri_map, CLASS_ORDER)
    fU = UC.uniformUsage(10, list(range(10)))
    fP = UC.classUsage(comp, set(CLASS_ORDER))
    assert np.allclose(fP, [uc["usage_P_classes"][c] for c in CLASS_ORDER], atol=1e-15, rtol=0)

    # gate A: w = 1 reproduces the usage_cert SGC values
    ref = uc["reference_codes"]["sgc_triplet_projection_stops_excluded"]
    W61_1 = PM.weightMatrix(pos61, (1, 1, 1))
    gA = {u: abs(PM.fullD(tri_lab, nbr61, W61_1, chem, f, 10) - ref[f"Dw_{u}"]) for u, f in (("U", fU), ("P", fP))}
    gateA = {"absdiff": gA, "pass": max(gA.values()) <= 1e-9}

    # gate B: kernel vs brute definition at the actual empirical weight vectors
    rng = np.random.default_rng(23)
    diffs = []
    for _, w, _ in WEIGHTS:
        for _ in range(5):
            f = rng.dirichlet(np.ones(10))
            for codons, nbr, pos, ww in ((c64, nbr64, pos64, w), (c61, nbr61, pos61, w), (c16, nbr16, pos16, w[:2])):
                lab = UC.randomStart(len(codons), 10, rng)
                d = PM.fullD(lab, nbr, PM.weightMatrix(pos, ww), chem, f, 10)
                diffs.append(abs(d - PM.bruteD(lab, codons, ww, chem, f, 10)) / max(1.0, d))
    gateB = {"n_cases": len(diffs), "max_reldiff": float(max(diffs)), "pass": bool(max(diffs) <= 1e-12)}

    # gate C: each of the three silent lifts is exact, Dpos3(lift_p(L2)) = Dpos2(L2; w without p);
    # and our p=3 lift equals usage_cert.liftLabels
    cd = []
    for _ in range(60):
        L2 = UC.randomStart(16, 10, rng)
        f = rng.dirichlet(np.ones(10))
        w = tuple(rng.uniform(0.3, 3.0, 3))
        for p in range(3):
            L3 = silentLift(L2, c16, c64, p)
            d2 = PM.fullD(L2, nbr16, PM.weightMatrix(pos16, dropPos(w, p)), chem, f, 10)
            d3 = PM.fullD(L3, nbr64, PM.weightMatrix(pos64, w), chem, f, 10)
            cd.append(abs(d2 - d3) / max(1.0, d2))
        assert np.array_equal(silentLift(L2, c16, c64, 2), UC.liftLabels(L2, c16, c64))
    gateC = {"n_cases": len(cd), "max_reldiff": float(max(cd)), "pass": max(cd) <= 1e-12}

    results = {"gateA": gateA, "gateB": gateB, "gateC": gateC, "proteome": prot_meta, "representatives": REPS,
               "weights": [{"set": s, "w": list(w), "source": src} for s, w, src in WEIGHTS],
               "config": {"n2": [N2_RANDOM, N2_WARM], "n3": [N3_RANDOM, N3_WARM], "workers": PM.N_WORKERS,
                          "n_iter_n2": "10000*N^2", "n_iter_n3": "2000*N^2", "tol": TOL},
               "sgc": [], "n2": [], "n3": [], "reductions": []}
    best_labs = {}
    print("gates:", gateA, gateB, gateC, flush=True)

    def dump():
        results["elapsed_seconds"] = round(time.time() - t_start, 1)
        json.dump(results, open(OUT_JSON + ".part", "w"), indent=1, default=lambda o: o.item())
        os.replace(OUT_JSON + ".part", OUT_JSON)
        json.dump(best_labs, open(OUT_LAB + ".part", "w"), indent=1)
        os.replace(OUT_LAB + ".part", OUT_LAB)

    if not (gateA["pass"] and gateB["pass"] and gateC["pass"]):
        results["aborted"] = "gate failure"
        dump()
        sys.exit("ABORT: gate failure")

    # SGC 10-class projection, stops excluded, under every triplet weight vector
    for usage, f in (("U", fU), ("P", fP)):
        for s, w, _ in WEIGHTS:
            results["sgc"].append({"set": s, "usage": usage, "w": list(w),
                                   "Dpos": PM.fullD(tri_lab, nbr61, PM.weightMatrix(pos61, w), chem, f, 10),
                                   "position_synonymy": positionSynonymy(tri_lab, nbr61, pos61, f, 10, 3)})
    dump()
    PM.anneal(UC.randomStart(16, 10, np.random.default_rng(0)), nbr16, PM.weightMatrix(pos16, (1, 1)),
              chem, fU, 10, 100, 1.0, 0)

    seed0 = 910000
    n2_pairs = []
    for _, w, _ in WEIGHTS:
        for p in range(3):
            if dropPos(w, p) not in n2_pairs:
                n2_pairs.append(dropPos(w, p))
    n3_times = []

    for usage, f in (("U", fU), ("P", fP)):
        best2 = {}
        prev = None
        for w2 in n2_pairs:
            job = dict(nbr=nbr16, W=PM.weightMatrix(pos16, w2), dist=chem, f=f, k=10, N=16, n_iter=10000 * 256,
                       seed0=seed0, warm_labs=[prev] if prev is not None else [],
                       warm_names=["prev_n2"] if prev is not None else [])
            seed0 += 100000
            res, lab = PM.runSearch(job, N2_RANDOM, N2_WARM if prev is not None else 0)
            res.update(usage=usage, n=2, w=list(w2))
            best2[w2] = {"D": res["best"], "lab": lab}
            prev = lab
            results["n2"].append(res)
            best_labs[f"n2_{usage}_w{w2}"] = {c16[i]: CLASS_ORDER[lab[i]] for i in range(16)}
            dump()
            print(f"n=2 {usage} w={w2}: best={res['best']:.9f} random {res['random']['hits']}/{N2_RANDOM} "
                  f"({res['seconds']}s, t={time.time() - t_start:.0f}s)", flush=True)

        for s, w, _ in WEIGHTS:
            W = PM.weightMatrix(pos64, w)
            lifts = [silentLift(best2[dropPos(w, p)]["lab"], c16, c64, p) for p in range(3)]
            nr, nw = N3_RANDOM, N3_WARM
            left = BUDGET_S - (time.time() - t_start)
            n3_left = 2 * len(WEIGHTS) - len(results["n3"])
            if n3_times and left < n3_left * max(n3_times):
                nr, nw = 128, 32
                results["reductions"].append(f"n3 {usage} {s}: {nr}+{nw}")
            job = dict(nbr=nbr64, W=W, dist=chem, f=f, k=10, N=64, n_iter=2000 * 64 * 64, seed0=seed0,
                       warm_labs=lifts, warm_names=[f"lift_silence_p{p + 1}" for p in range(3)])
            seed0 += 100000
            res, lab = PM.runSearch(job, nr, nw)
            n3_times.append(res["seconds"])
            D3 = res["best"]
            sil = [best2[dropPos(w, p)]["D"] for p in range(3)]
            lift_vals = [PM.fullD(L, nbr64, W, chem, f, 10) for L in lifts]
            mn = min(sil)
            syn = positionSynonymy(lab, nbr64, pos64, f, 10, 3)
            sgc = next(x for x in results["sgc"] if x["set"] == s and x["usage"] == usage)
            res.update(set=s, usage=usage, n=3, w=list(w), D3=D3, silent_D2=sil,
                       silent_lift_D3=lift_vals, lift_exact_maxdiff=max(abs(a - b) for a, b in zip(sil, lift_vals)),
                       D2_Hinv=best2[dropPos(w, 2)]["D"], D2_Hterm=best2[dropPos(w, 1)]["D"],
                       A_inv=1.0 - D3 / best2[dropPos(w, 2)]["D"], A_term=1.0 - D3 / best2[dropPos(w, 1)]["D"],
                       conjecture={"lhs_D3": D3, "rhs_min_p_D2": mn, "argmin_p": int(np.argmin(sil)) + 1,
                                   "gap": mn - D3, "equal": abs(mn - D3) <= TOL * max(1.0, abs(mn))},
                       certified=D3 <= min(lift_vals) + TOL * max(1.0, min(lift_vals)),
                       random_only_certified=res["random"]["best"] <= min(lift_vals) + TOL * max(1.0, min(lift_vals)),
                       position_synonymy=syn, silenced_position=silencedPosition(syn),
                       sgc_Dpos=sgc["Dpos"], sgc_excess_percent=100.0 * (sgc["Dpos"] - D3) / D3)
            if not res["certified"]:
                res["SEARCH_FAILURE"] = True
            results["n3"].append(res)
            best_labs[f"n3_{usage}_{s}"] = {c64[i]: CLASS_ORDER[lab[i]] for i in range(64)}
            dump()
            writeReport(results)
            wh = " ".join(f"{k}:{v['hits']}/{v['restarts']}" for k, v in res["warm"].items())
            print(f"n=3 {usage} {s} w={w}: D3={D3:.9f} minD2={mn:.9f} gap={mn - D3:.3e} random "
                  f"{res['random']['hits']}/{res['random']['restarts']} {wh} A_inv={res['A_inv']:+.5f} "
                  f"A_term={res['A_term']:+.5f} syn={[round(d['raw'], 3) for d in syn]} "
                  f"({res['seconds']}s, t={time.time() - t_start:.0f}s)", flush=True)

    results["finished"] = True
    dump()
    writeReport(results)


if __name__ == "__main__":
    main()
