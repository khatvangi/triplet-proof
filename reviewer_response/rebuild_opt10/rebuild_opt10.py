"""ten-class optima (usage_cert R2/R3 and position_empirical grid) under the rebuilt descriptors.

kernels, gates, search procedures and seeds schedules are imported read-only from
reviewer_response/usage_cert/ (UC) and reviewer_response/position_model/ (PM, the kernel module
position_empirical imports); the only change is the 10x10 class distance matrix, built from
descriptor_rebuild/NEW_VECTORS.csv (pc1..pc8, each z-scored ddof=0) at reps A,L,W,S,N,D,R,H,P,C.

order (budget): gate, (c), (a) U, (b) U, (a) P, (b) P.

run from repo root:  python reviewer_response/rebuild_opt10/rebuild_opt10.py
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
# numba reads NUMBA_CACHE_DIR once at import; pin it here so nothing is written into usage_cert/ etc.
os.environ["NUMBA_CACHE_DIR"] = os.path.join(HERE, "numba_cache")
sys.dont_write_bytecode = True
import numba.core.config  # noqa: E402,F401

import json  # noqa: E402
import time  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

sys.path.insert(0, os.path.join(HERE, "..", "position_model"))
import position_model as PM  # noqa: E402  (read-only import; brings in usage_cert as PM.UC)

assert numba.core.config.CACHE_DIR.startswith(HERE), numba.core.config.CACHE_DIR
UC, CLASS_ORDER, REPS = PM.UC, PM.CLASS_ORDER, PM.REPS
RR = os.path.join(PM.ROOT, "reviewer_response")
NEW_VECTORS = os.path.join(RR, "descriptor_rebuild", "NEW_VECTORS.csv")

UC.N_WORKERS = PM.N_WORKERS = 62
T0 = time.time()
BUDGET_S = 21 * 60          # wall budget for the searches (task budget 28 min incl. setup/report)
TOL = 1e-9
# (a): usage_cert counts; (b): machine is shared (load ~88/64), so start at the stated floor
A_N2, A_N3, A_WARM, A_FLOOR = 512, 512, 64, 256
B_N2, B_N2_WARM, B_N3, B_N3_WARM = 128, 16, 128, 16
WEIGHTS = [("L-Ec", (1.37, 1.0, 0.71)), ("L-Ec-LOO", (1.91, 1.0, 1.38)), ("L-Sc", (1.40, 1.0, 0.71)),
           ("L-Sc-LOO", (1.16, 1.0, 0.57)), ("K-str", (0.98, 1.0, 1.39)), ("K-int", (0.99, 1.0, 1.79)),
           ("UNIF", (1.0, 1.0, 1.0))]
GATE_TARGET = 20.017516292  # usage_cert R2 U n=2 optimum
GATE_SEED0 = 1000 + 100000 * 4  # usage_cert seed0 for order index 4 = (R2, n=2)
OUT_JSON = os.path.join(HERE, "RESULTS.json")
OUT_LAB = os.path.join(HERE, "best_labelings.json")


def loadNewVectors():
    """z-score pc1..pc8 across the 20 AAs (ddof=0, as arch_cert does for the old pca columns)."""
    df = pd.read_csv(NEW_VECTORS)
    x = df[[f"pc{i}" for i in range(1, 9)]].to_numpy(float)
    x = (x - x.mean(axis=0)) / x.std(axis=0, ddof=0)
    aa_vec = {aa: x[i] for i, aa in enumerate(df["aa"].str.upper())}
    return {cls: aa_vec[aa] for cls, aa in REPS.items()}


def silentLift(lab2, c16, c64, p):
    m = dict(zip(c16, lab2))
    return np.array([m[c[:p] + c[p + 1:]] for c in c64], dtype=np.int64)


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


def left():
    return BUDGET_S - (time.time() - T0)


def main():
    comp, prot_meta = UC.loadProteome()
    assert not prot_meta["fallback_used"], prot_meta
    tri = pd.read_csv(os.path.join(RR, "codon_table_export.csv"))
    tri["codon"] = tri["codon"].str.upper()
    tri_map = {r.codon: PM.RAA10[r.aa.upper()] for r in tri.itertuples() if not r.is_stop}
    c61, nbr61 = PM.buildSenseGraph(3, PM.STOPS)
    c64, nbr64 = PM.buildSenseGraph(3)
    c16, nbr16 = PM.buildSenseGraph(2)
    pos64, pos16 = PM.arcPositions(c64, nbr64), PM.arcPositions(c16, nbr16)
    pos61 = PM.arcPositions(c61, nbr61)
    tri_lab = PM.labelsFromMap(c61, tri_map, CLASS_ORDER)
    dbl = PM.build_doublet_table(tri, PM.RAA10)
    dbl_lab = PM.labelsFromMap(c16, dict(zip(dbl["codon"], dbl["aa"])), CLASS_ORDER)
    dbl_present = sorted(set(dbl_lab.tolist()))
    fU = UC.uniformUsage(10, list(range(10)))
    fP = UC.classUsage(comp, set(CLASS_ORDER))
    fU9 = UC.uniformUsage(10, dbl_present)
    fP9 = UC.classUsage(comp, {CLASS_ORDER[i] for i in dbl_present})
    uc = json.load(open(os.path.join(RR, "usage_cert", "RESULTS.json")))
    assert np.allclose(fP, [uc["usage_P_classes"][c] for c in CLASS_ORDER], atol=1e-15, rtol=0)

    chem_old = PM.distMatrix(CLASS_ORDER, PM.loadClassVectors())
    chem = PM.distMatrix(CLASS_ORDER, loadNewVectors())
    results = {"vectors": os.path.relpath(NEW_VECTORS, PM.ROOT), "representatives": REPS,
               "zscore": "pc1..pc8, ddof=0 across 20 AAs", "proteome": prot_meta,
               "dist_new": chem.tolist(), "dist_old": chem_old.tolist(),
               "config": {"a": [A_N2, A_N3, A_WARM, A_FLOOR], "b_n2": [B_N2, B_N2_WARM],
                          "b_n3": [B_N3, B_N3_WARM], "workers": PM.N_WORKERS, "budget_s": BUDGET_S},
               "reductions": [], "a": [], "b_n2": [], "b_n3": []}
    best_labs = {}

    def dump():
        results["elapsed_seconds"] = round(time.time() - T0, 1)
        json.dump(results, open(OUT_JSON + ".part", "w"), indent=1, default=lambda o: o.item())
        os.replace(OUT_JSON + ".part", OUT_JSON)
        json.dump(best_labs, open(OUT_LAB + ".part", "w"), indent=1)
        os.replace(OUT_LAB + ".part", OUT_LAB)

    # ── GATE: old vectors, one usage_cert seed, R2 U n=2 must reach 20.017516292 within 1e-6 ──
    # (the SGC reference with the old vectors must also equal usage_cert's value to 1e-9)
    sgc_old_U = UC.fullDw(tri_lab, nbr61, chem_old, fU, 10)
    gate = {"sgc_old_U": sgc_old_U, "sgc_old_U_absdiff":
            abs(sgc_old_U - uc["reference_codes"]["sgc_triplet_projection_stops_excluded"]["Dw_U"]), "tries": []}
    for r in range(16):  # walk usage_cert's seed list in order until one seed hits (usage_cert hit rate 272/512)
        seed = GATE_SEED0 + r
        lab = UC.randomStart(16, 10, np.random.default_rng(seed))
        UC.anneal(lab, nbr16, chem_old, fU, 10, 10000 * 256, 1.0, seed)
        d = UC.fullDw(lab, nbr16, chem_old, fU, 10)
        gate["tries"].append({"seed": seed, "Dw": d})
        if abs(d - GATE_TARGET) <= 1e-6:
            break
    gate["seed_used"] = gate["tries"][-1]["seed"]
    gate["Dw"] = gate["tries"][-1]["Dw"]
    gate["absdiff"] = abs(gate["Dw"] - GATE_TARGET)
    gate["pass"] = bool(gate["absdiff"] <= 1e-6 and gate["sgc_old_U_absdiff"] <= 1e-9)
    results["gate"] = gate
    print("gate:", {k: v for k, v in gate.items() if k != "tries"}, f"tries={len(gate['tries'])}", flush=True)
    if not gate["pass"]:
        results["aborted"] = "gate failure"
        dump()
        sys.exit("ABORT: gate failure")

    # ── (c) reference codes under the new vectors, equal weights ──
    results["c"] = {"sgc_triplet_projection_stops_excluded": {"U": UC.fullDw(tri_lab, nbr61, chem, fU, 10),
                                                               "P": UC.fullDw(tri_lab, nbr61, chem, fP, 10)},
                    "doublet_projection": {"U": UC.fullDw(dbl_lab, nbr16, chem, fU9, 10),
                                           "P": UC.fullDw(dbl_lab, nbr16, chem, fP9, 10),
                                           "missing": ["X8"], "note": "f renormalised over the 9 present classes"},
                    "sgc_Dpos_by_weight": {u: {s: PM.fullD(tri_lab, nbr61, PM.weightMatrix(pos61, w), chem, f, 10)
                                               for s, w in WEIGHTS} for u, f in (("U", fU), ("P", fP))}}
    print("(c):", results["c"]["sgc_triplet_projection_stops_excluded"], results["c"]["doublet_projection"], flush=True)
    dump()

    # warm the jit in the parent so forked workers inherit compiled code
    UC.anneal(UC.randomStart(16, 10, np.random.default_rng(0)), nbr16, chem, fU, 10, 100, 1.0, 0)
    PM.anneal(UC.randomStart(16, 10, np.random.default_rng(0)), nbr16, PM.weightMatrix(pos16, (1, 1)),
              chem, fU, 10, 100, 1.0, 0)

    def runA(usage, f):
        """usage_cert R2 (U) / R3 (P): k=10, equal weights, free sizes, all 4^n codons sense."""
        name = "R2" if usage == "U" else "R3"
        seed0 = 7_000_000 + (0 if usage == "U" else 1_000_000)
        job2 = dict(nbr=nbr16, dist=chem, f=f, k=10, N=16, n_iter=10000 * 256, classes=CLASS_ORDER,
                    seed0=seed0, warm_lab=None)
        r2, lab2 = UC.runSearch(job2, A_N2, 0)
        n3 = A_N3
        if left() < 0.5 * BUDGET_S and usage == "P":
            n3 = A_FLOOR
            results["reductions"].append(f"(a) {usage} n=3 random {A_N3}->{n3}")
        lift = UC.liftLabels(lab2, c16, c64)
        job3 = dict(nbr=nbr64, dist=chem, f=f, k=10, N=64, n_iter=2000 * 64 * 64, classes=CLASS_ORDER,
                    seed0=seed0 + 100000, warm_lab=lift)
        r3, lab3 = UC.runSearch(job3, n3, A_WARM)
        lift_val = UC.fullDw(lift, nbr64, chem, f, 10)
        D2, D3 = r2["best_Dw"], r3["best_Dw"]
        row = {"scenario": name, "usage": usage, "D2": D2, "D3": D3, "lift_Dw": lift_val,
               "lift_exact_check": abs(lift_val - D2), "advantage": D2 - D3,
               "advantage_percent": 100.0 * (D2 - D3) / D2,
               "certified": bool(D3 <= lift_val + TOL), "random_only_certified":
               bool(r3["random"]["best"] <= lift_val + TOL), "n2": r2, "n3": r3,
               "n3_position_synonymy": positionSynonymy(lab3, nbr64, pos64, 3)}
        best_labs[f"a_{name}_n2"] = {c16[i]: CLASS_ORDER[lab2[i]] for i in range(16)}
        best_labs[f"a_{name}_n3"] = {c64[i]: CLASS_ORDER[lab3[i]] for i in range(64)}
        results["a"].append(row)
        sgc = results["c"]["sgc_triplet_projection_stops_excluded"][usage]
        results["c"][f"sgc_excess_over_n3_opt_{usage}"] = {"abs": sgc - D3, "percent": 100.0 * (sgc - D3) / D3}
        dump()
        print(f"(a) {name} {usage}: D2*={D2:.9f} ({r2['random']['hits']}/{r2['random']['restarts']}) "
              f"D3*={D3:.9f} (rand {r3['random']['hits']}/{r3['random']['restarts']}, warm "
              f"{r3['warm']['hits']}/{r3['warm']['restarts']}) lift={lift_val:.9f} cert={row['certified']} "
              f"t={time.time() - T0:.0f}s", flush=True)

    n2_pairs = []
    for _, w in WEIGHTS:
        for p in range(3):
            if dropPos(w, p) not in n2_pairs:
                n2_pairs.append(dropPos(w, p))
    n3_times = []

    def runB(usage, f):
        """position_empirical grid: n=2 on every dropped-position pair, n=3 at each weight set."""
        seed0 = 9_100_000 + (0 if usage == "U" else 5_000_000)
        best2, prev = {}, None
        for w2 in n2_pairs:
            job = dict(nbr=nbr16, W=PM.weightMatrix(pos16, w2), dist=chem, f=f, k=10, N=16, n_iter=10000 * 256,
                       seed0=seed0, warm_labs=[prev] if prev is not None else [],
                       warm_names=["prev_n2"] if prev is not None else [])
            seed0 += 100000
            res, lab = PM.runSearch(job, B_N2, B_N2_WARM if prev is not None else 0)
            res.update(usage=usage, n=2, w=list(w2))
            best2[w2] = {"D": res["best"], "lab": lab}
            prev = lab
            results["b_n2"].append(res)
            best_labs[f"b_n2_{usage}_w{w2}"] = {c16[i]: CLASS_ORDER[lab[i]] for i in range(16)}
            dump()
            print(f"(b) n=2 {usage} w={w2}: {res['best']:.9f} random {res['random']['hits']}/{B_N2} "
                  f"({res['seconds']}s, t={time.time() - T0:.0f}s)", flush=True)
        for s, w in WEIGHTS:
            W = PM.weightMatrix(pos64, w)
            lifts = [silentLift(best2[dropPos(w, p)]["lab"], c16, c64, p) for p in range(3)]
            job = dict(nbr=nbr64, W=W, dist=chem, f=f, k=10, N=64, n_iter=2000 * 64 * 64, seed0=seed0,
                       warm_labs=lifts, warm_names=[f"lift_silence_p{p + 1}" for p in range(3)])
            seed0 += 100000
            res, lab = PM.runSearch(job, B_N3, B_N3_WARM)
            n3_times.append(res["seconds"])
            D3 = res["best"]
            sil = [best2[dropPos(w, p)]["D"] for p in range(3)]
            lift_vals = [PM.fullD(L, nbr64, W, chem, f, 10) for L in lifts]
            mn = min(sil)
            syn = positionSynonymy(lab, nbr64, pos64, 3)
            Dinv, Dterm = best2[dropPos(w, 2)]["D"], best2[dropPos(w, 1)]["D"]
            res.update(set=s, usage=usage, n=3, w=list(w), D3=D3, silent_D2=sil, silent_lift_D3=lift_vals,
                       lift_exact_maxdiff=max(abs(a - b) for a, b in zip(sil, lift_vals)),
                       D2_Hinv=Dinv, D2_Hterm=Dterm, A_inv=1.0 - D3 / Dinv, A_term=1.0 - D3 / Dterm,
                       conjecture={"rhs_min_p_D2": mn, "argmin_p": int(np.argmin(sil)) + 1, "gap": mn - D3,
                                   "equal": bool(abs(mn - D3) <= TOL * max(1.0, abs(mn)))},
                       certified=bool(D3 <= min(lift_vals) + TOL * max(1.0, min(lift_vals))),
                       position_synonymy_raw=syn, silenced_position=silencedPosition(syn),
                       sgc_Dpos=results["c"]["sgc_Dpos_by_weight"][usage][s])
            results["b_n3"].append(res)
            best_labs[f"b_n3_{usage}_{s}"] = {c64[i]: CLASS_ORDER[lab[i]] for i in range(64)}
            dump()
            wh = " ".join(f"{k}:{v['hits']}/{v['restarts']}" for k, v in res["warm"].items())
            print(f"(b) n=3 {usage} {s}: D3={D3:.9f} minD2={mn:.9f} gap={mn - D3:.3e} rand "
                  f"{res['random']['hits']}/{res['random']['restarts']} {wh} A_inv={res['A_inv']:+.5f} "
                  f"A_term={res['A_term']:+.5f} syn={[round(x, 3) for x in syn]} "
                  f"({res['seconds']}s, t={time.time() - T0:.0f}s)", flush=True)

    runA("U", fU)
    runB("U", fU)
    runA("P", fP)
    runB("P", fP)
    results["finished"] = True
    dump()


if __name__ == "__main__":
    main()
