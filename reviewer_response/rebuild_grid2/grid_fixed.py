"""(A) rebuild_grid n=3 cells rerun with the fixed certificate.

reuses rebuild_grid/rebuild_grid.py read-only (imported: loadNewVectors, silentLift, dropPos,
positionSynonymy, silencedPosition) and through it position_model (kernels, runSearch, graph).
n=2 optima are taken from rebuild_grid's finished run (same script, same seeds); each doublet
label is re-scored here and must reproduce its stored D2* to 1e-9.

changes vs rebuild_grid:
  - warm starts for n=3 = the three silent lifts (drop p1, drop p2, drop p3), each built with
    silentLift at the dropped position, plus the previous n=3 cell (same neighbour rule).
    rebuild_grid's 'lift_Hterm' was UC.liftLabels (silent at p3, not p2), so it was never the
    H-term lift; the drop-p1 lift was never a warm start at all.
  - certificate: D3 <= min over the three lift values (not just the H-inv lift).
random restarts use rebuild_grid's seed schedule (7000 + 100000*cell_index), so the random part
reproduces rebuild_grid's random restarts exactly when the restart count is the same.

order: U r1=3 (all six r3, including the two ABOVE cells), then P (all 12). cells are passed on argv.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.dont_write_bytecode = True
sys.path.insert(0, os.path.join(HERE, "..", "rebuild_grid"))
import rebuild_grid as RG  # noqa: E402  (read-only reuse; module-level code only sets paths)
import numba.core.config  # noqa: E402
os.environ["NUMBA_CACHE_DIR"] = os.path.join(HERE, "numba_cache")
numba.core.config.reload_config()

import json  # noqa: E402
import time  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

PM, UC, CLASS_ORDER = RG.PM, RG.UC, RG.CLASS_ORDER
PM.N_WORKERS = 62
TOL = 1e-9
OUT = os.path.join(HERE, "RESULTS.json")
OUT_LAB = os.path.join(HERE, "best_labelings_grid.json")
N_RANDOM = int(os.environ.get("G_RANDOM", 128))
N_WARM = int(os.environ.get("G_WARM", 32))
assert N_RANDOM >= 96 and N_WARM >= 16


def cellIndex(usage, r1, r3):
    """rebuild_grid's cell counter: U n2 (12), U n3 (12), P n2 (12), P n3 (12)."""
    j = PM.R1S.index(r1) * len(PM.R3S) + PM.R3S.index(r3)
    return (12 + j) if usage == "U" else (36 + j)


def updateResults(key, value):
    r = json.load(open(OUT)) if os.path.exists(OUT) else {}
    r[key] = value
    json.dump(r, open(OUT + ".part", "w"), indent=1, default=lambda o: o.item() if hasattr(o, "item") else str(o))
    os.replace(OUT + ".part", OUT)


def main(cells):
    T0 = time.time()
    comp, prot_meta = UC.loadProteome()
    assert not prot_meta["fallback_used"], prot_meta
    c64, nbr64 = PM.buildSenseGraph(3)
    c16, nbr16 = PM.buildSenseGraph(2)
    pos64, pos16 = PM.arcPositions(c64, nbr64), PM.arcPositions(c16, nbr16)
    fU = UC.uniformUsage(10, list(range(10)))
    fP = UC.classUsage(comp, set(CLASS_ORDER))
    chem = PM.distMatrix(CLASS_ORDER, RG.loadNewVectors())
    rg = json.load(open(os.path.join(HERE, "..", "rebuild_grid", "RESULTS.json")))
    rg_lab = json.load(open(os.path.join(HERE, "..", "rebuild_grid", "best_labelings.json")))
    assert rg.get("finished") and not rg["reductions"], "rebuild_grid run must be finished and uncut"
    assert np.allclose(np.array(rg["dist_new"]), chem, atol=1e-12, rtol=0)

    # doublet optima from rebuild_grid, re-scored
    best2 = {}
    for x in rg["n2"]:
        u, w = x["usage"], tuple(x["w"])
        f = fU if u == "U" else fP
        m = rg_lab[f"n2_{u}_w{w}"]
        lab = np.array([CLASS_ORDER.index(m[c]) for c in c16], dtype=np.int64)
        d = PM.fullD(lab, nbr16, PM.weightMatrix(pos16, w), chem, f, 10)
        assert abs(d - x["best"]) <= 1e-9 * max(1.0, d), (u, w, d, x["best"])
        best2[(u, w)] = {"D": x["best"], "lab": lab}
    rg3 = {(x["usage"], x["r1"], x["r3"]): x for x in rg["n3"]}

    PM.anneal(UC.randomStart(16, 10, np.random.default_rng(0)), nbr16, PM.weightMatrix(pos16, (1, 1)),
              chem, fU, 10, 100, 1.0, 0)

    state = json.load(open(OUT)).get("grid", {}) if os.path.exists(OUT) else {}
    state.setdefault("config", {"N_RANDOM": N_RANDOM, "N_WARM": N_WARM, "workers": PM.N_WORKERS,
                                "seed_schedule": "rebuild_grid: 7000 + 100000*cell_index",
                                "warm": "silent lifts drop p1/p2/p3 + previous n=3 cell (as run here)"})
    state.setdefault("n3", [])
    labs = json.load(open(OUT_LAB)) if os.path.exists(OUT_LAB) else {}
    prev3 = None
    for usage, r1, r3 in cells:
        f = fU if usage == "U" else fP
        w = (r1, 1.0, r3)
        W = PM.weightMatrix(pos64, w)
        d2_drop, lifts = [], []
        for p in range(3):
            wd = RG.dropPos(w, p)
            key = (usage, wd) if (usage, wd) in best2 else (usage, wd[::-1])
            lp = RG.silentLift(best2[key]["lab"], c16, c64, p, flip=key[1] != wd)
            d2_drop.append(best2[key]["D"])
            lifts.append(lp)
        lift_vals = [PM.fullD(lp, nbr64, W, chem, f, 10) for lp in lifts]
        for p in range(3):  # a silent lift scores exactly its doublet optimum
            assert abs(lift_vals[p] - d2_drop[p]) <= 1e-9 * max(1.0, d2_drop[p]), (w, p, lift_vals, d2_drop)
        warm, names, seen = [], [], set()
        for p in range(3):
            if lifts[p].tobytes() not in seen:
                seen.add(lifts[p].tobytes()); warm.append(lifts[p]); names.append(f"lift_drop_p{p + 1}")
        # neighbour rule of rebuild_grid: the previous cell in its (r1, r3) order; this run's code if it
        # was just searched here, else rebuild_grid's best labeling of that cell
        ci = cellIndex(usage, r1, r3) - (12 if usage == "U" else 36)
        pw = (PM.R1S[(ci - 1) // 6], 1.0, PM.R3S[(ci - 1) % 6]) if ci > 0 else None
        prev = prev3 if prev3 is not None and prev3[0] == usage and prev3[1] == pw else None
        if prev is None:
            if ci > 0:
                pr1, pr3 = PM.R1S[(ci - 1) // 6], PM.R3S[(ci - 1) % 6]
                m = rg_lab.get(f"n3_{usage}_w{(pr1, 1.0, pr3)}")
                if m is not None:
                    prev = (usage, (pr1, 1.0, pr3), np.array([CLASS_ORDER.index(m[c]) for c in c64], dtype=np.int64),
                            "rebuild_grid")
        if prev is not None and prev[2].tobytes() not in seen:
            warm.append(prev[2]); names.append(f"n3_neighbour{prev[1]}_{prev[3]}")
        job = dict(nbr=nbr64, W=W, dist=chem, f=f, k=10, N=64, n_iter=2000 * 64 * 64,
                   seed0=7000 + 100000 * cellIndex(usage, r1, r3), warm_labs=warm, warm_names=names)
        nw = max(N_WARM, 2 * len(warm))
        res, lab = PM.runSearch(job, N_RANDOM, nw)
        D3s = res["best"]
        min_lift = min(lift_vals)
        cert = D3s <= min_lift + TOL * max(1.0, min_lift)
        D3 = D3s if cert else min_lift
        if not cert:
            lab = lifts[int(np.argmin(lift_vals))]
        d3min2 = min(d2_drop)
        rel = TOL * max(1.0, d3min2)
        verdict = "yes" if abs(D3 - d3min2) <= rel else ("below" if D3 < d3min2 else "ABOVE")
        argmin = [p + 1 for p in range(3) if abs(d2_drop[p] - d3min2) <= rel]
        syn = RG.positionSynonymy(lab, nbr64, pos64, 3)
        d2_inv, d2_term = best2[(usage, (r1, 1.0))]["D"], d2_drop[1]
        old = rg3[(usage, r1, r3)]
        rec = dict(res, usage=usage, r1=r1, r3=r3, w=list(w), restarts=[N_RANDOM, nw],
                   seed0=job["seed0"], warm_names=names,
                   D2_inv=d2_inv, D2_term=d2_term, D3_search=D3s, D3=D3, A_inv=1 - D3 / d2_inv,
                   A_term=1 - D3 / d2_term, D2_drop_p1_p2_p3=d2_drop, lift_values=lift_vals,
                   min_D2=d3min2, argmin_positions=argmin, vs_min_D2=verdict,
                   certified_all_three_lifts=bool(cert), synonymy_by_position=syn,
                   silenced_position=RG.silencedPosition(syn),
                   rebuild_grid_D3=old["best"], rebuild_grid_vs_min_D2=(
                       "yes" if abs(old["best"] - d3min2) <= rel else ("below" if old["best"] < d3min2 else "ABOVE")),
                   D3_change_vs_rebuild_grid=D3 - old["best"])
        state["n3"] = [x for x in state["n3"] if not (x["usage"] == usage and x["r1"] == r1 and x["r3"] == r3)]
        state["n3"].append(rec)
        state["elapsed_seconds_last_call"] = round(time.time() - T0, 1)
        updateResults("grid", state)
        labs[f"n3_{usage}_w{w}"] = {c64[i]: CLASS_ORDER[lab[i]] for i in range(64)}
        json.dump(labs, open(OUT_LAB + ".part", "w"), indent=1)
        os.replace(OUT_LAB + ".part", OUT_LAB)
        prev3 = (usage, w, lab, "this run")
        wh = " ".join(f"{k}:{v['hits']}/{v['restarts']}" for k, v in res.get("warm", {}).items())
        print(f"n=3 {usage} w={w}: search={D3s:.9f} D3*={D3:.9f} lifts={[round(v, 9) for v in lift_vals]} "
              f"{verdict} cert3={cert} sil={rec['silenced_position']} rand {res['random']['hits']}/{N_RANDOM} {wh} "
              f"old={old['best']:.9f} ({res['seconds']}s, t={time.time() - T0:.0f}s)", flush=True)


if __name__ == "__main__":
    spec = sys.argv[1]  # "U3" or "P"
    if spec == "U3":
        cells = [("U", 3.0, r3) for r3 in PM.R3S]
    elif spec == "P":
        # informative cells first (rebuild_grid: 3,1 ABOVE; 3,2 and 1,1 below), then the rest in grid order
        first = [("P", 3.0, 1.0), ("P", 3.0, 2.0), ("P", 1.0, 1.0)]
        cells = first + [("P", r1, r3) for r1 in PM.R1S for r3 in PM.R3S if ("P", r1, r3) not in first]
    else:
        cells = [(spec.split(",")[0], float(spec.split(",")[1]), float(spec.split(",")[2]))]
    main(cells)
