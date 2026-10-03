"""(1) stop-codon variant with the NEW vectors: finish/ part (1) repeated exactly
(usage_cert kernels, 2000*N^2 schedule, 512 random + 64 warm, UC warm temperature ladder),
k=10 RAA10 reps, uniform weights, free sizes, n=3 on the 61 sense codons (stops dropped),
usage U and P. the only change vs finish/ is the class distance matrix (rebuild_opt10's
loadNewVectors) and the warm source (silent lift of rebuild_opt10's NEW n=2 optima).
two independent seed batches. writes part1.json incrementally.
run from repo root:  python reviewer_response/rebuild_close/part1_stop.py"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
os.environ["NUMBA_CACHE_DIR"] = os.path.join(HERE, "numba_cache")
sys.dont_write_bytecode = True
import numba.core.config  # noqa: E402,F401

import json  # noqa: E402
import time  # noqa: E402

import numpy as np  # noqa: E402

RR = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(RR, "position_model"))
import position_model as PM  # noqa: E402  (read-only; brings usage_cert as PM.UC)
import pandas as pd  # noqa: E402

# position_model may repoint the numba cache on import; pin it back to rebuild_close/
os.environ["NUMBA_CACHE_DIR"] = os.path.join(HERE, "numba_cache")
numba.core.config.reload_config()
UC, CLASS_ORDER = PM.UC, PM.CLASS_ORDER


def loadNewVectors():
    """copy of rebuild_opt10.loadNewVectors: pc1..pc8 z-scored ddof=0, at the RAA10 reps."""
    df = pd.read_csv(os.path.join(RR, "descriptor_rebuild", "NEW_VECTORS.csv"))
    x = df[[f"pc{i}" for i in range(1, 9)]].to_numpy(float)
    x = (x - x.mean(axis=0)) / x.std(axis=0, ddof=0)
    aa_vec = {aa: x[i] for i, aa in enumerate(df["aa"].str.upper())}
    return {cls: aa_vec[aa] for cls, aa in PM.REPS.items()}
UC.N_WORKERS = 62
D2_OPT = {"U": 21.344573500, "P": 19.187429770}
SEEDS = {"batch1": 11_000_000, "batch2": 13_000_000}  # disjoint from usage_cert/finish/rebuild_opt10
OUT = os.path.join(HERE, "part1.json")


def main():
    t0 = time.time()
    chem = PM.distMatrix(CLASS_ORDER, loadNewVectors())
    ro = json.load(open(os.path.join(RR, "rebuild_opt10", "RESULTS.json")))
    assert np.allclose(chem, np.array(ro["dist_new"]), atol=1e-12, rtol=0)
    ro_lab = json.load(open(os.path.join(RR, "rebuild_opt10", "best_labelings.json")))
    comp, prot = UC.loadProteome()
    assert not prot["fallback_used"]
    fU = UC.uniformUsage(10, list(range(10)))
    fP = UC.classUsage(comp, set(CLASS_ORDER))
    c61, nbr61 = PM.buildSenseGraph(3, PM.STOPS)
    c16, nbr16 = PM.buildSenseGraph(2)
    assert len(c61) == 61 and nbr61.max() < 61
    UC.anneal(UC.randomStart(16, 10, np.random.default_rng(0)), nbr16, chem, fU, 10, 100, 1.0, 0)
    out = {"config": {"N_RANDOM": UC.N_RANDOM, "N_WARM": UC.N_WARM, "n_iter": 2000 * 61 * 61,
                      "seeds": SEEDS, "workers": UC.N_WORKERS}, "runs": []}
    for batch, s0 in SEEDS.items():
        for idx, (usage, f, key) in enumerate((("U", fU, "a_R2_n2"), ("P", fP, "a_R3_n2"))):
            lab2 = np.array([CLASS_ORDER.index(ro_lab[key][c]) for c in c16], dtype=np.int64)
            d2 = UC.fullDw(lab2, nbr16, chem, f, 10)
            assert abs(d2 - D2_OPT[usage]) < 1e-8, (usage, d2)
            warm = UC.liftLabels(lab2, c16, c61)
            job = dict(nbr=nbr61, dist=chem, f=f, k=10, N=61, n_iter=2000 * 61 * 61, classes=CLASS_ORDER,
                       seed0=s0 + 100000 * idx, warm_lab=warm)
            res, lab = UC.runSearch(job, UC.N_RANDOM, UC.N_WARM)
            res.update(batch=batch, usage=usage, D2_opt=D2_OPT[usage], warm_lift_Dw=UC.fullDw(warm, nbr61, chem, f, 10),
                       gap=D2_OPT[usage] - res["best_Dw"],
                       gap_percent=100 * (D2_OPT[usage] - res["best_Dw"]) / D2_OPT[usage],
                       best_lab={c61[i]: CLASS_ORDER[lab[i]] for i in range(61)})
            out["runs"].append(res)
            out["elapsed_s"] = round(time.time() - t0, 1)
            json.dump(out, open(OUT + ".part", "w"), indent=1, default=lambda o: o.item())
            os.replace(OUT + ".part", OUT)
            print(f"{batch} {usage}: best={res['best_Dw']:.9f} rand {res['random']['hits']}/{res['random']['restarts']} "
                  f"warm {res['warm']['hits']}/{res['warm']['restarts']} gap%={res['gap_percent']:.3f} "
                  f"lift={res['warm_lift_Dw']:.6f} t={time.time() - t0:.0f}s", flush=True)
    out["finished"] = True
    json.dump(out, open(OUT, "w"), indent=1, default=lambda o: o.item())


if __name__ == "__main__":
    main()
