"""(3) search part: k=10 U n=2 optimum (usage_cert R2 procedure, 256 random restarts) with the
r=0.94 vectors, then the n=3 lift check: silent lift of the n=2 optimum onto 64 codons must score
exactly D2, and an n=3 search (128 random + 16 warm from the lift, UC ladder) must not go below it.
run from repo root:  python reviewer_response/rebuild_close/part3_search.py"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
os.environ["NUMBA_CACHE_DIR"] = os.path.join(HERE, "numba_cache")
sys.dont_write_bytecode = True
import numba.core.config  # noqa: E402,F401
import json  # noqa: E402
import time  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

RR = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(RR, "position_model"))
import position_model as PM  # noqa: E402
os.environ["NUMBA_CACHE_DIR"] = os.path.join(HERE, "numba_cache")
numba.core.config.reload_config()
UC, CLASS_ORDER = PM.UC, PM.CLASS_ORDER
UC.N_WORKERS = 62
OUT = os.path.join(HERE, "part3_search.json")


def main():
    t0 = time.time()
    df = pd.read_csv(os.path.join(HERE, "VECTORS_r094.csv"))
    x = df[[c for c in df.columns if c.startswith("pc")]].to_numpy(float)
    x = (x - x.mean(0)) / x.std(0, ddof=0)
    vec = {a: x[i] for i, a in enumerate(df["aa"].str.upper())}
    chem = PM.distMatrix(CLASS_ORDER, {c: vec[a] for c, a in PM.REPS.items()})
    fU = UC.uniformUsage(10, list(range(10)))
    c64, nbr64 = PM.buildSenseGraph(3)
    c16, nbr16 = PM.buildSenseGraph(2)
    UC.anneal(UC.randomStart(16, 10, np.random.default_rng(0)), nbr16, chem, fU, 10, 100, 1.0, 0)
    job2 = dict(nbr=nbr16, dist=chem, f=fU, k=10, N=16, n_iter=10000 * 256, classes=CLASS_ORDER,
                seed0=15_000_000, warm_lab=None)
    r2, lab2 = UC.runSearch(job2, 256, 0)
    print("n2", r2["best_Dw"], r2["random"]["hits"], round(time.time() - t0), flush=True)
    lift = UC.liftLabels(lab2, c16, c64)
    lift_val = UC.fullDw(lift, nbr64, chem, fU, 10)
    out = {"n2": r2, "n2_lab": {c16[i]: CLASS_ORDER[lab2[i]] for i in range(16)}, "lift_Dw": lift_val,
           "lift_exact_absdiff": abs(lift_val - r2["best_Dw"])}
    json.dump(out, open(OUT, "w"), indent=1, default=lambda o: o.item())
    job3 = dict(nbr=nbr64, dist=chem, f=fU, k=10, N=64, n_iter=2000 * 64 * 64, classes=CLASS_ORDER,
                seed0=15_100_000, warm_lab=lift)
    UC_NW = UC.N_WARM
    UC.N_WARM = 16  # the warm temperature ladder is spread over N_WARM restarts
    r3, lab3 = UC.runSearch(job3, 128, 16)
    UC.N_WARM = UC_NW
    out.update(n3=r3, D3=r3["best_Dw"], certified=bool(r3["best_Dw"] <= lift_val + 1e-9),
               below_lift=bool(r3["best_Dw"] < lift_val - 1e-9), seconds=round(time.time() - t0, 1))
    json.dump(out, open(OUT, "w"), indent=1, default=lambda o: o.item())
    print("n3", r3["best_Dw"], "lift", lift_val, out["certified"], round(time.time() - t0), flush=True)


if __name__ == "__main__":
    main()
