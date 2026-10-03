"""(3) correlation-threshold sensitivity, statistics part (no searches).
rebuilds the matrix from descriptor_rebuild/descriptors_raw.csv with build_descriptors.py's own
pipeline (z ddof=0, ordered pruning vs earlier kept, SVD, k97, each pc z-scored ddof=0), threshold
0.94 instead of 0.95; gate: threshold 0.95 must reproduce NEW_VECTORS.csv distances to 1e-9.
then runs stats.py's own functions (read-only import): codon shuffle (1e6), block null A 9! and
C 10! exhaustive, B 1e6 draws (also for the primary, same seed = its first 1e6 draws), part3, Landerer.
run from repo root:  python reviewer_response/rebuild_close/part3_stats.py"""
import json
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
os.environ["NUMBA_CACHE_DIR"] = os.path.join(HERE, "numba_cache")
sys.dont_write_bytecode = True
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

RR = os.path.dirname(HERE)
DR = os.path.join(RR, "descriptor_rebuild")
sys.path.insert(0, DR)
import stats as ST  # noqa: E402

OUT = os.path.join(HERE, "part3.json")


def buildMatrix(thr):
    raw = pd.read_csv(os.path.join(DR, "descriptors_raw.csv")).set_index("aa")
    Z = (raw - raw.mean()) / raw.std(ddof=0)
    C = Z.corr()
    kept, dropped = [], []
    for c in Z.columns:
        hit = [k for k in kept if abs(C.loc[c, k]) >= thr]
        if hit:
            dropped.append({"dropped": c, "by": hit[0], "abs_r": float(abs(C.loc[c, hit[0]]))})
        else:
            kept.append(c)
    Xz = Z[kept].to_numpy()
    U, S, Vt = np.linalg.svd(Xz, full_matrices=False)
    var = S ** 2 / (Xz.shape[0] - 1)
    cum = np.cumsum(var / var.sum())
    k97 = int(np.searchsorted(cum, 0.97) + 1)
    pcs = (Xz @ Vt.T)[:, :k97]
    resc = (pcs - pcs.mean(0)) / pcs.std(0, ddof=0)
    df = pd.DataFrame(resc, index=raw.index, columns=[f"pc{i + 1}" for i in range(k97)])
    return df, dict(threshold=thr, kept=kept, n_kept=len(kept), dropped=dropped, k97=k97,
                    cum_var=cum.tolist(), cum_var_at_k97=float(cum[k97 - 1]), cum_var_at_8=float(cum[7]))


def distOf(df):
    X = df.to_numpy()
    return np.sqrt(((X[:, None] - X[None]) ** 2).sum(-1))


def main():
    t0 = time.time()
    df95, m95 = buildMatrix(0.95)
    prim = pd.read_csv(os.path.join(DR, "NEW_VECTORS.csv")).set_index("aa")
    gate = float(np.abs(distOf(df95) - distOf(prim)).max())
    assert gate < 1e-9, gate
    df94, m94 = buildMatrix(0.94)
    path94 = os.path.join(HERE, "VECTORS_r094.csv")
    df94.to_csv(path94)
    res = {"gate_0.95_reproduces_primary_maxabs_dist_diff": gate, "build_0.95": m95, "build_0.94": m94}
    print("gate", gate, "dropped", m94["dropped"], "k97", m94["k97"], flush=True)

    def dump():
        res["elapsed_s"] = round(time.time() - t0, 1)
        json.dump(res, open(OUT + ".part", "w"), indent=1)
        os.replace(OUT + ".part", OUT)
    dump()
    tabs = ST.tables()
    for tag, path in (("primary", os.path.join(DR, "NEW_VECTORS.csv")), ("r094", path94)):
        vec, cols = ST.loadVectors(path, "pc", True)
        conds = {"A": ST.condition(tabs["A"], {c: vec[a] for c, a in ST.REPS.items()}),
                 "B": ST.condition(tabs["B"], vec),
                 "C": ST.condition(tabs["C"], {c: vec[a] for c, a in ST.REPS.items()})}
        r = {"columns": cols}
        if tag == "r094":
            r["s1_codon_shuffle"] = {k: ST.codonShuffle(c) for k, c in conds.items()}
            print(tag, "s1", round(time.time() - t0), flush=True)
            r["s2_block_null"] = {"A": ST.blockNull(conds["A"], True), "C": ST.blockNull(conds["C"], True),
                                  "B": ST.blockNull(conds["B"], False, n=1_000_000)}
        else:
            r["s2_block_null"] = {"B_1e6": ST.blockNull(conds["B"], False, n=1_000_000)}
        print(tag, "s2", round(time.time() - t0), flush=True)
        res[tag] = r
        dump()
        if tag == "r094":
            r["s3"] = ST.part3and1(vec, do_part1=False)
            print(tag, "s3", round(time.time() - t0), flush=True)
            dump()
            ST.ds.EB = os.path.join(HERE, "decode_r094")
            os.makedirs(ST.ds.EB, exist_ok=True)
            r["s4_landerer"] = ST.part3and1(vec, do_part1=True)["part1"]
            print(tag, "s4", round(time.time() - t0), flush=True)
            dump()
    res["finished"] = True
    dump()


if __name__ == "__main__":
    main()
