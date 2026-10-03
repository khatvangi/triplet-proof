"""(2) exact closed-form block-null mean of D and E for the SGC with NEW vectors, vs the sampled
mean in descriptor_rebuild (B: 1e7 draws) and the exhaustive means (A 9!, C 10!) as a positive control.
under a uniform permutation P of labels among blocks, for an edge between distinct blocks a!=b,
E[d(P(a),P(b))] = mean of d over ordered pairs of distinct labels (same for d^2); within-block edges
contribute 0. so E[D] = m_between * dbar / N and E[E] = m_between * d2bar / 2, m_between = directed
between-block edges. (4) all raw-descriptor pairs with |r| >= 0.93.
run from repo root:  python reviewer_response/rebuild_close/part2_4.py"""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
os.environ["NUMBA_CACHE_DIR"] = os.path.join(HERE, "numba_cache")
sys.dont_write_bytecode = True
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

RR = os.path.dirname(HERE)
DR = os.path.join(RR, "descriptor_rebuild")
sys.path.insert(0, DR)
import stats as ST  # noqa: E402  (read-only import)


def closedForm(cond):
    lab, E, dist = cond["lab"], cond["E"], cond["dist"]
    k = len(cond["classes"])
    between = int((lab[E[:, 0]] != lab[E[:, 1]]).sum())
    off = ~np.eye(k, dtype=bool)
    dbar, d2bar = dist[off].mean(), (dist[off] ** 2).mean()
    return {"directed_edges": int(len(E)), "between_block_directed_edges": between, "N": cond["N"], "k": k,
            "mean_d_distinct_ordered": float(dbar), "mean_d2_distinct_ordered": float(d2bar),
            "D_mean": float(between * dbar / cond["N"]), "E_mean": float(between * d2bar / 2.0)}


def main():
    vec, cols = ST.loadVectors(os.path.join(DR, "NEW_VECTORS.csv"), "pc", True)
    tabs = ST.tables()
    conds = {"A": ST.condition(tabs["A"], {c: vec[a] for c, a in ST.REPS.items()}),
             "B": ST.condition(tabs["B"], vec),
             "C": ST.condition(tabs["C"], {c: vec[a] for c, a in ST.REPS.items()})}
    run = json.load(open(os.path.join(DR, "runs", "new.json")))
    out = {"vectors": "descriptor_rebuild/NEW_VECTORS.csv (pc1..pc8, rescaled as stats.py)", "part2": {}}
    for c, cond in conds.items():
        cf = closedForm(cond)
        # sgc reproduces the run's value -> same edge set / distances
        assert abs(ST.DE(cond["lab"], cond, cond["dist"])[0] - run["sgc_D"][c]) < 1e-9
        bn = run["s2_block_null"][c]
        row = {"closed_form": cf, "exhaustive": bn["exhaustive"], "n_draws": bn["n_draws"]}
        for m in ("D", "E"):
            s = bn[m]
            se = s["null_sd"] / np.sqrt(s["n"])
            diff = s["null_mean"] - cf[f"{m}_mean"]
            row[m] = {"sampled_mean": s["null_mean"], "closed_form_mean": cf[f"{m}_mean"], "diff": float(diff),
                      "null_sd": s["null_sd"], "se": float(se), "diff_in_se": float(diff / se),
                      "rel_diff": float(diff / cf[f"{m}_mean"])}
        out["part2"][c] = row
        print(c, {m: (row[m]["sampled_mean"], row[m]["closed_form_mean"], row[m]["diff_in_se"]) for m in "DE"})
    out["part2"]["stats_py_draws"] = ("re-drawn: blockNull() builds fresh permutations with "
                                      "np.random.default_rng(SEED=20260924) (B) or itertools.permutations (A, C); "
                                      "no archived draws are read")
    # (4)
    raw = pd.read_csv(os.path.join(DR, "descriptors_raw.csv")).set_index("aa")
    R = np.corrcoef(raw.to_numpy(float).T)
    names = list(raw.columns)
    pairs = [{"a": names[i], "b": names[j], "r": float(R[i, j]), "abs_r": float(abs(R[i, j]))}
             for i in range(len(names)) for j in range(i + 1, len(names)) if abs(R[i, j]) >= 0.93]
    out["part4"] = {"n_descriptors": len(names), "pairs_abs_r_ge_0.93": sorted(pairs, key=lambda p: -p["abs_r"]),
                    "max_abs_r": float(np.abs(R - np.eye(len(names))).max())}
    for p in out["part4"]["pairs_abs_r_ge_0.93"]:
        print(p)
    json.dump(out, open(os.path.join(HERE, "part2_4.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
