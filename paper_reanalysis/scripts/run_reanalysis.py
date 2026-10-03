#!/usr/bin/env python3
"""the frozen grid: 2 nulls x 5 property representations x 2 objectives.

both nulls run through ONE scoring path, so the contrast cannot be an artefact
of two implementations.  gated on reproducing the published SGC values.

    python scripts/run_reanalysis.py --n 1000000
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
PROJ = HERE.parent
sys.path.insert(0, str(PROJ.parent))

from paper_reanalysis.src.reanalysis import (NULLS, Space, positiveControl,
                                             propertySets, score, zstats)

OUT = PROJ / "results"
BATCH = 2000
SEED = 20260923


def sweep(sp, props, null_name, n, rng, batch=BATCH):
    """accumulate E and D for every property set over one null."""
    keys = list(props)
    acc = {k: {"E": np.empty(n), "D": np.empty(n)} for k in keys}
    made = 0
    while made < n:
        b = min(batch, n - made)
        L = NULLS[null_name](sp, b, rng)
        for k in keys:
            E, D = score(props[k], L, sp)
            acc[k]["E"][made:made + b] = E
            acc[k]["D"][made:made + b] = D
        made += b
    return acc


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=1_000_000)
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)

    sp = Space()
    props = propertySets(sp)
    pc = positiveControl(sp, props)
    print("POSITIVE CONTROL:", json.dumps(pc), flush=True)
    if not pc["passed"]:
        print("positive control FAILED -- stopping, per the analysis plan.")
        sys.exit(1)

    obs = {k: dict(zip(("E", "D"),
                       [float(v[0]) for v in score(props[k], sp.labels[None, :], sp)]))
           for k in props}

    rows = []
    t0 = time.time()
    for null_name in ("codon", "block"):
        rng = np.random.default_rng(SEED)
        acc = sweep(sp, props, null_name, args.n, rng)
        for k in props:
            for objective in ("E", "D"):
                st = zstats(obs[k][objective], acc[k][objective])
                st.update({"null": null_name, "properties": k,
                           "objective": objective})
                rows.append(st)
        print(f"  {null_name} null done  [{time.time()-t0:.0f}s]", flush=True)

    t = pd.DataFrame(rows)[["null", "properties", "objective", "sgc",
                            "null_mean", "null_sd", "z", "n_better", "n_tie",
                            "n", "p_laplace", "null_min", "zero_variance"]]
    t.to_csv(OUT / "grid.csv", index=False)

    pd.set_option("display.width", 200)
    for null_name in ("codon", "block"):
        print(f"\n=== null = {null_name} ===")
        print(t[t.null == null_name][["properties", "objective", "sgc",
                                      "null_mean", "null_sd", "z", "n_better",
                                      "p_laplace", "zero_variance"]]
              .to_string(index=False, float_format=lambda v: f"{v:,.4f}"))

    # cross-validation against the preregistered 10^7 block-null run
    ref = {"E": -1.9153, "D": -2.5369}
    chk = {}
    for objective in ("E", "D"):
        got = float(t[(t.null == "block") & (t.properties == "pc8_whitened")
                      & (t.objective == objective)].z.iloc[0])
        chk[objective] = {"this_run": got, "block_null_D6_10e7": ref[objective],
                          "abs_diff": abs(got - ref[objective])}
    print("\n=== cross-validation vs the preregistered 10^7 block null ===")
    print(json.dumps(chk, indent=2))
    (OUT / "crossvalidation.json").write_text(json.dumps(chk, indent=2) + "\n")
    print(f"\nelapsed {time.time()-t0:.0f}s -> {OUT}")


if __name__ == "__main__":
    main()
