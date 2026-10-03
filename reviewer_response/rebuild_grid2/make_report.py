"""REPORT.md for rebuild_grid2; also merges opt20P_RESULTS.json into RESULTS.json.  run: python make_report.py"""
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
RR = os.path.join(HERE, "..")
TOL = 1e-9
R = json.load(open(os.path.join(HERE, "RESULTS.json"))) if os.path.exists(os.path.join(HERE, "RESULTS.json")) else {}
if os.path.exists(os.path.join(HERE, "opt20P_RESULTS.json")):
    R["opt20P"] = json.load(open(os.path.join(HERE, "opt20P_RESULTS.json")))
    json.dump(R, open(os.path.join(HERE, "RESULTS.json.part"), "w"), indent=1)
    os.replace(os.path.join(HERE, "RESULTS.json.part"), os.path.join(HERE, "RESULTS.json"))
RG = json.load(open(os.path.join(RR, "rebuild_grid", "RESULTS.json")))
PMR = json.load(open(os.path.join(RR, "position_model", "RESULTS.json")))
new = {(x["usage"], x["r1"], x["r3"]): x for x in R.get("grid", {}).get("n3", [])}
rg = {(x["usage"], x["r1"], x["r3"]): x for x in RG["n3"]}


def verdict(d3, m):
    return "yes" if abs(d3 - m) <= TOL * max(1.0, m) else ("below" if d3 < m else "ABOVE")


L = ["# rebuild_grid2: completed rebuilt-matrix grids (fixed certificate) + rebuild_opt20 usage P", "",
     "Read-only reuse: `rebuild_grid/rebuild_grid.py` (imported by `grid_fixed.py`) and "
     "`rebuild_opt20/run.py` (copied to `opt20P.py`, edits marked `rebuild_grid2`). NEW vectors = "
     "`descriptor_rebuild/NEW_VECTORS.csv`.", "",
     "## (A) ten-class position-weight grid, NEW vectors", "",
     "**Provenance.** The `rebuild_grid` run kept going after its REPORT.md was written. It finished all 48 "
     "cells (24 n=2, 24 n=3) with no restart cuts (`rebuild_grid/RESULTS.json`, finished = True). The n=2 "
     "optima used here come from that run. Each doublet label was re-scored and matched its stored D2* to 1e-9. "
     "The n=3 cells listed below as `rerun` were searched again in `grid_fixed.py`, with two fixes:", "",
     "1. **Certificate:** D3* is checked against all three silent lifts (drop p1, drop p2, drop p3). Each lift "
     "is built at its own dropped position, so its value under w equals its doublet optimum exactly (asserted).",
     "2. **Warm starts:** the three silent lifts, plus the previous n=3 cell (rebuild_grid's neighbour rule). "
     "rebuild_grid had only two lift warm starts. Its `lift_Hinv` is the drop-p3 lift. Its `lift_Hterm` "
     "came from `UC.liftLabels`, which is silent at p3, not p2, so it was never the H-term lift. The drop-p1 "
     "lift was missing, and that is the lift that wins whenever r1 = 3 and r3 ≤ 2.", "",
     "Random restarts reuse rebuild_grid's seed schedule (7000 + 100000·cell_index). With 128 random restarts, "
     "the random part reproduces rebuild_grid's random restarts exactly. U r1=1 rows were not rerun; they are "
     "rebuild_grid's values, re-checked here against all three lift values it stored.", "",
     "Columns: weights (r1, 1, r3). D2*(H-inv) = D2*(r1, 1), which drops p3. D2*(H-term) = D2*(r1, r3), which "
     "drops p2. A = 1 − D3*/D2*. 'vs min D2' compares D3* with min_p D2*(w∖p). 'yes' = equal (1e-9 rel), "
     "'below' = the triplet code beats every silent lift, 'ABOVE' = search failure. argmin p is in parentheses. "
     "'sil' = codon position(s) whose arcs are all synonymous in the best n=3 code. 'cert3' = D3*(search) ≤ "
     "min of the three lift values.", ""]
for usage in ("U", "P"):
    L += [f"### usage {usage}", "",
          "| r1 | r3 | D2*(H-inv) | D2*(H-term) | D3* | A_inv | A_term | sil | D3* vs min D2 (argmin p) | cert3 | rand hits | restarts | rebuild_grid D3* (vs min D2) | source |",
          "|:-:|:-:|--:|--:|--:|--:|--:|:-:|:-:|:-:|:-:|:-:|--:|:-:|"]
    for r1 in (1.0, 3.0):
        for r3 in (1.0, 2.0, 5.0, 10.0, 30.0, 100.0):
            o = rg[(usage, r1, r3)]
            om = o["drop"]["min_D2"]
            if (usage, r1, r3) in new:
                x = new[(usage, r1, r3)]
                L.append(f"| {r1:g} | {r3:g} | {x['D2_inv']:.6f} | {x['D2_term']:.6f} | {x['D3']:.9f} | "
                         f"{x['A_inv']:+.4f} | {x['A_term']:+.4f} | {x['silenced_position']} | "
                         f"{x['vs_min_D2']} ({'/'.join(map(str, x['argmin_positions']))}) | "
                         f"{'yes' if x['certified_all_three_lifts'] else 'NO'} | {x['random']['hits']}/"
                         f"{x['random']['restarts']} | {x['restarts'][0]}+{x['restarts'][1]} | "
                         f"{o['best']:.9f} ({verdict(o['best'], om)}) | rerun |")
            elif usage == "U" and r1 == 1.0:
                lv = o["drop"]["lift_value_under_w"]
                c3 = o["best"] <= min(lv) + TOL * max(1.0, min(lv))
                L.append(f"| {r1:g} | {r3:g} | {o['H_inv']['D2']:.6f} | {o['H_term']['D2']:.6f} | {o['best']:.9f} | "
                         f"{o['H_inv']['A']:+.4f} | {o['H_term']['A']:+.4f} | {o['silenced_position']} | "
                         f"{verdict(o['best'], om)} ({'/'.join(map(str, o['drop']['argmin_positions']))}) | "
                         f"{'yes' if c3 else 'NO'} | {o['random']['hits']}/{o['random']['restarts']} | 128+{o['restarts'][1]} | "
                         f"{o['best']:.9f} ({verdict(o['best'], om)}) | rebuild_grid |")
            else:
                # not rerun: rebuild_grid's search, certificate applied against all three lifts (no new search)
                lv = o["drop"]["lift_value_under_w"]
                c3 = o["best"] <= min(lv) + TOL * max(1.0, min(lv))
                d3 = o["best"] if c3 else min(lv)
                sil = o["silenced_position"] if c3 else f"{1 + lv.index(min(lv))} (lift)"
                L.append(f"| {r1:g} | {r3:g} | {o['H_inv']['D2']:.6f} | {o['H_term']['D2']:.6f} | {d3:.9f} | "
                         f"{1 - d3 / o['H_inv']['D2']:+.4f} | {1 - d3 / o['H_term']['D2']:+.4f} | {sil} | "
                         f"{verdict(d3, om)} ({'/'.join(map(str, o['drop']['argmin_positions']))}) | "
                         f"{'yes' if c3 else 'NO → lift'} | {o['random']['hits']}/{o['random']['restarts']} | "
                         f"128+{o['restarts'][1]} | {o['best']:.9f} ({verdict(o['best'], om)}) | rebuild_grid, 3-lift cert only (NOT rerun) |")
    L.append("")

# old-vector comparison for the A columns (position_model), for reference
L += ["### OLD vectors (position_model / finish), for comparison", "",
      "See `rebuild_grid/REPORT.md` (the 'old' side of each 'old → new' column). The OLD D3* values are "
      "unchanged by this rerun.", ""]

O = R.get("opt20P")
L += ["## (B) rebuild_opt20, usage P: n=3 vs n=4 optima at k = 20, OLD vs NEW matrix", ""]
if O:
    N4 = json.load(open(os.path.join(RR, "n4", "RESULTS.json")))

    def fmtw(w):
        return "(" + ", ".join(f"{v:g}" for v in w) + ")"

    def pooled(cells):
        p = {}
        for c in cells:
            key = (c["usage"], tuple(c["w4"]))
            p[key] = min(p.get(key, float("inf")), c["D4"])
        return p
    po, pn = pooled(N4["cells"]), pooled(O["cells"])
    old = {(c["usage"], c["weight_set"], c["hypothesis"]): c for c in N4["cells"]}
    gA = O["gateA_sgc_old"]
    L += [f"Gate A (old vectors, SGC D): {gA['Dpos']!r} vs n4 {gA['target']!r}, |diff| {gA['absdiff']:.1e} -> "
          f"**{'PASS' if gA['pass'] else 'FAIL'}**. Gate B was not rerun; rebuild_opt20 passed it with identical code. "
          f"Search: n=3 two independent batches of {O['config']['N3_BATCHES'][0]} random restarts (SA "
          f"{O['config']['n_iter_n3']}). n=4 random + warm from the triplet lift and all four deletion lifts (SA "
          f"{O['config']['n_iter_n4']}). Restart reductions: {O['reductions'] or 'none'}. The seed counter resumes "
          f"where rebuild_opt20's U block stopped. Finished: {O.get('finished', False)}; elapsed "
          f"{O.get('elapsed_seconds')} s.", "",
          "| weights | hyp | OLD D3* | OLD D4* | OLD G | NEW D3* | NEW D4* | NEW G | NEW cert | NEW rand hits | D4*=min del |",
          "|---|---|---|---|---|---|---|---|---|---|---|"]
    for c in O["cells"]:
        o = old.get(("P", c["weight_set"], c["hypothesis"]))
        nD4 = pn[("P", tuple(c["w4"]))]
        if o:
            oD4 = po[("P", tuple(o["w4"]))]
            os_ = f"{o['D3']:.6f} | {oD4:.6f} | {1 - oD4 / o['D3']:+.4f}"
        else:
            os_ = "— | — | —"
        L.append(f"| {c['weight_set']} {fmtw(c['w3'])} | {c['hypothesis']} | {os_} | {c['D3']:.6f} | {nD4:.6f} | "
                 f"{1 - nD4 / c['D3']:+.4f} | {'yes' if c['certified'] else 'NO'} | {c['random']['hits']}/"
                 f"{c['random']['restarts']} | {'yes' if c['D4_equals_min_deletion'] else 'no'} |")
    L += ["", "### n=3 batch agreement (NEW vectors, usage P)", "",
          "| usage | w sorted | batch A | batch B | agree (1e-6) | hits A / B | D3* |", "|---|---|---|---|---|---|---|"]
    for x in O["n3"]:
        L.append(f"| {x['usage']} | {fmtw(x['w_sorted'])} | {x['batchA_best']:.9f} | {x['batchB_best']:.9f} | "
                 f"{'yes' if x['batches_agree'] else 'NO'} | {x['batchA']['hits']}/{x['batchB']['hits']} | {x['best']:.9f} |")
    na = sum(x["batches_agree"] for x in O["n3"])
    L += ["", f"{na}/{len(O['n3'])} n=3 searches agree across batches.", ""]
else:
    L += ["(not run)", ""]
open(os.path.join(HERE, "REPORT.md"), "w").write("\n".join(L) + "\n")
print("\n".join(L))
