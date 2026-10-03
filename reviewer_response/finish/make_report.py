"""REPORT.md from finish/RESULTS.json (two short tables)."""
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
R = json.load(open(os.path.join(HERE, "RESULTS.json")))

L = ["# finish: R4 stop-codon variant and missing position_model P cells", ""]
L += ["## (1) R4: n=3 on 61 sense codons (stops excluded), uniform weights, k=10", "",
      "| usage | best Dw (this run) | hits random / warm | usage_cert R4 | best over both | n=2 optimum | gap | gap % | syn p1 / p2 / p3 (raw) |",
      "|:-:|--:|:-:|--:|--:|--:|--:|--:|:-:|"]
for r in R["part1_R4"]:
    syn = " / ".join(f"{s['raw']:.3f}" for s in r["synonymy_by_position"])
    L.append(f"| {r['usage']} | {r['best_Dw']:.9f} | {r['random']['hits']}/{r['random']['restarts']} / "
             f"{r['warm']['hits']}/{r['warm']['restarts']} | {r['usage_cert_R4_best']:.9f} | {r['best_over_both_runs']:.9f} | "
             f"{r['Dw_n2_opt']:.9f} | {r['gap']:.6f} | {r['gap_percent']:.2f} | {syn} |")
L += ["", "gap = n=2 optimum − best n=3 (61-sense) Dw; positive = the triplet code does better. synonymy = "
      "fraction of sense→sense arcs at that position whose endpoints share a class, in the best code from this run. "
      "Usage-weighted synonymy is in RESULTS.json.", "",
      "caveats: this run uses usage_cert's kernels and schedule on disjoint seeds. For U it did NOT reproduce usage_cert's "
      "best (18.866 vs 18.857; hit rates 4/512 and 1/512), so the U optimum is not converged and 18.857 is only an upper "
      "bound. For P both runs landed on the same value. No warm start from the n=2 lift reached the best in either usage. "
      "The synonymy profile belongs to one best code. Other codes with equal Dw may spread synonymy across the positions "
      "differently, so do not read the p1-vs-p3 split as unique.", ""]

g = R.get("part2_gate") or {}
L += ["## (2) position_model, usage P, n=3 cells that were missing", "",
      f"missing before this run: {R.get('part2_missing')}. reproduction gate (re-running cell (3,1,2) with its original "
      f"seeds): **{'PASS' if g.get('pass') else 'FAIL'}** (best {g.get('best_new')} vs {g.get('best_orig')}).", "",
      "| r1 | r3 | D3 best | hits random | warm hits | A (H-inv) | A (H-term) | certified | syn3 raw | SGC excess % |",
      "|:-:|:-:|--:|:-:|:-:|--:|--:|:-:|--:|--:|"]
for r in R["part2_cells"]:
    wh = ", ".join(f"{v['hits']}/{v['restarts']}" for v in r["warm"].values())
    L.append(f"| {r['r1']:g} | {r['r3']:g} | {r['best']:.9f} | {r['random']['hits']}/{r['random']['restarts']} | {wh} | "
             f"{r['H_inv']['A']:+.5f} | {r['H_term']['A']:+.5f} | {r['H_inv']['certified']} | "
             f"{r['pos3_synonymy']['raw']:.3f} | {r['sgc']['excess_percent']:.1f} |")
L += ["", "warm hits are in the order lift H-inv, lift H-term, previous n=3 cell."]
open(os.path.join(HERE, "REPORT.md"), "w").write("\n".join(L) + "\n")
print("\n".join(L))
