"""REPORT.md from RESULTS.json (one table per usage).  run: python reviewer_response/n4/make_report.py"""
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
R = json.load(open(os.path.join(HERE, "RESULTS.json")))


def syn(ps):
    return "/".join(f"{p['raw']:.2f}" for p in ps)


def fmtw(w):
    return "(" + ", ".join(f"{x:g}" for x in w) + ")"


L = ["# Triplet vs quadruplet at the full 20-residue alphabet", "",
     "Dpos (position_model kernel), k = 20 residue vectors (pca1..pca8 z-scored across the 20), "
     "H(n,4) with all 4^n codons sense, free class sizes (>= 1 codon each).", "",
     f"- gate 1 (SGC, 61 sense codons, w = 1, f = n_a/N): Dpos = {R['gate1']['Dpos']!r} vs published "
     f"{R['gate1']['published']!r}, |diff| = {R['gate1']['absdiff']:.1e} -> **{'PASS' if R['gate1']['pass'] else 'FAIL'}**",
     f"- gate 2 (silent lift preserves Dpos, any w4 / any insertion position, 100 random labelings): "
     f"max rel diff {R['gate2']['max_reldiff_any_insert_position']:.1e} -> **{'PASS' if R['gate2']['pass'] else 'FAIL'}**",
     f"- gate 3 (kernel vs brute force at n = 3, 4; position-permutation invariance): "
     f"{R['gate3']['max_reldiff_brute']:.1e} / {R['gate3']['max_reldiff_permutation']:.1e} -> "
     f"**{'PASS' if R['gate3']['pass'] else 'FAIL'}**",
     f"- search: n=3 {R['config']['N3_RANDOM']} random restarts, SA {R['config']['n_iter_n3']}; n=4 "
     f"{R['config']['N4_RANDOM']} random + {R['config']['N4_WARM']} warm (split over the triplet-lift and 4 "
     f"deletion-lift sources), SA {R['config']['n_iter_n4']}; recolor SA + steepest descent (recolor + swap).",
     f"- restart reductions: {R['reductions'] or 'none'}", "",
     "Columns: G = 1 - D4*/D3*. `D4 = min_p D3*(w\\p)` asks whether the best quadruplet code is just a "
     "triplet code (for the weights left after deleting position p) with that position made silent. "
     "`syn` = per-position fraction of synonymous arcs (raw) of the best code. `rand hits` = random "
     "restarts reaching D4*. `cert` = D4* <= every lift value.", ""]

# cells with identical (usage, w4) are the same optimisation problem; searches that disagree there
# show the n=4 search is not converged. the table reports the pooled minimum and flags the spread.
pool = {}
for c in R["cells"]:
    pool.setdefault((c["usage"], tuple(c["w4"])), []).append(c["D4"])
spread = {key: max(v) - min(v) for key, v in pool.items() if len(v) > 1}
L += ["## Search consistency (identical quadruplet weights searched more than once)", ""]
for key, s in spread.items():
    L.append(f"- usage {key[0]}, w4 = {fmtw(key[1])}: {len(pool[key])} searches, D4 values "
             f"{', '.join(f'{x:.6f}' for x in pool[key])}; spread {s:.6f}"
             f"{' **(NOT CONVERGED)**' if s > 1e-9 else ''}")
L += ["", "D4* in the tables below is the minimum over all searches of the same (usage, w4); it is an "
      "upper bound found by search, not a proven optimum. G is exact only if both D3* and D4* are optimal.", ""]

for tag in ("U", "P"):
    cells = [dict(c) for c in R["cells"] if c["usage"] == tag]
    for c in cells:
        c["D4_own"] = c["D4"]
        c["D4"] = min(pool[(c["usage"], tuple(c["w4"]))])
        c["G"] = 1.0 - c["D4"] / c["D3"]
        md = min(d["D3"] for d in c["deletions"])
        c["D4_equals_min_deletion"] = abs(c["D4"] - md) <= 1e-9 * max(1.0, md)
        c["D4_below_min_deletion_by"] = md - c["D4"]
    if not cells:
        continue
    L += [f"## Usage {tag}", "",
          "| weights | hyp | w4 | D3* | D4* | G | min_p D3*(w\\p) [p] | D4 = min del | syn n3 | syn n4 | rand hits | best from | cert |",
          "|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for c in cells:
        md = min(c["deletions"], key=lambda d: d["D3"])
        L.append(f"| {c['weight_set']} {fmtw(c['w3'])} | {c['hypothesis']} | {fmtw(c['w4'])} | {c['D3']:.6f} | "
                 f"{c['D4']:.6f} | {c['G']:+.5f} | {md['D3']:.6f} [{md['p']}] | "
                 f"{'yes' if c['D4_equals_min_deletion'] else 'NO (-%.6f)' % c['D4_below_min_deletion_by']} | "
                 f"{syn(c['pos_synonymy_n3'])} | {syn(c['pos_synonymy_n4'])} | "
                 f"{c['random']['hits']}/{c['random']['restarts']} | {c['best_from']} | "
                 f"{'yes' if c['certified'] else '**SEARCH FAILURE**'} |")
    L.append("")

b = R.get("blind")
if b:
    L += ["## Chemistry-blind (dist = 1), U usage, uniform weights, k = 20", "",
          f"- D3* (n=3, 64 codons) = {b['D3']:.6f}", f"- D4* (n=4, 256 codons) = {b['D4']:.6f} "
          f"(certified vs lifts: {b['certified']})",
          f"- analytic lower bound 3*log4(20) = {b['lower_bound_3log4k']:.6f} (any n)", ""]

L += ["## n=3 searches (sorted weights; codes for other orderings are position permutations)", "",
      "| usage | w sorted | D3* | random hits | syn per position |", "|---|---|---|---|---|"]
for s in R["n3"]:
    L.append(f"| {s['usage']} | {fmtw(s['w_sorted'])} | {s['best']:.9f} | {s['random']['hits']}/"
             f"{s['random']['restarts']} | {syn(s['pos_synonymy'])} |")
L.append("")
n3_hits = [s["random"]["hits"] for s in R["n3"] if s["usage"] != "U-blind"]
L += ["## Caveats", "",
      f"- n=3 chemistry searches reached their best value in {min(n3_hits)}-{max(n3_hits)} of 256 random "
      "restarts, so D3* is also a search upper bound. Where D4* < min_p D3*(w\\p) by a small margin, the "
      "gap may be real or may come from a suboptimal D3*; the two cannot be told apart from this run.",
      "- no n=4 random restart reached D4*; every n=4 best came from an annealed warm start.",
      "- `syn n4` is taken from the cell's own best code, which is not the pooled D4* code when a duplicate "
      "(usage, w4) search did better.", ""]
L.append(f"elapsed {R.get('elapsed_seconds')} s; finished = {R.get('finished', False)}")
open(os.path.join(HERE, "REPORT.md"), "w").write("\n".join(L) + "\n")
print("\n".join(L))
