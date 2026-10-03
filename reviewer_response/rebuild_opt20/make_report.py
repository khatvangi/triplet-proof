"""REPORT.md: OLD (n4/RESULTS.json) vs NEW (rebuild_opt20/RESULTS.json).  run: python make_report.py"""
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
R = json.load(open(os.path.join(HERE, "RESULTS.json")))
O = json.load(open(os.path.join(HERE, "..", "n4", "RESULTS.json")))


def fmtw(w):
    return "(" + ", ".join(f"{x:g}" for x in w) + ")"


def pooled(cells):
    """D4* = min over every search of the same (usage, w4), as in n4/REPORT.md."""
    p = {}
    for c in cells:
        key = (c["usage"], tuple(c["w4"]))
        p[key] = min(p.get(key, float("inf")), c["D4"])
    return p


po, pn = pooled(O["cells"]), pooled(R["cells"])
old = {(c["usage"], c["weight_set"], c["hypothesis"]): c for c in O["cells"]}
gA, gB = R["gateA_sgc_old"], R["gateB_n3_old"]
L = ["# n=3 vs n=4 optima at k = 20: OLD vs rebuilt descriptor matrix", "",
     "NEW = `descriptor_rebuild/NEW_VECTORS.csv` (8 PCs, each z-scored). OLD = `n4/RESULTS.json`. "
     "Kernels, graph, usage, search and lift/certification logic are imported unchanged from `n4/n4.py`. "
     "G = 1 - D4*/D3*. D4* is pooled over searches of the same (usage, w4), as in n4.", "",
     "## Gates (old vectors)", "",
     f"- SGC D = {gA['Dpos']!r} vs n4 gate 1 {gA['target']!r}, |diff| = {gA['absdiff']:.1e} -> "
     f"**{'PASS' if gA['pass'] else 'FAIL'}**",
     f"- n=3 U UNIF, n4's seed block, 256 restarts: best {gB['best']!r} vs n4 {gB['n4_best']!r} "
     f"(|diff| {gB['absdiff']:.1e}; n4 spread best..p50 = {gB['n4_best']:.6f}..{gB['n4_p50']:.6f}), "
     f"hits {gB['hits']}/256 -> **{'PASS' if gB['pass'] else 'FAIL'}**",
     f"- info: SGC D with NEW vectors = {R.get('sgc_D_new_info', float('nan')):.6f}", "",
     f"Search: n=3 two independent batches of {R['config']['N3_BATCHES'][0]} random restarts (SA "
     f"{R['config']['n_iter_n3']}); n=4 {R['config']['N4_RANDOM']} random + {R['config']['N4_WARM']} warm "
     f"(SA {R['config']['n_iter_n4']}). Restart reductions: {R['reductions'] or 'none'}.", ""]

for tag in ("U", "P"):
    rows = [c for c in R["cells"] if c["usage"] == tag]
    L += [f"## Usage {tag}", "",
          "| weights | hyp | OLD D3* | OLD D4* | OLD G | NEW D3* | NEW D4* | NEW G | NEW cert | NEW rand hits |",
          "|---|---|---|---|---|---|---|---|---|---|"]
    if not rows:
        L.append("| (not run) | | | | | | | | | |")
    for c in rows:
        o = old[(tag, c["weight_set"], c["hypothesis"])]
        oD4, nD4 = po[(tag, tuple(o["w4"]))], pn[(tag, tuple(c["w4"]))]
        L.append(f"| {c['weight_set']} {fmtw(c['w3'])} | {c['hypothesis']} | {o['D3']:.6f} | {oD4:.6f} | "
                 f"{1 - oD4 / o['D3']:+.4f} | {c['D3']:.6f} | {nD4:.6f} | {1 - nD4 / c['D3']:+.4f} | "
                 f"{'yes' if c['certified'] else 'NO'} | {c['random']['hits']}/{c['random']['restarts']} |")
    L.append("")

L += ["## n=3 batch agreement (NEW vectors)", "",
      "| usage | w sorted | batch A | batch B | agree (1e-6) | hits A / B | D3* |", "|---|---|---|---|---|---|---|"]
for x in R["n3"]:
    L.append(f"| {x['usage']} | {fmtw(x['w_sorted'])} | {x['batchA_best']:.9f} | {x['batchB_best']:.9f} | "
             f"{'yes' if x['batches_agree'] else 'NO'} | {x['batchA']['hits']}/{x['batchB']['hits']} | {x['best']:.9f} |")
na = sum(x["batches_agree"] for x in R["n3"])
L += ["", f"{na}/{len(R['n3'])} n=3 searches agree across batches. Hits are counted against the "
      "pooled best, so a batch that missed it shows 0.", ""]
L += ["## Conclusion (usage U only; P was not run)", "",
      "- **Q-inv: unchanged.** G = 0 in all 6 NEW cells. In every one the best n=4 code is the silent "
      "position-4 lift of the best triplet code. OLD was also 0 in all non-UNIF cells.",
      "- **Q-mon at in vitro rates: unchanged in kind, slightly larger.** K-str 2.9% -> 3.0%, K-int "
      "4.1% -> 5.6%. That is still 'a few percent'. Graded weights: w3=2 5.7% -> 6.3%, w3=5 17.9% -> 20.0%, "
      "w3=10 29.9% -> 33.7%. In every Q-mon cell D4* equals min_p D3*(w\\p). So G here is a ratio of two "
      "n=3 optima: the quadruplet gain is just deleting the heavy third position.",
      "- **UNIF: OLD G = +1.3%, NEW G = 0.** OLD found an n=4 code 0.30 below every lift. NEW found none, "
      "with only 96+16 restarts. This cell is search-limited, not a matrix effect that can be claimed.",
      "- Absolute D values are not comparable across matrices (the distance matrix changed; SGC D is "
      "23.92 OLD vs 23.84 NEW). Only G is comparable.", "",
      "## Caveats", "",
      "- **The n=3 search is not converged.** 0/17 two-batch pairs agree to 1e-6. Every batch best is hit "
      "by exactly 1 of 256 restarts, and the batch-to-batch gap is 0.01-0.23 (up to 0.7% of D3). Every D3* "
      "and D4* is a search upper bound, and G carries an uncertainty of that order. This is enough to "
      "move the UNIF cell and the smaller Q-mon gains by about +/-1 percentage point. It is not enough to "
      "turn Q-inv's 0 into a gain, or the w3=5/10 gains into nothing.",
      "- n=4 restarts were cut to the floor (96 random + 16 warm) in all 12 cells because of the budget. "
      "No random n=4 restart reached D4*; every best came from a warm start.",
      "- Usage P was not run: the U block took 1074 s on a machine shared with other jobs "
      "(load average ~87 on 64 cores), and the run was stopped before P to stay near the 28-minute budget.",
      ""]
L += [f"elapsed {R['elapsed_seconds']} s; finished = {R.get('finished', False)}"]
open(os.path.join(HERE, "REPORT.md"), "w").write("\n".join(L) + "\n")
