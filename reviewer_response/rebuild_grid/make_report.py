"""REPORT.md: OLD (position_model + finish for the missing P cells) vs NEW (rebuild_grid) grid."""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.dont_write_bytecode = True
RR = os.path.join(HERE, "..")
sys.path.insert(0, os.path.join(RR, "arch_cert"))
import numpy as np  # noqa: E402
from arch_cert import CLASS_ORDER, buildSenseGraph  # noqa: E402

TOL = 1e-9
pm = json.load(open(os.path.join(RR, "position_model", "RESULTS.json")))
fin = json.load(open(os.path.join(RR, "finish", "RESULTS.json")))
labs_old = {**json.load(open(os.path.join(RR, "position_model", "best_labelings.json"))),
            **json.load(open(os.path.join(RR, "finish", "best_labelings.json")))}
R = json.load(open(os.path.join(HERE, "RESULTS.json")))
labs_new = json.load(open(os.path.join(HERE, "best_labelings.json")))

c64, nbr64 = buildSenseGraph(3)


def silenced(m):
    """positions (1-based) whose every arc is synonymous in labeling m."""
    lab = np.array([CLASS_ORDER.index(m[c]) for c in c64])
    out = []
    for p in range(3):
        ok = all(lab[v] == lab[i] for i in range(64) for v in nbr64[i] if v >= 0 and c64[i][p] != c64[v][p])
        if ok:
            out.append(p + 1)
    return "none" if not out else "/".join(map(str, out))


def cellView(n2, n3, labs, usage, r1, r3):
    c = next((x for x in n3 if x["usage"] == usage and x["r1"] == r1 and x["r3"] == r3), None)
    if c is None:
        return None
    d2 = {tuple(x["w"]): x["best"] for x in n2 if x["usage"] == usage}
    drops = [d2[(1.0, r3)], d2[(r1, r3)], d2[(r1, 1.0)]]  # drop p1, p2, p3
    mn = min(drops)
    D3 = c["best"]
    rel = "yes" if abs(D3 - mn) <= TOL * max(1.0, mn) else ("below" if D3 < mn else "ABOVE")
    return {"rel": rel, "gap": D3 - mn, "mn": mn, "D2inv": c["H_inv"]["D2"], "D2term": c["H_term"]["D2"], "D3": D3, "Ainv": c["H_inv"]["A"],
            "Aterm": c["H_term"]["A"], "sil": silenced(labs[f"n3_{usage}_w{(r1, 1.0, r3)}"]),
            "eq": abs(D3 - mn) <= TOL * max(1.0, mn), "argmin": "/".join(str(p + 1) for p in range(3)
                                                                         if abs(drops[p] - mn) <= TOL * max(1.0, mn)),
            "cert": c["H_inv"]["certified"], "cut": c.get("restarts_cut", False),
            "hits": f"{c['random']['hits']}/{c['random']['restarts']}"}


old_n3 = pm["n3"] + fin["part2_cells"]
L = ["# rebuild_grid: ten-class position-weight grid, OLD vs NEW descriptors", "",
     f"NEW vectors: `{R['vectors']}` (z-scored: {R['zscore']}), reps {R['representatives']}. Machinery imported read-only from "
     "`position_model/` (kernels, runSearch, 128 random + max(32, 2·sources) warm restarts, warm sources, seed "
     "schedule 7000+100000·cell, certification). Only the 10×10 distance matrix changed.", "",
     "OLD = `position_model/RESULTS.json`, plus `finish/RESULTS.json` for the four P cells at r1=3, r3∈{5,10,30,100}.", ""]
g = R["gate"]
L += [f"**Gate** (OLD vectors, {g['cell']}, original seeds and warm start): best {g['best_new']!r} vs original "
      f"{g['best_orig']!r}, |Δ| = {g['absdiff']:.1e}; random hits {g['random_new']['hits']}/128 vs "
      f"{g['random_orig']['hits']}/128 → **{'PASS' if g['pass'] else 'FAIL'}**.", ""]
if R.get("reductions"):
    L += [f"Restart reductions (to the 96+16 floor): {len(R['reductions'])} cells: " + "; ".join(R["reductions"]), ""]
else:
    L += ["Restart reductions: none (every cell ran the full 128 + 32).", ""]
if not R.get("finished"):
    L += ["**RUN NOT FINISHED**: the cells missing from the NEW columns were not run.", ""]

L += ["## Findings", "",
      "- Every qualitative feature of the grid survives the descriptor rebuild. At r1=1, H-inv gives A=0 in 11 of "
      "the 12 cells (all except P (1,1,1)): the best triplet code is the p3-silenced lift. H-term A rises with r3 (U: 0.305→0.976; P: 0.268→0.973). "
      "At r1=3, r3≥5, A_inv=0 with p3 silenced. The only cells with A_inv>0 are r1=3, r3∈{1,2}, where p1 is "
      "silenced (or, in P (3,2), no position is), and P (1,1,1).",
      "- The only cells where the triplet beats every silent lift ('below') are P (1,1,1) (A=+0.0157 new vs +0.0037 old) and "
      "P (3,1,2) (26.082 vs min D2 26.221). The same two cells were 'below' with the old vectors.",
      "- NEW search failures: in 3 cells (U (3,1,1), U (3,1,2), P (3,1,1)) the search ended ABOVE the drop-p1 silent lift. "
      "position_model's H-inv certificate passes all three anyway, because it never checks the p1 lift. With the old vectors "
      "these cells reached the lift value. Their true D3* equals min D2 or lies below it; the corrected A values are listed below. "
      "The r1=3, r3=1/2 cells need a rerun with a drop-p1 lift warm start before anyone quotes them.",
      "- All numbers shift by only a few percent (D2*/D3* up 6.6–8.5% in U; in P they move between −0.6% and +1.9%). The SGC-excess column changes by at most "
      "about 60 percentage points, at r3=100 where the excess is already about 700–1600%.", ""]
L += ["Weights: triplet (r1, 1, r3). H-inv doublet (r1, 1), which drops p3. H-term doublet (r1, r3), which drops p2. "
      "A = 1 − D3*/D2*. min D2 = min over p of D2*(w minus p), where dropping p1 gives (1, r3). "
      "'=min' means D3* equals min D2 to 1e-9 relative. 'sil' lists the codon positions whose arcs are all synonymous in the best n=3 code. "
      "'cert' is position_model's certificate: D3* ≤ value of the H-inv silent lift. "
      "In the min-D2 column, 'yes' means D3* equals min D2. 'below' means the triplet code beats every silent lift. "
      "'ABOVE' means the search ended higher than a known silent lift: a **search failure**. That lift is a better "
      "n=3 code, so the true D3* is ≤ min D2. position_model's certificate misses this because it checks only the H-inv (drop-p3) lift.", ""]
for usage in ("U", "P"):
    L += [f"## usage {usage}", "",
          "| r1 | r3 | D2*(H-inv) old → new | D2*(H-term) old → new | D3* old → new | A_inv old → new | "
          "A_term old → new | sil old → new | D3*=min D2 (argmin p) old → new | cert old / new | rand hits new |",
          "|:-:|:-:|--:|--:|--:|--:|--:|:-:|:-:|:-:|:-:|"]
    for r1 in (1.0, 3.0):
        for r3 in (1.0, 2.0, 5.0, 10.0, 30.0, 100.0):
            o = cellView(pm["n2"], old_n3, labs_old, usage, r1, r3)
            n = cellView(R["n2"], R["n3"], labs_new, usage, r1, r3)
            f = lambda v, k, fmt: "—" if v is None else format(v[k], fmt)  # noqa: E731
            yn = lambda v: "—" if v is None else f"{v['rel']} ({v['argmin']})"  # noqa: E731
            L.append(f"| {r1:g} | {r3:g} | {f(o,'D2inv','.6f')} → {f(n,'D2inv','.6f')} | "
                     f"{f(o,'D2term','.6f')} → {f(n,'D2term','.6f')} | {f(o,'D3','.6f')} → {f(n,'D3','.6f')} | "
                     f"{f(o,'Ainv','+.4f')} → {f(n,'Ainv','+.4f')} | {f(o,'Aterm','+.4f')} → {f(n,'Aterm','+.4f')} | "
                     f"{f(o,'sil','s')} → {f(n,'sil','s')} | {yn(o)} → {yn(n)} | "
                     f"{f(o,'cert','')} / {f(n,'cert','')} | {f(n,'hits','s')}{' (cut)' if n and n['cut'] else ''} |")
    L += [""]

fails = [(u, x["r1"], x["r3"], x["best"], x["drop"]["min_D2"], x["drop"]["argmin_positions"],
          x["drop"]["lift_value_under_w"]) for u in ("U", "P") for x in R["n3"]
         if x["usage"] == u and not x["drop"]["certified_vs_min_D2"]]
L += ["## NEW cells where the search ended above min D2 (search failures)", ""]
if not fails:
    L += ["none", ""]
for u, r1, r3, d3, mn, am, lv in fails:
    c = next(x for x in R["n3"] if x["usage"] == u and x["r1"] == r1 and x["r3"] == r3)
    L.append(f"- {u} (r1={r1:g}, r3={r3:g}): search D3* = {d3:.9f} > min D2 = {mn:.9f} (drop p{"/".join(map(str, am))}). "
             f"Silent-lift values under w for dropping p1/p2/p3: {', '.join(f'{v:.9f}' for v in lv)}. "
             f"Using the best known n=3 value (the lift): A_inv = {1 - mn / c['H_inv']['D2']:+.5f}, "
             f"A_term = {1 - mn / c['H_term']['D2']:+.5f}. Treat these as bounds; this cell is not certified.")
L += [""]
L += ["## SGC 10-class projection (61 sense codons) under each triplet weight", "",
      "| usage | r1 | r3 | SGC Dpos old | SGC Dpos new | excess over D3* old % | excess over D3* new % | new, over best known (min of D3*, min D2) % |",
      "|:-:|:-:|:-:|--:|--:|--:|--:|--:|"]
for s in R["sgc"]:
    old = next(x for x in pm["sgc"] if x["usage"] == s["usage"] and x["r1"] == s["r1"] and x["r3"] == s["r3"])
    assert abs(old["Dpos"] - s["Dpos_old_vectors"]) <= 1e-9 * old["Dpos"]
    o = next((x for x in old_n3 if x["usage"] == s["usage"] and x["r1"] == s["r1"] and x["r3"] == s["r3"]), None)
    n = next((x for x in R["n3"] if x["usage"] == s["usage"] and x["r1"] == s["r1"] and x["r3"] == s["r3"]), None)
    eo = "—" if o is None else f"{100 * (old['Dpos'] - o['best']) / o['best']:.1f}"
    en = "—" if n is None else f"{n['sgc']['excess_percent']:.1f}"
    bk = "—" if n is None else min(n["best"], n["drop"]["min_D2"])
    eb = "—" if n is None else f"{100 * (s['Dpos'] - bk) / bk:.1f}"
    L.append(f"| {s['usage']} | {s['r1']:g} | {s['r3']:g} | {old['Dpos']:.6f} | {s['Dpos']:.6f} | {eo} | {en} | {eb} |")
L += ["", "The SGC Dpos with OLD vectors recomputed here matches position_model's stored values to 1e-9 (asserted).", ""]
open(os.path.join(HERE, "REPORT.md"), "w").write("\n".join(L) + "\n")
print("\n".join(L))
