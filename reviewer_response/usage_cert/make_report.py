"""build REPORT.md from RESULTS.json (no numbers are typed by hand)."""
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
R = json.load(open(os.path.join(HERE, "RESULTS.json")))


def sizes(d):
    return " ".join(str(v) for v in sorted(d.values(), reverse=True))


def status(r):
    if r["n"] == 2:
        return ""
    if r["stops_excluded"]:
        return "report only; " + ("<= restricted lift" if r["below_restricted_lift"] else "ABOVE restricted lift")
    if r.get("SEARCH_FAILURE"):
        return f"**SEARCH FAILURE**, n=3 upper bound = lift = {r['lift_Dw']:.9f}"
    s = "certified (<= lift)"
    if not r["random_only_certified"]:
        s += "; random-only best ABOVE lift"
    return s


g1, g2, p = R["gate1"], R["gate2"], R["proteome"]
L = ["# Usage-weighted doublet vs triplet optima, with lift certification", ""]
L += ["Script: `usage_cert.py` (from repo root: `python reviewer_response/usage_cert/usage_cert.py`); "
      "raw numbers in `RESULTS.json`, best labelings in `best_labelings.json`, log in `run.log`.", ""]
L += ["Metric: Dw = sum_a f_a (1/n_a) sum_{c in block(a)} sum_{c' in N(c)} dist(a, label(c')). "
      "Primary scenarios use all 4^n codons as sense, free class sizes (each >= 1).", ""]
L += ["## Gates", "",
      f"- **Gate 1** (f_a = n_a/N_sense reproduces the old D): triplet projection, stops excluded "
      f"Dw = {g1['triplet_Dw']!r} vs {g1['triplet_expected']!r} (|diff| = {g1['triplet_absdiff']:.2e}); "
      f"doublet majority-vote projection Dw = {g1['doublet_Dw']!r} vs {g1['doublet_expected']!r} "
      f"(|diff| = {g1['doublet_absdiff']:.2e}). **{'PASS' if g1['pass'] else 'FAIL'}** (tol 1e-9).",
      f"- **Gate 2** (silent lift L3(xyz) = L2(xy) on all 64 codons): doublet projection, uniform f over its 9 classes: "
      f"Dw(L2) = {g2['projection']['doublet_uniform9_n2']!r}, Dw(L3) = {g2['projection']['lift_n3']!r}; "
      f"plus {g2['n_random_labelings']} random doublet labelings with {g2['random_f']}. "
      f"Max |Dw(L3) - Dw(L2)| = {g2['max_absdiff']:.2e} (floating-point summation order only). "
      f"**{'PASS' if g2['pass'] else 'FAIL'}** (tol {g2['tolerance']:.0e}).", ""]
L += ["## Proteome usage (P)", ""]
if p["fallback_used"]:
    L += [f"Source page unreachable or unparsable ({p.get('error')}); **used the hardcoded fallback composition**.", ""]
else:
    L += [f"UniProtKB/Swiss-Prot release **{p['release']}** of **{p['date']}**, section 6.1 "
          f"'Composition in percent for the complete database', from {p['source']} "
          f"(raw page saved as `relstat.html`). The 20 values sum to {p['sum_percent']}% (rounding); "
          "class frequencies are sums over RAA10 members, renormalised to 1.", ""]
L += ["| class | " + " | ".join(R["usage_P_classes"]) + " |",
      "|---|" + "---|" * len(R["usage_P_classes"]),
      "| f_P | " + " | ".join(f"{v:.4f}" for v in R["usage_P_classes"].values()) + " |", ""]

rc = R["reference_codes"]
L += ["## Reference codes", "",
      "| code | Dw (U) | Dw (P) |", "|---|---|---|",
      f"| doublet majority-vote projection (9 classes; X8 absent, f renormalised over the 9) | "
      f"{rc['doublet_projection']['Dw_U']:.6f} | {rc['doublet_projection']['Dw_P']:.6f} |",
      f"| SGC 10-class triplet projection, UAA/UAG/UGA excluded | "
      f"{rc['sgc_triplet_projection_stops_excluded']['Dw_U']:.6f} | "
      f"{rc['sgc_triplet_projection_stops_excluded']['Dw_P']:.6f} |", ""]

L += ["## Optima", "",
      "Hits = restarts ending within 1e-9 of the row's best (random and lift-warm-started counted separately); "
      "p5/p50 over random restarts. G = Dw*(n=2) - Dw*(n=3).", "",
      "| row | n | k | usage | best Dw | random hits | warm hits | p5 | p50 | class sizes (best) | lift bound | G | G % | check |",
      "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
for r in R["rows"]:
    warm = f"{r['warm']['hits']}/{r['warm']['restarts']}" if "warm" in r else "-"
    lift = f"{r['lift_Dw']:.6f}" if "lift_Dw" in r else "-"
    G = f"{r['gap_G']:.6f}" if "gap_G" in r else "-"
    Gp = f"{r['gap_G_percent']:.3f}" if "gap_G_percent" in r else "-"
    st = status(r)
    if "reached_analytic" in r:
        st = (st + "; " if st else "") + f"3log4k={r['analytic_3log4k']:.0f} " + \
             ("reached" if r["reached_analytic"] else "NOT reached")
    if r.get("random_restarts_cut"):
        st += "; random restarts cut for budget"
    L.append(f"| {r['row']} | {r['n']}{' (61 sense)' if r['stops_excluded'] else ''} | {r['k']} | {r['usage']} | "
             f"{r['best_Dw']:.9f} | {r['random']['hits']}/{r['random']['restarts']} | {warm} | "
             f"{r['random']['p5']:.6f} | {r['random']['p50']:.6f} | {sizes(r['class_sizes_best'])} | "
             f"{lift} | {G} | {Gp} | {st} |")
LAB = json.load(open(os.path.join(HERE, "best_labelings.json")))


def silentPositions(lab):
    """codon positions along which the labeling is constant (a lift from the other two)."""
    out = []
    for p_ in range(3):
        if all(lab[c] == lab[c[:p_] + b + c[p_ + 1:]] for c in lab for b in "ACGU"
               if c[:p_] + b + c[p_ + 1:] in lab):
            out.append(p_ + 1)
    return out


L += ["", "## Structure of the n=3 optima", "",
      "Silent position = a codon position along which the best labeling never changes class, i.e. the "
      "labeling is a lift of a doublet code on the other two positions (H(3,4) is symmetric under "
      "position permutation, so the lift need not be on position 3).", "",
      "| row | k | silent positions of best n=3 labeling |", "|---|---|---|"]
for r in R["rows"]:
    if r["n"] == 3:
        sp = silentPositions(LAB[f"{r['row']}_k{r['k']}_n3"])
        L.append(f"| {r['row']} | {r['k']} | {', '.join(map(str, sp)) if sp else 'none (not a lift)'} |")
L += ["", "## Notes", "",
      "- Search: single-codon recolor SA (every class kept nonempty), geometric cooling by 1e4, "
      "then steepest descent to a strict local optimum (improvement threshold 1e-12) under all legal single "
      "recolors and all pairwise label swaps. Schedule lengths (iterations): R1 1000*N^2; chemistry n=2 "
      "10000*N^2; chemistry n=3 2000*N^2 (chosen from a timing probe on the loaded machine, `probe.py`).",
      "- Warm starts (n=3 chemistry rows, 64 each): restart 0 is pure descent from the lift of the best n=2 "
      "labeling; restarts 1..63 anneal from the lift with initial temperature scaled 1e-3..1e-1. Because "
      "restart 0 can only descend, warm starts alone guarantee best <= lift; the `random-only` flag "
      "therefore reports whether the *random* restarts also reached or passed the lift bound.",
      "- Caveat: with recolor-only SA, a singleton class cannot move directly (recoloring it would empty "
      "the class); it can only grow first. Swaps enter only at the final descent. This may contribute to "
      "hit rates below 100% on chemistry rows (not tested separately).",
      "- R4 rows (stops excluded, 61 sense codons): the lift is not exact there; the `lift bound` column is the "
      "Dw of the lifted best n=2 labeling restricted to the 61 sense codons. G for R4 compares against the "
      "all-sense n=2 optimum and is not a certified gap.",
      "- A row whose best is hit by only a few random restarts (see hits) is a best-found value, not a "
      "certified optimum.",
      f"- Random restarts cut for budget: {'yes' if any(r.get('random_restarts_cut') for r in R['rows']) else 'no'}.",
      "", f"Total wall time: {R['elapsed_seconds']:.0f} s."]
open(os.path.join(HERE, "REPORT.md"), "w").write("\n".join(L) + "\n")
print("\n".join(L))
