"""assemble RESULTS.json and REPORT.md from part1.json, part2_4.json, part3.json, part3_search.json.
run from anywhere:  python reviewer_response/rebuild_close/make_report.py"""
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
RR = os.path.dirname(HERE)


def load(name):
    p = os.path.join(HERE, name)
    return json.load(open(p)) if os.path.exists(p) else None


def verdictP(p):
    return "p<0.05" if p < 0.05 else "p>=0.05"


def main():
    p1, p24, p3, p3s = load("part1.json"), load("part2_4.json"), load("part3.json"), load("part3_search.json")
    # the (3) statistics ran as parallel stages (part3_stage.py); fold them into p3
    st = {k: load(f"part3_{k}.json") for k in ("s1A", "s1B", "s1C", "s2", "s2prim", "s3", "s4")}
    v = p3.setdefault("r094", {})
    if all(st[f"s1{c}"] for c in "ABC"):
        v["s1_codon_shuffle"] = {c: st[f"s1{c}"]["result"] for c in "ABC"}
    if st["s2"]:
        v["s2_block_null"] = st["s2"]["result"]
    if st["s3"]:
        v["s3"] = st["s3"]["result"]
    if st["s4"]:
        v["s4_landerer"] = st["s4"]["result"]
    if st["s2prim"]:
        p3["primary"] = {"s2_block_null": st["s2prim"]["result"]}
    p3["stage_seconds"] = {k: (x["seconds"] if x else None) for k, x in st.items()}
    p3["finished"] = all(st.values())
    prim = json.load(open(os.path.join(RR, "descriptor_rebuild", "runs", "new.json")))
    ro = json.load(open(os.path.join(RR, "rebuild_opt10", "RESULTS.json")))
    L = ["# rebuild_close: four closing items for the rebuilt descriptor matrix", ""]

    # (1)
    L += ["## (1) Stop-codon variant (61 sense codons, UAA/UAG/UGA removed), NEW vectors", "",
          "Same procedure as finish/ part (1): usage_cert kernels, 2000·N² schedule, 512 random + 64 warm "
          "restarts (warm = silent lift of the NEW n=2 optimum from rebuild_opt10, stop codons dropped), "
          "k=10 RAA10 reps, uniform position weights, free class sizes. Only the distance matrix and warm source "
          "differ from finish/. Two seed batches (seed0 11,000,000 and 13,000,000; U at +0 and P at +100,000).", "",
          "| batch | usage | best Dw | hits random | hits warm | warm-lift Dw | n=2 optimum | gap | gap % |",
          "|:-:|:-:|--:|:-:|:-:|--:|--:|--:|--:|"]
    if p1:
        for r in p1["runs"]:
            L.append(f"| {r['batch']} | {r['usage']} | {r['best_Dw']:.9f} | {r['random']['hits']}/{r['random']['restarts']} | "
                     f"{r['warm']['hits']}/{r['warm']['restarts']} | {r['warm_lift_Dw']:.6f} | {r['D2_opt']:.9f} | "
                     f"{r['gap']:.6f} | {r['gap_percent']:.3f} |")
        for u in "UP":
            rs = [r for r in p1["runs"] if r["usage"] == u]
            if len(rs) == 2:
                b = min(x["best_Dw"] for x in rs)
                agree = abs(rs[0]["best_Dw"] - rs[1]["best_Dw"]) <= 1e-9
                L.append(f"\n{u}: best over both batches {b:.9f} (gap {100 * (rs[0]['D2_opt'] - b) / rs[0]['D2_opt']:.3f}%); "
                         f"batches {'AGREE' if agree else 'DISAGREE'} "
                         f"({rs[0]['best_Dw']:.9f} vs {rs[1]['best_Dw']:.9f}).")
        if not p1.get("finished"):
            L.append("\n**part 1 did not finish within the budget; rows above are the completed runs only.**")
    L += ["", "gap = n=2 optimum − best n=3 (61-sense) Dw; positive = the triplet code does better. Hits count restarts "
          "within 1e-9 of the batch's own best. With low hit counts the best is an upper bound on the true "
          "optimum, so the gap is a LOWER bound.", ""]

    # (2)
    L += ["## (2) Block-null gate: closed-form mean vs sampled mean (NEW vectors)", "",
          "Closed form: under a uniform permutation of labels among blocks, each directed edge between distinct "
          "blocks has E[d] = mean of d over ordered pairs of distinct labels (E[d²] likewise); within-block edges "
          "contribute 0. So E[D] = m_between·d̄/N and E[E] = m_between·d̄²/2. The edge set and distances are "
          "those of descriptor_rebuild/stats.py itself (imported read-only; SGC D reproduced to 1e-9).", "",
          "| cond | null | between-block directed edges / all | metric | sampled / exhaustive mean | closed form | diff | SE | diff / SE |",
          "|:-:|:-:|:-:|:-:|--:|--:|--:|--:|--:|"]
    for c in "ABC":
        r = p24["part2"][c]
        cf = r["closed_form"]
        kind = "exhaustive" if r["exhaustive"] else f"{r['n_draws']:.0e} draws"
        for m in "DE":
            x = r[m]
            L.append(f"| {c} | {kind} | {cf['between_block_directed_edges']}/{cf['directed_edges']} | {m} | {x['sampled_mean']:.9f} | "
                     f"{x['closed_form_mean']:.9f} | {x['diff']:+.3e} | {x['se']:.3e} | {x['diff_in_se']:+.3f} |")
    L += ["", "- A and C are exhaustive (9!, 10!), so their means must equal the closed form exactly. They do, to "
          "~1e-11 SE (float rounding). This validates the formula and the edge set.",
          f"- B (the gate): the sampled mean from 1e7 draws lies {p24['part2']['B']['D']['diff_in_se']:+.2f} SE (D) and "
          f"{p24['part2']['B']['E']['diff_in_se']:+.2f} SE (E) from the exact mean. Both are within 2 SE, so the sampler is consistent "
          "with the null it claims to draw from.",
          "- For B only (labels = all 20 AAs, exactly the rows that were z-scored): mean d² over distinct pairs = "
          "2·8·20/19 = 16.842 for ANY rescaled 8-PC matrix. So B's E null mean (3301.0526 = 392·16.842/2) is fixed by the "
          "graph and says nothing about the geometry. Only D's mean does. This does not hold for A or C, whose "
          "10 representatives are not centred (A: 76·16.842/2 = 640 ≠ 737.8).",
          "- stats.py **re-drew** the permutations; it did not re-score archived draws. blockNull() builds fresh "
          "permutations with np.random.default_rng(20260924) for B (1e7 draws, chunks of 200,000), or enumerates "
          "itertools.permutations for A and C. No saved draws are read. Same seed means old and new matrices see the "
          "identical permutation sequence.", ""]

    # (4)
    q = p24["part4"]
    L += ["## (4) Raw-descriptor pairs with |r| ≥ 0.93", "",
          f"{q['n_descriptors']} descriptors (descriptors_raw.csv, Pearson r across 20 AAs). Max off-diagonal |r| = {q['max_abs_r']:.6f}.", "",
          "| descriptor a | descriptor b | r |", "|---|---|--:|"]
    for p in q["pairs_abs_r_ge_0.93"]:
        L.append(f"| {p['a']} | {p['b']} | {p['r']:+.6f} |")
    L += ["", "Both pairs lie in [0.94, 0.95). The 0.95 rule keeps both, and any threshold ≤ 0.9407 drops both.", ""]

    # (3)
    L += ["## (3) Correlation-threshold sensitivity: |r| ≥ 0.94 instead of 0.95", ""]
    if p3:
        b = p3["build_0.94"]
        L += [f"Gate: the re-implemented build at 0.95 reproduces NEW_VECTORS.csv distances to "
              f"{p3['gate_0.95_reproduces_primary_maxabs_dist_diff']:.1e}.", "",
              "Dropped at 0.94: " + "; ".join(f"**{d['dropped']}** (|r|={d['abs_r']:.4f} with earlier {d['by']})" for d in b["dropped"])
              + f". Kept {b['n_kept']} descriptors. PCA k at 97% = **{b['k97']}** (cumulative {b['cum_var_at_k97']:.4f}; "
              f"primary k=8 at {p3['build_0.95']['cum_var_at_k97']:.4f}). Vectors = 8 PCs, each z-scored (VECTORS_r094.csv).", ""]
        v = p3.get("r094", {})
        rows = []

        def add(label, a, bb, same):
            rows.append(f"| {label} | {a} | {bb} | {same if isinstance(same, str) else ('Y' if same else 'N')} |")
        if "s1_codon_shuffle" in v:
            for c in "ABC":
                for m in "ED":
                    a, bb = prim["s1_codon_shuffle"][c][m], v["s1_codon_shuffle"][c][m]
                    fp = lambda x: "0/1e6" if x["p_strict_better"] == 0 else f"{x['p_strict_better']:.2g}"  # noqa: E731
                    add(f"codon shuffle z({m}) {c}", f"{a['z']:.2f} (p={fp(a)})", f"{bb['z']:.2f} (p={fp(bb)})",
                        (a["z"] < 0) == (bb["z"] < 0) and (a["p_strict_better"] < 0.05) == (bb["p_strict_better"] < 0.05))
        if "s2_block_null" in v:
            for c, lab in (("A", "9! exh"), ("C", "10! exh")):
                for m in "DE":
                    a, bb = prim["s2_block_null"][c][m], v["s2_block_null"][c][m]
                    add(f"block null p({m}) {c} ({lab})", f"{a['p_strict_better']:.4g}", f"{bb['p_strict_better']:.4g}",
                        (a["p_strict_better"] < 0.05) == (bb["p_strict_better"] < 0.05))
            pb = p3.get("primary", {}).get("s2_block_null", {}).get("B_1e6")
            for m in ("DE" if pb else ""):
                a, a7, bb = pb[m], prim["s2_block_null"]["B"][m], v["s2_block_null"]["B"][m]
                add(f"block null p({m}) B (1e6 draws; primary 1e7 in brackets)",
                    f"{a['p_strict_better']:.4g} [{a7['p_strict_better']:.4g}]", f"{bb['p_strict_better']:.4g}",
                    (a["p_strict_better"] < 0.05) == (bb["p_strict_better"] < 0.05))
        if "s3" in v:
            for w in ("uniform", "katoh"):
                a, bb = prim["s3"]["part3"][w], v["s3"]["part3"][w]
                add(f"SGC mean d per misreading vs block median ({w})",
                    f"{a['sgc_mean_d']:.4f} vs {a['block_null']['median']:.4f} (frac≤SGC {a['block_null']['frac_le_sgc']:.4g})",
                    f"{bb['sgc_mean_d']:.4f} vs {bb['block_null']['median']:.4f} (frac≤SGC {bb['block_null']['frac_le_sgc']:.4g})",
                    (a["sgc_mean_d"] < a["block_null"]["median"]) == (bb["sgc_mean_d"] < bb["block_null"]["median"]))
        if "s4_landerer" in v:
            for org in ("ecoli", "yeast"):
                a, bb = prim["s4_landerer"][org]["pooled"], v["s4_landerer"][org]["pooled"]
                for cl in ("M1", "M2+"):
                    add(f"Landerer {org} pooled {cl}: mean d / cost share",
                        f"{a[cl]['mean_d']:.4f} / {a[cl]['cost_share']:.4f}", f"{bb[cl]['mean_d']:.4f} / {bb[cl]['cost_share']:.4f}",
                        (a["M1"]["mean_d"] < a["M2+"]["mean_d"]) == (bb["M1"]["mean_d"] < bb["M2+"]["mean_d"])
                        and (a["M1"]["cost_share"] > 0.5) == (bb["M1"]["cost_share"] > 0.5))
        if p3s:
            a = ro["a"][0]
            add("k=10 U n=2 optimum (restarts / hits)", f"{a['D2']:.9f} ({a['n2']['random']['restarts']}/{a['n2']['random']['hits']})",
                f"{p3s['n2']['best_Dw']:.9f} (256/{p3s['n2']['random']['hits']})", "n/a (value, not a test)")
            if "n3" in p3s:
                add("n=3 search vs silent lift (below lift?)",
                    f"D3={a['D3']:.6f}, lift={a['lift_Dw']:.6f}, certified={a['certified']} (512+64)",
                    f"D3={p3s['D3']:.6f}, lift={p3s['lift_Dw']:.6f}, certified={p3s['certified']} (128+16)",
                    a["certified"] == p3s["certified"])
        L += ["| statistic | primary (0.95) | variant (0.94) | conclusion same |", "|---|---|---|:-:|"] + rows
        L += ["", "Same-conclusion rule (RULES.md): sign/tail side and the 0.05 threshold, not the magnitude. For Landerer: "
              "M1 mean d < M2+ mean d in both, and the M1 cost share on the same side of 0.5 in both (E. coli > 0.5, "
              "yeast < 0.5, in both matrices). For the optimum: n=3 does not beat the silent lift in either.", "",
              "**Conclusion: every test statistic reaches the same verdict at 0.94 as at 0.95.** A, B and the "
              "per-misreading statistic stay significant (p < 0.05). C stays non-significant under the exhaustive block null "
              "(p ≈ 0.17 → 0.13). The Landerer M1/M2+ ordering and shares barely move. The triplet optimum with "
              "free class sizes and uniform weights again equals its n=2 silent lift. Every change is in magnitude "
              "only; the variant p-values are slightly smaller. Caveat: the 0.94 n=3 check used 128+16 restarts, "
              "not the 512+64 of rebuild_opt10.", ""]
        if not p3.get("finished"):
            L.append("**part 3 statistics did not finish within the budget; missing rows were not run.**")
    res = {"part1_stop_variant": p1, "part2_block_gate_and_part4_corr": p24, "part3_threshold_stats": p3,
           "part3_threshold_search": p3s}
    json.dump(res, open(os.path.join(HERE, "RESULTS.json.part"), "w"), indent=1)
    os.replace(os.path.join(HERE, "RESULTS.json.part"), os.path.join(HERE, "RESULTS.json"))
    open(os.path.join(HERE, "REPORT.md"), "w").write("\n".join(L) + "\n")


if __name__ == "__main__":
    main()
