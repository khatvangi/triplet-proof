"""assemble RESULTS.json and REPORT.md from runs/*.json, build_meta.json, table_S2_new.csv."""
import hashlib, json, os
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
L = lambda t: json.load(open(os.path.join(HERE, "runs", f"{t}.json")))  # noqa: E731
TAGS = ["old", "old_unscaled", "new", "new_k8_unscaled"]
R = {t: L(t) for t in TAGS if os.path.exists(os.path.join(HERE, "runs", f"{t}.json"))}
meta = json.load(open(os.path.join(HERE, "build_meta.json")))
f = lambda x, n=4: "—" if x is None else (f"{x:.{n}g}" if isinstance(x, float) else str(x))  # noqa: E731

rows = []  # (statistic, old, new, same, conclusion)


def add(name, get, concl, fmt=lambda v: f(v)):
    try:
        o, n = get(R["old"]), get(R["new"])
    except KeyError:
        rows.append((name, "not run", "not run", "—", "")); return
    same = concl(o) == concl(n)
    rows.append((name, fmt(o), fmt(n), "Y" if same else "N", concl(n) if same else f"old: {concl(o)} / new: {concl(n)}"))


def tail(s):
    return "SGC in lower tail (p<0.05)" if s["p_le"] < 0.05 and s["z"] < 0 else "SGC not in lower tail (p>=0.05)"


def blk(s):
    return "SGC significantly better than block null (p<0.05)" if s["p_strict_better"] < 0.05 else "not significant vs block null (p>=0.05)"


zp = lambda s: f"z={s['z']:.3f}, p≤={s['p_le']:.2g}"  # noqa: E731
zb = lambda s: f"z={s['z']:.3f}, p={s['p_strict_better']:.4g}"  # noqa: E731
for c, lab in (("A", "A doublet proj."), ("B", "B SGC"), ("C", "C RAA10 triplet")):
    for m in ("E", "D"):
        add(f"(1) codon-shuffle {lab} {m}", lambda r, c=c, m=m: r["s1_codon_shuffle"][c][m], tail, zp)
    for m in ("E", "D"):
        add(f"(1) equidistant {lab} {m}", lambda r, c=c, m=m: r["s1_codon_shuffle"][c]["equidistant"][m], tail, zp)
for c, lab in (("A", "A 9! exhaustive"), ("C", "C 10! exhaustive"), ("B", "B 1e7 draws")):
    for m in ("E", "D"):
        add(f"(2) block null {lab} {m}", lambda r, c=c, m=m: r["s2_block_null"][c][m], blk, zb)
add("(3) mean d, 190 pairs", lambda r: r["s3"]["mean_d_190"], lambda v: "scale only (no conclusion)")
for w in ("uniform", "katoh"):
    add(f"(3) SGC mean d / misreading ({w})", lambda r, w=w: r["s3"]["part3"][w]["sgc_mean_d"], lambda v: "scale only")
    add(f"(3) SGC vs codon-shuffle median ({w})", lambda r, w=w: (r["s3"]["part3"][w]["sgc_mean_d"], r["s3"]["part3"][w]["codon_shuffle"]["median"]),
        lambda v: "SGC below codon-shuffle median" if v[0] < v[1] else "SGC above", lambda v: f"{v[0]:.4f} vs {v[1]:.4f}")
    add(f"(3) SGC vs block-null median ({w})", lambda r, w=w: (r["s3"]["part3"][w]["sgc_mean_d"], r["s3"]["part3"][w]["block_null"]["median"]),
        lambda v: "SGC below block-null median" if v[0] < v[1] else "SGC above", lambda v: f"{v[0]:.4f} vs {v[1]:.4f}")
    add(f"(3) best block code sampled / +swap descent ({w})", lambda r, w=w: (r["s3"]["part3"][w]["block_best_sampled"], r["s3"]["part3"][w]["block_best_after_swap_descent"], r["s3"]["part3"][w]["sgc_mean_d"]),
        lambda v: "block codes exist below SGC" if v[1] < v[2] else "SGC is best found", lambda v: f"{v[0]:.4f} / {v[1]:.4f}")
for org in ("ecoli", "yeast"):
    for sub in ("pooled", "no_IL_destination", "no_cys_origin"):
        add(f"(4) {org} {sub}: mean d M1 / M2+", lambda r, o=org, s=sub: (r["s4_landerer"][o][s]["M1"]["mean_d"], r["s4_landerer"][o][s]["M2+"]["mean_d"]),
            lambda v: "M2+ costlier per event than M1" if v[1] > v[0] else "M1 costlier per event", lambda v: f"{v[0]:.4f} / {v[1]:.4f}")
        add(f"(4) {org} {sub}: cost share M1 / M2+", lambda r, o=org, s=sub: (r["s4_landerer"][o][s]["M1"]["cost_share"], r["s4_landerer"][o][s]["M2+"]["cost_share"]),
            lambda v: "M1 carries majority of cost" if v[0] > 0.5 else "M2+ carries majority of cost", lambda v: f"{v[0]:.3f} / {v[1]:.3f}")
    add(f"(4) {org} per-dataset median M2+ cost share", lambda r, o=org: r["s4_landerer"][o]["per_dataset_summary"]["M2plus_cost_share_median"],
        lambda v: "M2+ majority" if v > 0.5 else "M1 majority")
for u in ("U", "P"):
    add(f"(5) Dw_{u} doublet proj. vs SGC 10-class proj.", lambda r, u=u: (r["s5_usage"]["doublet_projection"][f"Dw_{u}"], r["s5_usage"]["sgc_10class_projection"][f"Dw_{u}"]),
        lambda v: "doublet projection lower Dw than SGC 10-class projection" if v[0] < v[1] else "SGC 10-class projection lower Dw", lambda v: f"{v[0]:.4f} vs {v[1]:.4f}")

# old-vector replication targets for the 10! / 9! enumeration
rep = {}
if "s2_block_null" in R.get("old", {}):
    b = R["old"]["s2_block_null"]
    rep = {"C_pD": b["C"]["D"]["p_strict_better"], "C_pE": b["C"]["E"]["p_strict_better"], "A_pD": b["A"]["D"]["p_strict_better"],
           "targets": {"C_pD": 0.34, "C_pE": 0.32, "A_pD": 0.0068}}
    rep["pass"] = round(rep["C_pD"], 2) == 0.34 and round(rep["C_pE"], 2) == 0.32 and abs(rep["A_pD"] - 0.0068) < 5e-5

s2 = pd.read_csv(os.path.join(HERE, "table_S2_new.csv"))
OLD_KEPT = open(os.path.join(os.path.dirname(HERE), "aa_props_export.csv")).readline().strip().split(",")[1:21]
OLD_NAMES = ["hydro_kd", "vol_zam", "polarity_gr", "pI", "sc_pka", "charge_pH7", "donors", "acceptors", "arom", "ring_count",
             "rotb", "refractivity", "bulkiness", "cf_alpha", "cf_beta", "cf_turn", "acc_pref", "flex", "helical_moment",
             "beta_moment", "hphob_fr", "sc_mass"]
ta = ["| # | concept | OLD column (hand-entered, src/io/aa_props_lib.py) | OLD status | NEW source | NEW status |", "|---|---|---|---|---|---|"]
for i, (on, r) in enumerate(zip(OLD_NAMES, s2.itertuples()), 1):
    ost = "kept" if f"{on}_z" in OLD_KEPT else "dropped (|r|>=0.95)"
    ta.append(f"| {i} | {r.concept} | {on} | {ost} | {r.route}: {r.accession_or_rule if isinstance(r.accession_or_rule, str) else ''} | {r.status} |")

gate = json.load(open(os.path.join(HERE, "runs", "gate.json")))
res = {"rules_sha256": open(os.path.join(HERE, "RULES.sha256")).read().split()[0],
       "rules_sha256_now": hashlib.sha256(open(os.path.join(HERE, "RULES.md"), "rb").read()).hexdigest(),
       "gate": gate, "build": meta, "old_vector_block_replication": rep,
       "variants_note": "k97 = 8 for the new matrix, so the fixed-8-PC rescaled variant IS the primary (NEW_VECTORS.csv == NEW_VECTORS_k8_rescaled.csv); unscaled variants run separately",
       "runs": R}
json.dump(res, open(os.path.join(HERE, "RESULTS.json"), "w"), indent=1)

tb = ["| statistic | OLD | NEW | same conclusion | conclusion |", "|---|---|---|---|---|"] + [f"| {a} | {b} | {c} | {d} | {e} |" for a, b, c, d, e in rows]
# variant table: key z-scores for unscaled variants
tv = ["| statistic | old rescaled | old unscaled | new primary (k=8 rescaled) | new k=8 unscaled |", "|---|---|---|---|---|"]
def g(t, fn):
    try: return f"{fn(R[t]):.3f}"
    except KeyError: return "not run"
for c in "ABC":
    for m in ("E", "D"):
        tv.append(f"| (1) z {c} {m} | " + " | ".join(g(t, lambda r: r["s1_codon_shuffle"][c][m]["z"]) for t in TAGS) + " |")
        tv.append(f"| (2) block p {c} {m} | " + " | ".join(g(t, lambda r: r["s2_block_null"][c][m]["p_strict_better"]) for t in TAGS) + " |")
open(os.path.join(HERE, "tables_generated.md"), "w").write(
    "## Table A\n\n" + "\n".join(ta) + "\n\n## Table B\n\n" + "\n".join(tb) + "\n\n## Variants\n\n" + "\n".join(tv) +
    "\n\n## old-vector enumeration replication\n\n" + json.dumps(rep, indent=1) + "\n")
print("\n".join(tb)); print("\n".join(tv)); print(rep)
