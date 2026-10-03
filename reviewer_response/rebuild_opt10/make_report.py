"""REPORT.md for rebuild_opt10: OLD (usage_cert/, position_empirical/) vs NEW (RESULTS.json) for every number."""
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
RR = os.path.join(HERE, "..")
R = json.load(open(os.path.join(HERE, "RESULTS.json")))
UCR = json.load(open(os.path.join(RR, "usage_cert", "RESULTS.json")))
PE = json.load(open(os.path.join(RR, "position_empirical", "RESULTS.json")))


def ucRow(name, n):
    return next(r for r in UCR["rows"] if r["row"] == name and r["n"] == n)


def f6(x):
    return "—" if x is None else f"{x:.6f}"


def sign(x, eps=1e-9):
    return "+" if x > eps else ("0" if abs(x) <= eps else "−")


L = ["# Ten-class optima under the rebuilt descriptor matrix", "",
     f"Vectors: `{R['vectors']}` ({R['zscore']}), representatives "
     f"{', '.join(f'{k}={v}' for k, v in R['representatives'].items())}. Proteome usage P: Swiss-Prot "
     f"{R['proteome'].get('release')}. Kernels, annealing + steepest descent, lift certification and "
     "the position grid are imported unchanged from `usage_cert/usage_cert.py` and "
     "`position_model/position_model.py` (the module `position_empirical/` imports); only the class "
     "distance matrix changes. OLD = `usage_cert/RESULTS.json`, `position_empirical/RESULTS.json`.", "",
     "## Gate", "",
     f"Old vectors (aa_props_export.csv, pca1..8 z-scored), usage_cert seed list for R2 U n=2 "
     f"(seed0 = 401000). Seed {R['gate']['tries'][0]['seed']} gave {R['gate']['tries'][0]['Dw']:.9f}; "
     f"seed {R['gate']['seed_used']} gave {R['gate']['Dw']:.9f}, |diff| from 20.017516292 = "
     f"{R['gate']['absdiff']:.1e} → **{'PASS' if R['gate']['pass'] else 'FAIL'}**. "
     f"SGC projection (old vectors, U) reproduces usage_cert to {R['gate']['sgc_old_U_absdiff']:.1e}.", ""]

# ── (a) ──
L += ["## (a) usage_cert R2/R3: equal position weights, free class sizes, all 4^n codons sense", "",
      "| quantity | usage | OLD | NEW | conclusion unchanged? |", "|---|:-:|--:|--:|---|"]
for a in R["a"]:
    u, nm = a["usage"], a["scenario"]
    o2, o3 = ucRow(nm, 2), ucRow(nm, 3)
    oadv = o2["best_Dw"] - o3["best_Dw"]
    n2, n3 = a["n2"], a["n3"]
    L += [f"| D2* | {u} | {o2['best_Dw']:.9f} | {a['D2']:.9f} | value changes (new metric) |",
          f"| D3* | {u} | {o3['best_Dw']:.9f} | {a['D3']:.9f} | value changes (new metric) |",
          f"| triplet advantage D2*−D3* | {u} | {oadv:.6f} ({100 * oadv / o2['best_Dw']:.3f}%) | "
          f"{a['advantage']:.6f} ({a['advantage_percent']:.3f}%) | "
          f"{'yes' if sign(oadv) == sign(a['advantage']) else 'NO'} (OLD sign {sign(oadv)}, NEW sign {sign(a['advantage'])}) |",
          f"| hits n=2 random | {u} | {o2['random']['hits']}/{o2['random']['restarts']} | "
          f"{n2['random']['hits']}/{n2['random']['restarts']} | — |",
          f"| hits n=3 random / warm | {u} | {o3['random']['hits']}/{o3['random']['restarts']} / "
          f"{o3['warm']['hits']}/{o3['warm']['restarts']} | {n3['random']['hits']}/{n3['random']['restarts']} / "
          f"{n3['warm']['hits']}/{n3['warm']['restarts']} | — |",
          f"| certified D3* ≤ lift | {u} | {o3['certified']} | {a['certified']} | "
          f"{'yes' if o3['certified'] == a['certified'] else 'NO'} |",
          f"| n=3 syn p1/p2/p3 (best code) | {u} | — | "
          f"{'/'.join(f'{x:.3f}' for x in a['n3_position_synonymy'])} | — |"]
L += ["", "Manuscript statements under test: \"D2* = D3* = 20.0175 (U)\" and \"18.836 vs 18.767 (P)\".", ""]

# ── (c) ──
c = R["c"]
oref = UCR["reference_codes"]
L += ["## (c) reference codes, equal weights", "",
      "| quantity | usage | OLD | NEW | conclusion unchanged? |", "|---|:-:|--:|--:|---|"]
for u in ("U", "P"):
    a = next((x for x in R["a"] if x["usage"] == u), None)
    osgc = oref["sgc_triplet_projection_stops_excluded"][f"Dw_{u}"]
    odbl = oref["doublet_projection"][f"Dw_{u}"]
    o3 = ucRow("R2" if u == "U" else "R3", 3)["best_Dw"]
    L += [f"| SGC 10-class projection (61 sense) | {u} | {osgc:.6f} | {c['sgc_triplet_projection_stops_excluded'][u]:.6f} | — |",
          f"| doublet projection (9 classes) | {u} | {odbl:.6f} | {c['doublet_projection'][u]:.6f} | — |",
          f"| SGC − doublet projection | {u} | {osgc - odbl:+.6f} | "
          f"{c['sgc_triplet_projection_stops_excluded'][u] - c['doublet_projection'][u]:+.6f} | "
          f"{'yes' if sign(osgc - odbl) == sign(c['sgc_triplet_projection_stops_excluded'][u] - c['doublet_projection'][u]) else 'NO'} |"]
    if a:
        ex = c[f"sgc_excess_over_n3_opt_{u}"]
        L += [f"| SGC excess over n=3 optimum | {u} | {osgc - o3:.6f} ({100 * (osgc - o3) / o3:.2f}%) | "
              f"{ex['abs']:.6f} ({ex['percent']:.2f}%) | {'yes' if sign(ex['abs']) == sign(osgc - o3) else 'NO'} (SGC above optimum) |"]
L += [""]

# ── (b) ──
L += ["## (b) position_empirical grid", "",
      "H-inv doublet = (w1, w2); H-term = (w1, w3); A = 1 − D3*/D2*. Conjecture: D3* = min_p D2*(w minus p) to 1e-9. "
      "\"Unchanged\" means the same silenced position, the same conjecture verdict, and the same sign of A_inv and A_term.", "",
      "| set | usage | | D2*(H-inv) | D2*(H-term) | D3* | silenced | A_inv | A_term | min_p D2* − D3* | conj. | cert. | SGC Dpos | unchanged? |",
      "|---|:-:|:-:|--:|--:|--:|:-:|--:|--:|--:|:-:|:-:|--:|---|"]
for r in R["b_n3"]:
    o = next(x for x in PE["n3"] if x["set"] == r["set"] and x["usage"] == r["usage"])
    osyn = o["silenced_position"]
    same = (str(osyn) == str(r["silenced_position"]) and o["conjecture"]["equal"] == r["conjecture"]["equal"]
            and sign(o["A_inv"], 1e-9) == sign(r["A_inv"], 1e-9) and sign(o["A_term"], 1e-9) == sign(r["A_term"], 1e-9))
    for tag, x, sil in (("OLD", o, osyn), ("NEW", r, r["silenced_position"])):
        L.append(f"| {r['set']} {tuple(r['w'])} | {r['usage']} | {tag} | {x['D2_Hinv']:.6f} | {x['D2_Hterm']:.6f} | "
                 f"{x['D3']:.6f} | {sil} | {x['A_inv']:+.5f} | {x['A_term']:+.5f} | {x['conjecture']['gap']:.2e} | "
                 f"{'EQUAL' if x['conjecture']['equal'] else 'FALSE'} | {'yes' if x['certified'] else 'NO'} | "
                 f"{x['sgc_Dpos']:.6f} | {('yes' if same else 'NO') if tag == 'NEW' else ''} |")
L += ["", "Hits (NEW, n=3): " + "; ".join(
    f"{r['usage']} {r['set']} rand {r['random']['hits']}/{r['random']['restarts']} "
    + " ".join(f"{k[-2:]}:{v['hits']}/{v['restarts']}" for k, v in r["warm"].items()) for r in R["b_n3"]), ""]

L += ["## Restarts and reductions", "",
      f"(a): n=2 {R['config']['a'][0]} random; n=3 {R['config']['a'][1]} random + {R['config']['a'][2]} silent-lift warm "
      f"(usage_cert counts). (b): n=2 {R['config']['b_n2'][0]} random + {R['config']['b_n2'][1]} warm from the previous "
      f"n=2 optimum; n=3 {R['config']['b_n3'][0]} random + {R['config']['b_n3'][1]} warm split over the three silent lifts "
      "(ceil → 6 per lift, 18 total). This is the stated floor, **below position_empirical's 256 + 32 (n=2: 256)**, "
      "chosen up front because the machine was shared (load ≈ 88 on 64 cores). "
      f"Further reductions: {R['reductions'] or 'none'}. Elapsed {R.get('elapsed_seconds')} s"
      f"{'' if R.get('finished') else ' — **RUN NOT FINISHED; cells missing above were not run**'}.", "",
      "Manuscript Table 2 could not be located in `manuscript_JME.md` (no D2*/D3* table there); the OLD column "
      "for (b) is `position_empirical/RESULTS.json`, which is the source those values were drawn from.", ""]
open(os.path.join(HERE, "REPORT.md"), "w").write("\n".join(L))
