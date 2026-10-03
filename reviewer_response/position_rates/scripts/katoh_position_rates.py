#!/usr/bin/env python
"""
codon-position-resolved misdecoding from Katoh & Suga 2023 Supp Fig 8
(parsed by parse_katoh_s8.py -> katoh2023_s8_parsed.csv).

only combinations whose mismatch position is uniquely determined can be
assigned to a position. two definitions are reported:
  strict  : exactly one '+' among 1st/2nd/3rd and no '±' (from the position rows)
  authors : total == '1' (the authors' own class, Fig 4 row "1"),
            position = the row carrying '+'
per position: n tested combos, n with detected misincorporation (intensity>0),
frequency = positives / tested (a per-opportunity rate), mean intensity
(0-5 semi-quantitative), and share of positive events.
"""
import json
import math
import pandas as pd


def wilson(k, n, z=1.96):
    # 95% wilson score interval for a binomial proportion
    if n == 0:
        return None
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return [c - h, c + h]

df = pd.read_csv("katoh2023_s8_parsed.csv", dtype=str, keep_default_na=False)
df["intensity"] = df.intensity.astype(int)
P = ["p1", "p2", "p3"]
df["n_plus"] = (df[P] == "+").sum(axis=1)
df["n_pm"] = (df[P] == "±").sum(axis=1)


def posOf(r):
    hits = [i + 1 for i, p in enumerate(P) if r[p] == "+"]
    return hits[0] if len(hits) == 1 else None


def summarise(sub, label):
    sub = sub.assign(pos=sub.apply(posOf, axis=1))
    unassigned = int(sub.pos.isna().sum())
    sub = sub.dropna(subset=["pos"])
    out = {"definition": label, "n_unassignable_in_class": unassigned, "positions": {}}
    npos_all = int((sub.intensity > 0).sum())
    for k in (1, 2, 3):
        g = sub[sub.pos == k]
        n, npos = len(g), int((g.intensity > 0).sum())
        out["positions"][k] = dict(
            n_tested=n, n_misdecoded=npos,
            freq=npos / n if n else None,
            freq_wilson95=wilson(npos, n),
            mean_intensity=float(g.intensity.mean()) if n else None,
            sum_intensity=int(g.intensity.sum()),
            share_of_events=npos / npos_all if npos_all else None,
            intensity_counts={int(a): int(b) for a, b in g.intensity.value_counts().sort_index().items()},
        )
    f = {k: out["positions"][k]["freq"] for k in (1, 2, 3)}
    m = {k: out["positions"][k]["mean_intensity"] for k in (1, 2, 3)}
    # katz log-ratio 95% interval for a ratio of two binomial proportions
    def katz(k):
        a, n1 = out["positions"][k]["n_misdecoded"], out["positions"][k]["n_tested"]
        b, n2 = out["positions"][2]["n_misdecoded"], out["positions"][2]["n_tested"]
        se = math.sqrt(1 / a - 1 / n1 + 1 / b - 1 / n2)
        r = (a / n1) / (b / n2)
        return [r * math.exp(-1.96 * se), r * math.exp(1.96 * se)]
    out["ratio_freq_eps1_eps2_katz95"] = katz(1)
    out["ratio_freq_eps3_eps2_katz95"] = katz(3)
    out["ratio_freq_eps1_eps2"] = f[1] / f[2]
    out["ratio_freq_eps3_eps2"] = f[3] / f[2]
    out["ratio_meanint_eps1_eps2"] = m[1] / m[2]
    out["ratio_meanint_eps3_eps2"] = m[3] / m[2]
    return out


strict = summarise(df[(df.n_plus == 1) & (df.n_pm == 0)], "strict: one '+' and no '±' in position rows")
authors = summarise(df[df.total == "1"], "authors: total=='1'")

# rows where the printed total disagrees with the printed + rows in a way
# that cannot come from correlated '±' across isoacceptors
flag = df[((df.n_pm == 0) & (df.n_plus.astype(str) != df.total))]
res = dict(strict=strict, authors=authors,
           source_inconsistencies=flag[["codon", "aa", "intensity", "p1", "p2", "p3", "total"]].to_dict("records"))

# per-position split by codon second base, since the paper's rules hinge on it
bysecond = {}
s = df[(df.n_plus == 1) & (df.n_pm == 0)].copy()
s["pos"] = s.apply(posOf, axis=1)
for b, g in s.groupby(s.codon.str[1]):
    bysecond[b] = {int(k): dict(n_tested=len(h), n_misdecoded=int((h.intensity > 0).sum()))
                   for k, h in g.groupby("pos")}
res["strict_by_codon_second_base"] = bysecond

json.dump(res, open("katoh2023_position_summary.json", "w"), indent=2)
for r in (strict, authors):
    print(r["definition"], "| unassignable:", r["n_unassignable_in_class"])
    for k, v in r["positions"].items():
        print(f"  pos{k}: tested={v['n_tested']} pos={v['n_misdecoded']} freq={v['freq']:.3f} "
              f"meanI={v['mean_intensity']:.3f} share={v['share_of_events']:.3f} {v['intensity_counts']}")
    print(f"  eps1/eps2 freq={r['ratio_freq_eps1_eps2']:.3f} eps3/eps2 freq={r['ratio_freq_eps3_eps2']:.3f} "
          f"| meanI eps1/eps2={r['ratio_meanint_eps1_eps2']:.3f} eps3/eps2={r['ratio_meanint_eps3_eps2']:.3f}")
for r in (strict, authors):
    print(r["definition"], {k: [round(x, 3) for x in v["freq_wilson95"]] for k, v in r["positions"].items()})
for r in (strict, authors):
    print(r["definition"], "katz95 eps1/eps2", [round(x, 3) for x in r["ratio_freq_eps1_eps2_katz95"]],
          "eps3/eps2", [round(x, 3) for x in r["ratio_freq_eps3_eps2_katz95"]])
print("source inconsistencies:", res["source_inconsistencies"])
print("by 2nd base:", json.dumps(bysecond))
