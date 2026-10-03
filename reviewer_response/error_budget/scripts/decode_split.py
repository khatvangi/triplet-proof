#!/usr/bin/env python
"""
PART 1 + PART 3 of the error budget: decoding (misreading) vs charging (misacylation).

part 1: every landerer 2024 data s4 psm (codon c -> destination class b) is classified
  M1  = some codon of b is hamming-1 from c (one codon-anticodon mismatch explains it)
  M2+ = no codon of b is hamming-1 from c (candidate misacylation, or >=2 mismatches)
  i and l are one mass class ('I/L'); i<->l is invisible, stops are excluded.
  d(a,b) for destination 'I/L': the member that is hamming-1 from c if exactly one is,
  otherwise the mean of d(a,I) and d(a,L). a sensitivity drops all 'I/L' destinations.
part 3: sgc mean d over non-synonymous single-nucleotide sense neighbours, uniform and
  katoh (0.98, 1, 1.39) position weights, vs codon-shuffle null, block-preserving null
  (20 aa permuted over the sgc's 20 synonymous codon sets, ser = one 6-codon set), and the
  best block-preserving code in 1e6 samples (+ a swap-descent refinement from it).

run: python reviewer_response/error_budget/scripts/decode_split.py  (from repo root or anywhere)
"""
import glob
import json
import os

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
EB = os.path.dirname(HERE)
RR = os.path.dirname(EB)
S4 = os.path.join(RR, "position_rates/raw/landerer2024/S4")
OUT = os.path.join(EB, "part1_part3_results.json")

B = "TCAG"
AA = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODE = {a + b + c: AA[16 * i + 4 * j + k]
        for i, a in enumerate(B) for j, b in enumerate(B) for k, c in enumerate(B)}
AAS = sorted(set(AA) - {"*"})
KATOH = (0.98, 1.0, 1.39)  # katoh & suga 2023 supp fig 8, strict, from position_rates/POSITION_RATES.json
SEED = 20260923
NSAMP = 1_000_000


# ── distances ──

def loadDist():
    # euclidean distance between pca1..pca8 vectors z-scored across the 20 aa (as prior runs)
    t = pd.read_csv(os.path.join(RR, "aa_props_export.csv")).set_index("aa").loc[AAS]
    X = t[[f"pca{i}" for i in range(1, 9)]].to_numpy(float)
    X = (X - X.mean(0)) / X.std(0, ddof=0)
    D = np.sqrt(((X[:, None] - X[None]) ** 2).sum(-1))
    return pd.DataFrame(D, index=AAS, columns=AAS)


def ham1(c):
    return [c[:k] + b + c[k + 1:] for k in range(3) for b in B if b != c[k]]


def cls(a):
    return "I/L" if a in ("I", "L") else a


# ── part 1 ──

def classify(c, dest, D):
    o = CODE[c]
    near = {CODE[n] for n in ham1(c)}
    members = ["I", "L"] if dest == "I/L" else [dest]
    hit = [m for m in members if m in near]
    kind = "M1" if hit else "M2+"
    use = hit if len(hit) == 1 else members
    return kind, float(np.mean([D.loc[o, m] for m in use]))


def load(org):
    base = f"{S4}/{org}_filtered"
    ev, ex = [], []
    for f in sorted(glob.glob(f"{base}/*_substitution_errors.csv")):
        ds = os.path.basename(f).split("_")[0]
        s = pd.read_csv(f, usecols=["codon", "origin", "destination"])
        ev.append(s.assign(ds=ds))
        c = pd.read_csv(f"{base}/{ds}_codon_counts.csv", usecols=["codon", "base_count", "error_count"])
        assert c.error_count.sum() == len(s), ds
        ex.append(c.assign(ds=ds))
    return pd.concat(ev, ignore_index=True), pd.concat(ex, ignore_index=True)


def classStats(ev, expo_total):
    out = {}
    tot_d = ev.d.sum()
    for k in ("M1", "M2+"):
        e = ev[ev.kind == k]
        out[k] = dict(n_events_psm=int(len(e)),
                      summed_intensity="UNAVAILABLE (S4 substitution files carry no intensity/abundance column; PSM counts only)",
                      event_share=float(len(e) / len(ev)) if len(ev) else None,
                      mean_d=float(e.d.mean()) if len(e) else None,
                      median_d=float(e.d.median()) if len(e) else None,
                      sum_d=float(e.d.sum()),
                      cost_share=float(e.d.sum() / tot_d) if tot_d else None,
                      rate_per_codon_exposure=float(len(e) / expo_total) if expo_total else None,
                      d_per_codon_exposure=float(e.d.sum() / expo_total) if expo_total else None)
    return out


def part1(D, known_pairs):
    res = {}
    for org in ("ecoli", "yeast"):
        ev, ex = load(org)
        ok = ev.origin == ev.codon.map(CODE)
        n_bad = int((~ok).sum())
        ev = ev[ok].copy()
        ev = ev[ev.codon.map(CODE) != "*"]
        ev = ev[~((ev.origin.isin(["I", "L"])) & (ev.destination == "I/L"))]
        kd = [classify(c, b, D) for c, b in zip(ev.codon, ev.destination)]
        ev["kind"] = [k for k, _ in kd]
        ev["d"] = [d for _, d in kd]
        ex = ex[ex.codon.map(CODE) != "*"]
        expo_ds = ex.groupby("ds").base_count.sum()
        expo = int(expo_ds.sum())
        r = dict(n_psm_dropped_codon_origin_mismatch=n_bad, n_events=int(len(ev)), n_datasets=int(ev.ds.nunique()),
                 codon_exposure_total=expo, pooled=classStats(ev, expo))
        r["no_IL_destination"] = classStats(ev[ev.destination != "I/L"], expo)
        r["no_cys_origin"] = classStats(ev[ev.origin != "C"], int(ex[ex.codon.map(CODE) != "C"].base_count.sum()))
        # per dataset
        per = []
        for ds, g in ev.groupby("ds"):
            s = classStats(g, int(expo_ds.get(ds, 0)))
            per.append(dict(ds=ds, n=len(g), **{f"{k}_{f}": s[k][f] for k in ("M1", "M2+")
                                                  for f in ("n_events_psm", "mean_d", "median_d", "cost_share", "rate_per_codon_exposure")}))
        per = pd.DataFrame(per)
        per.to_csv(os.path.join(EB, f"part1_per_dataset_{org}.csv"), index=False)
        r["per_dataset_summary"] = dict(
            file=f"part1_per_dataset_{org}.csv",
            M2plus_cost_share_median=float(per["M2+_cost_share"].median()),
            M2plus_cost_share_iqr=[float(per["M2+_cost_share"].quantile(.25)), float(per["M2+_cost_share"].quantile(.75))],
            M2plus_event_share_median=float((per["M2+_n_events_psm"] / per.n).median()))
        # per codon rates by class
        bc = ex.groupby("codon").base_count.sum()
        pc = ev.groupby(["codon", "kind"]).size().unstack(fill_value=0).reindex(bc.index, fill_value=0)
        pc = pc.reindex(columns=["M1", "M2+"], fill_value=0)
        pc["base_count"] = bc
        pc["rate_M1"] = pc["M1"] / pc.base_count
        pc["rate_M2+"] = pc["M2+"] / pc.base_count
        pc.to_csv(os.path.join(EB, f"part1_per_codon_{org}.csv"))
        # top m2+ pairs
        m2 = ev[ev.kind == "M2+"]
        top = m2.groupby(["origin", "destination"]).agg(n=("d", "size"), d=("d", "mean")).sort_values("n", ascending=False).head(10)
        r["top10_M2plus_pairs"] = [dict(pair=f"{a}->{b}", n_psm=int(t.n), d=float(t.d),
                                        known_misacylation_pair=known_pairs.get((a, b), "no sourced row in PART 2"))
                                   for (a, b), t in top.iterrows()]
        res[org] = r
    return res


# ── part 3 ──

def sgcGraph():
    sense = [c for c in CODE if CODE[c] != "*"]
    idx = {c: i for i, c in enumerate(sense)}
    edges, pos = [], []
    for c in sense:
        for k in range(3):
            for b in B:
                if b == c[k]:
                    continue
                n = c[:k] + b + c[k + 1:]
                if CODE[n] != "*":
                    edges.append((idx[c], idx[n]))
                    pos.append(k)
    return sense, np.array(edges), np.array(pos)


def meanD(lab, E, w, Dm):
    # lab: (n, 61) aa indices. mean d over directed non-synonymous edges, edge weight w
    a, b = lab[:, E[:, 0]], lab[:, E[:, 1]]
    ns = a != b
    return (Dm[a, b] * w).sum(1) / (ns * w).sum(1)


def part3(D, rng):
    Dm = D.to_numpy()
    sense, E, pos = sgcGraph()
    sgc = np.array([AAS.index(CODE[c]) for c in sense])
    W = {"uniform": np.ones(len(E)), "katoh": np.array(KATOH)[pos]}
    blocks = [np.where(sgc == i)[0] for i in range(20)]
    out = {"n_sense": len(sense), "n_directed_sense_edges": int(len(E)),
           "n_nonsyn_directed_edges_sgc": int((sgc[E[:, 0]] != sgc[E[:, 1]]).sum())}
    CH = 50_000
    shuf, blk, best = {k: [] for k in W}, {k: [] for k in W}, {k: (np.inf, None) for k in W}
    for s in range(0, NSAMP, CH):
        # codon-shuffle null: permute the 61 sense labels (per-aa counts kept; aug not pinned)
        L = np.array([rng.permutation(sgc) for _ in range(CH)])
        # block null: permute the 20 aa over the sgc's 20 synonymous codon sets
        P = np.array([rng.permutation(20) for _ in range(CH)])
        LB = P[:, sgc]
        for k, w in W.items():
            shuf[k].append(meanD(L, E, w, Dm))
            v = meanD(LB, E, w, Dm)
            blk[k].append(v)
            i = int(v.argmin())
            if v[i] < best[k][0]:
                best[k] = (float(v[i]), P[i].copy())
    for k, w in W.items():
        sh, bl = np.concatenate(shuf[k]), np.concatenate(blk[k])
        v_sgc = float(meanD(sgc[None], E, w, Dm)[0])
        # swap descent from the best sampled block code (a refinement, not a certificate)
        p = best[k][1].copy()
        cur = best[k][0]
        improved = True
        while improved:
            improved = False
            for i in range(20):
                for j in range(i + 1, 20):
                    q = p.copy(); q[i], q[j] = q[j], q[i]
                    v = float(meanD(q[sgc][None], E, w, Dm)[0])
                    if v < cur - 1e-12:
                        p, cur, improved = q, v, True
        out[k] = dict(sgc_mean_d=v_sgc,
                      codon_shuffle=dict(median=float(np.median(sh)), mean=float(sh.mean()), q01=float(np.quantile(sh, .01)),
                                         frac_le_sgc=float((sh <= v_sgc).mean()), n=len(sh)),
                      block_null=dict(median=float(np.median(bl)), mean=float(bl.mean()), q01=float(np.quantile(bl, .01)),
                                      frac_le_sgc=float((bl <= v_sgc).mean()), n=len(bl)),
                      block_best_sampled=best[k][0],
                      block_best_sampled_code={AAS[b]: AAS[int(best[k][1][b])] for b in range(20)},
                      block_best_after_swap_descent=cur,
                      block_best_after_swap_descent_code={AAS[b]: AAS[int(p[b])] for b in range(20)})
    return out


def main():
    D = loadDist()
    iu = np.triu_indices(20, 1)
    known = {}
    p2 = os.path.join(EB, "part2_misacylation_table.json")
    if os.path.exists(p2):
        for r in json.load(open(p2)):
            if isinstance(r.get("rate_per_aminoacylation"), (int, float)):
                known.setdefault((r["cognate_aa"], r["noncognate_aa"]), f"yes: {r.get('enzyme')} ({r.get('doi')})")
            else:
                known.setdefault((r["cognate_aa"], r["noncognate_aa"]), "listed in PART 2 but rate UNAVAILABLE")
    res = dict(distance=dict(definition="euclidean on pca1..pca8 z-scored across 20 aa (aa_props_export.csv)",
                             mean_d_190_pairs=float(D.to_numpy()[iu].mean()),
                             median_d_190_pairs=float(np.median(D.to_numpy()[iu])),
                             matrix=D.round(6).to_dict()),
               part1=part1(D, known),
               part3=part3(D, np.random.default_rng(SEED)),
               seed=SEED, n_samples=NSAMP, katoh_weights=KATOH)
    json.dump(res, open(OUT, "w"), indent=1)
    print(json.dumps({k: v for k, v in res.items() if k != "distance"}, indent=1)[:6000])
    print("mean d 190:", res["distance"]["mean_d_190_pairs"])


if __name__ == "__main__":
    main()
