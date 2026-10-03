#!/usr/bin/env python
"""
codon-position assignment of the per-PSM misincorporation events in
Landerer, Poehls & Toth-Petroczy 2024 MBE msae048, Data S4
(raw/landerer2024/S4/{ecoli,yeast}_filtered/*_substitution_errors.csv),
with codon exposure from the matching *_codon_counts.csv (base_count).

this is our own assignment, not a number the authors report. rule:
  an event (source codon c, destination aa d) is assigned to position k if
  some codon of d differs from c at position k only (hamming-1, near-cognate).
  I and L are one class ('I/L') because the authors merge them (equal mass).
  if no hamming-1 codon of d exists -> 'multi' (needs >=2 mismatches).
  if hamming-1 codons of d exist at more than one position -> 'ambiguous'.
the '1st/2nd/3rd' here are mRNA codon positions (position 3 = wobble).

per-opportunity rate at position k:
  eps_k = N_k / sum_c B_c * m_k(c)
  N_k   = PSMs assigned to k
  B_c   = codon exposure (base_count summed over datasets, PSM units, same
          unit as the authors' Data S2 detection rate)
  m_k(c)= number of the 3 alternative bases at k that give a detectable
          non-synonymous sense change (stop excluded, I<->L excluded).
variant 'observed_pairs': m_k counts only destination classes that are seen
at least once from that origin aa in the organism (crude proxy for
detectability; PTM-masked pairs are removed by the authors upstream).
variant 'no_cys': drops cys-origin events and cys codons (C->S and C->A are
the two most common substitutions overall and may be sample-prep artefacts).
"""
import glob
import json
import os
import pandas as pd

B = "TCAG"
AA = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODE = {a + b + c: AA[16 * i + 4 * j + k]
        for i, a in enumerate(B) for j, b in enumerate(B) for k, c in enumerate(B)}


def cls(aa):
    # detection class: I and L are indistinguishable by mass
    return "I/L" if aa in ("I", "L") else aa


def neighbours(codon, k):
    # classes reachable by a single change at position k (0-based)
    return [cls(CODE[codon[:k] + b + codon[k + 1:]]) for b in B if b != codon[k]]


def assignPos(codon, dest):
    ks = [k + 1 for k in range(3) if dest in neighbours(codon, k)]
    if not ks:
        return "multi"
    return str(ks[0]) if len(ks) == 1 else "ambiguous:" + "/".join(map(str, ks))


def load(org):
    base = f"raw/landerer2024/S4/{org}_filtered"
    ev, ex = [], []
    for f in sorted(glob.glob(f"{base}/*_substitution_errors.csv")):
        ds = os.path.basename(f).split("_")[0]
        s = pd.read_csv(f, usecols=["Protein", "Localization_in_Protein", "codon", "origin", "destination"])
        ev.append(s.assign(ds=ds))
        c = pd.read_csv(f"{base}/{ds}_codon_counts.csv", usecols=["codon", "base_count", "error_count"])
        # the per-codon error_count must equal the psm rows in the substitution file
        assert c.error_count.sum() == len(s), ds
        ex.append(c.assign(ds=ds))
    return pd.concat(ev, ignore_index=True), pd.concat(ex, ignore_index=True)


def rates(ev, expo, pairs=None):
    out = {}
    npos = {k: int((ev.pos == k).sum()) for k in "123"}
    for k in "123":
        opp = 0
        for codon, bc in expo.items():
            o = cls(CODE[codon])
            m = [d for d in neighbours(codon, int(k) - 1) if d not in ("*", o)]
            if pairs is not None:
                m = [d for d in m if (o, d) in pairs]
            opp += bc * len(m)
        out[k] = dict(n_psm=npos[k], opportunities=int(opp), eps=npos[k] / opp if opp else None)
    out["eps1_over_eps2"] = out["1"]["eps"] / out["2"]["eps"]
    out["eps3_over_eps2"] = out["3"]["eps"] / out["2"]["eps"]
    return out


def summarise(ev, ex, label):
    # drop rows whose codon does not encode the stated origin aa (source
    # annotation inconsistency; 84 e. coli psm rows, 0 yeast)
    ok = ev.origin == ev.codon.map(CODE)
    n_bad = int((~ok).sum())
    ev = ev[ok].copy()
    ev["pos"] = [assignPos(c, d) for c, d in zip(ev.codon, ev.destination)]
    sites = ev.drop_duplicates(["ds", "Protein", "Localization_in_Protein", "destination"])
    # sense codons only (e. coli has one TGA occurrence with base_count 1, a selenoprotein)
    expo = ex[ex.codon.map(CODE) != "*"].groupby("codon").base_count.sum().to_dict()
    pairs = {(cls(o), d) for o, d in zip(ev.origin, ev.destination)}
    single_psm = ev[ev.pos.isin(list("123"))]
    single_site = sites[sites.pos.isin(list("123"))]
    res = dict(
        label=label,
        n_psm_dropped_codon_origin_mismatch=n_bad,
        n_psm=len(ev), n_unique_sites=len(sites), n_datasets=int(ev.ds.nunique()),
        class_counts_psm=ev.pos.value_counts().to_dict(),
        class_counts_sites=sites.pos.value_counts().to_dict(),
        share_single_mismatch_psm={k: float((single_psm.pos == k).mean()) for k in "123"},
        share_single_mismatch_sites={k: float((single_site.pos == k).mean()) for k in "123"},
        per_opportunity_all_pairs=rates(ev, expo),
        per_opportunity_observed_pairs=rates(ev, expo, pairs),
        top_pairs_by_position={k: single_psm[single_psm.pos == k]
                               .assign(p=lambda x: x.codon + ">" + x.destination)
                               .p.value_counts().head(8).to_dict() for k in "123"},
    )
    return res


def main():
    out = {}
    for org in ("ecoli", "yeast"):
        ev, ex = load(org)
        out[org] = summarise(ev, ex, f"{org}: all events")
        out[org + "_no_cys"] = summarise(ev[ev.origin != "C"], ex[ex.codon.map(CODE) != "C"],
                                         f"{org}: cys-origin events and cys codons removed")
    json.dump(out, open("landerer2024_position_summary.json", "w"), indent=2)
    for key, r in out.items():
        print(f"== {r['label']}  dropped={r['n_psm_dropped_codon_origin_mismatch']} psm={r['n_psm']} sites={r['n_unique_sites']} datasets={r['n_datasets']}")
        print("  classes psm  :", r["class_counts_psm"])
        print("  classes sites:", r["class_counts_sites"])
        print("  single-mismatch share psm  :", {k: round(v, 3) for k, v in r["share_single_mismatch_psm"].items()})
        print("  single-mismatch share sites:", {k: round(v, 3) for k, v in r["share_single_mismatch_sites"].items()})
        for v in ("per_opportunity_all_pairs", "per_opportunity_observed_pairs"):
            q = r[v]
            print(f"  {v}: " + " ".join(f"eps{k}={q[k]['eps']:.3e}(N={q[k]['n_psm']})" for k in "123")
                  + f"  eps1/eps2={q['eps1_over_eps2']:.3f} eps3/eps2={q['eps3_over_eps2']:.3f}")
        print("  top pairs:", r["top_pairs_by_position"])


if __name__ == "__main__":
    main()
