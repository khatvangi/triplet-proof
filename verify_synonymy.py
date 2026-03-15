#!/usr/bin/env python3
"""
verify_synonymy.py — compute position-specific synonymy for all three conditions.
resolves the discrepancy between manuscript versions.

two denominator conventions:
  "all-subs":       synonymous / (all substitutions including sense→stop)
  "sense-to-sense": synonymous / (sense→sense substitutions only)
"""
import os, sys
sys.path.insert(0, os.path.abspath(os.path.dirname(__file__)))

import pandas as pd
from collections import Counter

from src.io.codon_io import load_codon_table_csv

RNA_BASES = ("A", "U", "G", "C")

RAA10 = {
    "A": "X1", "G": "X1",
    "V": "X2", "L": "X2", "I": "X2", "M": "X2",
    "F": "X3", "W": "X3", "Y": "X3",
    "S": "X4", "T": "X4",
    "N": "X5", "Q": "X5",
    "D": "X6", "E": "X6",
    "K": "X7", "R": "X7",
    "H": "X8",
    "P": "X9",
    "C": "X10",
}


def majority_vote(labels):
    cnt = Counter(labels)
    maxc = max(cnt.values())
    return sorted(lab for lab, c in cnt.items() if c == maxc)[0]


def build_doublet_table(triplet_df, raa_map):
    """project triplet SGC -> doublet via majority vote on pos 3."""
    t2c = {}
    for _, row in triplet_df.iterrows():
        aa = str(row["aa"]).strip().upper()
        codon = str(row["codon"]).strip().upper()
        is_stop = str(row.get("is_stop", "")).strip().upper() in ("TRUE", "1", "YES", "*")
        t2c[codon] = "*" if (is_stop or aa in ("*", "STOP")) else raa_map.get(aa, aa)
    rows = []
    for b1 in RNA_BASES:
        for b2 in RNA_BASES:
            pf = b1 + b2
            classes = [t2c.get(pf + b3, "*") for b3 in RNA_BASES]
            non_stop = [c for c in classes if c != "*"]
            cls = majority_vote(non_stop) if non_stop else "*"
            rows.append({"codon": pf, "aa": cls, "is_stop": cls == "*"})
    return pd.DataFrame(rows)


def compute_synonymy(codon_df, label="code"):
    """compute position-specific synonymy for a codon table.

    reports BOTH denominator conventions:
      all-subs: synonymous / total (including sense→stop as non-synonymous)
      sense-only: synonymous / sense-to-sense only (excluding sense→stop from denom)
    """
    sense = codon_df[~codon_df["is_stop"].astype(bool)]
    c2a = dict(zip(sense["codon"].str.upper(), sense["aa"].str.upper()))

    # build full codon→aa map including stops
    all_codons = {}
    for _, row in codon_df.iterrows():
        c = str(row["codon"]).strip().upper()
        is_stop = str(row.get("is_stop", "")).strip().upper() in ("TRUE", "1", "YES", "*")
        if is_stop:
            all_codons[c] = "*"
        else:
            all_codons[c] = str(row["aa"]).strip().upper()

    codon_len = len(list(c2a.keys())[0])

    syn = [0] * codon_len
    nonsyn_sense = [0] * codon_len  # sense→different-sense
    to_stop = [0] * codon_len       # sense→stop

    for codon, aa in c2a.items():
        for pos in range(codon_len):
            for base in RNA_BASES:
                if base == codon[pos]:
                    continue
                mut = codon[:pos] + base + codon[pos + 1:]
                mut_label = all_codons.get(mut)
                if mut_label is None:
                    # codon not in table (shouldn't happen for standard bases)
                    continue
                if mut_label == "*":
                    to_stop[pos] += 1
                elif mut_label == aa:
                    syn[pos] += 1
                else:
                    nonsyn_sense[pos] += 1

    print(f"\n{'=' * 65}")
    print(f"  {label}")
    print(f"  {len(c2a)} sense codons, codon length = {codon_len}")
    n_classes = len(set(c2a.values()))
    print(f"  {n_classes} distinct classes")
    print(f"{'=' * 65}")

    print(f"\n  {'pos':<6} {'syn':>5} {'ns':>5} {'→stop':>6} "
          f"{'all-subs%':>10} {'sense-only%':>12}")
    print(f"  {'-' * 50}")

    total_syn = 0
    total_ns = 0
    total_stop = 0

    for pos in range(codon_len):
        s = syn[pos]
        ns = nonsyn_sense[pos]
        st = to_stop[pos]
        total_all = s + ns + st
        total_sense = s + ns

        frac_all = s / total_all if total_all > 0 else 0
        frac_sense = s / total_sense if total_sense > 0 else 0

        print(f"  {pos + 1:<6} {s:>5} {ns:>5} {st:>6} "
              f"{frac_all * 100:>9.1f}% {frac_sense * 100:>11.1f}%")

        total_syn += s
        total_ns += ns
        total_stop += st

    grand_all = total_syn + total_ns + total_stop
    grand_sense = total_syn + total_ns
    frac_all = total_syn / grand_all if grand_all > 0 else 0
    frac_sense = total_syn / grand_sense if grand_sense > 0 else 0

    print(f"  {'-' * 50}")
    print(f"  {'total':<6} {total_syn:>5} {total_ns:>5} {total_stop:>6} "
          f"{frac_all * 100:>9.1f}% {frac_sense * 100:>11.1f}%")


def main():
    triplet_df = load_codon_table_csv()

    # condition B: 20 AA, n=3 (standard SGC)
    compute_synonymy(triplet_df, "CONDITION B: standard genetic code (20 AA, n=3)")

    # condition C: 10 AA, n=3 (reduced alphabet on triplet)
    from run_publication_controls import build_triplet_10aa
    triplet_10aa_df = build_triplet_10aa(triplet_df, RAA10)
    compute_synonymy(triplet_10aa_df, "CONDITION C: 10-AA triplet (reduced alphabet, n=3)")

    # condition A: 10 AA, n=2 (doublet)
    doublet_df = build_doublet_table(triplet_df, RAA10)
    compute_synonymy(doublet_df, "CONDITION A: doublet projection (10 AA, n=2)")


if __name__ == "__main__":
    main()
