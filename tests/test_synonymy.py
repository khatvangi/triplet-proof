"""
test_synonymy.py — verify position-specific synonymy numbers used in manuscript.

these are the ground-truth values (all-subs denominator: includes sense→stop
in denominator, excludes them from numerator). the manuscript must match these.
"""
import os, sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

import pytest
from verify_synonymy import (
    compute_synonymy, build_doublet_table,
    RNA_BASES, RAA10,
)
from src.io.codon_io import load_codon_table_csv
from run_publication_controls import build_triplet_10aa


@pytest.fixture(scope="module")
def triplet_df():
    return load_codon_table_csv()


def _all_subs_synonymy(codon_df):
    """compute synonymy with all-subs denominator (matching manuscript convention)."""
    from collections import Counter
    sense = codon_df[~codon_df["is_stop"].astype(bool)]
    c2a = dict(zip(sense["codon"].str.upper(), sense["aa"].str.upper()))

    all_codons = {}
    for _, row in codon_df.iterrows():
        c = str(row["codon"]).strip().upper()
        is_stop = str(row.get("is_stop", "")).strip().upper() in ("TRUE", "1", "YES", "*")
        all_codons[c] = "*" if is_stop else str(row["aa"]).strip().upper()

    codon_len = len(list(c2a.keys())[0])
    syn = [0] * codon_len
    total = [0] * codon_len  # all subs (including sense→stop)

    for codon, aa in c2a.items():
        for pos in range(codon_len):
            for base in RNA_BASES:
                if base == codon[pos]:
                    continue
                mut = codon[:pos] + base + codon[pos + 1:]
                mut_label = all_codons.get(mut)
                if mut_label is None:
                    continue
                total[pos] += 1
                if mut_label == aa:
                    syn[pos] += 1

    return {f"pos{i+1}": syn[i] / total[i] if total[i] > 0 else 0 for i in range(codon_len)}


def test_synonymy_20aa_triplet(triplet_df):
    """condition B: standard genetic code (20 AA, n=3)."""
    result = _all_subs_synonymy(triplet_df)
    # pos 1: 4.4%, pos 2: 0.0%, pos 3: 68.9%
    assert abs(result["pos1"] - 0.0437) < 0.005, f"pos1={result['pos1']:.4f}"
    assert result["pos2"] == 0.0, f"pos2={result['pos2']:.4f}, expected 0.0"
    assert abs(result["pos3"] - 0.6885) < 0.005, f"pos3={result['pos3']:.4f}"


def test_synonymy_10aa_triplet(triplet_df):
    """condition C: reduced 10-AA alphabet on triplet codons."""
    triplet_10aa_df = build_triplet_10aa(triplet_df, RAA10)
    result = _all_subs_synonymy(triplet_10aa_df)
    # pos 1: 26.2%, pos 2: 10.9%, pos 3: 76.5%
    assert abs(result["pos1"] - 0.262) < 0.005, f"pos1={result['pos1']:.4f}"
    assert abs(result["pos2"] - 0.109) < 0.005, f"pos2={result['pos2']:.4f}"
    assert abs(result["pos3"] - 0.765) < 0.005, f"pos3={result['pos3']:.4f}"


def test_synonymy_doublet(triplet_df):
    """condition A: doublet projection (10 AA, n=2)."""
    doublet_df = build_doublet_table(triplet_df, RAA10)
    result = _all_subs_synonymy(doublet_df)
    # pos 1: 33.3%, pos 2: 8.3%
    assert abs(result["pos1"] - 0.333) < 0.005, f"pos1={result['pos1']:.4f}"
    assert abs(result["pos2"] - 0.083) < 0.005, f"pos2={result['pos2']:.4f}"


def test_doublet_class_count(triplet_df):
    """doublet has 9 distinct classes, not 10 (X8/His lost to tie-breaking)."""
    doublet_df = build_doublet_table(triplet_df, RAA10)
    sense = doublet_df[~doublet_df["is_stop"].astype(bool)]
    n_classes = len(set(sense["aa"].str.upper()))
    assert n_classes == 9, f"expected 9 classes, got {n_classes}"
