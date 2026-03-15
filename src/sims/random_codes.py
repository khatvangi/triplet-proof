# src/sims/random_codes.py
import random
import pandas as pd
from collections import Counter, defaultdict

def fully_random_code(df: pd.DataFrame, aa_list=None, keep_stops=True, seed=None):
    rng = random.Random(seed)
    sense = df[~df["is_stop"]].copy()
    aa_pool = aa_list or sorted(list(set([a for a in sense["aa"].dropna().unique() if a != "*"])))
    out = {}
    for _, r in df.iterrows():
        c = r["codon"]
        if keep_stops and r["is_stop"]:
            out[c] = None
        else:
            out[c] = rng.choice(aa_pool)
    return out

def degeneracy_preserving_shuffle(df: pd.DataFrame, seed=None):
    """
    Keep the number of codons per AA the same as the input df;
    shuffle AA labels among sense codons. Stops remain stops; AUG remains Met.
    """
    rng = random.Random(seed)
    sense = df[~df["is_stop"]].copy()

    # Preserve AUG->Met (start) to keep a realistic start context
    sense_nostart = sense[sense["codon"] != "AUG"].copy()
    labels = sense_nostart["aa"].tolist()
    rng.shuffle(labels)

    out = {}
    # first assign stops
    for _, r in df.iterrows():
        if r["is_stop"]:
            out[r["codon"]] = None

    # assign AUG -> Met
    if "AUG" in df["codon"].values:
        out["AUG"] = "M"

    # assign the rest
    i = 0
    for _, r in df.iterrows():
        c = r["codon"]
        if r["is_stop"] or c == "AUG":
            continue
        out[c] = labels[i]
        i += 1
    return out

