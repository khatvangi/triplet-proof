# src/io/codon_io.py
import pandas as pd
import os

def load_codon_table_csv(path: str = "codon_table.csv") -> pd.DataFrame:
    """
    Load codon→AA mapping from CSV.
    Required columns: codon, aa
    Optional: is_start, is_stop
    """
    if not os.path.exists(path):
        raise FileNotFoundError(f"codon_io: file not found: {path}")
    df = pd.read_csv(path)
    # basic cleanup
    df["codon"] = df["codon"].astype(str).str.upper().str.strip()
    df["aa"] = df["aa"].astype(str).str.strip().str.upper()

    if "is_stop" not in df.columns:
        df["is_stop"] = df["aa"].isin(["*", "STOP"])
    if "is_start" not in df.columns:
        df["is_start"] = df["codon"].eq("AUG")
    # Preserve explicit base columns when available, or derive them from codon.
    for i, col in enumerate(("b1", "b2", "b3"), start=0):
        if col not in df.columns:
            df[col] = df["codon"].str[i]

    base_cols = [c for c in ("b1", "b2", "b3") if c in df.columns]
    ordered = ["codon", "aa", "is_start", "is_stop"] + base_cols
    extra = [c for c in df.columns if c not in ordered]
    return df[ordered + extra]


def build_sgc_map(df: pd.DataFrame | None = None, path: str = "codon_table.csv") -> pd.DataFrame:
    """
    Convenience: if df is None, loads codon_table.csv from repo root.
    Returns a clean DataFrame with columns ['codon','aa','is_start','is_stop'].
    """
    if df is None:
        df = load_codon_table_csv(path)
    cols = ["codon", "aa", "is_start", "is_stop"]
    for c in cols:
        if c not in df.columns:
            raise ValueError(f"codon_io.build_sgc_map: missing column '{c}'")
    return df.copy()
