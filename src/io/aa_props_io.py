# src/io/aa_props_io.py
from __future__ import annotations
from typing import Dict, Iterable, List, Optional, Tuple
import numpy as np
import pandas as pd
import os

# ------------------------- internals -------------------------

_NUMERIC_EXCLUDE = {"aa"}

def _numeric_cols(df: pd.DataFrame) -> List[str]:
    return [c for c in df.columns
            if c not in _NUMERIC_EXCLUDE and pd.api.types.is_numeric_dtype(df[c])]

def _zscore(X: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    mu = X.mean(axis=0, keepdims=True)
    sd = X.std(axis=0, ddof=0, keepdims=True)
    sd[sd == 0.0] = 1.0
    return (X - mu) / sd, mu, sd

def _pca_by_svd(Xz: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    # Xz: n×d z-scored. SVD: Xz = U S V^T
    U, S, Vt = np.linalg.svd(Xz, full_matrices=False)
    # components are columns of V; project: Z = Xz @ V
    var = (S**2) / (Xz.shape[0] - 1)
    return Vt.T, var  # components (d×d), variances (d,)

def _drop_highly_correlated(df: pd.DataFrame, cols: List[str], thresh: float) -> List[str]:
    if len(cols) <= 1:
        return cols
    C = df[cols].corr().abs()
    keep = []
    dropped = set()
    for i, c in enumerate(cols):
        if c in dropped:
            continue
        keep.append(c)
        # drop any later column with high corr to c
        for j in range(i+1, len(cols)):
            cj = cols[j]
            if cj not in dropped and C.iloc[i, j] >= thresh:
                dropped.add(cj)
    return keep

def _pca_columns(df: pd.DataFrame) -> List[str]:
    pcs = [c for c in df.columns if c.lower().startswith("pca")]
    if not pcs:
        return []
    # sort pca1, pca2, ...
    def _key(c: str):
        s = ''.join(ch for ch in c if ch.isdigit())
        return (int(s) if s.isdigit() else 10**6, c.lower())
    return sorted(pcs, key=_key)

# ------------------------- public API -------------------------

def load_aa_props(
    path: Optional[str] = None,
    *,
    auto_select: bool = False,
    corr_thresh: float = 0.95,
    pca_var: Optional[float] = None,
    standardize_before_pca: bool = True
) -> Tuple[pd.DataFrame, Dict]:
    """
    Load AA properties and (optionally) auto-select + PCA compress.

    Parameters
    ----------
    path : str or None
        Defaults to './aa_props.csv' if None.
    auto_select : bool
        If True, drop highly correlated raw features (|r| >= corr_thresh).
    corr_thresh : float
        Correlation threshold for dropping redundant features.
    pca_var : float or None
        If set (e.g., 0.97), create pca1..pcaK columns explaining >= pca_var variance.
    standardize_before_pca : bool
        Z-score features before PCA.

    Returns
    -------
    (df, meta)
      df   : DataFrame with 'aa' and selected feature columns (raw and/or pca*)
      meta : dict with keys: 'raw_used', 'pca_k', 'var_explained', 'z_mu', 'z_sd'
    """
    if path is None:
        path = "aa_props.csv"
    if not os.path.exists(path):
        raise FileNotFoundError(f"aa_props_io: file not found: {path}")

    # read CSV/Parquet
    df = pd.read_parquet(path) if path.endswith(".parquet") else pd.read_csv(path)
    if "aa" not in df.columns:
        raise ValueError("aa_props_io: missing 'aa' column.")
    df = df.copy()
    df["aa"] = df["aa"].astype(str).str.strip().str.upper()

    # work on numeric raw features
    raw_cols = _numeric_cols(df)
    if not raw_cols:
        raise ValueError("aa_props_io: no numeric columns found (besides 'aa').")

    # optional correlation-based pruning
    used_raw = raw_cols
    if auto_select and len(raw_cols) > 1:
        used_raw = _drop_highly_correlated(df, raw_cols, corr_thresh)

    # Optionally build PCA columns
    meta = {
        "raw_used": used_raw,
        "pca_k": None,
        "var_explained": None,
        "z_mu": None,
        "z_sd": None,
    }

    if pca_var is not None:
        sub = df[["aa"] + used_raw].drop_duplicates("aa").set_index("aa").sort_index()
        X = sub.to_numpy(dtype=float)
        if standardize_before_pca:
            Xz, mu, sd = _zscore(X)
        else:
            Xz, mu, sd = X, np.zeros((1, X.shape[1])), np.ones((1, X.shape[1]))
        comps, var = _pca_by_svd(Xz)  # comps: d×d, var: d,
        total = var.sum()
        cum = np.cumsum(var / total)
        k = int(np.searchsorted(cum, pca_var) + 1)
        # project to k PCs
        Z = Xz @ comps[:, :k]
        # attach pca columns (aligned to sorted AAs)
        pca_cols = [f"pca{i+1}" for i in range(k)]
        pca_df = pd.DataFrame(Z, index=sub.index, columns=pca_cols).reset_index()
        # merge back into df on 'aa'
        df = df.merge(pca_df, on="aa", how="left")

        meta.update({
            "pca_k": k,
            "var_explained": float(cum[k-1]),
            "z_mu": mu.ravel().tolist(),
            "z_sd": sd.ravel().tolist(),
        })

    return df, meta


def prop_names(df: pd.DataFrame) -> List[str]:
    """Preferred feature names for downstream: PCA if present else numeric raw."""
    pcs = _pca_columns(df)
    if pcs:
        return pcs
    raw = _numeric_cols(df)
    if not raw:
        raise ValueError("aa_props_io: no usable numeric columns.")
    return raw


def aa_feature_columns(df: pd.DataFrame) -> List[str]:
    """Alias kept for backward compatibility."""
    return prop_names(df)


def aa_prop_vector_map(
    df: pd.DataFrame,
    feature_cols: Optional[Iterable[str]] = None,
    standardize: bool = True
) -> Dict[str, np.ndarray]:
    """
    Map one-letter AA -> feature vector (np.ndarray).
    If feature_cols is None, uses prop_names(df).
    """
    if feature_cols is None:
        feature_cols = prop_names(df)
    feature_cols = list(feature_cols)

    for c in feature_cols:
        if c not in df.columns:
            raise ValueError(f"aa_props_io: feature column '{c}' not found.")
        if not pd.api.types.is_numeric_dtype(df[c]):
            raise ValueError(f"aa_props_io: feature column '{c}' must be numeric.")

    sub = df[["aa"] + feature_cols].drop_duplicates("aa").set_index("aa").sort_index()
    X = sub.to_numpy(dtype=float)

    if standardize:
        X, _, _ = _zscore(X)

    aa_list = list(sub.index)
    return {aa_list[i]: X[i, :].astype(float, copy=True) for i in range(X.shape[0])}


def aa_prop_matrix(
    df: pd.DataFrame,
    feature_cols: Optional[Iterable[str]] = None,
    standardize: bool = True
) -> Tuple[List[str], np.ndarray]:
    """Return (aa_list, X) aligned to sorted AAs."""
    if feature_cols is None:
        feature_cols = prop_names(df)
    feature_cols = list(feature_cols)

    sub = df[["aa"] + feature_cols].drop_duplicates("aa").set_index("aa").sort_index()
    X = sub.to_numpy(dtype=float)

    if standardize:
        X, _, _ = _zscore(X)

    return list(sub.index), X

