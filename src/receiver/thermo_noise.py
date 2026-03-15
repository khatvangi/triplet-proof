# src/receiver/thermo_noise.py
from __future__ import annotations
import math
from typing import Dict, Iterable, List, Tuple
import numpy as np
import pandas as pd

# -------------------------
# Helpers
# -------------------------

def _l2(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.linalg.norm(a - b, ord=2))

def _vecmap_to_nd(aa_vec_map: Dict[str, Iterable[float]]) -> Dict[str, np.ndarray]:
    return {k: (v if isinstance(v, np.ndarray) else np.asarray(v, dtype=float)) for k, v in aa_vec_map.items()}

def _sense_codons(df: pd.DataFrame) -> List[str]:
    cols = set(df.columns)
    if "is_stop" not in cols:
        df = df.assign(is_stop=False)
    sense = df.loc[~df["is_stop"].astype(bool)].copy()
    sense["codon"] = sense["codon"].astype(str).str.upper()
    return sense["codon"].tolist()

def _hamming(a: str, b: str) -> int:
    return sum(1 for x, y in zip(a, b) if x != y)

def _neighbors_one_edit(codon: str, alphabet=("A","U","G","C")) -> List[Tuple[str, int]]:
    outs = []
    for i in range(len(codon)):
        for base in alphabet:
            if base != codon[i]:
                c2 = codon[:i] + base + codon[i+1:]
                outs.append((c2, i))  # (neighbor, position_of_mismatch)
    return outs

# Wobble weighting at position-3 (codon index 2)
# "strict": pos3 mismatches are much less penalized (Crick GU wobble)
# "expanded": even more forgiving (proxy for inosine/expanded wobble)
# "none": no special tolerance at pos3
def _wobble_scale(rule: str) -> float:
    rule = (rule or "strict").lower()
    if rule == "strict":
        return 0.4
    if rule == "expanded":
        return 0.25
    return 1.0  # "none" or unknown

# Positional asymmetry penalties (larger = rarer errors)
# If asymmetry=True, enforce pos2 > pos1 > pos3 via K.
def _position_penalties(asymmetry: bool, K: float, wobble_rule: str) -> Tuple[float,float,float]:
    # Base penalties (lower -> more likely; higher -> less likely)
    if not asymmetry:
        p1 = p2 = p3 = 1.0
    else:
        # Core asymmetry: pos2 hardest, pos1 medium, pos3 easiest
        # Example: K=3.0 → p2=3.0, p1=2.0, p3=1.0
        p2 = float(K)
        p1 = max(1.0, 0.5 * (1.0 + K))  # sits between p3 and p2
        p3 = 1.0

    # Apply wobble softness at pos3
    wob = _wobble_scale(wobble_rule)
    p3 *= wob
    # Return penalties for indices 0,1,2 → (pos1, pos2, pos3)
    return (p1, p2, p3)

# Convert penalties into probabilities via a Boltzmann softmax
def _normalize_boltz(weights: List[float], RT: float) -> List[float]:
    # weights are penalties; probability ~ exp(-penalty/RT)
    x = np.asarray([-w/RT for w in weights], dtype=float)
    x -= np.max(x)  # stabilize
    ex = np.exp(x)
    s = ex.sum()
    if s == 0.0 or not np.isfinite(s):
        # fallback to uniform if pathological
        return [1.0/len(weights)] * len(weights)
    return (ex / s).tolist()

# -------------------------
# Public API
# -------------------------

def expected_property_distortion_noise_multi(
    codon_df: pd.DataFrame,
    aa_vec_map: Dict[str, Iterable[float]],
    *,
    RT: float = 1.0,
    wobble_rules: str = "strict",
    asymmetry: bool = True,
    K: float = 3.0,
    alphabet: Tuple[str, ...] = ("A","U","G","C"),
) -> float:
    """
    Compute expected property distortion from near-cognate misdecoding,
    using a softmax over single-position mismatches.

    Parameters
    ----------
    codon_df : DataFrame with columns ['codon','aa','is_stop'] (STOP rows allowed; they are ignored)
    aa_vec_map : dict AA->vector (np.ndarray or list)
    RT : float
        Softmax temperature (lower RT -> sharper selection against penalties)
    wobble_rules : {'strict','expanded','none'}
        Controls tolerance at codon position-3 mismatches (index 2).
    asymmetry : bool
        If True, imposes pos2 > pos1 > pos3 penalty hierarchy.
    K : float
        Strength of pos2 penalty when asymmetry=True (typical 2.0–5.0).
    alphabet : tuple
        Nucleotide alphabet.

    Returns
    -------
    float : average expected distortion over sense codons.
    """
    # Normalize AA vectors to ndarray
    aa_vec = _vecmap_to_nd(aa_vec_map)

    # Build a map codon->AA vector for sense codons only
    sense_codons = _sense_codons(codon_df)
    # Look up AA letters for these codons
    # Expect an 'aa' column; filter NaN/STOP already done by _sense_codons via is_stop
    df_idx = codon_df.set_index(codon_df["codon"].str.upper())
    codon_to_vec = {}
    for c in sense_codons:
        aa = str(df_idx.loc[c, "aa"]).strip().upper()
        if aa in {"*", "STOP"} or aa not in aa_vec:
            # skip if unknown or STOP sneaks in
            continue
        codon_to_vec[c] = aa_vec[aa]

    if not codon_to_vec:
        return 0.0

    # Positional penalties
    p1, p2, p3 = _position_penalties(asymmetry, K, wobble_rules)
    pos_pen = (p1, p2, p3)

    total = 0.0
    n_used = 0

    # For each codon, gather its one-edit neighbors that are also sense codons
    for c, v_c in codon_to_vec.items():
        neighs = []
        pen = []
        for c2, pos in _neighbors_one_edit(c, alphabet=alphabet):
            v_n = codon_to_vec.get(c2)
            if v_n is None:
                continue
            neighs.append((c2, v_n))
            pen.append(pos_pen[pos])  # penalty by mismatch position

        if not neighs:
            continue

        # Convert penalties to probabilities
        probs = _normalize_boltz(pen, RT=RT)
        # Expected distortion at this codon
        e_local = 0.0
        for (_, v_n), p in zip(neighs, probs):
            e_local += p * _l2(v_c, v_n)

        total += e_local
        n_used += 1

    if n_used == 0:
        return 0.0
    return total / n_used


def expected_property_distortion_noise_multi_general(
    codon_df: pd.DataFrame,
    aa_vec_map: Dict[str, Iterable[float]],
    *,
    RT: float = 1.0,
    wobble_positions: set | None = None,
    asymmetry: bool = True,
    K: float = 3.0,
    wobble_strength: float = 0.4,
    alphabet: Tuple[str, ...] = ("A","U","G","C"),
) -> float:
    """
    Generalized noise distortion for arbitrary codon length n.

    Unlike the triplet-specific version, this accepts an explicit set of
    wobble position indices instead of a wobble_rules string, allowing it
    to work with doublets (n=2), quadruplets (n=4), or any codon length.

    Parameters
    ----------
    codon_df : DataFrame with columns ['codon','aa'] (and optionally 'is_stop')
    aa_vec_map : dict AA->vector
    RT : float
        Softmax temperature.
    wobble_positions : set of int or None
        0-indexed positions treated as wobble (low-penalty).
        e.g., {2} for triplet pos-3 wobble; set() for no wobble.
    asymmetry : bool
        If True, non-wobble positions get high penalty K.
    K : float
        Penalty for non-wobble positions when asymmetry=True.
    wobble_strength : float
        Multiplicative scaling on wobble-position penalty (lower = more error-prone).
    alphabet : tuple
        Nucleotide alphabet.

    Returns
    -------
    float : average expected distortion over sense codons.
    """
    if wobble_positions is None:
        wobble_positions = set()

    aa_vec = _vecmap_to_nd(aa_vec_map)
    sense_codons = _sense_codons(codon_df)
    df_idx = codon_df.set_index(codon_df["codon"].str.upper())

    codon_to_vec = {}
    for c in sense_codons:
        aa = str(df_idx.loc[c, "aa"]).strip().upper()
        if aa in {"*", "STOP"} or aa not in aa_vec:
            continue
        codon_to_vec[c] = aa_vec[aa]

    if not codon_to_vec:
        return 0.0

    # build per-position penalties from codon length
    n_pos = len(next(iter(codon_to_vec)))
    pos_pen = {}
    for i in range(n_pos):
        if i in wobble_positions:
            pos_pen[i] = 1.0 * wobble_strength
        elif asymmetry:
            pos_pen[i] = float(K)
        else:
            pos_pen[i] = 1.0

    total = 0.0
    n_used = 0

    for c, v_c in codon_to_vec.items():
        neighs = []
        pen = []
        for c2, pos in _neighbors_one_edit(c, alphabet=alphabet):
            v_n = codon_to_vec.get(c2)
            if v_n is None:
                continue
            neighs.append((c2, v_n))
            pen.append(pos_pen[pos])

        if not neighs:
            continue

        probs = _normalize_boltz(pen, RT=RT)
        e_local = 0.0
        for (_, v_n), p in zip(neighs, probs):
            e_local += p * _l2(v_c, v_n)

        total += e_local
        n_used += 1

    if n_used == 0:
        return 0.0
    return total / n_used

