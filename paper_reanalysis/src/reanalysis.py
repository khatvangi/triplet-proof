"""nulls, objectives and property representations for the reanalysis.

nothing here is new science.  E and D are transcribed from
run_publication_controls.py, and the block null is the Haig & Hurst /
Freeland & Hurst scheme the manuscript's Methods actually describe.  the point
of the module is to run BOTH nulls and BOTH property regimes through one code
path, so the comparison cannot be an artefact of two different implementations.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))

from src.io.aa_props_io import aa_prop_vector_map, load_aa_props   # noqa: E402
from src.io.codon_io import load_codon_table_csv                   # noqa: E402
from src.sims.codon_graph import build_graph_n                     # noqa: E402

PUBLISHED_E = 2984.675406924604
PUBLISHED_D = 23.922359743253224


class Space:
    """sense codons, their Hamming-1 edges, and the SGC labelling."""

    def __init__(self, codon_df=None, n_bases=3):
        df = load_codon_table_csv() if codon_df is None else codon_df
        sense = df[~df["is_stop"].astype(bool)]
        c2a = {c.upper(): a.upper() for c, a in zip(sense["codon"], sense["aa"])}
        self.codons = sorted(c2a)
        pos = {c: i for i, c in enumerate(self.codons)}
        G = build_graph_n(n_bases)
        e = [(pos[u], pos[v]) for u, v in G.edges() if u in pos and v in pos]
        e = np.array(sorted(e), dtype=np.int32)
        self.EI, self.EJ = e[:, 0], e[:, 1]
        self.aas = sorted(set(c2a.values()))
        self.labels = np.array([self.aas.index(c2a[c]) for c in self.codons],
                               dtype=np.int16)
        self.n_sense = len(self.codons)
        self.n_edges = len(self.EI)


def propertySets(sp: Space) -> dict:
    """the five frozen representations, all on the same amino-acid ordering."""
    aa_df, _ = load_aa_props(str(REPO / "data" / "processed" / "aa_props.parquet"))
    vm = aa_prop_vector_map(aa_df)                      # standardize=True -> WHITENED
    pc8_whitened = np.array([vm[a] for a in sp.aas])

    raw = pd.read_parquet(REPO / "data" / "processed" / "aa_props.parquet")
    raw = raw.set_index("aa").loc[list(sp.aas)]
    pcs = sorted([c for c in raw.columns if c.startswith("pca")],
                 key=lambda s: int(s[3:]))
    pc8_raw = raw[pcs].to_numpy(float)                  # pure rotation, unwhitened
    desc20_z = raw[[c for c in raw.columns if c.endswith("_z")]].to_numpy(float)

    prs = _loadPrs()
    prs1d = np.array([[prs[a]] for a in sp.aas], dtype=float)

    # chemistry removed: ||e_a - e_b|| = 1 for a != b, 0 for a == b
    blind = np.eye(len(sp.aas)) / np.sqrt(2.0)

    return {"pc8_whitened": pc8_whitened, "pc8_raw": pc8_raw,
            "desc20_z": desc20_z, "prs1d": prs1d, "blind": blind}


def _loadPrs() -> dict:
    p = REPO / "code_evolution" / "data" / "nwk2007" / "prs_woese.csv"
    t = pd.read_csv(p, comment="#")
    return dict(zip(t.iloc[:, 0].astype(str), t.iloc[:, 1].astype(float)))


# --- objectives -------------------------------------------------------------

def score(P: np.ndarray, L: np.ndarray, sp: Space):
    """(E, D) for a batch of label vectors L of shape (B, n_sense).

    E = sum over unordered sense-sense edges of ||f(u)-f(v)||^2          (= H2)
    D = 2 * sum of the unsquared distances / n_sense                     (= 2*H1/n)
    the second is exactly what run_publication_controls.py computes: it sums over
    ORDERED (sender -> mutant) pairs, which double-counts each edge, and divides
    once by n_sense.  no 1/|N(c)| factor is applied there despite the Methods.
    """
    L = np.atleast_2d(L)
    d = P[L[:, sp.EI]] - P[L[:, sp.EJ]]
    sq = np.einsum("bij,bij->bi", d, d)
    E = sq.sum(axis=1)
    D = 2.0 * np.sqrt(sq).sum(axis=1) / sp.n_sense
    return E, D


# --- the two nulls ----------------------------------------------------------

def codonNull(sp: Space, n: int, rng) -> np.ndarray:
    """permute labels among individual sense codons.  DESTROYS every block --
    this is what the code ran, and why the statistic counts synonymy."""
    perm = np.argsort(rng.random((n, sp.n_sense)), axis=1)
    return sp.labels[perm]


def blockNull(sp: Space, n: int, rng) -> np.ndarray:
    """permute the amino-acid labels over the fixed block partition:
    L = pi[L_sgc].  every draw relabels the SAME partition, so the set of
    synonymous edges is identical in every code and synonymy cannot contribute
    to the variance.  this is the Methods' description and the Haig & Hurst /
    Freeland & Hurst scheme."""
    k = len(sp.aas)
    perms = np.argsort(rng.random((n, k)), axis=1).astype(np.int16)
    return perms[:, sp.labels]


NULLS = {"codon": codonNull, "block": blockNull}


def zstats(obs: float, null: np.ndarray) -> dict:
    mu, sd = float(null.mean()), float(null.std(ddof=1))
    # a statistic can be exactly invariant under a null (the block null cannot move
    # a chemistry-blind score at all).  float arithmetic leaves ~1e-15 of noise, so
    # test the spread against the scale of the statistic, not against exact zero --
    # otherwise sd=7e-15 divides into a meaningless z of ~1.
    invariant = sd <= 1e-12 * max(1.0, abs(mu))
    n_better = int((null < obs).sum())
    n_tie = int(np.isclose(null, obs, rtol=1e-12, atol=0.0).sum())
    # ties are NOT evidence for the SGC: count them against it.  with every draw
    # tied this gives p = 1, which is the honest reading of an invariant statistic.
    return {"sgc": obs, "null_mean": mu, "null_sd": sd,
            "z": float("nan") if invariant else (obs - mu) / sd,
            "n_better": n_better, "n_tie": n_tie, "n": int(null.size),
            "p_laplace": (n_better + n_tie + 1) / (null.size + 1),
            "null_min": float(null.min()), "zero_variance": invariant}


def positiveControl(sp: Space, props: dict) -> dict:
    """gate: reproduce the published SGC values exactly, or stop."""
    E, D = score(props["pc8_whitened"], sp.labels[None, :], sp)
    dE, dD = abs(E[0] - PUBLISHED_E), abs(D[0] - PUBLISHED_D)
    return {"E": float(E[0]), "D": float(D[0]), "abs_err_E": float(dE),
            "abs_err_D": float(dD), "passed": bool(dE < 1e-8 and dD < 1e-10)}
