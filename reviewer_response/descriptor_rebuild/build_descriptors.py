"""rebuild the 20-aa descriptor matrix from sourced data (RULES.md), then pca.
outputs: descriptors_raw.csv, table_S2_new.csv, NEW_VECTORS.csv (primary: k97 pcs, rescaled),
NEW_VECTORS_k8_rescaled.csv, NEW_VECTORS_k97_unscaled.csv, NEW_VECTORS_k8_unscaled.csv, build_meta.json"""
import hashlib, json, os
import numpy as np
import pandas as pd
import rdkit
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors as rd

HERE = os.path.dirname(os.path.abspath(__file__))
AAS = sorted("ARNDCQEGHILKMFPSTWYV")
AAIDX_ORDER = "ARNDCQEGHILKMFPSTWYV"  # column order of the aaindex1 'I' lines

# free l-amino acids, neutral (non-zwitterion) form; pubchem-style canonical smiles
SMILES = {
    "G": "NCC(=O)O", "A": "C[C@@H](C(=O)O)N", "V": "CC(C)[C@@H](C(=O)O)N",
    "L": "CC(C)C[C@@H](C(=O)O)N", "I": "CC[C@H](C)[C@@H](C(=O)O)N", "M": "CSCC[C@@H](C(=O)O)N",
    "F": "c1ccc(cc1)C[C@@H](C(=O)O)N", "W": "c1ccc2c(c1)c(c[nH]2)C[C@@H](C(=O)O)N",
    "Y": "c1cc(ccc1C[C@@H](C(=O)O)N)O", "S": "C([C@@H](C(=O)O)N)O", "T": "C[C@H]([C@@H](C(=O)O)N)O",
    "C": "C([C@@H](C(=O)O)N)S", "N": "C([C@@H](C(=O)O)N)C(=O)N", "Q": "C(CC(=O)N)[C@@H](C(=O)O)N",
    "D": "C([C@@H](C(=O)O)N)C(=O)O", "E": "C(CC(=O)O)[C@@H](C(=O)O)N", "K": "C(CCN)C[C@@H](C(=O)O)N",
    "R": "C(C[C@@H](C(=O)O)N)CN=C(N)N", "H": "c1c([nH]cn1)C[C@@H](C(=O)O)N", "P": "C1C[C@H](NC1)C(=O)O",
}

# rdkit rules: counts on the free amino acid; the backbone contributes the same constant to
# every residue except where the side chain alters it (pro), and a constant offset is removed
# by the z-score, so these equal side-chain counts up to that offset.
RULES_B = {
    "acceptors": ("rdMolDescriptors.CalcNumHBA(free aa)", lambda m: rd.CalcNumHBA(m)),
    "aromaticity": ("1 if any aromatic atom in the molecule else 0", lambda m: int(any(a.GetIsAromatic() for a in m.GetAtoms()))),
    "ring_count": ("rdMolDescriptors.CalcNumRings(free aa)", lambda m: rd.CalcNumRings(m)),
    "rotatable_bonds": ("rdMolDescriptors.CalcNumRotatableBonds(free aa) (default strictness)", lambda m: rd.CalcNumRotatableBonds(m)),
}

# fixed concept order (RULES.md). route: ('a', accession, reason) | ('b', rule key, reason) | ('drop', None, reason)
CONCEPTS = [
    ("hydropathy_KD", "a", "KYTJ820101", "named source present in AAindex"),
    ("volume_Zamyatnin", "a", "BIGC670101", "Zamyatnin 1972 not in AAindex; Bigelow 1967 residue volume is the AAindex volume entry derived, like Zamyatnin's, from solution partial specific volumes"),
    ("polarity_Grantham", "a", "GRAR740102", "named source present in AAindex"),
    ("isoelectric_point", "a", "ZIMJ680104", "AAindex isoelectric point entry"),
    ("sidechain_pKa", "drop", None, "no AAindex entry for side-chain pKa (FASG760104/5, JOND750102, FAUJ880113 are alpha-amino/carboxyl pKs); not computable by RDKit; no primary table saved -> DROPPED"),
    ("net_charge_pH7", "a", "KLEP840101", "AAindex net charge entry"),
    ("hbond_donors", "a", "FAUJ880109", "AAindex entry titled 'Number of hydrogen bond donors'"),
    ("hbond_acceptors", "b", "acceptors", "no AAindex acceptor-count entry; structure rule"),
    ("aromaticity", "b", "aromaticity", "no AAindex aromaticity indicator; structure rule"),
    ("ring_count", "b", "ring_count", "no AAindex ring-count entry; structure rule"),
    ("rotatable_bonds", "b", "rotatable_bonds", "no AAindex rotatable-bond entry; structure rule"),
    ("molar_refractivity", "a", "MCMT640101", "AAindex refractivity entry (McMeekin 1964)"),
    ("bulkiness_Zimmerman", "a", "ZIMJ680102", "named source present in AAindex"),
    ("chou_fasman_alpha", "a", "CHOP780201", "named source present in AAindex"),
    ("chou_fasman_beta", "a", "CHOP780202", "named source present in AAindex"),
    ("chou_fasman_turn", "a", "CHOP780203", "named source present in AAindex"),
    ("solvent_accessibility", "a", "ROSG850102", "no named source; relative accessibility concept -> Rose 1985 mean fractional area loss (= mean 1 - relative accessibility)"),
    ("flexibility", "a", "VINM940101", "no named source; 'normalized flexibility index' -> Vihinen 1994 normalized flexibility (B-values), average"),
    ("helical_moment", "a", "EISD860102", "per-residue helical moment is not a tabulated property; closest AAindex 'moment' entry is the atom-based hydrophobic moment (Eisenberg-McLachlan 1986), which is NOT helix-specific"),
    ("strand_moment", "drop", None, "no AAindex strand-moment entry (EISD860102 already used and not strand-specific); no structure rule -> DROPPED"),
    ("hydrophobic_surface", "a", "ROSG850101", "no named source; 'hydrophobic surface tendency' -> Rose 1985 mean area buried on transfer"),
    ("sidechain_mass", "a", "FASG760101", "molecular weight (Fasman 1976) = side-chain mass + a constant; identical after z-scoring"),
]


def readAAindex(path):
    ent, cur, key = {}, None, None
    lines = open(path).read().splitlines()
    i = 0
    while i < len(lines):
        ln = lines[i]
        if ln.startswith("H "):
            cur = {"acc": ln[2:].strip()}
        elif ln.startswith("D "):
            cur["title"] = ln[2:].strip()
        elif ln.startswith("I "):
            vals = (lines[i + 1] + " " + lines[i + 2]).split()
            cur["values"] = {a: (float(v) if v != "NA" else np.nan) for a, v in zip(AAIDX_ORDER, vals)}
            ent[cur["acc"]] = cur
            i += 2
        i += 1
    return ent


def main():
    src = os.path.join(HERE, "sources", "aaindex1")
    idx = readAAindex(src)
    mols = {a: Chem.MolFromSmiles(s) for a, s in SMILES.items()}
    cols, rows = {}, []
    for name, route, key, reason in CONCEPTS:
        if route == "a":
            e = idx[key]
            v = pd.Series(e["values"]).loc[AAS]
            assert not v.isna().any(), (key, "has NA")
            cols[name] = v.to_numpy(float)
            rows.append(dict(concept=name, route="AAindex1", accession_or_rule=key, title=e["title"],
                             units="as in AAindex entry", status="candidate", reason=reason))
        elif route == "b":
            desc, fn = RULES_B[key]
            cols[name] = np.array([fn(mols[a]) for a in AAS], float)
            rows.append(dict(concept=name, route="RDKit", accession_or_rule=desc, title=f"RDKit {rdkit.__version__} on free-AA SMILES",
                             units="count" if key != "aromaticity" else "0/1", status="candidate", reason=reason))
        else:
            rows.append(dict(concept=name, route="none", accession_or_rule="", title="", units="",
                             status="dropped", reason=reason))
    raw = pd.DataFrame(cols, index=AAS)
    raw.index.name = "aa"
    raw.to_csv(os.path.join(HERE, "descriptors_raw.csv"))
    # z-score (ddof=0) then ordered correlation pruning |r| >= 0.95 vs an earlier kept descriptor
    Z = (raw - raw.mean()) / raw.std(ddof=0)
    C = Z.corr()
    kept = []
    tab = {r["concept"]: r for r in rows}
    for c in Z.columns:
        hit = [k for k in kept if abs(C.loc[c, k]) >= 0.95]
        if hit:
            tab[c]["status"] = "dropped"
            tab[c]["reason"] += f" | DROPPED by correlation: |r|={abs(C.loc[c, hit[0]]):.3f} with earlier {hit[0]}"
        else:
            kept.append(c)
            tab[c]["status"] = "kept"
    pd.DataFrame(rows).to_csv(os.path.join(HERE, "table_S2_new.csv"), index=False)
    Xz = Z[kept].to_numpy()
    U, S, Vt = np.linalg.svd(Xz, full_matrices=False)
    var = S ** 2 / (Xz.shape[0] - 1)
    cum = np.cumsum(var / var.sum())
    k97 = int(np.searchsorted(cum, 0.97) + 1)
    P = Xz @ Vt.T
    out = {}
    for tag, k in (("k97", k97), ("k8", 8)):
        pcs = P[:, :k]
        resc = (pcs - pcs.mean(0)) / pcs.std(0, ddof=0)
        for sc, M in (("rescaled", resc), ("unscaled", pcs)):
            df = pd.DataFrame(M, index=AAS, columns=[f"pc{i+1}" for i in range(k)])
            df.index.name = "aa"
            fn = "NEW_VECTORS.csv" if (tag, sc) == ("k97", "rescaled") else f"NEW_VECTORS_{tag}_{sc}.csv"
            df.to_csv(os.path.join(HERE, fn))
            out[f"{tag}_{sc}"] = fn
    meta = dict(rdkit_version=rdkit.__version__, smiles=SMILES,
                aaindex1_sha256=hashlib.sha256(open(src, "rb").read()).hexdigest(),
                n_candidates=len(cols), kept=kept, n_kept=len(kept), k97=k97,
                cum_var=cum.tolist(), cum_var_at_k97=float(cum[k97 - 1]), cum_var_at_8=float(cum[7]) if len(cum) >= 8 else None,
                files=out, raw_corr_abs_max_offdiag=float((C.abs() - np.eye(len(C))).to_numpy().max()))
    json.dump(meta, open(os.path.join(HERE, "build_meta.json"), "w"), indent=1)
    print(json.dumps({k: v for k, v in meta.items() if k not in ("smiles",)}, indent=1))
    print(pd.DataFrame(rows)[["concept", "accession_or_rule", "status"]].to_string())


if __name__ == "__main__":
    main()
