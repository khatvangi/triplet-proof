"""fast reruns (1)-(5) for one amino-acid vector matrix.
usage: python stats.py <tag> <vectors.csv> <colprefix> <rescale 0/1>
writes runs/<tag>.json. every draw uses fixed seeds so old/new matrices see identical nulls."""
import itertools, json, math, os, sys, time
os.environ.setdefault("NUMBA_CACHE_DIR", os.path.join(os.path.dirname(os.path.abspath(__file__)), "numba_cache"))
sys.dont_write_bytecode = True
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
RR = os.path.dirname(HERE)
ROOT = os.path.dirname(RR)
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(RR, "arch_cert"))
sys.path.insert(0, os.path.join(RR, "usage_cert"))
sys.path.insert(0, os.path.join(RR, "error_budget", "scripts"))
from run_publication_controls import RAA10, build_doublet_table  # noqa: E402
import usage_cert as uc  # noqa: E402
import decode_split as ds  # noqa: E402

SEED = 20260924
BASES = "ACGU"
REPS = dict(zip([f"X{i}" for i in range(1, 11)], "ALWSNDRHPC"))
GATE_B, GATE_A = 23.9223597, 18.9372988


def loadVectors(path, prefix, rescale):
    df = pd.read_csv(path)
    cols = [c for c in df.columns if c.startswith(prefix)]
    X = df[cols].to_numpy(float)
    if rescale:
        X = (X - X.mean(0)) / X.std(0, ddof=0)
    return {a: X[i] for i, a in enumerate(df["aa"].str.upper())}, cols


def tables():
    tri = pd.read_csv(os.path.join(RR, "codon_table_export.csv"))
    tri["codon"] = tri["codon"].str.upper()
    tri["is_stop"] = tri["is_stop"].astype(bool)
    sense = tri[~tri.is_stop]
    B = dict(zip(sense.codon, sense.aa.str.upper()))
    C = {c: RAA10[a] for c, a in B.items()}
    dbl = build_doublet_table(tri, RAA10)
    assert not dbl.is_stop.any()
    A = dict(zip(dbl.codon, dbl.aa))
    return {"A": A, "B": B, "C": C}


def condition(code, vec):
    """code: {codon: label}. returns labels (int), classes, dist, directed edges, N."""
    codons = sorted(code)
    classes = sorted(set(code.values()))
    idx = {c: i for i, c in enumerate(codons)}
    lab = np.array([classes.index(code[c]) for c in codons])
    E = []
    for c in codons:
        for p in range(len(c)):
            for b in BASES:
                if b != c[p]:
                    m = c[:p] + b + c[p + 1:]
                    if m in idx:
                        E.append((idx[c], idx[m]))
    E = np.array(E)
    V = np.array([vec[x] for x in classes])
    dist = np.sqrt(((V[:, None] - V[None]) ** 2).sum(-1))
    return dict(lab=lab, classes=classes, dist=dist, E=E, N=len(codons))


def DE(L, cond, dist):
    """D = sum over directed sense edges d / N_sense ; E = sum over undirected edges d^2 = half the directed sum."""
    a, b = L[..., cond["E"][:, 0]], L[..., cond["E"][:, 1]]
    d = dist[a, b]
    return d.sum(-1) / cond["N"], (d ** 2).sum(-1) / 2.0


def zstat(v, null):
    m, s = null.mean(), null.std()
    tol = 1e-9 * abs(v)
    return dict(sgc=float(v), null_mean=float(m), null_sd=float(s), z=float((v - m) / s) if s > 0 else None,
                p_strict_better=float((null < v - tol).mean()), p_le=float((null <= v + tol).mean()),
                median=float(np.median(null)), n=int(len(null)))


def codonShuffle(cond, n=1_000_000, chunk=100_000):
    rng = np.random.default_rng(SEED)
    k = len(cond["classes"])
    eq = 1.0 - np.eye(k)
    out = {m: [] for m in ("D", "E", "D_eq", "E_eq")}
    for _ in range(n // chunk):
        L = rng.permuted(np.tile(cond["lab"], (chunk, 1)), axis=1)
        d, e = DE(L, cond, cond["dist"])
        out["D"].append(d); out["E"].append(e)
        d, e = DE(L, cond, eq)
        out["D_eq"].append(d); out["E_eq"].append(e)
    out = {m: np.concatenate(v) for m, v in out.items()}
    d0, e0 = DE(cond["lab"], cond, cond["dist"])
    d1, e1 = DE(cond["lab"], cond, eq)
    return {"D": zstat(d0, out["D"]), "E": zstat(e0, out["E"]),
            "equidistant": {"D": zstat(d1, out["D_eq"]), "E": zstat(e1, out["E_eq"])}}


def blockNull(cond, exhaustive, n=10_000_000, chunk=200_000):
    k = len(cond["classes"])
    M = np.zeros((k, k))
    for u, v in cond["E"]:
        M[cond["lab"][u], cond["lab"][v]] += 1
    I, J = np.nonzero(M)
    w = M[I, J]
    dist = cond["dist"]

    def score(P):
        d = dist[P[:, I], P[:, J]]
        return (d * w).sum(1) / cond["N"], (d ** 2 * w).sum(1) / 2.0

    ident = np.arange(k)[None]
    d0, e0 = score(ident)
    Ds, Es = [], []
    if exhaustive:
        it = itertools.permutations(range(k))
        while True:
            P = np.array(list(itertools.islice(it, 500_000)), dtype=np.int64)
            if len(P) == 0:
                break
            d, e = score(P); Ds.append(d); Es.append(e)
    else:
        rng = np.random.default_rng(SEED)
        for _ in range(n // chunk):
            P = rng.permuted(np.tile(np.arange(k), (chunk, 1)), axis=1)
            d, e = score(P); Ds.append(d); Es.append(e)
    Dn, En = np.concatenate(Ds), np.concatenate(Es)
    return {"D": zstat(float(d0[0]), Dn), "E": zstat(float(e0[0]), En),
            "n_blocks": k, "exhaustive": exhaustive, "n_draws": int(len(Dn))}


def part3and1(vec, do_part1):
    aas = ds.AAS
    X = np.array([vec[a] for a in aas])
    Dm = pd.DataFrame(np.sqrt(((X[:, None] - X[None]) ** 2).sum(-1)), index=aas, columns=aas)
    iu = np.triu_indices(20, 1)
    r = {"mean_d_190": float(Dm.to_numpy()[iu].mean()), "median_d_190": float(np.median(Dm.to_numpy()[iu]))}
    r["part3"] = ds.part3(Dm, np.random.default_rng(ds.SEED))  # decode_split's own seed -> old values reproducible
    if do_part1:
        known = {}
        p2 = os.path.join(RR, "error_budget", "part2_misacylation_table.json")
        for x in json.load(open(p2)):
            if isinstance(x.get("rate_per_aminoacylation"), (int, float)):
                known.setdefault((x["cognate_aa"], x["noncognate_aa"]), f"yes: {x.get('enzyme')} ({x.get('doi')})")
            else:
                known.setdefault((x["cognate_aa"], x["noncognate_aa"]), "listed in PART 2 but rate UNAVAILABLE")
        r["part1"] = ds.part1(Dm, known)
    return r


def usage(vec, tabs):
    chem = uc.distMatrix(uc.CLASS_ORDER, {c: vec[a] for c, a in REPS.items()})
    comp, meta = uc.loadProteome()
    assert not meta["fallback_used"]
    c61, nbr61 = uc.buildSenseGraph(3, uc.STOPS)
    c16, nbr16 = uc.buildSenseGraph(2)
    tri_lab = uc.labelsFromMap(c61, tabs["C"], uc.CLASS_ORDER)
    dbl_lab = uc.labelsFromMap(c16, tabs["A"], uc.CLASS_ORDER)
    pres = sorted(set(dbl_lab.tolist()))
    fP = uc.classUsage(comp, set(uc.CLASS_ORDER))
    fP9 = uc.classUsage(comp, {uc.CLASS_ORDER[i] for i in pres})
    fU = uc.uniformUsage(10, list(range(10)))
    return {"proteome": meta,
            "doublet_projection": {"Dw_U": float(uc.fullDw(dbl_lab, nbr16, chem, uc.uniformUsage(10, pres), 10)),
                                   "Dw_P": float(uc.fullDw(dbl_lab, nbr16, chem, fP9, 10))},
            "sgc_10class_projection": {"Dw_U": float(uc.fullDw(tri_lab, nbr61, chem, fU, 10)),
                                       "Dw_P": float(uc.fullDw(tri_lab, nbr61, chem, fP, 10))}}


def main():
    tag, path, prefix, rescale = sys.argv[1], sys.argv[2], sys.argv[3], bool(int(sys.argv[4]))
    gate_only = len(sys.argv) > 5 and sys.argv[5] == "gate"
    t0 = time.time()
    vec, cols = loadVectors(path, prefix, rescale)
    tabs = tables()
    conds = {"A": condition(tabs["A"], {c: vec[a] for c, a in REPS.items()}),
             "B": condition(tabs["B"], vec),
             "C": condition(tabs["C"], {c: vec[a] for c, a in REPS.items()})}
    res = {"tag": tag, "vectors": path, "columns": cols, "rescaled": rescale, "seed": SEED,
           "n_sense": {k: c["N"] for k, c in conds.items()}, "n_classes": {k: len(c["classes"]) for k, c in conds.items()},
           "sgc_D": {k: float(DE(c["lab"], c, c["dist"])[0]) for k, c in conds.items()},
           "sgc_E": {k: float(DE(c["lab"], c, c["dist"])[1]) for k, c in conds.items()}}
    if gate_only:
        res["gate"] = {"B_D": res["sgc_D"]["B"], "B_absdiff": abs(res["sgc_D"]["B"] - GATE_B),
                       "A_D": res["sgc_D"]["A"], "A_absdiff": abs(res["sgc_D"]["A"] - GATE_A)}
        print(json.dumps(res["gate"], indent=1))
        return
    out = os.path.join(HERE, "runs", f"{tag}.json")
    dump = lambda: json.dump(res, open(out, "w"), indent=1)  # noqa: E731
    res["s1_codon_shuffle"] = {k: codonShuffle(c) for k, c in conds.items()}
    print(tag, "s1", round(time.time() - t0), flush=True); dump()
    res["s2_block_null"] = {"A": blockNull(conds["A"], True), "C": blockNull(conds["C"], True),
                            "B": blockNull(conds["B"], False)}
    print(tag, "s2", round(time.time() - t0), flush=True); dump()
    res["s3"] = part3and1(vec, do_part1=False)
    print(tag, "s3", round(time.time() - t0), flush=True); dump()
    res["s5_usage"] = usage(vec, tabs)
    print(tag, "s5", round(time.time() - t0), flush=True); dump()
    # part1 writes per-dataset/per-codon csvs to ds.EB: redirect into this directory
    ds.EB = os.path.join(HERE, "runs", f"decode_{tag}")
    os.makedirs(ds.EB, exist_ok=True)
    res["s4_landerer"] = part3and1(vec, do_part1=True)["part1"]
    res["seconds"] = round(time.time() - t0, 1)
    print(tag, "s4 done", res["seconds"], flush=True); dump()


if __name__ == "__main__":
    main()
