#!/usr/bin/env python
"""
dataset-level bootstrap for the landerer per-opportunity ratios.
psms cluster strongly by pride dataset, so resample whole datasets
(with replacement, 2000 draws, seed 0) and recompute
eps1/eps2 and eps3/eps2 (all-pairs opportunity definition, all events,
codon/origin-mismatch rows dropped). reports 2.5/50/97.5 percentiles and
the leave-largest-dataset-out point estimate.
"""
import json
import numpy as np
from landerer_position_rates import load, assignPos, neighbours, cls, CODE

rng = np.random.default_rng(0)
out = {}
for org in ("ecoli", "yeast"):
    ev, ex = load(org)
    ex = ex[ex.codon.map(CODE) != "*"]
    ev = ev[ev.origin == ev.codon.map(CODE)].copy()
    ev["pos"] = [assignPos(c, d) for c, d in zip(ev.codon, ev.destination)]
    # m_k(c) for all 61 sense codons
    mk = {c: [len([d for d in neighbours(c, k) if d not in ("*", cls(a))]) for k in range(3)]
          for c, a in CODE.items() if a != "*"}
    dss = sorted(ev.ds.unique())
    N = np.array([[int(((ev.ds == d) & (ev.pos == k)).sum()) for k in "123"] for d in dss], float)
    O = np.array([[sum(bc * mk[c][k] for c, bc in zip(g.codon, g.base_count))
                   for k in range(3)] for d, g in ((d, ex[ex.ds == d]) for d in dss)], float)
    def ratio(idx):
        e = N[idx].sum(0) / O[idx].sum(0)
        return e[0] / e[1], e[2] / e[1]
    boots = np.array([ratio(rng.integers(0, len(dss), len(dss))) for _ in range(2000)])
    big = int(np.argmax(N.sum(1)))
    keep = [i for i in range(len(dss)) if i != big]
    pt = ratio(np.arange(len(dss)))
    out[org] = dict(
        n_datasets=len(dss),
        point=dict(eps1_eps2=pt[0], eps3_eps2=pt[1]),
        boot_2p5_50_97p5=dict(eps1_eps2=np.percentile(boots[:, 0], [2.5, 50, 97.5]).tolist(),
                              eps3_eps2=np.percentile(boots[:, 1], [2.5, 50, 97.5]).tolist()),
        largest_dataset=dss[big], largest_dataset_share_of_single_mismatch_psm=float(N[big].sum() / N.sum()),
        leave_largest_out=dict(zip(("eps1_eps2", "eps3_eps2"), ratio(np.array(keep)))),
    )
    print(org, json.dumps(out[org], indent=1))
json.dump(out, open("landerer2024_bootstrap.json", "w"), indent=2)
