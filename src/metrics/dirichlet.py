import numpy as np

def dirichlet_energy(G, prop_map: dict) -> float:
    """
    Sum over edges of ||p(u) - p(v)||^2 for nodes with property vectors.
    Skips edges touching stops/None.
    """
    E = 0.0
    for u, v in G.edges():
        pu = prop_map.get(u, None)
        pv = prop_map.get(v, None)
        if pu is None or pv is None:
            continue
        duv = np.array(pu) - np.array(pv)
        E += float(np.dot(duv, duv))
    return E

