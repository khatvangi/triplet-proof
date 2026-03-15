import networkx as nx
from itertools import product

BASES = ["U","C","A","G"]

def hamming_distance(a: str, b: str) -> int:
    return sum(x!=y for x,y in zip(a,b))

def build_graph_n(n: int) -> nx.Graph:
    nodes = [''.join(p) for p in product(BASES, repeat=n)]
    G = nx.Graph()
    G.add_nodes_from(nodes)
    for i, u in enumerate(nodes):
        for v in nodes[i+1:]:
            if hamming_distance(u,v)==1:
                G.add_edge(u,v)
    return G

def label_nodes_with_aa(G: nx.Graph, codon_to_aa: dict, exclude_stops=True):
    aa_labels = {}
    for node in G.nodes:
        aa = codon_to_aa.get(node, None)
        if exclude_stops and (aa is None or aa == "*" or aa == ""):
            aa_labels[node] = None
        else:
            aa_labels[node] = aa
    nx.set_node_attributes(G, aa_labels, "aa")
    return G

