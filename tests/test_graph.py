from src.sims.codon_graph import build_graph_n, hamming_distance

def test_hamming():
    assert hamming_distance("AUG","AUG")==0
    assert hamming_distance("AUG","AUA")==1
    assert hamming_distance("AUG","AAG")==1
    assert hamming_distance("AUG","CUG")==1
    assert hamming_distance("AUG","CUA")==2

def test_graph_counts():
    G = build_graph_n(3)
    assert len(G.nodes) == 64
    # Each node in a 3-letter, 4-state Hamming graph has degree 3*3=9
    assert set(dict(G.degree()).values()) == {9}

