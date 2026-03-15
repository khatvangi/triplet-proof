from src.io.codon_io import load_codon_table_csv
from src.io.aa_props_io import load_aa_props, aa_prop_vector_map, prop_names
from src.metrics.mi import mutual_info_by_position


def test_load_codon_table_preserves_base_columns():
    codon_df = load_codon_table_csv()

    for col in ("b1", "b2", "b3"):
        assert col in codon_df.columns

    aug = codon_df.loc[codon_df["codon"] == "AUG"].iloc[0]
    assert (aug["b1"], aug["b2"], aug["b3"]) == ("A", "U", "G")


def test_phase1_mutual_info_path_runs():
    codon_df = load_codon_table_csv()
    aa_df, _ = load_aa_props("aa_props.csv")
    feat_cols = prop_names(aa_df)
    aa_map = aa_prop_vector_map(aa_df, feature_cols=feat_cols)

    mi = mutual_info_by_position(codon_df, aa_map, feat_cols)

    assert set(mi) == set(feat_cols)
    assert all(len(vals) == 3 for vals in mi.values())
