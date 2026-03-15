import numpy as np
import pandas as pd
from sklearn.feature_selection import mutual_info_regression

def mutual_info_by_position(codon_df: pd.DataFrame, aa_prop_vecs: dict, prop_names: list):
    """
    MI between base-position one-hot and each AA property.
    Drops stops. Returns {prop: [MI_pos1, MI_pos2, MI_pos3]} (means over one-hot features).
    """
    df = codon_df.copy()
    df = df[(~df["is_stop"]) & (df["aa"].notna())]

    X_rows = []
    keep_idx = []
    for idx, row in df.iterrows():
        aa = row["aa"]
        if aa not in aa_prop_vecs:
            continue
        b1 = [int(row["b1"]==b) for b in ["U","C","A","G"]]
        b2 = [int(row["b2"]==b) for b in ["U","C","A","G"]]
        b3 = [int(row["b3"]==b) for b in ["U","C","A","G"]]
        X_rows.append(b1 + b2 + b3)
        keep_idx.append(idx)

    X = np.array(X_rows)
    df_keep = df.loc[keep_idx].reset_index(drop=True)

    results = {}
    for p_i, pname in enumerate(prop_names):
        Y = np.array([aa_prop_vecs[df_keep.loc[i,"aa"]][p_i] for i in range(len(df_keep))])
        MIs = []
        for pos in range(3):
            cols = list(range(pos*4, (pos+1)*4))
            mi = mutual_info_regression(X[:, cols], Y, discrete_features=True, random_state=0)
            MIs.append(float(np.mean(mi)))
        results[pname] = MIs
    return results

