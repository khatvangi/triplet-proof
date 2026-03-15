# build_aa_props.py  (place at REPO ROOT)
# builds the processed parquet from the 22-property library (aa_props_lib.py)
# for the full PCA feature set used by phase-2 runners.
import os, sys, json
import pandas as pd
sys.path.insert(0, os.path.abspath(os.path.dirname(__file__)))

from src.io.aa_props_lib import AA, aa_property_columns
from src.io.aa_props_io import load_aa_props

def main():
    # build a CSV-like DataFrame from the 22-property library
    rows = [{"aa": aa, **props} for aa, props in AA.items()]
    lib_df = pd.DataFrame(rows)
    # write a temporary CSV so load_aa_props can process it
    tmp_csv = "data/processed/_aa_props_full.csv"
    os.makedirs("data/processed", exist_ok=True)
    lib_df.to_csv(tmp_csv, index=False)

    # load with auto-select + PCA (22 props → ~8 PCs at 97% variance)
    aa_df, meta = load_aa_props(tmp_csv, auto_select=True, corr_thresh=0.95, pca_var=0.97)
    aa_df.to_parquet("data/processed/aa_props.parquet", index=False)
    with open("data/processed/aa_props_meta.json", "w") as f:
        json.dump(meta, f, indent=2)

    cols = [c for c in aa_df.columns if c.startswith("pca")]
    print(f"Wrote data/processed/aa_props.parquet")
    print(f"Features after auto-select: {meta.get('raw_used')}")
    print(f"PCA columns: {cols}")
    print(f"Explained variance: {meta.get('var_explained')}")

if __name__ == "__main__":
    main()

