"""(3) statistics split into independent stages so they run in parallel on a loaded machine.
usage: python part3_stage.py <stage>   stage in s1A s1B s1C s2 s2prim s3 s4
uses VECTORS_r094.csv (written by part3_stats.py, whose gate reproduced the primary matrix)
and stats.py's own functions (read-only import). writes part3_<stage>.json."""
import json
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
os.environ["NUMBA_CACHE_DIR"] = os.path.join(HERE, "numba_cache")
sys.dont_write_bytecode = True
RR = os.path.dirname(HERE)
DR = os.path.join(RR, "descriptor_rebuild")
sys.path.insert(0, DR)
import stats as ST  # noqa: E402


def main():
    stage = sys.argv[1]
    t0 = time.time()
    path = os.path.join(DR, "NEW_VECTORS.csv") if stage == "s2prim" else os.path.join(HERE, "VECTORS_r094.csv")
    vec, cols = ST.loadVectors(path, "pc", True)
    tabs = ST.tables()
    mk = {"A": lambda: ST.condition(tabs["A"], {c: vec[a] for c, a in ST.REPS.items()}),
          "B": lambda: ST.condition(tabs["B"], vec),
          "C": lambda: ST.condition(tabs["C"], {c: vec[a] for c, a in ST.REPS.items()})}
    if stage.startswith("s1"):
        r = ST.codonShuffle(mk[stage[2]]())
    elif stage == "s2":
        r = {"A": ST.blockNull(mk["A"](), True), "C": ST.blockNull(mk["C"](), True),
             "B": ST.blockNull(mk["B"](), False, n=1_000_000)}
    elif stage == "s2prim":
        r = {"B_1e6": ST.blockNull(mk["B"](), False, n=1_000_000)}
    elif stage == "s3":
        r = ST.part3and1(vec, do_part1=False)
    elif stage == "s4":
        ST.ds.EB = os.path.join(HERE, "decode_r094")
        os.makedirs(ST.ds.EB, exist_ok=True)
        r = ST.part3and1(vec, do_part1=True)["part1"]
    json.dump({"stage": stage, "vectors": path, "result": r, "seconds": round(time.time() - t0, 1)},
              open(os.path.join(HERE, f"part3_{stage}.json"), "w"), indent=1)
    print(stage, "done", round(time.time() - t0, 1), flush=True)


if __name__ == "__main__":
    main()
