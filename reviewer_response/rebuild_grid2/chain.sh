#!/bin/bash
cd /storage/kiran-stuff/triplet-proof
while kill -0 4063750 2>/dev/null; do sleep 2; done
OPT20_BUDGET=1 python reviewer_response/rebuild_grid2/opt20P.py > reviewer_response/rebuild_grid2/run_opt20P.log 2>&1
G_RANDOM=96 G_WARM=16 python reviewer_response/rebuild_grid2/grid_fixed.py P > reviewer_response/rebuild_grid2/run_gridP.log 2>&1
