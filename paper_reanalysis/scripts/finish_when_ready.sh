#!/bin/bash
# waits for the grid to land, builds the figure + table, and commits the result
# so the work survives a lost connection.  isolated worktree, own branch.
cd /storage/kiran-stuff/triplet-proof-reanalysis || exit 1
R=paper_reanalysis/results
for _ in $(seq 1 720); do            # up to 6 h
  [ -f "$R/grid.csv" ] && break
  sleep 30
done
if [ ! -f "$R/grid.csv" ]; then
  echo "grid.csv never appeared" > "$R/STATUS.txt"; exit 1
fi
python paper_reanalysis/scripts/make_report.py > "$R/report_stdout.txt" 2>&1
{
  echo "# Reanalysis status"
  echo
  echo "grid completed: $(date -r $R/grid.csv '+%Y-%m-%d %H:%M:%S')"
  echo "report built:   $(date '+%Y-%m-%d %H:%M:%S')"
  echo
  echo '## z grid'; echo
  cat "$R/TABLE_z_grid.md" 2>/dev/null
  echo; echo '## cross-validation vs the preregistered 10^7 block null'; echo
  cat "$R/crossvalidation.json" 2>/dev/null
} > "$R/STATUS.md"
git add -A
git -c user.name=khatvangi -c user.email=kinnaa@gmail.com commit -q -F - <<'MSG'
Reanalysis grid: results, figure and table

Produced unattended by finish_when_ready.sh after run_reanalysis.py completed.
Two nulls x five property representations x two objectives, one code path,
gated on reproducing the published SGC values (E to 4.5e-13, D to 2.5e-14).

See results/STATUS.md for the z grid and the cross-validation against the
preregistered 10^7 block-null run.

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_014bb6ZYvqeY4F5KtUs1Xruz
MSG
echo done >> "$R/STATUS.txt"
