#!/bin/bash
# Make tasks claimable again.
#   requeue.sh <run_dir> failed   -> release every failed task
#   requeue.sh <run_dir> stale    -> release claims without done/failed whose owner job is no longer running
# Afterwards submit workers again (submit.sh); finished tasks are never re-run.
RUN="$1"; MODE="$2"
cd "$RUN/state" || exit 1
case "$MODE" in
  failed)
    for f in failed/*; do [ -e "$f" ] || continue; id=${f#failed/}
      mv "$f" "failed_history_${id}_$(date +%s)" 2>/dev/null || rm -f "$f"; rm -rf "claim/$id"; echo "released $id"; done
    mkdir -p history && mv failed_history_* history/ 2>/dev/null ;;
  stale)
    # Only RUNNING jobs count: a job requeued by Slurm (node failure) keeps its id but its
    # workers are dead, so its claims are stale.
    running=$(squeue -h -u "$USER" -t R -o '%A %i' | tr ' ' '\n' | sort -u)
    for d in claim/*; do [ -e "$d" ] || continue; id=${d#claim/}
      [ -e "done/$id" ] || [ -e "failed/$id" ] && continue
      job=$(cut -d_ -f1 < "$d/owner" 2>/dev/null)
      if ! grep -qx "$job" <<<"$running"; then rm -rf "$d"; echo "released stale $id (job $job)"; fi
    done ;;
  *) echo "mode must be failed|stale"; exit 1 ;;
esac
