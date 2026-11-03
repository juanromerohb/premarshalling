#!/bin/bash
# Submit queue workers as a Slurm job array.
#   submit.sh <run_dir> <queue> <n_jobs> <slots_per_job> <threads_per_slot> <mem_gb_per_slot> [extra sbatch args]
# Each array element is one single-node job with slots*threads CPUs running `slots` workers.
# Examples:
#   submit.sh $RUN single 48 1 1 2          # 48 one-core jobs for 1-thread methods
#   submit.sh $RUN model 7 5 4 17           # 7 jobs x (5 runs x 4 threads) = 140 cores
set -euo pipefail
RUN="$1"; QUEUE="$2"; N="$3"; SLOTS="$4"; THREADS="$5"; MEM="$6"; shift 6
mkdir -p "$RUN/joblogs"
sbatch --parsable -J "pm-$QUEUE" -p all -N 1 -n 1 -c $((SLOTS * THREADS)) --mem="$((SLOTS * MEM))G" -m block:block \
  --array="1-$N" -o "$RUN/joblogs/slurm-%A_%a.out" "$@" \
  "$HOME/premarshalling/code/cluster/job.sh" "$RUN" "$QUEUE" "$SLOTS" "$THREADS"
