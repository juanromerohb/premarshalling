#!/bin/bash
# Slurm job body: starts <slots> queue workers of <threads> cores each and waits for them.
#   sbatch [resources] job.sh <run_dir> <queue> <slots> <threads>
# Resources must provide slots*threads CPUs on ONE node (-N 1 -c slots*threads).
#SBATCH -p all
#SBATCH -N 1
#SBATCH -t 10-00:00:00
set -uo pipefail
RUN="$1"; QUEUE="$2"; SLOTS="$3"; THREADS="$4"
mkdir -p "$RUN/joblogs"
echo "job ${SLURM_JOB_ID:-?} array ${SLURM_ARRAY_TASK_ID:-0} on $(hostname -s) cpus=$(awk '/Cpus_allowed_list/ {print $2}' /proc/self/status) start $(date '+%F %T')"
# Run a private snapshot of the worker: editing worker.sh later never affects running jobs
# (bash reads scripts incrementally), and the exact script used stays recorded with the run.
WORKER="$RUN/joblogs/worker_${SLURM_JOB_ID:-x}.sh"
cp "${PM_WORKER:-$HOME/premarshalling/code/cluster/worker.sh}" "$WORKER"
for ((k = 0; k < SLOTS; k++)); do
  bash "$WORKER" "$RUN" "$QUEUE" "$k" "$THREADS" \
    > "$RUN/joblogs/worker_${SLURM_JOB_ID:-x}_${SLURM_ARRAY_TASK_ID:-0}_s${k}.log" 2>&1 &
  sleep 2
done
wait
echo "end $(date '+%F %T')"
