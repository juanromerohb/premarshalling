#!/bin/bash
# Queue worker: repeatedly claims the next free task of <queue> and runs it.
#   worker.sh <run_dir> <queue: single|model> <slot> <threads_per_slot>
# A task is claimed atomically with `mkdir state/claim/<id>` (atomic on NFS), so any
# number of workers on any nodes can share one task list. Each task leaves exactly one of
# state/done/<id> or state/failed/<id>. The worker exits when nothing is left to claim.
set -uo pipefail
RUN="$1"; QUEUE="$2"; SLOT="$3"; THREADS="$4"
source "$HOME/premarshalling/cluster/env.sh"

BENCH="$PM_CODE/target/release/bench"
TASKS="$RUN/tasks_${QUEUE}.tsv"
mkdir -p "$RUN"/state/{claim,done,failed} "$RUN"/results "$RUN"/logs "$RUN"/solver_logs "$RUN"/work

# This slot's cores: slice of the CPUs Slurm bound to the job.
ALLOWED=$(awk '/Cpus_allowed_list/ {print $2}' /proc/self/status)
CPUS=$(python3 - "$ALLOWED" "$SLOT" "$THREADS" <<'EOF'
import sys
spec, slot, n = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
cpus = []
for part in spec.split(","):
    a, _, b = part.partition("-")
    cpus += range(int(a), int(b or a) + 1)
def socket(c):
    try:
        return int(open(f"/sys/devices/system/cpu/cpu{c}/topology/physical_package_id").read())
    except OSError:
        return 0
# Group by socket so a multi-threaded slot stays on one NUMA node whenever possible.
cpus.sort(key=lambda c: (socket(c), c))
sel = cpus[slot * n:(slot + 1) * n]
if len(sel) != n:
    sys.exit(f"slot {slot} x {n} threads does not fit in CPUs {spec}")
print(",".join(map(str, sel)))
EOF
) || exit 1

export GUROBI_THREADS="$THREADS"
export CP_SAT_NUM_WORKERS="$THREADS"
export CP_OPTIMIZER_WORKERS="$THREADS"
export OMP_NUM_THREADS="$THREADS"
TAG="${SLURM_JOB_ID:-nojob}_${SLURM_ARRAY_TASK_ID:-0}_s${SLOT}@$(hostname -s)"

# On scancel / job end: kill the running task's whole process group, leave its claim
# (it shows up as stale and is released by `requeue.sh <run> stale`).
CHILD=""
trap 'if [ -n "$CHILD" ]; then kill -TERM -- -"$CHILD" 2>/dev/null; kill -TERM "$CHILD" 2>/dev/null; sleep 2; kill -KILL -- -"$CHILD" 2>/dev/null; fi; echo "[$(date "+%F %T")] worker $TAG terminated by signal"; exit 143' TERM INT
echo "[$(date '+%F %T')] worker $TAG queue=$QUEUE cpus=$CPUS"

while IFS=$'\t' read -r ID METHOD FAMILY INST TAU LIMIT; do
  [ -e "$RUN/state/claim/$ID" ] && continue
  mkdir "$RUN/state/claim/$ID" 2>/dev/null || continue
  echo "$TAG $(date '+%F %T')" > "$RUN/state/claim/$ID/owner"

  case "$METHOD" in
    GLCT)  ALGO=(--algo glct) ;;
    HLCT)  ALGO=(--algo hlct) ;;
    BB0)   ALGO=(--algo bblct --bblct-init-heuristic none --bblct-lb-pruning false --bblct-dominance false --bblct-lb-sorting false) ;;
    BB1)   ALGO=(--algo bblct --bblct-init-heuristic none --bblct-lb-pruning true  --bblct-dominance false --bblct-lb-sorting false) ;;
    BB2)   ALGO=(--algo bblct --bblct-init-heuristic none --bblct-lb-pruning true  --bblct-dominance true  --bblct-lb-sorting false) ;;
    BB3)   ALGO=(--algo bblct --bblct-init-heuristic none --bblct-lb-pruning true  --bblct-dominance true  --bblct-lb-sorting true) ;;
    BB4|BBLCT) ALGO=(--algo bblct --bblct-init-heuristic hlct --bblct-lb-pruning true --bblct-dominance true --bblct-lb-sorting true) ;;
    CP1)   ALGO=(--algo lct1-cpsat) ;;
    CP2)   ALGO=(--algo lct2-cpsat) ;;
    IPLCT) ALGO=(--algo iplct) ;;
    *) echo "unknown method $METHOD" > "$RUN/state/failed/$ID"; continue ;;
  esac
  case "$METHOD" in
    BB*) ALGO+=(--bblct-cpu-limit-sec "$LIMIT") ;;
    CP*) ALGO+=(--cp-time-limit "$LIMIT") ;;
    IPLCT) ALGO+=(--ip-time-limit "$LIMIT") ;;
  esac

  # Memory: 16 GB per run. Slurm here has no cgroup enforcement, so 1-thread Rust runs get a
  # hard address-space limit; CP-SAT and Gurobi enforce it internally (max_memory_in_mb=16384,
  # MemLimit=16) and get a looser address-space guard (thread arenas inflate virtual size).
  case "$QUEUE" in
    single) VMEM_KB=$((16 * 1024 * 1024)) ;;
    *)      VMEM_KB=$((48 * 1024 * 1024)) ;;
  esac
  # Hard wall-clock guard well beyond the internal limit (catches hangs, never hit normally).
  HARD=$(( LIMIT > 0 ? LIMIT + ${PM_HARD_EXTRA:-900} : 600 ))

  WORK="$RUN/work/$ID"
  rm -rf "$WORK"; mkdir -p "$WORK"
  ln -s "$PM_CODE/scripts" "$WORK/scripts"   # the Rust code calls scripts/solve_lct_cpsat.py relative to cwd
  CSV="$RUN/results/$ID.csv"; rm -f "$CSV"
  LOG="$RUN/logs/$ID.log"

  START=$(date +%s)
  # Run in the background and `wait` so the TERM trap can fire mid-run. `timeout` puts the run in
  # its own process group, which Slurm (proctrack/pgid) would not kill on scancel.
  (
    cd "$WORK" || exit 97
    ulimit -v "$VMEM_KB"
    exec timeout -k 60 "$HARD" taskset -c "$CPUS" /usr/bin/time -v \
      "$BENCH" --instance-file "$PM_CODE/$INST" --tau "$TAU" "${ALGO[@]}" --csv "$CSV"
  ) < /dev/null > "$LOG" 2>&1 &
  CHILD=$!
  wait "$CHILD"
  RC=$?
  CHILD=""
  END=$(date +%s)

  # Keep solver logs (small), drop models / temporary files.
  for f in "$WORK"/tmp-agent/iplct/*.log; do [ -e "$f" ] && cp "$f" "$RUN/solver_logs/$ID.gurobi.log"; done
  rm -rf "$WORK"

  ROWS=$( [ -f "$CSV" ] && tail -n +2 "$CSV" | wc -l || echo 0 )
  INFO="rc=$RC rows=$ROWS wall=$((END - START))s host=$(hostname -s) tag=$TAG"
  if [ "$RC" -eq 0 ] && [ "$ROWS" -eq 1 ]; then
    echo "$INFO" > "$RUN/state/done/$ID"
  else
    echo "$INFO" > "$RUN/state/failed/$ID"
    echo "[$(date '+%F %T')] FAILED $ID $INFO"
  fi
done < "$TASKS"

echo "[$(date '+%F %T')] worker $TAG: queue exhausted"
