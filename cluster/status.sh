#!/bin/bash
# Progress of a run directory: status.sh <run_dir>
RUN="$1"
count() { ls "$RUN/state/$1" 2>/dev/null | awk -F'__' '{print ($1 ~ /^(CP1|CP2|IPLCT)$/ ? "model" : "single")}' | sort | uniq -c; }
declare -A TOT DONE FAIL CLAIM
for q in single model; do TOT[$q]=$( [ -f "$RUN/tasks_$q.tsv" ] && wc -l < "$RUN/tasks_$q.tsv" || echo 0 ); DONE[$q]=0; FAIL[$q]=0; CLAIM[$q]=0; done
while read -r n q; do DONE[$q]=$n; done < <(count done)
while read -r n q; do FAIL[$q]=$n; done < <(count failed)
while read -r n q; do CLAIM[$q]=$n; done < <(count claim)
for q in single model; do
  [ "${TOT[$q]}" -gt 0 ] || continue
  echo "$q: total=${TOT[$q]} done=${DONE[$q]} failed=${FAIL[$q]} running=$((CLAIM[$q] - DONE[$q] - FAIL[$q])) pending=$((TOT[$q] - CLAIM[$q]))"
done
echo "--- per method (done/failed)"
for st in done failed; do
  ls "$RUN/state/$st" 2>/dev/null | awk -F'__' -v s="$st" '{c[$1]++} END {for (k in c) print k, s, c[k]}'
done | sort
echo "--- failures (last 10)"
for f in $(ls -t "$RUN"/state/failed 2>/dev/null | head -10); do echo "$f: $(cat "$RUN/state/failed/$f")"; done
echo "--- jobs"
squeue -u "$USER" -o '%.12i %.10j %.3t %.10M %.4C %.6m %R'
