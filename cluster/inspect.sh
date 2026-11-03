#!/bin/bash
# Summaries of a run's logs: inspect.sh <run_dir>
R="$1"; cd "$R/logs" || exit 1
echo "--- peak RSS (MB) top 5 / per method max"
for f in *.log; do
  kb=$(grep -m1 'Maximum resident' "$f" | awk '{print $6}')
  echo "${f%%__*} $((kb / 1024)) ${f%.log}"
done | sort -k2,2n | tail -5
for f in *.log; do
  kb=$(grep -m1 'Maximum resident' "$f" | awk '{print $6}')
  echo "${f%%__*} $((kb / 1024))"
done | awk '{if ($2 > m[$1]) m[$1] = $2} END {for (k in m) print k, m[k] " MB"}' | sort
echo "--- CPU usage (Percent of CPU) per method, mean"
for f in *.log; do
  p=$(grep -m1 'Percent of CPU' "$f" | awk '{print $NF}' | tr -d '%')
  echo "${f%%__*} ${p:-0}"
done | awk '{s[$1] += $2; n[$1]++} END {for (k in s) printf "%s %.0f%%\n", k, s[k] / n[k]}' | sort
echo "--- wall (time -v) vs elapsed_ms reported, worst overhead"
for f in *.log; do
  id=${f%.log}; csv="$R/results/$id.csv"; [ -f "$csv" ] || continue
  w=$(grep -m1 'Elapsed (wall' "$f" | awk '{print $NF}' | awk -F: '{if (NF == 3) print $1 * 3600 + $2 * 60 + $3; else print $1 * 60 + $2}')
  e=$(tail -1 "$csv" | cut -d, -f14)
  echo "$id $w $e"
done | awk '{d = $2 - $3 / 1000; print d, $0}' | sort -n | tail -4
echo "--- stop reasons per method"
for c in "$R"/results/*.csv; do tail -1 "$c" | cut -d, -f2,17; done | sort | uniq -c
echo "--- slot CPU lists"
grep -h 'cpus=' "$R"/joblogs/worker_*.log | awk '{print $3, $5}' | sort -u | head -20
echo "--- gurobi license lines"
grep -h -iE 'license|WLS' "$R"/solver_logs/*.log 2>/dev/null | grep -v '^Set parameter' | sort | uniq -c | head
echo "--- errors in logs"
grep -l -iE 'panic|error|killed|Traceback' *.log 2>/dev/null | head
