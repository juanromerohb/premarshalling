#!/bin/bash
# Probe a compute node (run through srun/sbatch).
echo "host: $(hostname)"
lscpu | grep -E 'Model name|^CPU\(s\)|Thread|Socket|NUMA node|MHz'
free -g | head -2
grep -m1 flags /proc/cpuinfo | tr ' ' '\n' | grep -E '^(avx|avx2|bmi1|bmi2|fma|popcnt)$' | tr '\n' ' '; echo
echo "governor: $(cat /sys/devices/system/cpu/cpu0/cpufreq/scaling_governor 2>/dev/null)"
echo "no_turbo: $(cat /sys/devices/system/cpu/intel_pstate/no_turbo 2>/dev/null)"
echo "affinity: $(taskset -cp $$)"
ls -d /opt/amaconda3 2>&1
df -h /tmp | tail -1
ulimit -a | grep -E 'virtual|open files|processes'
env | grep -E '^SLURM_(CPUS|MEM|JOB_ID|NTASKS)' | sort
