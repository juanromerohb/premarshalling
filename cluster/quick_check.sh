#!/bin/bash
# Run every algorithm once on a small instance (sanity check of the toolchain).
source "$HOME/premarshalling/cluster/env.sh"
export GUROBI_THREADS=4 CP_SAT_NUM_WORKERS=4
ls -la /usr/bin/time
W=$(mktemp -d "$HOME/premarshalling/work_quick.XXXX"); cd "$W"; ln -s "$PM_CODE/scripts" scripts
"$PM_CODE/target/release/bench" --instance-file "$PM_CODE/instances/CVsub/CVh5s3/h5s3n1.txt" --tau 10 \
  --algo glct,hlct,bblct,lct1-cpsat,lct2-cpsat,iplct --cp-time-limit 60 --ip-time-limit 60 \
  --bblct-cpu-limit-sec 60 --csv out.csv 2>&1 | tail -12
cat out.csv
head -5 tmp-agent/iplct/*.log
cd; rm -rf "$W"
