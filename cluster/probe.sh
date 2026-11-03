#!/bin/bash
# Quick environment probe (run on login node or inside a job).
echo "host: $(hostname)"
lscpu | grep -E 'Model name|^CPU\(s\)|Thread|Socket|NUMA node|MHz'
free -g | head -2
ulimit -a | grep -E 'virtual|max memory|open files|processes'
for pkg in ortools gurobipy; do
  curl -s "https://pypi.org/pypi/$pkg/json" | python3 -c 'import sys,json; d=json.load(sys.stdin)["info"]; print(d["name"], d["version"], d["requires_python"])'
done
curl -s https://static.rust-lang.org/dist/channel-rust-stable.toml | grep -m1 -A2 '^\[pkg.rustc\]'
