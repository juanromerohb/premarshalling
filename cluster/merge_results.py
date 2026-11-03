#!/usr/bin/env python3
"""Merge the per-task CSVs of a run into one all_summary.csv per method.

Usage: merge_results.py <run_dir> <out_dir>

Writes <out_dir>/<METHOD>/all_summary.csv with the same columns as the bench CSV and the
instance name prefixed by its family (CV_h5s3n1.txt, BF_BF10_h8s16n1.txt, ...), matching
the layout of paper-results/. CVsub instances are labelled CV (they are CV instances).
Also reports tasks that have no result, so nothing silently goes missing.
"""
import csv
import os
import sys
from collections import defaultdict


def main():
    run_dir, out_dir = sys.argv[1], sys.argv[2]
    tasks = []
    for q in ("single", "model"):
        path = os.path.join(run_dir, f"tasks_{q}.tsv")
        if os.path.exists(path):
            with open(path) as fh:
                tasks += [line.rstrip("\n").split("\t") for line in fh if line.strip()]

    rows = defaultdict(list)
    header = None
    missing = []
    for task_id, method, family, rel, tau, limit in tasks:
        path = os.path.join(run_dir, "results", f"{task_id}.csv")
        if not os.path.exists(os.path.join(run_dir, "state", "done", task_id)) or not os.path.exists(path):
            missing.append(task_id)
            continue
        with open(path, newline="") as fh:
            r = list(csv.reader(fh))
        header = r[0]
        for row in r[1:]:
            label = "CV" if family == "CVsub" else family
            row[0] = f"{label}_{row[0]}"
            rows[method].append(row)

    for method, rs in sorted(rows.items()):
        d = os.path.join(out_dir, method)
        os.makedirs(d, exist_ok=True)
        rs.sort(key=lambda r: (float(r[6]), r[0]))
        with open(os.path.join(d, "all_summary.csv"), "w", newline="") as fh:
            w = csv.writer(fh, lineterminator="\n")
            w.writerow(header)
            w.writerows(rs)
        print(f"{method}: {len(rs)} rows")

    print(f"missing results: {len(missing)}")
    for t in missing[:50]:
        print("  ", t)


if __name__ == "__main__":
    main()
