#!/usr/bin/env python3
"""Generate the task lists for the paper experiments.

One task = one (method, instance, tau) run = one `bench --instance-file` process.
Output: <run_dir>/tasks_single.tsv (1-thread methods) and <run_dir>/tasks_model.tsv
(4-thread CP-SAT / Gurobi models). Columns:
    task_id  method  family  instance_rel_path  tau_min  limit_sec

Tasks are sorted longest-expected-first (time limit, tau, bay size) so that the
dynamic queue behaves like LPT scheduling.

Usage: make_tasks.py <code_dir> <run_dir> [--only METHOD[,METHOD...]] [--families F,...]
                     [--taus 10,15,20] [--max-per-category N] [--instances-file FILE]
"""
import argparse
import os
import re
import sys

TAUS = [10, 15, 20]
ALL_FAMILIES = ["CV", "BF", "BZ", "EMM", "ZJY"]

# method -> (queue, families, limit_sec)
METHODS = {
    "GLCT": ("single", ALL_FAMILIES, 0),
    "HLCT": ("single", ALL_FAMILIES, 0),
    "BB0": ("single", ["CV", "BF"], 595),
    "BB1": ("single", ["CV", "BF"], 595),
    "BB2": ("single", ["CV", "BF"], 595),
    "BB3": ("single", ["CV", "BF"], 595),
    "BB4": ("single", ["CV", "BF"], 595),
    "BBLCT": ("single", ALL_FAMILIES, 3595),
    "CP1": ("model", ["CVsub", "ZJY"], 3595),
    "CP2": ("model", ["CVsub"], 3595),
    "IPLCT": ("model", ["CVsub", "ZJY"], 3595),
}

HS_RE = re.compile(r"h(\d+)s(\d+)")


def list_instances(code_dir, family):
    root = os.path.join(code_dir, "instances", family)
    out = []
    for dirpath, _, files in os.walk(root):
        for f in files:
            if f.endswith(".txt"):
                out.append(os.path.relpath(os.path.join(dirpath, f), code_dir))
    out.sort()
    return out


def category(rel):
    m = HS_RE.search(os.path.basename(rel))
    return (int(m.group(1)), int(m.group(2)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("code_dir")
    ap.add_argument("run_dir")
    ap.add_argument("--only", default="")
    ap.add_argument("--families", default="")
    ap.add_argument("--taus", default="")
    ap.add_argument("--max-per-category", type=int, default=0)
    ap.add_argument("--instances-file", default="",
                    help="restrict to these 'METHOD TAU REL_PATH' lines")
    ap.add_argument("--limit-override", type=int, default=0,
                    help="override every non-zero time limit (tests only)")
    args = ap.parse_args()

    methods = [m for m in METHODS if not args.only or m in args.only.split(",")]
    fam_filter = set(args.families.split(",")) if args.families else None
    taus = [int(t) for t in args.taus.split(",")] if args.taus else TAUS
    wanted = None
    if args.instances_file:
        wanted = set()
        with open(args.instances_file) as fh:
            for line in fh:
                if line.strip() and not line.startswith("#"):
                    m, t, p = line.split()
                    wanted.add((m, int(t), p))

    queues = {"single": [], "model": []}
    for method in methods:
        queue, families, limit = METHODS[method]
        if args.limit_override and limit:
            limit = args.limit_override
        for fam in families:
            if fam_filter and fam not in fam_filter:
                continue
            insts = list_instances(args.code_dir, fam)
            if args.max_per_category:
                per_cat = {}
                kept = []
                for rel in insts:
                    key = (os.path.dirname(rel), category(rel))
                    per_cat[key] = per_cat.get(key, 0) + 1
                    if per_cat[key] <= args.max_per_category:
                        kept.append(rel)
                insts = kept
            for tau in taus:
                for rel in insts:
                    if wanted is not None and (method, tau, rel) not in wanted:
                        continue
                    h, s = category(rel)
                    stem = os.path.splitext(rel.replace("instances/", "", 1))[0]
                    task_id = f"{method}__t{tau}__{stem.replace('/', '_')}"
                    queues[queue].append((limit, tau, h * s, task_id, method, fam, rel))

    os.makedirs(args.run_dir, exist_ok=True)
    for queue, rows in queues.items():
        rows.sort(key=lambda r: (-r[0], -r[1], -r[2], r[3]))
        path = os.path.join(args.run_dir, f"tasks_{queue}.tsv")
        with open(path, "w") as fh:
            for limit, tau, _, task_id, method, fam, rel in rows:
                fh.write(f"{task_id}\t{method}\t{fam}\t{rel}\t{tau}\t{limit}\n")
        print(f"{path}: {len(rows)} tasks", file=sys.stderr)


if __name__ == "__main__":
    main()
