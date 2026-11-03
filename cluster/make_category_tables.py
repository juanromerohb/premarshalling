#!/usr/bin/env python3
"""Write the per-category tables of a results folder.

Usage: make_category_tables.py <results_dir> [--include-full]

For each <results_dir>/<METHOD>/all_summary.csv, writes <METHOD>/category-tables/tau{10,15,20}m.csv
with one row per (family, H, S) and one TOTAL row per family. The reduction of inaccessible containers
is undefined when the initial layout is already fully accessible (initial_acc == c), so those instances
are left out of that mean; --include-full counts them as 0% instead.
"""
import csv
import os
import sys
from collections import defaultdict

FAMILIES = ["CV", "CVsub", "BF", "BZ", "EMM", "ZJY"]
MODELS = {"CP1", "CP2", "IPLCT"}  # run on the CV subset, labelled CV in all_summary.csv
HEADER = ["tau", "banco", "H", "S", "num_instancias_categoria", "num_optimos", "acc_final_pct_prom",
          "acc_aumento_pp_prom", "inacc_reduccion_pct_prom", "tiempo_ejecucion_ms_prom",
          "nodos_explorados_prom", "nodos_explorados_optimos_prom"]


def mean(values):
    return sum(values) / len(values) if values else 0.0


def row(tau, family, h, s, rows, include_full):
    final, gain, reduction, ms, nodes, nodes_opt = [], [], [], [], [], []
    for r in rows:
        acc, initial, c = int(r["acc"]), int(r["initial_acc"]), int(r["c"])
        final.append(100.0 * acc / c)
        gain.append(100.0 * (acc - initial) / c)
        if initial < c:
            reduction.append(100.0 * (acc - initial) / (c - initial))
        elif include_full:
            reduction.append(0.0)
        ms.append(float(r["elapsed_ms"]))
        n = float(r["nodes"] or 0)
        nodes.append(n)
        if r["proved_optimal"] == "true":
            nodes_opt.append(n)
    return [f"{tau:.6f}", family, h, s, len(rows), len(nodes_opt)] + [
        f"{mean(v):.6f}" for v in (final, gain, reduction, ms, nodes, nodes_opt)]


def main():
    results_dir = sys.argv[1]
    include_full = "--include-full" in sys.argv[2:]
    for method in sorted(os.listdir(results_dir)):
        path = os.path.join(results_dir, method, "all_summary.csv")
        if not os.path.exists(path):
            continue
        groups = defaultdict(lambda: defaultdict(list))
        with open(path, newline="") as fh:
            for r in csv.DictReader(fh):
                family = r["instance"].split("_")[0]
                if method in MODELS and family == "CV":
                    family = "CVsub"
                tau = int(float(r["tau"])) // 60
                groups[tau][(family, int(r["h"]), int(r["s"]))].append(r)
        out_dir = os.path.join(results_dir, method, "category-tables")
        os.makedirs(out_dir, exist_ok=True)
        for tau, cats in sorted(groups.items()):
            with open(os.path.join(out_dir, f"tau{tau}m.csv"), "w", newline="") as fh:
                w = csv.writer(fh, lineterminator="\n")
                w.writerow(HEADER)
                for family in FAMILIES:
                    keys = sorted(k for k in cats if k[0] == family)
                    for k in keys:
                        w.writerow(row(tau, family, k[1], k[2], cats[k], include_full))
                    if keys:
                        w.writerow(row(tau, family, "TOTAL", "TOTAL", [r for k in keys for r in cats[k]], include_full))


if __name__ == "__main__":
    main()
