#!/usr/bin/env python3
"""Recompute the data tables of the paper from merged results.

Usage: make_tables.py <paper.tex> <results_dir> <out.tex> [--include-full]

Takes each table from <paper.tex> as a template (caption, label, headers, row layout), replaces
its numeric cells with values computed from <results_dir>/<METHOD>/all_summary.csv, writes the
tables to <out.tex> and prints how many cells differ from the template. The manuscript itself is
never modified.

Reduction is the mean per-instance reduction of inaccessible containers. It is undefined when the
initial layout is already fully accessible (initial_acc == c), so those instances are left out of
that mean only; #Optimal and CPU time still count them. --include-full counts them as 0% instead,
which is how the previous campaign's tables were computed.
"""
import csv
import os
import re
import sys
from decimal import Decimal, ROUND_HALF_UP

LABELS = ["tab:bb_adjust_summary", "tab:exact_combined", "tab:bb_all", "tab:heuristics",
          "tab:bb_adjust", "tab:bb_adjust_15", "tab:bb_adjust_20"]
BB = ["BB0", "BB1", "BB2", "BB3", "BB4"]
METHODS = BB + ["BBLCT", "CP1", "CP2", "IPLCT", "GLCT", "HLCT"]
# The CV subset of the exact comparison (Table exact_combined).
EXACT_CV = {(5, 3), (5, 4), (5, 5), (5, 6), (5, 7), (5, 8), (6, 4), (6, 5)}


def load(results_dir):
    data = {}
    for m in METHODS:
        with open(os.path.join(results_dir, m, "all_summary.csv"), newline="") as fh:
            rows = list(csv.DictReader(fh))
        for r in rows:
            r["set"] = r["instance"].split("_")[0]
            r["tau_i"] = int(float(r["tau"]))
            r["hs"] = (int(r["h"]), int(r["s"]))
        data[m] = rows
    return data


def fmt(x, digits):
    """Rounds half up, as in the published tables (1.125 -> 1.13)."""
    return str(Decimal(repr(round(x, 9))).quantize(Decimal(1).scaleb(-digits), rounding=ROUND_HALF_UP))


class Stats:
    def __init__(self, data, include_full):
        self.data = data
        self.include_full = include_full

    def rows(self, method, tau, fam, hs=None, subset=None):
        out = [r for r in self.data[method] if r["tau_i"] == tau and r["set"] == fam]
        if hs is not None:
            out = [r for r in out if r["hs"] == hs]
        if subset is not None:
            out = [r for r in out if r["hs"] in subset]
        return out

    def optimal(self, rows):
        return str(sum(r["proved_optimal"] == "true" for r in rows))

    def reduction(self, rows):
        vals = []
        for r in rows:
            a, ia, c = int(r["acc"]), int(r["initial_acc"]), int(r["c"])
            if ia == c:
                if self.include_full:
                    vals.append(0.0)
                continue
            vals.append(100.0 * (a - ia) / (c - ia))
        return fmt(sum(vals) / len(vals), 2) if vals else "---"

    def seconds(self, rows):
        return fmt(sum(float(r["elapsed_ms"]) for r in rows) / len(rows) / 1000, 1)

    def millis(self, rows):
        return fmt(sum(float(r["elapsed_ms"]) for r in rows) / len(rows), 4)


def split_row(line):
    body = line.strip()
    assert body.endswith("\\\\")
    return [c.strip() for c in body[:-2].split("&")]


def is_total(cells):
    return any("Total" in c for c in cells)


def key_cells(cells, k):
    """Returns (hs, N) of a data row with k trailing value cells; hs is None for a Total row."""
    if is_total(cells):
        return None, cells[-k - 1]
    return (int(cells[-k - 3]), int(cells[-k - 2])), cells[-k - 1]


def fill(block, k, values):
    """Replaces the last k cells of each data row; values(line_state, cells) -> list of k strings."""
    out = []
    state = {"set": None, "tau": None}
    changed = total = 0
    for line in block.split("\n"):
        m = re.search(r"\\multirow\{\d+\}\{\*\}\{\$\\tau=(\d+)\$s\}", line)
        if m:
            state["tau"] = int(m[1])
        m = re.search(r"\\multirow\{2\}\{\*\}\{(\d+)\}", line)
        if m:
            state["tau"] = int(m[1])
        m = re.search(r"\\multirow\{\d+\}\{\*\}\{(CV|BF|BZ|EMM|ZJY)\}", line)
        if m:
            state["set"] = m[1]
        s = line.strip()
        if s.endswith("\\\\") and s.count("&") >= k and not re.search(r"\\textsuperscript|\\multicolumn\{[345]\}|Reduction|#Optimal|\\\( H \\\)", s):
            cells = split_row(line)
            new = values(state, cells)
            for a, b in zip(cells[-k:], new):
                total += 1
                changed += a != b
            indent = line[: len(line) - len(line.lstrip())]
            line = indent + " & ".join(cells[:-k] + new) + " \\\\"
        out.append(line)
    return "\n".join(out), changed, total


def build(block, label, st):
    if label == "tab:bb_all":
        def values(state, cells):
            hs, _ = key_cells(cells, 9)
            res = []
            for tau in (600, 900, 1200):
                r = st.rows("BBLCT", tau, state["set"], hs)
                res += [st.optimal(r), st.reduction(r), st.seconds(r)]
            return res
        return fill(block, 9, values)

    if label == "tab:heuristics":
        def values(state, cells):
            hs, _ = key_cells(cells, 15)
            res = []
            for tau in (600, 900, 1200):
                g = st.rows("GLCT", tau, state["set"], hs)
                h = st.rows("HLCT", tau, state["set"], hs)
                b = st.rows("BBLCT", tau, state["set"], hs)
                res += [st.reduction(g), st.reduction(h), st.reduction(b), st.millis(g), st.millis(h)]
            return res
        return fill(block, 15, values)

    if label in ("tab:bb_adjust", "tab:bb_adjust_15", "tab:bb_adjust_20"):
        tau = {"tab:bb_adjust": 600, "tab:bb_adjust_15": 900, "tab:bb_adjust_20": 1200}[label]

        def values(state, cells):
            hs, _ = key_cells(cells, 15)
            rs = [st.rows(m, tau, state["set"], hs) for m in BB]
            return [st.optimal(r) for r in rs] + [st.reduction(r) for r in rs] + [st.seconds(r) for r in rs]
        return fill(block, 15, values)

    if label == "tab:bb_adjust_summary":
        def values(state, cells):
            fam = cells[-17]
            rs = [st.rows(m, state["tau"], fam) for m in BB]
            return [st.optimal(r) for r in rs] + [st.reduction(r) for r in rs] + [st.seconds(r) for r in rs]
        return fill(block, 15, values)

    if label == "tab:exact_combined":
        def values(state, cells):
            hs, _ = key_cells(cells, 12)
            subset = EXACT_CV if state["set"] == "CV" else None
            res = []
            rs = []
            for m in ("CP1", "CP2", "IPLCT", "BBLCT"):
                if m == "CP2" and state["set"] == "ZJY":
                    rs.append(None)
                else:
                    rs.append(st.rows(m, state["tau"], state["set"], hs, subset))
            for fn in (st.optimal, st.reduction, st.seconds):
                res += ["---" if r is None else fn(r) for r in rs]
            return res
        return fill(block, 12, values)

    raise ValueError(label)


def main():
    tex_path, results_dir, out_path = sys.argv[1:4]
    include_full = "--include-full" in sys.argv[4:]
    with open(tex_path, encoding="utf-8") as fh:
        tex = fh.read()
    st = Stats(load(results_dir), include_full)
    blocks = {}
    for b in re.findall(r"\\begin\{table\}.*?\\end\{table\}", tex, flags=re.S):
        m = re.search(r"\\label\{(tab:[^}]+)\}", b)
        if m:
            blocks[m[1]] = b
    with open(out_path, "w", encoding="utf-8", newline="\n") as out:
        out.write(f"% Generated by cluster/make_tables.py from {os.path.basename(os.path.normpath(results_dir))}.\n")
        out.write("% Reduction means include fully accessible instances as 0%.\n\n" if include_full
                  else "% Reduction means exclude instances with initial_acc == c.\n\n")
        for label in LABELS:
            new, changed, total = build(blocks[label], label, st)
            print(f"{label}: {changed} of {total} cells differ from the template")
            out.write(new + "\n\n")


if __name__ == "__main__":
    main()
