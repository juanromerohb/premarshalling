# CPMP-LCT in Rust

Rust implementation of exact and heuristic methods for the **Container Premarshalling Problem under Limited Crane Time (CPMP-LCT)**.

This repository accompanies the paper:

> Juan Romero-del-Hombrebueno, Consuelo Parreño-Torres.  
> *A branch-and-bound algorithm for the container premarshalling problem under limited crane time*.

The codebase includes:

- `BBLCT`: branch-and-bound algorithm for CPMP-LCT
- `GLCT` and `HLCT`: constructive and hybrid heuristics
- `IPLCT`: integer programming model
- `LCT1` and `LCT2`: constraint programming models with CP Optimizer and CP-SAT backends
- Benchmark instances and paper result tables

## Repository Layout

- `src/`: Rust implementation of algorithms, models, and utilities
- `src/bin/bench.rs`: batch benchmark runner
- `src/bin/bound_table.rs`: compares UB(tau) and LB(n) at the initial layout with the BBLCT solutions
- `scripts/solve_lct_cpsat.py`: OR-Tools CP-SAT wrapper
- `scripts/make_category_tables.py`: builds the category-level tables from the summary CSV files
- `instances/`: benchmark instances used in the experiments; `instances/CVsub/` is the CV subset of the comparison between exact methods
- `paper-results/`: CSV summaries and category-level tables reported in the paper
- `paper-results/bound-analysis/`: per-run values and LaTeX table written by `bound_table`

## Requirements

Minimum:

- Rust toolchain (`cargo`)

Optional, depending on what you want to run:

- `python3` + `ortools` for CP-SAT runs
- IBM `cpoptimizer` on `PATH` for CP Optimizer runs
- `gurobi_cl` on `PATH` for IPLCT runs

The repository was developed for the experimental setup described in the paper, where the code was run with Rust 1.99.1, OR-Tools CP-SAT 9.15.6755 and Gurobi 13.0.3.

## Build

```bash
cargo build --release
```

## Run a Single Instance

The default binary runs `BBLCT` on one instance:

```bash
cargo run --release -- instances/CV/CVh5s3/h5s3n1.txt 15
```

Arguments:

- `instance_path`: path to a `.txt` instance
- `tau_minutes`: crane-time limit in minutes; default is `15`

Example output includes instance size, initial accessibility, final accessibility, move count, crane time, and explored nodes.

## Run Benchmarks

Use the benchmark runner for one instance, one directory, or complete benchmark campaigns.

Example: run only the Rust methods on one instance:

```bash
cargo run --release --bin bench -- \
  --instance-file instances/CV/CVh5s3/h5s3n1.txt \
  --algo glct,hlct,bblct \
  --tau 15
```

Example: run a benchmark directory and save a CSV summary:

```bash
cargo run --release --bin bench -- \
  --dir instances/CV \
  --algo bblct \
  --tau 15 \
  --csv results/cv_bblct_tau15.csv
```

Useful options:

- `--algo`: `bblct`, `glct`, `hlct`, `tgh`, `lct1-cpo`, `lct2-cpo`, `lct1-cpsat`, `lct2-cpsat`, `iplct`
- `--bblct-cpu-limit-sec`: internal time limit for BBLCT
- `--cp-time-limit`: time limit for CP models
- `--ip-time-limit`: time limit for IPLCT
- `--model-warm-start`: `none`, `glct`, `hlct`
- `--save-solutions`, `--save-results`, `--output`, `--csv`

Note: `--algo all` requires every external solver backend to be installed. If you only want the Rust-native methods, specify them explicitly.

## Instance Format

Instances are plain text matrices. Each column is a stack, `0` denotes an empty slot, and rows are listed from top to bottom.

Example:

```text
0  0  0
0  0  0
1  5  4
7  6  9
3  2  8
```

## Included Benchmarks

The repository contains benchmark families used in the paper, including:

- `CV`
- `BF`
- `BZ`
- `EMM`
- `ZJY`

## Reproducibility

The `paper-results/` directory contains summary tables and category-level CSV files generated for the computational study. The benchmark runner can be used to reproduce or extend those experiments, using the settings described in the paper.

The results come from a complete re-run after the corrections listed in `CHANGELOG.md`. The results of the previous version remain available at the `v1.0.0` tag. The average reduction of inaccessible containers excludes instances whose initial layout is already fully accessible, since it is undefined for them. The category tables of the BB variants only cover CV and BF, the families on which those variants were run.

To recompute the comparison between the bounds and the BBLCT solutions:

```bash
cargo run --release --bin bound_table -- \
  --results paper-results/BBLCT/all_summary.csv \
  --instances instances \
  --out-csv results/bounds_by_run.csv \
  --out-tex results/table_bounds.tex
```

Because CP and IP runs depend on external solver installations and licenses, reproducibility of the full comparison requires matching solver availability in addition to the Rust environment.

## Citation

If you use this repository, please cite the associated paper. A finalized bibliographic entry will be added once publication details are available. The code and results are archived in Zenodo at https://doi.org/10.5281/zenodo.19361481, which always resolves to the latest version.

## License

See `LICENSE.txt`.
