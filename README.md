# CPMP-LCT in Rust

Rust implementation of exact and heuristic methods for the **Container Premarshalling Problem under Limited Crane Time (CPMP-LCT)**.

This repository accompanies the paper:

> Juan Romero del Hombrebueno, Consuelo Parreño-Torres.  
> *A branch-and-bound algorithm for the container premarshalling problem under limited crane time*.

The codebase includes:

- `BBLCT`: branch-and-bound algorithm for CPMP-LCT
- `GLCT` and `HLCT`: constructive and hybrid heuristics
- `IPLCT`: integer programming model
- `LCT1` and `LCT2`: constraint programming models with CP Optimizer and CP-SAT backends
- A desktop GUI to inspect instances and solutions
- Benchmark instances and paper result tables

## Repository Layout

- `src/`: Rust implementation of algorithms, models, utilities, and GUI
- `src/bin/bench.rs`: batch benchmark runner
- `src/bin/gui.rs`: desktop GUI entrypoint
- `scripts/solve_lct_cpsat.py`: OR-Tools CP-SAT wrapper
- `instances/`: benchmark instances used in the experiments
- `paper-results/`: CSV summaries and category-level tables reported in the paper

## Requirements

Minimum:

- Rust toolchain (`cargo`)

Optional, depending on what you want to run:

- `python3` + `ortools` for CP-SAT runs
- IBM `cpoptimizer` on `PATH` for CP Optimizer runs
- `gurobi_cl` on `PATH` for IPLCT runs

The repository was developed for the experimental setup described in the paper, where the Rust code was run with Rust 1.93.1.

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

## Launch the GUI

This repository also includes a precompiled Windows executable:

- `gui.exe`: ready-to-use desktop GUI for Windows

If you are on Windows and only want to explore instances or run algorithms interactively, you can launch `gui.exe` directly without building the Rust project first.

### Windows Quick Start

1. Run `gui.exe`.
2. In `Input`, choose one of:
   - `Manual`: create an instance by editing the grid
   - `Single .txt file`: open one benchmark instance
   - `Folder of .txt files`: run a batch over a directory
3. In `Algorithm and Limits`, choose `BBLCT`, `HLCT`, or `GLCT`, and set `Tau (minutes)`.
   - If you choose `BBLCT`, also set `BBLCT CPU limit (seconds)`.
4. Click `Run`.

### Brief GUI Guide

- In `Manual` mode, fill the `Manual Grid Editor` with container group values.
  - `Tier 1` is the bottom tier.
  - `0` means an empty slot.
- In `Single .txt file` mode, select an instance file from `instances/`.
- In `Folder of .txt files` mode, the GUI processes all `.txt` instances in the chosen directory and shows a results table.
- During a `BBLCT` run, the `Cancel BBLCT` button stops the search early.
- In manual mode, the `Results` panel includes a step-by-step replay with `Prev`, `Next`, and the `Step` slider.
- In file/folder modes, structured outputs are written under `results/`.

### Build From Source

```bash
cargo run --release --bin gui
```

On Linux, the GUI is configured to use `x11` and software OpenGL by default. If startup fails, try:

```bash
WINIT_UNIX_BACKEND=x11 LIBGL_ALWAYS_SOFTWARE=1 cargo run --release --bin gui
```

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

The `paper-results/` directory contains summary tables and category-level CSV files generated for the computational study. The benchmark runner can be used to reproduce or extend those experiments.

Because CP and IP runs depend on external solver installations and licenses, reproducibility of the full comparison requires matching solver availability in addition to the Rust environment.

## Citation

If you use this repository, please cite the associated paper. A finalized bibliographic entry will be added once publication details are available.

## License

See `LICENSE.txt`.
