#!/bin/bash
#SBATCH -J pm-build
#SBATCH -p all
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem=16G
#SBATCH -t 60
#SBATCH -o build-%j.out
# Release build on a compute node (same CPU as the experiment nodes).
# Dependencies must be fetched beforehand on the login node: `cargo fetch`.
set -euo pipefail
source "$HOME/premarshalling/cluster/env.sh"
cd "$PM_CODE"

export RUSTFLAGS="-C target-cpu=native"
export CARGO_PROFILE_RELEASE_LTO=fat
export CARGO_PROFILE_RELEASE_CODEGEN_UNITS=1
export CARGO_PROFILE_RELEASE_PANIC=unwind

hostname
rustc --version
cargo build --release --offline --bin bench -j "${SLURM_CPUS_PER_TASK:-8}"
ls -la target/release/bench
