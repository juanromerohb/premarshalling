#!/bin/bash
# One-time installation of the toolchain in the user's home (NFS, visible from all nodes).
#   Rust (rustup, pinned toolchain), standalone Python + OR-Tools venv, Gurobi (gurobi_cl).
# Run on the login node (needs internet). Idempotent.
set -euo pipefail

BASE="$HOME/premarshalling"
RUST_VERSION="1.99.1"
PY_VERSION="3.12"
ORTOOLS_VERSION="9.15.6755"
GUROBI_VERSION="13.0.3"
GUROBI_DIR="$BASE/gurobi/gurobi${GUROBI_VERSION//./}"

mkdir -p "$BASE"

# --- Rust ---
if ! command -v "$HOME/.cargo/bin/rustup" >/dev/null 2>&1; then
  curl -sSf https://sh.rustup.rs | sh -s -- -y --profile minimal --default-toolchain "$RUST_VERSION" --no-modify-path
fi
"$HOME/.cargo/bin/rustup" toolchain install "$RUST_VERSION" --profile minimal
"$HOME/.cargo/bin/rustup" default "$RUST_VERSION"

# --- Python + OR-Tools ---
if [ ! -x "$HOME/.local/bin/uv" ]; then
  curl -LsSf https://astral.sh/uv/install.sh | env UV_NO_MODIFY_PATH=1 sh
fi
export UV_PYTHON_INSTALL_DIR="$BASE/python"
"$HOME/.local/bin/uv" python install "$PY_VERSION"
if [ ! -x "$BASE/venv/bin/python3" ]; then
  "$HOME/.local/bin/uv" venv --python "$PY_VERSION" "$BASE/venv"
fi
"$HOME/.local/bin/uv" pip install --python "$BASE/venv/bin/python3" "ortools==$ORTOOLS_VERSION"

# --- Gurobi ---
if [ ! -x "$GUROBI_DIR/linux64/bin/gurobi_cl" ]; then
  tmp=$(mktemp -d)
  curl -sSfL -o "$tmp/gurobi.tgz" "https://packages.gurobi.com/${GUROBI_VERSION%.*}/gurobi${GUROBI_VERSION}_linux64.tar.gz"
  mkdir -p "$GUROBI_DIR"
  tar -xzf "$tmp/gurobi.tgz" -C "$GUROBI_DIR" --strip-components=1
  rm -rf "$tmp"
fi

echo "--- installed ---"
"$HOME/.cargo/bin/rustc" --version
"$BASE/venv/bin/python3" -c 'import sys, ortools; print("python", sys.version.split()[0], "ortools", ortools.__version__)'
"$GUROBI_DIR/linux64/bin/gurobi_cl" --version | head -1
