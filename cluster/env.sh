# Environment shared by every job. Source it: `source ~/premarshalling/cluster/env.sh`
export PM_BASE="$HOME/premarshalling"
export PM_CODE="$PM_BASE/code"

# Rust
export PATH="$HOME/.cargo/bin:$PATH"

# Python venv with OR-Tools (the Rust code calls `python3` from PATH)
export PATH="$PM_BASE/venv/bin:$PATH"

# Gurobi (the Rust code calls `gurobi_cl` from PATH)
export GUROBI_HOME="$PM_BASE/gurobi/gurobi1303/linux64"
export PATH="$GUROBI_HOME/bin:$PATH"
export LD_LIBRARY_PATH="$GUROBI_HOME/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
export GRB_LICENSE_FILE="$PM_BASE/gurobi/gurobi.lic"
