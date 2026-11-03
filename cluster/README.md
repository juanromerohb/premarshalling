# Running the paper experiments on the UCLM cluster (Slurm)

Scripts used to run every computational experiment of the paper on the UCLM
supercomputing service (Slurm, nodes with 2x Intel Xeon E5-2660 v4,
28 cores, 128 GB RAM, Oracle Linux 9).

## Software (installed in the user's home by `install.sh`)

| Component | Version |
|-----------|---------|
| Rust      | 1.99.1 (`rustup`) |
| Python    | 3.12 (standalone, `uv`) |
| OR-Tools CP-SAT | 9.15.6755 |
| Gurobi (`gurobi_cl`) | 13.0.3, WLS academic license (`GRB_LICENSE_FILE`) |

`env.sh` puts everything on `PATH`. The license file is **not** part of the repository.

The scripts expect the repository at `~/premarshalling/code` and read their environment from
`~/premarshalling/cluster`, so link that folder after cloning:

```bash
ln -s ~/premarshalling/code/cluster ~/premarshalling/cluster
```

Build (on a compute node, so `target-cpu=native` matches the experiment CPUs):

```bash
cd ~/premarshalling/code && cargo fetch --target x86_64-unknown-linux-gnu
sbatch cluster/build.sh      # release + LTO + codegen-units=1 + target-cpu=native
```

## Experiment design

One task = one `(method, instance, tau)` = one `bench --instance-file` process, so a crash
or panic only affects that run and is recorded. `make_tasks.py` writes two task lists:

| Method | Families | Limit (s) | Queue |
|--------|----------|-----------|-------|
| GLCT, HLCT | CV, BF, BZ, EMM, ZJY | — | single (1 thread) |
| BB0–BB4 (component analysis) | CV, BF | 595 | single |
| BBLCT | CV, BF, BZ, EMM, ZJY | 3595 | single |
| CP1 (LCT1, CP-SAT) | CVsub, ZJY | 3595 | model (4 threads) |
| CP2 (LCT2, CP-SAT) | CVsub | 3595 | model |
| IPLCT (Gurobi) | CVsub, ZJY | 3595 | model |

All with tau = 10, 15, 20 min. Limits = 600/3600 s minus the 5 s margin described in the paper.
BB variants: BB0 = plain tree search; BB1 = +lower-bound pruning; BB2 = +dominance;
BB3 = +node ordering (LB sorting); BB4 = +HLCT initial incumbent (= BBLCT).

Per run: 16 GB memory (hard address-space limit for the Rust methods; `max_memory_in_mb=16384`
for CP-SAT and `MemLimit=16` for Gurobi), exactly 1 or 4 dedicated cores (`taskset`), solver
threads fixed through `GUROBI_THREADS` / `CP_SAT_NUM_WORKERS`, and a hard wall-clock guard
(`timeout`, limit + 900 s) that is never reached in normal runs.

## Execution

Workers claim tasks dynamically (`mkdir` on NFS is atomic), longest tasks first. Each Slurm
array element is a single-node job running several workers ("slots").

```bash
source ~/premarshalling/cluster/env.sh; cd $PM_CODE
RUN=~/premarshalling/runs/rerun_2026-10
python3 cluster/make_tasks.py $PM_CODE $RUN
bash cluster/submit.sh $RUN model 7 5 4 6     # 7 jobs x 5 slots x 4 threads, 6 GB/slot requested (peak seen: 2.4 GB)
bash cluster/submit.sh $RUN single 40 1 1 2   # 40 one-core jobs
bash cluster/status.sh $RUN                   # progress
bash cluster/inspect.sh $RUN                  # memory / CPU / stop-reason summary
```

Recovery: `requeue.sh $RUN stale` (claims of cancelled/dead jobs) and `requeue.sh $RUN failed`,
then submit workers again. Finished tasks are never re-run.

Gurobi WLS: tested with 24+ concurrent `gurobi_cl` processes on 4 nodes without license
errors; sessions are counted per machine, so keep the model jobs on at most ~8 nodes.

## Results

```bash
python3 cluster/merge_results.py $RUN $RUN/merged    # -> merged/<METHOD>/all_summary.csv
python3 cluster/make_category_tables.py $RUN/merged  # -> merged/<METHOD>/category-tables/
python3 cluster/make_tables.py <manuscript.tex> $RUN/merged tables.tex
```

`merged/` has the same layout and columns as `paper-results/`. `make_tables.py` takes each table of the
manuscript as a template, so it needs the LaTeX source of the paper, which is not part of this repository. When averaging the
reduction of inaccessible containers, exclude instances whose initial layout is already
fully accessible (`initial_acc == c`), for which the reduction is undefined.
