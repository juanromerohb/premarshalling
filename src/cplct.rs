//! Constraint-programming models for the CPMP-LCT.
//!
//! The public API in this module exposes two solver variants: `Lct1`, which
//! minimizes the number of inaccessible containers, and `Lct2`, which
//! compresses the final objective to inaccessible priority groups when
//! priorities are unique. The implementation is split into two layers:
//! - paper-facing model construction (`build_model_data`, `generate_cpo_model`,
//!   and the transition tables);
//! - engineering wrappers that call CP Optimizer or a CP-SAT helper script,
//!   optionally inject warm starts, and reconstruct the resulting move sequence.

use crate::bay::{Bay, Move};
use crate::crane::CraneParams;
use crate::heuristics::{glct, hlct};
use crate::solution::Solution;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::fmt::Write;
use std::fs;
use std::path::Path;
use std::process::Command;
use std::sync::atomic::{AtomicU64, Ordering};
use std::time::{SystemTime, UNIX_EPOCH};

static UNIQUE_STAMP_COUNTER: AtomicU64 = AtomicU64::new(0);

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
/// CPMP-LCT model variant to solve.
pub enum CplctModel {
    /// Original LCT1 objective: minimize inaccessible containers.
    Lct1,
    /// Reduced LCT2 objective: minimize inaccessible priority groups.
    Lct2,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
/// Solver backend used to optimize the selected model.
pub enum CplctBackend {
    /// IBM CP Optimizer through an emitted `.cpo` model.
    CpOptimizer,
    /// OR-Tools CP-SAT through the Python wrapper in `scripts/solve_lct_cpsat.py`.
    CpSat,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
/// Heuristic used to provide a warm start to the CP backend.
pub enum WarmStartHeuristic {
    Glct,
    Hlct,
}

#[derive(Clone, Copy, Debug)]
/// Configuration of a CPMP-LCT solve.
pub struct CplctConfig {
    pub model: CplctModel,
    pub backend: CplctBackend,
    pub tau: f64,
    pub time_limit_sec: f64,
    pub warm_start_heuristic: Option<WarmStartHeuristic>,
}

#[derive(Clone, Debug)]
/// Result returned by [`solve_cplct`].
pub struct CplctResult {
    pub solution: Solution,
    pub optimal: bool,
    pub status: String,
    pub time_seconds: f64,
    pub objective_inaccessible: usize,
}

#[derive(Clone, Debug)]
struct ModelData {
    s: usize,
    h: usize,
    p: usize,
    c: usize,
    k: usize,
    init_x: Vec<Vec<usize>>,
    group_count: Vec<usize>,
    table_b: Vec<[i32; 8]>,
    h0_0_s: Vec<f64>,
    h0_rs: Vec<Vec<f64>>,
    h1_sr: Vec<Vec<f64>>,
    v_load: Vec<f64>,
    v_unload: Vec<f64>,
}

#[derive(Clone, Debug)]
struct BackendSolveData {
    status: String,
    optimal: bool,
    time_seconds: f64,
    objective_inaccessible: Option<usize>,
    y_ones: Vec<[usize; 3]>,
    z_ones: Vec<[usize; 3]>,
}

#[derive(Clone, Debug)]
struct WarmStartData {
    y_ones: Vec<[usize; 3]>,
    z_ones: Vec<[usize; 3]>,
}

#[derive(Serialize)]
struct CpsatInput {
    model: String,
    s: usize,
    h: usize,
    p: usize,
    c: usize,
    k: usize,
    tau: f64,
    time_limit_sec: f64,
    init_x: Vec<Vec<usize>>,
    group_count: Vec<usize>,
    table_b: Vec<[i32; 8]>,
    h0_0_s: Vec<f64>,
    h0_rs: Vec<Vec<f64>>,
    h1_sr: Vec<Vec<f64>>,
    v_load: Vec<f64>,
    v_unload: Vec<f64>,
    warm_start_y_ones: Vec<[usize; 3]>,
    warm_start_z_ones: Vec<[usize; 3]>,
}

#[derive(Deserialize)]
struct CpsatOutput {
    status: String,
    optimal: bool,
    time_seconds: f64,
    objective: Option<usize>,
    y_ones: Vec<[usize; 3]>,
    z_ones: Vec<[usize; 3]>,
}

/// Solves a CPMP-LCT instance with the selected model / backend combination.
///
/// The returned solution is reconstructed as a concrete relocation sequence and
/// scored again with the internal crane-time model. For `Lct2`, the bay must
/// have unique priorities because the compact objective only matches `Lct1`
/// under that assumption.
pub fn solve_cplct(bay: &Bay, crane: &CraneParams, cfg: &CplctConfig) -> CplctResult {
    if cfg.model == CplctModel::Lct2 {
        validate_lct2_unique_priorities(bay);
    }

    let data = build_model_data(bay, crane, cfg.tau);
    let warm_start = cfg
        .warm_start_heuristic
        .map(|heuristic| build_warm_start_data(bay, crane, cfg.tau, heuristic, data.k));

    let backend_data = match cfg.backend {
        CplctBackend::CpOptimizer => solve_with_cpoptimizer(&data, cfg, warm_start.as_ref()),
        CplctBackend::CpSat => solve_with_cpsat(&data, cfg, warm_start.as_ref()),
    };

    let solution = reconstruct_solution(
        bay,
        crane,
        data.k,
        &backend_data.y_ones,
        &backend_data.z_ones,
    );

    let objective_inaccessible = backend_data
        .objective_inaccessible
        .unwrap_or_else(|| bay.c.saturating_sub(solution.acc));

    CplctResult {
        solution,
        optimal: backend_data.optimal,
        status: backend_data.status,
        time_seconds: backend_data.time_seconds,
        objective_inaccessible,
    }
}

/// Generates transition table B.
///
/// This is the 8-variable occupancy / move table used by the current LCT model
/// builder.
pub fn generate_table_b() -> Vec<[i32; 8]> {
    let mut tuples = Vec::new();

    for mask in 0..(1 << 8) {
        let mut v = [0i32; 8];
        for (i, item) in v.iter_mut().enumerate() {
            *item = ((mask >> i) & 1) as i32;
        }
        if valid_tuple_b(&v) {
            tuples.push(v);
        }
    }

    tuples.sort_unstable();
    assert_eq!(tuples.len(), 7, "table B must have 7 tuples");
    tuples
}

/// Generates transition table A.
///
/// The current CPMP-LCT implementation only needs table B.
pub fn generate_table_a() -> Vec<[i32; 12]> {
    let mut tuples = Vec::new();

    for mask in 0..(1 << 12) {
        let mut v = [0i32; 12];
        for (i, item) in v.iter_mut().enumerate() {
            *item = ((mask >> i) & 1) as i32;
        }
        if valid_tuple_a(&v) {
            tuples.push(v);
        }
    }

    tuples.sort_unstable();
    assert_eq!(tuples.len(), 16, "table A must have 16 tuples");
    tuples
}

fn validate_lct2_unique_priorities(bay: &Bay) {
    for p in 1..=bay.p {
        assert!(
            bay.group_count[p] == 1,
            "LCT2 requires unique priorities (group {} has multiplicity {})",
            p,
            bay.group_count[p]
        );
    }
}

fn build_model_data(bay: &Bay, crane: &CraneParams, tau: f64) -> ModelData {
    let s = bay.s;
    let h = bay.h;
    let p = bay.p;
    let c = bay.c;

    // The stage bound is the available time divided by the
    // fastest possible relocation, rounded up and clamped to at least one stage.
    let t_min = crane.t_min(h);
    let mut k = (tau / t_min).ceil() as usize;
    if k == 0 {
        k = 1;
    }

    let mut init_x = vec![vec![0usize; h + 1]; s + 1];
    for (si, row) in init_x.iter_mut().enumerate().skip(1) {
        for (ti, cell) in row.iter_mut().enumerate().skip(1) {
            *cell = bay.g(si - 1, ti - 1) as usize;
        }
    }

    let mut group_count = vec![0usize; p + 1];
    for (g, gc) in group_count.iter_mut().enumerate().take(p + 1).skip(1) {
        *gc = bay.group_count[g] as usize;
    }
    group_count[0] = s * h - c;

    let mut h0_0_s = vec![0.0f64; s + 1];
    for (si, item) in h0_0_s.iter_mut().enumerate().skip(1) {
        *item = crane.r_u(0, si);
    }

    let mut h0_rs = vec![vec![0.0f64; s + 1]; s + 1];
    let mut h1_sr = vec![vec![0.0f64; s + 1]; s + 1];
    for r in 1..=s {
        for st in 1..=s {
            h0_rs[r][st] = crane.r_u(r, st);
            h1_sr[r][st] = crane.r_l(r, st);
        }
    }

    let mut v_load = vec![0.0f64; h + 1];
    let mut v_unload = vec![0.0f64; h + 1];
    for t in 1..=h {
        v_load[t] = crane.v_l(t);
        v_unload[t] = crane.v_u(t);
    }

    ModelData {
        s,
        h,
        p,
        c,
        k,
        init_x,
        group_count,
        table_b: generate_table_b(),
        h0_0_s,
        h0_rs,
        h1_sr,
        v_load,
        v_unload,
    }
}

fn solve_with_cpoptimizer(
    data: &ModelData,
    cfg: &CplctConfig,
    warm_start: Option<&WarmStartData>,
) -> BackendSolveData {
    let out_dir = Path::new("tmp-agent/cplct");
    fs::create_dir_all(out_dir).expect("cannot create tmp-agent/cplct");

    let stamp = unique_stamp();
    let model_path = out_dir.join(format!(
        "lct_{:?}_{:?}_{}.cpo",
        cfg.model, cfg.backend, stamp
    ));

    let cpo = generate_cpo_model(data, cfg, warm_start);
    fs::write(&model_path, cpo).expect("cannot write CPO model file");

    let workers = std::env::var("CP_OPTIMIZER_WORKERS")
        .ok()
        .and_then(|v| v.parse::<usize>().ok())
        .filter(|&v| v >= 1)
        .or_else(|| {
            std::env::var("SLURM_CPUS_PER_TASK")
                .ok()
                .and_then(|v| v.parse::<usize>().ok())
                .filter(|&v| v >= 1)
        })
        .unwrap_or(1);

    let args = vec![
        "-c".to_string(),
        format!("read {}", model_path.to_string_lossy()),
        format!("set Workers {}", workers),
        format!("set timelimit {}", cfg.time_limit_sec),
        "optimize".to_string(),
        "display solution".to_string(),
        "quit".to_string(),
    ];

    let output = Command::new("cpoptimizer")
        .args(&args)
        .output()
        .expect("failed to run cpoptimizer (is it installed and on PATH?)");

    let mut log = String::from_utf8_lossy(&output.stdout).to_string();
    if !output.stderr.is_empty() {
        log.push('\n');
        log.push_str(&String::from_utf8_lossy(&output.stderr));
    }

    if !output.status.success() {
        panic!("cpoptimizer failed:\n{}", log);
    }

    let optimal = log.contains("Optimal solution found");
    let status = if optimal {
        "OPTIMAL".to_string()
    } else if log.contains("Infeasible problem") {
        "INFEASIBLE".to_string()
    } else if log.contains("No solution found") {
        "NO_SOLUTION".to_string()
    } else if log.contains("solution found") {
        "FEASIBLE".to_string()
    } else {
        "UNKNOWN".to_string()
    };

    let objective_inaccessible = parse_cpoptimizer_objective(&log);
    let time_seconds = parse_cpoptimizer_time(&log);
    let y_ones = parse_cpoptimizer_binary_ones(&log, 'y');
    let z_ones = parse_cpoptimizer_binary_ones(&log, 'z');

    BackendSolveData {
        status,
        optimal,
        time_seconds,
        objective_inaccessible,
        y_ones,
        z_ones,
    }
}

fn solve_with_cpsat(
    data: &ModelData,
    cfg: &CplctConfig,
    warm_start: Option<&WarmStartData>,
) -> BackendSolveData {
    let check = Command::new("python3")
        .args(["-c", "import ortools"])
        .output()
        .expect("failed to run python3");
    assert!(
        check.status.success(),
        "OR-Tools python package not found. Install with: pip install ortools"
    );

    let out_dir = Path::new("tmp-agent/cplct");
    fs::create_dir_all(out_dir).expect("cannot create tmp-agent/cplct");

    let stamp = unique_stamp();
    let input_path = out_dir.join(format!("lct_cpsat_input_{}.json", stamp));
    let output_path = out_dir.join(format!("lct_cpsat_output_{}.json", stamp));

    let input = CpsatInput {
        model: match cfg.model {
            CplctModel::Lct1 => "lct1".to_string(),
            CplctModel::Lct2 => "lct2".to_string(),
        },
        s: data.s,
        h: data.h,
        p: data.p,
        c: data.c,
        k: data.k,
        tau: cfg.tau,
        time_limit_sec: cfg.time_limit_sec,
        init_x: data.init_x.clone(),
        group_count: data.group_count.clone(),
        table_b: data.table_b.clone(),
        h0_0_s: data.h0_0_s.clone(),
        h0_rs: data.h0_rs.clone(),
        h1_sr: data.h1_sr.clone(),
        v_load: data.v_load.clone(),
        v_unload: data.v_unload.clone(),
        warm_start_y_ones: warm_start.map_or_else(Vec::new, |w| w.y_ones.clone()),
        warm_start_z_ones: warm_start.map_or_else(Vec::new, |w| w.z_ones.clone()),
    };

    let json = serde_json::to_string_pretty(&input).expect("failed to serialize CP-SAT input");
    fs::write(&input_path, json).expect("cannot write CP-SAT input file");

    let output = Command::new("python3")
        .arg("scripts/solve_lct_cpsat.py")
        .arg(&input_path)
        .arg(&output_path)
        .output()
        .expect("failed to run scripts/solve_lct_cpsat.py");

    if !output.status.success() {
        let stdout = String::from_utf8_lossy(&output.stdout);
        let stderr = String::from_utf8_lossy(&output.stderr);
        panic!(
            "CP-SAT wrapper failed.\nstdout:\n{}\nstderr:\n{}",
            stdout, stderr
        );
    }

    let out_json = fs::read_to_string(&output_path).expect("cannot read CP-SAT output file");
    let parsed: CpsatOutput =
        serde_json::from_str(&out_json).expect("failed to parse CP-SAT output JSON");

    BackendSolveData {
        status: parsed.status,
        optimal: parsed.optimal,
        time_seconds: parsed.time_seconds,
        objective_inaccessible: parsed.objective,
        y_ones: parsed.y_ones,
        z_ones: parsed.z_ones,
    }
}

fn generate_cpo_model(
    data: &ModelData,
    cfg: &CplctConfig,
    warm_start: Option<&WarmStartData>,
) -> String {
    let mut cpo = String::new();

    writeln!(cpo, "// CPMP-LCT {:?} {:?}", cfg.model, cfg.backend).unwrap();
    writeln!(
        cpo,
        "// S={} H={} P={} C={} K={} tau={:.6}",
        data.s, data.h, data.p, data.c, data.k, cfg.tau
    )
    .unwrap();

    // Decision variables: bay state, move origins / destinations, horizontal
    // movement indicators, and the model-specific final-stage objective variables.
    for s in 1..=data.s {
        for t in 1..=data.h {
            for k in 0..=data.k {
                writeln!(cpo, "x_{}_{}_{} = intVar(0..{});", s, t, k, data.p).unwrap();
                writeln!(cpo, "d_{}_{}_{} = intVar(0..1);", s, t, k).unwrap();
            }
        }
    }

    for s in 1..=data.s {
        for t in 1..=data.h {
            for k in 1..=data.k {
                writeln!(cpo, "y_{}_{}_{} = intVar(0..1);", s, t, k).unwrap();
                writeln!(cpo, "z_{}_{}_{} = intVar(0..1);", s, t, k).unwrap();
            }
        }
    }

    for s in 1..=data.s {
        for r in 1..=data.s {
            if s == r {
                continue;
            }
            for k in 1..=data.k {
                writeln!(cpo, "g_{}_{}_{} = intVar(0..1);", s, r, k).unwrap();
            }
            for k in 2..=data.k {
                writeln!(cpo, "u_{}_{}_{} = intVar(0..1);", s, r, k).unwrap();
            }
        }
    }

    match cfg.model {
        CplctModel::Lct1 => {
            for s in 1..=data.s {
                for t in 1..=data.h {
                    writeln!(cpo, "q_{}_{} = intVar(0..1);", s, t).unwrap();
                }
            }
        }
        CplctModel::Lct2 => {
            for p in 1..=data.p {
                writeln!(cpo, "qp_{} = intVar(0..1);", p).unwrap();
            }
        }
    }

    writeln!(cpo).unwrap();
    append_cpoptimizer_starting_point(&mut cpo, warm_start);
    writeln!(cpo).unwrap();

    // Objective: LCT1 counts inaccessible containers, LCT2 counts inaccessible priority groups.
    match cfg.model {
        CplctModel::Lct1 => {
            let mut vars = Vec::new();
            for s in 1..=data.s {
                for t in 1..=data.h {
                    vars.push(format!("q_{}_{}", s, t));
                }
            }
            writeln!(cpo, "minimize(sum([{}]));", vars.join(", ")).unwrap();
        }
        CplctModel::Lct2 => {
            let vars: Vec<String> = (1..=data.p).map(|p| format!("qp_{}", p)).collect();
            writeln!(cpo, "minimize(sum([{}]));", vars.join(", ")).unwrap();
        }
    }

    writeln!(cpo).unwrap();

    // Initial state and container-multiplicity constraints.
    for s in 1..=data.s {
        for t in 1..=data.h {
            let x0 = data.init_x[s][t];
            let d0 = if x0 > 0 { 1 } else { 0 };
            writeln!(cpo, "x_{}_{}_0 == {};", s, t, x0).unwrap();
            writeln!(cpo, "d_{}_{}_0 == {};", s, t, d0).unwrap();
        }
    }

    for k in 1..=data.k {
        let stage_vars: Vec<String> = (1..=data.s)
            .flat_map(|s| (1..=data.h).map(move |t| format!("x_{}_{}_{}", s, t, k)))
            .collect();
        for p in 0..=data.p {
            writeln!(
                cpo,
                "count([{}], {}) == {};",
                stage_vars.join(", "),
                p,
                data.group_count[p]
            )
            .unwrap();
        }
    }

    for k in 1..=data.k {
        for s in 1..=data.s {
            for t in 1..=data.h {
                writeln!(cpo, "(x_{}_{}_{} > 0) == d_{}_{}_{};", s, t, k, s, t, k).unwrap();
                writeln!(
                    cpo,
                    "((x_{}_{}_{} == x_{}_{}_{}) == (d_{}_{}_{} == d_{}_{}_{}));",
                    s,
                    t,
                    k,
                    s,
                    t,
                    k - 1,
                    s,
                    t,
                    k,
                    s,
                    t,
                    k - 1
                )
                .unwrap();
            }
        }
    }

    if data.k >= 2 {
        for k in 1..data.k {
            for s in 1..=data.s {
                for t in 1..=data.h {
                    writeln!(cpo, "y_{}_{}_{} <= d_{}_{}_{};", s, t, k, s, t, k + 1).unwrap();
                }
            }
        }
    }

    for k in 1..=data.k {
        let y_vars: Vec<String> = (1..=data.s)
            .flat_map(|s| (1..=data.h).map(move |t| format!("y_{}_{}_{}", s, t, k)))
            .collect();
        let z_vars: Vec<String> = (1..=data.s)
            .flat_map(|s| (1..=data.h).map(move |t| format!("z_{}_{}_{}", s, t, k)))
            .collect();
        writeln!(cpo, "sum([{}]) <= 1;", y_vars.join(", ")).unwrap();
        writeln!(cpo, "sum([{}]) <= 1;", z_vars.join(", ")).unwrap();
    }

    if data.k >= 2 {
        for k in 1..data.k {
            let z_curr: Vec<String> = (1..=data.s)
                .flat_map(|s| (1..=data.h).map(move |t| format!("z_{}_{}_{}", s, t, k)))
                .collect();
            let z_next: Vec<String> = (1..=data.s)
                .flat_map(|s| (1..=data.h).map(move |t| format!("z_{}_{}_{}", s, t, k + 1)))
                .collect();
            writeln!(
                cpo,
                "sum([{}]) <= sum([{}]);",
                z_next.join(", "),
                z_curr.join(", ")
            )
            .unwrap();
        }
    }

    // Stage-to-stage occupancy transitions are enforced with the paper's table B.
    let table_b_literal = table8_literal(&data.table_b);
    if data.h >= 2 {
        for k in 1..=data.k {
            for s in 1..=data.s {
                for t in 1..data.h {
                    let vars = vec![
                        format!("d_{}_{}_{}", s, t, k - 1),
                        format!("d_{}_{}_{}", s, t + 1, k - 1),
                        format!("d_{}_{}_{}", s, t, k),
                        format!("d_{}_{}_{}", s, t + 1, k),
                        format!("z_{}_{}_{}", s, t, k),
                        format!("z_{}_{}_{}", s, t + 1, k),
                        format!("y_{}_{}_{}", s, t, k),
                        format!("y_{}_{}_{}", s, t + 1, k),
                    ];
                    writeln!(
                        cpo,
                        "allowedAssignments([{}], {});",
                        vars.join(", "),
                        table_b_literal
                    )
                    .unwrap();
                }
            }
        }
    }

    // Link per-slot move indicators to stack-to-stack horizontal movement variables.
    for k in 1..=data.k {
        for s in 1..=data.s {
            for r in 1..=data.s {
                if s == r {
                    continue;
                }
                let mut vars = Vec::new();
                for t in 1..=data.h {
                    vars.push(format!("z_{}_{}_{}", s, t, k));
                    vars.push(format!("y_{}_{}_{}", r, t, k));
                }
                writeln!(
                    cpo,
                    "g_{}_{}_{} == (sum([{}]) == 2);",
                    s,
                    r,
                    k,
                    vars.join(", ")
                )
                .unwrap();

                if k >= 2 {
                    let mut vars_u = Vec::new();
                    for t in 1..=data.h {
                        vars_u.push(format!("y_{}_{}_{}", s, t, k - 1));
                        vars_u.push(format!("z_{}_{}_{}", r, t, k));
                    }
                    writeln!(
                        cpo,
                        "u_{}_{}_{} == (sum([{}]) == 2);",
                        s,
                        r,
                        k,
                        vars_u.join(", ")
                    )
                    .unwrap();
                }
            }
        }
    }

    let mut time_terms = Vec::new();

    for s in 1..=data.s {
        for r in 1..=data.s {
            if s == r {
                continue;
            }
            time_terms.push(format!("{:.6} * g_{}_{}_1", data.h0_0_s[s], s, r));
        }
    }

    for r in 1..=data.s {
        for s in 1..=data.s {
            if r == s {
                continue;
            }
            for k in 2..=data.k {
                time_terms.push(format!("{:.6} * u_{}_{}_{}", data.h0_rs[r][s], r, s, k));
            }
        }
    }

    for s in 1..=data.s {
        for r in 1..=data.s {
            if s == r {
                continue;
            }
            for k in 1..=data.k {
                time_terms.push(format!("{:.6} * g_{}_{}_{}", data.h1_sr[s][r], s, r, k));
            }
        }
    }

    for s in 1..=data.s {
        for t in 1..=data.h {
            for k in 1..=data.k {
                time_terms.push(format!("{:.6} * z_{}_{}_{}", data.v_load[t], s, t, k));
                time_terms.push(format!("{:.6} * y_{}_{}_{}", data.v_unload[t], s, t, k));
            }
        }
    }

    if time_terms.is_empty() {
        writeln!(cpo, "0 <= {:.6};", cfg.tau).unwrap();
    } else {
        writeln!(cpo, "{} <= {:.6};", time_terms.join(" + "), cfg.tau).unwrap();
    }

    let kf = data.k;

    // Final-stage accessibility objective constraints differ between LCT1 and LCT2.
    match cfg.model {
        CplctModel::Lct1 => {
            if data.h >= 2 {
                for s in 1..=data.s {
                    for t in 1..data.h {
                        let terms: Vec<String> = ((t + 1)..=data.h)
                            .map(|j| format!("(x_{}_{}_{} > x_{}_{}_{})", s, j, kf, s, t, kf))
                            .collect();
                        writeln!(
                            cpo,
                            "sum([{}]) <= {} * q_{}_{};",
                            terms.join(", "),
                            data.h,
                            s,
                            t
                        )
                        .unwrap();
                    }
                }
            }

            for s in 1..=data.s {
                for t in 1..=data.h {
                    let mut terms = Vec::new();
                    for i in 1..=data.s {
                        for j in 1..=data.h {
                            terms.push(format!(
                                "((x_{}_{}_{} > x_{}_{}_{}) * q_{}_{})",
                                s, t, kf, i, j, kf, i, j
                            ));
                        }
                    }
                    writeln!(
                        cpo,
                        "sum([{}]) <= {} * q_{}_{};",
                        terms.join(", "),
                        data.c,
                        s,
                        t
                    )
                    .unwrap();
                }
            }
        }
        CplctModel::Lct2 => {
            for p in 1..=data.p {
                let mut terms = Vec::new();
                for s in 1..=data.s {
                    for t in 1..=data.h {
                        if t >= data.h {
                            continue;
                        }
                        let blockers: Vec<String> = ((t + 1)..=data.h)
                            .map(|j| format!("(x_{}_{}_{} > x_{}_{}_{})", s, j, kf, s, t, kf))
                            .collect();
                        if blockers.is_empty() {
                            continue;
                        }
                        let inner = format!("sum([{}])", blockers.join(", "));
                        terms.push(format!("((x_{}_{}_{} == {}) * ({}))", s, t, kf, p, inner));
                    }
                }
                if terms.is_empty() {
                    writeln!(cpo, "0 <= {} * qp_{};", data.c, p).unwrap();
                } else {
                    writeln!(cpo, "sum([{}]) <= {} * qp_{};", terms.join(", "), data.c, p).unwrap();
                }
            }

            if data.p >= 2 {
                for p in 1..data.p {
                    writeln!(cpo, "qp_{} <= qp_{};", p, p + 1).unwrap();
                }
            }
        }
    }

    cpo
}

fn run_warm_start_heuristic(
    bay: &Bay,
    crane: &CraneParams,
    tau: f64,
    heuristic: WarmStartHeuristic,
) -> Solution {
    match heuristic {
        WarmStartHeuristic::Glct => {
            let mut bay_copy = *bay;
            glct::glct(&mut bay_copy, crane, tau)
        }
        WarmStartHeuristic::Hlct => hlct::hlct(bay, crane, tau),
    }
}

fn build_warm_start_data(
    bay: &Bay,
    crane: &CraneParams,
    tau: f64,
    heuristic: WarmStartHeuristic,
    k_limit: usize,
) -> WarmStartData {
    let warm_sol = run_warm_start_heuristic(bay, crane, tau, heuristic);
    assert!(
        warm_sol.len() <= k_limit,
        "warm-start solution has {} moves but model allows at most {} stages",
        warm_sol.len(),
        k_limit
    );

    let mut y_ones = Vec::with_capacity(warm_sol.len());
    let mut z_ones = Vec::with_capacity(warm_sol.len());

    // Warm starts are passed as stage / stack / tier triples matching the emitted model.
    for (idx, dm) in warm_sol.moves.iter().enumerate() {
        let stage = idx + 1;
        z_ones.push([stage, dm.src + 1, dm.src_tier + 1]);
        y_ones.push([stage, dm.dst + 1, dm.dst_tier + 1]);
    }

    WarmStartData { y_ones, z_ones }
}

fn append_cpoptimizer_starting_point(cpo: &mut String, warm_start: Option<&WarmStartData>) {
    let Some(warm_start) = warm_start else {
        return;
    };
    if warm_start.y_ones.is_empty() && warm_start.z_ones.is_empty() {
        return;
    }

    writeln!(cpo, "startingPoint {{").unwrap();
    for &[k, s, t] in &warm_start.y_ones {
        writeln!(cpo, "  y_{}_{}_{} = 1;", s, t, k).unwrap();
    }
    for &[k, s, t] in &warm_start.z_ones {
        writeln!(cpo, "  z_{}_{}_{} = 1;", s, t, k).unwrap();
    }
    writeln!(cpo, "}}").unwrap();
}

fn reconstruct_solution(
    bay: &Bay,
    crane: &CraneParams,
    k: usize,
    y_ones: &[[usize; 3]],
    z_ones: &[[usize; 3]],
) -> Solution {
    let mut src_by_stage: HashMap<usize, usize> = HashMap::new();
    let mut dst_by_stage: HashMap<usize, usize> = HashMap::new();

    for &[kk, s, _t] in z_ones {
        match src_by_stage.get(&kk) {
            None => {
                src_by_stage.insert(kk, s);
            }
            Some(&prev) => {
                assert!(
                    prev == s,
                    "invalid solver output: multiple source stacks in stage {}",
                    kk
                );
            }
        }
    }

    for &[kk, s, _t] in y_ones {
        match dst_by_stage.get(&kk) {
            None => {
                dst_by_stage.insert(kk, s);
            }
            Some(&prev) => {
                assert!(
                    prev == s,
                    "invalid solver output: multiple destination stacks in stage {}",
                    kk
                );
            }
        }
    }

    let mut bay_state = *bay;
    let mut sol = Solution::new();

    // Solver output is sparse, so reconstruction first compacts it into one
    // source stack and one destination stack per used stage.
    for stage in 1..=k {
        let src = src_by_stage.get(&stage).copied();
        let dst = dst_by_stage.get(&stage).copied();

        match (src, dst) {
            (None, None) => continue,
            (Some(s), Some(d)) => {
                assert!(
                    s != d,
                    "invalid solver output: source == destination at stage {}",
                    stage
                );
                let m = Move {
                    src: s - 1,
                    dst: d - 1,
                };
                let dm = bay_state.apply_move(m);
                let prev_dst = sol.last_dst_1based();
                let dt = crane.move_time(s, dm.src_tier + 1, d, dm.dst_tier + 1, prev_dst);
                sol.push(dm, dt);
            }
            _ => {
                panic!(
                    "invalid solver output: missing source or destination for stage {}",
                    stage
                );
            }
        }
    }

    sol.acc = bay_state.acc();
    sol
}

fn parse_cpoptimizer_binary_ones(log: &str, prefix: char) -> Vec<[usize; 3]> {
    let mut out = Vec::new();

    for line in log.lines() {
        let ln = line.trim();
        if !ln.starts_with(prefix) {
            continue;
        }
        if !ln.contains(" = intVar(") {
            continue;
        }

        let mut parts = ln.splitn(2, '=');
        let name = parts.next().unwrap().trim();
        let rhs = parts.next().unwrap_or("").trim();

        let start = rhs.find("intVar(").map(|i| i + 7).unwrap_or(0);
        let end = rhs[start..]
            .find(')')
            .map(|i| start + i)
            .unwrap_or(rhs.len());
        let val = rhs[start..end].trim();

        if val.contains("..") {
            panic!(
                "cpoptimizer solution has non-fixed domain for {}: intVar({})",
                name, val
            );
        }

        if val != "1" {
            continue;
        }

        let toks: Vec<&str> = name.split('_').collect();
        assert!(
            toks.len() == 4,
            "unexpected CP variable name format: {}",
            name
        );
        let s: usize = toks[1].parse().expect("invalid stack index in var name");
        let t: usize = toks[2].parse().expect("invalid tier index in var name");
        let k: usize = toks[3].parse().expect("invalid stage index in var name");
        out.push([k, s, t]);
    }

    out.sort_unstable();
    out
}

fn parse_cpoptimizer_objective(log: &str) -> Option<usize> {
    for line in log.lines() {
        let ln = line.trim();
        if ln.starts_with("! Best objective") {
            if let Some(colon) = ln.find(':') {
                let tail = ln[colon + 1..].trim();
                let token = tail.split_whitespace().next().unwrap_or("");
                if let Ok(v) = token.parse::<f64>() {
                    return Some(v.round() as usize);
                }
            }
        }
    }

    for line in log.lines() {
        let ln = line.trim();
        if let Some(pos) = ln.find("objective") {
            let tail = ln[pos + "objective".len()..].trim();
            let token = tail
                .trim_start_matches('=')
                .trim_end_matches('.')
                .split_whitespace()
                .next()
                .unwrap_or("");
            if let Ok(v) = token.parse::<f64>() {
                return Some(v.round() as usize);
            }
        }
    }

    None
}

fn parse_cpoptimizer_time(log: &str) -> f64 {
    for line in log.lines() {
        let ln = line.trim();
        if !ln.contains("Time spent in solve") {
            continue;
        }
        if let Some(colon) = ln.find(':') {
            let tail = ln[colon + 1..].trim();
            let token = tail.split_whitespace().next().unwrap_or("");
            let token = token.trim_end_matches('s');
            if let Ok(v) = token.parse::<f64>() {
                return v;
            }
        }
    }
    0.0
}

fn unique_stamp() -> String {
    let now = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .expect("clock before unix epoch");
    let nanos = now.as_nanos();
    let pid = std::process::id() as u64;
    let counter = UNIQUE_STAMP_COUNTER.fetch_add(1, Ordering::Relaxed);
    let random_tail: u64 = rand::random();

    format!("{:x}_{:x}_{:x}_{:x}", nanos, pid, counter, random_tail)
}

fn table8_literal(table: &[[i32; 8]]) -> String {
    let mut rows = Vec::with_capacity(table.len());
    for row in table {
        rows.push(format!(
            "[{}, {}, {}, {}, {}, {}, {}, {}]",
            row[0], row[1], row[2], row[3], row[4], row[5], row[6], row[7]
        ));
    }
    format!("[{}]", rows.join(", "))
}

fn valid_tuple_b(v: &[i32; 8]) -> bool {
    let d0 = v[0];
    let d1 = v[1];
    let d0n = v[2];
    let d1n = v[3];
    let z0 = v[4];
    let z1 = v[5];
    let y0 = v[6];
    let y1 = v[7];

    if d1 > d0 || d1n > d0n {
        return false;
    }

    if d0n != d0 + y0 - z0 || d1n != d1 + y1 - z1 {
        return false;
    }

    if z0 > d0 || z1 > d1 {
        return false;
    }

    if y0 > 1 - d0 || y1 > 1 - d1 {
        return false;
    }

    if (z0 == 1 && y0 == 1) || (z1 == 1 && y1 == 1) {
        return false;
    }

    if z0 == 1 && d1 == 1 {
        return false;
    }

    if y0 == 1 && d1n == 1 {
        return false;
    }

    if y1 == 1 && d0n == 0 {
        return false;
    }

    if z0 + z1 > 1 || y0 + y1 > 1 {
        return false;
    }

    true
}

fn valid_tuple_a(v: &[i32; 12]) -> bool {
    let d0 = v[0];
    let d1 = v[1];
    let d0n = v[2];
    let d1n = v[3];
    let z0 = v[4];
    let z1 = v[5];
    let y0 = v[6];
    let y1 = v[7];
    let w0 = v[8];
    let w1 = v[9];
    let w0n = v[10];
    let w1n = v[11];

    if w0 > d0 || w1 > d1 || w0n > d0n || w1n > d1n {
        return false;
    }

    if d1 > d0 || d1n > d0n {
        return false;
    }

    if d0 == 0 && (d1 != 0 || w0 != 0 || w1 != 0) {
        return false;
    }
    if d0 == 1 && w0 == 1 && d1 == 1 && w1 != 1 {
        return false;
    }

    if d0n == 0 && (d1n != 0 || w0n != 0 || w1n != 0) {
        return false;
    }
    if d0n == 1 && w0n == 1 && d1n == 1 && w1n != 1 {
        return false;
    }

    if d0n != d0 + y0 - z0 || d1n != d1 + y1 - z1 {
        return false;
    }

    if z0 > d0 || z1 > d1 {
        return false;
    }
    if y0 > 1 - d0 || y1 > 1 - d1 {
        return false;
    }
    if (z0 == 1 && y0 == 1) || (z1 == 1 && y1 == 1) {
        return false;
    }

    if z0 == 1 && d1 == 1 {
        return false;
    }
    if y0 == 1 && d1n == 1 {
        return false;
    }
    if y1 == 1 && d0n == 0 {
        return false;
    }

    if z0 + z1 > 1 || y0 + y1 > 1 {
        return false;
    }

    if d0 == 1 && d0n == 1 && z0 == 0 && y0 == 0 && w0n != w0 {
        return false;
    }
    if d1 == 1 && d1n == 1 && z1 == 0 && y1 == 0 && w1n != w1 {
        return false;
    }

    if z1 == 1 && d0 == 1 && d0n == 1 && w0n != w0 {
        return false;
    }

    if y1 == 1 && w0n == 1 && w1n != 1 {
        return false;
    }

    true
}
