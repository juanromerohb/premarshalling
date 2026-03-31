//! Structured result writers for benchmark and GUI runs.

use std::fmt::Write as FmtWrite;
use std::fs;
use std::io::Write;
use std::path::PathBuf;
use std::process::Command;

use crate::bay::Bay;
use crate::benchmark::{Algorithm, BenchResult};
use crate::display;
use crate::instance;

/// Machine metadata recorded alongside structured result files.
pub struct MachineInfo {
    pub cpu_model: String,
    pub ram_total_mb: u64,
    pub os_info: String,
    pub rust_version: String,
}

impl MachineInfo {
    /// Collects machine information from the local environment.
    pub fn collect() -> Self {
        MachineInfo {
            cpu_model: read_cpu_model(),
            ram_total_mb: read_ram_mb(),
            os_info: read_os_info(),
            rust_version: read_rust_version(),
        }
    }
}

fn read_cpu_model() -> String {
    fs::read_to_string("/proc/cpuinfo")
        .ok()
        .and_then(|contents| {
            contents
                .lines()
                .find(|l| l.starts_with("model name"))
                .map(|l| l.split_once(':').map_or(l, |(_, v)| v).trim().to_string())
        })
        .unwrap_or_else(|| "unknown".to_string())
}

fn read_ram_mb() -> u64 {
    fs::read_to_string("/proc/meminfo")
        .ok()
        .and_then(|contents| {
            contents
                .lines()
                .find(|l| l.starts_with("MemTotal"))
                .map(|l| {
                    let kb: u64 = l
                        .split_whitespace()
                        .nth(1)
                        .and_then(|v| v.parse().ok())
                        .unwrap_or(0);
                    kb / 1024
                })
        })
        .unwrap_or(0)
}

fn read_os_info() -> String {
    Command::new("uname")
        .args(["-srm"])
        .output()
        .ok()
        .map(|o| String::from_utf8_lossy(&o.stdout).trim().to_string())
        .unwrap_or_else(|| "unknown".to_string())
}

fn read_rust_version() -> String {
    Command::new("rustc")
        .arg("--version")
        .output()
        .ok()
        .map(|o| String::from_utf8_lossy(&o.stdout).trim().to_string())
        .unwrap_or_else(|| "unknown".to_string())
}

fn current_date_ymd() -> String {
    Command::new("date")
        .arg("+%Y%m%d")
        .output()
        .ok()
        .map(|o| String::from_utf8_lossy(&o.stdout).trim().to_string())
        .unwrap_or_else(|| "00000000".to_string())
}

fn current_datetime() -> String {
    Command::new("date")
        .arg("+%Y-%m-%d %H:%M:%S")
        .output()
        .ok()
        .map(|o| String::from_utf8_lossy(&o.stdout).trim().to_string())
        .unwrap_or_else(|| "unknown".to_string())
}

use crate::bblct::{BblctConfig, HeuristicChoice};

impl Algorithm {
    /// Returns the option tag used in structured result folder names.
    pub fn options_tag(&self, tau: f64) -> String {
        let mut parts = Vec::new();

        let tau_minutes = tau / 60.0;
        if (tau_minutes - tau_minutes.round()).abs() < 0.01 {
            parts.push(format!("tau{}m", tau_minutes.round() as u64));
        } else {
            parts.push(format!("tau{:.1}m", tau_minutes));
        }

        match self {
            Algorithm::Bblct { config } => {
                let def = BblctConfig::default();
                if config.initial_heuristic != def.initial_heuristic {
                    parts.push(format!("init-{}", heuristic_tag(&config.initial_heuristic)));
                }
                if config.use_dominance != def.use_dominance {
                    parts.push("no-dom".into());
                }
                if config.use_lb_pruning != def.use_lb_pruning {
                    parts.push("no-lbp".into());
                }
                if config.use_lb_sorting != def.use_lb_sorting {
                    parts.push("no-lbs".into());
                }
                if config.node_heuristic != def.node_heuristic {
                    parts.push(format!("nodeh-{}", heuristic_tag(&config.node_heuristic)));
                }
                if config.node_heuristic_every != def.node_heuristic_every {
                    parts.push(format!("nhe{}", config.node_heuristic_every));
                }
            }
            Algorithm::Lct1Cpo { time_limit_sec }
            | Algorithm::Lct2Cpo { time_limit_sec }
            | Algorithm::Lct1Cpsat { time_limit_sec }
            | Algorithm::Lct2Cpsat { time_limit_sec } => {
                if (*time_limit_sec - 3600.0).abs() > 0.1 {
                    parts.push(format!("tlim{}", *time_limit_sec as u64));
                }
            }
            Algorithm::Iplct { time_limit_sec } => {
                if (*time_limit_sec - 3600.0).abs() > 0.1 {
                    parts.push(format!("tlim{}", *time_limit_sec as u64));
                }
            }
            _ => {}
        }

        if parts.is_empty() {
            "def".to_string()
        } else {
            parts.join("_")
        }
    }

    /// Returns a human-readable description of the algorithm options.
    pub fn options_description(&self) -> String {
        match self {
            Algorithm::Bblct { config } => {
                format!(
                    "initial_heuristic={}, dominance={}, lb_pruning={}, lb_sorting={}, node_heuristic={}, node_heuristic_every={}",
                    heuristic_tag(&config.initial_heuristic),
                    config.use_dominance,
                    config.use_lb_pruning,
                    config.use_lb_sorting,
                    heuristic_tag(&config.node_heuristic),
                    config.node_heuristic_every,
                )
            }
            Algorithm::Lct1Cpo { time_limit_sec }
            | Algorithm::Lct2Cpo { time_limit_sec }
            | Algorithm::Lct1Cpsat { time_limit_sec }
            | Algorithm::Lct2Cpsat { time_limit_sec } => {
                format!("time_limit={}s", time_limit_sec)
            }
            Algorithm::Iplct { time_limit_sec } => {
                format!("time_limit={}s", time_limit_sec)
            }
            _ => "none".to_string(),
        }
    }
}

fn heuristic_tag(h: &HeuristicChoice) -> &'static str {
    match h {
        HeuristicChoice::None => "none",
        HeuristicChoice::Glct => "glct",
        HeuristicChoice::Hlct => "hlct",
    }
}

/// Builds the per-algorithm result folder name for a given run.
pub fn build_folder_name(algo: &Algorithm, tau: f64, run_tag: Option<&str>) -> String {
    let tag = algo.options_tag(tau);
    let date = current_date_ymd();
    if let Some(run_tag) = run_tag {
        let run_tag = run_tag.trim();
        if !run_tag.is_empty() {
            return format!("{}_{}_{}", tag, date, run_tag);
        }
    }
    format!("{}_{}", tag, date)
}

fn resolve_results_base(instance_name: &str, instance_dir: &str) -> PathBuf {
    let stem = instance_name.strip_suffix(".txt").unwrap_or(instance_name);

    let dir_path = std::path::Path::new(instance_dir);
    let components: Vec<&str> = dir_path
        .components()
        .filter_map(|c| c.as_os_str().to_str())
        .collect();

    let set_name = find_set_name(&components);

    if stem.starts_with("BF") {
        let bf_prefix = stem.split('_').next().unwrap_or(stem);
        let inner_stem = stem
            .strip_prefix(bf_prefix)
            .and_then(|s| s.strip_prefix('_'))
            .unwrap_or(stem);
        let _ = inner_stem;
        return PathBuf::from("results")
            .join("BF")
            .join(bf_prefix)
            .join(stem);
    }

    if let Some((h, s)) = instance::try_parse_hs(instance_name) {
        let category = format!("{}h{}s{}", set_name, h, s);
        PathBuf::from("results")
            .join(set_name)
            .join(category)
            .join(stem)
    } else {
        PathBuf::from("results").join(set_name).join(stem)
    }
}

fn find_set_name(components: &[&str]) -> String {
    for &comp in components.iter().rev() {
        match comp {
            "CV" | "BZ" | "ZJY" | "EMM" | "BF" => return comp.to_string(),
            _ => {
                if comp.starts_with("BF") && comp[2..].chars().all(|c| c.is_ascii_digit()) {
                    return "BF".to_string();
                }
            }
        }
    }
    components.last().unwrap_or(&"unknown").to_string()
}

fn write_technical(
    path: &std::path::Path,
    machine: &MachineInfo,
    algo: &Algorithm,
    datetime: &str,
) {
    let mut f = fs::File::create(path).expect("cannot create technical.txt");

    writeln!(f, "=== Technical Details ===").unwrap();
    writeln!(f, "Date/Time:     {}", datetime).unwrap();
    writeln!(f, "CPU:           {}", machine.cpu_model).unwrap();
    writeln!(f, "RAM:           {} MB", machine.ram_total_mb).unwrap();
    writeln!(f, "OS:            {}", machine.os_info).unwrap();
    writeln!(f, "Rust version:  {}", machine.rust_version).unwrap();

    let solver = match algo {
        Algorithm::Lct1Cpo { .. } | Algorithm::Lct2Cpo { .. } => detect_solver_version("cpo"),
        Algorithm::Lct1Cpsat { .. } | Algorithm::Lct2Cpsat { .. } => detect_solver_version("cpsat"),
        Algorithm::Iplct { .. } => detect_solver_version("gurobi"),
        _ => "N/A".to_string(),
    };
    writeln!(f, "Solver:        {}", solver).unwrap();
}

fn detect_solver_version(backend: &str) -> String {
    match backend {
        "cpo" => Command::new("cpoptimizer")
            .arg("-version")
            .output()
            .ok()
            .map(|o| {
                let out = String::from_utf8_lossy(&o.stdout);
                let err = String::from_utf8_lossy(&o.stderr);
                let combined = format!("{}{}", out.trim(), err.trim());
                if combined.is_empty() {
                    "IBM CP Optimizer (version unknown)".to_string()
                } else {
                    combined
                        .lines()
                        .next()
                        .unwrap_or("IBM CP Optimizer")
                        .to_string()
                }
            })
            .unwrap_or_else(|| "IBM CP Optimizer (not found)".to_string()),
        "cpsat" => Command::new("python3")
            .args([
                "-c",
                "from ortools.sat.python import cp_model; import ortools; print(f'OR-Tools CP-SAT {ortools.__version__}')",
            ])
            .output()
            .ok()
            .map(|o| {
                let out = String::from_utf8_lossy(&o.stdout).trim().to_string();
                if out.is_empty() {
                    "OR-Tools CP-SAT (version unknown)".to_string()
                } else {
                    out
                }
            })
            .unwrap_or_else(|| "OR-Tools CP-SAT (not found)".to_string()),
        "gurobi" => Command::new("gurobi_cl")
            .arg("--version")
            .output()
            .ok()
            .map(|o| {
                let out = String::from_utf8_lossy(&o.stdout).trim().to_string();
                out.lines()
                    .next()
                    .unwrap_or("Gurobi (version unknown)")
                    .to_string()
            })
            .unwrap_or_else(|| "Gurobi (not found)".to_string()),
        _ => "N/A".to_string(),
    }
}

fn write_analysis(
    path: &std::path::Path,
    result: &BenchResult,
    algo: &Algorithm,
    initial_acc: usize,
) {
    let c = result.c;
    let initial_inacc = c - initial_acc;
    let final_acc = result.acc;
    let final_inacc = c - final_acc;

    let acc_increase_pct = if initial_acc == 0 {
        if final_acc > 0 {
            f64::INFINITY
        } else {
            0.0
        }
    } else {
        (final_acc as f64 - initial_acc as f64) / initial_acc as f64 * 100.0
    };

    let inacc_reduction_pct = if initial_inacc == 0 {
        0.0
    } else {
        (initial_inacc as f64 - final_inacc as f64) / initial_inacc as f64 * 100.0
    };

    let mut f = fs::File::create(path).expect("cannot create analysis.txt");

    writeln!(f, "=== Analysis ===").unwrap();
    writeln!(f, "Algorithm: {}", algo.name()).unwrap();
    writeln!(f, "Options:   {}", algo.options_description()).unwrap();
    writeln!(f).unwrap();
    writeln!(
        f,
        "Instance:  {}  (S={}, H={}, C={}, P={})",
        result.instance_name, result.s, result.h, result.c, result.p
    )
    .unwrap();
    writeln!(f, "Tau:       {:.1}s", result.tau).unwrap();
    let proved_optimal = match result.proved_optimal {
        Some(true) => "yes",
        Some(false) => "no",
        None => "N/A",
    };
    writeln!(f, "Proved optimal: {}", proved_optimal).unwrap();
    if let Some(stop_reason) = result.stop_reason.as_deref() {
        writeln!(f, "Stop reason:    {}", stop_reason).unwrap();
    }
    writeln!(f).unwrap();

    writeln!(f, "Initial:").unwrap();
    writeln!(
        f,
        "  Accessible:    {}/{}  ({:.1}%)",
        initial_acc,
        c,
        if c > 0 {
            100.0 * initial_acc as f64 / c as f64
        } else {
            100.0
        }
    )
    .unwrap();
    writeln!(
        f,
        "  Inaccessible:  {}/{}  ({:.1}%)",
        initial_inacc,
        c,
        if c > 0 {
            100.0 * initial_inacc as f64 / c as f64
        } else {
            0.0
        }
    )
    .unwrap();
    writeln!(f).unwrap();

    writeln!(f, "Final:").unwrap();
    writeln!(
        f,
        "  Accessible:    {}/{}  ({:.1}%)",
        final_acc,
        c,
        if c > 0 {
            100.0 * final_acc as f64 / c as f64
        } else {
            100.0
        }
    )
    .unwrap();
    writeln!(
        f,
        "  Inaccessible:  {}/{}  ({:.1}%)",
        final_inacc,
        c,
        if c > 0 {
            100.0 * final_inacc as f64 / c as f64
        } else {
            0.0
        }
    )
    .unwrap();
    writeln!(f).unwrap();

    writeln!(f, "Improvement:").unwrap();
    if acc_increase_pct.is_infinite() {
        writeln!(
            f,
            "  Accessible increase:     inf  (from {} to {})",
            initial_acc, final_acc
        )
        .unwrap();
    } else {
        writeln!(
            f,
            "  Accessible increase:     {:.1}%  (from {} to {})",
            acc_increase_pct, initial_acc, final_acc
        )
        .unwrap();
    }
    writeln!(
        f,
        "  Inaccessible reduction:  {:.1}%  (from {} to {})",
        inacc_reduction_pct, initial_inacc, final_inacc
    )
    .unwrap();
    writeln!(f).unwrap();

    writeln!(f, "Crane time:  {:.2}s", result.crane_time).unwrap();
    writeln!(f, "CPU time:    {:.1}ms", result.elapsed_ms).unwrap();
    writeln!(
        f,
        "Nodes:       {}",
        result
            .nodes_explored
            .map_or("N/A".to_string(), |n| n.to_string())
    )
    .unwrap();
    writeln!(f).unwrap();

    let mut sol_str = String::new();
    for (i, dm) in result.solution.moves.iter().enumerate() {
        if i > 0 {
            sol_str.push(' ');
        }
        write!(sol_str, "({},{})", dm.src + 1, dm.dst + 1).unwrap();
    }
    writeln!(f, "Solution: {}", sol_str).unwrap();
}

fn write_replay(path: &std::path::Path, bay: &Bay, result: &BenchResult) {
    let replay = display::replay_solution(bay, &result.solution, bay.p);
    fs::write(path, replay).expect("cannot create replay.txt");
}
struct Stats {
    mean: f64,
    std: f64,
    min: f64,
    max: f64,
}

impl Stats {
    fn compute(data: &[f64]) -> Self {
        let n = data.len() as f64;
        if data.is_empty() {
            return Stats {
                mean: 0.0,
                std: 0.0,
                min: 0.0,
                max: 0.0,
            };
        }
        let mean = data.iter().sum::<f64>() / n;
        let variance = data.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / n;
        let std = variance.sqrt();
        let min = data.iter().cloned().fold(f64::INFINITY, f64::min);
        let max = data.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        Stats {
            mean,
            std,
            min,
            max,
        }
    }
}

/// Saves the structured output bundle for one solved instance.
pub fn save_structured_results(
    result: &BenchResult,
    bay: &Bay,
    algo: &Algorithm,
    instance_dir: &str,
    machine: &MachineInfo,
    initial_acc: usize,
    folder_name: &str,
) {
    let base = resolve_results_base(&result.instance_name, instance_dir);
    let dir = base.join(folder_name);
    fs::create_dir_all(&dir).expect("cannot create results directory");

    let datetime = current_datetime();

    write_technical(&dir.join("technical.txt"), machine, algo, &datetime);
    write_analysis(&dir.join("analysis.txt"), result, algo, initial_acc);
    write_replay(&dir.join("replay.txt"), bay, result);
}

fn resolve_category_path(first_instance_name: &str, instance_dir: &str) -> PathBuf {
    let base = resolve_results_base(first_instance_name, instance_dir);
    base.parent()
        .map(PathBuf::from)
        .unwrap_or_else(|| PathBuf::from("results"))
}

/// Writes an aggregate analysis file for one algorithm/category combination.
pub fn write_category_analysis(
    results: &[&BenchResult],
    algo: &Algorithm,
    instance_dir: &str,
    folder_name: &str,
) {
    if results.is_empty() {
        return;
    }

    let category_path = resolve_category_path(&results[0].instance_name, instance_dir);
    let analysis_dir = category_path.join("analysis");
    fs::create_dir_all(&analysis_dir).expect("cannot create analysis directory");

    let path = analysis_dir.join(format!("{}.txt", folder_name));
    let mut f = fs::File::create(&path).expect("cannot create category analysis file");

    let n = results.len();
    let first = results[0];

    let cat_label = category_path
        .file_name()
        .and_then(|s| s.to_str())
        .unwrap_or("unknown");

    let initial_acc_pcts: Vec<f64> = results
        .iter()
        .map(|r| {
            if r.c > 0 {
                100.0 * r.initial_acc as f64 / r.c as f64
            } else {
                100.0
            }
        })
        .collect();
    let final_acc_pcts: Vec<f64> = results
        .iter()
        .map(|r| {
            if r.c > 0 {
                100.0 * r.acc as f64 / r.c as f64
            } else {
                100.0
            }
        })
        .collect();
    let acc_increase_pp: Vec<f64> = initial_acc_pcts
        .iter()
        .zip(final_acc_pcts.iter())
        .map(|(i, f)| f - i)
        .collect();
    let inacc_reduction_pcts: Vec<f64> = results
        .iter()
        .map(|r| {
            let initial_inacc = r.c - r.initial_acc;
            let final_inacc = r.c - r.acc;
            if initial_inacc == 0 {
                0.0
            } else {
                100.0 * (initial_inacc - final_inacc) as f64 / initial_inacc as f64
            }
        })
        .collect();
    let moves: Vec<f64> = results.iter().map(|r| r.moves as f64).collect();
    let crane_times: Vec<f64> = results.iter().map(|r| r.crane_time).collect();
    let cpu_times: Vec<f64> = results.iter().map(|r| r.elapsed_ms).collect();
    let nodes: Vec<f64> = results
        .iter()
        .filter_map(|r| r.nodes_explored.map(|n| n as f64))
        .collect();
    let fully_accessible = results.iter().filter(|r| r.acc == r.c).count();

    let s_init = Stats::compute(&initial_acc_pcts);
    let s_final = Stats::compute(&final_acc_pcts);
    let s_increase = Stats::compute(&acc_increase_pp);
    let s_inacc_red = Stats::compute(&inacc_reduction_pcts);
    let s_moves = Stats::compute(&moves);
    let s_crane = Stats::compute(&crane_times);
    let s_cpu = Stats::compute(&cpu_times);
    let has_nodes = !nodes.is_empty();
    let s_nodes = Stats::compute(&nodes);

    let date_str = if folder_name.len() >= 8 {
        let raw = &folder_name[folder_name.len() - 8..];
        format!("{}-{}-{}", &raw[0..4], &raw[4..6], &raw[6..8])
    } else {
        "unknown".to_string()
    };

    writeln!(f, "=== Category Analysis: {} ===", cat_label).unwrap();
    writeln!(f, "Algorithm: {}", algo.name()).unwrap();
    writeln!(f, "Options:   {}", algo.options_description()).unwrap();
    writeln!(f, "Tau:       {:.1}s", first.tau).unwrap();
    writeln!(f, "Date:      {}", date_str).unwrap();
    writeln!(f, "Instances: {}", n).unwrap();
    writeln!(f).unwrap();

    let c_values: Vec<usize> = results.iter().map(|r| r.c).collect();
    let p_values: Vec<usize> = results.iter().map(|r| r.p).collect();
    let c_min = c_values.iter().min().unwrap();
    let c_max = c_values.iter().max().unwrap();
    let p_min = p_values.iter().min().unwrap();
    let p_max = p_values.iter().max().unwrap();
    if c_min == c_max && p_min == p_max {
        writeln!(
            f,
            "Dimensions: S={}, H={}, C={}, P={}",
            first.s, first.h, c_min, p_min
        )
        .unwrap();
    } else {
        writeln!(
            f,
            "Dimensions: S={}, H={}, C={}..{}, P={}..{}",
            first.s, first.h, c_min, c_max, p_min, p_max
        )
        .unwrap();
    }
    writeln!(f).unwrap();

    // Stats table
    writeln!(
        f,
        "{:<24} {:>8} {:>8} {:>8} {:>8}",
        "", "Mean", "Std", "Min", "Max"
    )
    .unwrap();
    writeln!(
        f,
        "{:<24} {:>8} {:>8} {:>8} {:>8}",
        "", "----", "---", "---", "---"
    )
    .unwrap();

    write_stat_row(&mut f, "Initial acc (%)", &s_init, 1);
    write_stat_row(&mut f, "Final acc (%)", &s_final, 1);
    write_stat_row(&mut f, "Acc increase (pp)", &s_increase, 1);
    write_stat_row(&mut f, "Inacc reduction (%)", &s_inacc_red, 1);
    write_stat_row(&mut f, "Moves", &s_moves, 1);
    write_stat_row(&mut f, "Crane time (s)", &s_crane, 1);
    write_stat_row(&mut f, "CPU time (ms)", &s_cpu, 1);

    if has_nodes {
        write_stat_row(&mut f, "Nodes explored", &s_nodes, 0);
    } else {
        writeln!(
            f,
            "{:<24} {:>8} {:>8} {:>8} {:>8}",
            "Nodes explored", "--", "--", "--", "--"
        )
        .unwrap();
    }
    writeln!(f).unwrap();

    writeln!(
        f,
        "Fully accessible:  {}/{} ({:.1}%)",
        fully_accessible,
        n,
        100.0 * fully_accessible as f64 / n as f64
    )
    .unwrap();
    writeln!(f).unwrap();

    writeln!(f, "Per-instance details:").unwrap();
    if has_nodes {
        writeln!(
            f,
            "  {:<24} {:>7} {:>7} {:>6} {:>10} {:>9} {:>10}",
            "Instance", "Init%", "Final%", "Moves", "Crane(s)", "CPU(ms)", "Nodes"
        )
        .unwrap();
    } else {
        writeln!(
            f,
            "  {:<24} {:>7} {:>7} {:>6} {:>10} {:>9}",
            "Instance", "Init%", "Final%", "Moves", "Crane(s)", "CPU(ms)"
        )
        .unwrap();
    }

    for r in results {
        let init_pct = if r.c > 0 {
            100.0 * r.initial_acc as f64 / r.c as f64
        } else {
            100.0
        };
        let final_pct = if r.c > 0 {
            100.0 * r.acc as f64 / r.c as f64
        } else {
            100.0
        };
        if has_nodes {
            writeln!(
                f,
                "  {:<24} {:>7.1} {:>7.1} {:>6} {:>10.2} {:>9.1} {:>10}",
                r.instance_name,
                init_pct,
                final_pct,
                r.moves,
                r.crane_time,
                r.elapsed_ms,
                r.nodes_explored.map_or("--".to_string(), |n| n.to_string()),
            )
            .unwrap();
        } else {
            writeln!(
                f,
                "  {:<24} {:>7.1} {:>7.1} {:>6} {:>10.2} {:>9.1}",
                r.instance_name, init_pct, final_pct, r.moves, r.crane_time, r.elapsed_ms,
            )
            .unwrap();
        }
    }
}

fn write_stat_row(f: &mut fs::File, label: &str, s: &Stats, decimals: usize) {
    match decimals {
        0 => writeln!(
            f,
            "{:<24} {:>8.0} {:>8.0} {:>8.0} {:>8.0}",
            label, s.mean, s.std, s.min, s.max
        )
        .unwrap(),
        _ => writeln!(
            f,
            "{:<24} {:>8.1} {:>8.1} {:>8.1} {:>8.1}",
            label, s.mean, s.std, s.min, s.max
        )
        .unwrap(),
    }
}
