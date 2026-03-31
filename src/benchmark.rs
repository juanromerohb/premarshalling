//! Benchmark runner and batch execution utilities.

use std::fmt;
use std::fs;
use std::hash::{Hash, Hasher};
use std::io::Write;
use std::path::Path;
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::Arc;
use std::thread;
use std::time::{Duration, Instant, SystemTime};

use crate::bay::Bay;
use crate::bblct::{self, BblctConfig};
use crate::cplct::{self, CplctBackend, CplctConfig, CplctModel};
use crate::crane::CraneParams;
use crate::display;
use crate::heuristics::{glct, hlct, tgh};
use crate::instance;
use crate::iplct;
use crate::results_output::{self, MachineInfo};
use crate::solution::Solution;

#[derive(Clone, Debug)]
/// Algorithm choices exposed by the CLI benchmark runner.
pub enum Algorithm {
    Bblct { config: BblctConfig },
    Glct,
    Hlct,
    Tgh,
    Lct1Cpo { time_limit_sec: f64 },
    Lct2Cpo { time_limit_sec: f64 },
    Lct1Cpsat { time_limit_sec: f64 },
    Lct2Cpsat { time_limit_sec: f64 },
    Iplct { time_limit_sec: f64 },
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
/// Warm-start heuristic option shared by the CP and IP models.
pub enum ModelWarmStartHeuristic {
    Glct,
    Hlct,
}

impl Algorithm {
    /// Returns the short display name used in outputs and reports.
    pub fn name(&self) -> &'static str {
        match self {
            Algorithm::Bblct { .. } => "BBLCT",
            Algorithm::Glct => "GLCT",
            Algorithm::Hlct => "HLCT",
            Algorithm::Tgh => "TGH",
            Algorithm::Lct1Cpo { .. } => "LCT1-CPO",
            Algorithm::Lct2Cpo { .. } => "LCT2-CPO",
            Algorithm::Lct1Cpsat { .. } => "LCT1-CPSAT",
            Algorithm::Lct2Cpsat { .. } => "LCT2-CPSAT",
            Algorithm::Iplct { .. } => "IPLCT",
        }
    }

    /// Parses an algorithm name into a configured benchmark algorithm.
    pub fn from_name(
        name: &str,
        cp_time_limit_sec: f64,
        ip_time_limit_sec: f64,
        bblct_config: &BblctConfig,
    ) -> Option<Self> {
        match name.to_lowercase().as_str() {
            "bblct" => Some(Algorithm::Bblct {
                config: bblct_config.clone(),
            }),
            "glct" => Some(Algorithm::Glct),
            "hlct" => Some(Algorithm::Hlct),
            "tgh" => Some(Algorithm::Tgh),
            "lct1-cpo" | "lct1_cpo" | "lct1cpo" => Some(Algorithm::Lct1Cpo {
                time_limit_sec: cp_time_limit_sec,
            }),
            "lct2-cpo" | "lct2_cpo" | "lct2cpo" => Some(Algorithm::Lct2Cpo {
                time_limit_sec: cp_time_limit_sec,
            }),
            "lct1-cpsat" | "lct1_cpsat" | "lct1cpsat" => Some(Algorithm::Lct1Cpsat {
                time_limit_sec: cp_time_limit_sec,
            }),
            "lct2-cpsat" | "lct2_cpsat" | "lct2cpsat" => Some(Algorithm::Lct2Cpsat {
                time_limit_sec: cp_time_limit_sec,
            }),
            "iplct" => Some(Algorithm::Iplct {
                time_limit_sec: ip_time_limit_sec,
            }),
            _ => None,
        }
    }

    /// All algorithms with default parameters.
    pub fn all(
        cp_time_limit_sec: f64,
        ip_time_limit_sec: f64,
        bblct_config: &BblctConfig,
    ) -> Vec<Self> {
        vec![
            Algorithm::Tgh,
            Algorithm::Hlct,
            Algorithm::Glct,
            Algorithm::Bblct {
                config: bblct_config.clone(),
            },
            Algorithm::Lct1Cpo {
                time_limit_sec: cp_time_limit_sec,
            },
            Algorithm::Lct2Cpo {
                time_limit_sec: cp_time_limit_sec,
            },
            Algorithm::Lct1Cpsat {
                time_limit_sec: cp_time_limit_sec,
            },
            Algorithm::Lct2Cpsat {
                time_limit_sec: cp_time_limit_sec,
            },
            Algorithm::Iplct {
                time_limit_sec: ip_time_limit_sec,
            },
        ]
    }
}

impl fmt::Display for Algorithm {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.name())
    }
}

/// Result of running one algorithm on one instance.
pub struct BenchResult {
    pub instance_name: String,
    pub algorithm: String,
    pub s: usize,
    pub h: usize,
    pub c: usize,
    pub p: usize,
    pub tau: f64,
    pub initial_acc: usize,
    pub acc: usize,
    pub moves: usize,
    pub crane_time: f64,
    pub elapsed_ms: f64,
    pub nodes_explored: Option<u64>,
    pub proved_optimal: Option<bool>,
    pub stop_reason: Option<String>,
    pub solution: Solution,
}

impl BenchResult {
    /// Returns the final accessibility percentage.
    pub fn acc_pct(&self) -> f64 {
        if self.c == 0 {
            100.0
        } else {
            100.0 * self.acc as f64 / self.c as f64
        }
    }
}

fn map_cplct_warm_start(
    warm_start: Option<ModelWarmStartHeuristic>,
) -> Option<cplct::WarmStartHeuristic> {
    match warm_start {
        Some(ModelWarmStartHeuristic::Glct) => Some(cplct::WarmStartHeuristic::Glct),
        Some(ModelWarmStartHeuristic::Hlct) => Some(cplct::WarmStartHeuristic::Hlct),
        None => None,
    }
}

fn map_iplct_warm_start(
    warm_start: Option<ModelWarmStartHeuristic>,
) -> Option<iplct::WarmStartHeuristic> {
    match warm_start {
        Some(ModelWarmStartHeuristic::Glct) => Some(iplct::WarmStartHeuristic::Glct),
        Some(ModelWarmStartHeuristic::Hlct) => Some(iplct::WarmStartHeuristic::Hlct),
        None => None,
    }
}

/// Runs one algorithm on one instance and collects benchmark metrics.
pub fn run_single(
    algo: &Algorithm,
    instance_name: &str,
    bay: &Bay,
    crane: &CraneParams,
    tau: f64,
    bblct_cpu_limit_sec: Option<f64>,
    model_warm_start: Option<ModelWarmStartHeuristic>,
) -> BenchResult {
    let initial_acc = bay.acc();
    let start = Instant::now();

    let (solution, nodes, proved_optimal, stop_reason) = match algo {
        Algorithm::Bblct { config } => {
            if let Some(limit_sec) = bblct_cpu_limit_sec {
                let limit = limit_sec.max(0.0);
                let cancel = Arc::new(AtomicBool::new(false));

                let timer_thread = if limit <= 0.0 {
                    cancel.store(true, Ordering::Relaxed);
                    None
                } else {
                    let timer_cancel = Arc::clone(&cancel);
                    Some(thread::spawn(move || {
                        let start = Instant::now();
                        while !timer_cancel.load(Ordering::Relaxed) {
                            if start.elapsed().as_secs_f64() >= limit {
                                timer_cancel.store(true, Ordering::Relaxed);
                                break;
                            }
                            thread::sleep(Duration::from_millis(8));
                        }
                    }))
                };

                let (sol, stats) =
                    bblct::bblct_cancellable(bay, crane, tau, config, Some(Arc::clone(&cancel)));
                cancel.store(true, Ordering::Relaxed);
                if let Some(timer_thread) = timer_thread {
                    let _ = timer_thread.join();
                }

                let (proved_optimal, stop_reason) = match stats.stop_reason {
                    bblct::BblctStopReason::Completed => {
                        (Some(true), Some("completed".to_string()))
                    }
                    bblct::BblctStopReason::Cancelled => {
                        (Some(false), Some("cpu_limit".to_string()))
                    }
                };
                (sol, Some(stats.nodes_explored), proved_optimal, stop_reason)
            } else {
                let (sol, stats) = bblct::bblct(bay, crane, tau, config);
                (
                    sol,
                    Some(stats.nodes_explored),
                    Some(true),
                    Some("completed".to_string()),
                )
            }
        }
        Algorithm::Glct => {
            let mut bay_copy = *bay;
            let sol = glct::glct(&mut bay_copy, crane, tau);
            (sol, None, None, None)
        }
        Algorithm::Hlct => {
            let sol = hlct::hlct(bay, crane, tau);
            (sol, None, None, None)
        }
        Algorithm::Tgh => {
            let mut bay_copy = *bay;
            let sol = tgh::tgh(&mut bay_copy, crane).unwrap_or_else(Solution::new);
            (sol, None, None, None)
        }
        Algorithm::Lct1Cpo { time_limit_sec } => {
            let c = CplctConfig {
                model: CplctModel::Lct1,
                backend: CplctBackend::CpOptimizer,
                tau,
                time_limit_sec: *time_limit_sec,
                warm_start_heuristic: map_cplct_warm_start(model_warm_start),
            };
            let res = cplct::solve_cplct(bay, crane, &c);
            let stop_reason = Some(format!("solver_status:{}", res.status.to_ascii_lowercase()));
            (res.solution, None, Some(res.optimal), stop_reason)
        }
        Algorithm::Lct2Cpo { time_limit_sec } => {
            let c = CplctConfig {
                model: CplctModel::Lct2,
                backend: CplctBackend::CpOptimizer,
                tau,
                time_limit_sec: *time_limit_sec,
                warm_start_heuristic: map_cplct_warm_start(model_warm_start),
            };
            let res = cplct::solve_cplct(bay, crane, &c);
            let stop_reason = Some(format!("solver_status:{}", res.status.to_ascii_lowercase()));
            (res.solution, None, Some(res.optimal), stop_reason)
        }
        Algorithm::Lct1Cpsat { time_limit_sec } => {
            let c = CplctConfig {
                model: CplctModel::Lct1,
                backend: CplctBackend::CpSat,
                tau,
                time_limit_sec: *time_limit_sec,
                warm_start_heuristic: map_cplct_warm_start(model_warm_start),
            };
            let res = cplct::solve_cplct(bay, crane, &c);
            let stop_reason = Some(format!("solver_status:{}", res.status.to_ascii_lowercase()));
            (res.solution, None, Some(res.optimal), stop_reason)
        }
        Algorithm::Lct2Cpsat { time_limit_sec } => {
            let c = CplctConfig {
                model: CplctModel::Lct2,
                backend: CplctBackend::CpSat,
                tau,
                time_limit_sec: *time_limit_sec,
                warm_start_heuristic: map_cplct_warm_start(model_warm_start),
            };
            let res = cplct::solve_cplct(bay, crane, &c);
            let stop_reason = Some(format!("solver_status:{}", res.status.to_ascii_lowercase()));
            (res.solution, None, Some(res.optimal), stop_reason)
        }
        Algorithm::Iplct { time_limit_sec } => {
            let lp_path = tmp_lp_path("iplct", instance_name);
            let res = iplct::generate_and_solve_with_solution_options(
                bay,
                crane,
                tau,
                &lp_path,
                *time_limit_sec,
                iplct::IplctSolveOptions {
                    warm_start: map_iplct_warm_start(model_warm_start),
                },
            );
            let stop_reason = Some(if res.optimal {
                "solver_status:optimal".to_string()
            } else {
                "solver_status:not_proved_optimal".to_string()
            });
            (res.solution, None, Some(res.optimal), stop_reason)
        }
    };

    let elapsed_ms = start.elapsed().as_secs_f64() * 1000.0;

    BenchResult {
        instance_name: instance_name.to_string(),
        algorithm: algo.name().to_string(),
        s: bay.s,
        h: bay.h,
        c: bay.c,
        p: bay.p,
        tau,
        initial_acc,
        acc: solution.acc,
        moves: solution.len(),
        crane_time: solution.crane_time,
        elapsed_ms,
        nodes_explored: nodes,
        proved_optimal,
        stop_reason,
        solution,
    }
}

fn seeded_hash(name: &str, salt: u64) -> u64 {
    let mut hasher = std::collections::hash_map::DefaultHasher::new();
    name.hash(&mut hasher);
    salt.hash(&mut hasher);
    hasher.finish()
}

fn tmp_lp_path(prefix: &str, instance_name: &str) -> String {
    let stamp = SystemTime::now()
        .duration_since(SystemTime::UNIX_EPOCH)
        .map(|d| d.as_nanos())
        .unwrap_or(0);
    let hash = seeded_hash(instance_name, 0x49504C4354_u64);
    let _ = fs::create_dir_all("tmp-agent/iplct");
    format!(
        "tmp-agent/iplct/{}_{}_{}_{}.lp",
        prefix,
        instance_name.replace('/', "_"),
        hash,
        stamp
    )
}

/// Configuration of a benchmark batch run.
pub struct BenchConfig {
    pub instance_dir: String,
    pub instance_file: Option<String>,
    pub algorithms: Vec<Algorithm>,
    pub tau: f64,
    pub max_per_category: Option<usize>,
    pub output_dir: Option<String>,
    pub csv_path: Option<String>,
    pub save_solutions: bool,
    pub save_results: bool,
    pub run_tag: Option<String>,
    pub skip_category_analysis: bool,
    pub bblct_cpu_limit_sec: Option<f64>,
    pub model_warm_start: Option<ModelWarmStartHeuristic>,
}

/// Runs the requested algorithms on the selected instances.
pub fn run_benchmark(config: &BenchConfig) -> Vec<BenchResult> {
    let mut work_items: Vec<(usize, usize, Vec<(String, Bay)>)> = Vec::new();
    if let Some(ref file_path) = config.instance_file {
        let p = Path::new(file_path);
        let bay = instance::read_instance(file_path);
        let name = instance_name_from_path(p);
        let (h, s) = instance::try_parse_hs(&name).unwrap_or((bay.h, bay.s));
        work_items.push((h, s, vec![(name, bay)]));
    } else {
        let categories = instance::list_categories(&config.instance_dir);
        if categories.is_empty() {
            eprintln!("No categories found in {}", config.instance_dir);
            return Vec::new();
        }
        for &(h, s) in &categories {
            work_items.push((h, s, instance::read_category(&config.instance_dir, h, s)));
        }
    }

    eprintln!(
        "Benchmark: {} | tau={:.0}s | algorithms: {} | categories: {}",
        config.instance_dir,
        config.tau,
        config
            .algorithms
            .iter()
            .map(|a| a.name())
            .collect::<Vec<_>>()
            .join(", "),
        work_items
            .iter()
            .map(|(h, s, _)| format!("H{}S{}", h, s))
            .collect::<Vec<_>>()
            .join(", "),
    );

    let mut all_results = Vec::new();

    let machine_info = if config.save_results {
        Some(MachineInfo::collect())
    } else {
        None
    };

    let folder_names: Vec<String> = if config.save_results {
        config
            .algorithms
            .iter()
            .map(|a| {
                let algo_name = a.name();
                let tag =
                    results_output::build_folder_name(a, config.tau, config.run_tag.as_deref());
                format!("{}_{}", algo_name, tag)
            })
            .collect()
    } else {
        Vec::new()
    };

    for (h, s, instances) in work_items {
        let limit = config.max_per_category.unwrap_or(instances.len());
        let cat_name = format!("H{}S{}", h, s);

        eprintln!(
            "\n--- {} ({} instances) ---",
            cat_name,
            instances.len().min(limit)
        );

        for (name, bay) in instances.iter().take(limit) {
            let crane = CraneParams::new(bay.h);

            for algo in &config.algorithms {
                let result = run_single(
                    algo,
                    name,
                    bay,
                    &crane,
                    config.tau,
                    config.bblct_cpu_limit_sec,
                    config.model_warm_start,
                );

                eprint!(
                    "  {:<24} {:<6} acc={:>3}/{:<3} ({:>5.1}%) moves={:<3} crane={:>8.2}s  {:.1}ms",
                    name,
                    algo.name(),
                    result.acc,
                    result.c,
                    result.acc_pct(),
                    result.moves,
                    result.crane_time,
                    result.elapsed_ms,
                );
                if let Some(nodes) = result.nodes_explored {
                    eprint!("  nodes={}", nodes);
                }
                eprintln!();

                if config.save_solutions {
                    if let Some(ref output_dir) = config.output_dir {
                        save_solution_file(&result, bay, output_dir);
                    }
                }

                if config.save_results {
                    let algo_idx = config
                        .algorithms
                        .iter()
                        .position(|a| std::ptr::eq(a, algo))
                        .unwrap();
                    results_output::save_structured_results(
                        &result,
                        bay,
                        algo,
                        &config.instance_dir,
                        machine_info.as_ref().unwrap(),
                        result.initial_acc,
                        &folder_names[algo_idx],
                    );
                }

                all_results.push(result);
            }
        }

        if config.save_results && !config.skip_category_analysis {
            let cat_start =
                all_results.len() - instances.len().min(limit) * config.algorithms.len();
            for (algo_idx, algo) in config.algorithms.iter().enumerate() {
                let cat_results: Vec<&BenchResult> = all_results[cat_start..]
                    .iter()
                    .skip(algo_idx)
                    .step_by(config.algorithms.len())
                    .collect();
                if !cat_results.is_empty() {
                    results_output::write_category_analysis(
                        &cat_results,
                        algo,
                        &config.instance_dir,
                        &folder_names[algo_idx],
                    );
                }
            }
        }
    }

    all_results
}

/// Saves the benchmark results as a CSV file.
pub fn save_results_csv(results: &[BenchResult], path: &str) {
    let parent = Path::new(path).parent();
    if let Some(dir) = parent {
        if !dir.as_os_str().is_empty() {
            fs::create_dir_all(dir).expect("cannot create output directory");
        }
    }

    let mut f = fs::File::create(path).expect("cannot create CSV file");
    writeln!(
        f,
        "instance,algorithm,s,h,c,p,tau,initial_acc,initial_acc_pct,acc,acc_pct,moves,crane_time,elapsed_ms,nodes,proved_optimal,stop_reason"
    )
    .unwrap();

    for r in results {
        writeln!(
            f,
            "{},{},{},{},{},{},{:.1},{},{:.1},{},{:.1},{},{:.2},{:.1},{},{},{}",
            r.instance_name,
            r.algorithm,
            r.s,
            r.h,
            r.c,
            r.p,
            r.tau,
            r.initial_acc,
            if r.c == 0 {
                100.0
            } else {
                100.0 * r.initial_acc as f64 / r.c as f64
            },
            r.acc,
            r.acc_pct(),
            r.moves,
            r.crane_time,
            r.elapsed_ms,
            r.nodes_explored.map_or(String::new(), |n| n.to_string()),
            r.proved_optimal.map_or(String::new(), |v| v.to_string()),
            r.stop_reason.as_deref().unwrap_or(""),
        )
        .unwrap();
    }

    eprintln!("CSV saved to {}", path);
}

fn instance_name_from_path(path: &Path) -> String {
    let filename = path
        .file_name()
        .and_then(|s| s.to_str())
        .unwrap_or("unknown.txt")
        .to_string();
    let parent = path
        .parent()
        .and_then(|p| p.file_name())
        .and_then(|s| s.to_str())
        .unwrap_or("");
    if parent.starts_with("BF") && parent[2..].chars().all(|c| c.is_ascii_digit()) {
        format!("{}_{}", parent, filename)
    } else {
        filename
    }
}

/// Saves one human-readable solution replay file.
pub fn save_solution_file(result: &BenchResult, bay: &Bay, output_dir: &str) {
    fs::create_dir_all(output_dir).expect("cannot create output directory");

    let stem = result
        .instance_name
        .strip_suffix(".txt")
        .unwrap_or(&result.instance_name);
    let filename = format!("{}_{}.txt", stem, result.algorithm);
    let path = format!("{}/{}", output_dir, filename);

    let mut f = fs::File::create(&path).expect("cannot create solution file");

    // Header
    writeln!(f, "# Instance: {}", result.instance_name).unwrap();
    writeln!(f, "# Algorithm: {}", result.algorithm).unwrap();
    writeln!(
        f,
        "# S={} H={} C={} P={}",
        result.s, result.h, result.c, result.p
    )
    .unwrap();
    writeln!(f, "# tau = {:.1}s", result.tau).unwrap();
    writeln!(
        f,
        "# Result: acc={}/{} ({:.1}%) moves={} crane_time={:.2}s",
        result.acc,
        result.c,
        result.acc_pct(),
        result.moves,
        result.crane_time,
    )
    .unwrap();
    writeln!(f, "# Elapsed: {:.1}ms", result.elapsed_ms).unwrap();
    writeln!(f).unwrap();

    let replay = display::replay_solution(bay, &result.solution, bay.p);
    write!(f, "{}", replay).unwrap();
}

/// Prints a compact summary table grouped by `(H, S)` category.
pub fn print_summary_table(results: &[BenchResult]) {
    if results.is_empty() {
        return;
    }

    let mut algo_names: Vec<String> = Vec::new();
    for r in results {
        if !algo_names.contains(&r.algorithm) {
            algo_names.push(r.algorithm.clone());
        }
    }

    let mut categories: Vec<(usize, usize)> = Vec::new();
    for r in results {
        let cat = (r.h, r.s);
        if !categories.contains(&cat) {
            categories.push(cat);
        }
    }
    categories.sort();

    println!("\n=== Summary by category ===\n");

    // Header
    print!("{:<10}", "Category");
    for algo in &algo_names {
        print!(
            "  {:>8} {:>6} {:>8} {:>10}",
            format!("{}_acc%", algo),
            format!("{}_mv", algo),
            format!("{}_ct", algo),
            format!("{}_ms", algo),
        );
    }
    println!();
    println!("{}", "-".repeat(10 + algo_names.len() * 36));

    for &(h, s) in &categories {
        print!("H{}S{:<6}", h, s);

        for algo in &algo_names {
            let cat_results: Vec<&BenchResult> = results
                .iter()
                .filter(|r| r.h == h && r.s == s && r.algorithm == *algo)
                .collect();

            if cat_results.is_empty() {
                print!("  {:>8} {:>6} {:>8} {:>10}", "-", "-", "-", "-");
                continue;
            }

            let n = cat_results.len() as f64;
            let avg_acc = cat_results.iter().map(|r| r.acc_pct()).sum::<f64>() / n;
            let avg_moves = cat_results.iter().map(|r| r.moves).sum::<usize>() as f64 / n;
            let avg_ct = cat_results.iter().map(|r| r.crane_time).sum::<f64>() / n;
            let avg_ms = cat_results.iter().map(|r| r.elapsed_ms).sum::<f64>() / n;

            print!(
                "  {:>7.1}% {:>6.1} {:>7.1}s {:>9.1}",
                avg_acc, avg_moves, avg_ct, avg_ms,
            );
        }
        println!();
    }

    // Grand totals
    println!("{}", "-".repeat(10 + algo_names.len() * 36));
    print!("{:<10}", "TOTAL");
    for algo in &algo_names {
        let algo_results: Vec<&BenchResult> =
            results.iter().filter(|r| r.algorithm == *algo).collect();

        if algo_results.is_empty() {
            print!("  {:>8} {:>6} {:>8} {:>10}", "-", "-", "-", "-");
            continue;
        }

        let n = algo_results.len() as f64;
        let avg_acc = algo_results.iter().map(|r| r.acc_pct()).sum::<f64>() / n;
        let avg_moves = algo_results.iter().map(|r| r.moves).sum::<usize>() as f64 / n;
        let avg_ct = algo_results.iter().map(|r| r.crane_time).sum::<f64>() / n;
        let avg_ms = algo_results.iter().map(|r| r.elapsed_ms).sum::<f64>() / n;

        print!(
            "  {:>7.1}% {:>6.1} {:>7.1}s {:>9.1}",
            avg_acc, avg_moves, avg_ct, avg_ms,
        );
    }
    println!("\n");
}
