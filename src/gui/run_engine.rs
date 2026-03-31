//! Background execution engine for the desktop GUI.
//!
//! The GUI uses this module to normalize manual/file/folder inputs, run the
//! selected algorithm on a worker thread, and stream progress events back to the
//! interface.

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::mpsc::{self, Receiver, Sender};
use std::sync::Arc;
use std::thread;
use std::time::{Duration, Instant};

use crate::bay::Bay;
use crate::bblct::{self, BblctStopReason};
use crate::benchmark::{Algorithm, BenchResult};
use crate::crane::CraneParams;
use crate::heuristics::{glct, hlct, tgh};
use crate::instance;
use crate::results_output::{self, MachineInfo};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
/// Source mode selected in the GUI.
pub enum InputMode {
    Manual,
    File,
    Folder,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
/// Algorithm choices exposed by the desktop GUI.
pub enum GuiAlgorithm {
    Hlct,
    Glct,
    Tgh,
    Bblct,
}

impl GuiAlgorithm {
    /// Returns the long label shown in the GUI selector.
    pub fn label(self) -> &'static str {
        match self {
            GuiAlgorithm::Hlct => "HLCT (fast, deterministic, strongest heuristic)",
            GuiAlgorithm::Glct => "GLCT (very fast, deterministic greedy)",
            GuiAlgorithm::Tgh => "TGH (very fast, deterministic CPMPCT)",
            GuiAlgorithm::Bblct => "BBLCT (variable speed, deterministic, optimal/near-optimal)",
        }
    }

    /// Returns the short algorithm name shown in summaries.
    pub fn short_name(self) -> &'static str {
        match self {
            GuiAlgorithm::Hlct => "HLCT",
            GuiAlgorithm::Glct => "GLCT",
            GuiAlgorithm::Tgh => "TGH",
            GuiAlgorithm::Bblct => "BBLCT",
        }
    }

    /// Converts the GUI choice into the corresponding benchmark algorithm.
    pub fn to_algorithm(self) -> Algorithm {
        match self {
            GuiAlgorithm::Hlct => Algorithm::Hlct,
            GuiAlgorithm::Glct => Algorithm::Glct,
            GuiAlgorithm::Tgh => Algorithm::Tgh,
            GuiAlgorithm::Bblct => Algorithm::Bblct {
                config: bblct::BblctConfig::default(),
            },
        }
    }
}

#[derive(Clone)]
/// Input selected by the GUI before a run starts.
pub enum InputSource {
    Manual { name: String, bay: Bay },
    File { path: PathBuf },
    Folder { path: PathBuf },
}

impl InputSource {
    /// Returns the GUI input mode represented by this source.
    pub fn mode(&self) -> InputMode {
        match self {
            InputSource::Manual { .. } => InputMode::Manual,
            InputSource::File { .. } => InputMode::File,
            InputSource::Folder { .. } => InputMode::Folder,
        }
    }
}

#[derive(Clone)]
/// Request submitted by the GUI to the worker thread.
pub struct RunRequest {
    pub source: InputSource,
    pub algorithm: GuiAlgorithm,
    pub tau_minutes: f64,
    pub bblct_cpu_limit_sec: f64,
}

/// Handle returned to the GUI for a running background task.
pub struct RunHandle {
    pub receiver: Receiver<WorkerEvent>,
    pub cancel_flag: Arc<AtomicBool>,
}

#[derive(Clone)]
/// Compact per-instance progress information sent back to the GUI.
pub struct ProgressSummary {
    pub instance_name: String,
    pub acc: usize,
    pub c: usize,
    pub moves: usize,
    pub crane_time: f64,
    pub elapsed_ms: f64,
    pub proved_optimal: Option<bool>,
}

/// Worker-thread messages consumed by the GUI event loop.
pub enum WorkerEvent {
    Started {
        total_instances: usize,
    },
    Progress {
        completed: usize,
        total: usize,
        summary: ProgressSummary,
    },
    Finished(RunOutcome),
    Failed(String),
}

/// Full outcome of one instance processed by the GUI worker.
pub struct InstanceRunOutcome {
    pub initial_bay: Bay,
    pub initial_acc: usize,
    pub result: BenchResult,
    pub proved_optimal: Option<bool>,
}

/// Final outcome of a GUI run, possibly covering multiple instances.
pub struct RunOutcome {
    pub mode: InputMode,
    pub entries: Vec<InstanceRunOutcome>,
    pub cancelled: bool,
    pub outputs_saved: bool,
    pub errors: Vec<String>,
}

/// Starts a worker thread for the requested GUI run.
pub fn start_run(request: RunRequest) -> RunHandle {
    let (tx, rx) = mpsc::channel::<WorkerEvent>();
    let cancel_flag = Arc::new(AtomicBool::new(false));
    let cancel_for_thread = Arc::clone(&cancel_flag);

    thread::spawn(
        move || match execute_run(request, &tx, &cancel_for_thread) {
            Ok(outcome) => {
                let _ = tx.send(WorkerEvent::Finished(outcome));
            }
            Err(err) => {
                let _ = tx.send(WorkerEvent::Failed(err));
            }
        },
    );

    RunHandle {
        receiver: rx,
        cancel_flag,
    }
}

fn execute_run(
    request: RunRequest,
    tx: &Sender<WorkerEvent>,
    global_cancel: &Arc<AtomicBool>,
) -> Result<RunOutcome, String> {
    let mode = request.source.mode();
    let algorithm = request.algorithm.to_algorithm();
    let tau = request.tau_minutes * 60.0;
    let folder_name = format!(
        "{}_{}",
        algorithm.name(),
        results_output::build_folder_name(&algorithm, tau, None)
    );
    let (instances, instance_dir_for_output, outputs_saved) = collect_instances(&request.source)?;
    let machine_info = if outputs_saved {
        Some(MachineInfo::collect())
    } else {
        None
    };

    tx.send(WorkerEvent::Started {
        total_instances: instances.len(),
    })
    .map_err(|e| format!("cannot send start event: {}", e))?;

    let mut entries = Vec::with_capacity(instances.len());
    let mut errors = Vec::new();
    let mut cancelled = false;

    for (idx, (instance_name, bay)) in instances.iter().enumerate() {
        if global_cancel.load(Ordering::Relaxed) {
            cancelled = true;
            break;
        }

        let entry = run_one_instance(
            request.algorithm,
            &algorithm,
            instance_name,
            bay,
            tau,
            request.bblct_cpu_limit_sec,
            global_cancel,
        );

        if outputs_saved {
            let save_res = std::panic::catch_unwind(|| {
                results_output::save_structured_results(
                    &entry.result,
                    bay,
                    &algorithm,
                    &instance_dir_for_output,
                    machine_info
                        .as_ref()
                        .expect("machine info must be available"),
                    entry.initial_acc,
                    &folder_name,
                );
            });
            if save_res.is_err() {
                errors.push(format!(
                    "failed to save structured output for {}",
                    entry.result.instance_name
                ));
            }
        }

        let summary = ProgressSummary {
            instance_name: entry.result.instance_name.clone(),
            acc: entry.result.acc,
            c: entry.result.c,
            moves: entry.result.moves,
            crane_time: entry.result.crane_time,
            elapsed_ms: entry.result.elapsed_ms,
            proved_optimal: entry.proved_optimal,
        };
        entries.push(entry);

        tx.send(WorkerEvent::Progress {
            completed: idx + 1,
            total: instances.len(),
            summary,
        })
        .map_err(|e| format!("cannot send progress event: {}", e))?;

        if global_cancel.load(Ordering::Relaxed) {
            cancelled = true;
            break;
        }
    }

    if outputs_saved && matches!(mode, InputMode::Folder) && !entries.is_empty() {
        let mut by_category: BTreeMap<(usize, usize), Vec<&BenchResult>> = BTreeMap::new();
        for entry in &entries {
            by_category
                .entry((entry.result.h, entry.result.s))
                .or_default()
                .push(&entry.result);
        }
        for results in by_category.values() {
            let write_res = std::panic::catch_unwind(|| {
                results_output::write_category_analysis(
                    results,
                    &algorithm,
                    &instance_dir_for_output,
                    &folder_name,
                );
            });
            if write_res.is_err() {
                errors.push("failed to write category analysis file".to_string());
            }
        }
    }

    Ok(RunOutcome {
        mode,
        entries,
        cancelled,
        outputs_saved,
        errors,
    })
}

fn collect_instances(source: &InputSource) -> Result<(Vec<(String, Bay)>, String, bool), String> {
    match source {
        InputSource::Manual { name, bay } => Ok((vec![(name.clone(), *bay)], String::new(), false)),
        InputSource::File { path } => {
            if !path.exists() {
                return Err(format!("file not found: {}", path.display()));
            }
            let bay = instance::read_instance(
                path.to_str()
                    .ok_or_else(|| format!("invalid UTF-8 path: {}", path.display()))?,
            );
            let name = display_name_for_path(path);
            let dir = path
                .parent()
                .unwrap_or_else(|| Path::new("."))
                .to_string_lossy()
                .to_string();
            Ok((vec![(name, bay)], dir, true))
        }
        InputSource::Folder { path } => {
            if !path.exists() {
                return Err(format!("folder not found: {}", path.display()));
            }
            if !path.is_dir() {
                return Err(format!("not a directory: {}", path.display()));
            }
            let dir = path.to_string_lossy().to_string();
            let instances = instance::read_all_instances(&dir);
            if instances.is_empty() {
                return Err(format!("no .txt instances found in {}", path.display()));
            }
            Ok((instances, dir, true))
        }
    }
}

fn run_one_instance(
    algo_choice: GuiAlgorithm,
    algorithm: &Algorithm,
    instance_name: &str,
    bay: &Bay,
    tau: f64,
    bblct_cpu_limit_sec: f64,
    global_cancel: &Arc<AtomicBool>,
) -> InstanceRunOutcome {
    let crane = CraneParams::new(bay.h);
    let initial_acc = bay.acc();
    let start = Instant::now();

    let (solution, nodes_explored, proved_optimal, stop_reason) = match algo_choice {
        GuiAlgorithm::Hlct => (hlct::hlct(bay, &crane, tau), None, None, None),
        GuiAlgorithm::Glct => {
            let mut bay_copy = *bay;
            (glct::glct(&mut bay_copy, &crane, tau), None, None, None)
        }
        GuiAlgorithm::Tgh => {
            let mut bay_copy = *bay;
            let sol =
                tgh::tgh(&mut bay_copy, &crane).unwrap_or_else(crate::solution::Solution::new);
            (sol, None, None, None)
        }
        GuiAlgorithm::Bblct => {
            let cancel = Arc::new(AtomicBool::new(false));
            let timer_cancel = Arc::clone(&cancel);
            let timer_stop = Arc::new(AtomicBool::new(false));
            let timer_stop_flag = Arc::clone(&timer_stop);
            let limit = bblct_cpu_limit_sec.max(0.0);
            let timer_thread = thread::spawn(move || {
                let start = Instant::now();
                while !timer_stop_flag.load(Ordering::Relaxed) {
                    if start.elapsed().as_secs_f64() >= limit {
                        timer_cancel.store(true, Ordering::Relaxed);
                        break;
                    }
                    thread::sleep(Duration::from_millis(16));
                }
            });

            let bridge_cancel = Arc::clone(&cancel);
            let bridge_global = Arc::clone(global_cancel);
            let bridge_stop = Arc::new(AtomicBool::new(false));
            let bridge_stop_flag = Arc::clone(&bridge_stop);
            let bridge_thread = thread::spawn(move || {
                while !bridge_cancel.load(Ordering::Relaxed) {
                    if bridge_stop_flag.load(Ordering::Relaxed) {
                        break;
                    }
                    if bridge_global.load(Ordering::Relaxed) {
                        bridge_cancel.store(true, Ordering::Relaxed);
                        break;
                    }
                    thread::sleep(Duration::from_millis(16));
                }
            });

            let (sol, stats) = bblct::bblct_cancellable(
                bay,
                &crane,
                tau,
                &bblct::BblctConfig::default(),
                Some(Arc::clone(&cancel)),
            );
            cancel.store(true, Ordering::Relaxed);
            timer_stop.store(true, Ordering::Relaxed);
            bridge_stop.store(true, Ordering::Relaxed);
            let _ = timer_thread.join();
            let _ = bridge_thread.join();
            (
                sol,
                Some(stats.nodes_explored),
                Some(stats.stop_reason == BblctStopReason::Completed),
                Some(
                    match stats.stop_reason {
                        BblctStopReason::Completed => "completed",
                        BblctStopReason::Cancelled => "cancelled",
                    }
                    .to_string(),
                ),
            )
        }
    };

    let elapsed_ms = start.elapsed().as_secs_f64() * 1000.0;
    let result = BenchResult {
        instance_name: instance_name.to_string(),
        algorithm: algorithm.name().to_string(),
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
        nodes_explored,
        proved_optimal,
        stop_reason,
        solution,
    };

    InstanceRunOutcome {
        initial_bay: *bay,
        initial_acc,
        result,
        proved_optimal,
    }
}

fn display_name_for_path(path: &Path) -> String {
    let filename = path
        .file_name()
        .and_then(|s| s.to_str())
        .unwrap_or("instance.txt")
        .to_string();

    let parent = path
        .parent()
        .and_then(|p| p.file_name())
        .and_then(|s| s.to_str())
        .unwrap_or("");

    if parent.starts_with("BF") && parent[2..].chars().all(|c| c.is_ascii_digit()) {
        return format!("{}_{}", parent, filename);
    }

    filename
}
