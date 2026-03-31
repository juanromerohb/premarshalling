//! Main desktop GUI application.

#[cfg(target_os = "linux")]
use std::io::ErrorKind;
use std::path::PathBuf;
#[cfg(target_os = "linux")]
use std::process::Command;
use std::sync::atomic::Ordering;
use std::sync::mpsc::{Receiver, TryRecvError};
use std::sync::Arc;
use std::time::Duration;

use eframe::egui;
use egui::{Color32, RichText};
use egui_extras::{Column, TableBuilder};

use crate::bay::{Bay, MAX_HEIGHT, MAX_STACKS};
use crate::gui::replay::{build_replay_frames, ReplayFrame};
use crate::gui::run_engine::{
    start_run, GuiAlgorithm, InputMode, InputSource, ProgressSummary, RunHandle, RunOutcome,
    RunRequest, WorkerEvent,
};
use crate::gui::theme;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum UiRunStatus {
    Idle,
    Running,
    Completed,
    Cancelled,
    Failed,
}

struct ActiveRun {
    receiver: Receiver<WorkerEvent>,
    cancel_flag: Arc<std::sync::atomic::AtomicBool>,
    total_instances: usize,
    completed: usize,
    latest: Option<ProgressSummary>,
}

/// Main application state for the desktop GUI.
pub struct CpmpGuiApp {
    input_mode: InputMode,
    algorithm: GuiAlgorithm,
    tau_minutes: f64,
    bblct_cpu_limit_sec: f64,

    manual_name: String,
    manual_s: usize,
    manual_h: usize,
    manual_grid: Vec<Vec<u8>>, // [stack][tier], tier 0 = bottom
    manual_validation_msg: Option<String>,

    file_path: String,
    folder_path: String,

    status: UiRunStatus,
    active_run: Option<ActiveRun>,
    running_algorithm: Option<GuiAlgorithm>,
    last_error: Option<String>,
    outcome: Option<RunOutcome>,

    replay_frames: Vec<ReplayFrame>,
    replay_step: usize,
}

impl CpmpGuiApp {
    /// Creates the desktop GUI application with its default state.
    pub fn new(cc: &eframe::CreationContext<'_>) -> Self {
        theme::apply(&cc.egui_ctx);

        let manual_s = 4usize;
        let manual_h = 4usize;
        Self {
            input_mode: InputMode::Manual,
            algorithm: GuiAlgorithm::Hlct,
            tau_minutes: 15.0,
            bblct_cpu_limit_sec: 60.0,
            manual_name: "manual_instance".to_string(),
            manual_s,
            manual_h,
            manual_grid: vec![vec![0; manual_h]; manual_s],
            manual_validation_msg: None,
            file_path: String::new(),
            folder_path: String::new(),
            status: UiRunStatus::Idle,
            active_run: None,
            running_algorithm: None,
            last_error: None,
            outcome: None,
            replay_frames: Vec::new(),
            replay_step: 0,
        }
    }

    fn ensure_manual_shape(&mut self) {
        self.manual_s = self.manual_s.clamp(1, MAX_STACKS);
        self.manual_h = self.manual_h.clamp(1, MAX_HEIGHT);

        if self.manual_grid.len() != self.manual_s {
            self.manual_grid
                .resize(self.manual_s, vec![0; self.manual_h]);
        }
        for stack in &mut self.manual_grid {
            if stack.len() != self.manual_h {
                stack.resize(self.manual_h, 0);
            }
        }
    }

    fn build_request(&self) -> Result<RunRequest, String> {
        if self.tau_minutes <= 0.0 {
            return Err("Tau must be > 0 minutes.".to_string());
        }

        let source = match self.input_mode {
            InputMode::Manual => {
                let bay = manual_grid_to_bay(&self.manual_grid, self.manual_s, self.manual_h)?;
                InputSource::Manual {
                    name: self.manual_name.trim().to_string(),
                    bay,
                }
            }
            InputMode::File => {
                let path = self.file_path.trim();
                if path.is_empty() {
                    return Err("Please choose an instance file.".to_string());
                }
                InputSource::File { path: path.into() }
            }
            InputMode::Folder => {
                let path = self.folder_path.trim();
                if path.is_empty() {
                    return Err("Please choose an instances folder.".to_string());
                }
                InputSource::Folder { path: path.into() }
            }
        };

        if matches!(self.algorithm, GuiAlgorithm::Bblct) && self.bblct_cpu_limit_sec <= 0.0 {
            return Err("BBLCT CPU limit must be > 0 seconds.".to_string());
        }

        Ok(RunRequest {
            source,
            algorithm: self.algorithm,
            tau_minutes: self.tau_minutes,
            bblct_cpu_limit_sec: self.bblct_cpu_limit_sec,
        })
    }

    fn start_requested_run(&mut self) {
        self.last_error = None;
        self.outcome = None;
        self.replay_frames.clear();
        self.replay_step = 0;

        match self.build_request() {
            Ok(request) => {
                let RunHandle {
                    receiver,
                    cancel_flag,
                } = start_run(request);
                self.active_run = Some(ActiveRun {
                    receiver,
                    cancel_flag,
                    total_instances: 0,
                    completed: 0,
                    latest: None,
                });
                self.running_algorithm = Some(self.algorithm);
                self.status = UiRunStatus::Running;
                self.manual_validation_msg = None;
            }
            Err(err) => {
                self.status = UiRunStatus::Failed;
                self.last_error = Some(err);
            }
        }
    }

    fn poll_worker_events(&mut self) {
        let mut finished: Option<RunOutcome> = None;
        let mut failed: Option<String> = None;

        if let Some(run) = &mut self.active_run {
            loop {
                match run.receiver.try_recv() {
                    Ok(WorkerEvent::Started { total_instances }) => {
                        run.total_instances = total_instances;
                    }
                    Ok(WorkerEvent::Progress {
                        completed,
                        total,
                        summary,
                    }) => {
                        run.completed = completed;
                        run.total_instances = total;
                        run.latest = Some(summary);
                    }
                    Ok(WorkerEvent::Finished(outcome)) => {
                        finished = Some(outcome);
                        break;
                    }
                    Ok(WorkerEvent::Failed(err)) => {
                        failed = Some(err);
                        break;
                    }
                    Err(TryRecvError::Empty) => {
                        break;
                    }
                    Err(TryRecvError::Disconnected) => {
                        failed = Some("Worker disconnected unexpectedly.".to_string());
                        break;
                    }
                }
            }
        }

        if let Some(err) = failed {
            self.status = UiRunStatus::Failed;
            self.last_error = Some(err);
            self.active_run = None;
            self.running_algorithm = None;
            return;
        }

        if let Some(outcome) = finished {
            self.status = if outcome.cancelled {
                UiRunStatus::Cancelled
            } else {
                UiRunStatus::Completed
            };

            if !outcome.errors.is_empty() {
                self.last_error = Some(outcome.errors.join(" | "));
            }

            if matches!(outcome.mode, InputMode::Manual) {
                if let Some(first) = outcome.entries.first() {
                    self.replay_frames =
                        build_replay_frames(&first.initial_bay, &first.result.solution);
                    self.replay_step = 0;
                }
            }

            self.outcome = Some(outcome);
            self.active_run = None;
            self.running_algorithm = None;
        }
    }

    fn draw_mode_and_paths(&mut self, ui: &mut egui::Ui) {
        theme::card_frame().show(ui, |ui| {
            ui.heading("Input");
            ui.horizontal(|ui| {
                ui.selectable_value(&mut self.input_mode, InputMode::Manual, "Manual");
                ui.selectable_value(&mut self.input_mode, InputMode::File, "Single .txt file");
                ui.selectable_value(
                    &mut self.input_mode,
                    InputMode::Folder,
                    "Folder of .txt files",
                );
            });
            ui.separator();

            match self.input_mode {
                InputMode::Manual => {
                    ui.horizontal(|ui| {
                        ui.label("Instance name");
                        ui.text_edit_singleline(&mut self.manual_name);
                    });
                    ui.horizontal(|ui| {
                        ui.label("Stacks (S)");
                        let resp_s = ui.add(
                            egui::DragValue::new(&mut self.manual_s).clamp_range(1..=MAX_STACKS),
                        );
                        ui.label("Height (H)");
                        let resp_h = ui.add(
                            egui::DragValue::new(&mut self.manual_h).clamp_range(1..=MAX_HEIGHT),
                        );
                        if resp_s.changed() || resp_h.changed() {
                            self.ensure_manual_shape();
                        }
                    });
                }
                InputMode::File => {
                    ui.horizontal(|ui| {
                        ui.text_edit_singleline(&mut self.file_path);
                        if ui.button("Browse").clicked() {
                            match pick_file_with_fallback() {
                                PickerResult::Selected(path) => {
                                    self.file_path = path.display().to_string();
                                    self.last_error = None;
                                }
                                PickerResult::Cancelled => {}
                                PickerResult::Unavailable => {
                                    self.last_error = Some(
                                        "Could not open file picker. Install zenity or kdialog."
                                            .to_string(),
                                    );
                                }
                            }
                        }
                    });
                }
                InputMode::Folder => {
                    ui.horizontal(|ui| {
                        ui.text_edit_singleline(&mut self.folder_path);
                        if ui.button("Browse").clicked() {
                            match pick_folder_with_fallback() {
                                PickerResult::Selected(path) => {
                                    self.folder_path = path.display().to_string();
                                    self.last_error = None;
                                }
                                PickerResult::Cancelled => {}
                                PickerResult::Unavailable => {
                                    self.last_error = Some(
                                        "Could not open folder picker. Install zenity or kdialog."
                                            .to_string(),
                                    );
                                }
                            }
                        }
                    });
                }
            }
        });
    }

    fn draw_algo_controls(&mut self, ui: &mut egui::Ui) {
        theme::card_frame().show(ui, |ui| {
            ui.heading("Algorithm and Limits");
            egui::ComboBox::from_id_source("algo-select")
                .selected_text(self.algorithm.short_name())
                .show_ui(ui, |ui| {
                    for algo in [GuiAlgorithm::Bblct, GuiAlgorithm::Hlct, GuiAlgorithm::Glct] {
                        ui.selectable_value(&mut self.algorithm, algo, algo.label());
                    }
                });

            ui.horizontal(|ui| {
                ui.label("Tau (minutes)");
                ui.add(
                    egui::DragValue::new(&mut self.tau_minutes)
                        .speed(0.5)
                        .clamp_range(0.1..=10_000.0),
                );
            });

            if matches!(self.algorithm, GuiAlgorithm::Bblct) {
                ui.horizontal(|ui| {
                    ui.label("BBLCT CPU limit (seconds)");
                    ui.add(
                        egui::DragValue::new(&mut self.bblct_cpu_limit_sec)
                            .speed(1.0)
                            .clamp_range(1.0..=1_000_000.0),
                    );
                });
            }
        });
    }

    fn draw_manual_editor(&mut self, ui: &mut egui::Ui) {
        if !matches!(self.input_mode, InputMode::Manual) {
            return;
        }

        self.ensure_manual_shape();
        theme::card_frame().show(ui, |ui| {
            ui.heading("Manual Grid Editor");
            ui.label("Set container groups by cell. Tier 1 is bottom. 0 means empty.");
            ui.separator();

            egui::Grid::new("manual-grid")
                .spacing(egui::vec2(6.0, 6.0))
                .show(ui, |ui| {
                    ui.label("");
                    for s in 0..self.manual_s {
                        ui.label(RichText::new(format!("S{}", s + 1)).strong());
                    }
                    ui.end_row();

                    for t in (0..self.manual_h).rev() {
                        ui.label(RichText::new(format!("T{}", t + 1)).color(Color32::LIGHT_GRAY));
                        for s in 0..self.manual_s {
                            let value = &mut self.manual_grid[s][t];
                            let bg = theme::group_color(*value);
                            egui::Frame::none()
                                .fill(bg)
                                .stroke(egui::Stroke::new(1.0, theme::BORDER))
                                .rounding(egui::Rounding::same(6.0))
                                .show(ui, |ui| {
                                    ui.set_min_size(egui::vec2(54.0, 24.0));
                                    let mut num = *value as i32;
                                    let resp = ui.add(
                                        egui::DragValue::new(&mut num)
                                            .clamp_range(0..=u8::MAX as i32)
                                            .speed(1),
                                    );
                                    if resp.changed() {
                                        *value = num as u8;
                                    }
                                });
                        }
                        ui.end_row();
                    }
                });

            match validate_manual_grid(&self.manual_grid, self.manual_s, self.manual_h) {
                Ok(()) => {
                    self.manual_validation_msg = None;
                }
                Err(err) => {
                    self.manual_validation_msg = Some(err);
                }
            }

            if let Some(msg) = &self.manual_validation_msg {
                ui.colored_label(theme::WARNING, msg);
            }
        });
    }

    fn draw_run_controls(&mut self, ui: &mut egui::Ui) {
        theme::card_frame().show(ui, |ui| {
            ui.heading("Run");

            let status_label = match self.status {
                UiRunStatus::Idle => "Idle",
                UiRunStatus::Running => "Running",
                UiRunStatus::Completed => "Completed",
                UiRunStatus::Cancelled => "Cancelled",
                UiRunStatus::Failed => "Failed",
            };

            ui.horizontal(|ui| {
                ui.label("Status");
                ui.colored_label(theme::chip_color(status_label), status_label);
            });

            ui.horizontal(|ui| {
                let run_enabled = !matches!(self.status, UiRunStatus::Running);
                if ui
                    .add_enabled(run_enabled, egui::Button::new("Run"))
                    .clicked()
                {
                    self.start_requested_run();
                }

                let can_cancel = matches!(self.status, UiRunStatus::Running)
                    && matches!(self.running_algorithm, Some(GuiAlgorithm::Bblct));
                if ui
                    .add_enabled(can_cancel, egui::Button::new("Cancel BBLCT"))
                    .clicked()
                {
                    if let Some(run) = &self.active_run {
                        run.cancel_flag.store(true, Ordering::Relaxed);
                    }
                }
            });

            if let Some(run) = &self.active_run {
                let frac = if run.total_instances == 0 {
                    0.0
                } else {
                    run.completed as f32 / run.total_instances as f32
                };
                ui.add(
                    egui::ProgressBar::new(frac)
                        .text(format!("{}/{}", run.completed, run.total_instances)),
                );
                if let Some(latest) = &run.latest {
                    ui.label(format!(
                        "Latest: {} | acc={}/{} | moves={} | crane={:.2}s | cpu={:.1}ms",
                        latest.instance_name,
                        latest.acc,
                        latest.c,
                        latest.moves,
                        latest.crane_time,
                        latest.elapsed_ms
                    ));
                    if let Some(opt) = latest.proved_optimal {
                        ui.label(if opt {
                            "BBLCT status: proved optimal"
                        } else {
                            "BBLCT status: not proven optimal"
                        });
                    }
                }
            }

            if let Some(err) = &self.last_error {
                ui.colored_label(theme::DANGER, err);
            }
        });
    }

    fn draw_outcome(&mut self, ui: &mut egui::Ui) {
        let Some(outcome) = self.outcome.take() else {
            return;
        };

        theme::card_frame().show(ui, |ui| {
            ui.heading("Results");
            if outcome.outputs_saved {
                ui.label("Structured outputs written under results/ (file/folder modes).");
            }
            if outcome.cancelled {
                ui.colored_label(
                    theme::WARNING,
                    "Run stopped before finishing all instances.",
                );
            }

            if matches!(outcome.mode, InputMode::Manual) {
                self.draw_manual_result(ui, &outcome);
            } else {
                self.draw_batch_table(ui, &outcome);
            }
        });

        self.outcome = Some(outcome);
    }

    fn draw_manual_result(&mut self, ui: &mut egui::Ui, outcome: &RunOutcome) {
        let Some(entry) = outcome.entries.first() else {
            ui.label("No result available.");
            return;
        };
        let initial_pct = if entry.result.c == 0 {
            100.0
        } else {
            100.0 * entry.initial_acc as f64 / entry.result.c as f64
        };
        let final_pct = if entry.result.c == 0 {
            100.0
        } else {
            100.0 * entry.result.acc as f64 / entry.result.c as f64
        };

        ui.label(format!("Instance: {}", entry.result.instance_name));
        ui.label(format!(
            "Accessible: {} -> {} ({:.1}% -> {:.1}%)",
            entry.initial_acc, entry.result.acc, initial_pct, final_pct
        ));
        ui.label(format!("Crane time: {:.2}s", entry.result.crane_time));
        ui.label(format!("Moves: {}", entry.result.moves));
        ui.label(format!("CPU time: {:.1}ms", entry.result.elapsed_ms));
        if let Some(opt) = entry.proved_optimal {
            ui.label(if opt {
                "Proved optimal: yes"
            } else {
                "Proved optimal: no (interrupted or timeout)"
            });
        }

        ui.separator();
        if self.replay_frames.is_empty() {
            ui.label("No replay frames available.");
            return;
        }

        ui.horizontal(|ui| {
            if ui.button("Prev").clicked() && self.replay_step > 0 {
                self.replay_step -= 1;
            }
            if ui.button("Next").clicked() && self.replay_step + 1 < self.replay_frames.len() {
                self.replay_step += 1;
            }
            ui.add(
                egui::Slider::new(&mut self.replay_step, 0..=self.replay_frames.len() - 1)
                    .text("Step"),
            );
        });

        let frame = &self.replay_frames[self.replay_step];
        if let Some(desc) = &frame.move_desc {
            ui.label(desc);
        } else {
            ui.label("Initial state");
        }
        if let Some(delta) = frame.delta_time {
            ui.label(format!(
                "Delta: +{:.2}s | cumulative: {:.2}s",
                delta, frame.cumulative_time
            ));
        } else {
            ui.label(format!(
                "Cumulative crane time: {:.2}s",
                frame.cumulative_time
            ));
        }
        ui.label(format!("Accessible containers: {}", frame.acc));
        draw_bay(ui, &frame.bay, "manual-replay");
    }

    fn draw_batch_table(&self, ui: &mut egui::Ui, outcome: &RunOutcome) {
        if outcome.entries.is_empty() {
            ui.label("No solved instances.");
            return;
        }

        TableBuilder::new(ui)
            .striped(true)
            .resizable(true)
            .column(Column::auto())
            .column(Column::auto())
            .column(Column::auto())
            .column(Column::auto())
            .column(Column::auto())
            .column(Column::remainder())
            .header(20.0, |mut header| {
                header.col(|ui| {
                    ui.label("Instance");
                });
                header.col(|ui| {
                    ui.label("Acc%");
                });
                header.col(|ui| {
                    ui.label("Moves");
                });
                header.col(|ui| {
                    ui.label("Crane(s)");
                });
                header.col(|ui| {
                    ui.label("CPU(ms)");
                });
                header.col(|ui| {
                    ui.label("Optimal?");
                });
            })
            .body(|body| {
                body.heterogeneous_rows(outcome.entries.iter().map(|_| 22.0), |mut row| {
                    let idx = row.index();
                    let r = &outcome.entries[idx];
                    row.col(|ui| {
                        ui.label(&r.result.instance_name);
                    });
                    row.col(|ui| {
                        ui.label(format!("{:.1}", r.result.acc_pct()));
                    });
                    row.col(|ui| {
                        ui.label(format!("{}", r.result.moves));
                    });
                    row.col(|ui| {
                        ui.label(format!("{:.2}", r.result.crane_time));
                    });
                    row.col(|ui| {
                        ui.label(format!("{:.1}", r.result.elapsed_ms));
                    });
                    row.col(|ui| {
                        ui.label(match r.proved_optimal {
                            Some(true) => "Yes",
                            Some(false) => "No",
                            None => "-",
                        });
                    });
                });
            });
    }

    fn input_panel_width(&self) -> f32 {
        match self.input_mode {
            InputMode::Manual => 460.0,
            InputMode::File | InputMode::Folder => 820.0,
        }
    }

    fn algo_panel_width(&self) -> f32 {
        520.0
    }

    fn manual_grid_panel_width(&self) -> f32 {
        let grid_width = 110.0 + self.manual_s as f32 * 72.0;
        grid_width.clamp(420.0, 1600.0)
    }

    fn run_panel_width(&self) -> f32 {
        if matches!(self.status, UiRunStatus::Running) {
            560.0
        } else {
            520.0
        }
    }

    fn results_panel_width(&self) -> f32 {
        let Some(outcome) = &self.outcome else {
            return 860.0;
        };

        if matches!(outcome.mode, InputMode::Manual) {
            if let Some(first) = outcome.entries.first() {
                let bay_width = 110.0 + first.result.s as f32 * 56.0;
                return bay_width.clamp(520.0, 1800.0);
            }
            860.0
        } else {
            980.0
        }
    }

    fn common_column_width(&self) -> f32 {
        let mut w = self
            .input_panel_width()
            .max(self.algo_panel_width())
            .max(self.run_panel_width());
        if matches!(self.input_mode, InputMode::Manual) {
            w = w.max(self.manual_grid_panel_width());
        }
        if self.outcome.is_some() {
            w = w.max(self.results_panel_width());
        }
        w
    }
}

impl eframe::App for CpmpGuiApp {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        self.poll_worker_events();
        if matches!(self.status, UiRunStatus::Running) {
            ctx.request_repaint_after(Duration::from_millis(40));
        }

        egui::TopBottomPanel::top("header").show(ctx, |ui| {
            ui.add_space(4.0);
            ui.heading(RichText::new("CPMP-LCT Desktop GUI").color(theme::ACCENT));
            ui.label("Manual / File / Folder solving with BBLCT, HLCT and GLCT.");
            ui.add_space(2.0);
        });

        egui::CentralPanel::default().show(ctx, |ui| {
            egui::ScrollArea::vertical()
                .auto_shrink([false, false])
                .show(ui, |ui| {
                    let col_width = self.common_column_width();
                    draw_centered_block(ui, col_width, |ui| self.draw_mode_and_paths(ui));
                    draw_centered_block(ui, col_width, |ui| self.draw_algo_controls(ui));
                    draw_centered_block(ui, col_width, |ui| self.draw_manual_editor(ui));
                    draw_centered_block(ui, col_width, |ui| self.draw_run_controls(ui));
                    draw_centered_block(ui, col_width, |ui| self.draw_outcome(ui));
                });
        });
    }
}

fn draw_centered_block(
    ui: &mut egui::Ui,
    desired_width: f32,
    add_contents: impl FnOnce(&mut egui::Ui),
) {
    let avail = ui.available_width();
    let width = desired_width.clamp(320.0, avail.max(320.0)).min(avail);
    ui.horizontal(|ui| {
        let side_space = ((ui.available_width() - width) * 0.5).max(0.0);
        if side_space > 0.0 {
            ui.add_space(side_space);
        }
        ui.allocate_ui_with_layout(
            egui::vec2(width, 0.0),
            egui::Layout::top_down(egui::Align::Min),
            |ui| add_contents(ui),
        );
    });
}

fn draw_bay(ui: &mut egui::Ui, bay: &Bay, id: &str) {
    egui::Grid::new(id)
        .spacing(egui::vec2(4.0, 4.0))
        .show(ui, |ui| {
            ui.label("");
            for s in 0..bay.s {
                ui.label(RichText::new(format!("S{}", s + 1)).strong());
            }
            ui.end_row();

            for t in (0..bay.h).rev() {
                ui.label(RichText::new(format!("T{}", t + 1)).color(Color32::LIGHT_GRAY));
                for s in 0..bay.s {
                    let g = bay.g(s, t);
                    let bg = theme::group_color(g);
                    let fg = theme::text_on(bg);
                    egui::Frame::none()
                        .fill(bg)
                        .stroke(egui::Stroke::new(1.0, theme::BORDER))
                        .rounding(egui::Rounding::same(5.0))
                        .show(ui, |ui| {
                            ui.set_min_size(egui::vec2(38.0, 22.0));
                            let text = if g == 0 {
                                "".to_string()
                            } else {
                                format!("{}", g)
                            };
                            ui.centered_and_justified(|ui| {
                                ui.label(RichText::new(text).color(fg).monospace());
                            });
                        });
                }
                ui.end_row();
            }
        });
}

enum PickerResult {
    Selected(PathBuf),
    Cancelled,
    Unavailable,
}

fn pick_file_with_fallback() -> PickerResult {
    #[cfg(target_os = "linux")]
    {
        match run_dialog_command("zenity", &["--file-selection", "--file-filter=*.txt"]) {
            PickerResult::Selected(path) => return PickerResult::Selected(path),
            PickerResult::Cancelled => return PickerResult::Cancelled,
            PickerResult::Unavailable => {}
        }
        match run_dialog_command("kdialog", &["--getopenfilename", ".", "*.txt"]) {
            PickerResult::Selected(path) => return PickerResult::Selected(path),
            PickerResult::Cancelled => return PickerResult::Cancelled,
            PickerResult::Unavailable => {}
        }
        return PickerResult::Unavailable;
    }

    #[cfg(not(target_os = "linux"))]
    {
        if let Some(path) = rfd::FileDialog::new()
            .add_filter("txt", &["txt"])
            .pick_file()
        {
            PickerResult::Selected(path)
        } else {
            PickerResult::Cancelled
        }
    }
}

#[cfg(target_os = "linux")]
fn pick_folder_with_fallback() -> PickerResult {
    match run_dialog_command("zenity", &["--file-selection", "--directory"]) {
        PickerResult::Selected(path) => return PickerResult::Selected(path),
        PickerResult::Cancelled => return PickerResult::Cancelled,
        PickerResult::Unavailable => {}
    }
    match run_dialog_command("kdialog", &["--getexistingdirectory", "."]) {
        PickerResult::Selected(path) => return PickerResult::Selected(path),
        PickerResult::Cancelled => return PickerResult::Cancelled,
        PickerResult::Unavailable => {}
    }
    PickerResult::Unavailable
}

#[cfg(not(target_os = "linux"))]
fn pick_folder_with_fallback() -> PickerResult {
    if let Some(path) = rfd::FileDialog::new().pick_folder() {
        PickerResult::Selected(path)
    } else {
        PickerResult::Cancelled
    }
}

#[cfg(target_os = "linux")]
fn run_dialog_command(cmd: &str, args: &[&str]) -> PickerResult {
    let output = match Command::new(cmd).args(args).output() {
        Ok(out) => out,
        Err(e) if e.kind() == ErrorKind::NotFound => return PickerResult::Unavailable,
        Err(_) => return PickerResult::Unavailable,
    };

    if !output.status.success() {
        return PickerResult::Cancelled;
    }

    let path = String::from_utf8_lossy(&output.stdout).trim().to_string();
    if path.is_empty() {
        PickerResult::Cancelled
    } else {
        PickerResult::Selected(PathBuf::from(path))
    }
}

/// Validates a manually entered GUI grid before converting it to a [`Bay`].
pub fn validate_manual_grid(grid: &[Vec<u8>], s: usize, h: usize) -> Result<(), String> {
    if s == 0 || s > MAX_STACKS {
        return Err(format!("S must be in 1..={}", MAX_STACKS));
    }
    if h == 0 || h > MAX_HEIGHT {
        return Err(format!("H must be in 1..={}", MAX_HEIGHT));
    }
    if grid.len() != s {
        return Err("Invalid grid width.".to_string());
    }
    for (si, stack) in grid.iter().enumerate() {
        if stack.len() != h {
            return Err(format!("Invalid stack height at S{}.", si + 1));
        }
        let mut seen_zero = false;
        for (t, &g) in stack.iter().enumerate() {
            if g == 0 {
                seen_zero = true;
            } else if seen_zero {
                return Err(format!(
                    "Floating container at S{}, T{} (tier above an empty slot).",
                    si + 1,
                    t + 1
                ));
            }
        }
    }
    Ok(())
}

/// Converts a validated manual GUI grid into a [`Bay`].
pub fn manual_grid_to_bay(grid: &[Vec<u8>], s: usize, h: usize) -> Result<Bay, String> {
    validate_manual_grid(grid, s, h)?;
    let mut stacks: Vec<Vec<u8>> = vec![Vec::new(); s];
    for si in 0..s {
        for t in 0..h {
            let g = grid[si][t];
            if g == 0 {
                break;
            }
            stacks[si].push(g);
        }
    }
    Ok(Bay::from_vecs(s, h, &stacks))
}
