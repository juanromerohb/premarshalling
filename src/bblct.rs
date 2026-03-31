//! Branch-and-bound solver for the CPMP-LCT.
//!
//! This module combines heuristic incumbents, dominance rules, and the
//! crane-time lower bound into a depth-first exact search.

use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::Arc;

use crate::bay::Bay;
use crate::crane::CraneParams;
use crate::heuristics::{glct, hlct};
use crate::lower_bound;
use crate::solution::Solution;

#[derive(Clone, Debug, PartialEq)]
/// Heuristic used to initialize or strengthen the branch-and-bound search.
pub enum HeuristicChoice {
    None,
    Glct,
    Hlct,
}

#[derive(Clone, Debug)]
/// Configuration of the BBLCT branch-and-bound search.
pub struct BblctConfig {
    pub initial_heuristic: HeuristicChoice,
    pub use_dominance: bool,
    pub use_lb_pruning: bool,
    pub use_lb_sorting: bool,
    pub node_heuristic: HeuristicChoice,
    pub node_heuristic_every: usize,
}

impl Default for BblctConfig {
    fn default() -> Self {
        Self {
            initial_heuristic: HeuristicChoice::Hlct,
            use_dominance: true,
            use_lb_pruning: true,
            use_lb_sorting: true,
            node_heuristic: HeuristicChoice::None,
            node_heuristic_every: 1,
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
/// Reason why the search stopped.
pub enum BblctStopReason {
    Completed,
    Cancelled,
}

#[derive(Debug)]
/// Basic search statistics collected by BBLCT.
pub struct BblctStats {
    pub nodes_explored: u64,
    pub nodes_pruned: u64,
    pub stop_reason: BblctStopReason,
}

impl Default for BblctStats {
    fn default() -> Self {
        Self {
            nodes_explored: 0,
            nodes_pruned: 0,
            stop_reason: BblctStopReason::Completed,
        }
    }
}

/// Runs the branch-and-bound solver to completion.
pub fn bblct(
    bay: &Bay,
    crane: &CraneParams,
    tau: f64,
    config: &BblctConfig,
) -> (Solution, BblctStats) {
    bblct_cancellable(bay, crane, tau, config, None)
}

/// Runs the branch-and-bound solver with an optional cancellation flag.
///
/// When `cancel` becomes `true`, the current search is interrupted and the best
/// incumbent found so far is returned together with a cancelled stop reason.
pub fn bblct_cancellable(
    bay: &Bay,
    crane: &CraneParams,
    tau: f64,
    config: &BblctConfig,
    cancel: Option<Arc<AtomicBool>>,
) -> (Solution, BblctStats) {
    let mut bay = *bay;
    let mut pi = Solution::new();
    let mut pi_best = run_heuristic(&config.initial_heuristic, &bay, crane, tau);
    let mut stats = BblctStats::default();

    let cancelled = branch(
        &mut bay,
        crane,
        tau,
        &mut pi,
        &mut pi_best,
        &mut stats,
        config,
        &cancel,
    );
    stats.stop_reason = if cancelled {
        BblctStopReason::Cancelled
    } else {
        BblctStopReason::Completed
    };

    (pi_best, stats)
}

fn run_heuristic(choice: &HeuristicChoice, bay: &Bay, crane: &CraneParams, tau: f64) -> Solution {
    match choice {
        HeuristicChoice::None => {
            let mut sol = Solution::new();
            sol.acc = bay.acc();
            sol
        }
        HeuristicChoice::Glct => {
            let mut bay_copy = *bay;
            glct::glct(&mut bay_copy, crane, tau)
        }
        HeuristicChoice::Hlct => hlct::hlct(bay, crane, tau),
    }
}

fn merge_with_retimed_head(prefix: &Solution, tail: &Solution, crane: &CraneParams) -> Solution {
    let mut combined = prefix.clone();
    if tail.is_empty() {
        return combined;
    }

    combined.reserve(tail.len());

    let dm0 = tail.moves[0];
    let time0 = crane.move_time(
        dm0.src + 1,
        dm0.src_tier + 1,
        dm0.dst + 1,
        dm0.dst_tier + 1,
        combined.last_dst_1based(),
    );
    combined.push(dm0, time0);

    for idx in 1..tail.len() {
        combined.push(tail.moves[idx], tail.time_delta(idx));
    }

    combined
}

fn try_node_heuristic(
    config: &BblctConfig,
    bay: &Bay,
    crane: &CraneParams,
    tau: f64,
    pi: &Solution,
    pi_best: &mut Solution,
) {
    if config.node_heuristic == HeuristicChoice::None {
        return;
    }
    let depth = pi.len();
    if config.node_heuristic_every == 0 || depth % config.node_heuristic_every != 0 {
        return;
    }

    let tau_remaining = tau - pi.crane_time;
    if tau_remaining <= 0.0 {
        return;
    }

    let hsol = run_heuristic(&config.node_heuristic, bay, crane, tau_remaining);
    let mut combined = merge_with_retimed_head(pi, &hsol, crane);
    combined.acc = hsol.acc;

    if combined.crane_time > tau {
        return;
    }

    if combined.acc > pi_best.acc
        || (combined.acc == pi_best.acc && combined.crane_time < pi_best.crane_time)
    {
        *pi_best = combined;
    }
}

fn branch(
    bay: &mut Bay,
    crane: &CraneParams,
    tau: f64,
    pi: &mut Solution,
    pi_best: &mut Solution,
    stats: &mut BblctStats,
    config: &BblctConfig,
    cancel: &Option<Arc<AtomicBool>>,
) -> bool {
    if stats.nodes_explored & 0x3FF == 0 {
        if let Some(flag) = cancel {
            if flag.load(Ordering::Relaxed) {
                return true;
            }
        }
    }

    stats.nodes_explored += 1;

    let current_acc = bay.acc();
    if current_acc > pi_best.acc
        || (current_acc == pi_best.acc && pi.crane_time < pi_best.crane_time)
    {
        pi.acc = current_acc;
        *pi_best = pi.clone();
        pi_best.acc = current_acc;
    }

    if current_acc == bay.c {
        return false;
    }

    try_node_heuristic(config, bay, crane, tau, pi, pi_best);

    let mut children: Vec<(usize, usize, f64)> =
        Vec::with_capacity(bay.s * bay.s.saturating_sub(1));
    let prev_dst_1based = pi.last_dst_1based();
    let lb_target = pi_best.acc + 1;
    let last_move = pi.moves.last().copied();

    for src in 0..bay.s {
        if bay.is_empty(src) {
            continue;
        }
        for dst in 0..bay.s {
            if src == dst || bay.is_full(dst) {
                continue;
            }
            if let Some(est_time) = prune(
                bay,
                crane,
                tau,
                pi,
                src,
                dst,
                prev_dst_1based,
                lb_target,
                last_move,
                config,
            ) {
                children.push((src, dst, est_time));
            }
        }
    }

    if config.use_lb_sorting {
        children.sort_by(|a, b| a.2.partial_cmp(&b.2).unwrap());
    }

    for (src, dst, _) in children {
        let m = crate::bay::Move { src, dst };
        let dm = bay.apply_move(m);
        let time = crane.move_time(
            dm.src + 1,
            dm.src_tier + 1,
            dm.dst + 1,
            dm.dst_tier + 1,
            prev_dst_1based,
        );
        pi.push(dm, time);

        let cancelled = branch(bay, crane, tau, pi, pi_best, stats, config, cancel);

        pi.pop();
        bay.undo_move(dm);
        if cancelled {
            return true;
        }
    }

    false
}

fn prune(
    bay: &mut Bay,
    crane: &CraneParams,
    tau: f64,
    pi: &Solution,
    src: usize,
    dst: usize,
    prev_dst_1based: usize,
    lb_target: usize,
    last_move: Option<crate::bay::DetailedMove>,
    config: &BblctConfig,
) -> Option<f64> {
    if src == dst || bay.is_empty(src) || bay.is_full(dst) {
        return None;
    }

    if config.use_dominance {
        if let Some(last) = last_move {
            if src == last.dst {
                return None;
            }
            if dst == last.src && bay.top_group(src) == last.group {
                return None;
            }
        }
    }

    let m = crate::bay::Move { src, dst };
    let dm = bay.apply_move(m);
    let time = crane.move_time(
        dm.src + 1,
        dm.src_tier + 1,
        dm.dst + 1,
        dm.dst_tier + 1,
        prev_dst_1based,
    );
    let new_time = pi.crane_time + time;

    if new_time > tau {
        bay.undo_move(dm);
        return None;
    }

    let est_time = if config.use_lb_pruning {
        let lb = lower_bound::lower_bound(bay, crane, lb_target);
        let est = new_time + lb;

        bay.undo_move(dm);

        if est > tau {
            return None;
        }
        est
    } else {
        bay.undo_move(dm);
        new_time
    };

    Some(est_time)
}
