//! Replay-frame construction for the desktop GUI.

use crate::bay::{Bay, Move};
use crate::solution::Solution;

#[derive(Clone)]
/// Snapshot of one step in the GUI replay timeline.
pub struct ReplayFrame {
    pub step: usize,
    pub bay: Bay,
    pub acc: usize,
    pub move_desc: Option<String>,
    pub delta_time: Option<f64>,
    pub cumulative_time: f64,
}

/// Builds the full replay timeline from an initial bay and a solution.
pub fn build_replay_frames(initial: &Bay, solution: &Solution) -> Vec<ReplayFrame> {
    let mut frames = Vec::with_capacity(solution.len() + 1);
    let mut current = *initial;
    let mut cumulative_time = 0.0;

    frames.push(ReplayFrame {
        step: 0,
        bay: current,
        acc: current.acc(),
        move_desc: None,
        delta_time: None,
        cumulative_time: 0.0,
    });

    for (i, dm) in solution.moves.iter().enumerate() {
        let delta = solution.time_delta(i);
        cumulative_time += delta;

        current.apply_move(Move {
            src: dm.src,
            dst: dm.dst,
        });

        frames.push(ReplayFrame {
            step: i + 1,
            bay: current,
            acc: current.acc(),
            move_desc: Some(format!(
                "Move {}: S{} -> S{} (group {})",
                i + 1,
                dm.src + 1,
                dm.dst + 1,
                dm.group
            )),
            delta_time: Some(delta),
            cumulative_time,
        });
    }

    frames
}
