//! Text rendering helpers for bays and solution replays.

use std::fmt::Write;

use crate::bay::{Bay, Move};
use crate::solution::Solution;

/// Pretty-prints a bay using fixed-width stack cells.
pub fn pretty_bay(bay: &Bay) -> String {
    let width = digit_width(bay.p);
    let mut out = String::new();

    for t in (0..bay.h).rev() {
        for si in 0..bay.s {
            if si > 0 {
                out.push(' ');
            }
            if t < bay.he(si) {
                let g = bay.stacks[si][t];
                write!(out, "[{:>w$}]", g, w = width).unwrap();
            } else {
                write!(out, "[{:>w$}]", "", w = width).unwrap();
            }
        }
        out.push('\n');
    }

    if out.ends_with('\n') {
        out.pop();
    }
    out
}

/// Replays a solution step by step as a human-readable text trace.
///
/// `max_group` controls the display width and should be at least the maximum
/// group value that may appear during the replay.
pub fn replay_solution(bay: &Bay, sol: &Solution, max_group: usize) -> String {
    let mut out = String::new();
    let mut current = *bay;
    let width = digit_width(max_group.max(bay.p));

    writeln!(out, "=== Initial ===").unwrap();
    writeln!(out, "{}", pretty_bay_with_width(&current, width)).unwrap();
    writeln!(out, "acc = {}/{}", current.acc(), current.c).unwrap();

    let mut cumulative_time = 0.0;

    for (i, dm) in sol.moves.iter().enumerate() {
        let delta = sol.time_delta(i);
        cumulative_time += delta;

        current.apply_move(Move {
            src: dm.src,
            dst: dm.dst,
        });

        writeln!(out).unwrap();
        writeln!(
            out,
            "--- Move {}: S{} -> S{} (group {}) | +{:.2}s (total: {:.2}s) ---",
            i + 1,
            dm.src,
            dm.dst,
            dm.group,
            delta,
            cumulative_time,
        )
        .unwrap();
        writeln!(out, "{}", pretty_bay_with_width(&current, width)).unwrap();
        writeln!(out, "acc = {}/{}", current.acc(), current.c).unwrap();
    }

    writeln!(out).unwrap();
    writeln!(
        out,
        "=== Final: acc = {}/{} | moves = {} | crane time = {:.2}s ===",
        current.acc(),
        current.c,
        sol.len(),
        sol.crane_time,
    )
    .unwrap();

    out
}

fn pretty_bay_with_width(bay: &Bay, width: usize) -> String {
    let mut out = String::new();

    for t in (0..bay.h).rev() {
        for si in 0..bay.s {
            if si > 0 {
                out.push(' ');
            }
            if t < bay.he(si) {
                let g = bay.stacks[si][t];
                write!(out, "[{:>w$}]", g, w = width).unwrap();
            } else {
                write!(out, "[{:>w$}]", "", w = width).unwrap();
            }
        }
        out.push('\n');
    }

    if out.ends_with('\n') {
        out.pop();
    }
    out
}

fn digit_width(p: usize) -> usize {
    if p >= 100 {
        3
    } else {
        2
    }
}
