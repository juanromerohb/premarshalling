//! Hybrid heuristic for the CPMP-LCT.
//!
//! HLCT combines TGH prefixes on `p`-reduced bays with GLCT completions on the
//! original bay.

use crate::bay::{Bay, Move};
use crate::crane::CraneParams;
use crate::solution::Solution;

use super::glct;
use super::tgh;

#[inline]
fn is_better(acc: usize, crane_time: f64, best: &Solution) -> bool {
    acc > best.acc || (acc == best.acc && crane_time < best.crane_time)
}

#[inline]
fn combined_time(prefix: &Solution, tail: &Solution, crane: &CraneParams) -> f64 {
    if tail.is_empty() {
        return prefix.crane_time;
    }

    let dm0 = tail.moves[0];
    let first_delta = crane.move_time(
        dm0.src + 1,
        dm0.src_tier + 1,
        dm0.dst + 1,
        dm0.dst_tier + 1,
        prefix.last_dst_1based(),
    );

    let mut total = prefix.crane_time + first_delta;
    for idx in 1..tail.len() {
        total += tail.time_delta(idx);
    }
    total
}

/// Runs the hybrid limited-crane-time heuristic.
///
/// The input bay is not mutated; every candidate prefix is replayed on copies.
pub fn hlct(bay: &Bay, crane: &CraneParams, tau: f64) -> Solution {
    // Baseline incumbent from GLCT on the original bay.
    let mut bay_base = *bay;
    let mut best = glct::glct(&mut bay_base, crane, tau);

    for p in 1..=bay.p as u8 {
        let reduced = bay.p_reduction(p);
        let mut bay_reduced = reduced;

        let Some(pi_p) = tgh::tgh(&mut bay_reduced, crane) else {
            continue;
        };

        let mut bay_prefix = *bay;
        let mut prefix_sol = Solution::new();
        prefix_sol.reserve(pi_p.len());

        for (i, dm) in pi_p.moves.iter().enumerate() {
            let actual_dm = bay_prefix.apply_move(Move {
                src: dm.src,
                dst: dm.dst,
            });
            prefix_sol.push(actual_dm, pi_p.time_delta(i));

            if prefix_sol.crane_time > tau {
                break;
            }

            let remaining = tau - prefix_sol.crane_time;
            let mut bay_candidate = bay_prefix;
            let tail = if remaining > 0.0 {
                glct::glct(&mut bay_candidate, crane, remaining)
            } else {
                Solution::new()
            };

            let candidate_acc = bay_candidate.acc();
            let candidate_time = combined_time(&prefix_sol, &tail, crane);

            if !is_better(candidate_acc, candidate_time, &best) {
                continue;
            }

            let mut candidate = prefix_sol.clone();
            if !tail.is_empty() {
                candidate.reserve(tail.len());

                let dm0 = tail.moves[0];
                let time0 = crane.move_time(
                    dm0.src + 1,
                    dm0.src_tier + 1,
                    dm0.dst + 1,
                    dm0.dst_tier + 1,
                    prefix_sol.last_dst_1based(),
                );
                candidate.push(dm0, time0);

                for idx in 1..tail.len() {
                    candidate.push(tail.moves[idx], tail.time_delta(idx));
                }
            }
            candidate.acc = candidate_acc;
            best = candidate;
        }
    }

    best
}
