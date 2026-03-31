//! Integer-programming model for the CPMP-LCT.
//!
//! This module emits an LP formulation, optionally creates Gurobi warm starts
//! from GLCT or HLCT, and reconstructs concrete move sequences from Gurobi
//! solutions.

use crate::bay::Bay;
use crate::bay::Move;
use crate::crane::CraneParams;
use crate::heuristics::{glct, hlct};
use crate::solution::Solution;
use std::collections::BTreeMap;
use std::fmt::Write;
use std::fs;

/// Result returned by the IPLCT solve helpers.
pub struct IplctResult {
    pub acc: usize,
    pub optimal: bool,
    pub time_seconds: f64,
    pub solution: Solution,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
/// Heuristic used to seed Gurobi with a warm-start solution.
pub enum WarmStartHeuristic {
    Glct,
    Hlct,
}

#[derive(Clone, Copy, Debug, Default)]
/// Optional parameters for the IPLCT solve wrappers.
pub struct IplctSolveOptions {
    pub warm_start: Option<WarmStartHeuristic>,
}

/// Generates the LP model text and writes it to `output_path`.
pub fn generate_lp(bay: &Bay, crane: &CraneParams, tau: f64, output_path: &str) -> String {
    let ns = bay.s;
    let nh = bay.h;
    let np = bay.p;
    let nt = compute_t_max(bay, crane, tau);

    let mut lp = String::new();

    writeln!(lp, "\\ IPLCT model for CPMP-LCT").unwrap();
    writeln!(lp, "\\ S={}, H={}, P={}, T={}, tau={}", ns, nh, np, nt, tau).unwrap();
    writeln!(lp).unwrap();

    writeln!(lp, "Maximize").unwrap();
    write!(lp, " obj:").unwrap();
    for p in 1..=np {
        if p > 1 {
            write!(lp, " +").unwrap();
        }
        write!(lp, " a_{}", p).unwrap();
    }
    writeln!(lp).unwrap();
    writeln!(lp).unwrap();

    writeln!(lp, "Subject To").unwrap();

    writeln!(lp, "\\ [iplct:one-move-from]").unwrap();
    for t in 1..nt {
        write!(lp, " omf_{}: ", t).unwrap();
        let mut first = true;
        for s in 1..=ns {
            for h in 1..=nh {
                for p in 1..=np {
                    if !first {
                        write!(lp, " + ").unwrap();
                    }
                    write!(lp, "w_{}_{}_{}_{}", s, h, p, t).unwrap();
                    first = false;
                }
            }
        }
        writeln!(lp, " <= 1").unwrap();
    }

    writeln!(lp, "\\ [iplct:one-move-to]").unwrap();
    for t in 1..nt {
        write!(lp, " omt_{}: ", t).unwrap();
        let mut first = true;
        for s in 1..=ns {
            for h in 1..=nh {
                for p in 1..=np {
                    if !first {
                        write!(lp, " + ").unwrap();
                    }
                    write!(lp, "z_{}_{}_{}_{}", s, h, p, t).unwrap();
                    first = false;
                }
            }
        }
        writeln!(lp, " <= 1").unwrap();
    }

    writeln!(lp, "\\ [iplct:flow]").unwrap();
    for p in 1..=np {
        for t in 1..nt {
            write!(lp, " flow_{}_{}: ", p, t).unwrap();
            let mut first = true;
            for s in 1..=ns {
                for h in 1..=nh {
                    if !first {
                        write!(lp, " + ").unwrap();
                    }
                    write!(lp, "w_{}_{}_{}_{}", s, h, p, t).unwrap();
                    first = false;
                }
            }
            for s in 1..=ns {
                for h in 1..=nh {
                    write!(lp, " - z_{}_{}_{}_{}", s, h, p, t).unwrap();
                }
            }
            writeln!(lp, " = 0").unwrap();
        }
    }

    writeln!(lp, "\\ [iplct:consecutive-moves]").unwrap();
    if nt >= 3 {
        for t in 1..=(nt - 2) {
            write!(lp, " cm_{}: ", t).unwrap();
            let mut first = true;
            for s in 1..=ns {
                for h in 1..=nh {
                    for p in 1..=np {
                        if !first {
                            write!(lp, " + ").unwrap();
                        }
                        write!(lp, "w_{}_{}_{}_{}", s, h, p, t + 1).unwrap();
                        first = false;
                    }
                }
            }
            for s in 1..=ns {
                for h in 1..=nh {
                    for p in 1..=np {
                        write!(lp, " - w_{}_{}_{}_{}", s, h, p, t).unwrap();
                    }
                }
            }
            writeln!(lp, " <= 0").unwrap();
        }
    }

    writeln!(lp, "\\ [iplct:position]").unwrap();
    for s in 1..=ns {
        for h in 1..=nh {
            for p in 1..=np {
                for t in 1..nt {
                    writeln!(
                        lp,
                        " pos_{}_{}_{}_{}: x_{}_{}_{}_{}  + z_{}_{}_{}_{} - x_{}_{}_{}_{}  - w_{}_{}_{}_{} = 0",
                        s, h, p, t,
                        s, h, p, t,
                        s, h, p, t,
                        s, h, p, t + 1,
                        s, h, p, t,
                    )
                    .unwrap();
                }
            }
        }
    }

    writeln!(lp, "\\ [iplct:height-1]").unwrap();
    for s in 1..=ns {
        for t in 1..nt {
            write!(lp, " h1_{}_{}: ", s, t).unwrap();
            let mut first = true;
            for p in 1..=np {
                if !first {
                    write!(lp, " + ").unwrap();
                }
                write!(lp, "x_{}_{}_{}_{}", s, 1, p, t).unwrap();
                first = false;
            }
            for p in 1..=np {
                write!(lp, " + z_{}_{}_{}_{}", s, 1, p, t).unwrap();
            }
            writeln!(lp, " <= 1").unwrap();
        }
    }

    writeln!(lp, "\\ [iplct:height-h]").unwrap();
    for s in 1..=ns {
        for h in 1..nh {
            for t in 1..nt {
                write!(lp, " hh_{}_{}_{}: ", s, h, t).unwrap();
                let mut first = true;
                for p in 1..=np {
                    if !first {
                        write!(lp, " + ").unwrap();
                    }
                    write!(lp, "x_{}_{}_{}_{}", s, h + 1, p, t).unwrap();
                    first = false;
                }
                for p in 1..=np {
                    write!(lp, " + w_{}_{}_{}_{}", s, h, p, t).unwrap();
                }
                for p in 1..=np {
                    write!(lp, " + z_{}_{}_{}_{}", s, h + 1, p, t).unwrap();
                }
                for p in 1..=np {
                    write!(lp, " - x_{}_{}_{}_{}", s, h, p, t).unwrap();
                }
                writeln!(lp, " <= 0").unwrap();
            }
        }
    }

    writeln!(lp, "\\ [iplct:link-l]").unwrap();
    for s in 1..=ns {
        for k in 1..=ns {
            if k == s {
                continue;
            }
            for t in 1..nt {
                write!(lp, " ll_{}_{}_{}: ", s, k, t).unwrap();
                let mut first = true;
                for h in 1..=nh {
                    for p in 1..=np {
                        if !first {
                            write!(lp, " + ").unwrap();
                        }
                        write!(lp, "w_{}_{}_{}_{}", s, h, p, t).unwrap();
                        first = false;
                    }
                }
                for h in 1..=nh {
                    for p in 1..=np {
                        write!(lp, " + z_{}_{}_{}_{}", k, h, p, t).unwrap();
                    }
                }
                write!(lp, " - l_{}_{}_{}", s, k, t).unwrap();
                writeln!(lp, " <= 1").unwrap();
            }
        }
    }

    writeln!(lp, "\\ [iplct:link-u]").unwrap();
    for s in 1..=ns {
        for k in 1..=ns {
            if k == s {
                continue;
            }
            for t in 2..nt {
                write!(lp, " lu_{}_{}_{}: ", s, k, t).unwrap();
                let mut first = true;
                for h in 1..=nh {
                    for p in 1..=np {
                        if !first {
                            write!(lp, " + ").unwrap();
                        }
                        write!(lp, "z_{}_{}_{}_{}", s, h, p, t - 1).unwrap();
                        first = false;
                    }
                }
                for h in 1..=nh {
                    for p in 1..=np {
                        write!(lp, " + w_{}_{}_{}_{}", k, h, p, t).unwrap();
                    }
                }
                write!(lp, " - u_{}_{}_{}", s, k, t).unwrap();
                writeln!(lp, " <= 1").unwrap();
            }
        }
    }

    writeln!(lp, "\\ [iplct:strengthening]").unwrap();
    for s in 1..=ns {
        for h in 1..=nh {
            for t in 1..nt {
                write!(lp, " str_{}_{}_{}: ", s, h, t).unwrap();
                let mut first = true;
                for p in 1..=np {
                    if !first {
                        write!(lp, " + ").unwrap();
                    }
                    write!(lp, "x_{}_{}_{}_{}", s, h, p, t).unwrap();
                    first = false;
                }
                if h < nh {
                    for p in 1..=np {
                        write!(lp, " + x_{}_{}_{}_{}", s, h + 1, p, t).unwrap();
                    }
                }
                for p in 1..=np {
                    write!(lp, " + w_{}_{}_{}_{}", s, h, p, t).unwrap();
                }
                if h < nh {
                    for p in 1..=np {
                        write!(lp, " + z_{}_{}_{}_{}", s, h + 1, p, t).unwrap();
                    }
                }
                for p in 1..=np {
                    write!(lp, " + z_{}_{}_{}_{}", s, h, p, t).unwrap();
                }
                writeln!(lp, " <= 2").unwrap();
            }
        }
    }

    writeln!(lp, "\\ [iplct:avoid-transitive]").unwrap();
    if nt >= 3 {
        for s in 1..=ns {
            for t in 1..=(nt - 2) {
                write!(lp, " at_{}_{}: ", s, t).unwrap();
                let mut first = true;
                for h in 1..=nh {
                    for p in 1..=np {
                        if !first {
                            write!(lp, " + ").unwrap();
                        }
                        write!(lp, "z_{}_{}_{}_{}", s, h, p, t).unwrap();
                        first = false;
                    }
                }
                for h in 1..=nh {
                    for p in 1..=np {
                        write!(lp, " + w_{}_{}_{}_{}", s, h, p, t + 1).unwrap();
                    }
                }
                writeln!(lp, " <= 1").unwrap();
            }
        }
    }

    writeln!(lp, "\\ [iplct:avoid-same]").unwrap();
    if nt >= 3 {
        for s in 1..=ns {
            for p in 1..=np {
                for t in 1..=(nt - 2) {
                    write!(lp, " as_{}_{}_{}: ", s, p, t).unwrap();
                    let mut first = true;
                    for h in 1..=nh {
                        if !first {
                            write!(lp, " + ").unwrap();
                        }
                        write!(lp, "w_{}_{}_{}_{}", s, h, p, t).unwrap();
                        first = false;
                    }
                    for h in 1..=nh {
                        write!(lp, " + z_{}_{}_{}_{}", s, h, p, t + 1).unwrap();
                    }
                    writeln!(lp, " <= 1").unwrap();
                }
            }
        }
    }

    writeln!(lp, "\\ [iplct:lowest-end-from]").unwrap();
    if nt >= 2 {
        for s in 1..=ns {
            for h in 1..=nh {
                writeln!(lp, " lef_{}_{}: w_{}_{}_{}_{} = 0", s, h, s, h, 1, nt - 1).unwrap();
            }
        }
    }

    writeln!(lp, "\\ [iplct:lowest-end-to]").unwrap();
    if nt >= 2 {
        for s in 1..=ns {
            for h in 1..=nh {
                writeln!(lp, " let_{}_{}: z_{}_{}_{}_{} = 0", s, h, s, h, 1, nt - 1).unwrap();
            }
        }
    }

    writeln!(lp, "\\ [iplct:bottom-end]").unwrap();
    if nt >= 2 {
        for s in 1..=ns {
            for p in 2..=np {
                writeln!(lp, " be_{}_{}: w_{}_{}_{}_{} = 0", s, p, s, 1, p, nt - 1).unwrap();
            }
        }
    }

    writeln!(lp, "\\ [iplct:crane-time-limit]").unwrap();
    write!(lp, " ctl:").unwrap();
    let mut first = true;

    for s in 1..=ns {
        let ru_0_s = crane.r_u(0, s);
        for k in 1..=ns {
            if k == s {
                continue;
            }
            if !first {
                write!(lp, " +").unwrap();
            }
            write!(lp, " {:.6} l_{}_{}_{}", ru_0_s, s, k, 1).unwrap();
            first = false;
        }
    }

    for s in 1..=ns {
        for k in 1..=ns {
            if k == s {
                continue;
            }
            let rl = crane.r_l(s, k);
            for t in 1..nt {
                write!(lp, " + {:.6} l_{}_{}_{}", rl, s, k, t).unwrap();
            }
        }
    }

    for s in 1..=ns {
        for k in 1..=ns {
            if k == s {
                continue;
            }
            let ru = crane.r_u(s, k);
            for t in 2..nt {
                write!(lp, " + {:.6} u_{}_{}_{}", ru, s, k, t).unwrap();
            }
        }
    }

    for h in 1..=nh {
        let vl = crane.v_l(h);
        for s in 1..=ns {
            for p in 1..=np {
                for t in 1..nt {
                    write!(lp, " + {:.6} w_{}_{}_{}_{}", vl, s, h, p, t).unwrap();
                }
            }
        }
    }

    for h in 1..=nh {
        let vu = crane.v_u(h);
        for s in 1..=ns {
            for p in 1..=np {
                for t in 1..nt {
                    write!(lp, " + {:.6} z_{}_{}_{}_{}", vu, s, h, p, t).unwrap();
                }
            }
        }
    }

    writeln!(lp, " <= {:.6}", tau).unwrap();

    writeln!(lp, "\\ [iplct:rel-access]").unwrap();
    for s in 1..=ns {
        for h in 1..=nh {
            for p in 1..=np {
                writeln!(
                    lp,
                    " ra_{}_{}_{}: r_{}_{}_{} - x_{}_{}_{}_{} <= 0",
                    s, h, p, s, h, p, s, h, p, nt
                )
                .unwrap();
            }
        }
    }

    writeln!(lp, "\\ [iplct:rel-access-prop]").unwrap();
    for s in 1..=ns {
        for h in 1..nh {
            write!(lp, " rap_{}_{}: ", s, h).unwrap();
            let mut first = true;
            for p in 1..=np {
                if !first {
                    write!(lp, " + ").unwrap();
                }
                write!(lp, "x_{}_{}_{}_{}", s, h + 1, p, nt).unwrap();
                first = false;
            }
            for p in 1..=np {
                write!(lp, " + r_{}_{}_{}", s, h, p).unwrap();
            }
            for p in 1..=np {
                write!(lp, " - r_{}_{}_{}", s, h + 1, p).unwrap();
            }
            writeln!(lp, " <= 1").unwrap();
        }
    }

    writeln!(lp, "\\ [iplct:rel-access-block]").unwrap();
    for s in 1..=ns {
        for h in 1..nh {
            for p in 1..np {
                write!(lp, " rab_{}_{}_{}: r_{}_{}_{}", s, h, p, s, h, p).unwrap();
                for q in (p + 1)..=np {
                    write!(lp, " + x_{}_{}_{}_{}", s, h + 1, q, nt).unwrap();
                }
                writeln!(lp, " <= 1").unwrap();
            }
        }
    }

    writeln!(lp, "\\ [iplct:accessible-count]").unwrap();
    for p in 1..=np {
        write!(lp, " ac_{}: a_{}", p, p).unwrap();
        for s in 1..=ns {
            for h in 1..=nh {
                write!(lp, " - r_{}_{}_{}", s, h, p).unwrap();
            }
        }
        writeln!(lp, " <= 0").unwrap();
    }

    writeln!(lp, "\\ [iplct:group-lower]").unwrap();
    for p in 1..np {
        let cp = bay.group_count[p] as usize;
        writeln!(lp, " gl_{}: {} y_{} - a_{} <= 0", p, cp, p, p).unwrap();
    }

    writeln!(lp, "\\ [iplct:group-monotone]").unwrap();
    for p in 1..np {
        writeln!(lp, " gm_{}: y_{} - y_{} <= 0", p, p + 1, p).unwrap();
    }

    writeln!(lp, "\\ [iplct:group-upper]").unwrap();
    for p in 1..np {
        let cp1 = bay.group_count[p + 1] as usize;
        writeln!(lp, " gu_{}: a_{} - {} y_{} <= 0", p, p + 1, cp1, p).unwrap();
    }

    writeln!(lp).unwrap();

    writeln!(lp, "Bounds").unwrap();

    for s in 1..=ns {
        for h in 1..=nh {
            for p in 1..=np {
                let has = (h - 1) < bay.he(s - 1) && bay.stacks[s - 1][h - 1] == p as u8;
                let val = if has { 1 } else { 0 };
                writeln!(lp, " x_{}_{}_{}_{} = {}", s, h, p, 1, val).unwrap();
            }
        }
    }

    for p in 1..=np {
        let cp = bay.group_count[p] as usize;
        writeln!(lp, " 0 <= a_{} <= {}", p, cp).unwrap();
    }

    writeln!(lp).unwrap();

    writeln!(lp, "Binary").unwrap();

    for s in 1..=ns {
        for h in 1..=nh {
            for p in 1..=np {
                for t in 1..=nt {
                    writeln!(lp, " x_{}_{}_{}_{}", s, h, p, t).unwrap();
                }
            }
        }
    }
    for s in 1..=ns {
        for h in 1..=nh {
            for p in 1..=np {
                for t in 1..nt {
                    writeln!(lp, " w_{}_{}_{}_{}", s, h, p, t).unwrap();
                }
            }
        }
    }
    for s in 1..=ns {
        for h in 1..=nh {
            for p in 1..=np {
                for t in 1..nt {
                    writeln!(lp, " z_{}_{}_{}_{}", s, h, p, t).unwrap();
                }
            }
        }
    }
    for s in 1..=ns {
        for k in 1..=ns {
            if k == s {
                continue;
            }
            for t in 1..nt {
                writeln!(lp, " l_{}_{}_{}", s, k, t).unwrap();
            }
        }
    }
    for s in 1..=ns {
        for k in 1..=ns {
            if k == s {
                continue;
            }
            for t in 2..nt {
                writeln!(lp, " u_{}_{}_{}", s, k, t).unwrap();
            }
        }
    }
    for s in 1..=ns {
        for h in 1..=nh {
            for p in 1..=np {
                writeln!(lp, " r_{}_{}_{}", s, h, p).unwrap();
            }
        }
    }
    for p in 1..=np {
        writeln!(lp, " y_{}", p).unwrap();
    }

    writeln!(lp).unwrap();

    writeln!(lp, "General").unwrap();
    for p in 1..=np {
        writeln!(lp, " a_{}", p).unwrap();
    }
    writeln!(lp).unwrap();

    writeln!(lp, "End").unwrap();

    fs::write(output_path, &lp).expect("cannot write .lp file");
    output_path.to_string()
}

fn compute_t_max(bay: &Bay, crane: &CraneParams, tau: f64) -> usize {
    let t_min_h = crane.t_min(bay.h);
    (tau / t_min_h).floor() as usize + 1
}

fn count_rel_accessible_group(bay: &Bay, group: usize) -> usize {
    let mut count = 0usize;
    for s in 0..bay.s {
        for t in 0..bay.he(s) {
            if bay.stacks[s][t] as usize == group && !bay.is_blocked(s, t) {
                count += 1;
            }
        }
    }
    count
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

fn is_model_accessible(p_star: usize, group: usize, rel_accessible: bool, p_total: usize) -> bool {
    if p_star > p_total {
        return true;
    }
    group < p_star || (group == p_star && rel_accessible)
}

fn write_gurobi_warm_start_sol(
    bay: &Bay,
    crane: &CraneParams,
    tau: f64,
    lp_path: &str,
    heuristic: WarmStartHeuristic,
) -> String {
    let nt = compute_t_max(bay, crane, tau);
    let warm_sol = run_warm_start_heuristic(bay, crane, tau, heuristic);
    assert!(
        warm_sol.len() <= nt.saturating_sub(1),
        "warm-start solution has {} moves but IPLCT allows at most {} moves",
        warm_sol.len(),
        nt.saturating_sub(1)
    );

    let mut bay_final = *bay;
    for dm in &warm_sol.moves {
        let applied = bay_final.apply_move(Move {
            src: dm.src,
            dst: dm.dst,
        });
        assert_eq!(applied.group, dm.group, "warm-start move group mismatch");
    }

    let p_star = bay_final.p_star();
    let mut a_vals = vec![0usize; bay_final.p + 1];
    if p_star > bay_final.p {
        for (p, val) in a_vals.iter_mut().enumerate().take(bay_final.p + 1).skip(1) {
            *val = bay_final.group_count[p] as usize;
        }
    } else {
        for (p, val) in a_vals.iter_mut().enumerate().take(p_star).skip(1) {
            *val = bay_final.group_count[p] as usize;
        }
        a_vals[p_star] = count_rel_accessible_group(&bay_final, p_star);
    }

    let mut vars = BTreeMap::<String, i64>::new();

    for (idx, dm) in warm_sol.moves.iter().enumerate() {
        let t = idx + 1;
        let src = dm.src + 1;
        let dst = dm.dst + 1;
        let src_tier = dm.src_tier + 1;
        let dst_tier = dm.dst_tier + 1;
        let group = dm.group as usize;

        vars.insert(format!("w_{}_{}_{}_{}", src, src_tier, group, t), 1);
        vars.insert(format!("z_{}_{}_{}_{}", dst, dst_tier, group, t), 1);
        vars.insert(format!("l_{}_{}_{}", src, dst, t), 1);

        if t >= 2 {
            let prev_dst = warm_sol.moves[idx - 1].dst + 1;
            vars.insert(format!("u_{}_{}_{}", prev_dst, src, t), 1);
        }
    }

    for s in 0..bay_final.s {
        for h in 0..bay_final.he(s) {
            let g = bay_final.stacks[s][h] as usize;
            if g == 0 {
                continue;
            }
            let rel_accessible = !bay_final.is_blocked(s, h);
            if is_model_accessible(p_star, g, rel_accessible, bay_final.p) {
                vars.insert(format!("r_{}_{}_{}", s + 1, h + 1, g), 1);
            }
        }
    }

    for (p, &a) in a_vals.iter().enumerate().take(bay_final.p + 1).skip(1) {
        if a > 0 {
            vars.insert(format!("a_{}", p), a as i64);
        }
    }

    if bay_final.p >= 2 {
        for (p, &a) in a_vals.iter().enumerate().take(bay_final.p).skip(1) {
            if a == bay_final.group_count[p] as usize {
                vars.insert(format!("y_{}", p), 1);
            }
        }
    }

    let warm_path = format!(
        "{}.warmstart.sol",
        lp_path.strip_suffix(".lp").unwrap_or(lp_path)
    );
    let mut out = String::new();
    writeln!(out, "# Solution for model").unwrap();
    writeln!(out, "# Warm start heuristic = {:?}", heuristic).unwrap();
    for (name, value) in vars {
        writeln!(out, "{} {}", name, value).unwrap();
    }
    fs::write(&warm_path, out).expect("cannot write Gurobi warm-start file");
    warm_path
}

/// Generates the LP and solves it with Gurobi.
pub fn generate_and_solve(
    bay: &Bay,
    crane: &CraneParams,
    tau: f64,
    lp_path: &str,
    time_limit: f64,
) -> IplctResult {
    generate_and_solve_options(
        bay,
        crane,
        tau,
        lp_path,
        time_limit,
        IplctSolveOptions::default(),
    )
}

/// Generates the LP, solves it with Gurobi, and reconstructs the move sequence.
pub fn generate_and_solve_with_solution(
    bay: &Bay,
    crane: &CraneParams,
    tau: f64,
    lp_path: &str,
    time_limit: f64,
) -> IplctResult {
    generate_and_solve_with_solution_options(
        bay,
        crane,
        tau,
        lp_path,
        time_limit,
        IplctSolveOptions::default(),
    )
}

/// Generates and solves the LP using the provided solve options.
pub fn generate_and_solve_options(
    bay: &Bay,
    crane: &CraneParams,
    tau: f64,
    lp_path: &str,
    time_limit: f64,
    options: IplctSolveOptions,
) -> IplctResult {
    generate_and_solve_with_solution_options(bay, crane, tau, lp_path, time_limit, options)
}

/// Generates, solves, and reconstructs the LP solution using the provided options.
pub fn generate_and_solve_with_solution_options(
    bay: &Bay,
    crane: &CraneParams,
    tau: f64,
    lp_path: &str,
    time_limit: f64,
    options: IplctSolveOptions,
) -> IplctResult {
    generate_lp(bay, crane, tau, lp_path);
    let warm_start_path = options
        .warm_start
        .map(|heuristic| write_gurobi_warm_start_sol(bay, crane, tau, lp_path, heuristic));
    solve_gurobi_with_solution_impl(bay, crane, lp_path, time_limit, warm_start_path.as_deref())
}

/// Solves an already generated LP with Gurobi and returns only the objective value.
pub fn solve_gurobi(lp_path: &str, time_limit: f64) -> IplctResult {
    let (log, sol) = run_gurobi(lp_path, time_limit, None);
    let acc = parse_objective_from_sol(&sol).unwrap_or(0);
    let mut solution = Solution::new();
    solution.acc = acc;
    IplctResult {
        acc,
        optimal: log.contains("Optimal solution found"),
        time_seconds: parse_time_seconds(&log),
        solution,
    }
}

/// Solves an already generated LP with Gurobi and reconstructs the move sequence.
pub fn solve_gurobi_with_solution(
    bay: &Bay,
    crane: &CraneParams,
    lp_path: &str,
    time_limit: f64,
) -> IplctResult {
    solve_gurobi_with_solution_impl(bay, crane, lp_path, time_limit, None)
}

fn solve_gurobi_with_solution_impl(
    bay: &Bay,
    crane: &CraneParams,
    lp_path: &str,
    time_limit: f64,
    warm_start_path: Option<&str>,
) -> IplctResult {
    let (log, sol_text) = run_gurobi(lp_path, time_limit, warm_start_path);
    let optimal = log.contains("Optimal solution found");
    let time_seconds = parse_time_seconds(&log);
    let acc_obj = parse_objective_from_sol(&sol_text).unwrap_or(0);

    let solution = reconstruct_solution_from_sol(bay, crane, &sol_text);
    if acc_obj != 0 && solution.acc != acc_obj {
        panic!(
            "IPLCT solution inconsistency: objective acc={} but reconstructed acc={}",
            acc_obj, solution.acc
        );
    }

    IplctResult {
        acc: solution.acc,
        optimal,
        time_seconds,
        solution,
    }
}

fn run_gurobi(lp_path: &str, time_limit: f64, warm_start_path: Option<&str>) -> (String, String) {
    use std::process::Command;

    let sol_path = format!("{}.sol", lp_path.strip_suffix(".lp").unwrap_or(lp_path));
    let log_path = format!("{}.log", lp_path.strip_suffix(".lp").unwrap_or(lp_path));
    let threads = std::env::var("GUROBI_THREADS")
        .ok()
        .and_then(|v| v.parse::<usize>().ok())
        .filter(|&v| v >= 1)
        .or_else(|| {
            std::env::var("IPLCT_GUROBI_THREADS")
                .ok()
                .and_then(|v| v.parse::<usize>().ok())
                .filter(|&v| v >= 1)
        })
        .or_else(|| {
            std::env::var("CP_OPTIMIZER_WORKERS")
                .ok()
                .and_then(|v| v.parse::<usize>().ok())
                .filter(|&v| v >= 1)
        })
        .or_else(|| {
            std::env::var("SLURM_CPUS_PER_TASK")
                .ok()
                .and_then(|v| v.parse::<usize>().ok())
                .filter(|&v| v >= 1)
        })
        .unwrap_or(1);

    let mut args = vec![
        format!("ResultFile={}", sol_path),
        format!("LogFile={}", log_path),
        format!("Threads={}", threads),
        "MemLimit=16".to_string(),
    ];
    if time_limit.is_finite() {
        args.push(format!("TimeLimit={}", time_limit));
    }
    if let Some(path) = warm_start_path {
        args.push(format!("InputFile={}", path));
    }
    args.push(lp_path.to_string());

    let output = Command::new("gurobi_cl")
        .args(&args)
        .output()
        .expect("failed to run gurobi_cl — is Gurobi installed and on PATH?");

    let mut log = String::from_utf8_lossy(&output.stdout).to_string();
    if !output.stderr.is_empty() {
        log.push('\n');
        log.push_str(&String::from_utf8_lossy(&output.stderr));
    }
    let sol_text = fs::read_to_string(sol_path).unwrap_or_default();
    (log, sol_text)
}

fn parse_time_seconds(log: &str) -> f64 {
    log.lines()
        .find(|l| l.contains("Explored") && l.contains("seconds"))
        .and_then(|l| {
            let words: Vec<&str> = l.split_whitespace().collect();
            words.iter().position(|&w| w == "seconds").and_then(|i| {
                if i > 0 {
                    words[i - 1].parse::<f64>().ok()
                } else {
                    None
                }
            })
        })
        .unwrap_or(0.0)
}

fn parse_objective_from_sol(sol_text: &str) -> Option<usize> {
    sol_text
        .lines()
        .find(|l| l.starts_with("# Objective value"))
        .and_then(|l| {
            l.split('=')
                .nth(1)
                .and_then(|v| v.trim().parse::<f64>().ok())
                .map(|v| v.round() as usize)
        })
}

fn parse_move_vars(sol_text: &str, prefix: &str) -> BTreeMap<usize, (usize, usize, usize)> {
    let mut out = BTreeMap::new();
    for line in sol_text.lines() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let mut parts = line.split_whitespace();
        let name = parts.next().unwrap_or("");
        let value = parts
            .next()
            .and_then(|v| v.parse::<f64>().ok())
            .unwrap_or(0.0);
        if value < 0.5 || !name.starts_with(prefix) {
            continue;
        }
        let tokens: Vec<&str> = name.split('_').collect();
        if tokens.len() != 5 {
            continue;
        }
        let s = tokens[1].parse::<usize>().ok();
        let h = tokens[2].parse::<usize>().ok();
        let p = tokens[3].parse::<usize>().ok();
        let t = tokens[4].parse::<usize>().ok();
        if let (Some(s), Some(h), Some(p), Some(t)) = (s, h, p, t) {
            out.insert(t, (s, h, p));
        }
    }
    out
}

fn reconstruct_solution_from_sol(bay: &Bay, crane: &CraneParams, sol_text: &str) -> Solution {
    let ws = parse_move_vars(sol_text, "w_");
    let zs = parse_move_vars(sol_text, "z_");

    let mut bay_cur = *bay;
    let mut solution = Solution::new();
    let mut prev_dst_1based = 0usize;

    let max_t = ws.keys().chain(zs.keys()).copied().max().unwrap_or(0usize);

    for t in 1..=max_t {
        let w = ws.get(&t).copied();
        let z = zs.get(&t).copied();
        match (w, z) {
            (None, None) => continue,
            (Some(_), None) | (None, Some(_)) => {
                panic!("invalid IPLCT .sol: incomplete move at t={}", t);
            }
            (Some((ws, wh, wp)), Some((zs, zh, zp))) => {
                if ws == zs {
                    continue;
                }

                let src = ws.saturating_sub(1);
                let dst = zs.saturating_sub(1);
                if bay_cur.is_empty(src) {
                    panic!("invalid IPLCT .sol: source stack {} empty at t={}", ws, t);
                }
                if bay_cur.is_full(dst) {
                    panic!(
                        "invalid IPLCT .sol: destination stack {} full at t={}",
                        zs, t
                    );
                }
                let dm = bay_cur.apply_move(Move { src, dst });
                if dm.group as usize != wp || dm.group as usize != zp {
                    panic!(
                        "invalid IPLCT .sol: moved group {} but vars encode {} -> {} at t={}",
                        dm.group, wp, zp, t
                    );
                }
                if dm.src_tier + 1 != wh || dm.dst_tier + 1 != zh {
                    panic!(
                        "invalid IPLCT .sol: tier mismatch at t={} (src {}!={}, dst {}!={})",
                        t,
                        dm.src_tier + 1,
                        wh,
                        dm.dst_tier + 1,
                        zh
                    );
                }

                let delta = crane.move_time(
                    src + 1,
                    dm.src_tier + 1,
                    dst + 1,
                    dm.dst_tier + 1,
                    prev_dst_1based,
                );
                solution.push(dm, delta);
                prev_dst_1based = dst + 1;
            }
        }
    }

    solution.acc = bay_cur.acc();
    solution
}
