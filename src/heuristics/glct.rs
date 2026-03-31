//! Greedy heuristic for the CPMP-LCT.
//!
//! GLCT is an accessibility-oriented constructive heuristic.
//! It focuses on exposing low-priority blocked containers while
//! respecting the remaining crane-time budget.

use crate::bay::{Bay, Group, Move, MAX_STACKS};
use crate::crane::CraneParams;
use crate::solution::Solution;

/// Runs GLCT in place on `bay` with the remaining crane-time budget `tau_remaining`.
pub fn glct(bay: &mut Bay, crane: &CraneParams, tau_remaining: f64) -> Solution {
    let mut sol = Solution::new();
    let mut remaining = tau_remaining;

    while remaining > 0.0 {
        let Some(p) = bay.min_blocked_group() else {
            break;
        };

        let mut bl_stacks = bay.stacks_with_blocked(p);
        while !bl_stacks.is_empty() {
            let (st, n_blocking) = select_stack_agh(bay, p, &bl_stacks);

            let mut success = true;
            for _ in 0..n_blocking {
                if !try_relocate_agh(bay, crane, &mut sol, p, st, &mut remaining) {
                    remove_stack(&mut bl_stacks, st);
                    success = false;
                    break;
                }
            }

            if success && !stack_has_blocked_group(bay, st, p) {
                remove_stack(&mut bl_stacks, st);
            }
        }

        bl_stacks = bay.stacks_with_blocked(p);
        if !bl_stacks.is_empty() {
            let st = select_stack_agh(bay, p, &bl_stacks).0;
            if !try_reorder_agh(bay, crane, &mut sol, p, st, &mut remaining) {
                break;
            }
        }
    }

    sol.acc = bay.acc();
    sol
}

fn select_stack_agh(bay: &Bay, p: Group, bl_stacks: &[usize]) -> (usize, usize) {
    let mut best_s = bl_stacks[0];
    let mut best_n = usize::MAX;
    let mut best_min_gv = 0usize;
    let mut best_h = 0usize;
    let mut best_dist = usize::MAX;

    for &s in bl_stacks {
        let c_tier = find_topmost_blocked(bay, s, p).expect("blocked stack must contain blocked p");
        let t = find_lowest_blocker(bay, s, c_tier, p);
        let n_blocking = bay.he(s) - t;

        let min_gv_blocking = (t..bay.he(s))
            .map(|h| bay.stacks[s][h] as usize)
            .min()
            .unwrap_or(usize::MAX);
        let height = bay.he(s);
        let dist = distance_to_center_twice(s, bay.s);

        let better = n_blocking < best_n
            || (n_blocking == best_n && min_gv_blocking > best_min_gv)
            || (n_blocking == best_n && min_gv_blocking == best_min_gv && height > best_h)
            || (n_blocking == best_n
                && min_gv_blocking == best_min_gv
                && height == best_h
                && dist < best_dist)
            || (n_blocking == best_n
                && min_gv_blocking == best_min_gv
                && height == best_h
                && dist == best_dist
                && s < best_s);

        if better {
            best_s = s;
            best_n = n_blocking;
            best_min_gv = min_gv_blocking;
            best_h = height;
            best_dist = dist;
        }
    }

    (best_s, best_n)
}

fn try_relocate_agh(
    bay: &mut Bay,
    crane: &CraneParams,
    sol: &mut Solution,
    p: Group,
    st: usize,
    remaining: &mut f64,
) -> bool {
    if bay.is_empty(st) {
        return false;
    }

    let moving_gv = bay.top_group(st);

    let mut phi = Vec::with_capacity(bay.s);
    for dst in 0..bay.s {
        if dst == st || bay.is_full(dst) {
            continue;
        }
        if estimate_move_time(bay, crane, sol, st, dst) > *remaining {
            continue;
        }
        if !will_block_agh(bay, p, dst, moving_gv, bay.he(dst)) {
            phi.push(dst);
        }
    }

    if !phi.is_empty() {
        let to = select_phi_dst_agh(bay, crane, sol, st, &phi);
        return do_move_with_budget(bay, crane, sol, remaining, st, to);
    }

    let mut psi = Vec::with_capacity(bay.s);
    for k in 0..bay.s {
        if k == st || bay.is_empty(k) {
            continue;
        }

        let ctime_sk = estimate_move_time_after_pop(bay, crane, sol, st, k);
        if ctime_sk > *remaining {
            continue;
        }

        let top_after_pop = bay.he(k).saturating_sub(1);
        if will_block_agh(bay, p, k, moving_gv, top_after_pop) {
            continue;
        }

        let aux_gv = bay.top_group(k);
        let mut can_be_moved = false;
        for q in 0..bay.s {
            if q == st || q == k || bay.is_full(q) {
                continue;
            }
            let t_kq = estimate_move_time(bay, crane, sol, k, q);
            if ctime_sk + t_kq > *remaining {
                continue;
            }
            if !will_block_agh(bay, p, q, aux_gv, bay.he(q)) {
                can_be_moved = true;
                break;
            }
        }

        if can_be_moved {
            psi.push(k);
        }
    }

    if psi.is_empty() {
        return false;
    }

    let to = match select_psi_stack_agh(bay, crane, sol, p, st, &psi) {
        Some(v) => v,
        None => return false,
    };
    let ctime_sk = estimate_move_time_after_pop(bay, crane, sol, st, to);
    let aux_gv = bay.top_group(to);

    let mut psi2 = Vec::with_capacity(bay.s);
    for q in 0..bay.s {
        if q == st || q == to || bay.is_full(q) {
            continue;
        }
        let t_toq = estimate_move_time(bay, crane, sol, to, q);
        if ctime_sk + t_toq > *remaining {
            continue;
        }
        if !will_block_agh(bay, p, q, aux_gv, bay.he(q)) {
            psi2.push(q);
        }
    }

    let Some(to2) = select_psi2_dst_agh(bay, crane, sol, to, &psi2) else {
        return false;
    };

    let checkpoint = sol.len();
    if !do_move_with_budget(bay, crane, sol, remaining, to, to2) {
        return false;
    }
    if !do_move_with_budget(bay, crane, sol, remaining, st, to) {
        rollback_to_len(bay, sol, remaining, checkpoint);
        return false;
    }

    true
}

fn select_phi_dst_agh(
    bay: &Bay,
    crane: &CraneParams,
    sol: &Solution,
    st: usize,
    phi: &[usize],
) -> usize {
    let mut ordered = Vec::with_capacity(phi.len());
    for &s in phi {
        if !stack_has_blocked(bay, s) {
            ordered.push(s);
        }
    }
    let mut cand = if ordered.is_empty() {
        phi.to_vec()
    } else {
        ordered
    };

    if cand.len() > 1 {
        let mut min_gv = [0usize; MAX_STACKS];
        let mut max_min_gv = 0usize;
        for &s in &cand {
            let v = min_group_full_or_sentinel(bay, s);
            min_gv[s] = v;
            if v > max_min_gv {
                max_min_gv = v;
            }
        }
        cand.retain(|&s| min_gv[s] == max_min_gv);
    }

    if cand.len() > 1 {
        let mut ctime = [0.0f64; MAX_STACKS];
        let mut min_ct = f64::INFINITY;
        for &s in &cand {
            let t = estimate_move_time(bay, crane, sol, st, s);
            ctime[s] = t;
            if t < min_ct {
                min_ct = t;
            }
        }
        cand.retain(|&s| ctime[s] == min_ct);
    }

    if cand.len() > 1 {
        let mut heights = [0usize; MAX_STACKS];
        let mut min_he = usize::MAX;
        for &s in &cand {
            let he = bay.he(s);
            heights[s] = he;
            if he < min_he {
                min_he = he;
            }
        }
        cand.retain(|&s| heights[s] == min_he);
    }

    cand.sort_unstable();
    cand[0]
}

fn select_psi_stack_agh(
    bay: &Bay,
    crane: &CraneParams,
    sol: &Solution,
    p: Group,
    st: usize,
    psi: &[usize],
) -> Option<usize> {
    let mut cand = psi.to_vec();

    if cand.len() > 1 {
        let mut ordered_after_pop = Vec::with_capacity(cand.len());
        for &s in &cand {
            let top_after_pop = bay.he(s).saturating_sub(1);
            if top_after_pop == 0 || n_clean_in_prefix(bay, s, top_after_pop) >= top_after_pop {
                ordered_after_pop.push(s);
            }
        }
        if !ordered_after_pop.is_empty() {
            cand = ordered_after_pop;
        }
    }

    if cand.len() > 1 {
        let mut min_gv_after_pop = [0usize; MAX_STACKS];
        let mut max_min_gv = 0usize;
        for &s in &cand {
            let top_after_pop = bay.he(s).saturating_sub(1);
            let v = min_group_in_prefix_or_sentinel(bay, s, top_after_pop);
            min_gv_after_pop[s] = v;
            if v > max_min_gv {
                max_min_gv = v;
            }
        }
        cand.retain(|&s| min_gv_after_pop[s] == max_min_gv);
    }

    if cand.len() > 1 {
        let mut top_gv = [0u8; MAX_STACKS];
        let mut max_top_gv = 0u8;
        for &s in &cand {
            let gv = bay.top_group(s);
            top_gv[s] = gv;
            if gv > max_top_gv {
                max_top_gv = gv;
            }
        }
        cand.retain(|&s| top_gv[s] == max_top_gv);
    }

    if cand.len() > 1 {
        let mut ctime_after_pop = [0.0f64; MAX_STACKS];
        let mut min_ct = f64::INFINITY;
        for &s in &cand {
            let t = estimate_move_time_after_pop(bay, crane, sol, st, s);
            ctime_after_pop[s] = t;
            if t < min_ct {
                min_ct = t;
            }
        }
        cand.retain(|&s| ctime_after_pop[s] == min_ct);
    }

    if cand.len() > 1 {
        let mut dist_to_center = [0usize; MAX_STACKS];
        let mut min_dist = usize::MAX;
        for &s in &cand {
            let d = distance_to_center_twice(s, bay.s);
            dist_to_center[s] = d;
            if d < min_dist {
                min_dist = d;
            }
        }
        cand.retain(|&s| dist_to_center[s] == min_dist);
    }

    cand.sort_unstable();
    cand.first().copied().filter(|_| p > 0)
}

fn select_psi2_dst_agh(
    bay: &Bay,
    crane: &CraneParams,
    sol: &Solution,
    to: usize,
    psi2: &[usize],
) -> Option<usize> {
    if psi2.is_empty() {
        return None;
    }

    let mut ordered = Vec::with_capacity(psi2.len());
    for &s in psi2 {
        if !stack_has_blocked(bay, s) {
            ordered.push(s);
        }
    }
    let mut cand = if ordered.is_empty() {
        psi2.to_vec()
    } else {
        ordered
    };
    if cand.is_empty() {
        return None;
    }

    if cand.len() > 1 {
        let mut min_gv = [0usize; MAX_STACKS];
        let mut max_min_gv = 0usize;
        for &s in &cand {
            let v = min_group_full_or_sentinel(bay, s);
            min_gv[s] = v;
            if v > max_min_gv {
                max_min_gv = v;
            }
        }
        cand.retain(|&s| min_gv[s] == max_min_gv);
    }

    if cand.len() > 1 {
        let mut ctime = [0.0f64; MAX_STACKS];
        let mut min_ct = f64::INFINITY;
        for &s in &cand {
            let t = estimate_move_time(bay, crane, sol, to, s);
            ctime[s] = t;
            if t < min_ct {
                min_ct = t;
            }
        }
        cand.retain(|&s| ctime[s] == min_ct);
    }

    if cand.len() > 1 {
        let mut heights = [0usize; MAX_STACKS];
        let mut min_he = usize::MAX;
        for &s in &cand {
            let he = bay.he(s);
            heights[s] = he;
            if he < min_he {
                min_he = he;
            }
        }
        cand.retain(|&s| heights[s] == min_he);
    }

    cand.sort_unstable();
    cand.first().copied()
}

fn try_reorder_agh(
    bay: &mut Bay,
    crane: &CraneParams,
    sol: &mut Solution,
    p: Group,
    st: usize,
    remaining: &mut f64,
) -> bool {
    let checkpoint = sol.len();

    let mut has_non_full = false;
    let mut max_he = 0usize;
    for s in 0..bay.s {
        if s == st || bay.is_full(s) {
            continue;
        }
        has_non_full = true;
        let he = bay.he(s);
        if he > max_he {
            max_he = he;
        }
    }
    if !has_non_full {
        return false;
    }

    let mut res = usize::MAX;
    let mut best_dist = usize::MAX;
    for s in 0..bay.s {
        if s == st || bay.is_full(s) || bay.he(s) != max_he {
            continue;
        }
        let dist = distance_between_stacks(s, st);
        if dist < best_dist || (dist == best_dist && s < res) {
            best_dist = dist;
            res = s;
        }
    }
    if res == usize::MAX {
        return false;
    }

    let mut moved_to = [0usize; MAX_STACKS];

    for _ in 0..bay.h {
        let mut to = usize::MAX;
        let mut min_ct = f64::INFINITY;
        for s in 0..bay.s {
            if s == st || s == res || bay.is_full(s) {
                continue;
            }
            let ct = estimate_move_time(bay, crane, sol, st, s);
            if ct < min_ct || (ct == min_ct && s < to) {
                min_ct = ct;
                to = s;
            }
        }
        if to == usize::MAX {
            rollback_to_len(bay, sol, remaining, checkpoint);
            return false;
        }

        if !do_move_with_budget(bay, crane, sol, remaining, st, to) {
            rollback_to_len(bay, sol, remaining, checkpoint);
            return false;
        }
        moved_to[to] += 1;

        if bay.is_empty(st) {
            rollback_to_len(bay, sol, remaining, checkpoint);
            return false;
        }
        if bay.top_group(st) == p {
            break;
        }
    }

    if !do_move_with_budget(bay, crane, sol, remaining, st, res) {
        rollback_to_len(bay, sol, remaining, checkpoint);
        return false;
    }

    let mut from_st: Vec<usize> = (0..bay.s).filter(|&s| moved_to[s] > 0).collect();
    if from_st.is_empty() {
        rollback_to_len(bay, sol, remaining, checkpoint);
        return false;
    }

    let min_dist = from_st
        .iter()
        .map(|&s| distance_between_stacks(s, st))
        .min()
        .unwrap();
    from_st.retain(|&s| distance_between_stacks(s, st) == min_dist);
    from_st.sort_unstable();
    let from = from_st[0];

    if !do_move_with_budget(bay, crane, sol, remaining, from, st) {
        rollback_to_len(bay, sol, remaining, checkpoint);
        return false;
    }
    moved_to[from] -= 1;

    for s in 0..bay.s {
        for _ in 0..moved_to[s] {
            if !do_move_with_budget(bay, crane, sol, remaining, s, st) {
                rollback_to_len(bay, sol, remaining, checkpoint);
                return false;
            }
        }
    }

    let is_blocking = (0..bay.he(res).saturating_sub(1)).any(|h| bay.stacks[res][h] < p);
    if is_blocking && !do_move_with_budget(bay, crane, sol, remaining, res, st) {
        rollback_to_len(bay, sol, remaining, checkpoint);
        return false;
    }

    true
}

fn find_topmost_blocked(bay: &Bay, s: usize, p: Group) -> Option<usize> {
    for t in (0..bay.he(s)).rev() {
        if bay.stacks[s][t] == p && bay.is_blocked(s, t) {
            return Some(t);
        }
    }
    None
}

fn find_lowest_blocker(bay: &Bay, s: usize, c_tier: usize, p: Group) -> usize {
    for t in (c_tier + 1)..bay.he(s) {
        if bay.stacks[s][t] > p {
            return t;
        }
    }
    bay.he(s)
}

fn will_block_agh(bay: &Bay, p: Group, dst: usize, moving_gv: Group, top_exclusive: usize) -> bool {
    if top_exclusive == 0 {
        return false;
    }
    let start = top_clean_start_in_prefix(bay, dst, top_exclusive);
    for h in start..top_exclusive {
        let gv = bay.stacks[dst][h];
        if gv <= p && gv < moving_gv {
            return true;
        }
    }
    false
}

fn top_clean_start_in_prefix(bay: &Bay, s: usize, top_exclusive: usize) -> usize {
    let mut t = top_exclusive;
    while t > 0 && !is_blocked_in_prefix(bay, s, t - 1, top_exclusive) {
        t -= 1;
    }
    t
}

fn is_blocked_in_prefix(bay: &Bay, s: usize, t: usize, top_exclusive: usize) -> bool {
    let g = bay.stacks[s][t];
    for h in (t + 1)..top_exclusive {
        if bay.stacks[s][h] > g {
            return true;
        }
    }
    false
}

fn n_clean_in_prefix(bay: &Bay, s: usize, top_exclusive: usize) -> usize {
    if top_exclusive == 0 {
        0
    } else {
        top_exclusive - top_clean_start_in_prefix(bay, s, top_exclusive)
    }
}

fn min_group_in_prefix_or_sentinel(bay: &Bay, s: usize, top_exclusive: usize) -> usize {
    if top_exclusive == 0 {
        bay.p + 1
    } else {
        (0..top_exclusive)
            .map(|h| bay.stacks[s][h] as usize)
            .min()
            .unwrap_or(bay.p + 1)
    }
}

fn min_group_full_or_sentinel(bay: &Bay, s: usize) -> usize {
    min_group_in_prefix_or_sentinel(bay, s, bay.he(s))
}

fn stack_has_blocked_group(bay: &Bay, s: usize, p: Group) -> bool {
    (0..bay.he(s)).any(|h| bay.stacks[s][h] == p && bay.is_blocked(s, h))
}

fn stack_has_blocked(bay: &Bay, s: usize) -> bool {
    (0..bay.he(s)).any(|h| bay.is_blocked(s, h))
}

fn remove_stack(v: &mut Vec<usize>, s: usize) {
    if let Some(i) = v.iter().position(|&x| x == s) {
        v.remove(i);
    }
}

fn estimate_move_time_after_pop(
    bay: &Bay,
    crane: &CraneParams,
    sol: &Solution,
    src: usize,
    dst_after_pop: usize,
) -> f64 {
    let src_tier = bay.he(src);
    let dst_tier_after_pop = bay.he(dst_after_pop);
    crane.move_time(
        src + 1,
        src_tier,
        dst_after_pop + 1,
        dst_tier_after_pop,
        sol.last_dst_1based(),
    )
}

fn estimate_move_time(
    bay: &Bay,
    crane: &CraneParams,
    sol: &Solution,
    src: usize,
    dst: usize,
) -> f64 {
    let src_tier = bay.he(src);
    let dst_tier = bay.he(dst) + 1;
    crane.move_time(src + 1, src_tier, dst + 1, dst_tier, sol.last_dst_1based())
}

fn distance_to_center_twice(stack: usize, s: usize) -> usize {
    let center2 = (s as isize) - 1;
    ((2 * stack as isize) - center2).unsigned_abs()
}

fn distance_between_stacks(a: usize, b: usize) -> usize {
    a.abs_diff(b)
}

fn do_move_with_budget(
    bay: &mut Bay,
    crane: &CraneParams,
    sol: &mut Solution,
    remaining: &mut f64,
    from: usize,
    to: usize,
) -> bool {
    let time = estimate_move_time(bay, crane, sol, from, to);
    if time > *remaining {
        return false;
    }
    do_move(bay, crane, sol, from, to);
    *remaining -= time;
    true
}

fn rollback_to_len(bay: &mut Bay, sol: &mut Solution, remaining: &mut f64, len: usize) {
    while sol.len() > len {
        if let Some((dm, dt)) = sol.pop() {
            bay.undo_move(dm);
            *remaining += dt;
        }
    }
}

fn do_move(bay: &mut Bay, crane: &CraneParams, sol: &mut Solution, from: usize, to: usize) {
    let dm = bay.apply_move(Move { src: from, dst: to });
    let time = crane.move_time(
        dm.src + 1,
        dm.src_tier + 1,
        dm.dst + 1,
        dm.dst_tier + 1,
        sol.last_dst_1based(),
    );
    sol.push(dm, time);
}
