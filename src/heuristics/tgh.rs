//! Target-guided heuristic (TGH) for the classical CPMP.
//!
//! This module implements the constructive heuristic described in
//! Wang et al. (2015). The algorithm fixes containers in
//! descending group order, repeatedly selecting a target container / target
//! stack pair and executing a "giant move" that preserves already fixed
//! containers. The deterministic variant is used directly by HLCT; the
//! randomized selector is used to diversify target-pair choices.

use rand::Rng;

use crate::bay::{Bay, Group, Move, MAX_STACKS};
use crate::crane::CraneParams;
use crate::solution::Solution;

/// Mutable state carried by the TGH construction.
///
/// `fixed_height[s]` is the number of bottom tiers of stack `s` that have
/// already been fixed and should stay undisturbed except for the paper's
/// special temporary-displacement scenarios.
struct TghState {
    fixed_height: [usize; MAX_STACKS],
}

impl TghState {
    fn new() -> Self {
        TghState {
            fixed_height: [0; MAX_STACKS],
        }
    }

    fn uf(&self, bay: &Bay, s: usize) -> usize {
        bay.he(s) - self.fixed_height[s]
    }

    fn is_fixed(&self, s: usize, t: usize) -> bool {
        t < self.fixed_height[s]
    }

    fn fix(&mut self, s: usize) {
        self.fixed_height[s] += 1;
    }

    fn unfix(&mut self, s: usize) {
        debug_assert!(self.fixed_height[s] > 0);
        self.fixed_height[s] -= 1;
    }
}

fn is_clean(bay: &Bay, state: &TghState, s: usize) -> bool {
    let base = state.fixed_height[s];
    let top = bay.he(s);
    for t in (base + 1)..top {
        if bay.stacks[s][t] > bay.stacks[s][t - 1] {
            return false;
        }
    }
    true
}

fn is_dirty(bay: &Bay, state: &TghState, s: usize) -> bool {
    !is_clean(bay, state, s)
}

fn ldg(bay: &Bay, state: &TghState, s: usize) -> Group {
    let base = state.fixed_height[s];
    let top = bay.he(s);
    let mut max_dirty: Group = 0;
    let mut min_below: Group = Group::MAX;
    for t in base..top {
        let g = bay.stacks[s][t];
        if g > min_below {
            if g > max_dirty {
                max_dirty = g;
            }
        }
        if g < min_below {
            min_below = g;
        }
    }
    max_dirty
}

fn can_accommodate(bay: &Bay, state: &TghState, s: usize, g: Group) -> bool {
    if bay.is_full(s) {
        return false;
    }
    if !is_clean(bay, state, s) {
        return false;
    }
    if bay.is_empty(s) {
        return true;
    }
    bay.top_group(s) >= g
}

fn is_container_dirty(bay: &Bay, state: &TghState, s: usize, t: usize) -> bool {
    let base = state.fixed_height[s];
    if t < base {
        return false;
    }
    for i in (base + 1)..=t {
        if bay.stacks[s][i] > bay.stacks[s][i - 1] {
            return true;
        }
    }
    false
}

fn do_move(bay: &mut Bay, crane: &CraneParams, sol: &mut Solution, from: usize, to: usize) {
    let m = Move { src: from, dst: to };
    let dm = bay.apply_move(m);
    let time = crane.move_time(
        dm.src + 1,
        dm.src_tier + 1,
        dm.dst + 1,
        dm.dst_tier + 1,
        sol.last_dst_1based(),
    );
    sol.push(dm, time);
}

fn select_relocation_dst(
    bay: &Bay,
    state: &TghState,
    eligible: &[usize],
    c_group: Group,
    target_group: Group,
) -> usize {
    let mut clean = [false; MAX_STACKS];
    let mut empty = [false; MAX_STACKS];
    let mut top = [0; MAX_STACKS];
    let mut dirty_ldg = [0; MAX_STACKS];

    for &s in eligible {
        let is_empty_s = bay.is_empty(s);
        let is_clean_s = is_clean(bay, state, s);
        empty[s] = is_empty_s;
        clean[s] = is_clean_s;
        if !is_empty_s {
            top[s] = bay.top_group(s);
        }
        if !is_clean_s {
            dirty_ldg[s] = ldg(bay, state, s);
        }
    }

    let mut best: Option<(Group, usize)> = None;
    // Preference 1: clean stack that can already accommodate the container.
    for &s in eligible {
        if !bay.is_full(s) && clean[s] && (empty[s] || top[s] >= c_group) {
            let tg = if empty[s] { Group::MAX } else { top[s] };
            if best.is_none() || tg < best.unwrap().0 {
                best = Some((tg, s));
            }
        }
    }
    if let Some((_, s)) = best {
        return s;
    }

    let mut best2: Option<(Group, usize)> = None;
    // Preference 2: dirty stack whose largest dirty group does not exceed the container.
    for &s in eligible {
        if !clean[s] {
            let l = dirty_ldg[s];
            if l <= c_group && (best2.is_none() || l > best2.unwrap().0) {
                best2 = Some((l, s));
            }
        }
    }
    if let Some((_, s)) = best2 {
        return s;
    }

    let mut best3: Option<(Group, usize)> = None;
    // Preference 3: dirty stack with an intermediate largest dirty group.
    for &s in eligible {
        if !clean[s] {
            let l = dirty_ldg[s];
            if l > c_group && l < target_group {
                if best3.is_none() || l < best3.unwrap().0 {
                    best3 = Some((l, s));
                }
            }
        }
    }
    if let Some((_, s)) = best3 {
        return s;
    }

    let mut best4: Option<(Group, usize)> = None;
    // Preference 4: clean stack with the largest top group.
    for &s in eligible {
        if clean[s] && !empty[s] {
            let tg = top[s];
            if best4.is_none() || tg > best4.unwrap().0 {
                best4 = Some((tg, s));
            }
        }
    }
    if let Some((_, s)) = best4 {
        return s;
    }

    // Preference 5: dirty stack whose largest dirty group matches the current target group.
    for &s in eligible {
        if !clean[s] && dirty_ldg[s] == target_group {
            return s;
        }
    }

    // Final fallback: any empty stack, otherwise the first eligible one.
    for &s in eligible {
        if empty[s] {
            return s;
        }
    }

    eligible[0]
}

fn try_fulfillment(
    bay: &mut Bay,
    crane: &CraneParams,
    sol: &mut Solution,
    state: &TghState,
    from: usize,
    dst: usize,
    c_group: Group,
    prohibited: &[usize],
) -> bool {
    if bay.is_empty(dst) {
        return false;
    }
    let dst_top = bay.top_group(dst);
    let prohibited_mask = exclusion_mask(prohibited);

    let mut best: Option<(Group, usize)> = None;
    for s in 0..bay.s {
        if s == dst || s == from || bay.is_empty(s) || prohibited_mask[s] {
            continue;
        }
        let top_g = bay.top_group(s);
        if top_g > c_group && top_g <= dst_top {
            let t = bay.he(s) - 1;
            if is_container_dirty(bay, state, s, t) {
                if best.is_none() || top_g > best.unwrap().0 {
                    best = Some((top_g, s));
                }
            }
        }
    }

    if let Some((_, src_stack)) = best {
        // Fulfillment performs at most one auxiliary move, then relocation is reconsidered.
        do_move(bay, crane, sol, src_stack, dst);
        return true;
    }
    false
}

fn relocate_container(
    bay: &mut Bay,
    crane: &CraneParams,
    sol: &mut Solution,
    state: &TghState,
    from: usize,
    target_group: Group,
    prohibited: &[usize],
) {
    let prohibited_mask = exclusion_mask(prohibited);
    let mut eligible = [0usize; MAX_STACKS];

    loop {
        let c_group = bay.top_group(from);
        let mut n_eligible = 0usize;
        for s in 0..bay.s {
            if s != from && !bay.is_full(s) && !prohibited_mask[s] {
                eligible[n_eligible] = s;
                n_eligible += 1;
            }
        }

        let dst = select_relocation_dst(bay, state, &eligible[..n_eligible], c_group, target_group);

        if can_accommodate(bay, state, dst, c_group)
            && try_fulfillment(bay, crane, sol, state, from, dst, c_group, prohibited)
        {
            continue;
        }

        do_move(bay, crane, sol, from, dst);
        break;
    }
}

#[inline]
fn exclusion_mask(exclude: &[usize]) -> [bool; MAX_STACKS] {
    let mut mask = [false; MAX_STACKS];
    for &s in exclude {
        mask[s] = true;
    }
    mask
}

fn find_non_full(bay: &Bay, exclude: &[usize]) -> usize {
    let exclude_mask = exclusion_mask(exclude);
    for s in 0..bay.s {
        if !exclude_mask[s] && !bay.is_full(s) {
            return s;
        }
    }
    panic!("no non-full stack available");
}

fn find_stack_with_empty_slot(bay: &Bay, exclude: &[usize]) -> usize {
    find_non_full(bay, exclude)
}

fn select_temp_stack(bay: &Bay, state: &TghState, s_star: usize) -> usize {
    let mut best: Option<(usize, usize, usize)> = None;
    for s in 0..bay.s {
        if s == s_star || bay.is_full(s) {
            continue;
        }
        let height = bay.he(s);
        let dirty_bonus: usize = if is_dirty(bay, state, s) { 1 } else { 0 };
        let key = (height, dirty_bonus);
        if best.is_none() || key > (best.unwrap().0, best.unwrap().1) {
            best = Some((height, dirty_bonus, s));
        }
    }
    best.expect("no temp stack available").2
}

fn select_scenario3_case1_container(
    bay: &Bay,
    state: &TghState,
    ts: usize,
    s_star: usize,
) -> usize {
    let ts_is_clean = is_clean(bay, state, ts);

    if ts_is_clean {
        let ts_top = if bay.is_empty(ts) {
            Group::MAX
        } else {
            bay.top_group(ts)
        };

        let mut best: Option<(Group, usize)> = None;
        for s in 0..bay.s {
            if s == ts || s == s_star || bay.is_empty(s) {
                continue;
            }
            let g = bay.top_group(s);
            let t = bay.he(s) - 1;
            if is_container_dirty(bay, state, s, t) && g <= ts_top && !bay.is_full(ts) {
                if best.is_none() || g > best.unwrap().0 {
                    best = Some((g, s));
                }
            }
        }
        if let Some((_, s)) = best {
            return s;
        }

        let mut best2: Option<(Group, usize)> = None;
        for s in 0..bay.s {
            if s == ts || s == s_star || bay.is_empty(s) {
                continue;
            }
            let g = bay.top_group(s);
            let t = bay.he(s) - 1;
            if !is_container_dirty(bay, state, s, t) && g <= ts_top && !bay.is_full(ts) {
                if best2.is_none() || g > best2.unwrap().0 {
                    best2 = Some((g, s));
                }
            }
        }
        if let Some((_, s)) = best2 {
            return s;
        }

        let mut best3: Option<(Group, usize)> = None;
        for s in 0..bay.s {
            if s == ts || s == s_star || bay.is_empty(s) {
                continue;
            }
            let g = bay.top_group(s);
            let t = bay.he(s) - 1;
            if is_container_dirty(bay, state, s, t) {
                if best3.is_none() || g < best3.unwrap().0 {
                    best3 = Some((g, s));
                }
            }
        }
        if let Some((_, s)) = best3 {
            return s;
        }

        let mut best4: Option<(Group, usize)> = None;
        for s in 0..bay.s {
            if s == ts || s == s_star || bay.is_empty(s) {
                continue;
            }
            let g = bay.top_group(s);
            if best4.is_none() || g < best4.unwrap().0 {
                best4 = Some((g, s));
            }
        }
        best4.unwrap().1
    } else {
        let mut best: Option<(Group, usize)> = None;
        for s in 0..bay.s {
            if s == ts || s == s_star || bay.is_empty(s) {
                continue;
            }
            let g = bay.top_group(s);
            let t = bay.he(s) - 1;
            if is_container_dirty(bay, state, s, t) {
                if best.is_none() || g < best.unwrap().0 {
                    best = Some((g, s));
                }
            }
        }
        if let Some((_, s)) = best {
            return s;
        }

        let mut best2: Option<(Group, usize)> = None;
        for s in 0..bay.s {
            if s == ts || s == s_star || bay.is_empty(s) {
                continue;
            }
            let g = bay.top_group(s);
            if best2.is_none() || g < best2.unwrap().0 {
                best2 = Some((g, s));
            }
        }
        best2.unwrap().1
    }
}

fn select_scenario4_case2_container(
    bay: &Bay,
    state: &TghState,
    s_star: usize,
    sn_c: usize,
) -> usize {
    let mut best: Option<(Group, usize)> = None;
    for s in 0..bay.s {
        if s == s_star || s == sn_c || bay.is_empty(s) {
            continue;
        }
        let g = bay.top_group(s);
        let t = bay.he(s) - 1;
        if is_container_dirty(bay, state, s, t) {
            if best.is_none() || g < best.unwrap().0 {
                best = Some((g, s));
            }
        }
    }
    if let Some((_, s)) = best {
        return s;
    }

    let mut best2: Option<(Group, usize)> = None;
    for s in 0..bay.s {
        if s == s_star || s == sn_c || bay.is_empty(s) {
            continue;
        }
        let g = bay.top_group(s);
        if best2.is_none() || g < best2.unwrap().0 {
            best2 = Some((g, s));
        }
    }
    best2.unwrap().1
}

fn restore_displaced_fixed_container(
    bay: &mut Bay,
    crane: &CraneParams,
    sol: &mut Solution,
    state: &mut TghState,
    cs: usize,
    c_stack: usize,
    above: usize,
    current_group: Group,
) {
    for _ in 0..above {
        relocate_container(bay, crane, sol, state, c_stack, current_group, &[cs]);
    }

    do_move(bay, crane, sol, c_stack, cs);
    state.fix(cs);
}

fn giant_move_case1(
    bay: &mut Bay,
    crane: &CraneParams,
    sol: &mut Solution,
    state: &mut TghState,
    c_tier: usize,
    s_star: usize,
    current_group: Group,
) {
    let above_count = bay.he(s_star) - c_tier - 1;
    for _ in 0..above_count {
        relocate_container(bay, crane, sol, state, s_star, current_group, &[s_star]);
    }

    let ts = select_temp_stack(bay, state, s_star);

    let uf_s = state.uf(bay, s_star);

    let nslot = bay.empty_slots_except(&[ts, s_star]);

    // Case 1, scenario 1: enough external slots to clear the unfixed part directly.
    if nslot >= uf_s.saturating_sub(1) {
        do_move(bay, crane, sol, s_star, ts);
        let to_relocate = state.uf(bay, s_star);
        for _ in 0..to_relocate {
            relocate_container(bay, crane, sol, state, s_star, current_group, &[s_star, ts]);
        }
        do_move(bay, crane, sol, ts, s_star);
    // Case 1, scenario 2: one temporary slot is needed to shuttle c* away and back.
    } else if nslot > 0 {
        do_move(bay, crane, sol, s_star, ts);
        for _ in 0..(nslot - 1) {
            relocate_container(bay, crane, sol, state, s_star, current_group, &[s_star, ts]);
        }
        let temp_slot = find_stack_with_empty_slot(bay, &[ts, s_star]);
        do_move(bay, crane, sol, ts, temp_slot);
        let remaining = state.uf(bay, s_star);
        for _ in 0..remaining {
            do_move(bay, crane, sol, s_star, ts);
        }
        do_move(bay, crane, sol, temp_slot, s_star);
    // Case 1, scenario 3: no external slot is free, so another top container must be displaced.
    } else {
        let cs = select_scenario3_case1_container(bay, state, ts, s_star);
        let c_tier_before = bay.he(cs) - 1;
        let was_fixed = state.is_fixed(cs, c_tier_before);
        if was_fixed {
            state.unfix(cs);
        }

        do_move(bay, crane, sol, cs, ts);
        do_move(bay, crane, sol, s_star, cs);
        let remaining = state.uf(bay, s_star);
        for _ in 0..remaining {
            do_move(bay, crane, sol, s_star, ts);
        }
        do_move(bay, crane, sol, cs, s_star);

        if was_fixed {
            restore_displaced_fixed_container(
                bay,
                crane,
                sol,
                state,
                cs,
                ts,
                remaining,
                current_group,
            );
        }
    }
}

fn giant_move_case2(
    bay: &mut Bay,
    crane: &CraneParams,
    sol: &mut Solution,
    state: &mut TghState,
    c_star: (usize, usize),
    s_star: usize,
    current_group: Group,
) {
    let sn_c = c_star.0;
    let o_c = bay.he(sn_c) - c_star.1 - 1;
    let uf_s = state.uf(bay, s_star);

    if o_c == 0 && uf_s == 0 {
        do_move(bay, crane, sol, sn_c, s_star);
        return;
    }

    let f = uf_s + o_c;
    let nslot = bay.empty_slots_except(&[s_star, sn_c]);

    // Case 2 follows the four scenarios, depending on free slots in other stacks.
    if nslot >= f {
        giant_move_case2_s1(
            bay,
            crane,
            sol,
            state,
            sn_c,
            s_star,
            o_c,
            uf_s,
            current_group,
        );
    } else if nslot >= o_c + 1 {
        giant_move_case2_s2(
            bay,
            crane,
            sol,
            state,
            sn_c,
            s_star,
            o_c,
            nslot,
            current_group,
        );
    } else if nslot >= 1 {
        giant_move_case2_s3(
            bay,
            crane,
            sol,
            state,
            sn_c,
            s_star,
            o_c,
            nslot,
            current_group,
        );
    } else {
        giant_move_case2_s4(bay, crane, sol, state, sn_c, s_star, o_c, current_group);
    }
}

fn giant_move_case2_s1(
    bay: &mut Bay,
    crane: &CraneParams,
    sol: &mut Solution,
    state: &TghState,
    sn_c: usize,
    s_star: usize,
    mut from_sn: usize,
    mut from_s: usize,
    current_group: Group,
) {
    let prohibited = [s_star, sn_c];

    while from_sn > 0 || from_s > 0 {
        let from = pick_relocation_source(bay, sn_c, s_star, from_sn, from_s);
        relocate_container(bay, crane, sol, state, from, current_group, &prohibited);
        if from == sn_c {
            from_sn -= 1;
        } else {
            from_s -= 1;
        }
    }

    do_move(bay, crane, sol, sn_c, s_star);
}

fn giant_move_case2_s2(
    bay: &mut Bay,
    crane: &CraneParams,
    sol: &mut Solution,
    state: &TghState,
    sn_c: usize,
    s_star: usize,
    o_c: usize,
    nslot: usize,
    current_group: Group,
) {
    let prohibited = [s_star, sn_c];
    let mut from_sn = o_c;
    let mut from_s = nslot - 1 - o_c;

    for _ in 0..(nslot - 1) {
        let from = pick_relocation_source(bay, sn_c, s_star, from_sn, from_s);
        relocate_container(bay, crane, sol, state, from, current_group, &prohibited);
        if from == sn_c {
            from_sn -= 1;
        } else {
            from_s -= 1;
        }
    }

    let temp_slot = find_stack_with_empty_slot(bay, &[s_star, sn_c]);
    do_move(bay, crane, sol, sn_c, temp_slot);

    let remaining = state.uf(bay, s_star);
    for _ in 0..remaining {
        do_move(bay, crane, sol, s_star, sn_c);
    }

    do_move(bay, crane, sol, temp_slot, s_star);
}

fn giant_move_case2_s3(
    bay: &mut Bay,
    crane: &CraneParams,
    sol: &mut Solution,
    state: &TghState,
    sn_c: usize,
    s_star: usize,
    o_c: usize,
    nslot: usize,
    current_group: Group,
) {
    let prohibited = [s_star, sn_c];

    for _ in 0..(nslot - 1) {
        relocate_container(bay, crane, sol, state, sn_c, current_group, &prohibited);
    }

    let remaining_above = o_c - (nslot - 1);
    for _ in 0..remaining_above {
        do_move(bay, crane, sol, sn_c, s_star);
    }

    let temp_slot = find_stack_with_empty_slot(bay, &prohibited);
    do_move(bay, crane, sol, sn_c, temp_slot);

    let total_uf = state.uf(bay, s_star);
    for _ in 0..total_uf {
        do_move(bay, crane, sol, s_star, sn_c);
    }

    do_move(bay, crane, sol, temp_slot, s_star);
}

fn giant_move_case2_s4(
    bay: &mut Bay,
    crane: &CraneParams,
    sol: &mut Solution,
    state: &mut TghState,
    sn_c: usize,
    s_star: usize,
    o_c: usize,
    current_group: Group,
) {
    let cs = select_scenario4_case2_container(bay, state, s_star, sn_c);
    let c_tier_before = bay.he(cs) - 1;
    let was_fixed = state.is_fixed(cs, c_tier_before);
    if was_fixed {
        state.unfix(cs);
    }

    do_move(bay, crane, sol, cs, s_star);

    for _ in 0..o_c {
        do_move(bay, crane, sol, sn_c, s_star);
    }

    do_move(bay, crane, sol, sn_c, cs);

    let total_uf = state.uf(bay, s_star);
    for _ in 0..total_uf {
        do_move(bay, crane, sol, s_star, sn_c);
    }

    // Move c* to s*
    do_move(bay, crane, sol, cs, s_star);

    if was_fixed {
        debug_assert!(total_uf >= o_c + 1);
        let above = total_uf - (o_c + 1);
        restore_displaced_fixed_container(bay, crane, sol, state, cs, sn_c, above, current_group);
    }
}

fn pick_relocation_source(
    bay: &Bay,
    sn_c: usize,
    s_star: usize,
    from_sn: usize,
    from_s: usize,
) -> usize {
    if from_sn == 0 {
        s_star
    } else if from_s == 0 {
        sn_c
    } else {
        let g_sn = bay.top_group(sn_c);
        let g_s = bay.top_group(s_star);
        if g_sn > g_s {
            sn_c
        } else {
            s_star
        }
    }
}

fn giant_move(
    bay: &mut Bay,
    crane: &CraneParams,
    sol: &mut Solution,
    state: &mut TghState,
    c_star: (usize, usize),
    s_star: usize,
    current_group: Group,
) {
    if c_star.0 == s_star {
        giant_move_case1(bay, crane, sol, state, c_star.1, s_star, current_group);
    } else {
        giant_move_case2(bay, crane, sol, state, c_star, s_star, current_group);
    }
}

fn mark_clean_fixed(bay: &Bay, state: &mut TghState, g: Group) {
    for s in 0..bay.s {
        while state.fixed_height[s] < bay.he(s) && bay.stacks[s][state.fixed_height[s]] == g {
            state.fix(s);
        }
    }
}

fn fill_unfixed_with_group(
    bay: &Bay,
    state: &TghState,
    g: Group,
    result: &mut Vec<(usize, usize)>,
) {
    result.clear();
    for s in 0..bay.s {
        for t in state.fixed_height[s]..bay.he(s) {
            if bay.stacks[s][t] == g {
                result.push((s, t));
            }
        }
    }
}

fn fill_candidate_stacks(
    bay: &Bay,
    state: &TghState,
    candidates: &mut [usize; MAX_STACKS],
    no_reloc: &mut [usize; MAX_STACKS],
) -> Option<(usize, usize)> {
    let total_empty_slots: usize = (0..bay.s).map(|i| bay.h - bay.he(i)).sum();

    let mut n_candidates = 0usize;
    let mut n_no_reloc = 0usize;
    for s in 0..bay.s {
        if state.fixed_height[s] >= bay.h {
            continue;
        }
        let uf_s = state.uf(bay, s);
        let empty_here = bay.h - bay.he(s);
        let empty_others = total_empty_slots - empty_here;
        if uf_s <= empty_others {
            candidates[n_candidates] = s;
            n_candidates += 1;
            if uf_s == 0 {
                no_reloc[n_no_reloc] = s;
                n_no_reloc += 1;
            }
        }
    }
    if n_candidates == 0 {
        return None;
    }

    Some((n_candidates, n_no_reloc))
}

fn select_target_pair(
    bay: &Bay,
    state: &TghState,
    con_list: &[(usize, usize)],
    stk_list: &[usize],
) -> ((usize, usize), usize) {
    let mut best: Option<((usize, usize), usize, (usize, usize))> = None;
    let mut stack_uf = [0usize; MAX_STACKS];
    let mut stack_fixed = [0usize; MAX_STACKS];

    for &s in stk_list {
        stack_uf[s] = state.uf(bay, s);
        stack_fixed[s] = state.fixed_height[s];
    }

    for &(cs, ct) in con_list {
        let o_c = bay.he(cs) - ct - 1;
        for &s in stk_list {
            let uf_s = stack_uf[s];
            let f = if cs == s { uf_s } else { uf_s + o_c };
            let fixed_count = stack_fixed[s];
            let eval = (f, fixed_count);

            if best.is_none() || eval < best.unwrap().2 {
                best = Some(((cs, ct), s, eval));
            }
        }
    }

    let (c_star, s_star, _) = best.expect("no target pair found");
    (c_star, s_star)
}

fn select_target_pair_randomized<R: Rng>(
    _bay: &Bay,
    _state: &TghState,
    con_list: &[(usize, usize)],
    stk_list: &[usize],
    rng: &mut R,
) -> ((usize, usize), usize) {
    debug_assert!(!con_list.is_empty());
    debug_assert!(!stk_list.is_empty());

    let n_con = con_list.len();
    let n_stk = stk_list.len();
    let idx = rng.gen_range(0..(n_con * n_stk));

    let c_idx = idx / n_stk;
    let s_idx = idx % n_stk;
    (con_list[c_idx], stk_list[s_idx])
}

fn run_tgh_with_selector<F>(
    bay: &mut Bay,
    crane: &CraneParams,
    mut select_pair: F,
) -> Option<Solution>
where
    F: FnMut(&Bay, &TghState, &[(usize, usize)], &[usize]) -> ((usize, usize), usize),
{
    let mut sol = Solution::new();
    let mut state = TghState::new();
    let mut con_list: Vec<(usize, usize)> = Vec::new();
    let mut candidate_stacks = [0usize; MAX_STACKS];
    let mut no_reloc_stacks = [0usize; MAX_STACKS];

    // Containers are fixed group by group, from largest label to smallest.
    for g in (1..=bay.p as Group).rev() {
        mark_clean_fixed(bay, &mut state, g);

        fill_unfixed_with_group(bay, &state, g, &mut con_list);

        while !con_list.is_empty() {
            let (n_candidates, n_no_reloc) =
                fill_candidate_stacks(bay, &state, &mut candidate_stacks, &mut no_reloc_stacks)?;
            let stk_list: &[usize] = if n_no_reloc > 0 {
                &no_reloc_stacks[..n_no_reloc]
            } else {
                &candidate_stacks[..n_candidates]
            };

            let (c_star, s_star) = select_pair(bay, &state, &con_list, stk_list);

            // A giant move fixes one target container while preserving the fixed prefixes.
            giant_move(bay, crane, &mut sol, &mut state, c_star, s_star, g);

            state.fix(s_star);
            mark_clean_fixed(bay, &mut state, g);

            fill_unfixed_with_group(bay, &state, g, &mut con_list);
        }
    }

    sol.acc = bay.acc();
    Some(sol)
}

/// Runs the deterministic target-guided heuristic in place on `bay`.
///
/// The bay is mutated into the final TGH layout. The heuristic returns `None`
/// when no candidate target stack exists, which in practice signals an
/// unsolvable layout under the TGH assumptions.
pub fn tgh(bay: &mut Bay, crane: &CraneParams) -> Option<Solution> {
    run_tgh_with_selector(bay, crane, select_target_pair)
}

/// Runs a randomized variant of TGH by sampling target container / stack pairs.
///
/// The relocation logic is unchanged; only the target-pair selector is
/// randomized. The bay is mutated in place in the same way as [`tgh`].
pub fn tgh_randomized<R: Rng>(bay: &mut Bay, crane: &CraneParams, rng: &mut R) -> Option<Solution> {
    run_tgh_with_selector(bay, crane, |bay, state, con_list, stk_list| {
        select_target_pair_randomized(bay, state, con_list, stk_list, rng)
    })
}
