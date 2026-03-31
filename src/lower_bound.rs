//! Crane-time lower bound for the CPMP-LCT.
//!
//! The implementation follows the lower-bound construction described in the
//! main paper and operates directly on the current bay layout.

use crate::bay::{Bay, Group, MAX_HEIGHT, MAX_STACKS};
use crate::crane::CraneParams;

/// Returns the lower bound on crane time required to reach `n` accessible containers.
pub fn lower_bound(bay: &Bay, crane: &CraneParams, n: usize) -> f64 {
    if n == 0 {
        return 0.0;
    }
    if n > bay.c {
        return f64::INFINITY;
    }
    if bay.acc() >= n {
        return 0.0;
    }

    let pn = compute_pn(bay, n);
    let n1 = if pn == 0 {
        0
    } else {
        bay.cum_count[pn] as usize
    };
    let n2 = n - n1;
    let pn_group = pn as Group;
    let bh_pn_by_stack = compute_bh_pn_by_stack(bay, pn_group);

    let lb1 = lb1_value_with_bh(bay, crane, &bh_pn_by_stack);
    let lb2 = if n2 > 0 {
        lb2_value_with_bh(bay, crane, pn_group, n2, &bh_pn_by_stack)
    } else {
        0.0
    };

    lb1 + lb2
}

fn compute_pn(bay: &Bay, n: usize) -> usize {
    if bay.p == 0 {
        return 0;
    }

    let mut lo = 1usize;
    let mut hi = bay.p;
    let mut pn = 0usize;
    while lo <= hi {
        let mid = lo + (hi - lo) / 2;
        if bay.cum_count[mid] as usize <= n {
            pn = mid;
            lo = mid + 1;
        } else {
            hi = mid - 1;
        }
    }
    pn
}

fn compute_bh_pn_by_stack(bay: &Bay, pn: Group) -> [usize; MAX_STACKS] {
    let mut bh_pn_by_stack = [0usize; MAX_STACKS];
    for s in 0..bay.s {
        bh_pn_by_stack[s] = bay.bh_pn(s, pn);
    }
    bh_pn_by_stack
}

fn lb1_value_with_bh(bay: &Bay, crane: &CraneParams, bh_pn_by_stack: &[usize]) -> f64 {
    let mut total = 0.0;
    for s in 0..bay.s {
        let bh = bh_pn_by_stack[s];
        let he = bay.he(s);
        if bh < he {
            total += crane.t_min_sum(bh + 1, he);
        }
    }
    total
}

fn lb2_value_with_bh(
    bay: &Bay,
    crane: &CraneParams,
    pn: Group,
    n2: usize,
    bh_pn_by_stack: &[usize],
) -> f64 {
    let pn1 = pn + 1;
    if pn1 as usize > bay.p {
        return 0.0;
    }

    let pn1_count = bay.group_count[pn1 as usize] as usize;
    if pn1_count == 0 {
        return 0.0;
    }

    let mut re_hist = [0usize; MAX_HEIGHT + 1];
    let mut total_re_values = 0usize;
    let mut h_max = 0usize;

    for s in 0..bay.s {
        let bh_pn_s = bh_pn_by_stack[s];
        let he_s = bay.he(s);
        if he_s == 0 {
            continue;
        }

        let mut pre_re = 0usize;
        let mut stack_has_positive_re = false;
        let mut has_gt_above_to_bh = [false; MAX_HEIGHT];
        if bh_pn_s > 0 {
            let mut seen_gt = false;
            for t in (0..bh_pn_s).rev() {
                has_gt_above_to_bh[t] = seen_gt;
                if bay.stacks[s][t] > pn1 {
                    seen_gt = true;
                }
            }
        }
        for t in (0..he_s).rev() {
            if bay.stacks[s][t] != pn1 {
                continue;
            }

            let re = if t < bh_pn_s {
                if has_gt_above_to_bh[t] {
                    (bh_pn_s - (t + 1)).saturating_sub(pre_re)
                } else {
                    0
                }
            } else {
                0
            };

            if re > 0 {
                stack_has_positive_re = true;
            }

            debug_assert!(re <= MAX_HEIGHT);
            re_hist[re] += 1;
            total_re_values += 1;
            pre_re += re;
        }

        if stack_has_positive_re {
            let stack_h_max = bh_pn_s;
            if stack_h_max > h_max {
                h_max = stack_h_max;
            }
        }
    }

    if total_re_values == 0 || h_max == 0 {
        return 0.0;
    }

    let take = n2.min(total_re_values);
    if take == 0 {
        return 0.0;
    }

    let mut remaining = take;
    let mut sum_re = 0usize;
    for (re, &count) in re_hist.iter().enumerate() {
        if remaining == 0 {
            break;
        }
        if count == 0 {
            continue;
        }
        let picked = remaining.min(count);
        sum_re += re * picked;
        remaining -= picked;
    }

    if sum_re == 0 {
        return 0.0;
    }

    crane.t_min(h_max) * sum_re as f64
}

/// Returns the sorted `re` values used by the second lower-bound component.
pub fn compute_re_values(bay: &Bay, pn: Group) -> Vec<usize> {
    let mut re_values: Vec<usize> = compute_re_data(bay, pn)
        .into_iter()
        .map(|(re, _)| re)
        .collect();
    re_values.sort_unstable();
    re_values
}

fn compute_re_data(bay: &Bay, pn: Group) -> Vec<(usize, usize)> {
    let pn1 = pn + 1;
    if pn1 as usize > bay.p {
        return Vec::new();
    }

    let bh_pn_by_stack = compute_bh_pn_by_stack(bay, pn);
    let mut re_data = Vec::with_capacity(bay.group_count[pn1 as usize] as usize);

    for s in 0..bay.s {
        let bh_pn_s = bh_pn_by_stack[s];
        let he_s = bay.he(s);
        if he_s == 0 {
            continue;
        }

        let mut pre_re = 0usize;
        for t in (0..he_s).rev() {
            if bay.stacks[s][t] == pn1 {
                let re = if t < bh_pn_s {
                    let has_gt_pn1 = ((t + 1)..bh_pn_s).any(|u| bay.stacks[s][u] > pn1);
                    if has_gt_pn1 {
                        (bh_pn_s - (t + 1)).saturating_sub(pre_re)
                    } else {
                        0
                    }
                } else {
                    0
                };
                re_data.push((re, he_s));
                pre_re += re;
            }
        }
    }

    re_data.sort_unstable_by_key(|(re, _)| *re);
    re_data
}
