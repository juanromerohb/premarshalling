//! Crane-time model used by the heuristics and exact methods.
//!
//! This module implements the RTG travel-time computation. The public
//! methods expose the derived quantities used throughout the crate: horizontal
//! and vertical travel times, a conservative per-move lower bound, and the full
//! time of a relocation when the previous crane position is known.

use crate::bay::{MAX_HEIGHT, MAX_STACKS};

#[derive(Clone, Debug)]
/// RTG crane parameters together with small lookup tables for common distances.
///
/// Cached entries cover the crate's fixed maximum bay dimensions so repeated
/// move-time evaluations stay allocation-free.
pub struct CraneParams {
    pub ch: f64,
    pub cw: f64,
    pub vsep: f64,
    pub hsep: f64,
    pub msep: f64,
    pub v_max_vl: f64,
    pub v_max_vu: f64,
    pub v_max_rl: f64,
    pub v_max_ru: f64,
    pub d_v_max: f64,
    pub d_r_max: f64,
    pub max_h: usize,
    v_l_cache: [f64; MAX_HEIGHT + 1],
    v_u_cache: [f64; MAX_HEIGHT + 1],
    r_u_cache: [[f64; MAX_STACKS + 1]; MAX_STACKS + 1],
    r_l_cache: [[f64; MAX_STACKS + 1]; MAX_STACKS + 1],
    t_min_cache: [f64; MAX_HEIGHT + 1],
    t_min_prefix_cache: [f64; MAX_HEIGHT + 1],
}

impl CraneParams {
    /// Builds the crane model for bays of height at most `max_h`.
    pub fn new(max_h: usize) -> Self {
        debug_assert!(max_h <= MAX_HEIGHT);

        let mut crane = Self {
            ch: 2.591,
            cw: 2.438,
            vsep: 2.0,
            hsep: 1.0,
            msep: 0.3,
            v_max_vl: 0.50,
            v_max_vu: 1.00,
            v_max_rl: 7.0 / 6.0,
            v_max_ru: 13.0 / 6.0,
            d_v_max: 2.65,
            d_r_max: 2.50,
            max_h,
            v_l_cache: [0.0; MAX_HEIGHT + 1],
            v_u_cache: [0.0; MAX_HEIGHT + 1],
            r_u_cache: [[0.0; MAX_STACKS + 1]; MAX_STACKS + 1],
            r_l_cache: [[0.0; MAX_STACKS + 1]; MAX_STACKS + 1],
            t_min_cache: [0.0; MAX_HEIGHT + 1],
            t_min_prefix_cache: [0.0; MAX_HEIGHT + 1],
        };

        for h in 1..=max_h {
            crane.v_l_cache[h] = crane.compute_v_l(h);
            crane.v_u_cache[h] = crane.compute_v_u(h);
        }
        for s in 0..=MAX_STACKS {
            for k in 1..=MAX_STACKS {
                crane.r_u_cache[s][k] = crane.compute_r_u(s, k);
            }
        }
        for s in 1..=MAX_STACKS {
            for k in 1..=MAX_STACKS {
                crane.r_l_cache[s][k] = crane.compute_r_l(s, k);
            }
        }

        let ru_12 = crane.r_u_cache[1][2];
        let rl_12 = crane.r_l_cache[1][2];
        let vu_h = crane.v_u_cache[max_h];
        let mut prefix = 0.0;
        for h in 1..=max_h {
            let t = ru_12 + crane.v_l_cache[h] + rl_12 + vu_h;
            crane.t_min_cache[h] = t;
            prefix += t;
            crane.t_min_prefix_cache[h] = prefix;
        }

        crane
    }

    /// Returns the travel time for distance `x` under the paper's piecewise motion model.
    ///
    /// The crane accelerates up to distance `d_max`, cruises if the trip is
    /// long enough, and then decelerates symmetrically. The function is shared
    /// by horizontal and vertical loaded / unloaded movements.
    pub fn travel_time(&self, x: f64, v_max: f64, d_max: f64) -> f64 {
        if x < 0.0 {
            return 0.0;
        }
        if x < 2.0 * d_max {
            2.0 * (2.0 * d_max * x).sqrt() / v_max
        } else {
            (x + 2.0 * d_max) / v_max
        }
    }

    /// Returns the vertical distance from the top level to tier `h`.
    pub fn d_v(&self, h: usize) -> f64 {
        self.vsep + (self.max_h - h + 1) as f64 * self.ch
    }

    /// Returns the horizontal distance between stacks `s` and `k`.
    ///
    /// Stack `0` denotes the crane's initial safety position, so
    /// `d_r(0, k)` uses the dedicated truck-to-first-stack geometry.
    pub fn d_r(&self, s: usize, k: usize) -> f64 {
        if s == 0 {
            (k - 1) as f64 * self.msep + k as f64 * self.cw + self.hsep
        } else {
            ((s as isize - k as isize).unsigned_abs() as f64) * (self.msep + self.cw)
        }
    }

    #[inline]
    /// Returns the unloaded horizontal travel time `r^u(s, k)`.
    pub fn r_u(&self, s: usize, k: usize) -> f64 {
        if s <= MAX_STACKS && k <= MAX_STACKS {
            self.r_u_cache[s][k]
        } else {
            self.compute_r_u(s, k)
        }
    }

    #[inline]
    /// Returns the loaded horizontal travel time `r^l(s, k)`.
    pub fn r_l(&self, s: usize, k: usize) -> f64 {
        if s <= MAX_STACKS && k <= MAX_STACKS {
            self.r_l_cache[s][k]
        } else {
            self.compute_r_l(s, k)
        }
    }

    /// Returns the twistlock time for lifting from tier `h`.
    pub fn twistlock(&self, h: usize) -> f64 {
        5.0 * (self.max_h - h + 1) as f64
    }

    #[inline]
    /// Returns the loading time `v^l(h)`.
    ///
    /// This includes descending unloaded, twistlock handling, and ascending loaded.
    pub fn v_l(&self, h: usize) -> f64 {
        if h >= 1 && h <= self.max_h {
            self.v_l_cache[h]
        } else {
            self.compute_v_l(h)
        }
    }

    #[inline]
    /// Returns the unloading time `v^u(h)`.
    pub fn v_u(&self, h: usize) -> f64 {
        if h >= 1 && h <= self.max_h {
            self.v_u_cache[h]
        } else {
            self.compute_v_u(h)
        }
    }

    #[inline]
    fn compute_r_u(&self, s: usize, k: usize) -> f64 {
        self.travel_time(self.d_r(s, k), self.v_max_ru, self.d_r_max)
    }

    #[inline]
    fn compute_r_l(&self, s: usize, k: usize) -> f64 {
        self.travel_time(self.d_r(s, k), self.v_max_rl, self.d_r_max)
    }

    #[inline]
    fn compute_v_l(&self, h: usize) -> f64 {
        let dv = self.d_v(h);
        self.travel_time(dv, self.v_max_vu, self.d_v_max)
            + self.twistlock(h)
            + self.travel_time(dv, self.v_max_vl, self.d_v_max)
    }

    #[inline]
    fn compute_v_u(&self, h: usize) -> f64 {
        let dv = self.d_v(h);
        self.travel_time(dv, self.v_max_vl, self.d_v_max)
            + self.travel_time(dv, self.v_max_vu, self.d_v_max)
    }

    #[inline]
    /// Returns the conservative lower bound `t_min(h)` used by the LCT models.
    ///
    /// It assumes the source and destination stacks are adjacent and the
    /// destination tier is the fastest possible one.
    pub fn t_min(&self, h: usize) -> f64 {
        if h >= 1 && h <= self.max_h {
            self.t_min_cache[h]
        } else {
            self.r_u(1, 2) + self.v_l(h) + self.r_l(1, 2) + self.v_u(self.max_h)
        }
    }

    #[inline]
    /// Returns the sum of `t_min(h)` over tiers `from..=to`.
    pub fn t_min_sum(&self, from: usize, to: usize) -> f64 {
        if from > to {
            return 0.0;
        }
        if from >= 1 && to <= self.max_h {
            let left = if from > 1 {
                self.t_min_prefix_cache[from - 1]
            } else {
                0.0
            };
            self.t_min_prefix_cache[to] - left
        } else {
            let mut total = 0.0;
            for h in from..=to {
                total += self.t_min(h);
            }
            total
        }
    }

    #[inline]
    /// Returns the crane time of one relocation.
    ///
    /// All stack and tier arguments use the paper's one-based notation. The
    /// `prev_dst` argument is `0` for the first move and otherwise denotes the
    /// previous destination stack in one-based indexing.
    pub fn move_time(
        &self,
        src: usize,
        src_tier: usize,
        dst: usize,
        dst_tier: usize,
        prev_dst: usize,
    ) -> f64 {
        let r_u = if prev_dst <= MAX_STACKS && src <= MAX_STACKS {
            self.r_u_cache[prev_dst][src]
        } else {
            self.compute_r_u(prev_dst, src)
        };
        let r_l = if src <= MAX_STACKS && dst <= MAX_STACKS {
            self.r_l_cache[src][dst]
        } else {
            self.compute_r_l(src, dst)
        };
        let v_l = if src_tier >= 1 && src_tier <= self.max_h {
            self.v_l_cache[src_tier]
        } else {
            self.compute_v_l(src_tier)
        };
        let v_u = if dst_tier >= 1 && dst_tier <= self.max_h {
            self.v_u_cache[dst_tier]
        } else {
            self.compute_v_u(dst_tier)
        };
        r_u + v_l + r_l + v_u
    }
}
