/// Maximum number of stacks supported by the fixed-size bay representation.
pub const MAX_STACKS: usize = 20;
/// Maximum stack height supported by the fixed-size bay representation.
pub const MAX_HEIGHT: usize = 15;
/// Maximum number of containers in a bay.
pub const MAX_CONTAINERS: usize = MAX_STACKS * MAX_HEIGHT;
/// Maximum number of priority groups supported by the fixed-size counters.
pub const MAX_GROUPS: usize = MAX_CONTAINERS + 1;

/// Container priority / group identifier.
///
/// Lower values correspond to earlier retrieval.
pub type Group = u8;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
/// A relocation identified only by source and destination stacks.
///
/// Stack indices are zero-based in the in-memory bay representation.
pub struct Move {
    pub src: usize,
    pub dst: usize,
}

#[derive(Clone, Copy, Debug)]
/// A relocation enriched with the concrete source and destination tiers.
///
/// This is the move representation stored in [`crate::solution::Solution`].
pub struct DetailedMove {
    pub src: usize,
    pub src_tier: usize,
    pub dst: usize,
    pub dst_tier: usize,
    pub group: Group,
}

#[derive(Clone, Copy, Debug)]
/// Compact bay state used by the exact and heuristic algorithms.
pub struct Bay {
    pub s: usize,
    pub h: usize,
    pub stacks: [[Group; MAX_HEIGHT]; MAX_STACKS],
    pub heights: [u8; MAX_STACKS],
    pub p: usize,
    pub c: usize,
    pub group_count: [u16; MAX_GROUPS],
    pub cum_count: [u16; MAX_GROUPS],
}

impl Bay {
    /// Builds a bay from stack vectors listed from bottom tier to top tier.
    pub fn from_vecs(s: usize, h: usize, input: &[Vec<Group>]) -> Self {
        let mut bay = Self::empty(s, h);

        for (si, st) in input.iter().enumerate() {
            for (ti, &g) in st.iter().enumerate() {
                bay.stacks[si][ti] = g;
            }
            bay.heights[si] = st.len() as u8;
        }

        bay.recompute_counts();
        bay
    }

    /// Creates an empty bay with the requested number of stacks and height.
    pub fn empty(s: usize, h: usize) -> Self {
        Self {
            s,
            h,
            stacks: [[0; MAX_HEIGHT]; MAX_STACKS],
            heights: [0; MAX_STACKS],
            p: 0,
            c: 0,
            group_count: [0; MAX_GROUPS],
            cum_count: [0; MAX_GROUPS],
        }
    }

    fn recompute_counts(&mut self) {
        self.group_count = [0; MAX_GROUPS];
        self.c = 0;
        self.p = 0;

        for si in 0..self.s {
            for ti in 0..self.heights[si] as usize {
                let g = self.stacks[si][ti];
                if g > 0 {
                    self.c += 1;
                    self.group_count[g as usize] += 1;
                    if (g as usize) > self.p {
                        self.p = g as usize;
                    }
                }
            }
        }

        self.cum_count = [0; MAX_GROUPS];
        for i in 1..=self.p {
            self.cum_count[i] = self.cum_count[i - 1] + self.group_count[i];
        }
    }

    #[inline]
    /// Returns the height of stack `s`.
    pub fn he(&self, s: usize) -> usize {
        self.heights[s] as usize
    }

    #[inline]
    /// Returns the group stored at stack `s`, tier `t`, or `0` if the tier is empty.
    pub fn g(&self, s: usize, t: usize) -> Group {
        if t < self.he(s) {
            self.stacks[s][t]
        } else {
            0
        }
    }

    #[inline]
    /// Returns whether stack `s` contains no containers.
    pub fn is_empty(&self, s: usize) -> bool {
        self.heights[s] == 0
    }

    #[inline]
    /// Returns whether stack `s` is at capacity.
    pub fn is_full(&self, s: usize) -> bool {
        self.heights[s] as usize >= self.h
    }

    #[inline]
    /// Returns the topmost group of stack `s`.
    pub fn top_group(&self, s: usize) -> Group {
        debug_assert!(self.heights[s] > 0, "stack is empty");
        self.stacks[s][self.heights[s] as usize - 1]
    }

    #[inline]
    /// Applies a relocation and returns the fully specified move that was executed.
    pub fn apply_move(&mut self, m: Move) -> DetailedMove {
        debug_assert!(self.heights[m.src] > 0, "source stack is empty");
        let src_tier = self.heights[m.src] as usize - 1;
        let group = self.stacks[m.src][src_tier];
        self.heights[m.src] -= 1;

        let dst_tier = self.heights[m.dst] as usize;
        debug_assert!(dst_tier < self.h, "destination stack is full");
        self.stacks[m.dst][dst_tier] = group;
        self.heights[m.dst] += 1;

        DetailedMove {
            src: m.src,
            src_tier,
            dst: m.dst,
            dst_tier,
            group,
        }
    }

    #[inline]
    /// Reverts a previously applied detailed move.
    pub fn undo_move(&mut self, dm: DetailedMove) {
        self.heights[dm.dst] -= 1;
        self.heights[dm.src] += 1;
        self.stacks[dm.src][dm.src_tier] = dm.group;
    }

    /// Returns whether the container at `(s, t)` is blocked by a larger group above it.
    ///
    /// This is the local blocking predicate used throughout the CPMP-LCT code.
    pub fn is_blocked(&self, s: usize, t: usize) -> bool {
        let g = self.stacks[s][t];
        if g == 0 {
            return false;
        }
        let he = self.he(s);
        for h in (t + 1)..he {
            if self.stacks[s][h] > g {
                return true;
            }
        }
        false
    }

    #[inline]
    /// Returns whether the container is locally reachable inside its own stack.
    ///
    /// The method does not enforce the global accessibility definition,
    /// which also depends on earlier-priority containers.
    pub fn is_relatively_accessible(&self, s: usize, t: usize) -> bool {
        !self.is_blocked(s, t)
    }

    /// Returns whether the container at `(s, t)` is "misoverlaid".
    ///
    /// A container is misoverlaid if it sits above a smaller-priority container
    /// or above any container that is already blocked.
    pub fn is_misoverlaid(&self, s: usize, t: usize) -> bool {
        let g = self.stacks[s][t];
        if g == 0 {
            return false;
        }
        for h in 0..t {
            if self.stacks[s][h] > 0 && self.stacks[s][h] < g {
                return true;
            }
        }
        for h in 0..t {
            if self.is_blocked(s, h) {
                return true;
            }
        }
        false
    }

    /// Returns the tier of the first blocker above `(s, t)`, or the stack height if none exists.
    pub fn bh(&self, s: usize, t: usize) -> usize {
        let g = self.stacks[s][t];
        let he = self.he(s);
        for h in (t + 1)..he {
            if self.stacks[s][h] > g {
                return h;
            }
        }
        he
    }

    /// Returns the earliest blocking tier among blocked containers with group at most `pn`.
    pub fn bh_pn(&self, s: usize, pn: Group) -> usize {
        let he = self.he(s);
        let mut min_bh = he;
        for t in 0..he {
            let g = self.stacks[s][t];
            if g > 0 && g <= pn && self.is_blocked(s, t) {
                let b = self.bh(s, t);
                if b < min_bh {
                    min_bh = b;
                }
            }
        }
        min_bh
    }

    /// Returns the smallest group that is still blocked in the current layout.
    ///
    /// If every container is locally unblocked, this returns `p + 1`.
    pub fn p_star(&self) -> usize {
        for p in 1..=self.p {
            for si in 0..self.s {
                let he = self.he(si);
                for t in 0..he {
                    if self.stacks[si][t] == p as Group && self.is_blocked(si, t) {
                        return p;
                    }
                }
            }
        }
        self.p + 1
    }

    /// Returns the number of accessible containers in the CPMP-LCT sense.
    ///
    /// All groups strictly smaller than [`Self::p_star`] are fully accessible;
    /// in the critical group itself, only locally unblocked containers count.
    pub fn acc(&self) -> usize {
        let ps = self.p_star();
        if ps > self.p {
            return self.c;
        }
        let full_groups = if ps > 0 {
            self.cum_count[ps - 1] as usize
        } else {
            0
        };
        let mut n_ps = 0usize;
        for si in 0..self.s {
            for t in 0..self.he(si) {
                if self.stacks[si][t] == ps as Group && !self.is_blocked(si, t) {
                    n_ps += 1;
                }
            }
        }
        full_groups + n_ps
    }

    /// Returns the stacks that contain at least one blocked container of group `p`.
    pub fn stacks_with_blocked(&self, p: Group) -> Vec<usize> {
        let mut result = Vec::new();
        for si in 0..self.s {
            for t in 0..self.he(si) {
                if self.stacks[si][t] == p && self.is_blocked(si, t) {
                    result.push(si);
                    break;
                }
            }
        }
        result
    }

    /// Returns the stacks that contain at least one misoverlaid container of group `p`.
    pub fn stacks_with_misoverlaid(&self, p: Group) -> Vec<usize> {
        let mut result = Vec::new();
        for si in 0..self.s {
            for t in 0..self.he(si) {
                if self.stacks[si][t] == p && self.is_misoverlaid(si, t) {
                    result.push(si);
                    break;
                }
            }
        }
        result
    }

    /// Returns the largest group that still has a misoverlaid container.
    pub fn max_misoverlaid_group(&self) -> Option<Group> {
        for p in (1..=self.p).rev() {
            if !self.stacks_with_misoverlaid(p as Group).is_empty() {
                return Some(p as Group);
            }
        }
        None
    }

    /// Returns the smallest group that still has a blocked container.
    pub fn min_blocked_group(&self) -> Option<Group> {
        for p in 1..=self.p {
            if !self.stacks_with_blocked(p as Group).is_empty() {
                return Some(p as Group);
            }
        }
        None
    }

    /// Returns the `p`-reduced bay used by HLCT.
    ///
    /// All groups strictly larger than `p` are merged into `p + 1`.
    pub fn p_reduction(&self, p: Group) -> Bay {
        let mut reduced = *self;
        for si in 0..reduced.s {
            for ti in 0..reduced.he(si) {
                if reduced.stacks[si][ti] > p {
                    reduced.stacks[si][ti] = p + 1;
                }
            }
        }
        reduced.recompute_counts();
        reduced
    }

    /// Counts empty slots outside the stacks listed in `exclude`.
    pub fn empty_slots_except(&self, exclude: &[usize]) -> usize {
        let mut count = 0;
        for si in 0..self.s {
            if !exclude.contains(&si) {
                count += self.h - self.he(si);
            }
        }
        count
    }

    /// Counts locally accessible containers of the given group.
    pub fn count_relatively_accessible(&self, p: Group) -> usize {
        let mut count = 0;
        for si in 0..self.s {
            for t in 0..self.he(si) {
                if self.stacks[si][t] == p && !self.is_blocked(si, t) {
                    count += 1;
                }
            }
        }
        count
    }

    /// Returns whether every stack is clean, i.e. non-increasing from bottom to top.
    pub fn is_sorted(&self) -> bool {
        for si in 0..self.s {
            for t in 1..self.he(si) {
                if self.stacks[si][t] > self.stacks[si][t - 1] {
                    return false;
                }
            }
        }
        true
    }

    /// Pretty-prints the bay as rows from top tier to bottom tier.
    pub fn display(&self) -> String {
        let mut lines = Vec::new();
        for t in (0..self.h).rev() {
            let mut row = Vec::new();
            for si in 0..self.s {
                if t < self.he(si) {
                    row.push(format!("{:3}", self.stacks[si][t]));
                } else {
                    row.push("  .".to_string());
                }
            }
            lines.push(row.join(" "));
        }
        lines.join("\n")
    }
}
