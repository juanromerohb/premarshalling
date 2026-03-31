use crate::bay::DetailedMove;

#[derive(Clone, Debug)]
/// Sequence of executed relocations together with incremental crane-time data.
///
/// The solution stores both the detailed moves and the per-move time deltas so
/// exact and heuristic procedures can extend, rollback, and score sequences
/// without recomputing travel times from scratch.
pub struct Solution {
    pub moves: Vec<DetailedMove>,
    time_deltas: Vec<f64>,
    pub crane_time: f64,
    pub acc: usize,
}

impl Solution {
    /// Creates an empty solution.
    pub fn new() -> Self {
        Self {
            moves: Vec::new(),
            time_deltas: Vec::new(),
            crane_time: 0.0,
            acc: 0,
        }
    }

    /// Reserves storage for additional moves and time deltas.
    pub fn reserve(&mut self, additional: usize) {
        self.moves.reserve(additional);
        self.time_deltas.reserve(additional);
    }

    /// Appends a move and adds its crane-time contribution.
    pub fn push(&mut self, dm: DetailedMove, time_delta: f64) {
        self.moves.push(dm);
        self.time_deltas.push(time_delta);
        self.crane_time += time_delta;
    }

    /// Removes the last move, if any, and subtracts its crane-time contribution.
    pub fn pop(&mut self) -> Option<(DetailedMove, f64)> {
        if let Some(dm) = self.moves.pop() {
            let delta = self.time_deltas.pop().unwrap();
            self.crane_time -= delta;
            Some((dm, delta))
        } else {
            None
        }
    }

    /// Returns the previous destination stack in one-based notation.
    ///
    /// `0` denotes the crane's initial safety position.
    pub fn last_dst_1based(&self) -> usize {
        self.moves.last().map_or(0, |dm| dm.dst + 1)
    }

    /// Returns the number of relocations in the solution.
    pub fn len(&self) -> usize {
        self.moves.len()
    }

    /// Returns whether the solution contains no relocations.
    pub fn is_empty(&self) -> bool {
        self.moves.is_empty()
    }

    /// Returns the crane-time delta of move `i`.
    pub fn time_delta(&self, i: usize) -> f64 {
        self.time_deltas[i]
    }
}
