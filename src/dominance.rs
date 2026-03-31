//! Simple dominance checks used by BBLCT and related search code.

use crate::bay::Bay;
use crate::solution::Solution;

/// Returns whether moving from `src` would immediately undo the last destination.
pub fn violates_transitive(sol: &Solution, src: usize) -> bool {
    if let Some(last) = sol.moves.last() {
        src == last.dst
    } else {
        false
    }
}

/// Returns whether the move would recreate the same-group back-and-forth pattern.
pub fn violates_same_group(sol: &Solution, bay: &Bay, src: usize, dst: usize) -> bool {
    if let Some(last) = sol.moves.last() {
        if dst == last.src && !bay.is_empty(src) {
            bay.top_group(src) == last.group
        } else {
            false
        }
    } else {
        false
    }
}
