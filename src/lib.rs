//! Algorithms and tooling for the container premarshalling problem under limited crane time.
//!
//! The crate contains constructive heuristics, exact methods, model generators,
//! instance readers, reporting utilities, and a desktop GUI for running them.

pub mod bay;
pub mod bblct;
pub mod benchmark;
pub mod cplct;
pub mod crane;
pub mod display;
pub mod dominance;
pub mod gui;
pub mod heuristics;
pub mod instance;
pub mod iplct;
pub mod lower_bound;
pub mod results_output;
pub mod solution;
