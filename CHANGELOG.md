# Changelog

## v1.1.0 (2026-11-03)

This version corrects two implementation errors and replaces all the computational results.

### Fixed

- Lower bound: when computing the second component, the marginal relocations `re(c)` counted every container above `c` instead of starting at `b_c`, the lowest container that blocks `c`. This overestimated the bound on instances with repeated group values, so BBLCT could prune better solutions and report wrong optimal values.
- HLCT: the crane time of the GLCT completion was computed as if the crane started from its initial position instead of the end of the TGH prefix, so a solution could exceed the crane-time limit. Candidates whose total crane time exceeds the limit are now discarded.

### Changed

- `paper-results/` now holds the results of a complete re-run with Rust 1.99.1, OR-Tools CP-SAT 9.15.6755 and Gurobi 13.0.3. The previous results remain available at the `v1.0.0` tag.
- The average reduction of inaccessible containers excludes the 34 base instances whose initial layout is already fully accessible (1 in CV and 33 in EMM), since it is undefined for them.

### Added

- `scripts/make_category_tables.py`: builds the category-level tables from the summary CSV files.
- `src/bin/bound_table.rs` and `paper-results/bound-analysis/`: comparison of the bounds at the initial layout with the BBLCT solutions.
- `.zenodo.json`: archive metadata.

### Removed

- The desktop GUI, `gui.exe` and their dependencies. The results of the paper do not depend on them; the last version with the GUI is `v1.0.0`.
