//! Compares the crane-time lower bound at the initial layout with the BBLCT solutions.
//!
//! For each run of a BBLCT summary CSV, computes UB(tau) and LB(n) at the initial layout, where
//! n is the number of accessible containers of the solution, and writes the per-run values and a
//! LaTeX table of their averages. It also reports any run with n > UB(tau) or LB(n) > t, which
//! would contradict the validity of the bound.
//!
//! Usage: bound_table --results <all_summary.csv> --instances <dir> --out-csv <csv> --out-tex <tex>

use lct_rust::bay::Bay;
use lct_rust::crane::CraneParams;
use lct_rust::instance;
use lct_rust::lower_bound::lower_bound;
use std::collections::{BTreeMap, HashMap};
use std::fmt::Write as _;
use std::fs;
use std::path::Path;

const FAMILIES: [&str; 5] = ["CV", "BF", "BZ", "EMM", "ZJY"];
const TAUS: [u32; 3] = [600, 900, 1200];
/// The crane times of the summary CSV are rounded to two decimals.
const TIME_TOLERANCE: f64 = 0.01;

struct Run {
    instance: String,
    family: String,
    h: usize,
    s: usize,
    c: usize,
    tau: u32,
    initial_acc: usize,
    ub: usize,
    n: usize,
    lb_n: f64,
    t: f64,
    proved_optimal: bool,
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let results = arg_value(&args, "--results");
    let instances = arg_value(&args, "--instances");
    let out_csv = arg_value(&args, "--out-csv");
    let out_tex = arg_value(&args, "--out-tex");

    let paths = instance_paths(&instances);
    let runs = compute_runs(&results, &paths);

    report(&runs);
    fs::write(&out_csv, runs_csv(&runs)).expect("cannot write CSV");
    fs::write(&out_tex, latex_table(&runs)).expect("cannot write LaTeX table");
}

fn arg_value(args: &[String], name: &str) -> String {
    let pos = args.iter().position(|a| a == name);
    match pos.and_then(|i| args.get(i + 1)) {
        Some(value) => value.clone(),
        None => panic!("missing argument {}", name),
    }
}

/// Maps the instance names of the summary CSV (CV_h5s3n1.txt, BF_BF10_h8s16n1.txt) to their files.
fn instance_paths(root: &str) -> HashMap<String, String> {
    let mut paths = HashMap::new();
    for family in FAMILIES {
        let family_dir = Path::new(root).join(family);
        for sub in fs::read_dir(&family_dir).expect("cannot read instance directory") {
            let sub = sub.expect("cannot read instance directory").path();
            let sub_name = sub.file_name().unwrap().to_string_lossy().to_string();
            for file in fs::read_dir(&sub).expect("cannot read instance directory") {
                let file = file.expect("cannot read instance directory").path();
                let file_name = file.file_name().unwrap().to_string_lossy().to_string();
                let key = if family == "BF" {
                    format!("{}_{}_{}", family, sub_name, file_name)
                } else {
                    format!("{}_{}", family, file_name)
                };
                paths.insert(key, file.to_string_lossy().to_string());
            }
        }
    }
    paths
}

fn compute_runs(results: &str, paths: &HashMap<String, String>) -> Vec<Run> {
    let content = fs::read_to_string(results).expect("cannot read results CSV");
    let mut lines = content.lines();
    let header: Vec<&str> = lines.next().expect("empty results CSV").split(',').collect();
    let column = |name: &str| {
        header
            .iter()
            .position(|h| *h == name)
            .unwrap_or_else(|| panic!("missing column {}", name))
    };
    let (col_instance, col_tau, col_acc) = (column("instance"), column("tau"), column("acc"));
    let (col_time, col_optimal) = (column("crane_time"), column("proved_optimal"));
    let (col_c, col_initial) = (column("c"), column("initial_acc"));

    let mut runs = Vec::new();
    for line in lines {
        let f: Vec<&str> = line.split(',').collect();
        let name = f[col_instance];
        let path = paths
            .get(name)
            .unwrap_or_else(|| panic!("no instance file for {}", name));
        let bay = instance::read_instance(path);
        let crane = CraneParams::new(bay.h);
        let tau: f64 = f[col_tau].parse().expect("invalid tau");
        let n: usize = f[col_acc].parse().expect("invalid acc");

        assert_eq!(bay.c.to_string(), f[col_c], "c differs for {}", name);
        assert_eq!(
            bay.acc().to_string(),
            f[col_initial],
            "initial_acc differs for {}",
            name
        );

        runs.push(Run {
            instance: name.to_string(),
            family: name.split('_').next().unwrap().to_string(),
            h: bay.h,
            s: bay.s,
            c: bay.c,
            tau: tau as u32,
            initial_acc: bay.acc(),
            ub: upper_bound(&bay, &crane, tau),
            n,
            lb_n: lower_bound(&bay, &crane, n),
            t: f[col_time].parse().expect("invalid crane time"),
            proved_optimal: f[col_optimal] == "true",
        });
    }
    runs
}

/// UB(tau) = max { m in {0, ..., C} : LB(m) <= tau }.
fn upper_bound(bay: &Bay, crane: &CraneParams, tau: f64) -> usize {
    (0..=bay.c)
        .rev()
        .find(|&m| lower_bound(bay, crane, m) <= tau)
        .unwrap_or(0)
}

fn report(runs: &[Run]) {
    let n_above_ub = runs.iter().filter(|r| r.n > r.ub).count();
    let lb_above_t = runs
        .iter()
        .filter(|r| r.lb_n > r.t + TIME_TOLERANCE)
        .count();
    eprintln!("runs: {}", runs.len());
    eprintln!("n > UB(tau): {}", n_above_ub);
    eprintln!("LB(n) > t: {}", lb_above_t);
    for r in runs.iter().filter(|r| r.n > r.ub || r.lb_n > r.t + TIME_TOLERANCE) {
        eprintln!(
            "  {} tau={} UB={} n={} LB(n)={:.2} t={:.2}",
            r.instance, r.tau, r.ub, r.n, r.lb_n, r.t
        );
    }
}

fn runs_csv(runs: &[Run]) -> String {
    let mut out = String::from("instance,h,s,c,tau,initial_acc,ub,n,lb_n,t,proved_optimal\n");
    for r in runs {
        writeln!(
            out,
            "{},{},{},{},{},{},{},{},{:.2},{:.2},{}",
            r.instance, r.h, r.s, r.c, r.tau, r.initial_acc, r.ub, r.n, r.lb_n, r.t, r.proved_optimal
        )
        .unwrap();
    }
    out
}

fn mean(values: impl Iterator<Item = f64>) -> f64 {
    let (mut sum, mut count) = (0.0, 0usize);
    for value in values {
        sum += value;
        count += 1;
    }
    sum / count as f64
}

/// Returns the row cells: number of instances, then UB(tau), n, LB(n) and t for each tau.
fn row_cells(runs: &[&Run]) -> String {
    let mut cells = vec![(runs.len() / TAUS.len()).to_string()];
    for tau in TAUS {
        let at_tau = || runs.iter().filter(move |r| r.tau == tau);
        cells.push(format!("{:.2}", mean(at_tau().map(|r| r.ub as f64))));
        cells.push(format!("{:.2}", mean(at_tau().map(|r| r.n as f64))));
        cells.push(format!("{:.1}", mean(at_tau().map(|r| r.lb_n))));
        cells.push(format!("{:.1}", mean(at_tau().map(|r| r.t))));
    }
    cells.join(" & ")
}

fn latex_table(runs: &[Run]) -> String {
    let mut out = String::new();
    out.push_str(
        r"\begin{table}[p]
  \centering
  \caption{Comparison of the bounds at the initial layout with the BBLCT solutions over the five benchmark families. For each crane-time budget \( \tau \in \{600,900,1200\} \) seconds, the table reports the average upper bound on the number of accessible containers, \( UB(\tau) \); the average number of accessible containers in the optimal or best solution found, \( n \); the average lower bound on the crane time required to attain them, \( \operatorname{LB}(n) \), in seconds; and the average crane time of that solution, \( t \), in seconds.}
  \label{tab:bounds}
  \resizebox{\linewidth}{!}{
  \begin{tabular}{@{} c r r r @{\hspace{\bbGroupSep}} *{4}{r} @{\hspace{\bbGroupSep}} *{4}{r} @{\hspace{\bbGroupSep}} *{4}{r} @{}}
    \toprule
    & & & &
    \multicolumn{4}{c@{\hspace{\bbGroupSep}}}{\( \tau = 600 \) s} &
    \multicolumn{4}{c@{\hspace{\bbGroupSep}}}{\( \tau = 900 \) s} &
    \multicolumn{4}{c}{\( \tau = 1200 \) s} \\
    \cmidrule(l{0pt}r{\bbGroupSep}){5-8} \cmidrule(l{0pt}r{\bbGroupSep}){9-12} \cmidrule(l{0pt}r{0pt}){13-16}
    Set & \( H \) & \( S \) & \# &
    \( UB(\tau) \) & \( n \) & \( \operatorname{LB}(n) \) & \( t \) &
    \( UB(\tau) \) & \( n \) & \( \operatorname{LB}(n) \) & \( t \) &
    \( UB(\tau) \) & \( n \) & \( \operatorname{LB}(n) \) & \( t \) \\
    \midrule
",
    );

    for (i, family) in FAMILIES.iter().enumerate() {
        let of_family: Vec<&Run> = runs.iter().filter(|r| r.family == *family).collect();
        let mut categories: BTreeMap<(usize, usize), Vec<&Run>> = BTreeMap::new();
        for r in &of_family {
            categories.entry((r.h, r.s)).or_default().push(r);
        }

        for (j, ((h, s), category)) in categories.iter().enumerate() {
            let label = if j == 0 {
                format!(r"\multirow{{{}}}{{*}}{{{}}}", categories.len() + 1, family)
            } else {
                String::new()
            };
            writeln!(out, r"    {} & {} & {} & {} \\", label, h, s, row_cells(category)).unwrap();
        }
        out.push_str("    \\cmidrule{2-16}\n");
        writeln!(
            out,
            r"     & \multicolumn{{2}}{{c}}{{Total}} & {} \\",
            row_cells(&of_family)
        )
        .unwrap();
        if i + 1 < FAMILIES.len() {
            out.push_str("    \\cmidrule{1-16}\n");
        }
    }

    out.push_str("    \\bottomrule\n  \\end{tabular}}\n\\end{table}\n");
    out
}
