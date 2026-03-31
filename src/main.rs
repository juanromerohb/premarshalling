//! Minimal CLI entrypoint that runs the default BBLCT configuration on one instance.

use lct_rust::bblct;
use lct_rust::crane::CraneParams;
use lct_rust::instance;

fn main() {
    let args: Vec<String> = std::env::args().collect();

    if args.len() < 2 {
        eprintln!("Usage: lct-rust <instance_path> [tau_minutes]");
        eprintln!(
            "  instance_path: path to an instance file (e.g., instances/CV/CVh5s3/h5s3n1.txt)"
        );
        eprintln!("  tau_minutes:   crane time limit in minutes (default: 15)");
        std::process::exit(1);
    }

    let path = &args[1];
    let tau_minutes: f64 = args
        .get(2)
        .map_or(15.0, |s| s.parse().expect("invalid tau"));
    let tau = tau_minutes * 60.0;
    let bay = instance::read_instance(path);
    let crane = CraneParams::new(bay.h);

    println!("Instance: {}", path);
    println!("S={}, H={}, C={}, P={}", bay.s, bay.h, bay.c, bay.p);
    println!("tau = {} min ({} s)", tau_minutes, tau);
    println!("Initial acc = {}", bay.acc());
    println!();
    println!("Bay:");
    println!("{}", bay.display());
    println!();

    println!("Running BBLCT...");
    let start = std::time::Instant::now();
    let (solution, stats) = bblct::bblct(&bay, &crane, tau, &bblct::BblctConfig::default());
    let elapsed = start.elapsed();

    println!("Done in {:.3} s", elapsed.as_secs_f64());
    println!(
        "Accessible: {} / {} ({:.1}%)",
        solution.acc,
        bay.c,
        100.0 * solution.acc as f64 / bay.c as f64
    );
    println!("Moves: {}", solution.len());
    println!("Crane time: {:.2} s", solution.crane_time);
    println!("Nodes explored: {}", stats.nodes_explored);
}
