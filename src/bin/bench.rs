//! CLI entrypoint for batch benchmark runs.

use lct_rust::bblct::{BblctConfig, HeuristicChoice};
use lct_rust::benchmark::{self, Algorithm, BenchConfig, ModelWarmStartHeuristic};

fn main() {
    let args: Vec<String> = std::env::args().collect();

    let mut instance_dir = "instances/CV".to_string();
    let mut instance_file: Option<String> = None;
    let mut algo_name = "all".to_string();
    let mut tau_minutes = 30.0f64;
    let mut max_per_category: Option<usize> = None;
    let mut bblct_cpu_limit_sec: Option<f64> = None;
    let mut cp_time_limit_sec: f64 = 3600.0;
    let mut ip_time_limit_sec: f64 = 3600.0;
    let mut model_warm_start: Option<ModelWarmStartHeuristic> = None;
    let mut output_dir: Option<String> = None;
    let mut csv_path: Option<String> = None;
    let mut save_solutions = false;
    let mut save_results = false;
    let mut run_tag: Option<String> = None;
    let mut skip_category_analysis = false;

    let mut bblct_config = BblctConfig::default();

    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "-d" | "--dir" => {
                i += 1;
                instance_dir = arg_value(&args, i, "--dir").to_string();
            }
            "--instance-file" => {
                i += 1;
                instance_file = Some(arg_value(&args, i, "--instance-file").to_string());
            }
            "-a" | "--algo" => {
                i += 1;
                algo_name = arg_value(&args, i, "--algo").to_string();
            }
            "-t" | "--tau" => {
                i += 1;
                tau_minutes = arg_value(&args, i, "--tau").parse().expect("invalid tau");
            }
            "-n" | "--max" => {
                i += 1;
                max_per_category = Some(arg_value(&args, i, "--max").parse().expect("invalid max"));
            }
            "--bblct-cpu-limit-sec" => {
                i += 1;
                bblct_cpu_limit_sec = Some(
                    arg_value(&args, i, "--bblct-cpu-limit-sec")
                        .parse()
                        .expect("invalid --bblct-cpu-limit-sec"),
                );
            }
            "--cp-time-limit" => {
                i += 1;
                cp_time_limit_sec = arg_value(&args, i, "--cp-time-limit")
                    .parse()
                    .expect("invalid cp-time-limit");
            }
            "--ip-time-limit" => {
                i += 1;
                ip_time_limit_sec = arg_value(&args, i, "--ip-time-limit")
                    .parse()
                    .expect("invalid ip-time-limit");
            }
            "--model-warm-start" => {
                i += 1;
                model_warm_start =
                    parse_model_warm_start(arg_value(&args, i, "--model-warm-start"));
            }
            "-o" | "--output" => {
                i += 1;
                output_dir = Some(arg_value(&args, i, "--output").to_string());
            }
            "--csv" => {
                i += 1;
                csv_path = Some(arg_value(&args, i, "--csv").to_string());
            }
            "--run-tag" => {
                i += 1;
                run_tag = Some(arg_value(&args, i, "--run-tag").to_string());
            }
            "--save-solutions" => {
                save_solutions = true;
            }
            "--save-results" => {
                save_results = true;
            }
            "--skip-category-analysis" => {
                skip_category_analysis = true;
            }
            "--bblct-init-heuristic" => {
                i += 1;
                bblct_config.initial_heuristic =
                    parse_heuristic(arg_value(&args, i, "--bblct-init-heuristic"));
            }
            "--bblct-node-heuristic" => {
                i += 1;
                bblct_config.node_heuristic =
                    parse_heuristic(arg_value(&args, i, "--bblct-node-heuristic"));
            }
            "--bblct-dominance" => {
                i += 1;
                bblct_config.use_dominance = parse_bool(arg_value(&args, i, "--bblct-dominance"));
            }
            "--bblct-lb-pruning" => {
                i += 1;
                bblct_config.use_lb_pruning = parse_bool(arg_value(&args, i, "--bblct-lb-pruning"));
            }
            "--bblct-lb-sorting" => {
                i += 1;
                bblct_config.use_lb_sorting = parse_bool(arg_value(&args, i, "--bblct-lb-sorting"));
            }
            "--bblct-node-heuristic-every" => {
                i += 1;
                bblct_config.node_heuristic_every =
                    arg_value(&args, i, "--bblct-node-heuristic-every")
                        .parse()
                        .expect("invalid --bblct-node-heuristic-every");
            }
            "-h" | "--help" => {
                print_usage();
                return;
            }
            other => {
                eprintln!("Unknown option: {}", other);
                print_usage();
                std::process::exit(1);
            }
        }
        i += 1;
    }

    let tau = tau_minutes * 60.0;

    let algorithms = if algo_name == "all" {
        Algorithm::all(cp_time_limit_sec, ip_time_limit_sec, &bblct_config)
    } else {
        algo_name
            .split(',')
            .map(|name| {
                Algorithm::from_name(
                    name.trim(),
                    cp_time_limit_sec,
                    ip_time_limit_sec,
                    &bblct_config,
                )
                .unwrap_or_else(|| panic!("Unknown algorithm: {}", name))
            })
            .collect()
    };

    if save_solutions && output_dir.is_none() {
        output_dir = Some("results".to_string());
    }

    let config = BenchConfig {
        instance_dir,
        instance_file,
        algorithms,
        tau,
        max_per_category,
        output_dir,
        csv_path: csv_path.clone(),
        save_solutions,
        save_results,
        run_tag,
        skip_category_analysis,
        bblct_cpu_limit_sec,
        model_warm_start,
    };

    let results = benchmark::run_benchmark(&config);
    benchmark::print_summary_table(&results);

    if let Some(ref path) = csv_path {
        benchmark::save_results_csv(&results, path);
    }
}

fn print_usage() {
    eprintln!("Usage: bench [options]");
    eprintln!();
    eprintln!("Options:");
    eprintln!("  -d, --dir <path>                 Instance directory (default: instances/CV)");
    eprintln!("      --instance-file <path>       Run a single instance file");
    eprintln!("  -a, --algo <name>                Algorithm list or all (default: all)");
    eprintln!("                                   Available: bblct, hlct, glct, tgh,");
    eprintln!("                                   lct1-cpo, lct2-cpo, lct1-cpsat, lct2-cpsat,");
    eprintln!("                                   iplct");
    eprintln!("  -t, --tau <minutes>              Crane time limit in minutes (default: 30)");
    eprintln!("  -n, --max <count>                Max instances per category (default: all)");
    eprintln!("      --bblct-cpu-limit-sec <sec>  Internal wall-time limit for BBLCT (optional)");
    eprintln!("      --cp-time-limit <sec>        Time limit for CP models (default: 3600)");
    eprintln!("      --ip-time-limit <sec>        Time limit for IPLCT model (default: 3600)");
    eprintln!("      --model-warm-start <h>       none|glct|hlct (applies to CP/IP models)");
    eprintln!("  -o, --output <dir>               Output directory for solution files");
    eprintln!("      --csv <path>                 Save summary CSV to this path");
    eprintln!("      --save-solutions             Save individual human-readable solution files");
    eprintln!("      --save-results               Save structured results");
    eprintln!("      --run-tag <tag>              Append run tag to result folders");
    eprintln!("      --skip-category-analysis     Disable category-level analysis generation");
    eprintln!("      --bblct-init-heuristic <h>   none|glct|hlct");
    eprintln!("      --bblct-node-heuristic <h>   none|glct|hlct");
    eprintln!("      --bblct-dominance <bool>     true|false");
    eprintln!("      --bblct-lb-pruning <bool>    true|false");
    eprintln!("      --bblct-lb-sorting <bool>    true|false");
    eprintln!("      --bblct-node-heuristic-every <N>");
    eprintln!("  -h, --help                       Show this help");
}

fn arg_value<'a>(args: &'a [String], i: usize, flag: &str) -> &'a str {
    args.get(i)
        .map(|s| s.as_str())
        .unwrap_or_else(|| panic!("missing value for {}", flag))
}

fn parse_bool(value: &str) -> bool {
    match value.to_ascii_lowercase().as_str() {
        "1" | "true" | "yes" | "y" => true,
        "0" | "false" | "no" | "n" => false,
        _ => panic!("invalid boolean value: {}", value),
    }
}

fn parse_heuristic(value: &str) -> HeuristicChoice {
    match value.to_ascii_lowercase().as_str() {
        "none" => HeuristicChoice::None,
        "glct" => HeuristicChoice::Glct,
        "hlct" => HeuristicChoice::Hlct,
        _ => panic!("invalid heuristic: {} (expected none|glct|hlct)", value),
    }
}

fn parse_model_warm_start(value: &str) -> Option<ModelWarmStartHeuristic> {
    match value.to_ascii_lowercase().as_str() {
        "none" => None,
        "glct" => Some(ModelWarmStartHeuristic::Glct),
        "hlct" => Some(ModelWarmStartHeuristic::Hlct),
        _ => panic!(
            "invalid --model-warm-start: {} (expected none|glct|hlct)",
            value
        ),
    }
}
