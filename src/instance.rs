//! Instance loading and benchmark-dataset helpers.

use std::collections::BTreeSet;
use std::fs;
use std::path::{Path, PathBuf};

use crate::bay::{Bay, Group};
use crate::crane::CraneParams;

/// Fully specified instance including its crane-time limit and crane parameters.
pub struct Instance {
    pub name: String,
    pub bay: Bay,
    pub tau: f64,
    pub crane: CraneParams,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
/// Named benchmark families supported by the repository layout.
pub enum Benchmark {
    CV,
    BZ,
    ZJY,
    EMM,
}

impl Benchmark {
    /// Returns the default on-disk directory for the benchmark family.
    pub fn default_dir(&self) -> &'static str {
        match self {
            Benchmark::CV => "instances/CV",
            Benchmark::BZ => "instances/BZ",
            Benchmark::ZJY => "instances/ZJY",
            Benchmark::EMM => "instances/EMM",
        }
    }
}

/// Reads one text instance file into a [`Bay`].
pub fn read_instance(path: &str) -> Bay {
    let content = fs::read_to_string(path).expect("cannot read instance file");
    let lines: Vec<&str> = content.lines().collect();
    let h = lines.len();

    let rows: Vec<Vec<u8>> = lines
        .iter()
        .map(|line| {
            line.split_whitespace()
                .map(|x| x.parse::<u8>().expect("invalid number"))
                .collect()
        })
        .collect();

    let s = rows[0].len();

    let mut stacks: Vec<Vec<Group>> = vec![Vec::new(); s];
    for file_row in (0..h).rev() {
        for col in 0..s {
            let g = rows[file_row][col];
            if g > 0 {
                stacks[col].push(g);
            }
        }
    }

    Bay::from_vecs(s, h, &stacks)
}

/// Reads all instances from `dir` whose filenames encode category `(h, s)`.
pub fn read_category(dir: &str, h: usize, s: usize) -> Vec<(String, Bay)> {
    let pattern = format!("h{}s{}", h, s);
    let mut instances = Vec::new();
    let mut paths = Vec::new();
    collect_txt_paths(Path::new(dir), &mut paths);

    for path in paths {
        let filename = path.file_name().unwrap().to_string_lossy().to_string();
        if filename.contains(&pattern) {
            let path_str = path.to_string_lossy().to_string();
            let bay = read_instance(&path_str);
            let name = instance_display_name(&path);
            instances.push((name, bay));
        }
    }

    instances.sort_by(|a, b| a.0.cmp(&b.0));
    instances
}

/// Reads every `.txt` instance found recursively under `dir`.
pub fn read_all_instances(dir: &str) -> Vec<(String, Bay)> {
    let mut instances = Vec::new();
    let mut paths = Vec::new();
    collect_txt_paths(Path::new(dir), &mut paths);

    for path in paths {
        let name = instance_display_name(&path);
        let path_str = path.to_string_lossy().to_string();
        let bay = read_instance(&path_str);
        instances.push((name, bay));
    }

    instances.sort_by(|a, b| a.0.cmp(&b.0));
    instances
}

/// Lists all `(H, S)` categories detected under `dir`.
pub fn list_categories(dir: &str) -> Vec<(usize, usize)> {
    let mut categories = BTreeSet::new();
    let mut paths = Vec::new();
    collect_txt_paths(Path::new(dir), &mut paths);

    for path in paths {
        let name = path.file_name().unwrap().to_string_lossy().to_string();
        if let Some((h, s)) = try_parse_hs(&name) {
            categories.insert((h, s));
        }
    }

    categories.into_iter().collect()
}

fn collect_txt_paths(dir: &Path, out: &mut Vec<PathBuf>) {
    let entries = fs::read_dir(dir).expect("cannot read instance directory");
    for entry in entries {
        let entry = entry.expect("cannot read dir entry");
        let path = entry.path();
        if path.is_dir() {
            collect_txt_paths(&path, out);
        } else if path
            .extension()
            .and_then(|ext| ext.to_str())
            .is_some_and(|ext| ext.eq_ignore_ascii_case("txt"))
        {
            out.push(path);
        }
    }
}

fn instance_display_name(path: &Path) -> String {
    let filename = path.file_name().unwrap().to_string_lossy().to_string();

    let parent = path
        .parent()
        .and_then(|p| p.file_name())
        .and_then(|s| s.to_str())
        .unwrap_or("");

    if parent.starts_with("BF") && parent[2..].chars().all(|c| c.is_ascii_digit()) {
        return format!("{}_{}", parent, filename);
    }

    filename
}

/// Parses the `(H, S)` category encoded in a benchmark filename.
pub fn parse_hs_from_filename(name: &str) -> (usize, usize) {
    try_parse_hs(name).expect("cannot parse H and S from filename")
}

/// Tries to parse the `(H, S)` category encoded in a benchmark filename.
pub fn try_parse_hs(name: &str) -> Option<(usize, usize)> {
    let stem = Path::new(name).file_stem()?.to_str()?;
    let h_pos = stem.rfind('h')? + 1;
    let rest = &stem[h_pos..];
    let s_pos_rel = rest.find('s')?;
    let h: usize = rest[..s_pos_rel].parse().ok()?;
    let after_s = &rest[s_pos_rel + 1..];
    let s_len = after_s.chars().take_while(|c| c.is_ascii_digit()).count();
    if s_len == 0 {
        return None;
    }
    let s: usize = after_s[..s_len].parse().ok()?;
    Some((h, s))
}

/// Parses the `(P, F, H, S)` metadata encoded in a BZ benchmark filename.
pub fn parse_bz_params(name: &str) -> (usize, usize, usize, usize) {
    let stem = Path::new(name)
        .file_stem()
        .expect("no stem")
        .to_str()
        .expect("invalid name");

    if let Some(rest) = stem.strip_prefix('h') {
        let (h_str, rest) = rest.split_once('s').expect("invalid BZ name (missing 's')");
        let (s_str, rest) = rest.split_once('f').expect("invalid BZ name (missing 'f')");
        let (f_str, rest) = rest.split_once('p').expect("invalid BZ name (missing 'p')");
        let (p_str, _) = rest.split_once('n').expect("invalid BZ name (missing 'n')");

        let h: usize = h_str.parse().expect("invalid H");
        let s: usize = s_str.parse().expect("invalid S");
        let f: usize = f_str.parse().expect("invalid F");
        let p: usize = p_str.parse().expect("invalid P");
        return (p, f, h, s);
    }

    if let Some(rest) = stem.strip_prefix('p') {
        let (p_str, rest) = rest.split_once('f').expect("invalid BZ name (missing 'f')");
        let (f_str, rest) = rest.split_once('h').expect("invalid BZ name (missing 'h')");
        let (h_str, rest) = rest.split_once('s').expect("invalid BZ name (missing 's')");
        let (s_str, _) = rest.split_once('n').expect("invalid BZ name (missing 'n')");

        let p: usize = p_str.parse().expect("invalid P");
        let f: usize = f_str.parse().expect("invalid F");
        let h: usize = h_str.parse().expect("invalid H");
        let s: usize = s_str.parse().expect("invalid S");
        return (p, f, h, s);
    }

    panic!("invalid BZ filename format")
}
