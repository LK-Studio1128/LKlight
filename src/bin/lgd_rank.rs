/// lgd_rank: Reads all swarm GSO output files, sorts solutions by scoring and
/// luciferin, writes ranking files (rank_by_scoring.list, rank_by_luciferin.list).

use std::env;
use std::fs;
use std::io::{BufRead, BufReader, Write};
use std::path::PathBuf;

#[derive(Debug, Clone)]
struct Solution {
    id_swarm: u32,
    id_glowworm: u32,
    luciferin: f64,
    scoring: f64,
    coords: String,
    pdb_file: String,
}

fn parse_gso_file(path: &str, swarm_id: u32) -> Vec<Solution> {
    let file = match fs::File::open(path) {
        Ok(f) => f,
        Err(_) => return Vec::new(),
    };
    let reader = BufReader::new(file);
    let mut solutions = Vec::new();
    let mut id_glowworm = 0u32;

    for line in reader.lines() {
        let line = line.unwrap_or_default();
        let trimmed = line.trim();
        if trimmed.starts_with('#') || trimmed.is_empty() {
            continue;
        }
        // Format: (tx, ty, tz, qw, qx, qy, qz[, ...])  RecID LigID Luciferin NNeigh VRange Scoring
        if let Some(paren_end) = trimmed.find(')') {
            let coords = &trimmed[..=paren_end];
            let rest = trimmed[paren_end + 1..].split_whitespace().collect::<Vec<_>>();
            if rest.len() >= 6 {
                let luciferin: f64 = rest[2].parse().unwrap_or(0.0);
                let scoring: f64 = rest[5].parse().unwrap_or(0.0);
                solutions.push(Solution {
                    id_swarm: swarm_id,
                    id_glowworm,
                    luciferin,
                    scoring,
                    coords: coords.to_string(),
                    pdb_file: format!("lightdock_{}.pdb", id_glowworm),
                });
            }
        }
        id_glowworm += 1;
    }
    solutions
}

fn write_ranking(solutions: &[Solution], filename: &str, sort_by_scoring: bool) {
    let mut sorted = solutions.to_vec();
    if sort_by_scoring {
        sorted.sort_by(|a, b| b.scoring.partial_cmp(&a.scoring).unwrap_or(std::cmp::Ordering::Equal));
    } else {
        sorted.sort_by(|a, b| b.luciferin.partial_cmp(&a.luciferin).unwrap_or(std::cmp::Ordering::Equal));
    }
    let mut file = fs::File::create(filename).expect("Cannot create ranking file");
    writeln!(file, "# Rank Swarm GlowwormID Luciferin Scoring PDBFile Coordinates").unwrap();
    for (rank, sol) in sorted.iter().enumerate() {
        writeln!(
            file,
            "{:5} {:4} {:5} {:12.6} {:12.6} swarm_{}/{} {}",
            rank + 1,
            sol.id_swarm,
            sol.id_glowworm,
            sol.luciferin,
            sol.scoring,
            sol.id_swarm,
            sol.pdb_file,
            sol.coords
        ).unwrap();
    }
    println!("Wrote {} solutions to {}", sorted.len(), filename);
}

fn main() {
    let args: Vec<String> = env::args().collect();
    if args.len() < 3 {
        eprintln!("Usage: {} <num_swarms> <steps> [--result_file <file>]", args[0]);
        eprintln!("  Reads swarm_N/gso_STEPS.out for N=0..num_swarms-1");
        eprintln!("  Writes rank_by_scoring.list and rank_by_luciferin.list");
        std::process::exit(1);
    }

    let num_swarms: u32 = args[1].parse().expect("num_swarms must be integer");
    let steps: u32 = args[2].parse().expect("steps must be integer");

    let result_file_override: Option<String> = {
        let mut v = None;
        let mut i = 3;
        while i < args.len() {
            if (args[i] == "--result_file" || args[i] == "-f") && i + 1 < args.len() {
                v = Some(args[i + 1].clone());
                i += 2;
            } else {
                i += 1;
            }
        }
        v
    };

    let mut all_solutions: Vec<Solution> = Vec::new();

    for swarm_id in 0..num_swarms {
        let result_filename = match &result_file_override {
            Some(name) => format!("swarm_{}/{}", swarm_id, name),
            None => format!("swarm_{}/gso_{}.out", swarm_id, steps),
        };

        // Check for cluster representatives file
        let cluster_file = format!("swarm_{}/cluster_representatives.file", swarm_id);
        let cluster_ids: Option<Vec<u32>> = if PathBuf::from(&cluster_file).exists() {
            let content = fs::read_to_string(&cluster_file).unwrap_or_default();
            let ids = content.lines()
                .filter_map(|l| l.split(':').nth(3).and_then(|s| s.parse::<u32>().ok()))
                .collect();
            Some(ids)
        } else {
            None
        };

        let solutions = parse_gso_file(&result_filename, swarm_id);
        if solutions.is_empty() {
            eprintln!("Warning: {} not found or empty, skipping", result_filename);
            continue;
        }

        for sol in solutions {
            if let Some(ref ids) = cluster_ids {
                if ids.contains(&sol.id_glowworm) {
                    all_solutions.push(sol);
                }
            } else {
                all_solutions.push(sol);
            }
        }
    }

    println!("Total solutions collected: {}", all_solutions.len());
    write_ranking(&all_solutions, "rank_by_scoring.list", true);
    write_ranking(&all_solutions, "rank_by_luciferin.list", false);
    println!("Done.");
}
