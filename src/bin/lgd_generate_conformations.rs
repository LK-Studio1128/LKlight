/// lgd_generate_conformations: Given a LightDock GSO output file, applies the
/// docking transformations to generate final PDB structure files.

use rayon::prelude::*;
use std::env;
use std::fs;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

// Quaternion for rotation
#[derive(Clone, Copy, Debug)]
struct Quaternion {
    w: f64,
    x: f64,
    y: f64,
    z: f64,
}

impl Quaternion {
    fn new(w: f64, x: f64, y: f64, z: f64) -> Self {
        Quaternion { w, x, y, z }
    }

    fn rotate(&self, v: [f64; 3]) -> [f64; 3] {
        let (w, qx, qy, qz) = (self.w, self.x, self.y, self.z);
        let (vx, vy, vz) = (v[0], v[1], v[2]);
        let tx = 2.0 * (qy * vz - qz * vy);
        let ty = 2.0 * (qz * vx - qx * vz);
        let tz = 2.0 * (qx * vy - qy * vx);
        [
            vx + w * tx + qy * tz - qz * ty,
            vy + w * ty + qz * tx - qx * tz,
            vz + w * tz + qx * ty - qy * tx,
        ]
    }
}

// Minimal PDB atom record
#[derive(Clone, Debug)]
struct Atom {
    serial: u32,
    name: String,
    alt_loc: char,
    res_name: String,
    chain_id: char,
    res_seq: i32,
    i_code: char,
    x: f64,
    y: f64,
    z: f64,
    occupancy: f64,
    temp_factor: f64,
    element: String,
    charge: String,
    is_hetatm: bool,
}

fn parse_pdb(filename: &str) -> Vec<Atom> {
    let file = fs::File::open(filename).expect(&format!("Cannot open PDB: {}", filename));
    let reader = BufReader::new(file);
    let mut atoms = Vec::new();
    for line in reader.lines() {
        let line = line.unwrap_or_default();
        if line.len() < 27 {
            continue;
        }
        let record = &line[..6];
        let is_atom = record.trim() == "ATOM";
        let is_hetatm = record.trim() == "HETATM";
        if !is_atom && !is_hetatm {
            continue;
        }
        let serial: u32 = line[6..11].trim().parse().unwrap_or(0);
        let name = line[12..16].trim().to_string();
        let alt_loc = line.chars().nth(16).unwrap_or(' ');
        let res_name = line[17..21].trim().to_string();
        let chain_id = line.chars().nth(21).unwrap_or('A');
        let res_seq: i32 = line[22..26].trim().parse().unwrap_or(0);
        let i_code = line.chars().nth(26).unwrap_or(' ');
        let x: f64 = line[30..38].trim().parse().unwrap_or(0.0);
        let y: f64 = line[38..46].trim().parse().unwrap_or(0.0);
        let z: f64 = line[46..54].trim().parse().unwrap_or(0.0);
        let occupancy: f64 = if line.len() > 60 { line[54..60].trim().parse().unwrap_or(1.0) } else { 1.0 };
        let temp_factor: f64 = if line.len() > 66 { line[60..66].trim().parse().unwrap_or(0.0) } else { 0.0 };
        let element = if line.len() > 78 { line[76..78].trim().to_string() } else { String::new() };
        let charge = if line.len() > 80 { line[78..80].trim().to_string() } else { String::new() };
        atoms.push(Atom {
            serial,
            name,
            alt_loc,
            res_name,
            chain_id,
            res_seq,
            i_code,
            x, y, z,
            occupancy,
            temp_factor,
            element,
            charge,
            is_hetatm,
        });
    }
    atoms
}

fn write_pdb(atoms: &[Atom], filename: &str) {
    let file = fs::File::create(filename).expect(&format!("Cannot create {}", filename));
    let mut file = BufWriter::new(file);
    for atom in atoms {
        let record = if atom.is_hetatm { "HETATM" } else { "ATOM  " };
        writeln!(
            file,
            "{}{:5} {:<4}{}{:<4}{}{:4}{}   {:8.3}{:8.3}{:8.3}{:6.2}{:6.2}          {:>2}{:>2}",
            record,
            atom.serial,
            atom.name,
            atom.alt_loc,
            atom.res_name,
            atom.chain_id,
            atom.res_seq,
            atom.i_code,
            atom.x, atom.y, atom.z,
            atom.occupancy,
            atom.temp_factor,
            atom.element,
            atom.charge,
        ).unwrap();
    }
    writeln!(file, "END").unwrap();
}

// Parse setup.json for ANM configuration
fn read_setup_anm(setup_file: &str) -> (bool, usize, usize) {
    let content = fs::read_to_string(setup_file).unwrap_or_default();
    let use_anm = content.contains("\"use_anm\": true");
    let anm_rec = extract_usize_field(&content, "anm_rec").unwrap_or(10);
    let anm_lig = extract_usize_field(&content, "anm_lig").unwrap_or(10);
    (use_anm, anm_rec, anm_lig)
}

fn extract_usize_field(content: &str, field: &str) -> Option<usize> {
    let pattern = format!("\"{}\":", field);
    let pos = content.find(&pattern)?;
    let rest = &content[pos + pattern.len()..];
    let trimmed = rest.trim_start();
    let end = trimmed.find(|c: char| !c.is_ascii_digit())?;
    trimmed[..end].parse().ok()
}

#[derive(Debug)]
struct GSOEntry {
    translation: [f64; 3],
    rotation: Quaternion,
    rec_nmodes: Vec<f64>,
    lig_nmodes: Vec<f64>,
    luciferin: f64,
    scoring: f64,
}

/// Parse GSO output file. Lines look like:
/// (tx, ty, tz, qw, qx, qy, qz[, nm...])  RecID LigID Luciferin NNeigh VRange Score
fn parse_gso_output(filename: &str, num_anm_rec: usize, num_anm_lig: usize) -> Vec<GSOEntry> {
    let file = fs::File::open(filename).expect(&format!("Cannot open GSO file: {}", filename));
    let reader = BufReader::new(file);
    let mut entries = Vec::new();

    for line in reader.lines() {
        let line = line.unwrap_or_default();
        let trimmed = line.trim();
        if trimmed.starts_with('#') || trimmed.is_empty() {
            continue;
        }
        if !trimmed.starts_with('(') {
            continue;
        }
        let paren_end = match trimmed.find(')') {
            Some(p) => p,
            None => continue,
        };
        let inner = &trimmed[1..paren_end];
        let values: Vec<f64> = inner.split(',')
            .filter_map(|s| s.trim().parse::<f64>().ok())
            .collect();
        if values.len() < 7 {
            continue;
        }
        let translation = [values[0], values[1], values[2]];
        let rotation = Quaternion::new(values[3], values[4], values[5], values[6]);

        let mut rec_nmodes = Vec::new();
        let mut lig_nmodes = Vec::new();
        if num_anm_rec > 0 && values.len() >= 7 + num_anm_rec {
            rec_nmodes = values[7..7 + num_anm_rec].to_vec();
        }
        if num_anm_lig > 0 && values.len() >= 7 + num_anm_rec + num_anm_lig {
            lig_nmodes = values[7 + num_anm_rec..7 + num_anm_rec + num_anm_lig].to_vec();
        }

        // Parse rest of line
        let rest = trimmed[paren_end + 1..].split_whitespace().collect::<Vec<_>>();
        let luciferin: f64 = rest.get(2).and_then(|s| s.parse().ok()).unwrap_or(0.0);
        let scoring: f64 = rest.get(5).and_then(|s| s.parse().ok()).unwrap_or(0.0);

        entries.push(GSOEntry { translation, rotation, rec_nmodes, lig_nmodes, luciferin, scoring });
    }
    entries
}

/// Parse initial positions .dat file (no parentheses, space-separated)
fn parse_dat_file(filename: &str, num_anm_rec: usize, num_anm_lig: usize) -> Vec<GSOEntry> {
    let content = fs::read_to_string(filename).expect(&format!("Cannot read {}", filename));
    let mut entries = Vec::new();
    for line in content.lines() {
        let trimmed = line.trim();
        if trimmed.starts_with('#') || trimmed.is_empty() {
            continue;
        }
        let values: Vec<f64> = trimmed.split_whitespace()
            .filter_map(|s| s.parse::<f64>().ok())
            .collect();
        if values.len() < 7 {
            continue;
        }
        let translation = [values[0], values[1], values[2]];
        let rotation = Quaternion::new(values[3], values[4], values[5], values[6]);
        let rec_nmodes = if num_anm_rec > 0 && values.len() >= 7 + num_anm_rec {
            values[7..7 + num_anm_rec].to_vec()
        } else { Vec::new() };
        let lig_nmodes = if num_anm_lig > 0 && values.len() >= 7 + num_anm_rec + num_anm_lig {
            values[7 + num_anm_rec..7 + num_anm_rec + num_anm_lig].to_vec()
        } else { Vec::new() };
        entries.push(GSOEntry { translation, rotation, rec_nmodes, lig_nmodes, luciferin: 0.0, scoring: 0.0 });
    }
    entries
}

fn main() {
    let args: Vec<String> = env::args().collect();
    if args.len() < 5 {
        eprintln!("Usage: {} <lightdock_receptor.pdb> <lightdock_ligand.pdb> <gso_output_or_dat> <num_conformations> [--setup setup.json]",
                  args[0]);
        std::process::exit(1);
    }

    let receptor_file = &args[1];
    let ligand_file = &args[2];
    let output_file = &args[3];
    let num_conformations: usize = args[4].parse().expect("num_conformations must be integer");

    let setup_file = {
        let mut s = None;
        let mut i = 5;
        while i < args.len() {
            if (args[i] == "--setup" || args[i] == "-s") && i + 1 < args.len() {
                s = Some(args[i + 1].clone());
                i += 2;
            } else {
                i += 1;
            }
        }
        s
    };

    let (use_anm, num_anm_rec, num_anm_lig) = match &setup_file {
        Some(sf) => read_setup_anm(sf),
        None => (false, 0, 0),
    };

    println!("Reading receptor: {}", receptor_file);
    let receptor_atoms = parse_pdb(receptor_file);
    println!("  {} atoms", receptor_atoms.len());

    println!("Reading ligand: {}", ligand_file);
    let ligand_atoms = parse_pdb(ligand_file);
    println!("  {} atoms", ligand_atoms.len());

    // Parse entries
    let entries = if output_file.ends_with(".dat") {
        parse_dat_file(output_file, if use_anm { num_anm_rec } else { 0 }, if use_anm { num_anm_lig } else { 0 })
    } else {
        parse_gso_output(output_file, if use_anm { num_anm_rec } else { 0 }, if use_anm { num_anm_lig } else { 0 })
    };

    let n = num_conformations.min(entries.len());
    if n < num_conformations {
        eprintln!("Warning: only {} entries found (requested {}), clipping", n, num_conformations);
    }

    let destination = Path::new(output_file).parent().unwrap_or(Path::new("."));

    // Parallel: each conformation is independent (clone + rotate/translate + write).
    // Safe — this binary is invoked as a one-shot post-processing step, not nested in GSO rayon.
    entries.par_iter().take(n).enumerate().for_each(|(i, entry)| {
        let mut ligand_pose = ligand_atoms.clone();
        // Rotate then translate ligand
        for atom in &mut ligand_pose {
            let rotated = entry.rotation.rotate([atom.x, atom.y, atom.z]);
            atom.x = rotated[0] + entry.translation[0];
            atom.y = rotated[1] + entry.translation[1];
            atom.z = rotated[2] + entry.translation[2];
        }

        // Combine receptor + ligand into one PDB file
        let out_path = destination.join(format!("lightdock_{}.pdb", i));
        let out_str = out_path.to_str().unwrap();

        let mut all_atoms = receptor_atoms.clone();
        // Re-serial ligand atoms
        let rec_count = receptor_atoms.len() as u32;
        for (j, atom) in ligand_pose.iter_mut().enumerate() {
            atom.serial = rec_count + j as u32 + 1;
        }
        all_atoms.extend(ligand_pose);
        write_pdb(&all_atoms, out_str);
    });
    println!("Generated {} conformations", n);
}
