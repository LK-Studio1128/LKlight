/// lklight-setup: Rust port of lightdock3_setup.py
///
/// Sets up a LKlight simulation:
///   1. Reads receptor and ligand PDB files
///   2. Translates both structures to origin
///   3. Saves lightdock_<receptor>.pdb and lightdock_<ligand>.pdb
///   4. Generates N swarm positions on receptor surface (Fibonacci sphere)
///   5. Populates each swarm with M random glowworm positions
///   6. Writes initial_positions_N.dat files and swarm_N/ directories
///   7. Writes setup.json for use with lklight

use pdbtbx::PDB;
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::env;
use std::fs;
use std::io::Write;
use std::path::Path;

// ─── Setup file structure (must match LKlight's SetupFile) ──────────────────

#[derive(Serialize, Deserialize, Debug)]
struct SetupFile {
    seed: Option<u64>,
    anm_seed: u64,
    ftdock_file: Option<String>,
    noh: bool,
    anm_rec: usize,
    anm_lig: usize,
    swarms: u32,
    starting_points_seed: u32,
    verbose_parser: bool,
    noxt: bool,
    now: bool,
    restraints: Option<String>,
    use_anm: bool,
    glowworms: u32,
    membrane: bool,
    receptor_pdb: String,
    ligand_pdb: String,
    receptor_restraints: Option<HashMap<String, Vec<String>>>,
    ligand_restraints: Option<HashMap<String, Vec<String>>>,
}

// ─── CLI argument parser ──────────────────────────────────────────────────────

struct Args {
    receptor_pdb: String,
    ligand_pdb: String,
    swarms: u32,
    glowworms: u32,
    seed: u64,
    anm_seed: u64,
    use_anm: bool,
    anm_rec: usize,
    anm_lig: usize,
    noxt: bool,
    noh: bool,
    now: bool,
    swarm_radius: f64,
    restraints: Option<String>,
}

fn print_usage(prog: &str) {
    eprintln!(
        "Usage: {} receptor.pdb ligand.pdb [OPTIONS]

Options:
  -s, --swarms N          Number of swarms (default: 400)
  -g, --glowworms N       Glowworms per swarm (default: 200)
  --seed N                Random seed (default: 324324)
  --anm-seed N            ANM random seed (default: 324324)
  --anm                   Enable ANM (default: false)
  --anm-rec N             ANM modes for receptor (default: 10)
  --anm-lig N             ANM modes for ligand (default: 10)
  --noxt                  Remove OXT atoms
  --noh                   Remove hydrogen atoms
  --now                   Remove water molecules
  --swarm-radius R        Swarm surface radius in Angstroms (default: 10.0)
  --restraints FILE       Restraints file",
        prog
    );
}

fn parse_args() -> Args {
    let raw: Vec<String> = env::args().collect();
    if raw.len() < 3 {
        print_usage(&raw[0]);
        std::process::exit(1);
    }
    let mut args = Args {
        receptor_pdb: raw[1].clone(),
        ligand_pdb: raw[2].clone(),
        swarms: 400,
        glowworms: 200,
        seed: 324_324,
        anm_seed: 324_324,
        use_anm: false,
        anm_rec: 10,
        anm_lig: 10,
        noxt: false,
        noh: false,
        now: false,
        swarm_radius: 10.0,
        restraints: None,
    };
    let mut i = 3;
    while i < raw.len() {
        match raw[i].as_str() {
            "-s" | "--swarms" => { args.swarms = raw[i+1].parse().unwrap(); i += 2; }
            "-g" | "--glowworms" => { args.glowworms = raw[i+1].parse().unwrap(); i += 2; }
            "--seed" => { args.seed = raw[i+1].parse().unwrap(); i += 2; }
            "--anm-seed" => { args.anm_seed = raw[i+1].parse().unwrap(); i += 2; }
            "--anm" => { args.use_anm = true; i += 1; }
            "--anm-rec" => { args.anm_rec = raw[i+1].parse().unwrap(); i += 2; }
            "--anm-lig" => { args.anm_lig = raw[i+1].parse().unwrap(); i += 2; }
            "--noxt" => { args.noxt = true; i += 1; }
            "--noh" => { args.noh = true; i += 1; }
            "--now" => { args.now = true; i += 1; }
            "--swarm-radius" => { args.swarm_radius = raw[i+1].parse().unwrap(); i += 2; }
            "--restraints" => { args.restraints = Some(raw[i+1].clone()); i += 2; }
            _ => { eprintln!("Unknown option: {}", raw[i]); i += 1; }
        }
    }
    args
}

// ─── Geometry helpers ─────────────────────────────────────────────────────────

fn center_of_mass(pdb: &PDB) -> [f64; 3] {
    let mut cx = 0.0_f64;
    let mut cy = 0.0_f64;
    let mut cz = 0.0_f64;
    let mut count = 0usize;
    for atom in pdb.atoms() {
        cx += atom.x();
        cy += atom.y();
        cz += atom.z();
        count += 1;
    }
    if count == 0 { return [0.0, 0.0, 0.0]; }
    [cx / count as f64, cy / count as f64, cz / count as f64]
}

fn bounding_radius(pdb: &PDB, center: &[f64; 3]) -> f64 {
    let mut max_r2 = 0.0_f64;
    for atom in pdb.atoms() {
        let dx = atom.x() - center[0];
        let dy = atom.y() - center[1];
        let dz = atom.z() - center[2];
        let r2 = dx*dx + dy*dy + dz*dz;
        if r2 > max_r2 { max_r2 = r2; }
    }
    max_r2.sqrt()
}

fn translate_to_origin(pdb: &mut PDB, center: &[f64; 3]) {
    for atom in pdb.atoms_mut() {
        let _ = atom.set_x(atom.x() - center[0]);
        let _ = atom.set_y(atom.y() - center[1]);
        let _ = atom.set_z(atom.z() - center[2]);
    }
}

/// Fibonacci sphere: N uniformly distributed points on sphere of given radius.
fn fibonacci_sphere(n: u32, radius: f64) -> Vec<[f64; 3]> {
    let golden = std::f64::consts::PI * (3.0 - 5.0_f64.sqrt());
    let mut points = Vec::with_capacity(n as usize);
    for i in 0..n {
        let y = 1.0 - (i as f64 / (n as f64 - 1.0).max(1.0)) * 2.0;
        let r = (1.0 - y * y).max(0.0).sqrt();
        let theta = golden * i as f64;
        points.push([
            radius * r * theta.cos(),
            radius * y,
            radius * r * theta.sin(),
        ]);
    }
    points
}

/// Generate a uniformly random unit quaternion using Shoemake's method.
fn random_quaternion(rng: &mut StdRng) -> [f64; 4] {
    let u1: f64 = rng.gen();
    let u2: f64 = rng.gen::<f64>() * 2.0 * std::f64::consts::PI;
    let u3: f64 = rng.gen::<f64>() * 2.0 * std::f64::consts::PI;
    let s1 = (1.0 - u1).sqrt();
    let s2 = u1.sqrt();
    [
        s1 * u2.sin(),
        s1 * u2.cos(),
        s2 * u3.sin(),
        s2 * u3.cos(),
    ]
}

// ─── Main ──────────────────────────────────────────────────────────────────────

fn main() {
    let args = parse_args();

    // ── Read structures ──────────────────────────────────────────────────────
    println!("Reading receptor: {}", args.receptor_pdb);
    let (mut receptor, _) = pdbtbx::open(&args.receptor_pdb, pdbtbx::StrictnessLevel::Loose)
        .expect("Failed to open receptor PDB");
    println!("  {} atoms", receptor.atom_count());

    println!("Reading ligand: {}", args.ligand_pdb);
    let (mut ligand, _) = pdbtbx::open(&args.ligand_pdb, pdbtbx::StrictnessLevel::Loose)
        .expect("Failed to open ligand PDB");
    println!("  {} atoms", ligand.atom_count());

    // ── Translate to origin ──────────────────────────────────────────────────
    let rec_center = center_of_mass(&receptor);
    let lig_center = center_of_mass(&ligand);

    println!("Receptor center of mass: [{:.3}, {:.3}, {:.3}]",
             rec_center[0], rec_center[1], rec_center[2]);
    println!("Ligand center of mass:   [{:.3}, {:.3}, {:.3}]",
             lig_center[0], lig_center[1], lig_center[2]);

    translate_to_origin(&mut receptor, &rec_center);
    translate_to_origin(&mut ligand, &lig_center);

    // ── Save lightdock_ PDB files ────────────────────────────────────────────
    let rec_basename = Path::new(&args.receptor_pdb)
        .file_name().unwrap().to_str().unwrap();
    let lig_basename = Path::new(&args.ligand_pdb)
        .file_name().unwrap().to_str().unwrap();

    let rec_out = format!("lightdock_{}", rec_basename);
    let lig_out = format!("lightdock_{}", lig_basename);

    pdbtbx::save(&receptor, &rec_out, pdbtbx::StrictnessLevel::Loose)
        .expect("Failed to save receptor");
    println!("Saved receptor to: {}", rec_out);

    pdbtbx::save(&ligand, &lig_out, pdbtbx::StrictnessLevel::Loose)
        .expect("Failed to save ligand");
    println!("Saved ligand to: {}", lig_out);

    // ── Generate swarm positions ─────────────────────────────────────────────
    let rec_radius = bounding_radius(&receptor, &[0.0, 0.0, 0.0]);
    let sphere_radius = rec_radius + args.swarm_radius;
    println!("Receptor bounding radius: {:.2} Å, swarm sphere radius: {:.2} Å",
             rec_radius, sphere_radius);

    let swarm_centers = fibonacci_sphere(args.swarms, sphere_radius);
    let actual_swarms = swarm_centers.len() as u32;

    // ── Seed RNG ─────────────────────────────────────────────────────────────
    let mut rng: StdRng = SeedableRng::seed_from_u64(args.seed);

    // ── For each swarm: create directory, write initial positions ────────────
    for (swarm_id, center) in swarm_centers.iter().enumerate() {
        let swarm_dir = format!("swarm_{}", swarm_id);
        fs::create_dir_all(&swarm_dir).expect("Cannot create swarm directory");

        let pos_file = format!("initial_positions_{}.dat", swarm_id);
        let mut file = fs::File::create(&pos_file).expect("Cannot create positions file");

        for _ in 0..args.glowworms {
            // Small random displacement from swarm center (±1 Å)
            let dx: f64 = (rng.gen::<f64>() - 0.5) * 2.0;
            let dy: f64 = (rng.gen::<f64>() - 0.5) * 2.0;
            let dz: f64 = (rng.gen::<f64>() - 0.5) * 2.0;
            let tx = center[0] + dx;
            let ty = center[1] + dy;
            let tz = center[2] + dz;

            let q = random_quaternion(&mut rng);
            // w x y z
            write!(file, "{:.6} {:.6} {:.6} {:.6} {:.6} {:.6} {:.6}",
                   tx, ty, tz, q[3], q[0], q[1], q[2]).unwrap();

            // ANM extents (zero-initialized)
            if args.use_anm {
                for _ in 0..args.anm_rec {
                    write!(file, " 0.000000").unwrap();
                }
                for _ in 0..args.anm_lig {
                    write!(file, " 0.000000").unwrap();
                }
            }
            writeln!(file).unwrap();
        }
    }
    println!("Created {} swarms × {} glowworms = {} initial positions",
             actual_swarms, args.glowworms, actual_swarms * args.glowworms);

    // ── Write setup.json ─────────────────────────────────────────────────────
    let setup = SetupFile {
        seed: Some(args.seed),
        anm_seed: args.anm_seed,
        ftdock_file: None,
        noh: args.noh,
        anm_rec: if args.use_anm { args.anm_rec } else { 0 },
        anm_lig: if args.use_anm { args.anm_lig } else { 0 },
        swarms: actual_swarms,
        starting_points_seed: args.seed as u32,
        verbose_parser: false,
        noxt: args.noxt,
        now: args.now,
        restraints: args.restraints.clone(),
        use_anm: args.use_anm,
        glowworms: args.glowworms,
        membrane: false,
        receptor_pdb: rec_basename.to_string(),
        ligand_pdb: lig_basename.to_string(),
        receptor_restraints: None,
        ligand_restraints: None,
    };

    let json = serde_json::to_string_pretty(&setup).expect("JSON serialization failed");
    fs::write("setup.json", &json).expect("Cannot write setup.json");
    println!("Setup written to setup.json");
    println!("LKlight setup OK");
    println!("");
    println!("Next step:");
    println!("  lklight run setup.json initial_positions_0.dat 100 dfire");
}
