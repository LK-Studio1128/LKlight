/// simulator.rs — core LightDock simulation logic as a reusable library module.
/// Both `lightdock-rust` (original binary) and the new unified `lightdock` binary
/// call into these public functions, avoiding code duplication.

use super::constants::{DEFAULT_LIGHTDOCK_PREFIX, DEFAULT_LIG_NM_FILE, DEFAULT_REC_NM_FILE, DEFAULT_SEED};
use super::cpydock::CPYDOCK;
use super::ddna::DDNA;
use super::dfire::DFIRE;
use super::dfire2::DFIRE2;
use super::dna::DNA;
use super::mj3h::MJ3h;
use super::pisa::PISA;
use super::pydock::PYDOCK;
use super::sd::SD;
use super::scoring::{Method, Score};
use super::sipper::SIPPER;
use super::tobi::TOBI;
use super::vdw::VDW;
use super::GSO;
use npyz::NpyFile;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::error::Error;
use std::fs;
use std::fs::File;
use std::io::BufReader;
use std::path::Path;

// ─── Setup file structure ─────────────────────────────────────────────────────

#[derive(Serialize, Deserialize, Debug)]
pub struct SetupFile {
    pub seed: Option<u64>,
    pub anm_seed: u64,
    pub ftdock_file: Option<String>,
    pub noh: bool,
    pub anm_rec: usize,
    pub anm_lig: usize,
    pub swarms: u32,
    pub starting_points_seed: u32,
    pub verbose_parser: bool,
    pub noxt: bool,
    pub now: bool,
    pub restraints: Option<String>,
    pub use_anm: bool,
    pub glowworms: u32,
    pub membrane: bool,
    pub receptor_pdb: String,
    pub ligand_pdb: String,
    pub receptor_restraints: Option<HashMap<String, Vec<String>>>,
    pub ligand_restraints: Option<HashMap<String, Vec<String>>>,
}

pub fn read_setup_from_file<P: AsRef<Path>>(path: P) -> Result<SetupFile, Box<dyn Error>> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    Ok(serde_json::from_reader(reader)?)
}

pub fn parse_input_coordinates(swarm_filename: &str) -> Vec<Vec<f64>> {
    let contents = fs::read_to_string(swarm_filename).expect("Error reading input file");
    let mut positions: Vec<Vec<f64>> = Vec::new();
    for s in contents.lines() {
        let s = s.trim();
        if s.is_empty() { continue; }
        let position: Vec<f64> = s.split_whitespace()
            .filter_map(|tok| tok.parse::<f64>().ok())
            .collect();
        if !position.is_empty() {
            positions.push(position);
        }
    }
    positions
}

pub fn parse_method(method_str: &str) -> Option<Method> {
    match method_str.to_lowercase().as_str() {
        "ddna" => Some(Method::DDNA),
        "dfire" | "fastdfire" => Some(Method::DFIRE),
        "dfire2" => Some(Method::DFIRE2),
        "dna" => Some(Method::DNA),
        "mj3h" => Some(Method::MJ3H),
        "pydock" => Some(Method::PYDOCK),
        "cpydock" => Some(Method::CPYDOCK),
        "sd" => Some(Method::SD),
        "pisa" => Some(Method::PISA),
        "sipper" => Some(Method::SIPPER),
        "tobi" => Some(Method::TOBI),
        "vdw" => Some(Method::VDW),
        _ => None,
    }
}

pub fn parse_swarm_id(path: &Path) -> Option<i32> {
    path.file_name()
        .and_then(|s| s.to_str())
        .and_then(|s| s.strip_prefix("initial_positions_"))
        .and_then(|s| s.strip_suffix(".dat"))
        .and_then(|s| s.parse::<i32>().ok())
}

pub fn simulate(
    simulation_path: &str,
    setup: &SetupFile,
    swarm_filename: &str,
    steps: u32,
    method: Method,
) {
    let seed: u64 = setup.seed.unwrap_or(DEFAULT_SEED);

    println!("Reading starting positions from {:?}", swarm_filename);
    let file_path = Path::new(swarm_filename);
    let swarm_id = parse_swarm_id(file_path).expect("Could not parse swarm from swarm filename");
    println!("Swarm ID {:?}", swarm_id);
    let swarm_directory = format!("swarm_{}", swarm_id);

    if !fs::metadata(&swarm_directory)
        .map(|m| m.is_dir())
        .unwrap_or(false)
    {
        eprintln!("Output directory does not exist for swarm {:?}, creating it", swarm_id);
        fs::create_dir(&swarm_directory).expect("Error creating directory");
    }

    println!("Writing to swarm dir {:?}", swarm_directory);
    let positions = parse_input_coordinates(swarm_filename);

    let rec_basename = format!("{}{}", DEFAULT_LIGHTDOCK_PREFIX, setup.receptor_pdb);
    let receptor_filename = if simulation_path.is_empty() {
        rec_basename.clone()
    } else {
        Path::new(simulation_path).join(&rec_basename).to_string_lossy().into_owned()
    };
    println!("Reading receptor: {}", receptor_filename);
    let (receptor, _) = pdbtbx::open(&receptor_filename, pdbtbx::StrictnessLevel::Medium).unwrap();

    let lig_basename = format!("{}{}", DEFAULT_LIGHTDOCK_PREFIX, setup.ligand_pdb);
    let ligand_filename = if simulation_path.is_empty() {
        lig_basename.clone()
    } else {
        Path::new(simulation_path).join(&lig_basename).to_string_lossy().into_owned()
    };
    println!("Reading ligand: {}", ligand_filename);
    let (ligand, _) = pdbtbx::open(&ligand_filename, pdbtbx::StrictnessLevel::Medium).unwrap();

    // ANM data
    let mut rec_nm: Vec<f64> = Vec::new();
    let mut lig_nm: Vec<f64> = Vec::new();
    if setup.use_anm {
        if setup.anm_rec > 0 {
            let bytes = fs::read(DEFAULT_REC_NM_FILE)
                .unwrap_or_else(|e| panic!("Error reading rec ANM [{:?}]: {:?}", DEFAULT_REC_NM_FILE, e));
            let reader = NpyFile::new(&bytes[..]).unwrap();
            rec_nm = reader.into_vec::<f64>().unwrap();
            // n_atoms is implied by nm.len / (3 * n_modes) — no strict check needed
        }
        if setup.anm_lig > 0 {
            let bytes = fs::read(DEFAULT_LIG_NM_FILE)
                .unwrap_or_else(|e| panic!("Error reading lig ANM [{:?}]: {:?}", DEFAULT_LIG_NM_FILE, e));
            let reader = NpyFile::new(&bytes[..]).unwrap();
            lig_nm = reader.into_vec::<f64>().unwrap();
        }
    }

    // Restraints
    let rec_active  = restraints_or_empty(&setup.receptor_restraints, "active");
    let rec_passive = restraints_or_empty(&setup.receptor_restraints, "passive");
    let lig_active  = restraints_or_empty(&setup.ligand_restraints,   "active");
    let lig_passive = restraints_or_empty(&setup.ligand_restraints,   "passive");

    println!("Loading {:?} scoring function", method);
    let scoring: Box<dyn Score> = match method {
        Method::DFIRE => DFIRE::new(receptor, rec_active, rec_passive, rec_nm, setup.anm_rec,
                                    ligand, lig_active, lig_passive, lig_nm, setup.anm_lig, setup.use_anm),
        Method::DFIRE2 => DFIRE2::new(receptor, rec_active, rec_passive, rec_nm, setup.anm_rec,
                                      ligand, lig_active, lig_passive, lig_nm, setup.anm_lig, setup.use_anm),
        Method::DNA   => DNA::new(receptor, rec_active, rec_passive, rec_nm, setup.anm_rec,
                                  ligand, lig_active, lig_passive, lig_nm, setup.anm_lig, setup.use_anm),
        Method::MJ3H  => MJ3h::new(receptor, rec_active, rec_passive, ligand, lig_active, lig_passive),
        Method::PYDOCK => PYDOCK::new(receptor, rec_active, rec_passive, rec_nm, setup.anm_rec,
                                      ligand, lig_active, lig_passive, lig_nm, setup.anm_lig, setup.use_anm),
        Method::CPYDOCK => CPYDOCK::new(receptor, rec_active, rec_passive, rec_nm, setup.anm_rec,
                                       ligand, lig_active, lig_passive, lig_nm, setup.anm_lig, setup.use_anm),
        Method::SD => SD::new(receptor, rec_active, rec_passive, rec_nm, setup.anm_rec,
                              ligand, lig_active, lig_passive, lig_nm, setup.anm_lig, setup.use_anm),
        Method::DDNA  => DDNA::new(receptor, rec_active, rec_passive, rec_nm, setup.anm_rec,
                                   ligand, lig_active, lig_passive, lig_nm, setup.anm_lig, setup.use_anm),
        Method::PISA  => PISA::new(receptor, rec_active, rec_passive, rec_nm, setup.anm_rec,
                                   ligand, lig_active, lig_passive, lig_nm, setup.anm_lig, setup.use_anm),
        Method::SIPPER => SIPPER::new(receptor, rec_active, rec_passive, ligand, lig_active, lig_passive),
        Method::TOBI   => TOBI::new(receptor, rec_active, rec_passive, ligand, lig_active, lig_passive),
        Method::VDW   => VDW::new(receptor, rec_active, rec_passive, rec_nm, setup.anm_rec,
                                  ligand, lig_active, lig_passive, lig_nm, setup.anm_lig, setup.use_anm),
    };

    println!("Creating GSO with {} glowworms", positions.len());
    let mut gso = GSO::new(
        &positions, seed, &scoring,
        setup.use_anm, setup.anm_rec, setup.anm_lig,
        swarm_directory,
    );
    println!("Starting optimization ({} steps)", steps);
    gso.run(steps);
}

/// Parse a LightDock restraints file.
/// Format per line:
///   R <chain>.<resname>.<resnum>     → receptor active
///   R <chain>.<resname>.<resnum> P   → receptor passive
///   L <chain>.<resname>.<resnum>     → ligand active
///   L <chain>.<resname>.<resnum> P   → ligand passive
/// Returns (receptor_restraints, ligand_restraints) each being
/// HashMap with keys "active" and "passive".
pub fn parse_restraints_file(path: &str) -> (HashMap<String, Vec<String>>, HashMap<String, Vec<String>>) {
    let content = match fs::read_to_string(path) {
        Ok(c) => c,
        Err(e) => { eprintln!("Warning: cannot read restraints file {}: {}", path, e); return empty_restraints(); }
    };
    let mut rec: HashMap<String, Vec<String>> = [
        ("active".into(), vec![]),
        ("passive".into(), vec![]),
        ("blocked".into(), vec![]),
    ].into();
    let mut lig: HashMap<String, Vec<String>> = [
        ("active".into(), vec![]),
        ("passive".into(), vec![]),
        ("blocked".into(), vec![]),
    ].into();
    for line in content.lines() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') { continue; }
        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() < 2 { continue; }
        let mol    = parts[0];
        let res_id = parts[1];
        let tag    = parts.get(2).map(|&s| s).unwrap_or("");
        let key = match tag { "P" | "p" => "passive", "B" | "b" => "blocked", _ => "active" };
        match mol {
            "R" | "r" => rec.get_mut(key).unwrap().push(res_id.to_string()),
            "L" | "l" => lig.get_mut(key).unwrap().push(res_id.to_string()),
            _ => {}
        }
    }
    (rec, lig)
}

/// Standalone scoring: open two PDB files (rec + lig), build the requested
/// scoring function, and evaluate at identity pose (tx=ty=tz=0, q=identity).
/// Both PDB files should have their molecules at the desired relative positions.
pub fn score_pdb(
    rec_path: &str,
    lig_path: &str,
    tx: f64, ty: f64, tz: f64,
    qw: f64, qx: f64, qy: f64, qz: f64,
    method: Method,
) -> f64 {

    // Open with Medium strictness; padded PDBs from setup will be fine.
    let (receptor, _) = pdbtbx::open(rec_path, pdbtbx::StrictnessLevel::Loose)
        .unwrap_or_else(|_| panic!("score_pdb: cannot open {}", rec_path));
    let (ligand, _)   = pdbtbx::open(lig_path, pdbtbx::StrictnessLevel::Loose)
        .unwrap_or_else(|_| panic!("score_pdb: cannot open {}", lig_path));

    let none_vec: Vec<String> = vec![];
    let none_f:   Vec<f64>    = vec![];
    let scoring: Box<dyn Score> = match method {
        Method::DFIRE  => DFIRE::new(receptor,  none_vec.clone(), none_vec.clone(), none_f.clone(), 0,
                                     ligand,    none_vec.clone(), none_vec.clone(), none_f.clone(), 0, false),
        Method::DFIRE2 => DFIRE2::new(receptor, none_vec.clone(), none_vec.clone(), none_f.clone(), 0,
                                      ligand,   none_vec.clone(), none_vec.clone(), none_f.clone(), 0, false),
        Method::DNA    => DNA::new(receptor,    none_vec.clone(), none_vec.clone(), none_f.clone(), 0,
                                   ligand,      none_vec.clone(), none_vec.clone(), none_f.clone(), 0, false),
        Method::MJ3H   => MJ3h::new(receptor,   none_vec.clone(), none_vec.clone(),
                                    ligand,      none_vec.clone(), none_vec.clone()),
        Method::PYDOCK => PYDOCK::new(receptor, none_vec.clone(), none_vec.clone(), none_f.clone(), 0,
                                      ligand,   none_vec.clone(), none_vec.clone(), none_f.clone(), 0, false),
        Method::CPYDOCK => CPYDOCK::new(receptor, none_vec.clone(), none_vec.clone(), none_f.clone(), 0,
                                        ligand,  none_vec.clone(), none_vec.clone(), none_f.clone(), 0, false),
        Method::SD => SD::new(receptor, none_vec.clone(), none_vec.clone(), none_f.clone(), 0,
                              ligand,   none_vec.clone(), none_vec.clone(), none_f.clone(), 0, false),
        Method::DDNA   => DDNA::new(receptor,   none_vec.clone(), none_vec.clone(), none_f.clone(), 0,
                                    ligand,     none_vec.clone(), none_vec.clone(), none_f.clone(), 0, false),
        Method::PISA   => PISA::new(receptor,   none_vec.clone(), none_vec.clone(), none_f.clone(), 0,
                                    ligand,     none_vec.clone(), none_vec.clone(), none_f.clone(), 0, false),
        Method::SIPPER => SIPPER::new(receptor, none_vec.clone(), none_vec.clone(),
                                      ligand,   none_vec.clone(), none_vec.clone()),
        Method::TOBI   => TOBI::new(receptor,   none_vec.clone(), none_vec.clone(),
                                    ligand,     none_vec.clone(), none_vec.clone()),
        Method::VDW    => VDW::new(receptor,    none_vec.clone(), none_vec.clone(), none_f.clone(), 0,
                                   ligand,      none_vec.clone(), none_vec.clone(), none_f.clone(), 0, false),
    };
    let rot = super::qt::Quaternion::new(qw, qx, qy, qz);
    scoring.energy(&[tx, ty, tz], &rot, &[], &[])
}

fn empty_restraints() -> (HashMap<String, Vec<String>>, HashMap<String, Vec<String>>) {
    let e = || [("active".into(), vec![]), ("passive".into(), vec![]), ("blocked".into(), vec![])].into();
    (e(), e())
}

fn restraints_or_empty(
    opt: &Option<HashMap<String, Vec<String>>>,
    key: &str,
) -> Vec<String> {
    opt.as_ref().map(|m| m.get(key).cloned().unwrap_or_default()).unwrap_or_default()
}
