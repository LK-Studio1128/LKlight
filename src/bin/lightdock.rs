/// lklight — unified single-binary entry point.
///
/// Subcommands:
///   setup     <receptor.pdb> <ligand.pdb> [OPTIONS]
///   run       <setup.json> <initial_positions.dat> <steps> <method>
///   generate  <lightdock_rec.pdb> <lightdock_lig.pdb> <gso_output> <N>
///   cluster   <gso_output.out> [--cutoff 4.0]
///   rank      <num_swarms> <steps> [--filter-clusters]
///   top       <ranking_file> <N>
///   filter    <ranking_file> <restraints.list> [--receptor-cutoff %] [--ligand-cutoff %]
///   pipeline  <receptor.pdb> <ligand.pdb> <method> [OPTIONS]

extern crate npyz;
extern crate serde;
extern crate serde_json;

use lklight::anm::{build_atom_to_backbone, compute_anm, save_npy};
use lklight::simulator::{parse_method, parse_restraints_file, read_setup_from_file, simulate};
use rayon::prelude::*;
use std::collections::HashMap;
use std::env;
use std::fs;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;
use std::sync::{Arc, Mutex};
use std::thread;

const STACK_SIZE: usize = 8 * 1024 * 1024;

// ─── PDB atom (shared across generate / cluster) ──────────────────────────────

#[derive(Clone, Debug)]
struct Atom {
    serial: u32,
    name: String,
    res_name: String,
    chain_id: char,
    res_seq: i32,
    x: f64,
    y: f64,
    z: f64,
    occupancy: f64,
    temp_factor: f64,
    element: String,
    is_hetatm: bool,
}

fn parse_pdb(filename: &str) -> Vec<Atom> {
    let file = match fs::File::open(filename) {
        Ok(f) => f,
        Err(e) => { eprintln!("Cannot open {}: {}", filename, e); return Vec::new(); }
    };
    let mut atoms = Vec::new();
    for line in BufReader::new(file).lines() {
        let line = line.unwrap_or_default();
        if line.len() < 27 { continue; }
        let record = line[..6].trim();
        let is_hetatm = record == "HETATM";
        if record != "ATOM" && !is_hetatm { continue; }
        atoms.push(Atom {
            serial:     line[6..11].trim().parse().unwrap_or(0),
            name:       line[12..16].trim().to_string(),
            res_name:   line[17..21].trim().to_string(),
            chain_id:   line.chars().nth(21).unwrap_or('A'),
            res_seq:    line[22..26].trim().parse().unwrap_or(0),
            x:          line[30..38].trim().parse().unwrap_or(0.0),
            y:          line[38..46].trim().parse().unwrap_or(0.0),
            z:          line[46..54].trim().parse().unwrap_or(0.0),
            occupancy:  if line.len()>60 { line[54..60].trim().parse().unwrap_or(1.0) } else { 1.0 },
            temp_factor:if line.len()>66 { line[60..66].trim().parse().unwrap_or(0.0) } else { 0.0 },
            element:    if line.len()>78 { line[76..78].trim().to_string() } else { String::new() },
            is_hetatm,
        });
    }
    atoms
}

fn write_pdb(atoms: &[Atom], filename: &str) {
    let file = fs::File::create(filename).expect(&format!("Cannot create {}", filename));
    let mut writer = BufWriter::new(file);
    for a in atoms {
        writeln!(writer,
            "{}{:5} {:<4} {:<4}{}{:4}    {:8.3}{:8.3}{:8.3}{:6.2}{:6.2}          {:>2}",
            if a.is_hetatm { "HETATM" } else { "ATOM  " },
            a.serial, a.name, a.res_name, a.chain_id, a.res_seq,
            a.x, a.y, a.z, a.occupancy, a.temp_factor, a.element,
        ).unwrap();
    }
    writeln!(writer, "END").unwrap();
}

// ─── Quaternion (for generate conformation) ───────────────────────────────────

#[derive(Clone, Copy, Debug)]
struct Quaternion { w: f64, x: f64, y: f64, z: f64 }

impl Quaternion {
    fn new(w: f64, x: f64, y: f64, z: f64) -> Self { Quaternion { w, x, y, z } }
    fn rotate(&self, v: [f64; 3]) -> [f64; 3] {
        let (qw, qx, qy, qz) = (self.w, self.x, self.y, self.z);
        let (vx, vy, vz) = (v[0], v[1], v[2]);
        let tx = 2.0 * (qy * vz - qz * vy);
        let ty = 2.0 * (qz * vx - qx * vz);
        let tz = 2.0 * (qx * vy - qy * vx);
        [vx + qw * tx + qy * tz - qz * ty,
         vy + qw * ty + qz * tx - qx * tz,
         vz + qw * tz + qx * ty - qy * tx]
    }
}

// ─── GSO output entry ─────────────────────────────────────────────────────────

#[derive(Debug, Clone)]
struct GSOEntry {
    translation: [f64; 3],
    rotation: Quaternion,
    luciferin: f64,
    scoring: f64,
    glowworm_id: usize,
    all_extents: Vec<f64>,   // ANM extents: rec_ext_0..N, lig_ext_0..M
}

fn parse_gso_file(filename: &str) -> Vec<GSOEntry> {
    let file = match fs::File::open(filename) {
        Ok(f) => f,
        Err(_) => return Vec::new(),
    };
    let mut entries = Vec::new();
    let mut id = 0usize;
    for line in BufReader::new(file).lines() {
        let line = line.unwrap_or_default();
        let s = line.trim();
        if s.starts_with('#') || s.is_empty() { continue; }
        if let Some(pend) = s.find(')') {
            let vals: Vec<f64> = s[1..pend].split(',')
                .filter_map(|t| t.trim().parse::<f64>().ok()).collect();
            if vals.len() < 7 { id += 1; continue; }
            let rest: Vec<&str> = s[pend+1..].split_whitespace().collect();
            let luciferin = rest.get(2).and_then(|v| v.parse().ok()).unwrap_or(0.0);
            let scoring   = rest.get(5).and_then(|v| v.parse().ok()).unwrap_or(0.0);
            let all_extents = vals[7..].to_vec();  // ANM extents (may be empty)
            entries.push(GSOEntry {
                translation: [vals[0], vals[1], vals[2]],
                rotation: Quaternion::new(vals[3], vals[4], vals[5], vals[6]),
                luciferin, scoring,
                glowworm_id: id,
                all_extents,
            });
        }
        id += 1;
    }
    entries
}

// ═══════════════════════════════════════════════════════════════════════════════
// Subcommand: SETUP
// ═══════════════════════════════════════════════════════════════════════════════

/// Open a PDB file tolerantly: pad every ATOM/HETATM line to ≥80 chars before
/// parsing so that pdbtbx doesn't raise InvalidatingError for missing element
/// symbol / charge columns (cols 77-80), which are optional in older PDB files.
///
/// pdbtbx 0.11 additionally fails on free-text metadata records such as
/// "REMARK DATE:23-Dec-2018" (it tries to parse them as usize). On the first
/// parse failure we strip metadata records (REMARK/USER/HEADER/TITLE/...) and
/// retry; atom coordinates are never modified.
fn open_pdb_padded(path: &str) -> pdbtbx::PDB {
    use std::io::Cursor;
    let raw = fs::read_to_string(path)
        .unwrap_or_else(|e| panic!("Cannot read {}: {}", path, e));

    fn parse(raw: &str, path: &str) -> Option<pdbtbx::PDB> {
        let padded: String = raw.lines().map(|line| {
            let rec = &line[..line.len().min(6)];
            if rec.starts_with("ATOM") || rec.starts_with("HETATM") {
                if line.len() < 80 {
                    format!("{:<80}\n", line)   // pad with spaces on right
                } else {
                    format!("{}\n", line)
                }
            } else {
                format!("{}\n", line)
            }
        }).collect();
        let cursor = Cursor::new(padded.as_bytes().to_vec());
        let reader = std::io::BufReader::new(cursor);
        pdbtbx::open_pdb_raw(reader, pdbtbx::Context::show(path),
                             pdbtbx::StrictnessLevel::Loose)
            .ok().map(|(pdb, _errs)| pdb)
    }

    if let Some(pdb) = parse(&raw, path) {
        return pdb;
    }
    // pdbtbx 0.11 cannot parse free-text metadata lines; strip them and retry
    let stripped: String = raw.lines().filter(|line| {
        let rec = &line[..line.len().min(6)];
        !(rec.starts_with("REMARK") || rec.starts_with("USER")
          || rec.starts_with("HEADER") || rec.starts_with("TITLE")
          || rec.starts_with("COMPND") || rec.starts_with("SOURCE")
          || rec.starts_with("KEYWDS") || rec.starts_with("EXPDTA")
          || rec.starts_with("AUTHOR") || rec.starts_with("REVDAT")
          || rec.starts_with("JRNL") || rec.starts_with("FORMUL")
          || rec.starts_with("HET"))
    }).map(|l| format!("{}\n", l)).collect();
    match parse(&stripped, path) {
        Some(pdb) => pdb,
        None => panic!("Failed to parse PDB: {}", path),
    }
}

/// True if `atom` is a hydrogen (or deuterium).
///
/// Primary check uses the parsed element column; when that column is missing or
/// unreliable (common in legacy / tool-generated PDBs), we fall back to the atom
/// name. In PDB v3 the hydrogen name may carry a leading digit (e.g. `1H5'`,
/// `2HB`) or be a hydroxyl hydrogen (`HO5'`, `HO3'`, `HO'2`), so we inspect the
/// first *alphabetic* character. This is safe for protein/nucleic structures,
/// which contain no elements whose symbol starts with H other than hydrogen.
fn is_hydrogen_atom(atom: &pdbtbx::Atom) -> bool {
    if atom.element() == Some(&pdbtbx::Element::H) {
        return true;
    }
    let name = atom.name().trim();
    matches!(
        name.chars().find(|c| c.is_ascii_alphabetic()),
        Some('H') | Some('h') | Some('D') | Some('d')
    )
}

/// True if `residue` is a water molecule (covers common PDB / force-field names).
fn is_water_residue(res_name: &str) -> bool {
    matches!(
        res_name.trim().to_ascii_uppercase().as_str(),
        "HOH" | "WAT" | "H2O" | "SOL" | "TIP" | "TIP3" | "TIP4" | "TIP5" | "T3P" | "T4P" | "DOD"
    )
}

/// Apply `setup`-stage atom filtering in place, mirroring LightDock semantics:
/// `--noh` strips hydrogens, `--noxt` strips terminal `OXT`, `--now` strips water.
///
/// Filtering here (before the structure is written to `lightdock_<name>.pdb`)
/// guarantees the `run`/scoring stage never sees these atoms. This is essential
/// for the AMBER-based DNA/dDNA scoring, whose atom-type table does not include
/// non-standard hydrogen names such as `HO5'` / `HO3'` and would otherwise abort.
fn apply_setup_atom_filters(pdb: &mut pdbtbx::PDB, noh: bool, noxt: bool, now: bool) -> usize {
    let before = pdb.atom_count();
    if now {
        pdb.remove_residues_by(|res| is_water_residue(res.name().unwrap_or("")));
    }
    if noh || noxt {
        pdb.remove_atoms_by(|atom| {
            (noh && is_hydrogen_atom(atom))
                || (noxt && atom.name().trim().eq_ignore_ascii_case("OXT"))
        });
    }
    if noh || noxt || now {
        pdb.remove_empty();
    }
    before.saturating_sub(pdb.atom_count())
}

// ─── ANM helper: extract backbone and all-atom data from a pdbtbx PDB ─────────
fn extract_anm_data(pdb: &pdbtbx::PDB, is_protein: bool)
    -> (Vec<[f64;3]>, Vec<[f64;3]>, Vec<(i32, char)>, Vec<(i32, char)>)
{
    let bb_name = if is_protein { "CA" } else { "P" };
    let alt_bb  = "C4'";          // DNA fallback
    let mut bb_coords:   Vec<[f64;3]>     = Vec::new();
    let mut all_coords:  Vec<[f64;3]>     = Vec::new();
    let mut bb_keys:     Vec<(i32, char)> = Vec::new();
    let mut all_keys:    Vec<(i32, char)> = Vec::new();

    for chain in pdb.chains() {
        let cid = chain.id().chars().next().unwrap_or('A');
        for res in chain.residues() {
            let rseq = res.serial_number();
            let key  = (rseq as i32, cid);
            let mut has_bb = false;
            for atom in res.atoms() {
                let nm = atom.name().trim();
                // backbone atom?
                if nm == bb_name || (!has_bb && nm == alt_bb) {
                    bb_coords.push([atom.x(), atom.y(), atom.z()]);
                    bb_keys.push(key);
                    has_bb = true;
                }
                // all heavy atoms (skip H/D)
                if !nm.starts_with('H') && !nm.starts_with('D') {
                    all_coords.push([atom.x(), atom.y(), atom.z()]);
                    all_keys.push(key);
                }
            }
        }
    }
    (bb_coords, all_coords, bb_keys, all_keys)
}

// ─── Restraints-biased sphere point selection ─────────────────────────────────
/// Given a list of candidate sphere points (already on the sphere surface at radius r),
/// return `n` of them concentrated toward `bias_dir` (unit vector).
/// Weight = 0.2 + 0.8 * max(cosθ, 0), higher points near bias_dir are preferred.
/// Returned vec is deterministic (sorted by weight desc, then by index for stability).
fn biased_sphere_selection(candidates: &[[f64;3]], bias_dir: [f64;3], n: usize) -> Vec<[f64;3]> {
    let norm = (bias_dir[0]*bias_dir[0] + bias_dir[1]*bias_dir[1] + bias_dir[2]*bias_dir[2]).sqrt();
    if norm < 1e-9 || candidates.is_empty() {
        return candidates.iter().take(n).copied().collect();
    }
    let d = [bias_dir[0]/norm, bias_dir[1]/norm, bias_dir[2]/norm];
    let r = {
        let p = candidates[0];
        (p[0]*p[0]+p[1]*p[1]+p[2]*p[2]).sqrt().max(1e-9)
    };
    let mut indexed: Vec<(usize, f64)> = candidates.iter().enumerate().map(|(i, p)| {
        let pr = (p[0]*p[0]+p[1]*p[1]+p[2]*p[2]).sqrt().max(1e-9);
        let cos_th = (p[0]*d[0] + p[1]*d[1] + p[2]*d[2]) / pr;
        let w = 0.2_f64 + 0.8 * cos_th.max(0.0);
        (i, w)
    }).collect();
    indexed.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));
    let _ = r;
    indexed.iter().take(n).map(|(i, _)| candidates[*i]).collect()
}

fn cmd_setup(args: &[String]) {
    if args.len() < 2 {
        eprintln!("Usage: lklight setup <receptor.pdb> <ligand.pdb> [OPTIONS]");
        eprintln!("  -s N   swarms (default 400)   -g N  glowworms (default 200)");
        eprintln!("  --anm  enable ANM  --anm-rec N  --anm-lig N  --anm-rmsd R");
        eprintln!("  --seed N  --swarm-radius R  --restraints FILE");
        return;
    }
    use pdbtbx::PDB;
    use rand::rngs::StdRng;
    use rand::{Rng, SeedableRng};
    use serde::{Deserialize, Serialize};

    #[derive(Serialize, Deserialize)]
    struct SetupJson {
        seed: Option<u64>, anm_seed: u64, ftdock_file: Option<String>,
        noh: bool, anm_rec: usize, anm_lig: usize, swarms: u32,
        starting_points_seed: u32, verbose_parser: bool, noxt: bool, now: bool,
        restraints: Option<String>, use_anm: bool, glowworms: u32, membrane: bool,
        receptor_pdb: String, ligand_pdb: String,
        receptor_restraints: Option<HashMap<String, Vec<String>>>,
        ligand_restraints:   Option<HashMap<String, Vec<String>>>,
    }

    let rec_file = &args[0];
    let lig_file = &args[1];
    let (mut swarms, mut glowworms, mut seed, mut use_anm, mut anm_rec, mut anm_lig,
         mut swarm_radius, mut noxt, mut noh, mut now) =
        (400u32, 200u32, 324_324u64, false, 10usize, 10usize, 3.0_f64, false, false, false);
    let mut restraints_file: Option<String> = None;
    let mut anm_rec_rmsd = 0.5_f64;
    let mut anm_lig_rmsd = 0.5_f64;
    let mut i = 2usize;
    while i < args.len() {
        match args[i].as_str() {
            "-s"|"--swarms"       => { swarms = args[i+1].parse().unwrap(); i+=2; }
            "-g"|"--glowworms"    => { glowworms = args[i+1].parse().unwrap(); i+=2; }
            "--seed"              => { seed = args[i+1].parse().unwrap(); i+=2; }
            "--anm"               => { use_anm = true; i+=1; }
            "--anm-rec"           => { anm_rec = args[i+1].parse().unwrap(); i+=2; }
            "--anm-lig"           => { anm_lig = args[i+1].parse().unwrap(); i+=2; }
            "--anm-rmsd"          => {
                let v: f64 = args[i+1].parse().unwrap();
                anm_rec_rmsd = v; anm_lig_rmsd = v; i+=2;
            }
            "--anm-rec-rmsd"      => { anm_rec_rmsd = args[i+1].parse().unwrap(); i+=2; }
            "--anm-lig-rmsd"      => { anm_lig_rmsd = args[i+1].parse().unwrap(); i+=2; }
            "--swarm-radius"      => { swarm_radius = args[i+1].parse().unwrap(); i+=2; }
            "--noxt"              => { noxt = true; i+=1; }
            "--noh"               => { noh = true; i+=1; }
            "--now"               => { now = true; i+=1; }
            "--restraints"        => { restraints_file = Some(args[i+1].clone()); i+=2; }
            _                     => { i+=1; }
        }
    }

    fn com(pdb: &PDB) -> [f64;3] {
        let (mut cx,mut cy,mut cz,mut n) = (0.,0.,0.,0usize);
        for a in pdb.atoms() { cx+=a.x(); cy+=a.y(); cz+=a.z(); n+=1; }
        if n==0 { return [0.,0.,0.]; }
        [cx/n as f64, cy/n as f64, cz/n as f64]
    }
    fn bound_r(pdb: &PDB) -> f64 {
        let mut mx = 0.0f64;
        for a in pdb.atoms() { let r=a.x()*a.x()+a.y()*a.y()+a.z()*a.z(); if r>mx{mx=r;} }
        mx.sqrt()
    }
    /// Mean atom distance from the (origin-centered) CoM. LightDock's reference
    /// points sit near the molecular surface, which is far closer to the mean
    /// radius than to the max radius; using `bound_r` puts initial poses ~max-r
    /// away from the receptor, from which blind docking cannot converge.
    fn avg_r(pdb: &PDB) -> f64 {
        let mut s = 0.0f64; let mut n = 0usize;
        for a in pdb.atoms() { s += (a.x()*a.x()+a.y()*a.y()+a.z()*a.z()).sqrt(); n += 1; }
        if n == 0 { 10.0 } else { s / n as f64 }
    }
    fn translate(pdb: &mut PDB, c: &[f64;3]) {
        for a in pdb.atoms_mut() {
            let _ = a.set_x(a.x()-c[0]);
            let _ = a.set_y(a.y()-c[1]);
            let _ = a.set_z(a.z()-c[2]);
        }
    }
    fn fib_sphere(n: u32, r: f64) -> Vec<[f64;3]> {
        let g = std::f64::consts::PI * (3.0 - 5.0_f64.sqrt());
        (0..n).map(|i| {
            let y = 1.0 - (i as f64 / (n as f64 - 1.0).max(1.0)) * 2.0;
            let rr = (1.0-y*y).max(0.).sqrt();
            let th = g * i as f64;
            [r*rr*th.cos(), r*y, r*rr*th.sin()]
        }).collect()
    }
    fn rand_quat(rng: &mut StdRng) -> [f64;4] {
        let u1: f64 = rng.gen();
        let u2: f64 = rng.gen::<f64>() * 2.0 * std::f64::consts::PI;
        let u3: f64 = rng.gen::<f64>() * 2.0 * std::f64::consts::PI;
        let s1=(1.-u1).sqrt(); let s2=u1.sqrt();
        [s1*u2.sin(), s1*u2.cos(), s2*u3.sin(), s2*u3.cos()]
    }

    println!("Reading receptor: {}", rec_file);
    let mut rec = open_pdb_padded(rec_file);
    println!("Reading ligand:   {}", lig_file);
    let mut lig = open_pdb_padded(lig_file);

    // ── Atom filtering (--noh / --noxt / --now) ──────────────────────────────
    // Must run before CoM/translation/ANM/save so the cleaned structure is what
    // every downstream stage (including scoring) consumes.
    if noh || noxt || now {
        let rec_removed = apply_setup_atom_filters(&mut rec, noh, noxt, now);
        let lig_removed = apply_setup_atom_filters(&mut lig, noh, noxt, now);
        println!(
            "Atom filters [noh={} noxt={} now={}]: removed {} receptor / {} ligand atoms",
            noh, noxt, now, rec_removed, lig_removed
        );
    }

    let rc = com(&rec); let lc = com(&lig);
    println!("Receptor CoM: [{:.3},{:.3},{:.3}]  Ligand CoM: [{:.3},{:.3},{:.3}]",
             rc[0],rc[1],rc[2], lc[0],lc[1],lc[2]);
    translate(&mut rec, &rc); translate(&mut lig, &lc);

    let rec_base = Path::new(rec_file).file_name().unwrap().to_str().unwrap();
    let lig_base = Path::new(lig_file).file_name().unwrap().to_str().unwrap();
    let rec_out = format!("lightdock_{}", rec_base);
    let lig_out = format!("lightdock_{}", lig_base);
    pdbtbx::save(&rec, &rec_out, pdbtbx::StrictnessLevel::Loose).expect("Save receptor failed");
    pdbtbx::save(&lig, &lig_out, pdbtbx::StrictnessLevel::Loose).expect("Save ligand failed");
    println!("Saved {} and {}", rec_out, lig_out);

    // ── Parse restraints file (needed for ANM and swarm bias) ────────────────
    let (rec_rst, lig_rst): (Option<HashMap<String,Vec<String>>>, Option<HashMap<String,Vec<String>>>) =
        if let Some(ref rf) = restraints_file {
            let (r, l) = parse_restraints_file(rf);
            let rec_active_n = r.get("active").map(|v| v.len()).unwrap_or(0);
            let lig_active_n = l.get("active").map(|v| v.len()).unwrap_or(0);
            println!("Restraints: {} receptor active, {} ligand active", rec_active_n, lig_active_n);
            (Some(r), Some(l))
        } else { (None, None) };

    // ── ANM computation ──────────────────────────────────────────────────────
    if use_anm {
        let anm_rec_n = if anm_rec > 0 { anm_rec } else { 10 };
        let anm_lig_n = if anm_lig > 0 { anm_lig } else { 10 };

        println!("Computing ANM for receptor ({} modes) …", anm_rec_n);
        let (bb_rec, all_rec, bb_keys_rec, all_keys_rec) = extract_anm_data(&rec, true);
        let a2b_rec = build_atom_to_backbone(&bb_keys_rec, &all_keys_rec);
        let modes_rec = compute_anm(&bb_rec, &a2b_rec, all_rec.len(), anm_rec_n, anm_rec_rmsd);
        save_npy("lightdock_rec.nm.npy", &modes_rec.data)
            .expect("Cannot save rec_nm.npy");
        println!("Saved lightdock_rec.nm.npy  ({} modes × {} atoms × 3, rmsd={:.3})",
            modes_rec.n_modes, modes_rec.n_atoms, anm_rec_rmsd);

        println!("Computing ANM for ligand ({} modes) …", anm_lig_n);
        let (bb_lig, all_lig, bb_keys_lig, all_keys_lig) = extract_anm_data(&lig, true);
        let a2b_lig = build_atom_to_backbone(&bb_keys_lig, &all_keys_lig);
        let modes_lig = compute_anm(&bb_lig, &a2b_lig, all_lig.len(), anm_lig_n, anm_lig_rmsd);
        save_npy("lightdock_lig.nm.npy", &modes_lig.data)
            .expect("Cannot save lig_nm.npy");
        println!("Saved lightdock_lig.nm.npy  ({} modes × {} atoms × 3)",
            modes_lig.n_modes, modes_lig.n_atoms);
    }

    // ── Swarm positions ──────────────────────────────────────────────────────
    // Initial poses sit near the molecular surface: mean radius + small offset.
    // (Formerly `bound_r` (max radius) + 10 Å, which placed the ligand too far
    // from the receptor for blind docking to converge.)
    let sphere_r = avg_r(&rec) + swarm_radius;

    // Build candidate pool: 5× more points than needed for biased selection
    let n_candidates = (swarms * 5).max(2000);
    let all_candidates = fib_sphere(n_candidates, sphere_r);

    // Compute bias direction from receptor restraint residues (if any)
    let bias_dir: Option<[f64;3]> = rec_rst.as_ref().and_then(|rr| {
        let active = rr.get("active")?;
        if active.is_empty() { return None; }
        // Find Cα coordinates of restraint residues in origin-centered receptor
        let mut cx = 0.0_f64; let mut cy = 0.0; let mut cz = 0.0; let mut cnt = 0usize;
        for chain in rec.chains() {
            let cid = chain.id().chars().next().unwrap_or('A');
            for res in chain.residues() {
                let res_id = format!("{}.{}.{}",
                    cid,
                    res.name().unwrap_or(""),
                    res.serial_number());
                if active.contains(&res_id) {
                    for atom in res.atoms() {
                        if atom.name().trim() == "CA" {
                            cx += atom.x(); cy += atom.y(); cz += atom.z(); cnt += 1;
                        }
                    }
                }
            }
        }
        if cnt == 0 { return None; }
        Some([cx/cnt as f64, cy/cnt as f64, cz/cnt as f64])
    });

    let centers: Vec<[f64;3]> = if let Some(dir) = bias_dir {
        println!("Restraints-biased swarm placement (direction [{:.2},{:.2},{:.2}])",
            dir[0], dir[1], dir[2]);
        biased_sphere_selection(&all_candidates, dir, swarms as usize)
    } else {
        all_candidates.into_iter().take(swarms as usize).collect()
    };

    let mut rng: StdRng = SeedableRng::seed_from_u64(seed);

    for (sid, center) in centers.iter().enumerate() {
        fs::create_dir_all(format!("swarm_{}", sid)).unwrap();
        let mut f = fs::File::create(format!("initial_positions_{}.dat", sid)).unwrap();
        for _ in 0..glowworms {
            let dx: f64 = (rng.gen::<f64>()-0.5)*2.;
            let dy: f64 = (rng.gen::<f64>()-0.5)*2.;
            let dz: f64 = (rng.gen::<f64>()-0.5)*2.;
            let q = rand_quat(&mut rng);
            write!(f, "{:.6} {:.6} {:.6} {:.6} {:.6} {:.6} {:.6}",
                   center[0]+dx, center[1]+dy, center[2]+dz, q[3],q[0],q[1],q[2]).unwrap();
            if use_anm {
                for _ in 0..(anm_rec+anm_lig) { write!(f, " 0.000000").unwrap(); }
            }
            writeln!(f).unwrap();
        }
    }

    let setup = SetupJson {
        seed: Some(seed), anm_seed: seed, ftdock_file: None,
        noh, anm_rec: if use_anm{anm_rec}else{0}, anm_lig: if use_anm{anm_lig}else{0},
        swarms: swarms as u32, starting_points_seed: seed as u32,
        verbose_parser: false, noxt, now,
        restraints: restraints_file,
        use_anm, glowworms, membrane: false,
        receptor_pdb: rec_base.to_string(), ligand_pdb: lig_base.to_string(),
        receptor_restraints: rec_rst, ligand_restraints: lig_rst,
    };
    fs::write("setup.json", serde_json::to_string_pretty(&setup).unwrap()).unwrap();
    println!("Created {} swarms × {} glowworms. Written setup.json.", swarms, glowworms);
}

// ═══════════════════════════════════════════════════════════════════════════════
// Subcommand: RUN  (delegates to library simulate())
// ═══════════════════════════════════════════════════════════════════════════════

fn cmd_run(args: &[String]) {
    if args.len() < 4 {
        eprintln!("Usage: lklight run <setup.json> <initial_positions.dat> <steps> <method>");
        return;
    }
    let steps: u32 = match args[2].parse() {
        Ok(n) => n,
        Err(_) => { eprintln!("steps must be integer"); return; }
    };
    let method = match parse_method(&args[3]) {
        Some(m) => m,
        None => { eprintln!("Unknown method: {}", args[3]); return; }
    };
    let setup = match read_setup_from_file(&args[0]) {
        Ok(s) => s,
        Err(e) => { eprintln!("Cannot read setup: {:?}", e); return; }
    };
    let sim_path = Path::new(&args[0]).parent().unwrap();
    simulate(sim_path.to_str().unwrap(), &setup, &args[1], steps, method);
}

// ═══════════════════════════════════════════════════════════════════════════════
// Subcommand: GENERATE  (apply rotation+translation, write PDB files)
// ═══════════════════════════════════════════════════════════════════════════════

/// Load a .npy file as a flat Vec<f64>.  Returns None if the file is missing or unreadable.
fn load_npy_f64(path: &str) -> Option<Vec<f64>> {
    if !Path::new(path).exists() { return None; }
    let bytes = fs::read(path).ok()?;
    let reader = npyz::NpyFile::new(&bytes[..]).ok()?;
    reader.into_vec::<f64>().ok()
}

/// Apply ANM deformation to a list of atoms in-place.
/// `modes`: flat (n_modes × n_atoms × 3), `extents`: one scalar per mode.
fn apply_anm(atoms: &mut [Atom], modes: &[f64], extents: &[f64]) {
    let n_atoms = atoms.len();
    if n_atoms == 0 || modes.is_empty() || extents.is_empty() { return; }
    for (nm, &ext) in extents.iter().enumerate() {
        if ext == 0.0 { continue; }
        let base = nm * n_atoms * 3;
        if base + n_atoms * 3 > modes.len() { break; }
        for a in 0..n_atoms {
            atoms[a].x += modes[base + a * 3    ] * ext;
            atoms[a].y += modes[base + a * 3 + 1] * ext;
            atoms[a].z += modes[base + a * 3 + 2] * ext;
        }
    }
}

fn cmd_generate(args: &[String]) {
    if args.len() < 4 {
        eprintln!("Usage: lklight generate <rec.pdb> <lig.pdb> <gso_output|.dat> <N>");
        eprintln!("  Reads setup.json and lightdock_rec.nm.npy / lightdock_lig.nm.npy");
        eprintln!("  from CWD if ANM was used during setup.");
        return;
    }
    let rec_atoms_orig = parse_pdb(&args[0]);
    let lig_atoms_orig = parse_pdb(&args[1]);
    let gso_file  = &args[2];
    let n: usize  = args[3].parse().expect("N must be integer");

    // Try loading ANM modes and setup parameters
    let setup     = read_setup_from_file("setup.json").ok();
    let n_anm_rec = setup.as_ref().map(|s| s.anm_rec).unwrap_or(0);
    let n_anm_lig = setup.as_ref().map(|s| s.anm_lig).unwrap_or(0);
    let modes_rec = load_npy_f64("lightdock_rec.nm.npy");
    let modes_lig = load_npy_f64("lightdock_lig.nm.npy");

    let entries = parse_gso_file(gso_file);
    let dest = Path::new(gso_file).parent().unwrap_or(Path::new("."));
    let take_n = n.min(entries.len());

    // Parallel: each conformation is independent (clone + ANM + rotate/translate + write).
    // Safe to parallelize here because cmd_generate is called as a one-shot subcommand
    // OUTSIDE the GSO/swarm.rs rayon chain — no nested rayon.
    entries
        .par_iter()
        .take(take_n)
        .enumerate()
        .for_each(|(i, e)| {
            // ── Receptor: clone + apply ANM
            let mut rec_pose = rec_atoms_orig.clone();
            if let Some(ref modes) = modes_rec {
                if n_anm_rec > 0 && !e.all_extents.is_empty() {
                    let rec_ext = &e.all_extents[..n_anm_rec.min(e.all_extents.len())];
                    apply_anm(&mut rec_pose, modes, rec_ext);
                }
            }

            // ── Ligand: clone + apply ANM + rotate + translate
            let mut lig_pose = lig_atoms_orig.clone();
            if let Some(ref modes) = modes_lig {
                if n_anm_lig > 0 && e.all_extents.len() > n_anm_rec {
                    let lig_ext = &e.all_extents[n_anm_rec..n_anm_rec + n_anm_lig.min(e.all_extents.len().saturating_sub(n_anm_rec))];
                    apply_anm(&mut lig_pose, modes, lig_ext);
                }
            }
            let rec_len = rec_pose.len() as u32;
            for (j, atom) in lig_pose.iter_mut().enumerate() {
                let r = e.rotation.rotate([atom.x, atom.y, atom.z]);
                atom.x = r[0] + e.translation[0];
                atom.y = r[1] + e.translation[1];
                atom.z = r[2] + e.translation[2];
                atom.serial = rec_len + j as u32 + 1;
            }
            let mut all = rec_pose;
            all.extend(lig_pose);
            write_pdb(&all, dest.join(format!("lightdock_{}.pdb", i)).to_str().unwrap());
        });
    println!("Generated {} conformations in {:?}", take_n, dest);
}

// ═══════════════════════════════════════════════════════════════════════════════
// Subcommand: CLUSTER  (BSAS algorithm with Cα/P RMSD) — NEW
// ═══════════════════════════════════════════════════════════════════════════════

fn backbone_coords(atoms: &[Atom]) -> Vec<[f64;3]> {
    atoms.iter()
        .filter(|a| a.name == "CA" || a.name == "P")
        .map(|a| [a.x, a.y, a.z])
        .collect()
}

fn rmsd(a: &[[f64;3]], b: &[[f64;3]]) -> f64 {
    let n = a.len().min(b.len());
    if n == 0 { return f64::MAX; }
    let sum: f64 = (0..n).map(|i| {
        let d0=a[i][0]-b[i][0]; let d1=a[i][1]-b[i][1]; let d2=a[i][2]-b[i][2];
        d0*d0 + d1*d1 + d2*d2
    }).sum();
    (sum / n as f64).sqrt()
}

fn cmd_cluster(args: &[String]) {
    if args.is_empty() {
        eprintln!("Usage: lklight cluster <gso_output.out> [--cutoff 4.0]");
        return;
    }
    let gso_file = &args[0];
    let cutoff: f64 = args.windows(2)
        .find(|w| w[0] == "--cutoff")
        .and_then(|w| w[1].parse().ok())
        .unwrap_or(4.0);

    let mut entries = parse_gso_file(gso_file);
    if entries.is_empty() {
        eprintln!("No entries found in {}", gso_file);
        return;
    }
    entries.sort_by(|a, b| b.scoring.partial_cmp(&a.scoring).unwrap_or(std::cmp::Ordering::Equal));

    let swarm_dir = Path::new(gso_file).parent().unwrap_or(Path::new("."));

    // Load backbone Cα/P coords; fall back to translation point if PDB missing
    let mut backbone: HashMap<usize, Vec<[f64;3]>> = HashMap::new();
    let mut pdb_found = 0usize;
    for e in &entries {
        let pdb_path = swarm_dir.join(format!("lightdock_{}.pdb", e.glowworm_id));
        if pdb_path.exists() {
            let atoms = parse_pdb(pdb_path.to_str().unwrap());
            let bc = backbone_coords(&atoms);
            if !bc.is_empty() { pdb_found += 1; backbone.insert(e.glowworm_id, bc); continue; }
        }
        // Fallback: use translation coordinate as a single pseudo-backbone point
        backbone.insert(e.glowworm_id, vec![e.translation]);
    }
    if pdb_found == 0 {
        println!("Note: no lightdock_N.pdb files found — clustering by translation distance");
    } else {
        println!("Loaded {}/{} PDB files for backbone RMSD", pdb_found, entries.len());
    }

    let mut clusters: Vec<Vec<usize>> = Vec::new();
    let mut rep_indices: Vec<usize> = Vec::new();

    'outer: for (ei, e) in entries.iter().enumerate() {
        let bc = backbone.get(&e.glowworm_id).unwrap();
        for (ci, &ri) in rep_indices.iter().enumerate() {
            let rep_bc = backbone.get(&entries[ri].glowworm_id).unwrap();
            if rmsd(bc, rep_bc) <= cutoff {
                clusters[ci].push(ei);
                continue 'outer;
            }
        }
        clusters.push(vec![ei]);
        rep_indices.push(ei);
    }

    let out_path = swarm_dir.join("cluster_representatives.file");
    let mut out = fs::File::create(&out_path).expect("Cannot create cluster file");
    for (ci, members) in clusters.iter().enumerate() {
        let rep = &entries[rep_indices[ci]];
        writeln!(out, "{}:{}:{:.6}:{}:lightdock_{}.pdb",
            ci, members.len(), rep.scoring, rep.glowworm_id, rep.glowworm_id).unwrap();
    }
    println!("BSAS: {} clusters from {} solutions (cutoff={:.1}Å) → {:?}",
             clusters.len(), entries.len(), cutoff, out_path);
}

// ═══════════════════════════════════════════════════════════════════════════════
// Subcommand: FILTER  (restraints satisfaction filtering)  — NEW
// ═══════════════════════════════════════════════════════════════════════════════

/// Check what fraction of `restraint_ids` (e.g. ["A.TYR.234"]) are satisfied:
/// at least one atom of that residue is within `cutoff` Å of any atom in `other_atoms`.
fn restraint_satisfaction(atoms: &[Atom], restraint_ids: &[String], other_atoms: &[Atom], cutoff: f64) -> f64 {
    if restraint_ids.is_empty() { return 1.0; }
    let cutoff2 = cutoff * cutoff;
    let mut satisfied = 0usize;
    for res_id in restraint_ids {
        let parts: Vec<&str> = res_id.splitn(3, '.').collect();
        if parts.len() < 3 { continue; }
        let chain = parts[0]; let _resname = parts[1];
        let resseq: i32 = parts[2].parse().unwrap_or(-1);
        let res_atoms: Vec<&Atom> = atoms.iter()
            .filter(|a| a.chain_id.to_string() == chain && a.res_seq == resseq)
            .collect();
        'res: for ra in &res_atoms {
            for oa in other_atoms {
                let d2 = (ra.x-oa.x).powi(2) + (ra.y-oa.y).powi(2) + (ra.z-oa.z).powi(2);
                if d2 <= cutoff2 { satisfied += 1; break 'res; }
            }
        }
    }
    satisfied as f64 / restraint_ids.len() as f64
}

fn cmd_filter(args: &[String]) {
    if args.len() < 2 {
        eprintln!("Usage: lklight filter <ranking_file> <restraints.list> [OPTIONS]");
        eprintln!("  --contact-cutoff A    contact distance in Å (default 5.0)");
        eprintln!("  --rec-cutoff P        min receptor satisfaction 0-1 (default 0.4)");
        eprintln!("  --lig-cutoff P        min ligand satisfaction 0-1 (default 0.4)");
        return;
    }
    let ranking_file  = &args[0];
    let rst_file      = &args[1];
    let contact_cutoff: f64 = args.windows(2).find(|w| w[0]=="--contact-cutoff")
        .and_then(|w| w[1].parse().ok()).unwrap_or(5.0);
    let rec_min: f64  = args.windows(2).find(|w| w[0]=="--rec-cutoff")
        .and_then(|w| w[1].parse().ok()).unwrap_or(0.4);
    let lig_min: f64  = args.windows(2).find(|w| w[0]=="--lig-cutoff")
        .and_then(|w| w[1].parse().ok()).unwrap_or(0.4);

    let (rec_rst, lig_rst) = parse_restraints_file(rst_file);
    let rec_active = rec_rst.get("active").cloned().unwrap_or_default();
    let lig_active = lig_rst.get("active").cloned().unwrap_or_default();

    let content = fs::read_to_string(ranking_file).expect("Cannot read ranking file");
    let out_file = "filtered.list";
    let mut out = fs::File::create(out_file).unwrap();
    let mut total = 0usize; let mut passed = 0usize;

    for line in content.lines() {
        if line.starts_with('#') { writeln!(out, "{}", line).unwrap(); continue; }
        let parts: Vec<&str> = line.split_whitespace().collect();
        // Expect PDB path as 7th field (rank swarm glow luciferin scoring pdb_path ...)
        if parts.len() < 6 { continue; }
        let pdb_path = parts[5];
        total += 1;
        if !Path::new(pdb_path).exists() { continue; }
        let atoms = parse_pdb(pdb_path);
        if atoms.is_empty() { continue; }

        // Split into receptor (chain of rec_rst) and ligand (chain of lig_rst) atoms
        let rec_chain = rec_active.first().and_then(|s| s.split('.').next()).unwrap_or("A");
        let lig_chain = lig_active.first().and_then(|s| s.split('.').next()).unwrap_or("B");
        let rec_atoms: Vec<Atom> = atoms.iter().filter(|a| a.chain_id.to_string() == rec_chain).cloned().collect();
        let lig_atoms: Vec<Atom> = atoms.iter().filter(|a| a.chain_id.to_string() == lig_chain).cloned().collect();

        let rec_sat = restraint_satisfaction(&rec_atoms, &rec_active, &lig_atoms, contact_cutoff);
        let lig_sat = restraint_satisfaction(&lig_atoms, &lig_active, &rec_atoms, contact_cutoff);

        if rec_sat >= rec_min && lig_sat >= lig_min {
            writeln!(out, "{}  rec_sat={:.2} lig_sat={:.2}", line, rec_sat, lig_sat).unwrap();
            passed += 1;
        }
    }
    println!("Filter complete: {}/{} solutions passed (rec≥{:.0}% lig≥{:.0}%)",
             passed, total, rec_min*100., lig_min*100.);
    println!("Written to {}", out_file);
}

// ═══════════════════════════════════════════════════════════════════════════════
// Subcommand: RANK
// ═══════════════════════════════════════════════════════════════════════════════

#[derive(Clone)]
struct Solution {
    swarm: u32, glowworm: u32, luciferin: f64, scoring: f64, coords: String,
}

fn collect_solutions(num_swarms: u32, steps: u32, filter_clusters: bool) -> Vec<Solution> {
    let mut all = Vec::new();
    for sid in 0..num_swarms {
        let gso_path = format!("swarm_{}/gso_{}.out", sid, steps);
        let cluster_path = format!("swarm_{}/cluster_representatives.file", sid);
        let cluster_ids: Option<Vec<u32>> = if filter_clusters && Path::new(&cluster_path).exists() {
            let c = fs::read_to_string(&cluster_path).unwrap_or_default();
            Some(c.lines().filter_map(|l| l.split(':').nth(3).and_then(|s| s.parse().ok())).collect())
        } else { None };

        let entries = parse_gso_file(&gso_path);
        for e in entries {
            if let Some(ref ids) = cluster_ids {
                if !ids.contains(&(e.glowworm_id as u32)) { continue; }
            }
            all.push(Solution {
                swarm: sid, glowworm: e.glowworm_id as u32,
                luciferin: e.luciferin, scoring: e.scoring,
                coords: format!("({:.3},{:.3},{:.3})", e.translation[0], e.translation[1], e.translation[2]),
            });
        }
    }
    all
}

fn write_rank_file(solutions: &[Solution], filename: &str, by_scoring: bool) {
    let mut sorted = solutions.to_vec();
    if by_scoring {
        sorted.sort_by(|a,b| b.scoring.partial_cmp(&a.scoring).unwrap_or(std::cmp::Ordering::Equal));
    } else {
        sorted.sort_by(|a,b| b.luciferin.partial_cmp(&a.luciferin).unwrap_or(std::cmp::Ordering::Equal));
    }
    let mut f = fs::File::create(filename).unwrap();
    writeln!(f, "# Rank Swarm Glowworm Luciferin Scoring PDB Coords").unwrap();
    for (rank, s) in sorted.iter().enumerate() {
        writeln!(f, "{:5} {:4} {:5} {:12.6} {:12.6} swarm_{}/lightdock_{}.pdb {}",
            rank+1, s.swarm, s.glowworm, s.luciferin, s.scoring,
            s.swarm, s.glowworm, s.coords).unwrap();
    }
    println!("Wrote {} solutions to {}", sorted.len(), filename);
}

fn cmd_rank(args: &[String]) {
    if args.len() < 2 {
        eprintln!("Usage: lklight rank <num_swarms> <steps> [OPTIONS]");
        eprintln!("  --filter-clusters   only keep cluster representatives");
        eprintln!("  --clashes-cutoff S  drop solutions with scoring < S");
        return;
    }
    let num_swarms: u32 = args[0].parse().expect("num_swarms must be integer");
    let steps: u32      = args[1].parse().expect("steps must be integer");
    let filter = args.iter().any(|a| a == "--filter-clusters");
    let clashes_cutoff: Option<f64> = args.windows(2)
        .find(|w| w[0] == "--clashes-cutoff")
        .and_then(|w| w[1].parse().ok());

    let mut solutions = collect_solutions(num_swarms, steps, filter);
    if let Some(cutoff) = clashes_cutoff {
        solutions.retain(|s| s.scoring >= cutoff);
        println!("After clashes cutoff ({:.3}): {} solutions remaining", cutoff, solutions.len());
    }
    println!("Collected {} solutions", solutions.len());
    // Python lgd_rank writes: ranking.list (default), rank_by_scoring, rank_by_luciferin, rank_by_rmsd
    write_rank_file(&solutions, "ranking.list",          true);
    write_rank_file(&solutions, "rank_by_scoring.list",  true);
    write_rank_file(&solutions, "rank_by_luciferin.list", false);
    // rank_by_rmsd: only meaningful when evaluation.dat is present; write empty header otherwise
    {
        let rmsd_avail = Path::new("evaluation.dat").exists();
        let out = "rank_by_rmsd.list";
        let mut f = fs::File::create(out).unwrap();
        writeln!(f, "# Rank Swarm Glowworm Luciferin Scoring RMSD PDB Coords").unwrap();
        if !rmsd_avail {
            writeln!(f, "# No evaluation.dat found — RMSD data not available").unwrap();
        } else {
            for (rank, s) in solutions.iter().enumerate() {
                writeln!(f, "{:5} {:4} {:5} {:12.6} {:12.6}   N/A swarm_{}/lightdock_{}.pdb {}",
                    rank+1, s.swarm, s.glowworm, s.luciferin, s.scoring,
                    s.swarm, s.glowworm, s.coords).unwrap();
            }
        }
        println!("Wrote {} to {}", solutions.len(), out);
    }
}

/// lgd_rank_swarm: write per-swarm rank_by_scoring.list inside each swarm_N folder
fn cmd_rank_swarm(args: &[String]) {
    if args.len() < 2 {
        eprintln!("Usage: lklight rank_swarm <num_swarms> <steps>");
        return;
    }
    let num_swarms: u32 = args[0].parse().expect("num_swarms must be integer");
    let steps: u32      = args[1].parse().expect("steps must be integer");
    let mut total = 0u32;
    for sid in 0..num_swarms {
        let gso = format!("swarm_{}/gso_{}.out", sid, steps);
        if !Path::new(&gso).exists() { continue; }
        let entries = parse_gso_file(&gso);
        if entries.is_empty() { continue; }
        let solutions: Vec<Solution> = entries.iter().map(|e| Solution {
            swarm: sid, glowworm: e.glowworm_id as u32,
            luciferin: e.luciferin, scoring: e.scoring,
            coords: format!("({:.3},{:.3},{:.3})",
                e.translation[0], e.translation[1], e.translation[2]),
        }).collect();
        let out = format!("swarm_{}/rank_by_scoring.list", sid);
        write_rank_file(&solutions, &out, true);
        total += 1;
    }
    println!("Per-swarm ranking written for {}/{} swarms", total, num_swarms);
}

// ═══════════════════════════════════════════════════════════════════════════════
// Subcommand: TOP  (select top-N from ranking file)  — NEW
// ═══════════════════════════════════════════════════════════════════════════════

fn cmd_top(args: &[String]) {
    if args.len() < 2 {
        eprintln!("Usage: lklight top <ranking_file> <N>");
        eprintln!("  Copies swarm_X/lightdock_Y.pdb → top_1.pdb … top_N.pdb");
        eprintln!("  Also writes top_N.list with the header+lines.");
        return;
    }
    let n: usize = args[1].parse().expect("N must be integer");
    let base = Path::new(&args[0]).parent().unwrap_or(Path::new("."));
    let content = fs::read_to_string(&args[0]).expect("Cannot read ranking file");
    let list_file = format!("top_{}.list", n);
    let mut list_out = fs::File::create(&list_file).unwrap();
    let mut count = 0usize;
    for line in content.lines() {
        if line.starts_with('#') { writeln!(list_out, "{}", line).unwrap(); continue; }
        if count >= n { break; }
        writeln!(list_out, "{}", line).unwrap();
        // Field index 5 is the PDB path: swarm_X/lightdock_Y.pdb
        let fields: Vec<&str> = line.split_whitespace().collect();
        if fields.len() >= 6 {
            let pdb_src = base.join(fields[5]);
            let pdb_dst = format!("top_{}.pdb", count + 1);
            if pdb_src.exists() {
                let _ = fs::copy(&pdb_src, &pdb_dst);
            } else {
                eprintln!("Warning: {} not found, skipping PDB copy", pdb_src.display());
            }
        }
        count += 1;
    }
    println!("Top {} → {}, top_1.pdb … top_{}.pdb", count, list_file, count);
}

// ═══════════════════════════════════════════════════════════════════════════════
// Subcommand: PIPELINE  (end-to-end automation)  — NEW
// ═══════════════════════════════════════════════════════════════════════════════

fn cmd_pipeline(args: &[String]) {
    if args.len() < 3 {
        eprintln!("Usage: lklight pipeline <receptor.pdb> <ligand.pdb> <method> [OPTIONS]");
        eprintln!("  -s N           swarms (default 400)");
        eprintln!("  -g N           glowworms (default 200)");
        eprintln!("  --steps N      simulation steps (default 100)");
        eprintln!("  --top N        top-N results (default 10)");
        eprintln!("  --cluster      cluster each swarm after generate");
        eprintln!("  --threads N    parallel swarm simulations (default 1)");
        eprintln!("  --restraints F restraints file");
        return;
    }
    let rec    = args[0].clone();
    let lig    = args[1].clone();
    let method = args[2].clone();
    let mut swarms  = 400u32;
    let mut glowworms = 200u32;
    let mut steps   = 100u32;
    let mut top_n   = 10usize;
    let mut do_cluster = false;
    let mut threads = 1usize;
    let mut restraints_file: Option<String> = None;
    let mut i = 3usize;
    while i < args.len() {
        match args[i].as_str() {
            "-s"|"--swarms"    => { swarms    = args[i+1].parse().unwrap(); i+=2; }
            "-g"|"--glowworms" => { glowworms = args[i+1].parse().unwrap(); i+=2; }
            "--steps"          => { steps     = args[i+1].parse().unwrap(); i+=2; }
            "--top"            => { top_n     = args[i+1].parse().unwrap(); i+=2; }
            "--cluster"        => { do_cluster = true; i+=1; }
            "--threads"        => { threads   = args[i+1].parse().unwrap_or(1); i+=2; }
            "--restraints"     => { restraints_file = Some(args[i+1].clone()); i+=2; }
            _                  => { i+=1; }
        }
    }

    println!("=== LKlight Pipeline ===");
    println!("Receptor={} Ligand={} Method={} Swarms={} Glowworms={} Steps={} Threads={}",
             rec, lig, method, swarms, glowworms, steps, threads);

    // Step 1: Setup
    println!("\n--- Step 1: Setup ---");
    let mut setup_args: Vec<String> = vec![
        rec.clone(), lig.clone(),
        "-s".into(), swarms.to_string(),
        "-g".into(), glowworms.to_string(),
    ];
    if let Some(ref rf) = restraints_file {
        setup_args.push("--restraints".into()); setup_args.push(rf.clone());
    }
    cmd_setup(&setup_args);

    // Step 2: Simulate each swarm (optionally parallel)
    println!("\n--- Step 2: Simulate {} swarms (threads={}) ---", swarms, threads);
    let setup = match read_setup_from_file("setup.json") {
        Ok(s) => s,
        Err(e) => { eprintln!("Cannot read setup.json: {:?}", e); return; }
    };
    let parsed_method = match parse_method(&method) {
        Some(m) => m,
        None => { eprintln!("Unknown method: {}", method); return; }
    };

    let setup_arc = Arc::new(setup);
    let slots = Arc::new(Mutex::new(threads));
    let mut handles = Vec::new();

    for sid in 0..swarms {
        let dat = format!("initial_positions_{}.dat", sid);
        if !Path::new(&dat).exists() {
            eprintln!("Warning: {} not found, skipping", dat);
            continue;
        }
        // Acquire semaphore slot
        loop {
            let mut s = slots.lock().unwrap();
            if *s > 0 { *s -= 1; break; }
            drop(s);
            std::thread::sleep(std::time::Duration::from_millis(50));
        }
        let setup_c = Arc::clone(&setup_arc);
        let slots_c = Arc::clone(&slots);
        let h = thread::Builder::new()
            .stack_size(STACK_SIZE)
            .spawn(move || {
                println!("  Swarm {}/{}", sid+1, swarms);
                simulate("", &*setup_c, &dat, steps, parsed_method);
                *slots_c.lock().unwrap() += 1;
            }).unwrap();
        handles.push(h);
    }
    for h in handles { let _ = h.join(); }

    // Step 3: Generate conformations per swarm
    println!("\n--- Step 3: Generate conformations ---");
    let rec_ld = format!("lightdock_{}", Path::new(&rec).file_name().unwrap().to_str().unwrap());
    let lig_ld = format!("lightdock_{}", Path::new(&lig).file_name().unwrap().to_str().unwrap());
    for sid in 0..swarms {
        let gso = format!("swarm_{}/gso_{}.out", sid, steps);
        if !Path::new(&gso).exists() { continue; }
        cmd_generate(&vec![rec_ld.clone(), lig_ld.clone(), gso, glowworms.to_string()]);
    }

    // Step 4: Cluster (optional)
    if do_cluster {
        println!("\n--- Step 4: Cluster ---");
        for sid in 0..swarms {
            let gso = format!("swarm_{}/gso_{}.out", sid, steps);
            if Path::new(&gso).exists() { cmd_cluster(&vec![gso]); }
        }
    }

    // Step 5: Rank
    println!("\n--- Step 5: Rank ---");
    let rank_args: Vec<String> = if do_cluster {
        vec![swarms.to_string(), steps.to_string(), "--filter-clusters".into()]
    } else {
        vec![swarms.to_string(), steps.to_string()]
    };
    cmd_rank(&rank_args);

    // Step 6: Top
    println!("\n--- Step 6: Top {} ---", top_n);
    cmd_top(&vec!["rank_by_scoring.list".into(), top_n.to_string()]);

    println!("\n=== Pipeline complete. top_{}.list ready ===", top_n);
}
/// lgd_gso_to_csv: convert a ranking/GSO output to CSV
fn cmd_gso_to_csv(args: &[String]) {
    if args.len() < 2 {
        eprintln!("Usage: lklight gso_to_csv <ranking_file> <output.csv> [--sep ',']");
        return;
    }
    let ranking = &args[0];
    let csv_out = &args[1];
    let sep = args.windows(2).find(|w| w[0] == "--sep")
        .map(|w| w[1].as_str()).unwrap_or(",");
    let content = fs::read_to_string(ranking).expect("Cannot read ranking file");
    let mut out  = fs::File::create(csv_out).unwrap();
    let mut count = 0usize;
    for line in content.lines() {
        if line.starts_with('#') {
            let hdr = line.trim_start_matches('#').trim();
            writeln!(out, "{}", hdr.split_whitespace().collect::<Vec<_>>().join(sep)).unwrap();
            continue;
        }
        let fields: Vec<&str> = line.split_whitespace().collect();
        if fields.is_empty() { continue; }
        writeln!(out, "{}", fields.join(sep)).unwrap();
        count += 1;
    }
    println!("gso_to_csv: {} lines → {}", count, csv_out);
}

/// move_anm: generate N conformations by randomly sampling ANM modes
fn cmd_move_anm(args: &[String]) {
    if args.len() < 3 {
        eprintln!("Usage: lklight move_anm <pdb_file> <n_modes> <n_confs> [--rmsd 1.5] [--seed N]");
        return;
    }
    let pdb_file  = &args[0];
    let n_modes: usize  = args[1].parse().expect("n_modes must be int");
    let n_confs: usize  = args[2].parse().expect("n_confs must be int");
    let rmsd: f64 = args.windows(2).find(|w| w[0]=="--rmsd")
        .and_then(|w| w[1].parse().ok()).unwrap_or(1.5);
    let seed: u64 = args.windows(2).find(|w| w[0]=="--seed")
        .and_then(|w| w[1].parse().ok()).unwrap_or(324_324);

    use lklight::anm::{build_atom_to_backbone, compute_anm, save_npy};
    use rand::rngs::StdRng;
    use rand::{Rng, SeedableRng};

    // Parse structure with pdbtbx (padded)
    let pdb = open_pdb_padded(pdb_file);
    let (bb_coords, all_coords, bb_keys, all_keys) = extract_anm_data(&pdb, true);
    let atom_to_bb = build_atom_to_backbone(&bb_keys, &all_keys);
    let modes = compute_anm(&bb_coords, &atom_to_bb, all_coords.len(), n_modes, rmsd);
    let n_all = modes.n_atoms;

    let atoms = parse_pdb(pdb_file);
    let n_atoms = atoms.len().min(n_all);
    let mut rng: StdRng = SeedableRng::seed_from_u64(seed);

    for c in 0..n_confs {
        let mut out_atoms = atoms.clone();
        // Random normal extents
        let extents: Vec<f64> = (0..n_modes).map(|_| rng.gen::<f64>() * 2.0 - 1.0).collect();
        for nm in 0..n_modes {
            let ext = extents[nm] * rmsd;
            let base = nm * n_all * 3;
            for a in 0..n_atoms {
                if base + a*3 + 2 < modes.data.len() {
                    out_atoms[a].x += modes.data[base + a*3    ] * ext;
                    out_atoms[a].y += modes.data[base + a*3 + 1] * ext;
                    out_atoms[a].z += modes.data[base + a*3 + 2] * ext;
                }
            }
        }
        let out_name = format!("anm_{}_{}", c+1, Path::new(pdb_file).file_name()
            .and_then(|n| n.to_str()).unwrap_or("out.pdb"));
        write_pdb(&out_atoms, &out_name);
    }
    let _ = save_npy;
    println!("move_anm: generated {} conformations from {}", n_confs, pdb_file);
}

/// score: compute the scoring function energy for a pre-positioned receptor + ligand pair.
/// Usage: lklight score <rec.pdb> <lig.pdb> <method>
///   [--tx X] [--ty Y] [--tz Z] [--qw W] [--qx X] [--qy Y] [--qz Z]
/// Default: translation=(0,0,0), rotation=identity (both molecules at current positions).
fn cmd_score(args: &[String]) {
    use lklight::simulator::score_pdb;
    if args.len() < 3 {
        eprintln!("Usage: lklight score <rec.pdb> <lig.pdb> <method> [OPTIONS]");
        eprintln!("  --tx X  --ty Y  --tz Z       ligand translation (default 0 0 0)");
        eprintln!("  --qw W  --qx X  --qy Y  --qz Z  ligand rotation quaternion (default identity)");
        eprintln!("  Methods: dfire fastdfire dfire2 dna mj3h pydock cpydock sd vdw pisa sipper tobi ddna");
        return;
    }
    let rec_path = &args[0];
    let lig_path = &args[1];
    let method   = match parse_method(&args[2]) {
        Some(m) => m,
        None => { eprintln!("Unknown method: {}", args[2]); return; }
    };
    let tx: f64 = args.windows(2).find(|w| w[0]=="--tx").and_then(|w| w[1].parse().ok()).unwrap_or(0.0);
    let ty: f64 = args.windows(2).find(|w| w[0]=="--ty").and_then(|w| w[1].parse().ok()).unwrap_or(0.0);
    let tz: f64 = args.windows(2).find(|w| w[0]=="--tz").and_then(|w| w[1].parse().ok()).unwrap_or(0.0);
    let qw: f64 = args.windows(2).find(|w| w[0]=="--qw").and_then(|w| w[1].parse().ok()).unwrap_or(1.0);
    let qx: f64 = args.windows(2).find(|w| w[0]=="--qx").and_then(|w| w[1].parse().ok()).unwrap_or(0.0);
    let qy: f64 = args.windows(2).find(|w| w[0]=="--qy").and_then(|w| w[1].parse().ok()).unwrap_or(0.0);
    let qz: f64 = args.windows(2).find(|w| w[0]=="--qz").and_then(|w| w[1].parse().ok()).unwrap_or(0.0);

    let energy = score_pdb(rec_path, lig_path, tx, ty, tz, qw, qx, qy, qz, method);
    println!("Score ({:?}): {:.6}", method, energy);
}

// =============================================================================
// Subcommand: DIAMETER  (lgd_calculate_diameter equivalent)
// =============================================================================

fn cmd_diameter(args: &[String]) {
    if args.is_empty() {
        eprintln!("Usage: lklight diameter <pdb_file>");
        return;
    }
    let atoms = parse_pdb(&args[0]);
    if atoms.is_empty() { eprintln!("No atoms in {}", args[0]); return; }
    let n = atoms.len();
    let mut max_d2 = 0.0f64;
    let mut pair = (0usize, 0usize);
    for i in 0..n {
        for j in (i+1)..n {
            let d2 = (atoms[i].x - atoms[j].x).powi(2)
                   + (atoms[i].y - atoms[j].y).powi(2)
                   + (atoms[i].z - atoms[j].z).powi(2);
            if d2 > max_d2 { max_d2 = d2; pair = (i, j); }
        }
    }
    let d = max_d2.sqrt();
    println!("Diameter of {}: {:.3} Å  (atoms {} and {})",
        args[0], d, pair.0, pair.1);
}

// =============================================================================
// Subcommand: TRAJECTORY  (lgd_generate_trajectory equivalent)
// =============================================================================

fn cmd_trajectory(args: &[String]) {
    if args.len() < 5 {
        eprintln!("Usage: lklight trajectory <rec.pdb> <lig.pdb> <swarm_id> <glowworm_id> <steps>");
        eprintln!("  Generates trajectory_<gid>_step_<N>.pdb for each step in swarm_<sid>/");
        eprintln!("  Reads setup.json and ANM .npy files from CWD if present.");
        return;
    }
    let rec_atoms_orig = parse_pdb(&args[0]);
    let lig_atoms_orig = parse_pdb(&args[1]);
    let swarm_id: u32  = args[2].parse().expect("swarm_id must be int");
    let glowworm_id: usize = args[3].parse().expect("glowworm_id must be int");
    let max_steps: u32 = args[4].parse().expect("steps must be int");

    let setup     = read_setup_from_file("setup.json").ok();
    let n_anm_rec = setup.as_ref().map(|s| s.anm_rec).unwrap_or(0);
    let n_anm_lig = setup.as_ref().map(|s| s.anm_lig).unwrap_or(0);
    let modes_rec = load_npy_f64("lightdock_rec.nm.npy");
    let modes_lig = load_npy_f64("lightdock_lig.nm.npy");

    let mut found = 0u32;
    for step in 0..=max_steps {
        let gso_path = format!("swarm_{}/gso_{}.out", swarm_id, step);
        if !Path::new(&gso_path).exists() { continue; }
        let entries = parse_gso_file(&gso_path);
        // find the glowworm entry by id
        let e = match entries.iter().find(|e| e.glowworm_id == glowworm_id) {
            Some(e) => e.clone(),
            None => continue,
        };
        let mut rec_pose = rec_atoms_orig.clone();
        if let Some(ref modes) = modes_rec {
            if n_anm_rec > 0 && !e.all_extents.is_empty() {
                apply_anm(&mut rec_pose, modes, &e.all_extents[..n_anm_rec.min(e.all_extents.len())]);
            }
        }
        let mut lig_pose = lig_atoms_orig.clone();
        if let Some(ref modes) = modes_lig {
            if n_anm_lig > 0 && e.all_extents.len() > n_anm_rec {
                let start = n_anm_rec;
                let end   = start + n_anm_lig.min(e.all_extents.len().saturating_sub(start));
                apply_anm(&mut lig_pose, modes, &e.all_extents[start..end]);
            }
        }
        let rec_len = rec_pose.len() as u32;
        for (j, atom) in lig_pose.iter_mut().enumerate() {
            let r = e.rotation.rotate([atom.x, atom.y, atom.z]);
            atom.x = r[0] + e.translation[0];
            atom.y = r[1] + e.translation[1];
            atom.z = r[2] + e.translation[2];
            atom.serial = rec_len + j as u32 + 1;
        }
        let mut all = rec_pose; all.extend(lig_pose);
        let out = format!("trajectory_{}_step_{}.pdb", glowworm_id, step);
        write_pdb(&all, &out);
        found += 1;
    }
    println!("trajectory: generated {} PDB frames for glowworm {} in swarm {}",
        found, glowworm_id, swarm_id);
}

// =============================================================================
// Subcommand: MAP_CONTACTS  (lgd_map_contacts equivalent)
// =============================================================================

fn cmd_map_contacts(args: &[String]) {
    if args.len() < 3 {
        eprintln!("Usage: lklight map_contacts <rec.pdb> <lig.pdb> <gso_file> [--cutoff 5.0]");
        eprintln!("  Writes rec_contacts.pdb with B-factors = contact frequency (0-1).");
        return;
    }
    let rec_atoms = parse_pdb(&args[0]);
    let lig_atoms = parse_pdb(&args[1]);
    let gso_file  = &args[2];
    let cutoff: f64 = args.windows(2).find(|w| w[0]=="--cutoff")
        .and_then(|w| w[1].parse().ok()).unwrap_or(5.0);
    let cutoff2 = cutoff * cutoff;

    let setup     = read_setup_from_file("setup.json").ok();
    let n_anm_rec = setup.as_ref().map(|s| s.anm_rec).unwrap_or(0);
    let n_anm_lig = setup.as_ref().map(|s| s.anm_lig).unwrap_or(0);
    let modes_rec = load_npy_f64("lightdock_rec.nm.npy");
    let modes_lig = load_npy_f64("lightdock_lig.nm.npy");

    let entries = parse_gso_file(gso_file);
    if entries.is_empty() { eprintln!("No entries in {}", gso_file); return; }

    // Contact counter per receptor atom index
    let mut contact_count = vec![0u32; rec_atoms.len()];

    for e in &entries {
        let mut rec_pose = rec_atoms.clone();
        if let Some(ref modes) = modes_rec {
            if n_anm_rec > 0 && !e.all_extents.is_empty() {
                apply_anm(&mut rec_pose, modes, &e.all_extents[..n_anm_rec.min(e.all_extents.len())]);
            }
        }
        let mut lig_pose = lig_atoms.clone();
        if let Some(ref modes) = modes_lig {
            if n_anm_lig > 0 && e.all_extents.len() > n_anm_rec {
                let start = n_anm_rec;
                let end   = start + n_anm_lig.min(e.all_extents.len().saturating_sub(start));
                apply_anm(&mut lig_pose, modes, &e.all_extents[start..end]);
            }
        }
        for la in &lig_pose {
            let lx = la.x + e.translation[0];
            let ly = la.y + e.translation[1];
            let lz = la.z + e.translation[2];
            // Apply rotation
            let lr = e.rotation.rotate([lx, ly, lz]);
            for (ri, ra) in rec_pose.iter().enumerate() {
                let d2 = (ra.x-lr[0]).powi(2) + (ra.y-lr[1]).powi(2) + (ra.z-lr[2]).powi(2);
                if d2 <= cutoff2 { contact_count[ri] += 1; }
            }
        }
    }
    let n_entries = entries.len() as f64;
    let mut out_atoms = rec_atoms.clone();
    for (i, a) in out_atoms.iter_mut().enumerate() {
        a.temp_factor = contact_count[i] as f64 / n_entries; // 0..1
    }
    let out_path = Path::new(gso_file).with_extension("").to_string_lossy().to_string() + "_contacts.pdb";
    write_pdb(&out_atoms, &out_path);
    println!("map_contacts: {} glowworms, cutoff={:.1}Å → {}", entries.len(), cutoff, out_path);
}

// =============================================================================
// Subcommand: REFERENCE_POINTS  (lgd_calculate_reference_points equivalent)
// =============================================================================

fn cmd_reference_points(args: &[String]) {
    if args.is_empty() {
        eprintln!("Usage: lklight reference_points <pdb_file> [--save]");
        eprintln!("  Computes center of mass + 6 PCA poles of the structure.");
        eprintln!("  With --save: writes lightdock_<pdb>.npy containing the center.");
        return;
    }
    let atoms = parse_pdb(&args[0]);
    if atoms.is_empty() { eprintln!("No atoms in {}", args[0]); return; }
    let n = atoms.len() as f64;

    // Centre of mass
    let cx = atoms.iter().map(|a| a.x).sum::<f64>() / n;
    let cy = atoms.iter().map(|a| a.y).sum::<f64>() / n;
    let cz = atoms.iter().map(|a| a.z).sum::<f64>() / n;
    println!("Center of mass: [{:.3}, {:.3}, {:.3}]", cx, cy, cz);

    // 3x3 covariance matrix
    let (mut cxx,mut cxy,mut cxz,mut cyy,mut cyz,mut czz) = (0.0f64,0.,0.,0.,0.,0.);
    for a in &atoms {
        let (dx,dy,dz) = (a.x-cx, a.y-cy, a.z-cz);
        cxx+=dx*dx; cxy+=dx*dy; cxz+=dx*dz;
        cyy+=dy*dy; cyz+=dy*dz; czz+=dz*dz;
    }
    let (cxx,cxy,cxz,cyy,cyz,czz) = (cxx/n, cxy/n, cxz/n, cyy/n, cyz/n, czz/n);

    // Analytical eigenvalues for 3x3 symmetric matrix via Jacobi iteration (3 sweeps sufficient)
    let mut mat = [[cxx,cxy,cxz],[cxy,cyy,cyz],[cxz,cyz,czz]];
    let mut evec = [[1.0f64,0.,0.],[0.,1.,0.],[0.,0.,1.]];
    for _ in 0..50 {
        for p in 0..3 { for q in (p+1)..3 {
            if mat[p][q].abs() < 1e-12 { continue; }
            let theta = 0.5 * (mat[q][q]-mat[p][p]) / mat[p][q];
            let t = if theta >= 0.0 { 1.0/(theta+(1.0+theta*theta).sqrt()) }
                    else { -1.0/(-theta+(1.0+theta*theta).sqrt()) };
            let c = 1.0/(1.0+t*t).sqrt(); let s = t*c;
            // update mat
            let app = mat[p][p]; let aqq = mat[q][q]; let apq = mat[p][q];
            mat[p][p] = app - t*apq; mat[q][q] = aqq + t*apq; mat[p][q] = 0.0; mat[q][p] = 0.0;
            for r in 0..3 { if r!=p && r!=q {
                let arp = mat[r][p]; let arq = mat[r][q];
                mat[r][p] = c*arp - s*arq; mat[p][r] = mat[r][p];
                mat[r][q] = s*arp + c*arq; mat[q][r] = mat[r][q];
            }}
            for r in 0..3 {
                let erp = evec[r][p]; let erq = evec[r][q];
                evec[r][p] = c*erp - s*erq;
                evec[r][q] = s*erp + c*erq;
            }
        }}
    }
    // eigenvalues on diagonal of mat, eigenvectors as columns of evec
    for i in 0..3 {
        let eval = mat[i][i];
        let radius = eval.max(0.0).sqrt();  // principal radius = sqrt(variance)
        let ev = [evec[0][i], evec[1][i], evec[2][i]];
        let norm = (ev[0]*ev[0]+ev[1]*ev[1]+ev[2]*ev[2]).sqrt().max(1e-12);
        let ev = [ev[0]/norm, ev[1]/norm, ev[2]/norm];
        println!("  PC{}: eigenvalue={:.4}  radius={:.3} Å  axis=[{:.4},{:.4},{:.4}]",
            i+1, eval, radius, ev[0], ev[1], ev[2]);
        println!("       pole+: [{:.3},{:.3},{:.3}]",
            cx+radius*ev[0], cy+radius*ev[1], cz+radius*ev[2]);
        println!("       pole-: [{:.3},{:.3},{:.3}]",
            cx-radius*ev[0], cy-radius*ev[1], cz-radius*ev[2]);
    }
    if args.iter().any(|a| a=="--save") {
        let base = Path::new(&args[0]).file_name().and_then(|n| n.to_str()).unwrap_or("struct");
        let out  = format!("lightdock_{}.npy", base);
        let data = vec![cx, cy, cz];
        let _ = lklight::anm::save_npy(&out, &data);
        println!("Saved center to {}", out);
    }
}

fn print_help(prog: &str) {
    eprintln!("LKlight — unified docking tool\n");
    eprintln!("Usage: {} <subcommand> [args]\n", prog);
    eprintln!("Subcommands:");
    eprintln!("  setup            <rec.pdb> <lig.pdb> [-s N] [-g N] [--anm] [--anm-rec-rmsd R] [--restraints F]");
    eprintln!("  run              <setup.json> <initial_positions.dat> <steps> <method>");
    eprintln!("  generate         <rec.pdb> <lig.pdb> <gso_output> <N>    (ANM-aware)");
    eprintln!("  cluster          <gso_output.out> [--cutoff 4.0]");
    eprintln!("  rank             <num_swarms> <steps> [--filter-clusters] [--clashes-cutoff S]");
    eprintln!("  rank_swarm       <num_swarms> <steps>                     (per-swarm rank)");
    eprintln!("  top              <ranking_file> <N>                       (writes top_K.pdb)");
    eprintln!("  filter           <ranking_file> <restraints.list> [--rec-cutoff 0.4]");
    eprintln!("  gso_to_csv       <ranking_file> <output.csv> [--sep ',']");
    eprintln!("  move_anm         <pdb> <n_modes> <n_confs> [--rmsd 1.5]");
    eprintln!("  score            <rec.pdb> <lig.pdb> <method> [--tx X --ty Y --tz Z --qw W ...]");
    eprintln!("  diameter         <pdb_file>");
    eprintln!("  trajectory       <rec.pdb> <lig.pdb> <swarm_id> <glowworm_id> <steps>");
    eprintln!("  map_contacts     <rec.pdb> <lig.pdb> <gso_file> [--cutoff 5.0]");
    eprintln!("  reference_points <pdb_file> [--save]");
    eprintln!("  pipeline         <rec.pdb> <lig.pdb> <method> [--threads N] [--restraints F]");
    eprintln!("\nMethods: dfire fastdfire dfire2 dna mj3h pydock cpydock sd vdw pisa sipper tobi ddna");
}

fn dispatch() {
    env_logger::init();
    let args: Vec<String> = env::args().collect();
    if args.len() < 2 { print_help(&args[0]); return; }
    let rest = &args[2..].to_vec();
    match args[1].as_str() {
        "setup"    => cmd_setup(rest),
        "run"      => cmd_run(rest),
        "generate" => cmd_generate(rest),
        "cluster"  => cmd_cluster(rest),
        "rank"     => cmd_rank(rest),
        "top"      => cmd_top(rest),
        "filter"     => cmd_filter(rest),
        "gso_to_csv"       => cmd_gso_to_csv(rest),
        "rank_swarm"       => cmd_rank_swarm(rest),
        "move_anm"         => cmd_move_anm(rest),
        "score"            => cmd_score(rest),
        "diameter"         => cmd_diameter(rest),
        "trajectory"       => cmd_trajectory(rest),
        "map_contacts"     => cmd_map_contacts(rest),
        "reference_points" => cmd_reference_points(rest),
        "pipeline"         => cmd_pipeline(rest),
        "--help"|"-h"|"help" => print_help(&args[0]),
        other => {
            eprintln!("Unknown subcommand: '{}'\n", other);
            print_help(&args[0]);
        }
    }
}

fn main() {
    let child = thread::Builder::new()
        .stack_size(STACK_SIZE)
        .spawn(dispatch)
        .unwrap();
    child.join().unwrap();
}
