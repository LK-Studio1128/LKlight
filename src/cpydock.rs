/// cpydock.rs — PyDock scoring with contact-based SASA desolvation energy.
///
/// Bit-faithful port of the Python cpydock scoring function
/// (lightdock/scoring/cpydock: driver.py + energy/c/cpydock.c + freesasa).
///
/// Energy = (elec + 0.1×vdw + solv) × -1.0
/// where solv = -(solv_rec + solv_lig).
///
/// Reference-equivalent behaviours reproduced exactly:
///  1. Reference SASA per atom is computed with freesasa's Lee & Richards
///     algorithm on the *unbound* monomer (probe 1.4 Å, 20 slices, radii =
///     desolvation radii), NOT a lookup table.
///  2. Desolvation energy/radius coefficients come from the asp_type
///     (residue, atom) table first, falling back to the AMBER-type table.
///  3. The hydrogen flag array is int64 in Python but read as uint32 by the C
///     extension; the resulting flag view is f(i) = (i even) ? hyd[i/2] : 0.
///     The min-distance update condition mirrors the C source:
///     `!f_rec[i] && !f_lig[j]`.
///  4. min distance initialised to HUGE_DISTANCE = 10000.0 and the
///     solvation window is SOLVATION_DISTANCE2 = 6.4*6.4.

use super::amber::{AMBER_TYPES, ELE_CHARGES, NT_ELE_CHARGES, VDW_CHARGES, VDW_RADII};
use super::constants::{INTERFACE_CUTOFF2, MEMBRANE_PENALTY_SCORE};
use super::lr_sasa::lee_richards_sasa;
use super::qt::{rot3_apply, Quaternion};
use super::scoring::{membrane_intersection, satisfied_restraints, Score};
use log::{info, warn};
use pdbtbx::PDB;
use std::cell::RefCell;
use std::collections::HashMap;

// ── Physical constants (from cpydock.c) ──────────────────────────────────────
const EPSILON: f64 = 4.0;
const FACTOR: f64 = 332.0;
const ELEC_MAX_CUTOFF: f64 = 1.0 * EPSILON / FACTOR;
const ELEC_MIN_CUTOFF: f64 = -1.0 * EPSILON / FACTOR;
const VDW_CUTOFF: f64 = 1.0;
const ELEC_DIST_CUTOFF2: f64 = 30.0 * 30.0;
const VDW_DIST_CUTOFF: f64 = 10.0;
const VDW_DIST_CUTOFF2: f64 = VDW_DIST_CUTOFF * VDW_DIST_CUTOFF;
const SOLVATION_DISTANCE2: f64 = 6.4 * 6.4; // 40.959999999999994 (as in C)
const VDW_WEIGHT: f64 = 0.1;

// ── Desolvation coefficients (solvation.py) ──────────────────────────────────

/// asp_type_charges / asp_type_radius indexed by ASP type
const ASP_CHARGES: [f64; 11] = [
    0.0, 0.01918, 0.1108, -0.0391, -0.12604, -0.06256, -0.04255, -0.03128, -0.06877, 0.02576,
    0.00506,
];
const ASP_RADII: [f64; 11] = [
    0.0, 1.95, 1.8, 1.7, 1.7, 1.7, 1.6, 1.4, 1.4, 2.0, 1.85,
];

/// radius_per_asp: fallback desolvation radius by AMBER type
fn radius_per_asp(amber_type: &str) -> f64 {
    match amber_type {
        "C" | "CD" | "CT" | "CY" | "CZ" => 1.95,
        "C*" | "CA" | "CB" | "CC" | "CK" | "CM" | "CN" | "CQ" | "CR" | "CV" | "CW" => 1.8,
        "N" | "N*" | "N2" | "N3" | "NA" | "NB" | "NC" | "NT" | "NY" => 1.7,
        "O" | "O2" | "O3" => 1.4,
        "OH" | "OS" | "OW" => 1.6,
        "S" => 1.85,
        "SH" => 2.0,
        _ => 0.0,
    }
}

/// asp_type table: (residue, atom) → ASP type index
fn asp_type(res_name: &str, atom_name: &str) -> Option<usize> {
    let t = match res_name {
        "ABU" => match atom_name { "C"|"CB"|"CA"|"CG" => 1, "O" => 7, "N" => 3, _ => return None },
        "AHP" => match atom_name { "C"|"CB"|"CA"|"CG"|"CE"|"CD"|"CZ" => 1, "O" => 7, "N" => 3, _ => return None },
        "AHX" => match atom_name { "C"|"CB"|"CA"|"CG"|"CE"|"CD" => 1, "O" => 7, "N" => 3, _ => return None },
        "ALA" => match atom_name { "CB"|"CA"|"C" => 1, "O" => 7, "N" => 3, _ => return None },
        "APE" => match atom_name { "C"|"CB"|"CA"|"CG"|"CD" => 1, "O" => 7, "N" => 3, _ => return None },
        "ARG" => match atom_name { "C"|"CB"|"CA"|"CG"|"CD"|"CZ" => 1, "NE" => 3, "O" => 7, "NH1"|"NH2" => 5, "N" => 3, _ => return None },
        "ASN" => match atom_name { "C"|"CB"|"CA"|"CG" => 1, "O"|"OD1" => 7, "N"|"ND2" => 3, _ => return None },
        "ASP" => match atom_name { "C"|"CB"|"CA"|"CG" => 1, "O" => 7, "N" => 3, "OD1"|"OD2" => 8, _ => return None },
        "CYS" => match atom_name { "C"|"CB"|"CA" => 1, "O" => 7, "N" => 3, "SG" => 9, _ => return None },
        "GLN" => match atom_name { "C"|"CB"|"CA"|"CG"|"CD" => 1, "O"|"OE1" => 7, "N"|"NE2" => 3, _ => return None },
        "GLU" => match atom_name { "C"|"CB"|"CA"|"CG"|"CD" => 1, "O" => 7, "OE1"|"OE2" => 8, "N" => 3, _ => return None },
        "GLY" => match atom_name { "CA"|"C" => 1, "O" => 7, "N" => 3, _ => return None },
        "HID" => match atom_name { "C"|"CB"|"CA" => 1, "CE1"|"CG"|"CD2" => 2, "O" => 7, "N"|"NE2" => 3, "ND1" => 4, _ => return None },
        "HIE" => match atom_name { "C"|"CB"|"CA" => 1, "CE1"|"CG"|"CD2" => 2, "O" => 7, "N"|"ND1" => 3, "NE2" => 4, _ => return None },
        "HIP" => match atom_name { "C"|"CB"|"CA" => 1, "CE1"|"CG"|"CD2" => 2, "O" => 7, "N" => 3, "ND1"|"NE2" => 4, _ => return None },
        "HIS" => match atom_name { "C"|"CB"|"CA" => 1, "CE1"|"CG"|"CD2" => 2, "O" => 7, "N" => 3, "ND1"|"NE2" => 4, _ => return None },
        "HSC" => match atom_name { "C"|"CB"|"CA" => 1, "CE1"|"CG"|"CD2" => 2, "O" => 7, "N" => 3, "ND1"|"NE2" => 4, _ => return None },
        "HSE" => match atom_name { "C"|"CB"|"CA"|"CG" => 1, "OD"|"O" => 7, "N" => 3, _ => return None },
        "ILE" => match atom_name { "C"|"CB"|"CA"|"CD1"|"CG1"|"CG2" => 1, "O" => 7, "N" => 3, _ => return None },
        "LEU" => match atom_name { "C"|"CB"|"CA"|"CG"|"CD1"|"CD2" => 1, "O" => 7, "N" => 3, _ => return None },
        "LYS" => match atom_name { "C"|"CB"|"CA"|"CG"|"CE"|"CD" => 1, "NZ" => 4, "O" => 7, "N" => 3, _ => return None },
        "MET" => match atom_name { "C"|"CB"|"CA"|"CG"|"CE" => 1, "N" => 3, "O" => 7, "SD" => 10, _ => return None },
        "PHE" => match atom_name { "C"|"CB"|"CA" => 1, "CE1"|"CG"|"CZ"|"CD1"|"CD2"|"CE2" => 2, "O" => 7, "N" => 3, _ => return None },
        "PRO" => match atom_name { "C"|"CB"|"CA"|"CG"|"CD" => 1, "O" => 7, "N" => 3, _ => return None },
        "SER" => match atom_name { "C"|"CB"|"CA" => 1, "OG" => 6, "O" => 7, "N" => 3, _ => return None },
        "THR" => match atom_name { "C"|"CB"|"CA"|"CG2" => 1, "OG1" => 6, "O" => 7, "N" => 3, _ => return None },
        "TRP" => match atom_name { "C"|"CB"|"CA" => 1, "CZ2"|"CG"|"CH2"|"CE2"|"CE3"|"CD1"|"CD2"|"CZ3" => 2, "O" => 7, "N"|"NE1" => 3, _ => return None },
        "TYR" => match atom_name { "C"|"CB"|"CA" => 1, "CE1"|"CG"|"CZ"|"CD1"|"CD2"|"CE2" => 2, "OH" => 6, "O" => 7, "N" => 3, _ => return None },
        "VAL" => match atom_name { "C"|"CB"|"CA"|"CG1"|"CG2" => 1, "O" => 7, "N" => 3, _ => return None },
        _ => return None,
    };
    Some(t)
}

/// Desolvation energy coefficient: asp_type first, then AMBER-type table
/// (mirrors solvation.get_solvation including the OXT/CYX/HSD/HSP/HSE/HID
/// residue-renaming special cases).
fn des_energy(res_name: &str, atom_name: &str, amber_type: &str) -> f64 {
    let (mut rr, mut aa) = (res_name.to_string(), atom_name.to_string());
    if aa == "OXT" { rr = "ASP".into(); aa = "OD1".into(); }
    if rr == "CYX" { rr = "CYS".into(); }
    if rr == "HSD" { rr = "HID".into(); }
    if rr == "HSP" { rr = "HIP".into(); }
    if rr == "HSE" { rr = "HIE".into(); }
    if rr == "HID" { rr = "HIP".into(); }
    if let Some(t) = asp_type(&rr, &aa) {
        return ASP_CHARGES[t];
    }
    match amber_type {
        "CT" | "C" | "CD" | "CZ" | "CY" => 0.01918,
        "CA" | "CB" | "CC" | "CK" | "CM" | "CN" | "CQ" | "CR" | "CV" | "CW" | "C*" => 0.1108,
        "N" | "N*" | "NA" | "NB" | "NC" | "NY" | "NT" | "N2" => -0.0391,
        "N3" => -0.12604,
        "O" => -0.03128,
        "O2" | "O3" => -0.06877,
        "OH" | "OW" | "OS" => -0.04255,
        "S" => 0.00506,
        "SH" => 0.02576,
        _ => 0.0,
    }
}

/// Desolvation radius: asp_type first, then radius_per_asp (same renames).
fn des_radius(res_name: &str, atom_name: &str, amber_type: &str) -> f64 {
    let (mut rr, mut aa) = (res_name.to_string(), atom_name.to_string());
    if aa == "OXT" { rr = "ASP".into(); aa = "OD1".into(); }
    if rr == "CYX" { rr = "CYS".into(); }
    if rr == "HSD" { rr = "HID".into(); }
    if rr == "HSP" { rr = "HIP".into(); }
    if rr == "HSE" { rr = "HIE".into(); }
    if rr == "HID" { rr = "HIP".into(); }
    if let Some(t) = asp_type(&rr, &aa) {
        return ASP_RADII[t];
    }
    radius_per_asp(amber_type)
}

// ── Model ─────────────────────────────────────────────────────────────────────

pub struct CPYDOCKDockingModel {
    pub atoms: Vec<usize>,
    pub coordinates: Vec<[f64; 3]>,
    pub membrane: Vec<usize>,
    pub active_restraints: HashMap<String, Vec<usize>>,
    pub passive_restraints: HashMap<String, Vec<usize>>,
    pub num_anm: usize,
    pub nmodes: Vec<f64>,
    pub vdw_radii: Vec<f64>,
    pub vdw_charges: Vec<f64>,
    pub sqrt_vdw_charges: Vec<f64>,
    pub ele_charges: Vec<f64>,
    pub des_energy: Vec<f64>,  // desolvation energy coefficient per atom
    pub asa: Vec<f64>,         // reference SASA per atom (-1.0 = hydrogen, excluded)
    pub hydrogens: Vec<i32>,   // 1 = heavy, 0 = hydrogen (int64 semantics, see C read bug)
}

impl CPYDOCKDockingModel {
    pub fn new(
        structure: &PDB,
        active_restraints: &[String],
        passive_restraints: &[String],
        nmodes: &[f64],
        num_anm: usize,
    ) -> CPYDOCKDockingModel {
        let mut model = CPYDOCKDockingModel {
            atoms: Vec::new(),
            coordinates: Vec::new(),
            membrane: Vec::new(),
            active_restraints: HashMap::new(),
            passive_restraints: HashMap::new(),
            nmodes: nmodes.to_owned(),
            num_anm,
            vdw_radii: Vec::new(),
            vdw_charges: Vec::new(),
            sqrt_vdw_charges: Vec::new(),
            ele_charges: Vec::new(),
            des_energy: Vec::new(),
            asa: Vec::new(),
            hydrogens: Vec::new(),
        };

        let mut atom_index: u64 = 0;
        for chain in structure.chains() {
            for residue in chain.residues() {
                let res_name = residue.name().unwrap_or("UNK");
                let mut res_id = format!("{}.{}.{}", chain.id(), res_name, residue.serial_number());
                if let Some(c) = residue.insertion_code() {
                    res_id.push_str(c);
                }

                for atom in residue.atoms() {
                    let rec_atom_type = format!("{}{}", res_name, atom.name());
                    if rec_atom_type == "MMBBJ" {
                        model.membrane.push(atom_index as usize);
                    }

                    if active_restraints.contains(&res_id) {
                        model.active_restraints
                            .entry(res_id.clone())
                            .or_default()
                            .push(atom_index as usize);
                    }
                    if passive_restraints.contains(&res_id) {
                        model.passive_restraints
                            .entry(res_id.clone())
                            .or_default()
                            .push(atom_index as usize);
                    }

                    let atom_name = atom.name().trim();
                    let mut atom_id = format!("{}-{}", res_name, atom_name);

                    // Never panics: unknown atoms fall back to a generic element type,
                    // then to a neutral carbon-like type.
                    let amber_type = match AMBER_TYPES.get(&*atom_id) {
                        Some(&t) => t,
                        _ => {
                            let h_id = format!("{}-H", res_name);
                            if (atom_name == "H1" || atom_name == "H2" || atom_name == "H3")
                                && AMBER_TYPES.contains_key(&*h_id)
                            {
                                atom_id = h_id;
                                AMBER_TYPES[&*atom_id]
                            } else {
                                let elem =
                                    atom_name.chars().next().unwrap_or('C').to_ascii_uppercase();
                                atom_id = format!("*-{}", elem);
                                match AMBER_TYPES.get(&*atom_id) {
                                    Some(&t) => t,
                                    _ => {
                                        warn!("CPYDOCK Warning: Atom [{:?}] not supported, using neutral fallback", atom_id);
                                        "C"
                                    }
                                }
                            }
                        }
                    };

                    let ele_charge = match ELE_CHARGES.get(&*atom_id) {
                        Some(&c) => c,
                        _ => match NT_ELE_CHARGES.get(&*atom_id) {
                            Some(&c) => c,
                            _ => 0.0,
                        },
                    };
                    model.ele_charges.push(ele_charge);

                    let vdw_charge = *VDW_CHARGES.get(amber_type).unwrap_or(&0.086);
                    model.vdw_charges.push(vdw_charge);
                    model.sqrt_vdw_charges.push(vdw_charge.sqrt());

                    let vdw_radius = *VDW_RADII.get(amber_type).unwrap_or(&1.908);
                    model.vdw_radii.push(vdw_radius);

                    // Hydrogen check: element column (PDB) like Python
                    let is_h = atom.element().map(|e| e.symbol() == "H")
                        .unwrap_or_else(|| atom_name.starts_with('H'));
                    let heavy = !is_h;
                    model.hydrogens.push(if heavy { 1 } else { 0 });

                    // Desolvation energy coefficient
                    model.des_energy.push(des_energy(res_name, atom_name, amber_type));

                    // Reference SASA filled below (needs all heavy atoms first)
                    model.asa.push(-1.0); // placeholder; overwritten for heavy atoms

                    model.coordinates.push([atom.x(), atom.y(), atom.z()]);
                    atom_index += 1;
                }
            }
        }

        // ── Reference SASA with freesasa Lee-Richards on the unbound monomer ──
        // Only heavy atoms participate (as in driver.py); radii = des_radii.
        // Collect heavy atoms (with their global index) + desolvation radii.
        let mut hc: Vec<[f64; 3]> = Vec::new();
        let mut heavy_idx: Vec<usize> = Vec::new();
        let mut names: Vec<(String, String)> = Vec::new();
        let mut global_i = 0usize;
        for chain in structure.chains() {
            for residue in chain.residues() {
                let res_name = residue.name().unwrap_or("UNK").to_string();
                for atom in residue.atoms() {
                    let atom_name = atom.name().trim().to_string();
                    let is_h = atom.element().map(|e| e.symbol() == "H")
                        .unwrap_or_else(|| atom_name.starts_with('H'));
                    if is_h { global_i += 1; continue; }
                    hc.push([atom.x(), atom.y(), atom.z()]);
                    heavy_idx.push(global_i);
                    names.push((res_name.clone(), atom_name));
                    global_i += 1;
                }
            }
        }
        // amber types for radius fallback (same logic as main loop)
        let mut hr: Vec<f64> = Vec::with_capacity(names.len());
        for (rn, an) in &names {
            let atom_id = format!("{}-{}", rn, an);
            let at = match AMBER_TYPES.get(&*atom_id) {
                Some(&t) => t.to_string(),
                _ => {
                    let h_id = format!("{}-H", rn);
                    if (an == "H1" || an == "H2" || an == "H3") && AMBER_TYPES.contains_key(&*h_id)
                    {
                        AMBER_TYPES[&*h_id].to_string()
                    } else {
                        let elem = an.chars().next().unwrap_or('C').to_ascii_uppercase();
                        let a2 = format!("*-{}", elem);
                        AMBER_TYPES.get(&*a2).map(|s| s.to_string()).unwrap_or_else(|| "C".to_string())
                    }
                }
            };
            hr.push(des_radius(rn, an, &at));
        }
        let areas = lee_richards_sasa(&hc, &hr);
        for (&pos, &area) in heavy_idx.iter().zip(areas.iter()) {
            model.asa[pos] = area;
        }

        info!("CPYDOCK atoms read: {}", atom_index);
        model
    }
}

// ── Scoring struct ────────────────────────────────────────────────────────────

pub struct CPYDOCK {
    pub receptor: CPYDOCKDockingModel,
    pub ligand: CPYDOCKDockingModel,
    pub use_anm: bool,
}

impl CPYDOCK {
    pub fn new(
        receptor: PDB,
        rec_active: Vec<String>,
        rec_passive: Vec<String>,
        rec_nmodes: Vec<f64>,
        rec_num_anm: usize,
        ligand: PDB,
        lig_active: Vec<String>,
        lig_passive: Vec<String>,
        lig_nmodes: Vec<f64>,
        lig_num_anm: usize,
        use_anm: bool,
    ) -> Box<dyn Score> {
        Box::new(CPYDOCK {
            receptor: CPYDOCKDockingModel::new(&receptor, &rec_active, &rec_passive, &rec_nmodes, rec_num_anm),
            ligand:   CPYDOCKDockingModel::new(&ligand,   &lig_active, &lig_passive, &lig_nmodes, lig_num_anm),
            use_anm,
        })
    }
}

impl Score for CPYDOCK {
    fn energy(
        &self,
        translation: &[f64],
        rotation: &Quaternion,
        rec_nmodes: &[f64],
        lig_nmodes: &[f64],
    ) -> f64 {
        thread_local! {
            static SCRATCH: RefCell<(
                Vec<[f64; 3]>, Vec<[f64; 3]>,
                Vec<usize>, Vec<usize>,
                Vec<f64>, Vec<f64>,  // min_dist2 buffers
            )> = RefCell::new((Vec::new(), Vec::new(), Vec::new(), Vec::new(), Vec::new(), Vec::new()));
        }
        let rot_mat = rotation.to_matrix();

        SCRATCH.with(|sc| {
            let mut sc = sc.borrow_mut();
            let (rec_c, lig_c, iface_r, iface_l, min_rec, min_lig) = &mut *sc;

            let rec_n = self.receptor.coordinates.len();
            let lig_n = self.ligand.coordinates.len();

            if rec_c.len() != rec_n { rec_c.resize(rec_n, [0.0; 3]); }
            if lig_c.len() != lig_n { lig_c.resize(lig_n, [0.0; 3]); }
            if iface_r.len() != rec_n { iface_r.resize(rec_n, 0); }
            if iface_l.len() != lig_n { iface_l.resize(lig_n, 0); }
            if min_rec.len() != rec_n { min_rec.resize(rec_n, 0.0); }
            if min_lig.len() != lig_n { min_lig.resize(lig_n, 0.0); }

            rec_c.copy_from_slice(&self.receptor.coordinates);
            lig_c.copy_from_slice(&self.ligand.coordinates);
            for v in iface_r.iter_mut() { *v = 0; }
            for v in iface_l.iter_mut() { *v = 0; }
            const HUGE: f64 = 10000.0; // HUGE_DISTANCE in cpydock.c
            for v in min_rec.iter_mut() { *v = HUGE; }
            for v in min_lig.iter_mut() { *v = HUGE; }

            // Derive true ANM atom counts from the mode-vector length
            let lig_nm_n = if self.ligand.num_anm > 0 {
                self.ligand.nmodes.len() / (3 * self.ligand.num_anm)
            } else { lig_n };
            let rec_nm_n = if self.receptor.num_anm > 0 {
                self.receptor.nmodes.len() / (3 * self.receptor.num_anm)
            } else { rec_n };

            // Apply ligand motion
            for (i_atom, coord) in lig_c.iter_mut().enumerate() {
                let r = rot3_apply(&rot_mat, *coord);
                coord[0] = r[0] + translation[0];
                coord[1] = r[1] + translation[1];
                coord[2] = r[2] + translation[2];
                if self.use_anm && self.ligand.num_anm > 0 && i_atom < lig_nm_n {
                    for i_nm in 0..self.ligand.num_anm {
                        let b = i_nm * lig_nm_n * 3 + i_atom * 3;
                        coord[0] += self.ligand.nmodes[b]     * lig_nmodes[i_nm];
                        coord[1] += self.ligand.nmodes[b + 1] * lig_nmodes[i_nm];
                        coord[2] += self.ligand.nmodes[b + 2] * lig_nmodes[i_nm];
                    }
                }
            }
            if self.use_anm && self.receptor.num_anm > 0 {
                for (i_atom, coord) in rec_c.iter_mut().enumerate() {
                    if i_atom >= rec_nm_n { break; }
                    for i_nm in 0..self.receptor.num_anm {
                        let b = i_nm * rec_nm_n * 3 + i_atom * 3;
                        coord[0] += self.receptor.nmodes[b]     * rec_nmodes[i_nm];
                        coord[1] += self.receptor.nmodes[b + 1] * rec_nmodes[i_nm];
                        coord[2] += self.receptor.nmodes[b + 2] * rec_nmodes[i_nm];
                    }
                }
            }

            // ── Single (i,j) main loop, order identical to cpydock.c ────────
            // C: for i: for j: min_dist; elec; vdw; interface
            let rec_ele  = &self.receptor.ele_charges;
            let lig_ele  = &self.ligand.ele_charges;
            let rec_vdwq = &self.receptor.vdw_charges;
            let lig_vdwq = &self.ligand.vdw_charges;
            let rec_vdwr = &self.receptor.vdw_radii;
            let lig_vdwr = &self.ligand.vdw_radii;
            let rec_hyd  = &self.receptor.hydrogens;
            let lig_hyd  = &self.ligand.hydrogens;

            // hydrogen flag as seen by the C binary (uint32 view of int64):
            // f(i) = (i even) ? hyd[i/2] : 0
            let rec_flag = |i: usize| -> bool {
                if i % 2 == 0 { rec_hyd[i / 2] != 0 } else { false }
            };
            let lig_flag = |j: usize| -> bool {
                if j % 2 == 0 { lig_hyd[j / 2] != 0 } else { false }
            };

            let mut total_elec = 0.0f64;
            let mut total_vdw = 0.0f64;
            for (i, ra) in rec_c.iter().enumerate() {
                let rx = ra[0]; let ry = ra[1]; let rz = ra[2];
                for (j, la) in lig_c.iter().enumerate() {
                    let dx = rx - la[0]; let dy = ry - la[1]; let dz = rz - la[2];
                    let d2 = dx*dx + dy*dy + dz*dz;

                    // min distance update (C: if(!flag_i && !flag_j))
                    if !rec_flag(i) && !lig_flag(j) {
                        if d2 < min_rec[i] { min_rec[i] = d2; }
                        if d2 < min_lig[j] { min_lig[j] = d2; }
                    }

                    // Electrostatics (C clamps with ifs, same as clamp)
                    if d2 <= ELEC_DIST_CUTOFF2 {
                        let ae = (rec_ele[i] * lig_ele[j] / d2)
                            .clamp(ELEC_MIN_CUTOFF, ELEC_MAX_CUTOFF);
                        total_elec += ae;
                    }

                    // Van der Waals (C: sqrt(a*b), pow(x,6)/pow(y,3))
                    if d2 <= VDW_DIST_CUTOFF2 {
                        let vdw_e = (rec_vdwq[i] * lig_vdwq[j]).sqrt();
                        let vdw_r = rec_vdwr[i] + lig_vdwr[j];
                        let p6 = vdw_r.powf(6.0) / d2.powf(3.0);
                        let k = vdw_e * (p6*p6 - 2.0*p6);
                        total_vdw += if k > VDW_CUTOFF { VDW_CUTOFF } else { k };
                    }

                    if d2 <= INTERFACE_CUTOFF2 {
                        iface_r[i] = 1;
                        iface_l[j] = 1;
                    }
                }
            }

            let total_elec_kcal = total_elec * FACTOR / EPSILON;

            // ── Desolvation (C solvation loop, same order & clamps) ────────
            let rec_asa = &self.receptor.asa;
            let lig_asa = &self.ligand.asa;
            let rec_des = &self.receptor.des_energy;
            let lig_des = &self.ligand.des_energy;

            let mut total_solv = 0.0_f64;
            for i in 0..rec_n {
                let mut solv_rec = 0.0;
                if min_rec[i] <= SOLVATION_DISTANCE2 && min_rec[i] > 0.0 && rec_asa[i] > 0.0 {
                    solv_rec = -10.0 * min_rec[i].sqrt() + 65.0;
                }
                if solv_rec > rec_asa[i] { solv_rec = rec_asa[i]; }
                total_solv += solv_rec * rec_des[i];
            }
            for j in 0..lig_n {
                let mut solv_lig = 0.0;
                if min_lig[j] <= SOLVATION_DISTANCE2 && min_lig[j] > 0.0 && lig_asa[j] > 0.0 {
                    solv_lig = -10.0 * min_lig[j].sqrt() + 65.0;
                }
                if solv_lig > lig_asa[j] { solv_lig = lig_asa[j]; }
                total_solv += solv_lig * lig_des[j];
            }

            // score = (elec + 0.1×vdw - solv_rec - solv_lig) × -1
            let score = (total_elec_kcal + VDW_WEIGHT * total_vdw - total_solv) * -1.0;

            let perc_r = satisfied_restraints(iface_r, &self.receptor.active_restraints);
            let perc_l = satisfied_restraints(iface_l, &self.ligand.active_restraints);
            let intersection = membrane_intersection(iface_r, &self.receptor.membrane);
            let penalty = if intersection > 0.0 { MEMBRANE_PENALTY_SCORE * intersection } else { 0.0 };

            score + perc_r * score + perc_l * score - penalty
        })
    }
}
