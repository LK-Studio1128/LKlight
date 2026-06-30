/// cpydock.rs — PyDock scoring with contact-based SASA desolvation energy.
///
/// Equivalent to the Python cpydock scoring function (cpydock/energy/c/cpydock.c).
/// Energy = (elec + 0.1 × vdw + solv) × -1.0
/// where solv = -(solv_rec + solv_lig) (desolvation is favourable when buried).
///
/// Desolvation formula per atom i (receptor):
///   d = min distance to any heavy ligand atom
///   if d ≤ 6.4 Å and asa[i] > 0:
///       solv = min(-10×d + 65, asa[i])
///       total_solv_rec += solv × des_energy[i]

use super::amber::{AMBER_TYPES, ELE_CHARGES, NT_ELE_CHARGES, VDW_CHARGES, VDW_RADII};
use super::constants::{INTERFACE_CUTOFF2, MEMBRANE_PENALTY_SCORE};
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
const SOLVATION_DISTANCE2: f64 = 6.4 * 6.4; // 40.96 Å²
const VDW_WEIGHT: f64 = 0.1;

// ── Desolvation energy coefficient by AMBER type (charges_per_asp_amber) ─────
fn des_energy_by_amber(amber_type: &str) -> f64 {
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

// ── Reference SASA per (residue, atom) from solvation.py reference_area ──────
fn reference_asa(res_name: &str, atom_name: &str) -> Option<f64> {
    let v = match res_name {
        "ALA" => match atom_name { "C"=>1.68,"CB"=>56.29,"CA"=>6.32,"O"=>23.84,"N"=>3.04, _=>return None },
        "ARG" => match atom_name { "C"=>1.78,"CB"=>22.95,"CA"=>6.98,"CG"=>30.21,"O"=>24.79,"CD"=>41.17,"N"=>4.16,"CZ"=>19.85,"NE"=>34.45,"NH1"=>60.12,"NH2"=>61.74, _=>return None },
        "ASN" => match atom_name { "C"=>3.74,"CB"=>28.23,"CA"=>5.0,"CG"=>13.03,"O"=>25.02,"N"=>4.86,"OD1"=>30.66,"ND2"=>53.44, _=>return None },
        "ASP" => match atom_name { "C"=>4.43,"CB"=>31.86,"CA"=>6.95,"CG"=>20.76,"O"=>23.38,"N"=>5.36,"OD1"=>23.23,"OD2"=>36.41, _=>return None },
        "CYS" => match atom_name { "C"=>3.61,"CB"=>33.83,"CA"=>7.23,"O"=>24.54,"N"=>4.49,"SG"=>64.53, _=>return None },
        "GLN" => match atom_name { "C"=>3.5,"CB"=>18.83,"CA"=>6.97,"CG"=>28.65,"O"=>24.86,"CD"=>15.46,"N"=>5.66,"NE2"=>52.85,"OE1"=>34.89, _=>return None },
        "GLU" => match atom_name { "C"=>2.7,"CB"=>23.63,"CA"=>7.1,"CG"=>32.13,"O"=>25.99,"CD"=>23.66,"OE2"=>34.58,"N"=>2.9,"OE1"=>29.9, _=>return None },
        "GLY" => match atom_name { "CA"=>40.04,"C"=>8.29,"O"=>27.15,"N"=>13.87, _=>return None },
        "HIS"|"HID"|"HIE"|"HIP"|"HSC"|"HSE" => match atom_name { "C"=>4.77,"CE1"=>52.36,"CB"=>26.63,"CA"=>4.88,"CG"=>4.01,"O"=>25.71,"N"=>2.88,"CD2"=>37.62,"ND1"=>19.57,"NE2"=>27.31, _=>return None },
        "ILE" => match atom_name { "C"=>2.52,"CB"=>7.86,"CA"=>6.22,"O"=>22.93,"N"=>4.52,"CD1"=>69.96,"CG1"=>23.65,"CG2"=>51.08, _=>return None },
        "LEU" => match atom_name { "C"=>2.05,"CB"=>24.49,"CA"=>6.34,"CG"=>10.99,"O"=>23.43,"CD1"=>64.19,"CD2"=>67.77,"N"=>4.31, _=>return None },
        "LYS" => match atom_name { "C"=>2.11,"CB"=>23.8,"CA"=>9.32,"CG"=>21.55,"CE"=>42.83,"CD"=>28.36,"NZ"=>59.79,"O"=>23.7,"N"=>6.55, _=>return None },
        "MET" => match atom_name { "C"=>2.73,"CB"=>26.63,"CA"=>6.79,"CG"=>28.25,"O"=>23.26,"SD"=>51.12,"CE"=>62.36,"N"=>4.27, _=>return None },
        "PHE" => match atom_name { "C"=>2.12,"CZ"=>37.78,"CB"=>21.04,"CA"=>5.31,"CG"=>3.18,"O"=>23.86,"N"=>3.51,"CE1"=>37.4,"CE2"=>37.5,"CD1"=>22.73,"CD2"=>22.7, _=>return None },
        "PRO" => match atom_name { "C"=>1.27,"CB"=>37.06,"CA"=>11.0,"CG"=>41.44,"O"=>18.32,"CD"=>31.15,"N"=>0.53, _=>return None },
        "SER" => match atom_name { "C"=>4.94,"OG"=>31.08,"CB"=>46.15,"CA"=>9.06,"O"=>23.7,"N"=>6.44, _=>return None },
        "THR" => match atom_name { "C"=>3.56,"CB"=>15.51,"CA"=>6.38,"OG1"=>30.13,"O"=>24.81,"N"=>4.58,"CG2"=>62.7, _=>return None },
        "TRP" => match atom_name { "C"=>4.13,"CZ2"=>37.81,"CB"=>25.13,"CA"=>2.3,"CG"=>3.01,"CH2"=>38.2,"O"=>24.82,"N"=>2.58,"CE2"=>6.66,"CE3"=>21.85,"CD1"=>32.84,"CD2"=>2.78,"CZ3"=>37.25,"NE1"=>25.66, _=>return None },
        "TYR" => match atom_name { "C"=>2.07,"CZ"=>17.63,"CB"=>22.44,"CA"=>5.54,"CG"=>2.8,"O"=>23.87,"OH"=>38.04,"N"=>3.67,"CE1"=>37.67,"CE2"=>37.62,"CD1"=>23.04,"CD2"=>23.12, _=>return None },
        "VAL" => match atom_name { "C"=>1.68,"CB"=>9.17,"CA"=>5.8,"O"=>23.84,"N"=>3.04,"CG1"=>58.93,"CG2"=>59.67, _=>return None },
        "ABU" => match atom_name { "C"=>1.68,"CB"=>9.17,"CA"=>5.8,"CG"=>59.67,"O"=>23.84,"N"=>3.04, _=>return None },
        _ => return None,
    };
    Some(v)
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
    pub is_heavy: Vec<bool>,   // true if not hydrogen
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
            is_heavy: Vec::new(),
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

                    // Hydrogen check: atom name starts with 'H'
                    let heavy = !atom_name.starts_with('H');
                    model.is_heavy.push(heavy);

                    // Desolvation energy coefficient
                    model.des_energy.push(des_energy_by_amber(amber_type));

                    // Reference SASA
                    let asa_val = if !heavy {
                        -1.0  // hydrogen → excluded
                    } else {
                        reference_asa(res_name, atom_name).unwrap_or(100.0)
                    };
                    model.asa.push(asa_val);

                    model.coordinates.push([atom.x(), atom.y(), atom.z()]);
                    atom_index += 1;
                }
            }
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
            const HUGE: f64 = 1.0e9;
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

            // ── Phase 1: parallel ELEC + VDW (receptor atoms in parallel) ────────
            let rec_ele  = &self.receptor.ele_charges;
            let lig_ele  = &self.ligand.ele_charges;
            let rec_svdw = &self.receptor.sqrt_vdw_charges;
            let lig_svdw = &self.ligand.sqrt_vdw_charges;
            let rec_vdwr = &self.receptor.vdw_radii;
            let lig_vdwr = &self.ligand.vdw_radii;
            let lig_slice: &[[f64; 3]] = lig_c.as_slice();

            let (total_elec_raw, total_vdw) = rec_c.iter().enumerate()
                .map(|(i, ra)| {
                    let rx = ra[0]; let ry = ra[1]; let rz = ra[2];
                    let mut ei = 0.0f64;
                    let mut vi = 0.0f64;
                    for (j, la) in lig_slice.iter().enumerate() {
                        let dx = rx - la[0]; let dy = ry - la[1]; let dz = rz - la[2];
                        let d2 = dx*dx + dy*dy + dz*dz;
                        if d2 <= ELEC_DIST_CUTOFF2 {
                            let ae = (rec_ele[i] * lig_ele[j] / d2)
                                .clamp(ELEC_MIN_CUTOFF, ELEC_MAX_CUTOFF);
                            ei += ae;
                        }
                        if d2 <= VDW_DIST_CUTOFF2 {
                            let vdw_e = rec_svdw[i] * lig_svdw[j];
                            let vdw_r = rec_vdwr[i] + lig_vdwr[j];
                            let p6 = vdw_r.powi(6) / d2.powi(3);
                            vi += (vdw_e * (p6*p6 - 2.0*p6)).min(VDW_CUTOFF);
                        }
                    }
                    (ei, vi)
                })
                .fold((0.0, 0.0), |(e1, v1), (e2, v2)| (e1+e2, v1+v2));

            let total_elec = total_elec_raw * FACTOR / EPSILON;

            // ── Phase 2: min_rec, min_lig, interface flags (sequential) ───────
            for (i, ra) in rec_c.iter().enumerate() {
                let rec_heavy = self.receptor.is_heavy[i];
                for (j, la) in lig_slice.iter().enumerate() {
                    let dx = ra[0]-la[0]; let dy = ra[1]-la[1]; let dz = ra[2]-la[2];
                    let d2 = dx*dx + dy*dy + dz*dz;
                    if rec_heavy && self.ligand.is_heavy[j] {
                        if d2 < min_rec[i] { min_rec[i] = d2; }
                        if d2 < min_lig[j] { min_lig[j] = d2; }
                    }
                    if d2 <= INTERFACE_CUTOFF2 {
                        iface_r[i] = 1;
                        iface_l[j] = 1;
                    }
                }
            }

            // ── Desolvation ─────────────────────────────────────────────────
            let mut total_solv = 0.0_f64;
            for i in 0..rec_n {
                let asa = self.receptor.asa[i];
                if asa > 0.0 && min_rec[i] <= SOLVATION_DISTANCE2 && min_rec[i] > 0.0 {
                    let d = min_rec[i].sqrt();
                    let solv = f64::min(-10.0 * d + 65.0, asa);
                    total_solv += solv * self.receptor.des_energy[i];
                }
            }
            for j in 0..lig_n {
                let asa = self.ligand.asa[j];
                if asa > 0.0 && min_lig[j] <= SOLVATION_DISTANCE2 && min_lig[j] > 0.0 {
                    let d = min_lig[j].sqrt();
                    let solv = f64::min(-10.0 * d + 65.0, asa);
                    total_solv += solv * self.ligand.des_energy[j];
                }
            }

            // score = (elec + 0.1×vdw - solv_rec - solv_lig) × -1
            let score = (total_elec + VDW_WEIGHT * total_vdw - total_solv) * -1.0;

            let perc_r = satisfied_restraints(iface_r, &self.receptor.active_restraints);
            let perc_l = satisfied_restraints(iface_l, &self.ligand.active_restraints);
            let intersection = membrane_intersection(iface_r, &self.receptor.membrane);
            let penalty = if intersection > 0.0 { MEMBRANE_PENALTY_SCORE * intersection } else { 0.0 };

            score + perc_r * score + perc_l * score - penalty
        })
    }
}
