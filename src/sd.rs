/// sd.rs — SwarmDock scoring function.
///
/// Reference: SwarmDock and the use of Normal Modes in Protein-Protein Docking
/// https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2996808/
///
/// Equivalent to the Python sd scoring function (sd/energy/c/sd.c).
/// Uses AMBER94 charges and VdW parameters (same as pydock/cpydock).
/// Key differences from cpydock:
///   - No desolvation energy
///   - Distance cutoff 9 Å with cubic switching function between 7-9 Å
///   - VdW cap = 5,000,000 (effectively unlimited per-accumulation)
///   - Elec converted per-pair (not summed then converted)
///   - VdW accumulated per receptor atom across all ligand atoms (clamped)

use super::amber::{AMBER_TYPES, ELE_CHARGES, NT_ELE_CHARGES, VDW_CHARGES, VDW_RADII};
use super::constants::{INTERFACE_CUTOFF2, MEMBRANE_PENALTY_SCORE};
use super::qt::{rot3_apply, Quaternion};
use super::scoring::{membrane_intersection, satisfied_restraints, Score};
use log::{info, warn};
use pdbtbx::PDB;
use std::cell::RefCell;
use std::collections::HashMap;

// ── Constants (from sd.c) ─────────────────────────────────────────────────────
const EPSILON: f64 = 4.0;
const FACTOR: f64 = 332.0;
const CUTON: f64 = 7.0;
const CUTOFF: f64 = 9.0;
const CUTON2: f64 = CUTON * CUTON;   // 49.0
const CUTOFF2: f64 = CUTOFF * CUTOFF; // 81.0
const VDW_CUTOFF: f64 = 5_000_000.0;
/// (CUTOFF² - CUTON²)³
const SWITCH_DENOM: f64 = (CUTOFF2 - CUTON2) * (CUTOFF2 - CUTON2) * (CUTOFF2 - CUTON2);

#[inline(always)]
fn switch_fn(d2: f64) -> f64 {
    let cd2 = CUTOFF2 - d2;
    cd2 * cd2 * (CUTOFF2 + 2.0 * d2 - 3.0 * CUTON2) / SWITCH_DENOM
}

// ── Model ─────────────────────────────────────────────────────────────────────

pub struct SDDockingModel {
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
}

impl SDDockingModel {
    pub fn new(
        structure: &PDB,
        active_restraints: &[String],
        passive_restraints: &[String],
        nmodes: &[f64],
        num_anm: usize,
    ) -> SDDockingModel {
        let mut model = SDDockingModel {
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
        };

        let mut atom_index: u64 = 0;
        for chain in structure.chains() {
            for residue in chain.residues() {
                let raw_res_name = residue.name().unwrap_or("UNK");
                // Match the Python reference: SDAdapter maps HIS -> HID before
                // looking up amber_types / charges ("HIS" itself never appears
                // as a key in the amber table).
                let res_name = if raw_res_name == "HIS" { "HID" } else { raw_res_name };
                let mut res_id = format!("{}.{}.{}", chain.id(), raw_res_name, residue.serial_number());
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
                                        warn!("SD Warning: Atom [{:?}] not supported, using neutral fallback", atom_id);
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

                    model.coordinates.push([atom.x(), atom.y(), atom.z()]);
                    atom_index += 1;
                }
            }
        }
        info!("SD atoms read: {}", atom_index);
        model
    }
}

// ── Scoring struct ────────────────────────────────────────────────────────────

pub struct SD {
    pub receptor: SDDockingModel,
    pub ligand: SDDockingModel,
    pub use_anm: bool,
}

impl SD {
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
        Box::new(SD {
            receptor: SDDockingModel::new(&receptor, &rec_active, &rec_passive, &rec_nmodes, rec_num_anm),
            ligand:   SDDockingModel::new(&ligand,   &lig_active, &lig_passive, &lig_nmodes, lig_num_anm),
            use_anm,
        })
    }
}

impl Score for SD {
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
                HashMap<(i32, i32, i32), Vec<usize>>,
            )> = RefCell::new((Vec::new(), Vec::new(), Vec::new(), Vec::new(), HashMap::new()));
        }
        let rot_mat = rotation.to_matrix();

        SCRATCH.with(|sc| {
            let mut sc = sc.borrow_mut();
            let (rec_c, lig_c, iface_r, iface_l, lig_grid) = &mut *sc;

            let rec_n = self.receptor.coordinates.len();
            let lig_n = self.ligand.coordinates.len();

            if rec_c.len() != rec_n { rec_c.resize(rec_n, [0.0; 3]); }
            if lig_c.len() != lig_n { lig_c.resize(lig_n, [0.0; 3]); }
            if iface_r.len() != rec_n { iface_r.resize(rec_n, 0); }
            if iface_l.len() != lig_n { iface_l.resize(lig_n, 0); }

            rec_c.copy_from_slice(&self.receptor.coordinates);
            lig_c.copy_from_slice(&self.ligand.coordinates);
            for v in iface_r.iter_mut() { *v = 0; }
            for v in iface_l.iter_mut() { *v = 0; }

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

            // ── Spatial grid for ligand (cell = CUTOFF = 9 Å) ────────────────
            lig_grid.clear();
            for (j, la) in lig_c.iter().enumerate() {
                let k = ((la[0] / CUTOFF).floor() as i32,
                         (la[1] / CUTOFF).floor() as i32,
                         (la[2] / CUTOFF).floor() as i32);
                lig_grid.entry(k).or_default().push(j);
            }

            // ── Phase 1: parallel energy (receptor outer, 9Å grid inner) ───
            let rec_ele  = &self.receptor.ele_charges;
            let lig_ele  = &self.ligand.ele_charges;
            let rec_svdw = &self.receptor.sqrt_vdw_charges;
            let lig_svdw = &self.ligand.sqrt_vdw_charges;
            let rec_vdwr = &self.receptor.vdw_radii;
            let lig_vdwr = &self.ligand.vdw_radii;
            let lig_slice: &[[f64; 3]] = lig_c.as_slice();
            let grid_ref: &HashMap<(i32,i32,i32), Vec<usize>> = &*lig_grid;

            let energy: f64 = rec_c.iter().enumerate()
                .map(|(i, ra)| {
                    let x1 = ra[0]; let y1 = ra[1]; let z1 = ra[2];
                    let cx = (x1 / CUTOFF).floor() as i32;
                    let cy = (y1 / CUTOFF).floor() as i32;
                    let cz = (z1 / CUTOFF).floor() as i32;
                    // Bit-exactness with the Python reference requires processing
                    // ligand atoms in ascending index order (floating-point
                    // summation is order-dependent).
                    let mut nb: Vec<usize> = Vec::new();
                    for dx in -1_i32..=1 {
                        for dy in -1_i32..=1 {
                            for dz in -1_i32..=1 {
                                if let Some(js) = grid_ref.get(&(cx+dx, cy+dy, cz+dz)) {
                                    nb.extend_from_slice(js);
                                }
                            }
                        }
                    }
                    nb.sort_unstable();
                    let mut atom_vdw = 0.0f64;
                    let mut ei = 0.0f64;
                    for &j in &nb {
                        let la = &lig_slice[j];
                        let d2 = (x1-la[0]).powi(2) + (y1-la[1]).powi(2) + (z1-la[2]).powi(2);
                        if d2 < CUTOFF2 {
                            let atom_elec = rec_ele[i] * lig_ele[j] / d2 * FACTOR / EPSILON;
                            let vdw_e = rec_svdw[i] * lig_svdw[j];
                            let vdw_r = rec_vdwr[i] + lig_vdwr[j];
                            let p6 = vdw_r.powi(6) / d2.powi(3);
                            atom_vdw += vdw_e * (p6*p6 - 2.0*p6);
                            if atom_vdw > VDW_CUTOFF { atom_vdw = VDW_CUTOFF; }
                            let pair_e = atom_elec + atom_vdw;
                            if d2 < CUTON2 { ei += pair_e; }
                            else { ei += pair_e * switch_fn(d2); }
                        }
                    }
                    ei
                })
                .sum();

            let score = energy * -1.0;

            // ── Phase 2: interface flags (sequential, INTERFACE_CUTOFF=3.9Å) ──
            for (i, ra) in rec_c.iter().enumerate() {
                let x1 = ra[0]; let y1 = ra[1]; let z1 = ra[2];
                let cx = (x1 / CUTOFF).floor() as i32;
                let cy = (y1 / CUTOFF).floor() as i32;
                let cz = (z1 / CUTOFF).floor() as i32;
                for dx in -1_i32..=1 {
                    for dy in -1_i32..=1 {
                        for dz in -1_i32..=1 {
                            if let Some(js) = lig_grid.get(&(cx+dx, cy+dy, cz+dz)) {
                                for &j in js {
                                    let la = &lig_c[j];
                                    let d2 = (x1-la[0]).powi(2) + (y1-la[1]).powi(2) + (z1-la[2]).powi(2);
                                    if d2 <= INTERFACE_CUTOFF2 {
                                        iface_r[i] = 1;
                                        iface_l[j] = 1;
                                    }
                                }
                            }
                        }
                    }
                }
            }

            let perc_r = satisfied_restraints(iface_r, &self.receptor.active_restraints);
            let perc_l = satisfied_restraints(iface_l, &self.ligand.active_restraints);
            let intersection = membrane_intersection(iface_r, &self.receptor.membrane);
            let penalty = if intersection > 0.0 { MEMBRANE_PENALTY_SCORE * intersection } else { 0.0 };

            score + perc_r * score + perc_l * score - penalty
        })
    }
}
