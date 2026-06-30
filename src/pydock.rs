use super::amber::{AMBER_TYPES, ELE_CHARGES, NT_ELE_CHARGES, VDW_CHARGES, VDW_RADII};
use super::constants::{INTERFACE_CUTOFF2, MEMBRANE_PENALTY_SCORE};
use super::qt::{rot3_apply, Quaternion};
use std::cell::RefCell;
use super::scoring::{membrane_intersection, satisfied_restraints, Score};
use pdbtbx::PDB;
use std::collections::HashMap;

use log::{info, warn};

const EPSILON: f64 = 4.0;
const FACTOR: f64 = 332.0;
const MAX_ES_CUTOFF: f64 = 1.0;
const MIN_ES_CUTOFF: f64 = -1.0;
const VDW_CUTOFF: f64 = 1.0;
const ELEC_DIST_CUTOFF: f64 = 30.0;
const ELEC_DIST_CUTOFF2: f64 = ELEC_DIST_CUTOFF * ELEC_DIST_CUTOFF;
const VDW_DIST_CUTOFF: f64 = 10.0;
const VDW_DIST_CUTOFF2: f64 = VDW_DIST_CUTOFF * VDW_DIST_CUTOFF;
const ELEC_MAX_CUTOFF: f64 = MAX_ES_CUTOFF * EPSILON / FACTOR;
const ELEC_MIN_CUTOFF: f64 = MIN_ES_CUTOFF * EPSILON / FACTOR;


pub struct PYDOCKDockingModel {
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

impl<'a> PYDOCKDockingModel {
    pub(crate) fn new(
        structure: &'a PDB,
        active_restraints: &'a [String],
        passive_restraints: &'a [String],
        nmodes: &[f64],
        num_anm: usize,
    ) -> PYDOCKDockingModel {
        let mut model = PYDOCKDockingModel {
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
                let res_name = match residue.name() {
                    Some(name) => name,
                    None => panic!("PDB Parsing Error: Residue name error"),
                };
                let mut res_id = format!("{}.{}.{}", chain.id(), res_name, residue.serial_number());
                if let Some(c) = residue.insertion_code() {
                    res_id.push_str(c);
                }

                for atom in residue.atoms() {
                    // Membrane beads MMB.BJ
                    let rec_atom_type = format!("{}{}", res_name, atom.name());
                    if rec_atom_type == "MMBBJ" {
                        model.membrane.push(atom_index as usize);
                    }

                    if active_restraints.contains(&res_id) {
                        match model.active_restraints.get_mut(&res_id) {
                            Some(atom_indexes) => {
                                atom_indexes.push(atom_index as usize);
                            }
                            None => {
                                model
                                    .active_restraints
                                    .insert(res_id.to_string(), vec![atom_index as usize]);
                            }
                        }
                    }

                    if passive_restraints.contains(&res_id) {
                        match model.passive_restraints.get_mut(&res_id) {
                            Some(atom_indexes) => {
                                atom_indexes.push(atom_index as usize);
                            }
                            None => {
                                model
                                    .passive_restraints
                                    .insert(res_id.to_string(), vec![atom_index as usize]);
                            }
                        }
                    }

                    let atom_name = atom.name().trim();
                    let mut atom_id = format!("{}-{}", res_name, atom_name);

                    // Calculate AMBER type (never panics: unknown atoms fall back to a
                    // generic element type, then to a neutral carbon-like type)
                    let amber_type = match AMBER_TYPES.get(&*atom_id) {
                        Some(&amber) => amber,
                        _ => {
                            let h_id = format!("{}-H", res_name);
                            if (atom_name == "H1" || atom_name == "H2" || atom_name == "H3")
                                && AMBER_TYPES.contains_key(&*h_id)
                            {
                                atom_id = h_id;
                                AMBER_TYPES[&*atom_id]
                            } else {
                                let atom_element =
                                    atom_name.chars().next().unwrap_or('C').to_ascii_uppercase();
                                atom_id = format!("*-{}", atom_element);
                                match AMBER_TYPES.get(&*atom_id) {
                                    Some(&amber) => amber,
                                    _ => {
                                        warn!(
                                            "PYDOCK Warning: Atom [{:?}] not supported, using neutral fallback",
                                            atom_id
                                        );
                                        "C"
                                    }
                                }
                            }
                        }
                    };

                    // Assign electrostatics charge (defaults to 0.0 for unknown atoms)
                    let ele_charge = match ELE_CHARGES.get(&*atom_id) {
                        Some(&charge) => charge,
                        _ => match NT_ELE_CHARGES.get(&*atom_id) {
                            Some(&charge) => charge,
                            _ => 0.0,
                        },
                    };
                    model.ele_charges.push(ele_charge);

                    // Assign VDW charge and radius (carbon-like defaults for unknown types)
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
        info!("Atoms read: {}", atom_index);
        model
    }
}

pub struct PYDOCK {
    pub receptor: PYDOCKDockingModel,
    pub ligand: PYDOCKDockingModel,
    pub use_anm: bool,
}

impl<'a> PYDOCK {
    pub fn new(
        receptor: PDB,
        rec_active_restraints: Vec<String>,
        rec_passive_restraints: Vec<String>,
        rec_nmodes: Vec<f64>,
        rec_num_anm: usize,
        ligand: PDB,
        lig_active_restraints: Vec<String>,
        lig_passive_restraints: Vec<String>,
        lig_nmodes: Vec<f64>,
        lig_num_anm: usize,
        use_anm: bool,
    ) -> Box<dyn Score + 'a> {
        let d = PYDOCK {
            receptor: PYDOCKDockingModel::new(
                &receptor,
                &rec_active_restraints,
                &rec_passive_restraints,
                &rec_nmodes,
                rec_num_anm,
            ),
            ligand: PYDOCKDockingModel::new(
                &ligand,
                &lig_active_restraints,
                &lig_passive_restraints,
                &lig_nmodes,
                lig_num_anm,
            ),
            use_anm,
        };
        Box::new(d)
    }
}

impl Score for PYDOCK {
    fn energy(
        &self,
        translation: &[f64],
        rotation: &Quaternion,
        rec_nmodes: &[f64],
        lig_nmodes: &[f64],
    ) -> f64 {
        thread_local! {
            static SCRATCH: RefCell<(Vec<[f64;3]>, Vec<[f64;3]>, Vec<usize>, Vec<usize>)> =
                RefCell::new((Vec::new(), Vec::new(), Vec::new(), Vec::new()));
        }
        let rot_mat = rotation.to_matrix();

        SCRATCH.with(|sc| {
        let mut sc = sc.borrow_mut();
        let (receptor_coordinates, ligand_coordinates, interface_receptor, interface_ligand) = &mut *sc;
        let rec_num_atoms = self.receptor.coordinates.len();
        let lig_num_atoms = self.ligand.coordinates.len();
        if receptor_coordinates.len() != rec_num_atoms { receptor_coordinates.resize(rec_num_atoms, [0.0;3]); }
        if ligand_coordinates.len() != lig_num_atoms { ligand_coordinates.resize(lig_num_atoms, [0.0;3]); }
        if interface_receptor.len() != rec_num_atoms { interface_receptor.resize(rec_num_atoms, 0); }
        if interface_ligand.len() != lig_num_atoms { interface_ligand.resize(lig_num_atoms, 0); }
        receptor_coordinates.copy_from_slice(&self.receptor.coordinates);
        ligand_coordinates.copy_from_slice(&self.ligand.coordinates);
        for v in interface_receptor.iter_mut() { *v = 0; }
        for v in interface_ligand.iter_mut() { *v = 0; }

        let lig_nm_n = if self.ligand.num_anm > 0 {
            self.ligand.nmodes.len() / (3 * self.ligand.num_anm)
        } else { lig_num_atoms };
        let rec_nm_n = if self.receptor.num_anm > 0 {
            self.receptor.nmodes.len() / (3 * self.receptor.num_anm)
        } else { rec_num_atoms };

        for (i_atom, coordinate) in ligand_coordinates.iter_mut().enumerate() {
            let r = rot3_apply(&rot_mat, *coordinate);
            coordinate[0] = r[0] + translation[0];
            coordinate[1] = r[1] + translation[1];
            coordinate[2] = r[2] + translation[2];
            if self.use_anm && self.ligand.num_anm > 0 && i_atom < lig_nm_n {
                for i_nm in 0..self.ligand.num_anm {
                    let b = i_nm * lig_nm_n * 3 + i_atom * 3;
                    coordinate[0] += self.ligand.nmodes[b]   * lig_nmodes[i_nm];
                    coordinate[1] += self.ligand.nmodes[b+1] * lig_nmodes[i_nm];
                    coordinate[2] += self.ligand.nmodes[b+2] * lig_nmodes[i_nm];
                }
            }
        }
        if self.use_anm && self.receptor.num_anm > 0 {
            for (i_atom, coordinate) in receptor_coordinates.iter_mut().enumerate() {
                if i_atom >= rec_nm_n { break; }
                for i_nm in 0..self.receptor.num_anm {
                    let b = i_nm * rec_nm_n * 3 + i_atom * 3;
                    coordinate[0] += self.receptor.nmodes[b]   * rec_nmodes[i_nm];
                    coordinate[1] += self.receptor.nmodes[b+1] * rec_nmodes[i_nm];
                    coordinate[2] += self.receptor.nmodes[b+2] * rec_nmodes[i_nm];
                }
            }
        }

        // ── Phase 1: parallel pairwise energy (receptor atoms in parallel) ───
        let rec_ele  = &self.receptor.ele_charges;
        let lig_ele  = &self.ligand.ele_charges;
        let rec_svdw = &self.receptor.sqrt_vdw_charges;
        let lig_svdw = &self.ligand.sqrt_vdw_charges;
        let rec_vdwr = &self.receptor.vdw_radii;
        let lig_vdwr = &self.ligand.vdw_radii;
        let lig_slice: &[[f64; 3]] = ligand_coordinates.as_slice();

        let (total_elec_raw, total_vdw) = receptor_coordinates.iter().enumerate()
            .map(|(i, ra)| {
                let rx = ra[0]; let ry = ra[1]; let rz = ra[2];
                let mut ei = 0.0f64;
                let mut vi = 0.0f64;
                for (j, la) in lig_slice.iter().enumerate() {
                    let dx = rx - la[0];
                    let dy = ry - la[1];
                    let dz = rz - la[2];
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
        let score = -(total_elec + total_vdw);

        // ── Phase 2: interface flags (INTERFACE_CUTOFF=3.9Å, fast sequential) ─
        for (i, ra) in receptor_coordinates.iter().enumerate() {
            for (j, la) in lig_slice.iter().enumerate() {
                let dx = ra[0]-la[0]; let dy = ra[1]-la[1]; let dz = ra[2]-la[2];
                if dx*dx + dy*dy + dz*dz <= INTERFACE_CUTOFF2 {
                    interface_receptor[i] = 1;
                    interface_ligand[j] = 1;
                }
            }
        }

        let perc_r = satisfied_restraints(interface_receptor, &self.receptor.active_restraints);
        let perc_l = satisfied_restraints(interface_ligand, &self.ligand.active_restraints);
        let intersection = membrane_intersection(interface_receptor, &self.receptor.membrane);
        let penalty = if intersection > 0.0 { MEMBRANE_PENALTY_SCORE * intersection } else { 0.0 };

        score + perc_r * score + perc_l * score - penalty
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::qt::Quaternion;
    use std::env;

    #[test]
    fn test_1azp() {
        let cargo_path = match env::var("CARGO_MANIFEST_DIR") {
            Ok(val) => val,
            Err(_) => String::from("."),
        };
        let test_path: String = format!("{}/tests/1azp", cargo_path);

        let receptor_filename: String = format!("{}/1azp_receptor.pdb", test_path);
        let (receptor, _errors) =
            pdbtbx::open(&receptor_filename, pdbtbx::StrictnessLevel::Strict).unwrap();

        let ligand_filename: String = format!("{}/1azp_ligand.pdb", test_path);
        let (ligand, _errors) =
            pdbtbx::open(&ligand_filename, pdbtbx::StrictnessLevel::Strict).unwrap();

        let scoring = PYDOCK::new(
            receptor,
            Vec::new(),
            Vec::new(),
            Vec::new(),
            0,
            ligand,
            Vec::new(),
            Vec::new(),
            Vec::new(),
            0,
            false,
        );

        let translation = vec![0., 0., 0.];
        let rotation = Quaternion::default();
        let energy = scoring.energy(&translation, &rotation, &Vec::new(), &Vec::new());
        assert!((energy - (-364.88126358158974)).abs() < 1e-8,
            "energy={energy} expected≈-364.88126358158974");
    }
}
