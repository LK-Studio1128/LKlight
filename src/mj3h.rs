use super::constants::{INTERFACE_CUTOFF2, MEMBRANE_PENALTY_SCORE};
use super::qt::{rot3_apply, Quaternion};
use super::scoring::{membrane_intersection, satisfied_restraints, Score};
use pdbtbx::PDB;
use std::cell::RefCell;
use std::collections::HashMap;

lazy_static::lazy_static! {
    static ref MJ3H_POTENTIALS: Vec<Vec<f64>> = {
        let raw = include_bytes!("../data/MJ_potentials.dat");
        let contents = std::str::from_utf8(raw).expect("MJ_potentials.dat UTF-8 error");
        let mut matrix = vec![vec![0.0f64; NUM_RESIDUE_TYPES]; NUM_RESIDUE_TYPES];
        let mut in_mj3h = false;
        let mut row = 0usize;
        for line in contents.lines() {
            let trimmed = line.trim();
            if trimmed.starts_with(">MJ3h") { in_mj3h = true; row = 0; continue; }
            if trimmed.starts_with('<') { in_mj3h = false; continue; }
            if in_mj3h && row < NUM_RESIDUE_TYPES && !trimmed.is_empty() {
                let vals: Vec<f64> = trimmed.split_whitespace()
                    .filter_map(|s| s.parse::<f64>().ok()).collect();
                if vals.len() == NUM_RESIDUE_TYPES { matrix[row] = vals; row += 1; }
            }
        }
        assert_eq!(row, NUM_RESIDUE_TYPES, "MJ3h matrix incomplete");
        matrix
    };
}

// MJ3h: Miyazawa-Jernigan residue-level contact potential
// Uses side-chain centroid (all heavy atoms except backbone N, CA, C, O)
// Distance cutoffs (squared): min 6.2 Å² (2.49 Å), max 42.25 Å² (6.5 Å)
// Penalization = 3.0 for clashing contacts (dist² < 6.2)
// energy *= -1.0  (more negative = better)
// Ref: Miyazawa & Jernigan, J Mol Biol 1996;256:623–644

const NUM_RESIDUE_TYPES: usize = 20;
const MAX_DIST_SQ: f64 = 42.25;  // 6.5^2
const MIN_DIST_SQ: f64 = 6.2;    // ~2.49^2
const PENALIZATION: f64 = 3.0;

// Residue name → row/column index in the 20×20 MJ3h matrix
// Order matches the Python MJPotential.residues dict
fn res_to_idx(name: &str) -> Option<usize> {
    match name {
        "LEU" => Some(0),  "PHE" => Some(1),  "ILE" => Some(2),  "MET" => Some(3),
        "VAL" => Some(4),  "TRP" => Some(5),  "CYS" => Some(6),  "TYR" => Some(7),
        "HIS" => Some(8),  "ALA" => Some(9),  "THR" => Some(10), "GLY" => Some(11),
        "PRO" => Some(12), "ARG" => Some(13), "GLN" => Some(14), "SER" => Some(15),
        "ASN" => Some(16), "GLU" => Some(17), "ASP" => Some(18), "LYS" => Some(19),
        _ => None,
    }
}

// Atoms EXCLUDED from side-chain centroid.
// Faithful port of the Python reference (mj3h/driver.py):
//   not_considered_atoms = ["O", "C", "N", "H"]  (exact name match)
// Note CA *participates* in the centroid; only these four names are skipped.
fn is_excluded(name: &str) -> bool {
    matches!(name, "O" | "C" | "N" | "H")
}

pub struct MJ3hDockingModel {
    pub residue_types: Vec<usize>,
    // Per residue: list of sidechain heavy-atom coordinates (for on-the-fly centroid)
    pub residue_atoms: Vec<Vec<[f64; 3]>>,
    pub active_restraints: HashMap<String, Vec<usize>>,
    pub passive_restraints: HashMap<String, Vec<usize>>,
    pub membrane: Vec<usize>,
}

impl MJ3hDockingModel {
    fn new(
        structure: &PDB,
        active_restraints: &[String],
        passive_restraints: &[String],
    ) -> MJ3hDockingModel {
        let mut model = MJ3hDockingModel {
            residue_types: Vec::new(),
            residue_atoms: Vec::new(),
            active_restraints: HashMap::new(),
            passive_restraints: HashMap::new(),
            membrane: Vec::new(),
        };

        let mut res_index: usize = 0;
        for chain in structure.chains() {
            for residue in chain.residues() {
                let res_name = match residue.name() {
                    Some(name) => name,
                    None => continue,
                };
                let res_idx = match res_to_idx(res_name) {
                    Some(idx) => idx,
                    None => continue,  // skip non-standard residues
                };

                // Collect centroid-contributing atoms (exact-name exclusion of O/C/N/H)
                let contributing: Vec<[f64; 3]> = residue
                    .atoms()
                    .filter(|a| !is_excluded(a.name()))
                    .map(|a| [a.x(), a.y(), a.z()])
                    .collect();

                // Python: if count == 0 the residue contributes nothing and is removed
                let atoms = contributing;

                if atoms.is_empty() {
                    continue;  // skip if no usable atoms at all
                }

                let mut res_id = format!("{}.{}.{}", chain.id(), res_name, residue.serial_number());
                if let Some(c) = residue.insertion_code() {
                    res_id.push_str(c);
                }

                if active_restraints.contains(&res_id) {
                    model.active_restraints
                        .insert(res_id.clone(), vec![res_index]);
                }
                if passive_restraints.contains(&res_id) {
                    model.passive_restraints
                        .insert(res_id.clone(), vec![res_index]);
                }

                model.residue_types.push(res_idx);
                model.residue_atoms.push(atoms);
                res_index += 1;
            }
        }
        model
    }

    // Compute centroid for residue i (no transformation — use for receptor)
    pub fn centroid(&self, i: usize) -> [f64; 3] {
        let atoms = &self.residue_atoms[i];
        let n = atoms.len() as f64;
        let cx = atoms.iter().map(|a| a[0]).sum::<f64>() / n;
        let cy = atoms.iter().map(|a| a[1]).sum::<f64>() / n;
        let cz = atoms.iter().map(|a| a[2]).sum::<f64>() / n;
        [cx, cy, cz]
    }

    // Compute centroid for residue i after applying rotation + translation
    pub fn centroid_transformed(&self, i: usize, rot_mat: &[[f64;3];3], translation: &[f64]) -> [f64; 3] {
        let atoms = &self.residue_atoms[i];
        let n = atoms.len() as f64;
        let mut cx = 0.0_f64;
        let mut cy = 0.0_f64;
        let mut cz = 0.0_f64;
        for atom in atoms {
            let r = rot3_apply(rot_mat, *atom);
            cx += r[0] + translation[0];
            cy += r[1] + translation[1];
            cz += r[2] + translation[2];
        }
        [cx / n, cy / n, cz / n]
    }
}

pub struct MJ3h {
    pub receptor: MJ3hDockingModel,
    pub ligand: MJ3hDockingModel,
}

impl MJ3h {
    pub fn new(
        receptor: PDB,
        rec_active_restraints: Vec<String>,
        rec_passive_restraints: Vec<String>,
        ligand: PDB,
        lig_active_restraints: Vec<String>,
        lig_passive_restraints: Vec<String>,
    ) -> Box<dyn Score> {
        let _ = &*MJ3H_POTENTIALS; // ensure pre-loaded
        Box::new(MJ3h {
            receptor: MJ3hDockingModel::new(&receptor, &rec_active_restraints, &rec_passive_restraints),
            ligand:   MJ3hDockingModel::new(&ligand,   &lig_active_restraints, &lig_passive_restraints),
        })
    }
}

impl Score for MJ3h {
    fn energy(
        &self,
        translation: &[f64],
        rotation: &Quaternion,
        _rec_nmodes: &[f64],
        _lig_nmodes: &[f64],
    ) -> f64 {
        thread_local! {
            static SCRATCH: RefCell<(Vec<[f64;3]>, Vec<[f64;3]>, Vec<usize>, Vec<usize>)> =
                RefCell::new((Vec::new(), Vec::new(), Vec::new(), Vec::new()));
        }
        let pot = &*MJ3H_POTENTIALS;
        let rot_mat = rotation.to_matrix();

        SCRATCH.with(|sc| {
        let mut sc = sc.borrow_mut();
        let (rec_c, lig_c, iface_r, iface_l) = &mut *sc;
        let n_rec = self.receptor.residue_types.len();
        let n_lig = self.ligand.residue_types.len();
        if rec_c.len() != n_rec { rec_c.resize(n_rec, [0.0;3]); }
        if lig_c.len() != n_lig { lig_c.resize(n_lig, [0.0;3]); }
        if iface_r.len() != n_rec { iface_r.resize(n_rec, 0); }
        if iface_l.len() != n_lig { iface_l.resize(n_lig, 0); }
        for (i, c) in rec_c.iter_mut().enumerate() { *c = self.receptor.centroid(i); }
        for (j, c) in lig_c.iter_mut().enumerate() { *c = self.ligand.centroid_transformed(j, &rot_mat, translation); }
        for v in iface_r.iter_mut() { *v = 0; }
        for v in iface_l.iter_mut() { *v = 0; }

        let mut energy = 0.0f64;
        for (i, rc) in rec_c.iter().enumerate() {
            let ri = self.receptor.residue_types[i];
            for (j, lc) in lig_c.iter().enumerate() {
                let dx = rc[0]-lc[0]; let dy = rc[1]-lc[1]; let dz = rc[2]-lc[2];
                let dist_sq = dx*dx + dy*dy + dz*dz;
                if dist_sq < MAX_DIST_SQ {
                    if dist_sq < MIN_DIST_SQ { energy += PENALIZATION; }
                    else { let rj = self.ligand.residue_types[j]; energy += pot[ri][rj]; }
                    if dist_sq <= INTERFACE_CUTOFF2 { iface_r[i] = 1; iface_l[j] = 1; }
                }
            }
        }
        energy *= -1.0;
        let perc_r = satisfied_restraints(iface_r, &self.receptor.active_restraints);
        let perc_l = satisfied_restraints(iface_l, &self.ligand.active_restraints);
        let intersection = membrane_intersection(iface_r, &self.receptor.membrane);
        let penalty = if intersection > 0.0 { MEMBRANE_PENALTY_SCORE * intersection } else { 0.0 };
        energy + perc_r * energy + perc_l * energy - penalty
        })
    }
}
