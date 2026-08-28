use super::constants::{INTERFACE_CUTOFF2, MEMBRANE_PENALTY_SCORE};
use super::qt::{rot3_apply, Quaternion};
use super::scoring::{membrane_intersection, satisfied_restraints, Score};
use npyz::NpyFile;
use pdbtbx::PDB;
use std::cell::RefCell;
use std::collections::HashMap;

// DFIRE2: 167 atom types, 30 distance bins (0.5 Å per bin, max 15 Å)
// bin = (distance * 2.0) as usize   (must be < 30)
// energy = Σ potential[a*167*30 + b*30 + bin] / 100.0
// Ref: Yang & Zhou, Protein Science 17:1212-1219 (2008)

const ATOM_TYPES: usize = 167;
const DIST_BINS: usize = 30;
const MAX_DIST_SQ: f64 = 15.0 * 15.0; // bin 30 → dist 15 Å upper bound

macro_rules! hashmap {
    ($( $key: expr => $val: expr ),*) => {{
         let mut map = ::std::collections::HashMap::new();
         $( map.insert($key, $val); )*
         map
    }}
}

lazy_static! {
    // Potential loaded once from embedded .npy binary
    static ref DFIRE2_POTENTIAL: Vec<f64> = {
        let raw = include_bytes!("../data/dfire2_energies.npy");
        NpyFile::new(&raw[..])
            .expect("Cannot parse dfire2_energies.npy")
            .into_vec::<f64>()
            .expect("Cannot read DFIRE2 potentials")
    };

    // Maps "RESNAME ATOMNAME" (with space) → DFIRE2 atom type index (0–166)
    // Ported verbatim from lightdock-0.9.4/lightdock/scoring/dfire2/driver.py
    static ref ATOMNUMBER: HashMap<&'static str, usize> = hashmap![
        "ALA N" => 0,   "ALA CA" => 1,  "ALA C" => 2,   "ALA O" => 3,   "ALA CB" => 4,
        "CYS N" => 5,   "CYS CA" => 6,  "CYS C" => 7,   "CYS O" => 8,   "CYS CB" => 9,  "CYS SG" => 10,
        "ASP N" => 11,  "ASP CA" => 12, "ASP C" => 13,  "ASP O" => 14,  "ASP CB" => 15,
        "ASP CG" => 16, "ASP OD1" => 17,"ASP OD2" => 18,
        "GLU N" => 19,  "GLU CA" => 20, "GLU C" => 21,  "GLU O" => 22,  "GLU CB" => 23,
        "GLU CG" => 24, "GLU CD" => 25, "GLU OE1" => 26,"GLU OE2" => 27,
        "PHE N" => 28,  "PHE CA" => 29, "PHE C" => 30,  "PHE O" => 31,  "PHE CB" => 32,
        "PHE CG" => 33, "PHE CD1" => 34,"PHE CD2" => 35,"PHE CE1" => 36,"PHE CE2" => 37,
        "PHE CZ" => 38,
        "GLY N" => 39,  "GLY CA" => 40, "GLY C" => 41,  "GLY O" => 42,
        "HIS N" => 43,  "HIS CA" => 44, "HIS C" => 45,  "HIS O" => 46,  "HIS CB" => 47,
        "HIS CG" => 48, "HIS ND1" => 49,"HIS CD2" => 50,"HIS CE1" => 51,"HIS NE2" => 52,
        "ILE N" => 53,  "ILE CA" => 54, "ILE C" => 55,  "ILE O" => 56,  "ILE CB" => 57,
        "ILE CG1" => 58,"ILE CG2" => 59,"ILE CD1" => 60,
        "LYS N" => 61,  "LYS CA" => 62, "LYS C" => 63,  "LYS O" => 64,  "LYS CB" => 65,
        "LYS CG" => 66, "LYS CD" => 67, "LYS CE" => 68, "LYS NZ" => 69,
        "LEU N" => 70,  "LEU CA" => 71, "LEU C" => 72,  "LEU O" => 73,  "LEU CB" => 74,
        "LEU CG" => 75, "LEU CD1" => 76,"LEU CD2" => 77,
        "MET N" => 78,  "MET CA" => 79, "MET C" => 80,  "MET O" => 81,  "MET CB" => 82,
        "MET CG" => 83, "MET SD" => 84, "MET CE" => 85,
        "ASN N" => 86,  "ASN CA" => 87, "ASN C" => 88,  "ASN O" => 89,  "ASN CB" => 90,
        "ASN CG" => 91, "ASN OD1" => 92,"ASN ND2" => 93,
        "PRO N" => 94,  "PRO CA" => 95, "PRO C" => 96,  "PRO O" => 97,  "PRO CB" => 98,
        "PRO CG" => 99, "PRO CD" => 100,
        "GLN N" => 101, "GLN CA" => 102,"GLN C" => 103, "GLN O" => 104, "GLN CB" => 105,
        "GLN CG" => 106,"GLN CD" => 107,"GLN OE1" => 108,"GLN NE2" => 109,
        "ARG N" => 110, "ARG CA" => 111,"ARG C" => 112, "ARG O" => 113, "ARG CB" => 114,
        "ARG CG" => 115,"ARG CD" => 116,"ARG NE" => 117,"ARG CZ" => 118,
        "ARG NH1" => 119,"ARG NH2" => 120,
        "SER N" => 121, "SER CA" => 122,"SER C" => 123, "SER O" => 124, "SER CB" => 125,
        "SER OG" => 126,
        "THR N" => 127, "THR CA" => 128,"THR C" => 129, "THR O" => 130, "THR CB" => 131,
        "THR OG1" => 132,"THR CG2" => 133,
        "VAL N" => 134, "VAL CA" => 135,"VAL C" => 136, "VAL O" => 137, "VAL CB" => 138,
        "VAL CG1" => 139,"VAL CG2" => 140,
        "TRP N" => 141, "TRP CA" => 142,"TRP C" => 143, "TRP O" => 144, "TRP CB" => 145,
        "TRP CG" => 146,"TRP CD1" => 147,"TRP CD2" => 148,"TRP NE1" => 149,"TRP CE2" => 150,
        "TRP CE3" => 151,"TRP CZ2" => 152,"TRP CZ3" => 153,"TRP CH2" => 154,
        "TYR N" => 155, "TYR CA" => 156,"TYR C" => 157, "TYR O" => 158, "TYR CB" => 159,
        "TYR CG" => 160,"TYR CD1" => 161,"TYR CD2" => 162,"TYR CE1" => 163,"TYR CE2" => 164,
        "TYR CZ" => 165,"TYR OH" => 166
    ];
}

pub struct DFIRE2DockingModel {
    pub atom_indices: Vec<usize>,
    pub residue_numbers: Vec<i32>,
    pub coordinates: Vec<[f64; 3]>,
    pub membrane: Vec<usize>,
    pub active_restraints: HashMap<String, Vec<usize>>,
    pub passive_restraints: HashMap<String, Vec<usize>>,
    pub num_anm: usize,
    pub nmodes: Vec<f64>,
}

impl DFIRE2DockingModel {
    fn new(
        structure: &PDB,
        active_restraints: &[String],
        passive_restraints: &[String],
        nmodes: &[f64],
        num_anm: usize,
    ) -> DFIRE2DockingModel {
        let mut model = DFIRE2DockingModel {
            atom_indices: Vec::new(),
            residue_numbers: Vec::new(),
            coordinates: Vec::new(),
            membrane: Vec::new(),
            active_restraints: HashMap::new(),
            passive_restraints: HashMap::new(),
            nmodes: nmodes.to_owned(),
            num_anm,
        };

        let mut atom_index: usize = 0;
        for chain in structure.chains() {
            for residue in chain.residues() {
                let res_name = match residue.name() {
                    Some(name) => name,
                    None => panic!("DFIRE2: Residue name error"),
                };
                let mut res_id = format!("{}.{}.{}", chain.id(), res_name, residue.serial_number());
                if let Some(c) = residue.insertion_code() {
                    res_id.push_str(c);
                }

                for atom in residue.atoms() {
                    let atom_key = format!("{} {}", res_name, atom.name());

                    // Skip atoms not in DFIRE2 dict (H, OXT, etc.)
                    let dfire2_idx = match ATOMNUMBER.get(&atom_key[..]) {
                        Some(&idx) => idx,
                        None => continue,
                    };

                    // Membrane bead
                    if res_name == "MMB" && atom.name() == "BJ" {
                        model.membrane.push(atom_index);
                    }

                    if active_restraints.contains(&res_id) {
                        model.active_restraints
                            .entry(res_id.clone())
                            .or_insert_with(Vec::new)
                            .push(atom_index);
                    }
                    if passive_restraints.contains(&res_id) {
                        model.passive_restraints
                            .entry(res_id.clone())
                            .or_insert_with(Vec::new)
                            .push(atom_index);
                    }

                    model.atom_indices.push(dfire2_idx);
                    model.residue_numbers.push(residue.serial_number() as i32);
                    model.coordinates.push([atom.x(), atom.y(), atom.z()]);
                    atom_index += 1;
                }
            }
        }
        model
    }
}

pub struct DFIRE2 {
    pub receptor: DFIRE2DockingModel,
    pub ligand: DFIRE2DockingModel,
    pub use_anm: bool,
}

impl DFIRE2 {
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
    ) -> Box<dyn Score> {
        let _ = &*DFIRE2_POTENTIAL; // ensure pre-loaded
        let d = DFIRE2 {
            receptor: DFIRE2DockingModel::new(
                &receptor,
                &rec_active_restraints,
                &rec_passive_restraints,
                &rec_nmodes,
                rec_num_anm,
            ),
            ligand: DFIRE2DockingModel::new(
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

impl Score for DFIRE2 {
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
        let potential = &*DFIRE2_POTENTIAL;
        let rot_mat = rotation.to_matrix();

        SCRATCH.with(|sc| {
            let mut sc = sc.borrow_mut();
            let (rec_c, lig_c, iface_r, iface_l) = &mut *sc;
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

            let rec_nm_n = if self.receptor.num_anm > 0 {
                self.receptor.nmodes.len() / (3 * self.receptor.num_anm)
            } else { rec_n };
            let lig_nm_n = if self.ligand.num_anm > 0 {
                self.ligand.nmodes.len() / (3 * self.ligand.num_anm)
            } else { lig_n };

            if self.use_anm && self.receptor.num_anm > 0 {
                for (i, c) in rec_c.iter_mut().enumerate() {
                    if i >= rec_nm_n { break; }
                    for k in 0..self.receptor.num_anm {
                        let b = k * rec_nm_n * 3 + i * 3;
                        c[0] += self.receptor.nmodes[b]   * rec_nmodes[k];
                        c[1] += self.receptor.nmodes[b+1] * rec_nmodes[k];
                        c[2] += self.receptor.nmodes[b+2] * rec_nmodes[k];
                    }
                }
            }
            for (i, c) in lig_c.iter_mut().enumerate() {
                let r = rot3_apply(&rot_mat, *c);
                c[0] = r[0] + translation[0];
                c[1] = r[1] + translation[1];
                c[2] = r[2] + translation[2];
                if self.use_anm && self.ligand.num_anm > 0 && i < lig_nm_n {
                    for k in 0..self.ligand.num_anm {
                        let b = k * lig_nm_n * 3 + i * 3;
                        c[0] += self.ligand.nmodes[b]   * lig_nmodes[k];
                        c[1] += self.ligand.nmodes[b+1] * lig_nmodes[k];
                        c[2] += self.ligand.nmodes[b+2] * lig_nmodes[k];
                    }
                }
            }

            // ── Phase 1: score on concatenated complex, full i<j double loop ──
            // Faithful port of lightdock cdfire2.c: receptor atoms followed by
            // ligand atoms in ONE array; every pair i<j is evaluated except
            // pairs from the same residue (res_indexes[i] == res_indexes[j]).
            // Ligand residue numbers are offset by the last receptor residue
            // number so cross-molecule same-numbered residues are NOT excluded.
            // Includes intra-molecular contributions, which are rigid-motion
            // invariant — this fixes the constant ~169 offset vs. the Python
            // reference.
            let rec_idx = &self.receptor.atom_indices;
            let lig_idx = &self.ligand.atom_indices;
            let rec_res = &self.receptor.residue_numbers;
            let lig_res_offset = match (rec_res.last(), self.ligand.residue_numbers.first()) {
                (Some(&last), Some(&first)) => last - first,
                _ => 0,
            };
            let mut res_all: Vec<i32> = Vec::with_capacity(rec_n + lig_n);
            res_all.extend_from_slice(rec_res);
            for r in self.ligand.residue_numbers.iter() {
                res_all.push(r + lig_res_offset);
            }

            let mut coords: Vec<[f64; 3]> = Vec::with_capacity(rec_n + lig_n);
            coords.extend_from_slice(rec_c.as_slice());
            coords.extend_from_slice(lig_c.as_slice());
            let idx_all: Vec<usize> = rec_idx.iter().chain(lig_idx.iter()).copied().collect();

            let n_atoms = coords.len();
            let mut score_raw = 0.0f64;
            for i in 0..n_atoms {
                let ri = &coords[i];
                let atom_a = idx_all[i];
                let res_i = res_all[i];
                for j in (i + 1)..n_atoms {
                    if res_i == res_all[j] {
                        continue;
                    }
                    let lj = &coords[j];
                    let dx = ri[0] - lj[0];
                    let dy = ri[1] - lj[1];
                    let dz = ri[2] - lj[2];
                    let dist2 = dx * dx + dy * dy + dz * dz;
                    if dist2 <= MAX_DIST_SQ {
                        let dist = dist2.sqrt();
                        let bin = (dist * 2.0) as usize;
                        if bin < DIST_BINS {
                            score_raw += potential
                                [atom_a * ATOM_TYPES * DIST_BINS + idx_all[j] * DIST_BINS + bin];
                        }
                    }
                }
            }

            let score = score_raw / 100.0;

            // ── Phase 2: interface flags (sequential, INTERFACE_CUTOFF=3.9Å) ──
            for (i, ra) in rec_c.iter().enumerate() {
                for (j, la) in lig_c.iter().enumerate() {
                    let dx = ra[0]-la[0]; let dy = ra[1]-la[1]; let dz = ra[2]-la[2];
                    if dx*dx + dy*dy + dz*dz <= INTERFACE_CUTOFF2 {
                        iface_r[i] = 1;
                        iface_l[j] = 1;
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
