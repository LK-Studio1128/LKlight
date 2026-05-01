use super::constants::{INTERFACE_CUTOFF2, MEMBRANE_PENALTY_SCORE};
use super::qt::{rot3_apply, Quaternion};
use super::scoring::{membrane_intersection, satisfied_restraints, Score};
use pdbtbx::PDB;
use std::cell::RefCell;
use std::collections::HashMap;

macro_rules! hashmap {
    ($( $key: expr => $val: expr ),*) => {{
         let mut map = ::std::collections::HashMap::new();
         $( map.insert($key, $val); )*
         map
    }}
}

pub fn r3_to_numerical(residue_name: &str) -> usize {
    match residue_name {
        "ALA" => 0,
        "CYS" => 1,
        "ASP" => 2,
        "GLU" => 3,
        "PHE" => 4,
        "GLY" => 5,
        "HIS" => 6,
        "ILE" => 7,
        "LYS" => 8,
        "LEU" => 9,
        "MET" => 10,
        "ASN" => 11,
        "PRO" => 12,
        "GLN" => 13,
        "ARG" => 14,
        "SER" => 15,
        "THR" => 16,
        "VAL" => 17,
        "TRP" => 18,
        "TYR" => 19,
        "MMB" => 20,
        "MMY" => 0,
        _ => return 999  // unsupported residue — caller skips atom
    }
}

// DFIRE only uses 20 distance bins
const DIST_TO_BINS: &[usize] = &[
    1, 1, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 14, 15, 15, 16, 16, 17, 17, 18, 18, 19,
    19, 20, 20, 21, 21, 22, 22, 23, 23, 24, 24, 25, 25, 26, 26, 27, 27, 28, 28, 29, 29, 30, 30, 31,
    32,
];

lazy_static! {
    // Potential table loaded once at startup from embedded DCparams binary
    static ref DFIRE_POTENTIAL: Vec<f64> = {
        let raw = include_bytes!("../data/DCparams");
        std::str::from_utf8(raw)
            .expect("DCparams is not valid UTF-8")
            .split_whitespace()
            .take(169 * 169 * 20)
            .map(|s| s.parse::<f64>().expect("DCparams parse error"))
            .collect()
    };

    static ref ATOMNUMBER: HashMap<&'static str, usize> = hashmap![
        "ALAN" => 0, "ALACA" => 1, "ALAC" => 2, "ALAO" => 3, "ALACB" => 4,
        "CYSN" => 0, "CYSCA" => 1, "CYSC" => 2, "CYSO" => 3, "CYSCB" => 4, "CYSSG" => 5,
        "ASPN" => 0, "ASPCA" => 1, "ASPC" => 2, "ASPO" => 3, "ASPCB" => 4, "ASPCG" => 5, "ASPOD1" => 6, "ASPOD2" => 7,
        "GLUN" => 0, "GLUCA" => 1, "GLUC" => 2, "GLUO" => 3, "GLUCB" => 4, "GLUCG" => 5, "GLUCD" => 6, "GLUOE1" => 7, "GLUOE2" => 8,
        "PHEN" => 0, "PHECA" => 1, "PHEC" => 2, "PHEO" => 3, "PHECB" => 4, "PHECG" => 5, "PHECD1" => 6, "PHECD2" => 7, "PHECE1" => 8, "PHECE2" => 9, "PHECZ" => 10,
        "GLYN" => 0, "GLYCA" => 1, "GLYC" => 2, "GLYO" => 3,
        "HISN" => 0, "HISCA" => 1, "HISC" => 2, "HISO" => 3, "HISCB" => 4, "HISCG" => 5, "HISND1" => 6, "HISCD2" => 7, "HISCE1" => 8, "HISNE2" => 9,
        "ILEN" => 0, "ILECA" => 1, "ILEC" => 2, "ILEO" => 3, "ILECB" => 4, "ILECG1" => 5, "ILECG2" => 6, "ILECD1" => 7,
        "LYSN" => 0, "LYSCA" => 1, "LYSC" => 2, "LYSO" => 3, "LYSCB" => 4, "LYSCG" => 5, "LYSCD" => 6, "LYSCE" => 7, "LYSNZ" => 8,
        "LEUN" => 0, "LEUCA" => 1, "LEUC" => 2, "LEUO" => 3, "LEUCB" => 4, "LEUCG" => 5, "LEUCD1" => 6, "LEUCD2" => 7,
        "METN" => 0, "METCA" => 1, "METC" => 2, "METO" => 3, "METCB" => 4, "METCG" => 5, "METSD" => 6, "METCE" => 7,
        "ASNN" => 0, "ASNCA" => 1, "ASNC" => 2, "ASNO" => 3, "ASNCB" => 4, "ASNCG" => 5, "ASNOD1" => 6, "ASNND2" => 7,
        "PRON" => 0, "PROCA" => 1, "PROC" => 2, "PROO" => 3, "PROCB" => 4, "PROCG" => 5, "PROCD" => 6,
        "GLNN" => 0, "GLNCA" => 1, "GLNC" => 2, "GLNO" => 3, "GLNCB" => 4, "GLNCG" => 5, "GLNCD" => 6, "GLNOE1" => 7, "GLNNE2" => 8,
        "ARGN" => 0, "ARGCA" => 1, "ARGC" => 2, "ARGO" => 3, "ARGCB" => 4, "ARGCG" => 5, "ARGCD" => 6, "ARGNE" => 7, "ARGCZ" => 8, "ARGNH1" => 9, "ARGNH2" => 10,
        "SERN" => 0, "SERCA" => 1, "SERC" => 2, "SERO" => 3, "SERCB" => 4, "SEROG" => 5,
        "THRN" => 0, "THRCA" => 1, "THRC" => 2, "THRO" => 3, "THRCB" => 4, "THROG1" => 5, "THRCG2" => 6,
        "VALN" => 0, "VALCA" => 1, "VALC" => 2, "VALO" => 3, "VALCB" => 4, "VALCG1" => 5, "VALCG2" => 6,
        "TRPN" => 0, "TRPCA" => 1, "TRPC" => 2, "TRPO" => 3, "TRPCB" => 4, "TRPCG" => 5, "TRPCD1" => 6, "TRPCD2" => 7, "TRPCE2" => 8, "TRPNE1" => 9, "TRPCE3" => 10, "TRPCZ3" => 11, "TRPCH2" => 12, "TRPCZ2" => 13,
        "TYRN" => 0, "TYRCA" => 1, "TYRC" => 2, "TYRO" => 3, "TYRCB" => 4, "TYRCG" => 5, "TYRCD1" => 6, "TYRCD2" => 7, "TYRCE1" => 8, "TYRCE2" => 9, "TYRCZ" => 10, "TYROH" => 11,
        "MMBBJ" => 0, "MMYDU" => 0];

    // Atom type and residue translation matrix
    static ref ATOMRES: Vec<Vec<usize>> = vec![vec![74, 75, 76, 77, 78, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                                               vec![0, 1, 2, 3, 4, 5, 0, 0, 0, 0, 0, 0, 0, 0],
                                               vec![122, 123, 124, 125, 126, 127, 128, 129, 0, 0, 0, 0, 0, 0],
                                               vec![113, 114, 115, 116, 117, 118, 119, 120, 121, 0, 0, 0, 0, 0],
                                               vec![14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 0, 0, 0],
                                               vec![79, 80, 81, 82, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                                               vec![130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 0, 0, 0, 0],
                                               vec![25, 26, 27, 28, 29, 30, 31, 32, 0, 0, 0, 0, 0, 0],
                                               vec![151, 152, 153, 154, 155, 156, 157, 158, 159, 0, 0, 0, 0, 0],
                                               vec![33, 34, 35, 36, 37, 38, 39, 40, 0, 0, 0, 0, 0, 0],
                                               vec![6, 7, 8, 9, 10, 11, 12, 13, 0, 0, 0, 0, 0, 0],
                                               vec![105, 106, 107, 108, 109, 110, 111, 112, 0, 0, 0, 0, 0, 0],
                                               vec![160, 161, 162, 163, 164, 165, 166, 0, 0, 0, 0, 0, 0, 0],
                                               vec![96, 97, 98, 99, 100, 101, 102, 103, 104, 0, 0, 0, 0, 0],
                                               vec![140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 0, 0, 0],
                                               vec![90, 91, 92, 93, 94, 95, 0, 0, 0, 0, 0, 0, 0, 0],
                                               vec![83, 84, 85, 86, 87, 88, 89, 0, 0, 0, 0, 0, 0, 0],
                                               vec![41, 42, 43, 44, 45, 46, 47, 0, 0, 0, 0, 0, 0, 0],
                                               vec![48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61],
                                               vec![62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 0, 0],
                                               vec![167, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                                               vec![74, 75, 76, 77, 78, 0, 0, 0, 0, 0, 0, 0, 0, 0]];
}

pub struct DFIREDockingModel {
    pub atoms: Vec<usize>,
    pub coordinates: Vec<[f64; 3]>,
    pub membrane: Vec<usize>,
    pub active_restraints: HashMap<String, Vec<usize>>,
    pub passive_restraints: HashMap<String, Vec<usize>>,
    pub num_anm: usize,
    pub nmodes: Vec<f64>,
}

impl<'a> DFIREDockingModel {
    fn new(
        structure: &'a PDB,
        active_restraints: &'a [String],
        passive_restraints: &'a [String],
        nmodes: &[f64],
        num_anm: usize,
    ) -> DFIREDockingModel {
        let mut model = DFIREDockingModel {
            atoms: Vec::new(),
            coordinates: Vec::new(),
            membrane: Vec::new(),
            active_restraints: HashMap::new(),
            passive_restraints: HashMap::new(),
            nmodes: nmodes.to_owned(),
            num_anm,
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

                    let rnuma = r3_to_numerical(res_name);
                    if rnuma == 999 { continue; }  // unsupported residue — skip
                    let anuma = match ATOMNUMBER.get(&rec_atom_type[..]) {
                        Some(&a) => a,
                        _ => continue,  // unsupported atom type — skip
                    };
                    let atoma = ATOMRES[rnuma][anuma];
                    model.atoms.push(atoma);
                    model.coordinates.push([atom.x(), atom.y(), atom.z()]);
                    atom_index += 1;
                }
            }
        }
        model
    }
}

pub struct DFIRE {
    pub receptor: DFIREDockingModel,
    pub ligand: DFIREDockingModel,
    pub use_anm: bool,
}

impl<'a> DFIRE {
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
        let _ = &*DFIRE_POTENTIAL; // ensure potential is pre-loaded
        let d = DFIRE {
            receptor: DFIREDockingModel::new(
                &receptor,
                &rec_active_restraints,
                &rec_passive_restraints,
                &rec_nmodes,
                rec_num_anm,
            ),
            ligand: DFIREDockingModel::new(
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

    pub fn get_potential(&self, x: usize, y: usize, z: usize) -> f64 {
        DFIRE_POTENTIAL[x + 169 * (y + 20 * z)]
    }
}

impl Score for DFIRE {
    fn energy(
        &self,
        translation: &[f64],
        rotation: &Quaternion,
        rec_nmodes: &[f64],
        lig_nmodes: &[f64],
    ) -> f64 {
        // Thread-local scratch: avoids heap allocation on every call
        thread_local! {
            static SCRATCH: RefCell<(
                Vec<[f64; 3]>, // receptor coords
                Vec<[f64; 3]>, // ligand coords
                Vec<usize>,    // interface_receptor flags
                Vec<usize>,    // interface_ligand flags
            )> = RefCell::new((Vec::new(), Vec::new(), Vec::new(), Vec::new()));
        }

        let potential = &*DFIRE_POTENTIAL;
        let rot_mat = rotation.to_matrix(); // precompute once, reuse per atom

        SCRATCH.with(|sc| {
            let mut sc = sc.borrow_mut();
            let (rec_c, lig_c, iface_r, iface_l) = &mut *sc;

            let rec_n = self.receptor.coordinates.len();
            let lig_n = self.ligand.coordinates.len();

            // Grow scratch buffers if needed; otherwise just overwrite in place
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

            // Apply ANM to receptor (no rotation/translation)
            if self.use_anm && self.receptor.num_anm > 0 {
                for (i_atom, coord) in rec_c.iter_mut().enumerate() {
                    if i_atom >= rec_nm_n { break; }
                    for i_nm in 0..self.receptor.num_anm {
                        let base = i_nm * rec_nm_n * 3 + i_atom * 3;
                        coord[0] += self.receptor.nmodes[base]     * rec_nmodes[i_nm];
                        coord[1] += self.receptor.nmodes[base + 1] * rec_nmodes[i_nm];
                        coord[2] += self.receptor.nmodes[base + 2] * rec_nmodes[i_nm];
                    }
                }
            }

            // Rotate + translate + ANM for ligand
            for (i_atom, coord) in lig_c.iter_mut().enumerate() {
                let r = rot3_apply(&rot_mat, *coord);
                coord[0] = r[0] + translation[0];
                coord[1] = r[1] + translation[1];
                coord[2] = r[2] + translation[2];
                if self.use_anm && self.ligand.num_anm > 0 && i_atom < lig_nm_n {
                    for i_nm in 0..self.ligand.num_anm {
                        let base = i_nm * lig_nm_n * 3 + i_atom * 3;
                        coord[0] += self.ligand.nmodes[base]     * lig_nmodes[i_nm];
                        coord[1] += self.ligand.nmodes[base + 1] * lig_nmodes[i_nm];
                        coord[2] += self.ligand.nmodes[base + 2] * lig_nmodes[i_nm];
                    }
                }
            }

            // ── Phase 1: parallel score (receptor atoms outer, ligand inner) ──
            let rec_atoms = &self.receptor.atoms;
            let lig_atoms = &self.ligand.atoms;
            let lig_slice: &[[f64; 3]] = lig_c.as_slice();
            let lig_n_atoms = lig_slice.len();

            let score_raw: f64 = rec_c.iter().enumerate()
                .map(|(i, ra)| {
                    let rx = ra[0]; let ry = ra[1]; let rz = ra[2];
                    let atoma = rec_atoms[i];
                    let mut s = 0.0f64;
                    for j in 0..lig_n_atoms {
                        let la = &lig_slice[j];
                        let dx = rx - la[0]; let dy = ry - la[1]; let dz = rz - la[2];
                        let dist2 = dx*dx + dy*dy + dz*dz;
                        if dist2 <= 225.0 {
                            let atomb = lig_atoms[j];
                            let d = dist2.sqrt() * 2.0 - 1.0;
                            let dfire_bin = DIST_TO_BINS[d as usize] - 1;
                            s += potential[atoma * 169 * 20 + atomb * 20 + dfire_bin];
                        }
                    }
                    s
                })
                .sum();

            let score = (score_raw * 0.0157 - 4.7) * -1.0;

            // ── Phase 2: interface flags (sequential, INTERFACE_CUTOFF=3.9Å) ──
            for (i, ra) in rec_c.iter().enumerate() {
                for (j, la) in lig_slice.iter().enumerate() {
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::qt::Quaternion;
    use std::env;

    // #[test]
    // fn test_read_potentials() {
    //     let mut scoring = DFIRE {
    //         potential: Vec::with_capacity(168 * 168 * 20),
    //     };
    //     scoring.load_potentials();
    //     assert_eq!(scoring.potential[0], 10.0);
    //     assert_eq!(scoring.potential[2], -0.624030868);
    //     assert_eq!(scoring.potential[4998], -0.0458685914);
    //     assert_eq!(scoring.potential[168*168*20-1], 0.0);
    // }

    #[test]
    fn test_2oob() {
        let cargo_path = match env::var("CARGO_MANIFEST_DIR") {
            Ok(val) => val,
            Err(_) => String::from("."),
        };
        let test_path: String = format!("{}/tests/2oob", cargo_path);

        let receptor_filename: String = format!("{}/2oob_receptor.pdb", test_path);
        let (receptor, _errors) =
            pdbtbx::open(&receptor_filename, pdbtbx::StrictnessLevel::Strict).unwrap();

        let ligand_filename: String = format!("{}/2oob_ligand.pdb", test_path);
        let (ligand, _errors) =
            pdbtbx::open(&ligand_filename, pdbtbx::StrictnessLevel::Strict).unwrap();

        let scoring = DFIRE::new(
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
        assert!((energy - 16.7540569503498).abs() < 1e-9, "energy = {}", energy);
    }
}
