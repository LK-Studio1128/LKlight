use super::constants::INTERFACE_CUTOFF;
use super::qt::{rot3_apply, Quaternion};
use std::cell::RefCell;
use super::scoring::{satisfied_restraints, Score};
use pdbtbx::PDB;
use std::collections::HashMap;

// SIPPER: sequence-independent pair potential for protein-protein docking
// Ref: Viswanath et al., J Chem Inf Model 51, 2515-2527 (2011)
//
// Residue-level contact potential: one energy added per contacting residue pair
// Contact defined as any atom-atom distance² < 25.0 Å² (5 Å)
// Energy: total_energy = -1.0 * (total_sipper - 0.019 * total_oda)
// Interface cutoff: 3.9 Å (INTERFACE_CUTOFF²)
// ODA weighting is optional; defaults to 0.0 per residue when no .oda file

const DIST2_CUTOFF: f64 = 25.0;
const INTERFACE_CUTOFF2: f64 = INTERFACE_CUTOFF * INTERFACE_CUTOFF;
const ODA_WEIGHT: f64 = 0.019;

fn res_to_index(name: &str) -> Option<usize> {
    match name {
        "ALA" => Some(0),  "ARG" => Some(1),  "ASN" => Some(2),  "ASP" => Some(3),
        "CYS" => Some(4),  "GLN" => Some(5),  "GLU" => Some(6),  "GLY" => Some(7),
        "HIS" | "HIP" | "HID" => Some(8),
        "ILE" => Some(9),  "LEU" => Some(10), "LYS" => Some(11), "MET" => Some(12),
        "PHE" => Some(13), "PRO" => Some(14), "SER" => Some(15), "THR" => Some(16),
        "TRP" => Some(17), "TYR" => Some(18), "VAL" => Some(19),
        _ => None,
    }
}

// 20×20 SIPPER energy matrix (row-major, order matches res_to_index above)
#[rustfmt::skip]
const SIPPER_ENERGY: [[f64; 20]; 20] = [
    [ 0.0203,-0.0863,-0.4557,-0.2390,-0.1248,-0.1364,-0.3296,-0.1967, 0.0617,-0.1510, 0.3005,-0.5890, 0.3480, 0.3347,-0.3933,-0.0485,-0.5509,-0.1197,-0.0258, 0.2342],
    [-0.0863,-0.1238, 0.0758, 0.3463, 0.0516, 0.1506, 0.3716, 0.1986, 0.0733, 0.0189, 0.0798,-0.2429, 0.1691, 0.2749, 0.0556, 0.0307,-0.1641, 0.4066, 0.3020,-0.1588],
    [-0.4557, 0.0758, 0.1364,-0.0668,-0.7083, 0.1542,-0.0790,-0.3221,-0.0791,-0.3255,-0.2293,-0.2471,-0.0370,-0.1128, 0.0580,-0.2231, 0.1560, 0.1709, 0.2394,-0.2417],
    [-0.2390, 0.3463,-0.0668,-0.3435,-0.6855,-0.3219,-0.4705,-0.1295, 0.0531,-0.3258,-0.3965, 0.3560,-0.2127,-0.4031,-0.3132, 0.0224,-0.3750,-0.2714, 0.1343,-0.3264],
    [-0.1248, 0.0516,-0.7083,-0.6855,-0.8104, 0.2972,-0.1732,-0.0774, 0.2898,-0.0418, 0.2815, 0.2155, 0.1408, 0.1201,-1.4258, 0.2154,-0.7522, 0.4573, 0.0093,-0.2817],
    [-0.1364, 0.1506, 0.1542,-0.3219, 0.2972,-0.0158,-0.1981,-0.0469,-0.3700, 0.2238, 0.0074,-0.0123,-0.2583, 0.0009, 0.1196,-0.1580,-0.0776, 0.0481, 0.2801, 0.1851],
    [-0.3296, 0.3716,-0.0790,-0.4705,-0.1732,-0.1981,-0.3159,-0.4593,-0.0476,-0.0528,-0.3962, 0.1902, 0.0891,-0.2260,-0.2402,-0.1619,-0.0918,-0.1682,-0.0743,-0.2364],
    [-0.1967, 0.1986,-0.3221,-0.1295,-0.0774,-0.0469,-0.4593,-0.1134, 0.0774,-0.1675,-0.0427,-0.2603, 0.1048,-0.1419,-0.0669,-0.2623,-0.2991,-0.0405, 0.0623,-0.1274],
    [ 0.0617, 0.0733,-0.0791, 0.0531, 0.2898,-0.3700,-0.0476, 0.0774,-0.0065, 0.3230,-0.0958,-0.0730, 0.5929, 0.1252,-0.1376,-0.3582,-0.0609, 0.1252, 0.1903,-0.0374],
    [-0.1510, 0.0189,-0.3255,-0.3258,-0.0418, 0.2238,-0.0528,-0.1675, 0.3230, 0.1240, 0.4622,-0.2669, 0.0782, 0.4798,-0.0497,-0.0236,-0.1774, 0.2025, 0.1888, 0.2096],
    [ 0.3005, 0.0798,-0.2293,-0.3965, 0.2815, 0.0074,-0.3962,-0.0427,-0.0958, 0.4622, 0.3728, 0.0411, 0.5373, 0.4378, 0.1188,-0.1046,-0.2240, 0.1605,-0.0031, 0.4125],
    [-0.5890,-0.2429,-0.2471, 0.3560, 0.2155,-0.0123, 0.1902,-0.2603,-0.0730,-0.2669, 0.0411,-0.5909, 0.0855,-0.0428,-0.3767,-0.0221,-0.2437, 0.2437, 0.0405,-0.4158],
    [ 0.3480, 0.1691,-0.0370,-0.2127, 0.1408,-0.2583, 0.0891, 0.1048, 0.5929, 0.0782, 0.5373, 0.0855, 0.7834, 0.2749, 0.2739,-0.2994,-0.0021, 0.2643, 0.1597, 0.4304],
    [ 0.3347, 0.2749,-0.1128,-0.4031, 0.1201, 0.0009,-0.2260,-0.1419, 0.1252, 0.4798, 0.4378,-0.0428, 0.2749, 0.6080, 0.0835, 0.0307, 0.0467, 0.4828, 0.1903, 0.1886],
    [-0.3933, 0.0556, 0.0580,-0.3132,-1.4258, 0.1196,-0.2402,-0.0669,-0.1376,-0.0497, 0.1188,-0.3767, 0.2739, 0.0835,-0.1661, 0.1552,-0.0582, 0.5238, 0.3497,-0.2488],
    [-0.0485, 0.0307,-0.2231, 0.0224, 0.2154,-0.1580,-0.1619,-0.2623,-0.3582,-0.0236,-0.1046,-0.0221,-0.2994, 0.0307, 0.1552,-0.2341,-0.1813,-0.1190, 0.0475,-0.0925],
    [-0.5509,-0.1641, 0.1560,-0.3750,-0.7522,-0.0776,-0.0918,-0.2991,-0.0609,-0.1774,-0.2240,-0.2437,-0.0021, 0.0467,-0.0582,-0.1813,-0.1036, 0.0086, 0.0089, 0.1553],
    [-0.1197, 0.4066, 0.1709,-0.2714, 0.4573, 0.0481,-0.1682,-0.0405, 0.1252, 0.2025, 0.1605, 0.2437, 0.2643, 0.4828, 0.5238,-0.1190, 0.0086, 0.2568, 0.4536, 0.0098],
    [-0.0258, 0.3020, 0.2394, 0.1343, 0.0093, 0.2801,-0.0743, 0.0623, 0.1903, 0.1888,-0.0031, 0.0405, 0.1597, 0.1903, 0.3497, 0.0475, 0.0089, 0.4536,-0.0037, 0.1134],
    [ 0.2342,-0.1588,-0.2417,-0.3264,-0.2817, 0.1851,-0.2364,-0.1274,-0.0374, 0.2096, 0.4125,-0.4158, 0.4304, 0.1886,-0.2488,-0.0925, 0.1553, 0.0098, 0.1134, 0.1246],
];

pub struct SIPPERDockingModel {
    pub residue_types: Vec<usize>,
    pub residue_atom_ranges: Vec<(usize, usize)>,
    pub coordinates: Vec<[f64; 3]>,
    pub oda: Vec<f64>,
    pub active_restraints: HashMap<String, Vec<usize>>,
    pub passive_restraints: HashMap<String, Vec<usize>>,
}

impl SIPPERDockingModel {
    fn new(
        structure: &PDB,
        active_restraints: &[String],
        passive_restraints: &[String],
    ) -> SIPPERDockingModel {
        let mut model = SIPPERDockingModel {
            residue_types: Vec::new(),
            residue_atom_ranges: Vec::new(),
            coordinates: Vec::new(),
            oda: Vec::new(),
            active_restraints: HashMap::new(),
            passive_restraints: HashMap::new(),
        };

        for chain in structure.chains() {
            for residue in chain.residues() {
                let res_name = match residue.name() { Some(n) => n, None => continue };
                let res_idx = match res_to_index(res_name) {
                    Some(i) => i,
                    None => continue,
                };

                let atom_start = model.coordinates.len();
                for atom in residue.atoms() {
                    if atom.element().map(|e| e.symbol() == "H").unwrap_or(false) { continue; }
                    let aname = atom.name().trim();
                    if aname == "H" || aname.starts_with("H") { continue; }
                    model.coordinates.push([atom.x(), atom.y(), atom.z()]);
                }
                let atom_end = model.coordinates.len();
                if atom_end == atom_start { continue; }

                let res_index = model.residue_types.len();
                model.residue_types.push(res_idx);
                model.residue_atom_ranges.push((atom_start, atom_end));
                model.oda.push(0.0);

                let mut res_id = format!("{}.{}.{}", chain.id(), res_name, residue.serial_number());
                if let Some(c) = residue.insertion_code() { res_id.push_str(c); }
                let atom_indices: Vec<usize> = (atom_start..atom_end).collect();
                if active_restraints.contains(&res_id) {
                    model.active_restraints.insert(res_id.clone(), atom_indices.clone());
                }
                if passive_restraints.contains(&res_id) {
                    model.passive_restraints.insert(res_id.clone(), atom_indices);
                }
                let _ = res_index;
            }
        }
        model
    }
}

pub struct SIPPER {
    pub receptor: SIPPERDockingModel,
    pub ligand: SIPPERDockingModel,
    pub use_anm: bool,
}

impl<'a> SIPPER {
    pub fn new(
        receptor: PDB,
        rec_active_restraints: Vec<String>,
        rec_passive_restraints: Vec<String>,
        ligand: PDB,
        lig_active_restraints: Vec<String>,
        lig_passive_restraints: Vec<String>,
    ) -> Box<dyn Score + 'a> {
        Box::new(SIPPER {
            receptor: SIPPERDockingModel::new(&receptor, &rec_active_restraints, &rec_passive_restraints),
            ligand: SIPPERDockingModel::new(&ligand, &lig_active_restraints, &lig_passive_restraints),
            use_anm: false,
        })
    }
}

impl Score for SIPPER {
    fn energy(
        &self,
        translation: &[f64],
        rotation: &Quaternion,
        _rec_nmodes: &[f64],
        _lig_nmodes: &[f64],
    ) -> f64 {
        thread_local! {
            static SCRATCH: RefCell<(Vec<[f64;3]>, Vec<usize>, Vec<usize>)> =
                RefCell::new((Vec::new(), Vec::new(), Vec::new()));
        }
        let rot_mat = rotation.to_matrix();

        SCRATCH.with(|sc| {
        let mut sc = sc.borrow_mut();
        let (lig_c, iface_r, iface_l) = &mut *sc;
        let n_lig_atoms = self.ligand.coordinates.len();
        let n_rec_atoms = self.receptor.coordinates.len();
        if lig_c.len() != n_lig_atoms { lig_c.resize(n_lig_atoms, [0.0;3]); }
        if iface_r.len() != n_rec_atoms { iface_r.resize(n_rec_atoms, 0); }
        if iface_l.len() != n_lig_atoms { iface_l.resize(n_lig_atoms, 0); }
        for (c, src) in lig_c.iter_mut().zip(self.ligand.coordinates.iter()) {
            let r = rot3_apply(&rot_mat, *src);
            *c = [r[0]+translation[0], r[1]+translation[1], r[2]+translation[2]];
        }
        for v in iface_r.iter_mut() { *v = 0; }
        for v in iface_l.iter_mut() { *v = 0; }

        let rec_coords = &self.receptor.coordinates;
        let n_rec_res = self.receptor.residue_types.len();
        let n_lig_res = self.ligand.residue_types.len();
        let mut total_sipper = 0.0_f64;
        let mut total_oda = 0.0_f64;

        // The C extension reads the int64 `indexes` and `atoms_per_residue`
        // arrays through a uint32 pointer. For values < 2^32 the resulting
        // view is v[i] = (i even) ? arr[i/2] : 0, i.e. every other residue is
        // seen as type 0 (ALA) with zero atoms. We reproduce that here so the
        // residue-pair loop and the per-residue atom ranges match the binary.
        let eff_type = |res_types: &[usize], i: usize| -> usize {
            if i % 2 == 0 { res_types[i / 2] } else { 0 }
        };
        // atom count of residue i as seen by the C binary
        let rec_eff_cnt = |i: usize| -> usize {
            if i % 2 == 0 {
                self.receptor.residue_atom_ranges[i / 2].1 - self.receptor.residue_atom_ranges[i / 2].0
            } else { 0 }
        };
        let lig_eff_cnt = |j: usize| -> usize {
            if j % 2 == 0 {
                self.ligand.residue_atom_ranges[j / 2].1 - self.ligand.residue_atom_ranges[j / 2].0
            } else { 0 }
        };
        // running atom offset: because odd residues contribute 0 atoms, the
        // C offsets accumulate as sum of even counts only.
        let mut rec_off = Vec::with_capacity(n_rec_res);
        {
            let mut acc = 0usize;
            for i in 0..n_rec_res {
                rec_off.push(acc);
                acc += rec_eff_cnt(i);
            }
        }
        let mut lig_off = Vec::with_capacity(n_lig_res);
        {
            let mut acc = 0usize;
            for j in 0..n_lig_res {
                lig_off.push(acc);
                acc += lig_eff_cnt(j);
            }
        }

        for i in 0..n_rec_res {
            let ri = eff_type(&self.receptor.residue_types, i);
            let ra_start = rec_off[i];
            let ra_end = ra_start + rec_eff_cnt(i);

            for j in 0..n_lig_res {
                let rj = eff_type(&self.ligand.residue_types, j);
                let la_start = lig_off[j];
                let la_end = la_start + lig_eff_cnt(j);
                let mut contacted = false;

                for ai in ra_start..ra_end {
                    let rc = &rec_coords[ai];
                    for aj in la_start..la_end {
                        let lc = &lig_c[aj];
                        let dx = rc[0] - lc[0];
                        let dx2 = dx * dx;
                        if dx2 > DIST2_CUTOFF { continue; }
                        let dy = rc[1] - lc[1];
                        let dy2 = dy * dy;
                        if dy2 > DIST2_CUTOFF { continue; }
                        let dz = rc[2] - lc[2];
                        let dz2 = dz * dz;
                        if dz2 > DIST2_CUTOFF { continue; }
                        let dist2 = dx2 + dy2 + dz2;

                        if dist2 < DIST2_CUTOFF {
                            // C: break only the inner (atom_j) loop, so every
                            // contacting receptor atom adds the energy once.
                            total_sipper += SIPPER_ENERGY[ri][rj];
                            total_oda += self.receptor.oda[i] + self.ligand.oda[j];
                            contacted = true;
                            break;
                        }
                        if dist2 <= INTERFACE_CUTOFF2 {
                            iface_r[ai] = 1;
                            iface_l[aj] = 1;
                        }
                    }
                }
                let _ = contacted;
            }
        }

        let score = -1.0 * (total_sipper - ODA_WEIGHT * total_oda);
        let perc_r = satisfied_restraints(iface_r, &self.receptor.active_restraints);
        let perc_l = satisfied_restraints(iface_l, &self.ligand.active_restraints);
        score + perc_r * score + perc_l * score
        })
    }
}
