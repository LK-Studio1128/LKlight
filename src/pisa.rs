use super::constants::INTERFACE_CUTOFF;
use super::qt::{rot3_apply, Quaternion};
use super::scoring::{satisfied_restraints, Score};
use pdbtbx::PDB;
use std::cell::RefCell;
use std::collections::HashMap;

// PISA: atomic statistical potential for protein-protein docking
// Ref: Viswanath et al., Proteins 81(4), 2012
//
// 32 atom types (1-32 in Python, used as 0-31 here), 3 distance bins
// Spline bins: [2.0,3.0) [3.0,4.0) [4.0,4.5) [4.5,6.0) [6.0,8.0)
// Energy range: [2.0, 8.0] Å
// Parameters loaded from data/pisa.params (upper triangle, 1584 values)

const NUM_ATOM_TYPES: usize = 32;
const NUM_BINS: usize = 3;
const MIN_DIST: f64 = 2.0;
const MAX_DIST: f64 = 8.0;
const INTERFACE_CUTOFF2: f64 = INTERFACE_CUTOFF * INTERFACE_CUTOFF;

lazy_static::lazy_static! {
    static ref PISA_ENERGY_TABLE: Vec<f64> = {
        let raw = include_bytes!("../data/pisa.params");
        let contents = std::str::from_utf8(raw).expect("pisa.params UTF-8 error");
        let mut values = contents.lines().filter_map(|l| l.trim().parse::<f64>().ok());
        let mut energy = vec![0.0f64; NUM_ATOM_TYPES * NUM_ATOM_TYPES * NUM_BINS];
        for i in 0..NUM_ATOM_TYPES {
            for j in i..NUM_ATOM_TYPES {
                for r in 0..NUM_BINS {
                    let v = values.next().unwrap_or(0.0);
                    energy[i * NUM_ATOM_TYPES * NUM_BINS + j * NUM_BINS + r] = v;
                    energy[j * NUM_ATOM_TYPES * NUM_BINS + i * NUM_BINS + r] = v;
                }
            }
        }
        energy
    };
}

const R1_SPLINE: [f64; 5] = [2.0, 3.0, 4.0, 4.5, 6.0];
const R2_SPLINE: [f64; 5] = [3.0, 4.0, 4.5, 6.0, 8.0];

fn get_distance_to_bin(dist: f64) -> i32 {
    if dist >= MAX_DIST || dist < MIN_DIST { return -1; }
    let mut i = 4i32;
    while dist < R1_SPLINE[i as usize] || dist >= R2_SPLINE[i as usize] {
        i -= 1;
    }
    i
}

fn get_atom_type(atom_name: &str, res_name: &str) -> i32 {
    if res_name == "LYS" && atom_name == "NZ" { return 1; }
    if atom_name == "N" { return 2; }
    if atom_name == "C"
        || (res_name == "ASN" && atom_name == "CG")
        || (res_name == "GLN" && atom_name == "CD")
    { return 3; }
    if atom_name == "O"
        || (res_name == "ASN" && atom_name == "OD1")
        || (res_name == "GLN" && (atom_name == "OE1" || atom_name == "OE"))
    { return 4; }
    if atom_name == "CA" && res_name != "PRO" { return 5; }
    if (res_name == "ALA" && atom_name == "CB")
        || (res_name == "ILE" && (atom_name == "CG2" || atom_name == "CD1"))
        || (res_name == "LEU" && (atom_name == "CD1" || atom_name == "CD2" || atom_name == "CD" || atom_name == "CE"))
        || (res_name == "THR" && atom_name == "CG2")
        || (res_name == "VAL" && (atom_name == "CG1" || atom_name == "CG2"))
    { return 6; }
    if (res_name == "ARG" && (atom_name == "CB" || atom_name == "CG"))
        || (res_name == "ASN" && atom_name == "CB")
        || (res_name == "GLN" && (atom_name == "CB" || atom_name == "CG"))
        || (res_name == "GLU" && atom_name == "CB")
        || (res_name == "HIS" && atom_name == "CB")
        || (res_name == "ILE" && (atom_name == "CB" || atom_name == "CG1"))
        || (res_name == "LEU" && (atom_name == "CB" || atom_name == "CG"))
        || (res_name == "LYS" && (atom_name == "CB" || atom_name == "CG" || atom_name == "CD"))
        || (res_name == "MET" && atom_name == "CB")
        || (res_name == "PHE" && atom_name == "CB")
        || (res_name == "PRO" && (atom_name == "CB" || atom_name == "CG"))
        || (res_name == "TRP" && atom_name == "CB")
        || (res_name == "TYR" && atom_name == "CB")
        || (res_name == "VAL" && atom_name == "CB")
    { return 7; }
    if (res_name == "PHE" && matches!(atom_name, "CG"|"CD1"|"CD2"|"CE1"|"CE2"|"CZ"))
        || (res_name == "TRP" && matches!(atom_name, "CE3"|"CZ2"|"CZ3"|"CH2"))
        || (res_name == "TYR" && matches!(atom_name, "CG"|"CD1"|"CD2"|"CE1"|"CE2"|"CC"|"CD"|"CE"|"CH"))
    { return 8; }
    if res_name == "TYR" && (atom_name == "CZ" || atom_name == "CF") { return 9; }
    if (res_name == "SER" && atom_name == "OG")
        || (res_name == "THR" && atom_name == "OG1")
        || (res_name == "TYR" && atom_name == "OH")
    { return 10; }
    if res_name == "TRP" && (atom_name == "CG" || atom_name == "CD2") { return 11; }
    if res_name == "TRP" && (atom_name == "CD1" || atom_name == "CE2") { return 12; }
    if res_name == "TRP" && atom_name == "NE1" { return 13; }
    if res_name == "MET" && atom_name == "CG" { return 14; }
    if res_name == "MET" && matches!(atom_name, "SD"|"S"|"SE") { return 15; }
    if res_name == "MET" && atom_name == "CE" { return 16; }
    if res_name == "LYS" && (atom_name == "CE" || atom_name == "CZ") { return 17; }
    if (res_name == "SER" && atom_name == "CB") || (res_name == "THR" && atom_name == "CB") { return 18; }
    if res_name == "PRO" && (atom_name == "CD" || atom_name == "CA") { return 19; }
    if res_name == "CYS" && atom_name == "CB" { return 20; }
    if res_name == "CYS" && matches!(atom_name, "SG"|"S"|"SE") { return 21; }
    if res_name == "HIS" && (atom_name == "CG" || atom_name == "CD2") { return 22; }
    if res_name == "HIS" && (atom_name == "ND1" || atom_name == "NE2") { return 23; }
    if res_name == "HIS" && atom_name == "CE1" { return 24; }
    if res_name == "ARG" && atom_name == "CD" { return 25; }
    if res_name == "ARG" && atom_name == "NE" { return 26; }
    if res_name == "ARG" && atom_name == "CZ" { return 27; }
    if res_name == "ARG" && (atom_name == "NH1" || atom_name == "NH2") { return 28; }
    if (res_name == "ASN" && atom_name == "ND2") || (res_name == "GLN" && atom_name == "NE2") { return 29; }
    if (res_name == "ASP" && atom_name == "CB") || (res_name == "GLU" && atom_name == "CG") { return 30; }
    if (res_name == "ASP" && atom_name == "CG")
        || (res_name == "GLU" && matches!(atom_name, "CD"|"CD1"))
    { return 31; }
    if atom_name == "OXT"
        || (res_name == "ASP" && (atom_name == "OD1" || atom_name == "OD2"))
        || (res_name == "GLU" && matches!(atom_name, "OE1"|"OE2"|"OE11"|"OE21"|"OE"))
        || atom_name == "OX2"
    { return 32; }
    -1
}

pub struct PISADockingModel {
    pub coordinates: Vec<[f64; 3]>,
    pub atom_types: Vec<usize>,
    pub active_restraints: HashMap<String, Vec<usize>>,
    pub passive_restraints: HashMap<String, Vec<usize>>,
    pub num_anm: usize,
    pub nmodes: Vec<f64>,
}

impl PISADockingModel {
    fn new(
        structure: &PDB,
        active_restraints: &[String],
        passive_restraints: &[String],
        nmodes: &[f64],
        num_anm: usize,
    ) -> PISADockingModel {
        let mut model = PISADockingModel {
            coordinates: Vec::new(),
            atom_types: Vec::new(),
            active_restraints: HashMap::new(),
            passive_restraints: HashMap::new(),
            num_anm,
            nmodes: nmodes.to_owned(),
        };

        for chain in structure.chains() {
            for residue in chain.residues() {
                let res_name = match residue.name() { Some(n) => n, None => continue };
                let mut res_id = format!("{}.{}.{}", chain.id(), res_name, residue.serial_number());
                if let Some(c) = residue.insertion_code() { res_id.push_str(c); }

                let mut res_atom_indices: Vec<usize> = Vec::new();

                for atom in residue.atoms() {
                    let aname = atom.name().trim();
                    let atype = get_atom_type(aname, res_name);
                    if atype == -1 { continue; }
                    let idx = model.coordinates.len();
                    model.coordinates.push([atom.x(), atom.y(), atom.z()]);
                    model.atom_types.push((atype - 1) as usize);
                    res_atom_indices.push(idx);
                }

                if !res_atom_indices.is_empty() {
                    if active_restraints.contains(&res_id) {
                        model.active_restraints.insert(res_id.clone(), res_atom_indices.clone());
                    }
                    if passive_restraints.contains(&res_id) {
                        model.passive_restraints.insert(res_id.clone(), res_atom_indices);
                    }
                }
            }
        }
        model
    }
}

pub struct PISA {
    pub receptor: PISADockingModel,
    pub ligand: PISADockingModel,
    pub use_anm: bool,
}

impl<'a> PISA {
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
        let _ = &*PISA_ENERGY_TABLE; // ensure pre-loaded
        let p = PISA {
            receptor: PISADockingModel::new(
                &receptor, &rec_active_restraints, &rec_passive_restraints, &rec_nmodes, rec_num_anm,
            ),
            ligand: PISADockingModel::new(
                &ligand, &lig_active_restraints, &lig_passive_restraints, &lig_nmodes, lig_num_anm,
            ),
            use_anm,
        };
        Box::new(p)
    }
}

impl Score for PISA {
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
        let energy_table = &*PISA_ENERGY_TABLE;
        let rot_mat = rotation.to_matrix();
        let max_dist2 = MAX_DIST * MAX_DIST;

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

            let mut num_contacts = vec![0.0f64; NUM_ATOM_TYPES * NUM_ATOM_TYPES * NUM_BINS];

            const CELL: f64 = MAX_DIST;
            let mut grid: HashMap<(i32,i32,i32), Vec<usize>> = HashMap::with_capacity(rec_n);
            for (i, c) in rec_c.iter().enumerate() {
                let k = ((c[0]/CELL).floor() as i32, (c[1]/CELL).floor() as i32, (c[2]/CELL).floor() as i32);
                grid.entry(k).or_default().push(i);
            }

            for (j, lc) in lig_c.iter().enumerate() {
                let cx = (lc[0]/CELL).floor() as i32;
                let cy = (lc[1]/CELL).floor() as i32;
                let cz = (lc[2]/CELL).floor() as i32;
                let jt0 = self.ligand.atom_types[j];
                for dx in -1..=1i32 { for dy in -1..=1i32 { for dz in -1..=1i32 {
                    if let Some(cells) = grid.get(&(cx+dx, cy+dy, cz+dz)) {
                        for &i in cells {
                            let rc = &rec_c[i];
                            let ddx = rc[0]-lc[0]; let ddy = rc[1]-lc[1]; let ddz = rc[2]-lc[2];
                            let dist2 = ddx*ddx + ddy*ddy + ddz*ddz;
                            if dist2 > max_dist2 { continue; }
                            let dist = dist2.sqrt();
                            if dist2 <= INTERFACE_CUTOFF2 { iface_r[i] = 1; iface_l[j] = 1; }
                            let r = get_distance_to_bin(dist);
                            if r < 0 { continue; }
                            let it0 = self.receptor.atom_types[i];
                            let (it, jt) = if jt0 < it0 { (jt0, it0) } else { (it0, jt0) };
                            let base = it * NUM_ATOM_TYPES * NUM_BINS + jt * NUM_BINS;
                            match r {
                                0 => num_contacts[base]     += 1.0,
                                1 => { num_contacts[base]   += 4.0-dist; num_contacts[base+1] += dist-3.0; }
                                2 => num_contacts[base+1]   += 1.0,
                                3 => { num_contacts[base+1] += 4.0-(dist/3.0); num_contacts[base+2] += (dist/3.0)-3.0; }
                                4 => num_contacts[base+2]   += 1.0,
                                _ => {}
                            }
                        }
                    }
                }}}
            }

            let mut total_energy = 0.0f64;
            for i in 0..NUM_ATOM_TYPES {
                for j in 0..NUM_ATOM_TYPES {
                    for r in 0..NUM_BINS {
                        let nc = num_contacts[i * NUM_ATOM_TYPES * NUM_BINS + j * NUM_BINS + r];
                        if nc != 0.0 {
                            total_energy += nc * energy_table[i * NUM_ATOM_TYPES * NUM_BINS + j * NUM_BINS + r];
                        }
                    }
                }
            }

            let score = total_energy * -1.0;
            let perc_r = satisfied_restraints(iface_r, &self.receptor.active_restraints);
            let perc_l = satisfied_restraints(iface_l, &self.ligand.active_restraints);
            score + perc_r * score + perc_l * score
        })
    }
}
