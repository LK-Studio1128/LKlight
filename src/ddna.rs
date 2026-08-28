use super::constants::INTERFACE_CUTOFF;
use super::qt::{rot3_apply, Quaternion};
use super::scoring::{satisfied_restraints, Score};
use pdbtbx::PDB;
use std::cell::RefCell;
use std::collections::HashMap;

// DDNA: distance-dependent knowledge-based potential for protein-DNA docking
// Ref: Zhang et al., J Med Chem 48, 2325-2335 (2005)
//
// 19 atom types, 21 distance bins, 3-D potential table [21][20][20]
// Distance binning: jj = MAP[int(d * 2.0)], used when 0 < jj <= 20
// Energy: (sum * 0.0021297 - 5.4738) * -1.0
// Score += restraint bias

const DIST_CUTOFF: f64 = 12.0;
const INTERFACE_CUTOFF2: f64 = INTERFACE_CUTOFF * INTERFACE_CUTOFF;
const ENERGY_SCALE: f64 = 0.0021297;
const ENERGY_OFFSET: f64 = -5.4738;

const NUM_BINS: usize = 21;
const TABLE_SIZE: usize = NUM_BINS * 20 * 20; // matches Python: 21*20*20

lazy_static::lazy_static! {
    static ref DDNA_DIST_MAP: Vec<i32> = {
        let mut m = vec![-1i32; 700];
        for i in 0..700usize {
            m[i] = if i == 0 { -1 }          // Python _createmap loops range(1,50), map[0] stays -1
                   else if i < 4 { 1 }
                   else if i < 16 { (i as i32) - 2 }
                   else if i < 50 { (i as i32 / 2) + 6 }
                   else { -1 };
        }
        m
    };
    static ref DDNA_POTENTIALS: Vec<f64> = {
        let raw = include_bytes!("../data/fort.21_xscore_noH_Met");
        let contents = std::str::from_utf8(raw).expect("DDNA potentials UTF-8 error");
        let mut pot = vec![0.0f64; TABLE_SIZE];
        let d1 = 20usize * 20;
        let d2 = 20usize;
        for line in contents.lines() {
            let fields: Vec<&str> = line.split_whitespace().collect();
            if fields.len() == 6 {
                let energy: f64 = match fields[2].parse() { Ok(v) => v, Err(_) => continue };
                let m: usize    = match fields[3].parse() { Ok(v) => v, Err(_) => continue };
                let i: usize    = match fields[4].parse() { Ok(v) => v, Err(_) => continue };
                let j: usize    = match fields[5].parse() { Ok(v) => v, Err(_) => continue };
                if m < NUM_BINS && i < 20 && j < 20 {
                    pot[m * d1 + i * d2 + j] = energy;
                    pot[m * d1 + j * d2 + i] = energy;
                }
            }
        }
        pot
    };
}

// Python atom_types list order (index = ddna atom type)
fn atom_type_index(t: &str) -> Option<usize> {
    match t {
        "C.3"   => Some(0),
        "C.2"   => Some(1),
        "C.ar"  => Some(2),
        "C.cat" => Some(3),
        "N.4"   => Some(4),
        "N.2"   => Some(5),
        "N.ar"  => Some(6),
        "N.am"  => Some(7),
        "N.pl3" => Some(8),
        "O.3"   => Some(9),
        "O.2"   => Some(10),
        "O.co2" => Some(11),
        "S.3"   => Some(12),
        "S.o2"  => Some(13),
        "P.3"   => Some(14),
        "F"     => Some(15),
        "Cl"    => Some(16),
        "Br"    => Some(17),
        "Met"   => Some(18),
        _ => None,
    }
}

// Map residue+atom to DDNA atom type string
fn atom_map(res: &str, atom: &str) -> Option<&'static str> {
    match (res, atom) {
        // CYS
        ("CYS","N")=>"N.am", ("CYS","CA")=>"C.3", ("CYS","C")=>"C.2", ("CYS","O")=>"O.2",
        ("CYS","CB")=>"C.3", ("CYS","SG")=>"S.3",
        // MET
        ("MET","N")=>"N.am", ("MET","CA")=>"C.3", ("MET","C")=>"C.2", ("MET","O")=>"O.2",
        ("MET","CB")=>"C.3", ("MET","CG")=>"C.3", ("MET","SD")=>"S.3", ("MET","CE")=>"C.3",
        // PHE
        ("PHE","N")=>"N.am", ("PHE","CA")=>"C.3", ("PHE","C")=>"C.2", ("PHE","O")=>"O.2",
        ("PHE","CB")=>"C.3", ("PHE","CG")=>"C.ar", ("PHE","CD1")=>"C.ar", ("PHE","CD2")=>"C.ar",
        ("PHE","CE1")=>"C.ar", ("PHE","CE2")=>"C.ar", ("PHE","CZ")=>"C.ar",
        // ILE
        ("ILE","N")=>"N.am", ("ILE","CA")=>"C.3", ("ILE","C")=>"C.2", ("ILE","O")=>"O.2",
        ("ILE","CB")=>"C.3", ("ILE","CG1")=>"C.3", ("ILE","CG2")=>"C.3",
        ("ILE","CD")=>"C.3", ("ILE","CD1")=>"C.3",
        // LEU
        ("LEU","N")=>"N.am", ("LEU","CA")=>"C.3", ("LEU","C")=>"C.2", ("LEU","O")=>"O.2",
        ("LEU","CB")=>"C.3", ("LEU","CG")=>"C.3", ("LEU","CD1")=>"C.3", ("LEU","CD2")=>"C.3",
        // VAL
        ("VAL","N")=>"N.am", ("VAL","CA")=>"C.3", ("VAL","C")=>"C.2", ("VAL","O")=>"O.2",
        ("VAL","CB")=>"C.3", ("VAL","CG1")=>"C.3", ("VAL","CG2")=>"C.3",
        // TRP
        ("TRP","N")=>"N.am", ("TRP","CA")=>"C.3", ("TRP","C")=>"C.2", ("TRP","O")=>"O.2",
        ("TRP","CB")=>"C.3", ("TRP","CG")=>"C.2", ("TRP","CD1")=>"C.2", ("TRP","CD2")=>"C.ar",
        ("TRP","NE1")=>"N.pl3", ("TRP","CE2")=>"C.ar", ("TRP","CE3")=>"C.ar",
        ("TRP","CZ2")=>"C.ar", ("TRP","CZ3")=>"C.ar", ("TRP","CH2")=>"C.ar",
        // TYR
        ("TYR","N")=>"N.am", ("TYR","CA")=>"C.3", ("TYR","C")=>"C.2", ("TYR","O")=>"O.2",
        ("TYR","CB")=>"C.3", ("TYR","CG")=>"C.ar", ("TYR","CD1")=>"C.ar", ("TYR","CD2")=>"C.ar",
        ("TYR","CE1")=>"C.ar", ("TYR","CE2")=>"C.ar", ("TYR","CZ")=>"C.ar", ("TYR","OH")=>"O.3",
        // ALA
        ("ALA","N")=>"N.am", ("ALA","CA")=>"C.3", ("ALA","C")=>"C.2", ("ALA","O")=>"O.2", ("ALA","CB")=>"C.3",
        // GLY
        ("GLY","N")=>"N.am", ("GLY","CA")=>"C.3", ("GLY","C")=>"C.2", ("GLY","O")=>"O.2",
        // THR
        ("THR","N")=>"N.am", ("THR","CA")=>"C.3", ("THR","C")=>"C.2", ("THR","O")=>"O.2",
        ("THR","CB")=>"C.3", ("THR","OG1")=>"O.3", ("THR","CG2")=>"C.3",
        // SER
        ("SER","N")=>"N.am", ("SER","CA")=>"C.3", ("SER","C")=>"C.2", ("SER","O")=>"O.2",
        ("SER","CB")=>"C.3", ("SER","OG")=>"O.3",
        // GLN
        ("GLN","N")=>"N.am", ("GLN","CA")=>"C.3", ("GLN","C")=>"C.2", ("GLN","O")=>"O.2",
        ("GLN","CB")=>"C.3", ("GLN","CG")=>"C.3", ("GLN","CD")=>"C.2", ("GLN","OE1")=>"O.2", ("GLN","NE2")=>"N.am",
        // ASN
        ("ASN","N")=>"N.am", ("ASN","CA")=>"C.3", ("ASN","C")=>"C.2", ("ASN","O")=>"O.2",
        ("ASN","CB")=>"C.3", ("ASN","CG")=>"C.2", ("ASN","OD1")=>"O.2", ("ASN","ND2")=>"N.am",
        // GLU
        ("GLU","N")=>"N.am", ("GLU","CA")=>"C.3", ("GLU","C")=>"C.2", ("GLU","O")=>"O.2",
        ("GLU","CB")=>"C.3", ("GLU","CG")=>"C.3", ("GLU","CD")=>"C.2", ("GLU","OE1")=>"O.co2", ("GLU","OE2")=>"O.co2",
        // ASP
        ("ASP","N")=>"N.am", ("ASP","CA")=>"C.3", ("ASP","C")=>"C.2", ("ASP","O")=>"O.2",
        ("ASP","CB")=>"C.3", ("ASP","CG")=>"C.2", ("ASP","OD1")=>"O.co2", ("ASP","OD2")=>"O.co2",
        // HIS
        ("HIS","N")=>"N.am", ("HIS","CA")=>"C.3", ("HIS","C")=>"C.2", ("HIS","O")=>"O.2",
        ("HIS","CB")=>"C.3", ("HIS","CG")=>"C.2", ("HIS","ND1")=>"N.pl3", ("HIS","CD2")=>"C.2",
        ("HIS","CE1")=>"C.2", ("HIS","NE2")=>"N.2",
        // ARG
        ("ARG","N")=>"N.am", ("ARG","CA")=>"C.3", ("ARG","C")=>"C.2", ("ARG","O")=>"O.2",
        ("ARG","CB")=>"C.3", ("ARG","CG")=>"C.3", ("ARG","CD")=>"C.3", ("ARG","NE")=>"N.pl3",
        ("ARG","CZ")=>"C.cat", ("ARG","NH1")=>"N.pl3", ("ARG","NH2")=>"N.pl3",
        // LYS
        ("LYS","N")=>"N.am", ("LYS","CA")=>"C.3", ("LYS","C")=>"C.2", ("LYS","O")=>"O.2",
        ("LYS","CB")=>"C.3", ("LYS","CG")=>"C.3", ("LYS","CD")=>"C.3", ("LYS","CE")=>"C.3", ("LYS","NZ")=>"N.4",
        // PRO
        ("PRO","N")=>"N.am", ("PRO","CA")=>"C.3", ("PRO","C")=>"C.2", ("PRO","O")=>"O.2",
        ("PRO","CB")=>"C.3", ("PRO","CG")=>"C.3", ("PRO","CD")=>"C.3",
        // DNA/RNA nucleotides (T, DT, A, DA, G, DG, C, DC) - matches the
        // Python DDNAAdapter atom_map table exactly (driver.py lines 204-389)
        ("T","P")|("DT","P")|("A","P")|("DA","P")|("G","P")|("DG","P")|("C","P")|("DC","P") => "P.3",
        ("T","O1P")|("T","O2P")|("DT","O1P")|("DT","O2P")|("A","O1P")|("A","O2P")|("DA","O1P")|("DA","O2P")|("G","O1P")|("G","O2P")|("DG","O1P")|("DG","O2P")|("C","O1P")|("C","O2P")|("DC","O1P")|("DC","O2P") => "O.co2",
        ("T","O5*")|("T","O4*")|("T","O3*")|("DT","O5'")|("DT","O4'")|("DT","O3'")|("A","O5*")|("A","O4*")|("A","O3*")|("A","O2*")|("DA","O5'")|("DA","O4'")|("DA","O3'")|("DA","O2'")|("G","O5*")|("G","O4*")|("G","O3*")|("G","O2*")|("DG","O5'")|("DG","O4'")|("DG","O3'")|("DG","O2'")|("C","O5*")|("C","O4*")|("C","O3*")|("C","O2*")|("DC","O5'")|("DC","O4'")|("DC","O3'")|("DC","O2'") => "O.3",
        ("T","C5*")|("T","C4*")|("T","C3*")|("T","C2*")|("T","C1*")|("T","C5M")|("DT","C5'")|("DT","C4'")|("DT","C3'")|("DT","C2'")|("DT","C1'")|("DT","C7")|("A","C5*")|("A","C4*")|("A","C3*")|("A","C2*")|("A","C1*")|("DA","C5'")|("DA","C4'")|("DA","C3'")|("DA","C2'")|("DA","C1'")|("G","C5*")|("G","C4*")|("G","C3*")|("G","C2*")|("G","C1*")|("DG","C5'")|("DG","C4'")|("DG","C3'")|("DG","C2'")|("DG","C1'")|("C","C5*")|("C","C4*")|("C","C3*")|("C","C2*")|("C","C1*")|("DC","C5'")|("DC","C4'")|("DC","C3'")|("DC","C2'")|("DC","C1'") => "C.3",
        ("T","N1")|("T","N3")|("DT","N1")|("DT","N3")|("G","N1")|("DG","N1")|("C","N1")|("DC","N1") => "N.am",
        ("T","C2")|("T","C4")|("T","C5")|("T","C6")|("DT","C2")|("DT","C4")|("DT","C5")|("DT","C6")|("A","C8")|("DA","C8")|("G","C8")|("G","C5")|("G","C6")|("G","C2")|("G","C4")|("DG","C8")|("DG","C5")|("DG","C6")|("DG","C2")|("DG","C4")|("C","C2")|("C","C4")|("C","C5")|("C","C6")|("DC","C2")|("DC","C4")|("DC","C5")|("DC","C6") => "C.2",
        ("T","O2")|("T","O4")|("DT","O2")|("DT","O4")|("G","O6")|("DG","O6")|("C","O2")|("DC","O2") => "O.2",
        ("A","N9")|("A","N6")|("DA","N9")|("DA","N6")|("G","N9")|("G","N2")|("DG","N9")|("DG","N2")|("C","N4")|("DC","N4") => "N.pl3",
        ("A","N7")|("DA","N7")|("G","N7")|("G","N3")|("DG","N7")|("DG","N3")|("C","N3")|("DC","N3") => "N.2",
        ("A","C5")|("A","C6")|("A","C2")|("A","C4")|("DA","C5")|("DA","C6")|("DA","C2")|("DA","C4") => "C.ar",
        ("A","N1")|("A","N3")|("DA","N1")|("DA","N3") => "N.ar",
        _ => return None,
    }.into()
}

pub struct DDNADockingModel {
    pub coordinates: Vec<[f64; 3]>,
    pub atom_types: Vec<usize>,
    pub active_restraints: HashMap<String, Vec<usize>>,
    pub passive_restraints: HashMap<String, Vec<usize>>,
    pub num_anm: usize,
    pub nmodes: Vec<f64>,
}

impl DDNADockingModel {
    fn new(
        structure: &PDB,
        active_restraints: &[String],
        passive_restraints: &[String],
        nmodes: &[f64],
        num_anm: usize,
    ) -> DDNADockingModel {
        let mut model = DDNADockingModel {
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
                    if atom.element().map(|e| e.symbol() == "H").unwrap_or(false) { continue; }
                    if aname == "H" || aname.starts_with("H") { continue; }

                    let type_str = match atom_map(res_name, aname) {
                        Some(t) => t,
                        None => continue,
                    };
                    let type_idx = match atom_type_index(type_str) {
                        Some(i) => i,
                        None => continue,
                    };
                    let idx = model.coordinates.len();
                    model.coordinates.push([atom.x(), atom.y(), atom.z()]);
                    model.atom_types.push(type_idx);
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

pub struct DDNA {
    pub receptor: DDNADockingModel,
    pub ligand: DDNADockingModel,
    pub use_anm: bool,
}

impl<'a> DDNA {
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
        let _ = &*DDNA_POTENTIALS; // ensure pre-loaded
        let d = DDNA {
            receptor: DDNADockingModel::new(
                &receptor,
                &rec_active_restraints,
                &rec_passive_restraints,
                &rec_nmodes,
                rec_num_anm,
            ),
            ligand: DDNADockingModel::new(
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

impl Score for DDNA {
    fn energy(
        &self,
        translation: &[f64],
        rotation: &Quaternion,
        rec_nmodes: &[f64],
        lig_nmodes: &[f64],
    ) -> f64 {
        thread_local! {
            static SCRATCH: RefCell<(Vec<[f64;3]>, Vec<[f64;3]>, Vec<usize>, Vec<usize>, HashMap<(i32,i32,i32),Vec<usize>>)> =
                RefCell::new((Vec::new(), Vec::new(), Vec::new(), Vec::new(), HashMap::new()));
        }
        let potentials = &*DDNA_POTENTIALS;
        let dist_map   = &*DDNA_DIST_MAP;
        let rot_mat    = rotation.to_matrix();
        let cutoff2    = DIST_CUTOFF * DIST_CUTOFF;

        SCRATCH.with(|sc| {
            let mut sc = sc.borrow_mut();
            let (rec_c, lig_c, iface_r, iface_l, grid) = &mut *sc;
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

            const CELL: f64 = DIST_CUTOFF; // cell == cutoff → 27-cell check sufficient
            grid.clear();
            for (i, c) in rec_c.iter().enumerate() {
                let k = ((c[0]/CELL).floor() as i32, (c[1]/CELL).floor() as i32, (c[2]/CELL).floor() as i32);
                grid.entry(k).or_default().push(i);
            }

            let mut energy_sum = 0.0f64;
            for (j, lc) in lig_c.iter().enumerate() {
                let cx = (lc[0]/CELL).floor() as i32;
                let cy = (lc[1]/CELL).floor() as i32;
                let cz = (lc[2]/CELL).floor() as i32;
                let lat = self.ligand.atom_types[j];
                for dx in -1..=1i32 { for dy in -1..=1i32 { for dz in -1..=1i32 {
                    if let Some(cells) = grid.get(&(cx+dx, cy+dy, cz+dz)) {
                        for &i in cells {
                            let rc = &rec_c[i];
                            let dist2 = (rc[0]-lc[0])*(rc[0]-lc[0])
                                      + (rc[1]-lc[1])*(rc[1]-lc[1])
                                      + (rc[2]-lc[2])*(rc[2]-lc[2]);
                            if dist2 <= cutoff2 {
                                let d = dist2.sqrt();
                                // The Cython source declares `cdef unsigned
                                // int d` and assigns the distance, truncating
                                // it to an integer before binning and before
                                // the interface test. Replicate that.
                                let d_int = d as u32;
                                let d_int_f = d_int as f64;
                                if d_int_f <= INTERFACE_CUTOFF {
                                    iface_r[i] = 1;
                                    iface_l[j] = 1;
                                }
                                let bin_key = (d_int_f * 2.0) as usize;
                                if bin_key < 700 {
                                    let jj = dist_map[bin_key];
                                    if jj > 0 && jj <= 20 {
                                        let rat = self.receptor.atom_types[i];
                                        let idx = (jj as usize) * 20 * 20 + rat * 20 + lat;
                                        let mut u = potentials[idx];
                                        if u < -5.0 { u = 0.0; }
                                        energy_sum += u;
                                    }
                                }
                            }
                        }
                    }
                }}}
            }

            let score = (energy_sum * ENERGY_SCALE + ENERGY_OFFSET) * -1.0;
            let perc_r = satisfied_restraints(iface_r, &self.receptor.active_restraints);
            let perc_l = satisfied_restraints(iface_l, &self.ligand.active_restraints);
            score + perc_r * score + perc_l * score
        })
    }
}
