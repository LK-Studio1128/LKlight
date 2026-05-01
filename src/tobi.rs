use super::constants::INTERFACE_CUTOFF;
use super::qt::{rot3_apply, Quaternion};
use std::cell::RefCell;
use super::scoring::{satisfied_restraints, Score};
use pdbtbx::PDB;
use std::collections::HashMap;

// TOBI: 22×22 residue-level contact potential
// Indices 0-19: standard amino acids (ALA..VAL), 20: NH backbone, 21: OC backbone
// Ref: Tobi & Elber, Proteins 2000;41:40-56

const NUM_TYPES: usize = 22;

fn res_to_idx(name: &str) -> Option<usize> {
    match name {
        "ALA" => Some(0),  "ARG" => Some(1),  "ASN" => Some(2),  "ASP" => Some(3),
        "CYS" => Some(4),  "GLN" => Some(5),  "GLU" => Some(6),  "GLY" => Some(7),
        "HIS" => Some(8),  "ILE" => Some(9),  "LEU" => Some(10), "LYS" => Some(11),
        "MET" => Some(12), "PHE" => Some(13), "PRO" => Some(14), "SER" => Some(15),
        "THR" => Some(16), "TRP" => Some(17), "TYR" => Some(18), "VAL" => Some(19),
        _ => None,
    }
}

// Atom type classification for sidechain centroid (Python atom_types flattened)
fn tobi_atom_type(res_name: &str, atom_name: &str) -> bool {
    let key = format!("{}{}", res_name, atom_name);
    matches!(key.as_str(),
        "ALACB"|"ARGCB"|"ASNCB"|"ASPCB"|"CYSCB"|"GLNCB"|"GLUCB"|"HISCB"|"ILECB"|
        "LEUCB"|"LYSCB"|"METCB"|"PHECB"|"PROCB"|"PROCG"|"PROCD"|"THRCB"|"TRPCB"|
        "TYRCB"|"VALCB"|"LYSCE"|"LYSNZ"|"LYSCD"|"ASPCG"|"ASPOD1"|"ASPOD2"|"GLUCD"|
        "GLUOE1"|"GLUOE2"|"ARGCZ"|"ARGNH1"|"ARGNH2"|"ASNCG"|"ASNOD1"|"ASNND2"|
        "GLNCD"|"GLNOE1"|"GLNNE2"|"ARGCD"|"ARGNE"|"SERCB"|"SEROG"|"THROG1"|"TYROH"|
        "HISCG"|"HISND1"|"HISCD2"|"HISCE1"|"HISNE2"|"TRPNE1"|"TYRCE1"|"TYRCE2"|
        "TYRCZ"|"ARGCG"|"GLNCG"|"GLUCG"|"ILECG1"|"LEUCG"|"LYSCG"|"METCG"|"METSD"|
        "PHECG"|"PHECD1"|"PHECD2"|"PHECE1"|"PHECE2"|"PHECZ"|"THRCG2"|"TRPCG"|
        "TRPCD1"|"TRPCD2"|"TRPCE2"|"TRPCE3"|"TRPCZ2"|"TRPCZ3"|"TRPCH2"|"TYRCG"|
        "TYRCD1"|"TYRCD2"|"ILECG2"|"ILECD1"|"ILECD"|"LEUCD1"|"LEUCD2"|"METCE"|
        "VALCG1"|"VALCG2"|"CYSSG"|"N"|"CA"|"C"|"O"|"GLYCA"
    )
}

// Step-1 potential table (22×22, row-major, rows/cols: ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL NH OC)
#[rustfmt::skip]
const TOBI_SC_1: [[f64; NUM_TYPES]; NUM_TYPES] = [
    [-0.59, 0.76,-0.02,-0.76,-0.53,-0.71, 0.95,-0.64,-0.62,-0.90,-0.20, 1.34,-0.75,-1.35,-0.54, 0.02,-0.06,-1.42,-0.88,-1.61, 0.09, 0.27],
    [ 0.76, 0.07, 0.26,-1.59,-0.35, 0.21,-2.01, 0.07, 0.17, 1.03,-0.60, 1.40,-0.52,-1.01,-1.19, 0.37,-0.48,-2.47,-0.96,-0.41, 1.37,-1.33],
    [-0.02, 0.26,-0.34,-1.24,-0.37,-0.25,-0.83, 0.24,-0.45,-0.48, 0.02,-0.61,-0.58, 0.93, 0.27,-0.32,-0.58,-1.27,-1.60, 0.26, 0.19,-0.42],
    [-0.76,-1.59,-1.24, 0.96,-0.47,-0.43, 0.18, 0.12,-0.91, 1.44, 0.39,-2.15,-0.85, 0.71,-0.07,-1.16,-0.70,-1.10,-0.37,-1.12,-0.15, 0.38],
    [-0.53,-0.35,-0.37,-0.47,10.00, 1.55, 1.38, 0.28,-2.26,-0.20,-0.27,-1.03,-0.36,-0.35, 1.81, 0.18, 3.04, 3.27, 1.14,-0.51,-0.63, 1.40],
    [-0.71, 0.21,-0.25,-0.43, 1.55, 0.04, 0.08,-0.65,-0.15,-1.61,-0.21,-0.64,-1.02, 0.67,-0.30,-0.47,-0.41,-1.54,-1.11,-0.78, 0.51,-0.44],
    [ 0.95,-2.01,-0.83, 0.18, 1.38, 0.08, 1.69,-0.62,-0.90,-0.29, 0.46,-1.69,-1.65,-0.95,-0.56,-0.93,-0.10,-1.09,-0.75,-0.42, 0.07,-0.04],
    [-0.64, 0.07, 0.24, 0.12, 0.28,-0.65,-0.62,-0.19,-1.26,-0.19,-0.59,-0.25,-0.69, 0.75, 0.75, 0.40, 0.13,-0.35,-0.70,-1.44,-0.76, 0.60],
    [-0.62, 0.17,-0.45,-0.91,-2.26,-0.15,-0.90,-1.26, 2.41,-1.84,-1.85,-0.08,-1.77,-2.61,-0.62,-0.33,-1.56, 0.38,-1.77,-0.39,-0.01,-0.05],
    [-0.90, 1.03,-0.48, 1.44,-0.20,-1.61,-0.29,-0.19,-1.84,-0.93,-1.54,-0.21,-0.89,-2.24, 0.43,-0.64,-0.50,-2.84,-1.19,-0.85, 1.36,-0.81],
    [-0.20,-0.60, 0.02, 0.39,-0.27,-0.21, 0.46,-0.59,-1.85,-1.54,-1.92,-1.20,-1.60,-1.92,-1.33,-0.60,-0.68,-0.78,-1.33,-1.52, 0.16,-0.67],
    [ 1.34, 1.40,-0.61,-2.15,-1.03,-0.64,-1.69,-0.25,-0.08,-0.21,-1.20, 0.40,-0.06,-0.34,-0.12,-0.67, 0.61,-1.67,-1.22,-0.91, 0.86,-0.71],
    [-0.75,-0.52,-0.58,-0.85,-0.36,-1.02,-1.65,-0.69,-1.77,-0.89,-1.60,-0.06,10.00,-2.98,-2.10,-1.30,-0.38,-0.04,-2.35,-2.25, 0.68,-0.68],
    [-1.35,-1.01, 0.93, 0.71,-0.35, 0.67,-0.95, 0.75,-2.61,-2.24,-1.92,-0.34,-2.98,-2.24,-1.70, 0.57,-1.50, 0.36,-1.66, 0.16, 0.03,-0.21],
    [-0.54,-1.19, 0.27,-0.07, 1.81,-0.30,-0.56, 0.75,-0.62, 0.43,-1.33,-0.12,-2.10,-1.70,-2.64,-0.95,-0.16,-1.63,-1.40,-1.23, 0.54,-0.15],
    [ 0.02, 0.37,-0.32,-1.16, 0.18,-0.47,-0.93, 0.40,-0.33,-0.64,-0.60,-0.67,-1.30, 0.57,-0.95, 0.27,-0.23,-0.89,-0.88, 0.42, 0.02,-0.21],
    [-0.06,-0.48,-0.58,-0.70, 3.04,-0.41,-0.10, 0.13,-1.56,-0.50,-0.68, 0.61,-0.38,-1.50,-0.16,-0.23, 0.60, 0.38,-0.14,-0.04, 0.07,-0.03],
    [-1.42,-2.47,-1.27,-1.10, 3.27,-1.54,-1.09,-0.35, 0.38,-2.84,-0.78,-1.67,-0.04, 0.36,-1.63,-0.89, 0.38, 1.40,-2.23, 0.32,-0.44,-0.46],
    [-0.88,-0.96,-1.60,-0.37, 1.14,-1.11,-0.75,-0.70,-1.77,-1.19,-1.33,-1.22,-2.35,-1.66,-1.40,-0.88,-0.14,-2.23,-1.63,-0.63, 0.29,-0.20],
    [-1.61,-0.41, 0.26,-1.12,-0.51,-0.78,-0.42,-1.44,-0.39,-0.85,-1.52,-0.91,-2.25, 0.16,-1.23, 0.42,-0.04, 0.32,-0.63,-1.45, 0.31,-0.33],
    [ 0.09, 1.37, 0.19,-0.15,-0.63, 0.51, 0.07,-0.76,-0.01, 1.36, 0.16, 0.86, 0.68, 0.03, 0.54, 0.02, 0.07,-0.44, 0.29, 0.31, 0.63,-0.02],
    [ 0.27,-1.33,-0.42, 0.38, 1.40,-0.44,-0.04, 0.60,-0.05,-0.81,-0.67,-0.71,-0.68,-0.21,-0.15,-0.21,-0.03,-0.46,-0.20,-0.33,-0.02,-0.09],
];

#[rustfmt::skip]
const TOBI_SC_2: [[f64; NUM_TYPES]; NUM_TYPES] = [
    [-0.58,-0.40, 0.02, 0.85, 0.54, 0.04,-0.05, 1.56, 0.41,-0.24, 0.12, 0.94,-0.07,-0.15, 1.75, 0.16,-0.36, 0.56, 0.42,-0.17,-0.25, 0.09],
    [-0.40, 0.72, 0.37,-0.56, 0.09, 0.48,-0.64,-0.31, 1.42, 0.16, 0.38, 1.50,-0.42, 1.11, 0.70, 0.43, 0.36,-1.14,-0.20,-0.97,-0.07,-0.19],
    [ 0.02, 0.37, 1.42,-0.06, 1.13, 0.20,-0.43,-0.04, 0.81, 0.47, 0.38,-0.29, 2.71,-0.64,-0.18, 0.33, 0.04, 0.25, 0.24, 0.84,-0.20, 0.24],
    [ 0.85,-0.56,-0.06, 0.22,-0.25,-0.14,-0.18,-0.09,-0.50, 0.23, 0.81,-0.24, 0.20,-1.10, 1.26,-0.27,-0.21, 0.57,-0.42,-0.67, 0.39, 0.30],
    [ 0.54, 0.09, 1.13,-0.25, 1.10,-0.66, 0.59,-0.13, 6.43, 0.12, 0.86, 0.33, 0.55,-0.20, 1.54, 0.84, 1.76, 2.67, 1.37, 0.96,-1.20,-0.21],
    [ 0.04, 0.48, 0.20,-0.14,-0.66, 1.00, 0.09,-0.11, 0.85, 0.42,-0.17, 0.21, 0.34, 0.05, 0.73,-0.45, 0.27,-0.47, 0.60, 0.08, 0.00, 0.00],
    [-0.05,-0.64,-0.43,-0.18, 0.59, 0.09, 0.41, 0.22,-0.12, 1.18, 0.33,-1.10,-0.26, 0.65, 1.02,-0.04,-0.51, 0.85, 0.45, 1.43, 0.37,-0.03],
    [ 1.56,-0.31,-0.04,-0.09,-0.13,-0.11, 0.22,-0.17, 0.51, 0.28, 0.96, 0.71,-0.74, 0.45, 0.06, 0.53, 0.42,-0.15,-0.32, 0.67,-0.39, 0.15],
    [ 0.41, 1.42, 0.81,-0.50, 6.43, 0.85,-0.12, 0.51,-0.19,-0.75,-1.22, 0.78, 4.03, 0.98, 0.07, 0.39, 1.14, 1.54,-0.06,-0.51, 0.00, 0.24],
    [-0.24, 0.16, 0.47, 0.23, 0.12, 0.42, 1.18, 0.28,-0.75,-1.08,-1.28, 0.26,-1.22,-0.78, 1.53, 0.82, 0.54, 0.30,-0.96, 1.55, 0.25,-0.56],
    [ 0.12, 0.38, 0.38, 0.81, 0.86,-0.17, 0.33, 0.96,-1.22,-1.28, 0.26, 0.62,-0.12, 0.37,-0.91, 0.25, 0.00,-0.01, 0.08,-0.54, 0.10, 0.15],
    [ 0.94, 1.50,-0.29,-0.24, 0.33, 0.21,-1.10, 0.71, 0.78, 0.26, 0.62, 1.74, 0.78, 0.22,-0.40, 0.00,-0.13,-0.15, 0.41, 0.93, 0.19, 0.11],
    [-0.07,-0.42, 2.71, 0.20, 0.55, 0.34,-0.26,-0.74, 4.03,-1.22,-0.12, 0.78,10.00, 0.74, 0.81, 0.12, 0.00, 0.30, 1.16, 6.00, 0.05, 0.43],
    [-0.15, 1.11,-0.64,-1.10,-0.20, 0.05, 0.65, 0.45, 0.98,-0.78, 0.37, 0.22, 0.74, 1.23, 0.83,-0.13, 0.35,-0.17,-0.66,-0.57,-0.38,-0.18],
    [ 1.75, 0.70,-0.18, 1.26, 1.54, 0.73, 1.02, 0.06, 0.07, 1.53,-0.91,-0.40, 0.81, 0.83,-0.56,-0.31, 0.06, 1.90,-0.18, 0.52,-0.20, 0.74],
    [ 0.16, 0.43, 0.33,-0.27, 0.84,-0.45,-0.04, 0.53, 0.39, 0.82, 0.25, 0.00, 0.12,-0.13,-0.31,-0.76,-0.46,-0.62, 0.02, 0.37,-0.09, 0.23],
    [-0.36, 0.36, 0.04,-0.21, 1.76, 0.27,-0.51, 0.42, 1.14, 0.54, 0.00,-0.13, 0.00, 0.35, 0.06,-0.46, 0.01,-0.71,-0.38,-0.26,-0.36, 0.28],
    [ 0.56,-1.14, 0.25, 0.57, 2.67,-0.47, 0.85,-0.15, 1.54, 0.30,-0.01,-0.15, 0.30,-0.17, 1.90,-0.62,-0.71,-1.69,-0.59, 0.16, 0.14,-0.29],
    [ 0.42,-0.20, 0.24,-0.42, 1.37, 0.60, 0.45,-0.32,-0.06,-0.96, 0.08, 0.41, 1.16,-0.66,-0.18, 0.02,-0.38,-0.59, 2.07,-0.41,-0.31,-0.08],
    [-0.17,-0.97, 0.84,-0.67, 0.96, 0.08, 1.43, 0.67,-0.51, 1.55,-0.54, 0.93, 6.00,-0.57, 0.52, 0.37,-0.26, 0.16,-0.41, 1.30,-0.04, 0.45],
    [-0.25,-0.07,-0.20, 0.39,-1.20, 0.00, 0.37,-0.39, 0.00, 0.25, 0.10, 0.19, 0.05,-0.38,-0.20,-0.09,-0.36, 0.14,-0.31,-0.04,-0.18, 0.09],
    [ 0.09,-0.19, 0.24, 0.30,-0.21, 0.00,-0.03, 0.15, 0.24,-0.56, 0.15, 0.11, 0.43,-0.18, 0.74, 0.23, 0.28,-0.29,-0.08, 0.45, 0.09,-0.24],
];

pub struct TOBIDockingModel {
    pub coordinates: Vec<[f64; 3]>,
    pub tobi_types: Vec<usize>,
    pub active_restraints: HashMap<String, Vec<usize>>,
    pub passive_restraints: HashMap<String, Vec<usize>>,
}

impl TOBIDockingModel {
    fn new(
        structure: &PDB,
        active_restraints: &[String],
        passive_restraints: &[String],
    ) -> TOBIDockingModel {
        let mut model = TOBIDockingModel {
            coordinates: Vec::new(),
            tobi_types: Vec::new(),
            active_restraints: HashMap::new(),
            passive_restraints: HashMap::new(),
        };

        for chain in structure.chains() {
            for residue in chain.residues() {
                let res_name = match residue.name() {
                    Some(n) => n,
                    None => continue,
                };
                let res_idx = match res_to_idx(res_name) {
                    Some(idx) => idx,
                    None => continue,
                };

                let mut n_coord: Option<[f64; 3]> = None;
                let mut o_coord: Option<[f64; 3]> = None;
                let mut cx = 0.0_f64;
                let mut cy = 0.0_f64;
                let mut cz = 0.0_f64;
                let mut count = 0usize;

                for atom in residue.atoms() {
                    if atom.name() == "H" || atom.element().map(|e| e.symbol() == "H").unwrap_or(false) {
                        continue;
                    }
                    let aname = atom.name().trim();
                    if aname == "N" {
                        n_coord = Some([atom.x(), atom.y(), atom.z()]);
                    }
                    if aname == "O" {
                        o_coord = Some([atom.x(), atom.y(), atom.z()]);
                    }
                    if tobi_atom_type(res_name, aname) {
                        cx += atom.x();
                        cy += atom.y();
                        cz += atom.z();
                        count += 1;
                    }
                }

                if count == 0 {
                    continue;
                }

                let point_idx_start = model.coordinates.len();

                if let Some(nc) = n_coord {
                    model.coordinates.push(nc);
                    model.tobi_types.push(20);
                }
                if let Some(oc) = o_coord {
                    model.coordinates.push(oc);
                    model.tobi_types.push(21);
                }

                let centroid = [cx / count as f64, cy / count as f64, cz / count as f64];
                model.coordinates.push(centroid);
                model.tobi_types.push(res_idx);

                let mut res_id = format!("{}.{}.{}", chain.id(), res_name, residue.serial_number());
                if let Some(c) = residue.insertion_code() {
                    res_id.push_str(c);
                }
                let point_count = model.coordinates.len() - point_idx_start;
                let point_indices: Vec<usize> = (point_idx_start..point_idx_start + point_count).collect();

                if active_restraints.contains(&res_id) {
                    model.active_restraints.insert(res_id.clone(), point_indices.clone());
                }
                if passive_restraints.contains(&res_id) {
                    model.passive_restraints.insert(res_id.clone(), point_indices);
                }
            }
        }
        model
    }
}

pub struct TOBI {
    pub receptor: TOBIDockingModel,
    pub ligand: TOBIDockingModel,
}

impl<'a> TOBI {
    pub fn new(
        receptor: PDB,
        rec_active_restraints: Vec<String>,
        rec_passive_restraints: Vec<String>,
        ligand: PDB,
        lig_active_restraints: Vec<String>,
        lig_passive_restraints: Vec<String>,
    ) -> Box<dyn Score + 'a> {
        Box::new(TOBI {
            receptor: TOBIDockingModel::new(&receptor, &rec_active_restraints, &rec_passive_restraints),
            ligand: TOBIDockingModel::new(&ligand, &lig_active_restraints, &lig_passive_restraints),
        })
    }
}

impl Score for TOBI {
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
        let interface_cutoff2 = INTERFACE_CUTOFF * INTERFACE_CUTOFF;

        SCRATCH.with(|sc| {
        let mut sc = sc.borrow_mut();
        let (lig_c, iface_r, iface_l) = &mut *sc;
        let n_lig = self.ligand.coordinates.len();
        let n_rec = self.receptor.coordinates.len();
        if lig_c.len() != n_lig { lig_c.resize(n_lig, [0.0;3]); }
        if iface_r.len() != n_rec { iface_r.resize(n_rec, 0); }
        if iface_l.len() != n_lig { iface_l.resize(n_lig, 0); }
        for (c, src) in lig_c.iter_mut().zip(self.ligand.coordinates.iter()) {
            let r = rot3_apply(&rot_mat, *src);
            *c = [r[0]+translation[0], r[1]+translation[1], r[2]+translation[2]];
        }
        for v in iface_r.iter_mut() { *v = 0; }
        for v in iface_l.iter_mut() { *v = 0; }

        let rec_coords = &self.receptor.coordinates;
        let mut energy = 0.0_f64;

        const CELL: f64 = 8.0;
        let mut grid: std::collections::HashMap<(i32,i32,i32), Vec<usize>> =
            std::collections::HashMap::with_capacity(n_rec);
        for (i, c) in rec_coords.iter().enumerate() {
            let k = ((c[0]/CELL).floor() as i32, (c[1]/CELL).floor() as i32, (c[2]/CELL).floor() as i32);
            grid.entry(k).or_default().push(i);
        }

        for (j, lc) in lig_c.iter().enumerate() {
            let cx = (lc[0]/CELL).floor() as i32;
            let cy = (lc[1]/CELL).floor() as i32;
            let cz = (lc[2]/CELL).floor() as i32;
            let lt = self.ligand.tobi_types[j];
            for ddx in -1..=1i32 { for ddy in -1..=1i32 { for ddz in -1..=1i32 {
                if let Some(cells) = grid.get(&(cx+ddx, cy+ddy, cz+ddz)) {
                    for &i in cells {
                        let rc = &rec_coords[i];
                        let dx = rc[0]-lc[0]; let dy = rc[1]-lc[1]; let dz = rc[2]-lc[2];
                        let d2 = dx*dx + dy*dy + dz*dz;
                        if d2 <= interface_cutoff2 { iface_r[i] = 1; iface_l[j] = 1; }
                        if d2 <= 64.0 {
                            let rt = self.receptor.tobi_types[i];
                            let rec_bb = rt >= 20;
                            let lig_bb = lt >= 20;
                            if rec_bb && lig_bb {
                                if d2 <= 20.25 { energy += TOBI_SC_1[rt][lt]; }
                                else if d2 <= 36.0  { energy += TOBI_SC_2[rt][lt]; }
                            } else if rec_bb || lig_bb {
                                if d2 <= 30.25 { energy += TOBI_SC_1[rt][lt]; }
                                else if d2 <= 49.0  { energy += TOBI_SC_2[rt][lt]; }
                            } else {
                                if d2 <= 42.25 { energy += TOBI_SC_1[rt][lt]; }
                                else if d2 <= 64.0  { energy += TOBI_SC_2[rt][lt]; }
                            }
                        }
                    }
                }
            }}}
        }

        let score = energy * -1.0;
        let perc_r = satisfied_restraints(iface_r, &self.receptor.active_restraints);
        let perc_l = satisfied_restraints(iface_l, &self.ligand.active_restraints);
        score + perc_r * score + perc_l * score
        })
    }
}
