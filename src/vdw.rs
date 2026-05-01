use super::constants::INTERFACE_CUTOFF2;
use super::pydock::PYDOCKDockingModel;
use super::qt::{rot3_apply, Quaternion};
use std::cell::RefCell;
use super::scoring::{satisfied_restraints, Score};
use pdbtbx::PDB;

const VDW_CUTOFF: f64 = 1.0;
const VDW_DIST_CUTOFF: f64 = 10.0;
const VDW_DIST_CUTOFF2: f64 = VDW_DIST_CUTOFF * VDW_DIST_CUTOFF;

pub struct VDW {
    pub receptor: PYDOCKDockingModel,
    pub ligand: PYDOCKDockingModel,
    pub use_anm: bool,
}

impl<'a> VDW {
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
        let d = VDW {
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

impl Score for VDW {
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
        let (rec_c, lig_c, iface_r, iface_l) = &mut *sc;
        let rec_n = self.receptor.coordinates.len();
        let lig_n = self.ligand.coordinates.len();
        if rec_c.len() != rec_n { rec_c.resize(rec_n, [0.0;3]); }
        if lig_c.len() != lig_n { lig_c.resize(lig_n, [0.0;3]); }
        if iface_r.len() != rec_n { iface_r.resize(rec_n, 0); }
        if iface_l.len() != lig_n { iface_l.resize(lig_n, 0); }
        rec_c.copy_from_slice(&self.receptor.coordinates);
        lig_c.copy_from_slice(&self.ligand.coordinates);
        for v in iface_r.iter_mut() { *v = 0; }
        for v in iface_l.iter_mut() { *v = 0; }
        let rec_num_atoms = rec_n;
        let lig_num_atoms = lig_n;

        let lig_nm_n = if self.ligand.num_anm > 0 {
            self.ligand.nmodes.len() / (3 * self.ligand.num_anm)
        } else { lig_num_atoms };
        let rec_nm_n = if self.receptor.num_anm > 0 {
            self.receptor.nmodes.len() / (3 * self.receptor.num_anm)
        } else { rec_num_atoms };

        for (i_atom, coordinate) in lig_c.iter_mut().enumerate() {
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
            for (i_atom, coordinate) in rec_c.iter_mut().enumerate() {
                if i_atom >= rec_nm_n { break; }
                for i_nm in 0..self.receptor.num_anm {
                    let b = i_nm * rec_nm_n * 3 + i_atom * 3;
                    coordinate[0] += self.receptor.nmodes[b]   * rec_nmodes[i_nm];
                    coordinate[1] += self.receptor.nmodes[b+1] * rec_nmodes[i_nm];
                    coordinate[2] += self.receptor.nmodes[b+2] * rec_nmodes[i_nm];
                }
            }
        }

        let mut total_vdw = 0.0;
        for (i, ra) in rec_c.iter().enumerate() {
            let x1 = ra[0];
            let y1 = ra[1];
            let z1 = ra[2];
            for (j, la) in lig_c.iter().enumerate() {
                let distance2 = (x1 - la[0]) * (x1 - la[0])
                    + (y1 - la[1]) * (y1 - la[1])
                    + (z1 - la[2]) * (z1 - la[2]);

                if distance2 <= VDW_DIST_CUTOFF2 {
                    let vdw_energy =
                        (self.receptor.vdw_charges[i] * self.ligand.vdw_charges[j]).sqrt();
                    let vdw_radius = self.receptor.vdw_radii[i] + self.ligand.vdw_radii[j];
                    let p6 = vdw_radius.powi(6) / distance2.powi(3);
                    let mut k = vdw_energy * (p6 * p6 - 2.0 * p6);
                    if k > VDW_CUTOFF {
                        k = VDW_CUTOFF;
                    }
                    total_vdw += k;
                }

                if distance2 <= INTERFACE_CUTOFF2 {
                    iface_r[i] = 1;
                    iface_l[j] = 1;
                }
            }
        }

        let score = total_vdw * -1.0;
        let perc_r = satisfied_restraints(iface_r, &self.receptor.active_restraints);
        let perc_l = satisfied_restraints(iface_l, &self.ligand.active_restraints);
        score + perc_r * score + perc_l * score
        })
    }
}
