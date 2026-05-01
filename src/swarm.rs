use super::glowworm::{distance_sq, Glowworm};
use super::qt::Quaternion;
use super::scoring::Score;
use rand::Rng;
use rayon::prelude::*;
use std::fs::File;
use std::io::{BufWriter, Error, Write};

pub struct Swarm<'a> {
    pub glowworms: Vec<Glowworm<'a>>,
    pos_scratch:   Vec<[f64; 3]>,
    rot_scratch:   Vec<Quaternion>,
}

impl<'a> Default for Swarm<'a> {
    fn default() -> Self {
        Swarm::new()
    }
}

impl<'a> Swarm<'a> {
    pub fn new() -> Self {
        Swarm {
            glowworms:   Vec::new(),
            pos_scratch: Vec::new(),
            rot_scratch: Vec::new(),
        }
    }

    pub fn add_glowworms(
        &mut self,
        positions: &[Vec<f64>],
        scoring: &'a Box<dyn Score>,
        use_anm: bool,
        rec_num_anm: usize,
        lig_num_anm: usize,
    ) {
        for (i, position) in positions.iter().enumerate() {
            // Translation component
            let translation = [position[0], position[1], position[2]];
            // Rotation component
            let rotation = Quaternion::new(position[3], position[4], position[5], position[6]);
            // ANM for receptor
            let mut rec_nmodes: Vec<f64> = Vec::new();
            if use_anm && rec_num_anm > 0 {
                for j in 7..7 + rec_num_anm {
                    rec_nmodes.push(positions[i][j]);
                }
            }
            // ANM for ligand
            let mut lig_nmodes: Vec<f64> = Vec::new();
            if use_anm && lig_num_anm > 0 {
                for j in 7 + rec_num_anm..positions[i].len() {
                    lig_nmodes.push(positions[i][j]);
                }
            }
            let glowworm = Glowworm::new(
                i as u32,
                translation,
                rotation,
                rec_nmodes,
                lig_nmodes,
                scoring,
                use_anm,
            );
            self.glowworms.push(glowworm);
        }
    }

    pub fn update_luciferin(&mut self) {
        self.glowworms.par_iter_mut().for_each(|gw| gw.compute_luciferin());
    }

    pub fn movement_phase(&mut self, rng: &mut rand::prelude::StdRng) {
        let n = self.glowworms.len();

        // ── G1: fill reusable scratch (no malloc after first call) ───────────
        self.pos_scratch.resize(n, [0.0; 3]);
        self.rot_scratch.resize(n, Quaternion::new(1.0, 0.0, 0.0, 0.0));
        for (i, gw) in self.glowworms.iter().enumerate() {
            self.pos_scratch[i] = gw.translation;
            self.rot_scratch[i] = gw.rotation;
        }
        // ANM still cloned (per-glowworm snapshot; empty when use_anm=false)
        let anm_recs: Vec<Vec<f64>> = self.glowworms.iter().map(|gw| gw.rec_nmodes.clone()).collect();
        let anm_ligs: Vec<Vec<f64>> = self.glowworms.iter().map(|gw| gw.lig_nmodes.clone()).collect();

        // ── Parallel neighbor search (sqrt-free, already parallelized) ───────
        let neighbors: Vec<Vec<u32>> = self.glowworms
            .par_iter()
            .map(|g1| {
                let vr2 = g1.vision_range * g1.vision_range;
                self.glowworms
                    .iter()
                    .filter(|g2| g2.id != g1.id && g1.luciferin < g2.luciferin
                                 && distance_sq(g1, g2) < vr2)
                    .map(|g2| g2.id)
                    .collect()
            })
            .collect();

        // ── G1: move neighbor lists (no clone) + compute probabilities ───────
        let luciferins: Vec<f64> = self.glowworms.iter().map(|gw| gw.luciferin).collect();
        for (gw, nbrs) in self.glowworms.iter_mut().zip(neighbors.into_iter()) {
            gw.neighbors = nbrs;
            gw.compute_probability_moving_toward_neighbor(&luciferins);
        }

        // ── G2: pre-generate randoms, then parallel move ─────────────────────
        let randoms: Vec<f64> = (0..n).map(|_| rng.gen()).collect();
        let (gws, pos_s, rot_s) = (&mut self.glowworms, &self.pos_scratch, &self.rot_scratch);
        gws.par_iter_mut()
            .zip(randoms.par_iter())
            .for_each(|(gw, &r)| {
                let nid = gw.select_random_neighbor(r) as usize;
                gw.move_towards(nid as u32, &pos_s[nid], &rot_s[nid], &anm_recs[nid], &anm_ligs[nid]);
                gw.update_vision_range();
            });
    }

    pub fn save(&mut self, step: u32, output_directory: &str) -> Result<(), Error> {
        let path = format!("{}/gso_{}.out", output_directory, step);
        let mut output = BufWriter::new(File::create(path)?);
        writeln!(
            output,
            "#Coordinates  RecID  LigID  Luciferin  Neighbor's number  Vision Range  Scoring"
        )?;
        for glowworm in self.glowworms.iter() {
            write!(
                output,
                "({:.7}, {:.7}, {:.7}, {:.7}, {:.7}, {:.7}, {:.7}",
                glowworm.translation[0],
                glowworm.translation[1],
                glowworm.translation[2],
                glowworm.rotation.w,
                glowworm.rotation.x,
                glowworm.rotation.y,
                glowworm.rotation.z
            )?;
            if glowworm.use_anm && !glowworm.rec_nmodes.is_empty() {
                for i in 0..glowworm.rec_nmodes.len() {
                    write!(output, ", {:.7}", glowworm.rec_nmodes[i])?;
                }
            }
            if glowworm.use_anm && !glowworm.lig_nmodes.is_empty() {
                for i in 0..glowworm.lig_nmodes.len() {
                    write!(output, ", {:.7}", glowworm.lig_nmodes[i])?;
                }
            }
            writeln!(
                output,
                ")    0    0   {:.8}  {:?} {:.3} {:.8}",
                glowworm.luciferin,
                glowworm.neighbors.len(),
                glowworm.vision_range,
                glowworm.scoring
            )?;
        }
        Ok(())
    }
}
