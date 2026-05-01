use super::constants::{DEFAULT_NMODES_STEP, DEFAULT_ROTATION_STEP, DEFAULT_TRANSLATION_STEP};
use super::qt::Quaternion;
use super::scoring::Score;
use std::f64;

pub struct Glowworm<'a> {
    pub id: u32,
    pub translation: [f64; 3],
    pub rotation: Quaternion,
    pub rec_nmodes: Vec<f64>,
    pub lig_nmodes: Vec<f64>,
    pub scoring_function: &'a Box<dyn Score>,
    pub rho: f64,
    pub gamma: f64,
    pub beta: f64,
    pub luciferin: f64,
    pub vision_range: f64,
    pub max_vision_range: f64,
    pub max_neighbors: u32,
    pub neighbors: Vec<u32>,
    pub probabilities: Vec<f64>,
    pub scoring: f64,
    pub moved: bool,
    pub step: u32,
    pub use_anm: bool,
}

impl<'a> Glowworm<'a> {
    pub fn new(
        id: u32,
        translation: [f64; 3],
        rotation: Quaternion,
        rec_nmodes: Vec<f64>,
        lig_nmodes: Vec<f64>,
        scoring_function: &'a Box<dyn Score>,
        use_anm: bool,
    ) -> Self {
        Glowworm {
            id,
            translation,
            rotation,
            rec_nmodes,
            lig_nmodes,
            scoring_function,
            rho: 0.5,
            gamma: 0.4,
            beta: 0.08,
            luciferin: 5.0,
            vision_range: 0.2,
            max_vision_range: 5.0,
            max_neighbors: 5,
            neighbors: Vec::new(),
            probabilities: Vec::new(),
            scoring: 0.0,
            moved: false,
            step: 0,
            use_anm,
        }
    }

    pub fn compute_luciferin(&mut self) {
        if self.moved || self.step == 0 {
            self.scoring = self.scoring_function.energy(
                &self.translation,
                &self.rotation,
                &self.rec_nmodes,
                &self.lig_nmodes,
            );
        }
        self.luciferin = (1.0 - self.rho) * self.luciferin + self.gamma * self.scoring;
        self.step += 1;
    }

    pub fn distance(&self, other: &Glowworm) -> f64 {
        let x1 = self.translation[0];
        let x2 = other.translation[0];
        let y1 = self.translation[1];
        let y2 = other.translation[1];
        let z1 = self.translation[2];
        let z2 = other.translation[2];
        ((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2) + (z1 - z2) * (z1 - z2)).sqrt()
    }

    pub fn is_neighbor(&self, other: &Glowworm) -> bool {
        if self.id != other.id && self.luciferin < other.luciferin {
            return distance_sq(self, other) < self.vision_range * self.vision_range;
        }
        false
    }

    pub fn update_vision_range(&mut self) {
        self.vision_range = (self.max_vision_range).min((0_f64).max(
            self.vision_range
                + self.beta * f64::from(self.max_neighbors as i32 - (self.neighbors.len() as i32)),
        ));
    }

    pub fn compute_probability_moving_toward_neighbor(&mut self, luciferins: &[f64]) {
        self.probabilities.clear();

        let mut total_sum: f64 = 0.0;
        let mut difference: f64;
        for neighbor_id in &self.neighbors {
            difference = luciferins[*neighbor_id as usize] - self.luciferin;
            self.probabilities.push(difference);
            total_sum += difference;
        }

        for i in 0..self.neighbors.len() {
            self.probabilities[i] /= total_sum;
        }
    }

    pub fn select_random_neighbor(&mut self, random_number: f64) -> u32 {
        if self.neighbors.is_empty() {
            return self.id;
        }

        let mut sum_probabilities: f64 = 0.0;
        let mut i: usize = 0;
        while sum_probabilities < random_number {
            sum_probabilities += self.probabilities[i];
            i += 1;
        }
        self.neighbors[i - 1]
    }

    pub fn move_towards(
        &mut self,
        other_id: u32,
        other_position: &[f64],
        other_rotation: &Quaternion,
        other_anm_rec: &[f64],
        other_anm_lig: &[f64],
    ) {
        self.moved = self.id != other_id;
        if self.id != other_id {
            // Translation component
            let dx = other_position[0] - self.translation[0];
            let dy = other_position[1] - self.translation[1];
            let dz = other_position[2] - self.translation[2];
            let norm: f64 = (dx * dx + dy * dy + dz * dz).sqrt();
            let coef: f64 = DEFAULT_TRANSLATION_STEP / norm;
            self.translation[0] += dx * coef;
            self.translation[1] += dy * coef;
            self.translation[2] += dz * coef;

            // Rotation component
            self.rotation = self.rotation.slerp(other_rotation, DEFAULT_ROTATION_STEP);
            self.rotation.normalize();

            // ANM component
            if self.use_anm && !self.rec_nmodes.is_empty() {
                let mut cum2 = 0.0f64;
                for (a, &b) in self.rec_nmodes.iter().zip(other_anm_rec.iter()) { cum2 += (b-a)*(b-a); }
                let coef = DEFAULT_NMODES_STEP / cum2.sqrt().max(1e-14);
                for (a, &b) in self.rec_nmodes.iter_mut().zip(other_anm_rec.iter()) { *a += (b - *a) * coef; }
            }
            if self.use_anm && !self.lig_nmodes.is_empty() {
                let mut cum2 = 0.0f64;
                for (a, &b) in self.lig_nmodes.iter().zip(other_anm_lig.iter()) { cum2 += (b-a)*(b-a); }
                let coef = DEFAULT_NMODES_STEP / cum2.sqrt().max(1e-14);
                for (a, &b) in self.lig_nmodes.iter_mut().zip(other_anm_lig.iter()) { *a += (b - *a) * coef; }
            }
        }
    }
}

/// Returns squared Euclidean distance between two glowworms (no sqrt — use for comparisons).
#[inline]
pub fn distance_sq(one: &Glowworm, two: &Glowworm) -> f64 {
    let dx = one.translation[0] - two.translation[0];
    let dy = one.translation[1] - two.translation[1];
    let dz = one.translation[2] - two.translation[2];
    dx * dx + dy * dy + dz * dz
}

/// Returns Euclidean distance between two glowworms.
#[inline]
pub fn distance(one: &Glowworm, two: &Glowworm) -> f64 {
    distance_sq(one, two).sqrt()
}
