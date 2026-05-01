/// anm.rs — Anisotropic Network Model (ANM) implementation in pure Rust.
///
/// Replicates the essential behaviour of ProDy's ANM / extendModel pipeline
/// used by the Python LightDock setup script:
///   1. Build 3N×3N Kirchhoff-style Hessian from Cα/P backbone coordinates.
///   2. Full symmetric eigendecomposition (nalgebra).
///   3. Sort modes by eigenvalue (ascending), skip 6 rigid-body modes.
///   4. Extend backbone modes to every heavy atom in the structure.
///   5. Scale modes so that a random sampling produces the desired RMSD.
///   6. Return as a flat Vec<f64> shaped (n_modes, n_atoms, 3).
///   7. Write to .npy (1-D float64 format readable by NumPy and npyz).

use nalgebra::DMatrix;
use std::io::{BufWriter, Write};

/// Default ANM parameters (match ProDy defaults).
const ANM_CUTOFF: f64 = 15.0;   // Å — contact cut-off
const ANM_GAMMA:  f64 = 1.0;    // spring constant

// ─── Public types ─────────────────────────────────────────────────────────────

/// Result of an ANM calculation.
pub struct AnmModes {
    /// Flat array shaped (n_modes, n_atoms, 3).
    /// Sized n_modes * n_atoms * 3 elements (f64).
    pub data: Vec<f64>,
    pub n_modes: usize,
    pub n_atoms: usize,
}

// ─── Core computation ─────────────────────────────────────────────────────────

/// Build the 3N × 3N Hessian matrix from N backbone (Cα / P) coordinates.
fn build_hessian(backbone: &[[f64; 3]]) -> Vec<f64> {
    let n = backbone.len();
    let dim = 3 * n;
    let cutoff2 = ANM_CUTOFF * ANM_CUTOFF;
    let mut h = vec![0.0f64; dim * dim];

    for i in 0..n {
        for j in (i + 1)..n {
            let dx = backbone[j][0] - backbone[i][0];
            let dy = backbone[j][1] - backbone[i][1];
            let dz = backbone[j][2] - backbone[i][2];
            let d2 = dx * dx + dy * dy + dz * dz;
            if d2 > cutoff2 { continue; }

            let coeff = ANM_GAMMA / d2;
            let r = [dx, dy, dz];

            for a in 0..3usize {
                for b in 0..3usize {
                    let val = coeff * r[a] * r[b];
                    // off-diagonal blocks
                    h[(3*i+a)*dim + (3*j+b)] -= val;
                    h[(3*j+a)*dim + (3*i+b)] -= val;
                    // diagonal blocks
                    h[(3*i+a)*dim + (3*i+b)] += val;
                    h[(3*j+a)*dim + (3*j+b)] += val;
                }
            }
        }
    }
    h
}

/// Compute ANM modes for a molecule.
///
/// * `backbone`       – Cα / P coordinates, N_bb atoms.
/// * `atom_to_bb`     – for every heavy atom (index in full structure),
///                      the index of its backbone atom in `backbone`.
/// * `n_atoms_total`  – total heavy atoms (used for sizing the output).
/// * `n_modes`        – number of non-trivial modes requested.
/// * `rmsd`           – target per-mode RMSD (Å) for scaling (default 0.5).
///
/// Returns `AnmModes` whose `data` is shaped `(n_modes, n_atoms_total, 3)`.
pub fn compute_anm(
    backbone:      &[[f64; 3]],
    atom_to_bb:    &[usize],
    n_atoms_total: usize,
    n_modes:       usize,
    rmsd:          f64,
) -> AnmModes {
    let n_bb  = backbone.len();
    let dim   = 3 * n_bb;
    let zeros = AnmModes {
        data:    vec![0.0; n_modes * n_atoms_total * 3],
        n_modes, n_atoms: n_atoms_total,
    };
    if n_bb < 7 { return zeros; }

    // ── Build Hessian ────────────────────────────────────────────────────────
    eprintln!("[ANM] Building {}×{} Hessian for {} backbone atoms …", dim, dim, n_bb);
    let h_data = build_hessian(backbone);
    let hess = DMatrix::from_row_slice(dim, dim, &h_data);

    // ── Eigendecompose ───────────────────────────────────────────────────────
    eprintln!("[ANM] Running symmetric eigendecomposition …");
    let eig = nalgebra::linalg::SymmetricEigen::new(hess);

    // Sort eigenvalue / eigenvector pairs ascending
    let mut order: Vec<usize> = (0..dim).collect();
    order.sort_by(|&a, &b| {
        eig.eigenvalues[a]
            .partial_cmp(&eig.eigenvalues[b])
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    // ── Select non-trivial modes (skip first 6) ──────────────────────────────
    let available  = dim.saturating_sub(6);
    let actual     = n_modes.min(available);

    // Variances = 1 / λ  (InvFreq) for modes 6..6+actual
    let variances: Vec<f64> = (0..actual).map(|m| {
        let ev = eig.eigenvalues[order[6 + m]];
        if ev.abs() < 1e-12 { 0.0 } else { 1.0 / ev }
    }).collect();

    // Scale factor per mode:
    //   coef ≈ sqrt(sum(variances))   (analytical approximation for the
    //          Python coef = mean(|z~N(0,1)|·var^0.5) statistic)
    let sum_var: f64 = variances.iter().sum();
    let coef = if sum_var > 1e-14 { sum_var.sqrt() } else { 1.0 };
    let base_scale = (n_atoms_total as f64).sqrt() * rmsd / coef;

    // ── Extend backbone modes to all atoms and scale ─────────────────────────
    eprintln!("[ANM] Extending {} modes to {} atoms …", actual, n_atoms_total);
    let mut out = vec![0.0f64; n_modes * n_atoms_total * 3];

    for m in 0..actual {
        let col_idx = order[6 + m];
        let var_sqrt = variances[m].sqrt();

        // Magnitude of this eigenvector (should be ≈ 1, but normalise anyway)
        let mag: f64 = {
            let col = eig.eigenvectors.column(col_idx);
            col.iter().map(|v| v * v).sum::<f64>().sqrt()
        };
        let scale = if mag < 1e-14 { 0.0 } else { base_scale * var_sqrt / mag };

        let col = eig.eigenvectors.column(col_idx);
        let base = m * n_atoms_total * 3;

        for a in 0..n_atoms_total {
            let bb = if a < atom_to_bb.len() { atom_to_bb[a] } else { 0 };
            for d in 0..3 {
                out[base + a * 3 + d] = col[3 * bb + d] * scale;
            }
        }
    }
    // modes `actual..n_modes` remain zero (padding)

    eprintln!("[ANM] Done ({} non-trivial modes, {} padded).",
        actual, n_modes.saturating_sub(actual));

    AnmModes { data: out, n_modes, n_atoms: n_atoms_total }
}

// ─── npy writer ───────────────────────────────────────────────────────────────

/// Write a 1-D f64 array in NumPy .npy v1.0 format (little-endian float64).
/// The resulting file is readable by `np.load(...)` and npyz's `NpyFile::new`.
pub fn save_npy(path: &str, data: &[f64]) -> std::io::Result<()> {
    let mut f = BufWriter::new(std::fs::File::create(path)?);

    // Magic: \x93NUMPY + version 1.0
    f.write_all(b"\x93NUMPY\x01\x00")?;

    // Header dict
    let header_raw = format!(
        "{{'descr': '<f8', 'fortran_order': False, 'shape': ({},), }}",
        data.len()
    );
    // Total header area (magic=6 + version=2 + header_len_field=2 + header)
    // must be a multiple of 64.
    // header_len_field(2) + header.len() must be rounded up to next mult-of-64
    // after accounting for the 10-byte prefix.
    let pad_to = {
        let prefix = 6 + 2 + 2; // magic + version + 2-byte len field
        let raw_total = prefix + header_raw.len() + 1; // +1 for trailing '\n'
        ((raw_total + 63) / 64) * 64
    };
    let pad_len = pad_to - (6 + 2 + 2) - 1; // header content length before '\n'
    let header = format!("{:<width$}\n", header_raw, width = pad_len);
    assert_eq!(header.len(), pad_len + 1);

    let hlen = (header.len() as u16).to_le_bytes();
    f.write_all(&hlen)?;
    f.write_all(header.as_bytes())?;

    // Data — write entire payload as a single byte slice
    let byte_slice: &[u8] = unsafe {
        std::slice::from_raw_parts(data.as_ptr() as *const u8, data.len() * 8)
    };
    f.write_all(byte_slice)?;
    Ok(())
}

// ─── Helper: build atom_to_backbone mapping ──────────────────────────────────
// These helpers take simple slices to avoid depending on pdbtbx in this module.

/// Given parallel slices of backbone-atom residue keys and all-atom residue keys,
/// build `atom_to_bb[i]` = index in backbone list of the Cα of atom i.
///
/// If a residue has no backbone atom, the previous backbone atom is used as fallback.
pub fn build_atom_to_backbone(
    bb_res_keys:  &[(i32, char)],   // (res_serial, chain) for each backbone atom
    all_res_keys: &[(i32, char)],   // same for every heavy atom in order
) -> Vec<usize> {
    use std::collections::HashMap;
    let bb_map: HashMap<(i32, char), usize> = bb_res_keys
        .iter().enumerate().map(|(i, &k)| (k, i)).collect();
    let mut last = 0usize;
    all_res_keys.iter().map(|key| {
        if let Some(&idx) = bb_map.get(key) { last = idx; }
        last
    }).collect()
}
