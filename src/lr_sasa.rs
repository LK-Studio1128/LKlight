/// lr_sasa.rs — Lee & Richards solvent-accessible surface area.
///
/// Faithful port of freesasa 2.2 (src/nb.c + src/sasa_lr.c) so that per-atom
/// areas match the Python reference (freesasa.calc) to the last bit:
///   - cell-list neighbour construction with the exact same traversal order,
///   - the same slice/arc arithmetic, insertion sort and summation order.
///
/// Algorithm (Lee & Richards, 1971): for each atom, sweep n_slices horizontal
/// planes through its sphere (radius = vdw + probe); on each plane the atom's
/// circle is partially buried by neighbour circles; the exposed arc length is
/// summed and multiplied by delta*Ri to give the SASA contribution.

const TWOPI: f64 = 2.0 * std::f64::consts::PI;

// ── neighbour lists (port of nb.c) ───────────────────────────────────────────

struct Nb {
    nn: Vec<usize>,          // number of neighbours per atom
    nb: Vec<Vec<usize>>,     // neighbour indices, in discovery order
    xyd: Vec<Vec<f64>>,      // sqrt(dx^2+dy^2) per neighbour
    xd: Vec<Vec<f64>>,       // dx = xj - xi (as seen from atom i)
    yd: Vec<Vec<f64>>,       // dy = yj - yi
}

/// Build the neighbour list exactly like freesasa_nb_new()/nb.c.
/// radii must already include the probe radius (R = atom_radius + probe).
fn build_neighbours(coords: &[[f64; 3]], radii: &[f64]) -> Nb {
    let n = coords.len();
    let mut nn = vec![0usize; n];
    let mut nb: Vec<Vec<usize>> = (0..n).map(|_| Vec::with_capacity(128)).collect();
    let mut xyd: Vec<Vec<f64>> = (0..n).map(|_| Vec::with_capacity(128)).collect();
    let mut xd: Vec<Vec<f64>> = (0..n).map(|_| Vec::with_capacity(128)).collect();
    let mut yd: Vec<Vec<f64>> = (0..n).map(|_| Vec::with_capacity(128)).collect();

    // cell size = 2 * max radius (with probe)
    let rmax = radii.iter().cloned().fold(0.0f64, f64::max);
    let d = 2.0 * rmax;
    assert!(d > 0.0);

    // bounds (cell_list_bounds)
    let mut xmin = coords[0][0];
    let mut xmax = coords[0][0];
    let mut ymin = coords[0][1];
    let mut ymax = coords[0][1];
    let mut zmin = coords[0][2];
    let mut zmax = coords[0][2];
    for c in &coords[1..] {
        xmin = xmin.min(c[0]); xmax = xmax.max(c[0]);
        ymin = ymin.min(c[1]); ymax = ymax.max(c[1]);
        zmin = zmin.min(c[2]); zmax = zmax.max(c[2]);
    }
    xmin -= d / 2.0; xmax += d / 2.0;
    ymin -= d / 2.0; ymax += d / 2.0;
    zmin -= d / 2.0; zmax += d / 2.0;
    let nx = ((xmax - xmin) / d).ceil() as usize;
    let ny = ((ymax - ymin) / d).ceil() as usize;
    let nz = ((zmax - zmin) / d).ceil() as usize;
    let n_cells = nx * ny * nz;

    // cells (fill_cells): atoms in increasing index order per cell
    let mut cell_atoms: Vec<Vec<usize>> = (0..n_cells).map(|_| Vec::new()).collect();
    let cell_of = |x: f64, y: f64, z: f64| -> usize {
        let ix = (((x - xmin) / d) as isize).clamp(0, nx as isize - 1) as usize;
        let iy = (((y - ymin) / d) as isize).clamp(0, ny as isize - 1) as usize;
        let iz = (((z - zmin) / d) as isize).clamp(0, nz as isize - 1) as usize;
        ix + nx * (iy + ny * iz)
    };
    for i in 0..n {
        let c = coords[i];
        cell_atoms[cell_of(c[0], c[1], c[2])].push(i);
    }

    // forward-neighbour cells of each cell (fill_nb): scalar product with
    // (1,1,1) must be non-negative; loops in the same (ix,iy,iz) order as C.
    let mut cell_nb: Vec<Vec<usize>> = (0..n_cells).map(|_| Vec::new()).collect();
    for icx in 0..nx {
        for icy in 0..ny {
            for icz in 0..nz {
                let ic = icx + nx * (icy + ny * icz);
                let xmin_i = if icx > 0 { icx - 1 } else { 0 };
                let xmax_i = if icx < nx - 1 { icx + 1 } else { icx };
                let ymin_i = if icy > 0 { icy - 1 } else { 0 };
                let ymax_i = if icy < ny - 1 { icy + 1 } else { icy };
                let zmin_i = if icz > 0 { icz - 1 } else { 0 };
                let zmax_i = if icz < nz - 1 { icz + 1 } else { icz };
                for i in xmin_i..=xmax_i {
                    for j in ymin_i..=ymax_i {
                        for k in zmin_i..=zmax_i {
                            if (i as isize - icx as isize)
                                + (j as isize - icy as isize)
                                + (k as isize - icz as isize)
                                >= 0
                            {
                                cell_nb[ic].push(i + nx * (j + ny * k));
                            }
                        }
                    }
                }
            }
        }
    }

    // register pairs (nb_fill_list + nb_calc_cell_pair + nb_add_pair)
    for ic in 0..n_cells {
        for &jc in &cell_nb[ic] {
            let a_i = &cell_atoms[ic];
            let a_j = &cell_atoms[jc];
            for (ii, &ia) in a_i.iter().enumerate() {
                let ri = radii[ia];
                let xi = coords[ia][0];
                let yi = coords[ia][1];
                let zi = coords[ia][2];
                let start = if ic == jc { ii + 1 } else { 0 };
                for &ja in &a_j[start..] {
                    let rj = radii[ja];
                    let xj = coords[ja][0];
                    let yj = coords[ja][1];
                    let zj = coords[ja][2];
                    let cut2 = (ri + rj) * (ri + rj);
                    let dx = xj - xi;
                    let dy = yj - yi;
                    let dz = zj - zi;
                    if dx * dx + dy * dy + dz * dz < cut2 {
                        // nb_add_pair
                        let nni = nn[ia];
                        let nnj = nn[ja];
                        nb[ia].push(ja);
                        nb[ja].push(ia);
                        let dxy = (dx * dx + dy * dy).sqrt();
                        xyd[ia].push(dxy);
                        xyd[ja].push(dxy);
                        xd[ia].push(dx);
                        xd[ja].push(-dx);
                        yd[ia].push(dy);
                        yd[ja].push(-dy);
                        nn[ia] = nni + 1;
                        nn[ja] = nnj + 1;
                    }
                }
            }
        }
    }

    Nb { nn, nb, xyd, xd, yd }
}

// ── SASA per atom (port of sasa_lr.c atom_area) ─────────────────────────────

/// Insertion sort of the arc intervals by start point (sort_arcs in C).
fn sort_arcs(arc: &mut [f64], n: usize) {
    let mut i = 2;
    while i < 2 * n {
        let mut tmp = [arc[i], arc[i + 1]];
        let mut j = i;
        while j > 0 && arc[j - 2] > tmp[0] {
            arc[j] = arc[j - 2];
            arc[j + 1] = arc[j - 1];
            j -= 2;
        }
        arc[j] = tmp[0];
        arc[j + 1] = tmp[1];
        i += 2;
    }
}

/// Sum of exposed arcs given buried intervals (exposed_arc_length in C).
fn exposed_arc_length(arc: &mut [f64], n: usize) -> f64 {
    if n == 0 {
        return TWOPI;
    }
    sort_arcs(arc, n);
    let mut sum = arc[0];
    let mut sup = arc[1];
    let mut i2 = 2;
    while i2 < 2 * n {
        if sup < arc[i2] {
            sum += arc[i2] - sup;
        }
        let tmp = arc[i2 + 1];
        if tmp > sup {
            sup = tmp;
        }
        i2 += 2;
    }
    sum + TWOPI - sup
}

/// Area of atom i (atom_area in C). `nb` contains the pre-built neighbour
/// lists; `coords`/`radii` (with probe) are the same arrays used to build it.
fn atom_area(i: usize, coords: &[[f64; 3]], radii: &[f64], nb: &Nb, nslices: usize) -> f64 {
    let nni = nb.nn[i];
    let zi = coords[i][2];
    let ri = radii[i];

    // scratch arrays for neighbours (z and R)
    let mut z_nb = vec![0.0f64; nni];
    let mut r_nb = vec![0.0f64; nni];
    for (j, &nbi) in nb.nb[i].iter().enumerate() {
        z_nb[j] = coords[nbi][2];
        r_nb[j] = radii[nbi];
    }

    let delta = 2.0 * ri / nslices as f64;
    let mut z = zi - ri - 0.5 * delta;
    let mut sasa = 0.0f64;
    // worst case: each neighbour contributes at most 2 arcs (4 doubles)
    let mut arc = vec![0.0f64; 4 * nni.max(1)];

    for _islice in 0..nslices {
        z += delta;
        let di = (zi - z).abs();
        let ri_p2 = ri * ri - di * di;
        if ri_p2 < 0.0 {
            continue; // round-off
        }
        let ri_p = ri_p2.sqrt();
        if ri_p <= 0.0 {
            continue; // more round-off
        }
        let mut n_arcs = 0usize;
        let mut is_buried = false;
        for j in 0..nni {
            let zj = z_nb[j];
            let dj = (zj - z).abs();
            let rj = r_nb[j];
            if dj < rj {
                let rj_p2 = rj * rj - dj * dj;
                let rj_p = rj_p2.sqrt();
                let dij = nb.xyd[i][j];
                if dij >= ri_p + rj_p {
                    continue; // not in contact
                }
                if dij + ri_p < rj_p {
                    is_buried = true; // circle i completely inside j
                    break;
                }
                if dij + rj_p < ri_p {
                    continue; // circle j completely inside i
                }
                let alpha = ((ri_p2 + dij * dij - rj_p2) / (2.0 * ri_p * dij)).acos();
                let beta = nb.yd[i][j].atan2(nb.xd[i][j]) + std::f64::consts::PI;
                let mut inf = beta - alpha;
                let mut sup = beta + alpha;
                if inf < 0.0 {
                    inf += TWOPI;
                }
                if sup > 2.0 * std::f64::consts::PI {
                    sup -= TWOPI;
                }
                let narc2 = 2 * n_arcs;
                if sup < inf {
                    arc[narc2] = 0.0;
                    arc[narc2 + 1] = sup;
                    arc[narc2 + 2] = inf;
                    arc[narc2 + 3] = TWOPI;
                    n_arcs += 2;
                } else {
                    arc[narc2] = inf;
                    arc[narc2 + 1] = sup;
                    n_arcs += 1;
                }
            }
        }
        if !is_buried {
            sasa += delta * ri * exposed_arc_length(&mut arc[..], n_arcs);
        }
    }
    sasa
}

/// Compute per-atom solvent-accessible surface area with the Lee & Richards
/// algorithm (single-threaded, deterministic, bit-identical to freesasa 2.2
/// with default parameters: probe = 1.4, n_slices = 20).
///
/// `atom_radii` are the bare van-der-Waals/solvation radii; the probe is added
/// internally, matching freesasa's `init_lr`.
pub fn lee_richards_sasa(coords: &[[f64; 3]], atom_radii: &[f64]) -> Vec<f64> {
    let n = coords.len();
    let probe = 1.4;
    let nslices = 20;
    if n == 0 {
        return Vec::new();
    }
    let radii: Vec<f64> = atom_radii.iter().map(|r| r + probe).collect();
    let nb = build_neighbours(coords, &radii);
    let mut sasa = vec![0.0f64; n];
    for i in 0..n {
        sasa[i] = atom_area(i, coords, &radii, &nb, nslices);
    }
    sasa
}
