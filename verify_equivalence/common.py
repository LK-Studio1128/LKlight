"""
common.py — shared utilities for the LKlight / LightDock equivalence benchmark.

Pure-python PDB parsing, GSO output parsing, CAPRI-class assessment and
top-N success-rate evaluation.  NumPy is used when available (Kabsch
receptor alignment); without it the evaluator falls back to a no-alignment
RMSD approximation (still consistent across engines, which is what the
equivalence test needs).
"""

import json
import math
import os
import re

try:
    import numpy as np
    HAVE_NUMPY = True
except Exception:  # pragma: no cover
    HAVE_NUMPY = False


# ─────────────────────────────────────────────────────────────────────────────
# PDB parsing
# ─────────────────────────────────────────────────────────────────────────────

def parse_pdb(path):
    """Return a list of atom dicts: chain, segid, resnum, resname, atomname, x, y, z."""
    atoms = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                if len(line) < 54:
                    continue
                try:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                except ValueError:
                    continue
                atoms.append({
                    "chain": (line[21:22] or "A").strip() or "A",
                    "segid": line[72:76].strip(),
                    "resnum": _resnum(line[22:26]),
                    "resname": line[17:20].strip(),
                    "atomname": line[12:16].strip(),
                    "x": x, "y": y, "z": z,
                })
    return atoms


def _resnum(s):
    try:
        return int(float(s.strip()))
    except ValueError:
        return -1


def read_setup(path):
    """Read a LightDock setup.json (same schema for the Python and Rust engines)."""
    with open(path) as fh:
        return json.load(fh)


def receptor_fingerprint(receptor_pdb):
    """{(chain, resnum, atomname)} set identifying receptor atoms by identity."""
    fp = set()
    for a in parse_pdb(receptor_pdb):
        fp.add((a["chain"], a["resnum"], a["atomname"]))
    return fp


def split_model(atoms, rec_fp):
    """Split a merged model PDB into receptor / ligand atom lists via fingerprint."""
    rec, lig = [], []
    for a in atoms:
        if (a["chain"], a["resnum"], a["atomname"]) in rec_fp:
            rec.append(a)
        else:
            lig.append(a)
    return rec, lig


def _fingerprint_hits(atoms, fp):
    return sum(1 for a in atoms if (a["chain"], a["resnum"], a["atomname"]) in fp)


def match_native_chains(native_atoms, rec_fp, lig_fp):
    """
    Guess which native chain(s) are receptor vs ligand by fingerprint hit ratio.
    Returns (rec_chains, lig_chains) as sets of (chain, segid) tuples.
    """
    groups = {}
    for a in native_atoms:
        key = (a["chain"], a["segid"])
        groups.setdefault(key, []).append(a)

    scored = []
    for key, grp in groups.items():
        if not grp:
            continue
        hr = _fingerprint_hits(grp, rec_fp) / len(grp)
        hl = _fingerprint_hits(grp, lig_fp) / len(grp)
        scored.append((key, hr, hl, len(grp)))

    scored.sort(key=lambda t: -max(t[1], t[2]))
    rec_keys, lig_keys = set(), set()
    rec_seen, lig_seen = 0, 0
    for key, hr, hl, n in scored:
        if hr >= hl and (rec_seen == 0 or hr > 0.5):
            rec_keys.add(key)
            rec_seen += n
        elif hl > hr and (lig_seen == 0 or hl > 0.5):
            lig_keys.add(key)
            lig_seen += n
        else:
            # ambiguous: assign to whichever side has seen less mass
            (rec_keys if rec_seen <= lig_seen else lig_keys).add(key)
            if rec_seen <= lig_seen:
                rec_seen += n
            else:
                lig_seen += n
    return rec_keys, lig_keys


def atoms_in_chains(atoms, keys):
    return [a for a in atoms if (a["chain"], a["segid"]) in keys]


# ─────────────────────────────────────────────────────────────────────────────
# GSO output parsing (mirrors LKlight's parse_gso_file)
# ─────────────────────────────────────────────────────────────────────────────

def parse_gso(filename):
    """
    Return [(glowworm_id, luciferin, scoring), ...] in file order.
    Line format: "(tx, ty, tz, q0, q1, q2, q3[, ext...])  <rest>"
    where rest[2] = luciferin and rest[5] = scoring (0-indexed, whitespace split).
    """
    out = []
    gid = 0
    with open(filename) as fh:
        for raw in fh:
            s = raw.strip()
            if not s or s.startswith("#"):
                gid += 1
                continue
            pend = s.find(")")
            if pend < 0:
                gid += 1
                continue
            rest = s[pend + 1:].split()
            luciferin = float(rest[2]) if len(rest) > 2 else 0.0
            scoring = float(rest[5]) if len(rest) > 5 else 0.0
            out.append((gid, luciferin, scoring))
            gid += 1
    return out


# ─────────────────────────────────────────────────────────────────────────────
# Geometry helpers
# ─────────────────────────────────────────────────────────────────────────────

def _coords(atoms):
    return [(a["x"], a["y"], a["z"]) for a in atoms]


def _pair_match(a_list, b_list):
    """Match atoms by (chain, resnum, atomname); returns (A_idx, B_idx) pairs."""
    bmap = {}
    for j, b in enumerate(b_list):
        bmap.setdefault((b["chain"], b["resnum"], b["atomname"]), []).append(j)
    pairs = []
    for i, a in enumerate(a_list):
        cand = bmap.get((a["chain"], a["resnum"], a["atomname"]), [])
        if cand:
            pairs.append((i, cand[0]))
    return pairs


def _mean(coords):
    n = len(coords)
    if n == 0:
        return (0.0, 0.0, 0.0)
    return (sum(c[0] for c in coords) / n,
            sum(c[1] for c in coords) / n,
            sum(c[2] for c in coords) / n)


def _rmsd_pairs(a_list, b_list, pairs):
    if not pairs:
        return float("inf")
    s = 0.0
    for i, j in pairs:
        d0 = a_list[i][0] - b_list[j][0]
        d1 = a_list[i][1] - b_list[j][1]
        d2 = a_list[i][2] - b_list[j][2]
        s += d0 * d0 + d1 * d1 + d2 * d2
    return math.sqrt(s / len(pairs))


def _kabsch_rot(p, q):
    """Rotation matrix aligning P onto Q (P, Q: n x 3 lists). Returns 3x3 list."""
    n = len(p)
    pc = _mean(p)
    qc = _mean(q)
    p0 = [(x - pc[0], y - pc[1], z - pc[2]) for (x, y, z) in p]
    q0 = [(x - qc[0], y - qc[1], z - qc[2]) for (x, y, z) in q]
    h = [[0.0, 0.0, 0.0] for _ in range(3)]
    for i in range(n):
        for r in range(3):
            for c in range(3):
                h[r][c] += p0[i][r] * q0[i][c]
    import numpy as _np
    u, _s, vt = _np.linalg.svd(_np.array(h))
    d = _np.identity(3)
    d[2][2] = _np.linalg.det(_np.dot(vt.T, u.T))
    rot = _np.dot(vt.T, _np.dot(d, u.T))
    return rot, pc, qc


def _apply_rot(points, rot):
    return [(
        rot[0][0] * x + rot[0][1] * y + rot[0][2] * z,
        rot[1][0] * x + rot[1][1] * y + rot[1][2] * z,
        rot[2][0] * x + rot[2][1] * y + rot[2][2] * z,
    ) for (x, y, z) in points]


# ─────────────────────────────────────────────────────────────────────────────
# CAPRI assessment
# ─────────────────────────────────────────────────────────────────────────────

CONTACT_CUTOFF = 5.0  # Angstrom, interface contact definition


def interface_contacts(rec_atoms, lig_atoms, cutoff=CONTACT_CUTOFF):
    """Set of (rec_idx, lig_idx) atom pairs within cutoff (numpy-accelerated)."""
    c2 = cutoff * cutoff
    if HAVE_NUMPY and len(rec_atoms) * len(lig_atoms) > 1e4:
        import numpy as _np
        R = _np.array([[a["x"], a["y"], a["z"]] for a in rec_atoms], dtype=_np.float64)
        L = _np.array([[a["x"], a["y"], a["z"]] for a in lig_atoms], dtype=_np.float64)
        contacts = set()
        step = 512
        for i0 in range(0, len(R), step):
            block = R[i0:i0 + step]
            d2 = _np.sum((block[:, None, :] - L[None, :, :]) ** 2, axis=2)
            ii, jj = _np.nonzero(d2 <= c2)
            for a, b in zip(ii + i0, jj):
                contacts.add((int(a), int(b)))
        return contacts
    contacts = set()
    for i, ra in enumerate(rec_atoms):
        for j, la in enumerate(lig_atoms):
            dx = ra["x"] - la["x"]
            dy = ra["y"] - la["y"]
            dz = ra["z"] - la["z"]
            if dx * dx + dy * dy + dz * dz <= c2:
                contacts.add((i, j))
    return contacts


def capri_assessment(model_atoms, native_atoms, rec_fp, lig_fp):
    """
    Evaluate one model against the native complex.
    Returns dict: fnat, lrmsd, rec_rmsd, capri_class (-1 incorrect … 3 high).
    """
    m_rec, m_lig = split_model(model_atoms, rec_fp)
    n_rec_keys, n_lig_keys = match_native_chains(native_atoms, rec_fp, lig_fp)
    n_rec = atoms_in_chains(native_atoms, n_rec_keys)
    n_lig = atoms_in_chains(native_atoms, n_lig_keys)

    if not m_rec or not m_lig or not n_rec or not n_lig:
        return {"fnat": 0.0, "lrmsd": float("inf"), "rec_rmsd": float("inf"),
                "capri_class": -1, "error": "chain mapping failed"}

    # L-RMSD: superpose the LIGAND (CAPRI definition), then RMSD.
    # Note: model conformations live in the LightDock reference-point frame,
    # so no frame pre-alignment is needed — ligand superposition handles it.
    lrmsd = float("inf")
    lp = _pair_match(m_lig, n_lig)
    lig_model = _coords(m_lig)
    lig_native = _coords(n_lig)
    if HAVE_NUMPY and len(lp) >= 3:
        p = [lig_model[i] for i, _j in lp]
        q = [lig_native[j] for _i, j in lp]
        rot, pc, qc = _kabsch_rot(p, q)
        aligned = []
        for (x, y, z) in lig_model:
            xx, yy, zz = x - pc[0], y - pc[1], z - pc[2]
            aligned.append((rot[0][0] * xx + rot[0][1] * yy + rot[0][2] * zz + qc[0],
                            rot[1][0] * xx + rot[1][1] * yy + rot[1][2] * zz + qc[1],
                            rot[2][0] * xx + rot[2][1] * yy + rot[2][2] * zz + qc[2]))
        lrmsd = _rmsd_pairs(aligned, lig_native, lp)
    elif lp:
        lrmsd = _rmsd_pairs(lig_model, lig_native, lp)

    # fraction of native contacts (contact pair sets are frame-invariant)
    c_model = interface_contacts(m_rec, m_lig)
    c_native = interface_contacts(n_rec, n_lig)
    fnat = 0.0
    if c_native:
        def res_pairs(contacts, rec, lig):
            return set(
                (rec[i]["resnum"], lig[j]["resnum"])
                for (i, j) in contacts
            )
        pm = res_pairs(c_model, m_rec, m_lig)
        pn = res_pairs(c_native, n_rec, n_lig)
        if pn:
            fnat = len(pm & pn) / len(pn)

    # receptor alignment quality (diagnostic; should be ~0 for rigid docking)
    rp = _pair_match(m_rec, n_rec)
    rec_rmsd = _rmsd_pairs(_coords(m_rec), _coords(n_rec), rp)

    capri = -1
    if fnat >= 0.5 and lrmsd <= 1.0:
        capri = 3
    elif fnat >= 0.3 and lrmsd <= 5.0:
        capri = 2
    elif fnat >= 0.1 and lrmsd <= 10.0:
        capri = 1
    return {"fnat": fnat, "lrmsd": lrmsd, "rec_rmsd": rec_rmsd,
            "capri_class": capri}


# ─────────────────────────────────────────────────────────────────────────────
# Top-N success-rate evaluation
# ─────────────────────────────────────────────────────────────────────────────

def evaluate_engine(case_dir, rec_pdb, lig_pdb, native_pdb, steps, swarms,
                    top_cutoffs=(1, 5, 10, 20, 50, 100)):
    """
    Run the shared evaluation for one engine's output in case_dir.
    Reads swarm_*/gso_{steps}.out, ranks by scoring (desc), picks top-N
    conformations (lightdock_{gid}.pdb) and scores them with CAPRI.

    Returns dict with per-model metrics and success rates.
    """
    rec_fp = receptor_fingerprint(rec_pdb)
    native_atoms = parse_pdb(native_pdb)

    all_sol = []  # (scoring, swarm, glowworm_id)
    for s in range(swarms):
        gso = os.path.join(case_dir, f"swarm_{s}", f"gso_{steps}.out")
        if not os.path.exists(gso):
            # engines save every 10 steps (and step 1); fall back to the last snapshot
            import glob
            hits = glob.glob(os.path.join(case_dir, f"swarm_{s}", "gso_*.out"))
            if hits:
                gso = max(hits, key=lambda p: int(os.path.basename(p)[4:-4]))
            else:
                continue
        for (gid, _luc, scoring) in parse_gso(gso):
            all_sol.append((scoring, s, gid))

    all_sol.sort(key=lambda t: -t[0])  # higher score first (LightDock convention)

    # only fully assess the top-N solutions that the success-rate cutoffs need;
    # the rest are counted but not scored (keeps big complexes tractable)
    assess_n = max(top_cutoffs)
    models = []
    lig_fp = receptor_fingerprint(lig_pdb)
    for rank, (scoring, s, gid) in enumerate(all_sol):
        if rank >= assess_n:
            break
        pdb_path = os.path.join(case_dir, f"swarm_{s}", f"lightdock_{gid}.pdb")
        if not os.path.exists(pdb_path):
            continue
        model_atoms = parse_pdb(pdb_path)
        res = capri_assessment(model_atoms, native_atoms, rec_fp, lig_fp)
        models.append({"rank": rank + 1, "swarm": s, "glowworm": gid,
                       "scoring": scoring, **res})

    success = {}
    for n in top_cutoffs:
        topn = [m for m in models if m["rank"] <= n]
        ok = sum(1 for m in topn if m["capri_class"] >= 1)  # acceptable or better
        den = max(1, min(n, len(topn)))
        success[n] = ok / den

    return {"success_rates": success, "models": models, "n_solutions": len(all_sol)}


def parse_official_list(list_path, top_cutoffs=(1, 5, 10, 20, 50, 100)):
    """
    Parse a LightDock official results .list file (columns:
    name iRMSD ligandRMSD Fnat score).  Used as the Python reference
    baseline without re-running the Python engine.

    Returns a dict compatible with evaluate_engine() (success_rates,
    n_solutions, models), so both engines produce identical per-model
    records for the summary tables and plots.
    """
    rows = []  # (score, irmsd, lrmsd, fnat, name)
    with open(list_path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.split()
            if len(f) < 5:
                continue
            try:
                rows.append((float(f[4]), float(f[1]), float(f[2]), float(f[3]), f[0]))
            except ValueError:
                continue
    rows.sort(key=lambda t: -t[0])

    def capri_class(fnat, lrmsd):
        if fnat >= 0.5 and lrmsd <= 1.0:
            return 3
        if fnat >= 0.3 and lrmsd <= 5.0:
            return 2
        if fnat >= 0.1 and lrmsd <= 10.0:
            return 1
        return -1

    # per-model records (rank, swarm, glowworm parsed from "swarm_N_M.pdb")
    models = []
    for rank, (_score, irmsd, lrmsd, fnat, name) in enumerate(rows):
        if rank >= max(top_cutoffs):
            break
        swarm = glowworm = -1
        base = os.path.basename(name)
        if base.startswith("swarm_"):
            parts = base[6:].split("_")
            try:
                swarm = int(parts[0])
                glowworm = int(parts[1].split(".")[0])
            except (ValueError, IndexError):
                pass
        models.append({"rank": rank + 1, "swarm": swarm, "glowworm": glowworm,
                       "scoring": float(rows[rank][0]),
                       "fnat": fnat, "lrmsd": lrmsd,
                       "rec_rmsd": irmsd, "capri_class": capri_class(fnat, lrmsd)})

    success = {}
    for n in top_cutoffs:
        topn = rows[:n]
        ok = sum(1 for _s, _i, lrmsd, fnat, _nm in topn
                 if capri_class(fnat, lrmsd) >= 1)
        success[n] = ok / n if n else 0.0
    return {"success_rates": success, "n_solutions": len(rows), "models": models}

