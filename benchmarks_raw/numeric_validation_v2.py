#!/usr/bin/env python3
"""Numeric validation v2 — pose-robust, pure-comparison design.

Applies rigid transforms to the ligand PDB *file* first (numpy), then
evaluates both engines at identity. Both engines see identical inputs ->
any diff is a real engine diff.

Poses: identity / translate(5,-3,2.5) / rotate60(axis=(1,1,0)) / transrot.
13 methods x 4 poses = 52 items (pydock python-side unavailable -> noted).
"""
import subprocess, math, os, sys

LKLIGHT = "/Users/luoxiaowen/Desktop/LKDock/byi/LKlight/target/release/LKlight"
OUT = "/Users/luoxiaowen/Desktop/LKDock/LKlight论文/benchmarks_raw/numeric_validation_raw.tsv"
WORK = "/tmp/numval"

# ── numpy/scipy/lightdock paths MUST be set before any numpy import ──────────
sys.path.insert(1, "/tmp/numval/deps_n")  # numpy 1.23.5 (cp39 arm64 wheel)
sys.path.append("/tmp/numval/deps")       # scipy 1.13.1
sys.path.append("/Users/luoxiaowen/mambaforge/envs/docking_env/lib/python3.9/site-packages")

import numpy as np

def load_atoms(path):
    atoms = []
    for line in open(path):
        if line.startswith(("ATOM", "HETATM")):
            atoms.append([
                float(line[30:38]), float(line[38:46]), float(line[46:54]),
            ])
    return np.array(atoms)

def write_transformed(src, dst, R=None, t=None):
    out = []
    for line in open(src):
        if line.startswith(("ATOM", "HETATM")):
            xyz = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
            if R is not None:
                xyz = R @ xyz
            if t is not None:
                xyz = xyz + t
            out.append(line[:30] + "%8.3f%8.3f%8.3f" % tuple(xyz) + line[54:])
        else:
            out.append(line)
    open(dst, "w").writelines(out)

def rot_matrix(axis, deg):
    axis = np.asarray(axis, dtype=float)
    axis = axis / np.linalg.norm(axis)
    a = math.radians(deg)
    c, s = math.cos(a), math.sin(a)
    K = np.array([[0, -axis[2], axis[1]],
                  [axis[2], 0, -axis[0]],
                  [-axis[1], axis[0], 0]])
    return np.eye(3) + s * K + (1 - c) * (K @ K)

SYSTEMS = {
    "2oob": ("/tmp/numval/2oob_receptor_noh.pdb", "/tmp/numval/2oob_ligand_noh.pdb"),
    "1azp": ("/tmp/numval/1azp_receptor_noh.pdb", "/tmp/numval/1azp_ligand_noh.pdb"),
}
PROTEIN_METHODS = ["dfire", "fastdfire", "dfire2", "mj3h", "cpydock",
                   "sd", "vdw", "pisa", "sipper", "tobi"]
DNA_METHODS = ["dna", "ddna"]

from lightdock.pdbutil.PDBIO import parse_complex_from_file
from lightdock.structure.complex import Complex
import importlib

_cache = {}
def get_parsed(path):
    # Pose files are rewritten per (pose, method): a path-based cache would
    # serve a stale file (e.g. 2oob ligand) to a later method (dna/1azp).
    # Only cache the fixed reference PDBs, never the pose_*_ligand.pdb files.
    if "pose_" in os.path.basename(path):
        atoms, residues, chains = parse_complex_from_file(path)
        return Complex(chains, atoms, structure_file_name=path)
    if path not in _cache:
        atoms, residues, chains = parse_complex_from_file(path)
        _cache[path] = Complex(chains, atoms, structure_file_name=path)
    return _cache[path]

def py_score(method, rec_path, lig_path):
    receptor = get_parsed(rec_path).clone()
    ligand = get_parsed(lig_path).clone()
    module = importlib.import_module("lightdock.scoring.%s.driver" % method)
    ScoringFunction = getattr(module, "DefinedScoringFunction")
    ModelAdapter = getattr(module, "DefinedModelAdapter")
    adapter = ModelAdapter(receptor, ligand)
    sf = ScoringFunction()
    return float(sf(
        adapter.receptor_model,
        adapter.receptor_model.coordinates[0],
        adapter.ligand_model,
        adapter.ligand_model.coordinates[0],
    ))

def rust_score(method, rec_path, lig_path):
    r = subprocess.run(
        [LKLIGHT, "score", rec_path, lig_path, method,
         "--tx", "0", "--ty", "0", "--tz", "0",
         "--qw", "1", "--qx", "0", "--qy", "0", "--qz", "0"],
        capture_output=True, text=True, timeout=300)
    if r.returncode != 0:
        return None, "rc=%d: %s" % (r.returncode, r.stderr.strip()[:150])
    for tok in reversed(r.stdout.replace(":", " ").replace(",", " ").split()):
        try:
            return float(tok), ""
        except ValueError:
            continue
    return None, "unparseable stdout"

def main():
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w") as f:
        f.write("pose\tmethod\tsystem\tpython_score\trust_score\tdiff_abs\tstatus\n")

        A = {}
        for m in PROTEIN_METHODS:
            A[m] = SYSTEMS["2oob"]
        for m in DNA_METHODS:
            A[m] = SYSTEMS["1azp"]

        poses = [
            ("identity", None, None),
            ("translate", None, np.array([5.0, -3.0, 2.5])),
            ("rotate", rot_matrix((1, 1, 0), 60), None),
            ("transrot", rot_matrix((1, 1, 0), 60), np.array([5.0, -3.0, 2.5])),
        ]

        tol = 1e-6
        total_ok = 0
        total_items = 0

        for pose_name, R, t in poses:
            for method, (recf, ligf) in A.items():
                lig_x = os.path.join(WORK, "pose_%s_ligand.pdb" % pose_name)
                write_transformed(ligf, lig_x, R=R, t=t)
                label = "%s/%s" % (pose_name, method)
                # python
                try:
                    ps = py_score(method, recf, lig_x)
                    perr = ""
                except Exception as e:
                    ps = None
                    perr = "%s: %s" % (type(e).__name__, str(e)[:80])
                # rust
                rs, rerr = rust_score(method, recf, lig_x)
                if ps is not None and rs is not None:
                    diff = abs(ps - rs)
                    status = "OK" if diff < tol else "MISMATCH"
                    if status == "OK":
                        total_ok += 1
                    total_items += 1
                    line = "%s\t%s\t%s\t%.10f\t%.6f\t%.3g\t%s\n" % (
                        pose_name, method, "2oob" if method in PROTEIN_METHODS else "1azp",
                        ps, rs, diff, status)
                else:
                    line = "%s\t%s\t%s\t%s\t%s\t\t%s\n" % (
                        pose_name, method, "2oob" if method in PROTEIN_METHODS else "1azp",
                        ps if ps is not None else "ERR:%s" % perr,
                        rs if rs is not None else "ERR:%s" % rerr, "ERROR")
                f.write(line)
                print("%-44s %s" % (label, line.strip().split("\t")[-1]), flush=True)

        f.write("# summary\t%d/%d OK (tol=%.0e)\n" % (total_ok, total_items, tol))

    print("\n===== %d/%d OK (tol=%.0e) =====" % (total_ok, total_items, tol))

if __name__ == "__main__":
    main()
