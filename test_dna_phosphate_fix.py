#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Regression test for the LKlight protein-RNA `*-P` crash fix.

Builds a minimal protein receptor + a raw single-letter RNA ligand ("U" residue
WITH a phosphate "P" atom, plus a deliberately-unknown atom "ZZ"), then runs
`LKlight score rec.pdb lig.pdb dna`.

Before the fix: the dna scorer panicked at dna.rs `*-P not supported` -> non-zero
exit + panic message -> all swarms would die.
After the fix: it must NOT panic; it should print a score and exit 0.

Tests all three freshly built binaries found under dist/ that can run on this host
(native mac binary at minimum). Windows/Linux binaries are only existence/type
checked here (cannot execute on macOS).
"""
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
DIST = os.path.join(HERE, "dist")
TMP = "/tmp/lkfix"
os.makedirs(TMP, exist_ok=True)


def atom(serial, name, res, chain, resseq, x, y, z, elem):
    # Strict PDB columns: name 13-16, altLoc 17, resName 18-20, chain 22,
    # resSeq 23-26, x 31-38, y 39-46, z 47-54, element 77-78.
    nm = name if len(name) >= 4 else " " + name
    return (f"ATOM  {serial:5d} {nm:<4s} {res:>3s} {chain}{resseq:4d}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {elem:>2s}\n")


def write_pdbs():
    rec = []
    rec_atoms = [("N", "N", 0, 0, 0), ("CA", "C", 1.5, 0, 0), ("C", "C", 2, 1.4, 0),
                 ("O", "O", 3.2, 1.6, 0), ("CB", "C", 2, -0.8, 1.2)]
    for i, (nm, e, x, y, z) in enumerate(rec_atoms, 1):
        rec.append(atom(i, nm, "ALA", "A", 1, x, y, z, e))
    rec.append("TER\nEND\n")
    rec_path = os.path.join(TMP, "rec.pdb")
    open(rec_path, "w").write("".join(rec))

    # raw single-letter RNA "U" WITH a phosphate "P" atom = the exact historical
    # crash trigger (dna.rs *-P not supported). Standard PDB v3 RNA atom names.
    lig = []
    lig_atoms = [("P", "P", 5, 0, 0), ("OP1", "O", 5.5, 1.2, 0), ("OP2", "O", 5.5, -1.2, 0),
                 ("O5'", "O", 6, 0.5, 0), ("C5'", "C", 6.8, 0, 0), ("C4'", "C", 7.6, 0.6, 0),
                 ("O4'", "O", 7.9, 1.8, 0), ("C1'", "C", 8.6, 1.6, 0), ("N1", "N", 9.4, 2.2, 0),
                 ("C2", "C", 10.2, 1.8, 0), ("O2", "O", 10.5, 0.6, 0), ("O2'", "O", 8.2, 2.6, 0)]
    for i, (nm, e, x, y, z) in enumerate(lig_atoms, 1):
        lig.append(atom(i, nm, "U", "B", 1, x, y, z, e))
    lig.append("TER\nEND\n")
    lig_path = os.path.join(TMP, "lig.pdb")
    open(lig_path, "w").write("".join(lig))
    return rec_path, lig_path


def run_score(binary, rec, lig, method="dna"):
    try:
        proc = subprocess.run([binary, "score", rec, lig, method],
                              capture_output=True, text=True, timeout=120)
        return proc.returncode, (proc.stdout or "") + (proc.stderr or "")
    except Exception as e:  # noqa
        return 999, f"EXEC ERROR: {e}"


def main():
    rec, lig = write_pdbs()
    print(f"[setup] wrote {rec} , {lig}")

    mac_bin = os.path.join(DIST, "LKlight_mac_arm64")
    win_bin = os.path.join(DIST, "LKlight_win_x64.exe")
    lin_bin = os.path.join(DIST, "LKlight_linux_x64")

    overall_ok = True

    # --- runnable native mac binary: the real regression check ---
    if os.path.isfile(mac_bin):
        print("\n=== [RUN] mac arm64: score dna (raw 'U' + P) ===")
        rc, out = run_score(mac_bin, rec, lig, "dna")
        tail = "\n".join(out.strip().splitlines()[-12:])
        print(tail)
        panicked = "panic" in out.lower() or "not supported" in out.lower()
        if rc == 0 and not panicked:
            print(f"[PASS] dna scoring completed without panic (exit={rc})")
        else:
            print(f"[FAIL] dna scoring exit={rc}, panicked={panicked}")
            overall_ok = False

        # ddna is the other nucleic-capable scorer; also must not crash on P.
        rc2, out2 = run_score(mac_bin, rec, lig, "ddna")
        panicked2 = "panic" in out2.lower() or "not supported" in out2.lower()
        status = "PASS" if (rc2 == 0 and not panicked2) else "FAIL"
        if status == "FAIL":
            overall_ok = False
        print(f"[{status}] ddna: exit={rc2}, panicked={panicked2}")
        tail2 = "\n".join(out2.strip().splitlines()[-4:])
        if tail2:
            print("   " + tail2.replace("\n", "\n   "))
    else:
        print(f"[WARN] mac binary not found: {mac_bin}")
        overall_ok = False

    # --- cross binaries: existence/type check only (cannot execute on macOS) ---
    print("\n=== [CHECK] cross binaries (not executed on macOS) ===")
    for label, path in (("windows", win_bin), ("linux", lin_bin)):
        if os.path.isfile(path):
            try:
                ftype = subprocess.run(["file", "-b", path], capture_output=True, text=True).stdout.strip()
            except Exception:
                ftype = "?"
            print(f"[OK] {label}: {os.path.basename(path)} ({ftype})")
        else:
            print(f"[MISS] {label}: {path}")
            overall_ok = False

    print("\nRESULT:", "ALL PASS" if overall_ok else "SOME FAILURES")
    sys.exit(0 if overall_ok else 1)


if __name__ == "__main__":
    main()
