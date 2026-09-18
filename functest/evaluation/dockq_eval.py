#!/usr/bin/env python3
"""Standard interface-quality indicators for LKlight decoys.

For every decoy, relative to a bound (native) reference complex:

  Fnat     fraction of native residue-residue contacts recovered
           (heavy-atom pairs within 5.0 A, one entry per residue pair)
  L-RMSD   ligand RMSD after least-squares superposition on the receptor
           backbone (ligand measured over all heavy atoms, which keeps the
           definition meaningful for nucleic-acid ligands)
  iRMSD    interface RMSD: interface residues are those with any heavy-atom
           pair within 10.0 A across the partner molecules in the native;
           those residues (both sides) are superposed as rigid bodies on
           their all-heavy-atom coordinates and the same atoms are used for
           the RMSD
  DockQ    continuous score of Basu & Wallner (2016),
           DockQ = (Fnat + 1/(1+(iRMSD/1.5)^2) + 1/(1+(L-RMSD/8.5)^2)) / 3
  CAPRI    discrete class from (Fnat, L-RMSD) after Mendez et al. (2003):
           High       Fnat >= 0.5 and L-RMSD <= 1.0 A
           Medium     Fnat >= 0.3 and L-RMSD <= 2.0 A
           Acceptable Fnat >= 0.1 and L-RMSD <= 4.0 A
           Incorrect  otherwise

Residue identity is keyed on (chain, resseq, icode), so multi-chain ligands
(for example a two-strand DNA duplex) are handled correctly.  The receptor is
one or more explicit chains; the ligand is every other chain in the file.
Water and hetero records are ignored.  No external docking library is used, so
the numbers are reproducible from the archived PDB files alone.

Usage
    dockq_eval.py --curve <receptor.pdb> <ligand.pdb> <decoy_prefix> <Nmax>
                  [--rec-chain A]
    dockq_eval.py <receptor.pdb> <ligand.pdb> <decoy.pdb> [...] [--rec-chain A]
"""
import argparse
import numpy as np

CONTACT_CUTOFF = 5.0      # A, heavy-atom, for Fnat
INTERFACE_CUTOFF = 10.0   # A, heavy-atom, for the interface definition
BACKBONE = {"N", "CA", "C", "O"}
ALT_OK = (" ", "A")


def parse_pdb(path):
    """Return a list of (chain, resseq, icode, atomname, xyz) heavy atoms."""
    atoms = []
    for ln in open(path, errors="ignore"):
        if not ln.startswith("ATOM  "):
            continue
        if ln[16] not in ALT_OK:
            continue
        elem = (ln[76:78].strip() or ln[12:16].strip()[0]).upper()
        if elem == "H":
            continue
        atoms.append((ln[21], int(ln[22:26]), ln[26], ln[12:16].strip(),
                      (float(ln[30:38]), float(ln[38:46]), float(ln[46:54]))))
    return atoms


def split_by_chains(atoms, rec_chains):
    rec = [a for a in atoms if a[0] in rec_chains]
    lig = [a for a in atoms if a[0] not in rec_chains]
    return rec, lig


def group_residues(atoms):
    res = {}
    for ch, rs, ic, an, xyz in atoms:
        res.setdefault((ch, rs, ic), []).append(xyz)
    return res


def residue_contacts(rec, lig, cutoff=CONTACT_CUTOFF):
    """Set of ((rec chain,resseq,icode), (lig chain,resseq,icode)) in contact."""
    cut2 = cutoff * cutoff
    rres, lres = group_residues(rec), group_residues(lig)
    larr = [(k, np.asarray(v)) for k, v in lres.items()]
    pairs = set()
    for rk, rxyz in rres.items():
        ra = np.asarray(rxyz)
        for lk, la in larr:
            d2 = ((ra[:, None, :] - la[None, :, :]) ** 2).sum(-1)
            if d2.min() <= cut2:
                pairs.add((rk, lk))
    return pairs


def interface_residues(rec, lig, cutoff=INTERFACE_CUTOFF):
    cut2 = cutoff * cutoff
    rres, lres = group_residues(rec), group_residues(lig)
    larr = [(k, np.asarray(v)) for k, v in lres.items()]
    ir, il = set(), set()
    for rk, rxyz in rres.items():
        ra = np.asarray(rxyz)
        for lk, la in larr:
            if ((ra[:, None, :] - la[None, :, :]) ** 2).sum(-1).min() <= cut2:
                ir.add(rk)
                il.add(lk)
    return ir, il


def atoms_of(atoms, wanted, backbone_only=False):
    out = []
    for ch, rs, ic, an, xyz in atoms:
        if (ch, rs, ic) not in wanted:
            continue
        if backbone_only and an not in BACKBONE:
            continue
        out.append(((ch, rs, ic, an), xyz))
    return out


def match(atoms_a, atoms_b, backbone_only=False, subset=None):
    idx = {}
    for ch, rs, ic, an, xyz in atoms_b:
        if backbone_only and an not in BACKBONE:
            continue
        if subset is not None and (ch, rs, ic) not in subset:
            continue
        idx[(ch, rs, ic, an)] = xyz
    P, Q = [], []
    for ch, rs, ic, an, xyz in atoms_a:
        if backbone_only and an not in BACKBONE:
            continue
        if subset is not None and (ch, rs, ic) not in subset:
            continue
        k = (ch, rs, ic, an)
        if k in idx:
            P.append(xyz)
            Q.append(idx[k])
    return np.asarray(P), np.asarray(Q)


def kabsch(P, Q):
    Pc, Qc = P.mean(0), Q.mean(0)
    V, _, Wt = np.linalg.svd((P - Pc).T @ (Q - Qc))
    d = np.sign(np.linalg.det(V @ Wt))
    R = V @ np.diag([1.0, 1.0, d]) @ Wt
    return lambda x: (R @ (np.asarray(x) - Pc).T).T + Qc


def rmsd(A, B):
    return float(np.sqrt(((A - B) ** 2).sum(1).mean()))


def capri_class(fnat, lrmsd):
    if fnat >= 0.5 and lrmsd <= 1.0:
        return "High"
    if fnat >= 0.3 and lrmsd <= 2.0:
        return "Medium"
    if fnat >= 0.1 and lrmsd <= 4.0:
        return "Acceptable"
    return "Incorrect"


class Reference:
    def __init__(self, rec_path, lig_path, rec_chains):
        self.rec_chains = rec_chains
        atoms = parse_pdb(rec_path) + parse_pdb(lig_path)
        self.rec, self.lig = split_by_chains(atoms, rec_chains)
        self.pairs = residue_contacts(self.rec, self.lig)
        self.iface_r, self.iface_l = interface_residues(self.rec, self.lig)
        if not self.pairs:
            raise SystemExit(
                "reference complex has no interface contacts — check that the "
                "receptor/ligand files are the bound complex, not a separated pair"
            )

    def evaluate(self, decoy_path):
        atoms = parse_pdb(decoy_path)
        dec_rec, dec_lig = split_by_chains(atoms, self.rec_chains)
        dec_pairs = residue_contacts(dec_rec, dec_lig)
        fnat = len(dec_pairs & self.pairs) / len(self.pairs)

        P, Q = match(dec_rec, self.rec, backbone_only=True)
        if len(P) < 3:
            P, Q = match(dec_rec, self.rec)
        f = kabsch(P, Q)
        PL, QL = match(dec_lig, self.lig)
        lrmsd = rmsd(f(PL), QL) if len(PL) else float("nan")

        Pi_r, Qi_r = match(dec_rec, self.rec, subset=self.iface_r)
        Pi_l, Qi_l = match(dec_lig, self.lig, subset=self.iface_l)
        Pi, Qi = np.vstack([Pi_r, Pi_l]), np.vstack([Qi_r, Qi_l])
        irmsd = rmsd(kabsch(Pi, Qi)(Pi), Qi) if len(Pi) >= 3 else float("nan")

        dockq = (fnat
                 + 1.0 / (1.0 + (irmsd / 1.5) ** 2)
                 + 1.0 / (1.0 + (lrmsd / 8.5) ** 2)) / 3.0
        return {"fnat": fnat, "lrmsd": lrmsd, "irmsd": irmsd, "dockq": dockq,
                "capri": capri_class(fnat, lrmsd),
                "n_native_pairs": len(self.pairs), "n_decoy_pairs": len(dec_pairs)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("receptor")
    ap.add_argument("ligand")
    ap.add_argument("decoys", nargs="+")
    ap.add_argument("--rec-chain", dest="rec_chains", action="append", default=None,
                    help="receptor chain id (repeatable); default A")
    ap.add_argument("--curve", action="store_true",
                    help="treat the last DECOYS value as Nmax and expand "
                         "<prefix>1.pdb .. <prefix>Nmax.pdb")
    args = ap.parse_args()
    rec_chains = set(args.rec_chains or ["A"])
    ref = Reference(args.receptor, args.ligand, rec_chains)
    print(f"# reference: {ref.n_native_pairs if hasattr(ref,'n_native_pairs') else len(ref.pairs)} "
          f"native inter-molecular residue contacts; "
          f"interface residues {len(ref.iface_r)}R + {len(ref.iface_l)}L; "
          f"receptor chains {sorted(rec_chains)}")

    if args.curve:
        prefix, nmax = args.decoys[0], int(args.decoys[1])
        print(f"{'top-N':>6} {'success':>8} {'bestDockQ':>10} {'bestFnat':>9} "
              f"{'bestL-RMSD':>11} {'bestiRMSD':>10} {'bestCAPRI':>10}")
        for n in (1, 5, 10, 20, 50, 100, 200, 500, 1000):
            if n > nmax:
                break
            rows = []
            for i in range(1, n + 1):
                try:
                    rows.append(ref.evaluate(f"{prefix}{i}.pdb"))
                except FileNotFoundError:
                    pass
            if not rows:
                continue
            best = max(rows, key=lambda r: r["dockq"])
            ok = sum(1 for r in rows if r["capri"] != "Incorrect")
            print(f"{n:>6} {('YES' if ok else 'no'):>8} {best['dockq']:>10.3f} "
                  f"{best['fnat']:>9.3f} {best['lrmsd']:>11.2f} "
                  f"{best['irmsd']:>10.2f} {best['capri']:>10}")
        return

    print(f"{'decoy':>28} {'Fnat':>7} {'L-RMSD':>8} {'iRMSD':>7} {'DockQ':>7}  CAPRI")
    for d in args.decoys:
        r = ref.evaluate(d)
        print(f"{d:>28} {r['fnat']:>7.3f} {r['lrmsd']:>8.2f} {r['irmsd']:>7.2f} "
              f"{r['dockq']:>7.3f}  {r['capri']}")


if __name__ == "__main__":
    main()
