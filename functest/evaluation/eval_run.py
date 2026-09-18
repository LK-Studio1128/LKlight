#!/usr/bin/env python3
"""Pool the decoys of a multi-swarm LKlight run and report CAPRI/DockQ curves.

Reads every `rank_s*.list` in a run directory (as written by successive
`lklight rank <swarm> <N>` calls), pools the decoys into one global ranking by
the reported scoring value, evaluates each decoy with `dockq_eval.py`, and
writes

    eval_metrics.csv   one row per decoy: score, Fnat, L-RMSD, iRMSD, DockQ, CAPRI
    eval_summary.json  pooled top-N curve, oracle curve, sampling ceiling

The "scoring-ranked" curve is what a real run reports; the "oracle" curve ranks
the same decoys by true L-RMSD and therefore isolates *sampling* from *scoring*
— if the oracle curve is also flat, the sampling budget never produced a
near-native pose.

Usage
    eval_run.py <run_dir> <receptor.pdb> <ligand.pdb> [--rec-chain A ...]
"""
import argparse
import csv
import glob
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import dockq_eval as de  # noqa: E402


def pooled_rows(run_dir):
    rows = []
    for f in sorted(glob.glob(os.path.join(run_dir, "rank_s*.list"))):
        for i, ln in enumerate(open(f, errors="ignore")):
            if i == 0 or not ln.strip():
                continue
            p = ln.split()
            try:
                rows.append((float(p[4]), p[5]))
            except (IndexError, ValueError):
                continue
    rows.sort(key=lambda r: r[0])
    # `lklight rank` appends cumulatively across swarms, so the rank_s*.list
    # files overlap heavily; keep one row per decoy (identical score) so that
    # the pooled set is the full unique ensemble, not a duplicated subset.
    seen, dedup = set(), []
    for score, path in rows:
        if path in seen:
            continue
        seen.add(path)
        dedup.append((score, path))
    return dedup


def curve(seq, ns):
    out = []
    for n in ns:
        sub = seq[:n]
        if not sub:
            continue
        best = max(sub, key=lambda r: r["dockq"])
        out.append(dict(n=n,
                        hits=sum(1 for r in sub if r["capri"] != "Incorrect"),
                        best_dockq=best["dockq"], best_fnat=best["fnat"],
                        best_lrmsd=best["lrmsd"], best_irmsd=best["irmsd"],
                        best_capri=best["capri"]))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run_dir")
    ap.add_argument("receptor")
    ap.add_argument("ligand")
    ap.add_argument("--rec-chain", dest="rec_chains", action="append", default=None)
    ap.add_argument("--max", type=int, default=2000)
    args = ap.parse_args()
    rec_chains = set(args.rec_chains or ["A"])

    rd = os.path.abspath(args.run_dir)
    ref = de.Reference(args.receptor, args.ligand, rec_chains)
    rows = pooled_rows(rd)[: args.max]

    metrics, cache = [], {}
    for score, path in rows:
        if path not in cache:
            try:
                cache[path] = ref.evaluate(os.path.join(rd, path))
            except FileNotFoundError:
                cache[path] = None
        if cache[path]:
            m = dict(score=score, path=path, **cache[path])
            metrics.append(m)

    with open(os.path.join(rd, "eval_metrics.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["score", "decoy", "fnat", "lrmsd", "irmsd", "dockq", "capri"])
        for m in metrics:
            w.writerow([f"{m['score']:.6f}", m["path"], f"{m['fnat']:.4f}",
                        f"{m['lrmsd']:.3f}", f"{m['irmsd']:.3f}",
                        f"{m['dockq']:.4f}", m["capri"]])

    ns = [1, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000]
    scored = curve(metrics, ns)
    oracle = curve(sorted(metrics, key=lambda r: r["lrmsd"]), ns)
    ceiling = min(metrics, key=lambda r: r["lrmsd"])
    summary = dict(
        run_dir=os.path.basename(rd),
        n_decoys=len(metrics),
        native_pairs=len(ref.pairs),
        interface_residues=dict(receptor=len(ref.iface_r), ligand=len(ref.iface_l)),
        acceptable_in_pool=sum(1 for m in metrics if m["capri"] != "Incorrect"),
        decoys_with_contacts=sum(1 for m in metrics if m["fnat"] > 0),
        sampling_ceiling=dict(lrmsd=ceiling["lrmsd"], fnat=ceiling["fnat"],
                              irmsd=ceiling["irmsd"], dockq=ceiling["dockq"]),
        scoring_ranked=scored, oracle=oracle,
    )
    with open(os.path.join(rd, "eval_summary.json"), "w") as fh:
        json.dump(summary, fh, indent=1)

    print(f"# {summary['run_dir']}: {len(metrics)} decoys, "
          f"{summary['acceptable_in_pool']} acceptable-or-better, "
          f"sampling ceiling L-RMSD {ceiling['lrmsd']:.2f} A "
          f"(Fnat {ceiling['fnat']:.3f}, DockQ {ceiling['dockq']:.3f})")
    print(f"{'top-N':>6} {'hits':>5} {'bestDockQ':>9} {'bestFnat':>9} "
          f"{'bestL-RMSD':>11} {'bestCAPRI':>10}   | {'hits':>5} {'bestL-RMSD':>11}   (oracle)")
    for s, o in zip(scored, oracle):
        print(f"{s['n']:>6} {s['hits']:>5} {s['best_dockq']:>9.3f} "
              f"{s['best_fnat']:>9.3f} {s['best_lrmsd']:>11.2f} {s['best_capri']:>10}   | "
              f"{o['hits']:>5} {o['best_lrmsd']:>11.2f}")


if __name__ == "__main__":
    main()
