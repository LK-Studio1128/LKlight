#!/usr/bin/env python3
"""Summarize equivalence benchmark results into paper-ready markdown tables.

Usage: python3 summarize_results.py <results_dir> [--top 1 5 10 20 50 100]
Reads <results_dir>/summary.tsv and prints markdown tables comparing engines.
"""
import argparse
import csv
import os
import sys


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("results_dir")
    ap.add_argument("--top", nargs="*", type=int, default=[1, 5, 10, 20, 50, 100])
    args = ap.parse_args()

    tsv = os.path.join(args.results_dir, "summary.tsv")
    rows = {}
    with open(tsv) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            rows.setdefault(r["case"], {})[r["engine"]] = r

    cases = sorted(rows)
    print(f"# Accuracy-equivalence summary ({len(cases)} cases)\n")
    print("| case | engine | n_solutions | "
          + " | ".join(f"top{n}" for n in args.top) + " | status |")
    print("|---|---|---|" + "---|" * len(args.top) + "---|")
    for c in cases:
        for eng in ("lk", "py"):
            r = rows[c].get(eng)
            if not r:
                continue
            vals = " | ".join(f"{r.get(f'top{n}', '-'):>6}" for n in args.top)
            print(f"| {c} | {eng} | {r.get('n_solutions', '-')} | {vals} | {r.get('status', 'ok')} |")

    # per-top-N mean over cases per engine
    print("\n## Mean success rate (across cases)\n")
    print("| top-N | LKlight | Python |")
    print("|---|---|---|")
    for n in args.top:
        means = {}
        for eng in ("lk", "py"):
            vals = []
            for c in cases:
                r = rows[c].get(eng)
                if r and r.get(f"top{n}", "-") not in ("-", ""):
                    try:
                        vals.append(float(r[f"top{n}"]))
                    except ValueError:
                        pass
            means[eng] = sum(vals) / len(vals) if vals else float("nan")
        print(f"| {n} | {means['lk']:.3f} | {means['py']:.3f} |")


if __name__ == "__main__":
    main()
