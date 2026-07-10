#!/usr/bin/env python3
"""Validate the C bigwig reader (bwdump) against oracle.json (from pybigtools).
Usage: check.py <bwdump-binary> <fixture-dir>   (dir has small.bw, many.bw, oracle.json)

For each query the oracle stores the ground-truth per-base signal (uncovered = 0)
as non-zero runs. We expand that to a per-base array over the window and compare
it to the same expansion of bwdump's reported intervals (clipped to the window).
An exact per-base match is format-independent (any section split decodes the
same), so bedGraph/varStep/fixedStep all validate through one path."""
import json
import subprocess
import sys

TOL = 1e-4
bwdump, d = sys.argv[1], sys.argv[2]
oracle = json.load(open(f"{d}/oracle.json"))
fails = total = 0


def expand(runs, s, e):
    a = [0.0] * (e - s)
    for rs, re, v in runs:
        for p in range(max(rs, s), min(re, e)):
            a[p - s] = v
    return a


for fn, spec in oracle["files"].items():
    for q in spec["queries"]:
        total += 1
        s, e = q["start"], q["end"]
        want = expand([(r[0], r[1], r[2]) for r in q["runs"]], s, e)
        out = subprocess.run([bwdump, f"{d}/{fn}", q["chrom"], str(s), str(e)],
                             capture_output=True, text=True)
        got_runs = []
        for ln in out.stdout.splitlines():
            p = ln.split("\t")
            got_runs.append((int(p[0]), int(p[1]), float(p[2])))
        got = expand(got_runs, s, e)
        if any(abs(a - b) > TOL for a, b in zip(want, got)):
            fails += 1
            diffs = [(s + i, want[i], got[i]) for i in range(len(want))
                     if abs(want[i] - got[i]) > TOL]
            print(f"FAIL {fn} {q['chrom']}:{s}-{e}  ({len(diffs)} base(s) differ)")
            print(f"  first diffs (pos, want, got): {diffs[:6]}")
            if out.stderr.strip():
                print(f"  stderr: {out.stderr.strip()}")

print(f"{total - fails}/{total} queries match" + (" -- ALL PASS" if not fails else f" -- {fails} FAILED"))
sys.exit(1 if fails else 0)
