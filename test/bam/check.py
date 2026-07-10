#!/usr/bin/env python3
"""Validate the C bamcov reader (bamdump) against oracle.json.
Usage: check.py <bamdump-binary> <fixture-dir>   (dir has synth.bam, oracle.json)

Per query the oracle stores the ground-truth per-base coverage (uncovered = 0) as
non-zero runs. Expand that and bamdump's reported runs to per-base arrays over the
window and compare -- an exact per-base match regardless of how runs are split."""
import json
import subprocess
import sys

bamdump, d = sys.argv[1], sys.argv[2]
oracle = json.load(open(f"{d}/oracle.json"))
fails = total = 0


def expand(runs, s, e):
    a = [0] * (e - s)
    for rs, re, v in runs:
        for p in range(max(rs, s), min(re, e)):
            a[p - s] = v
    return a


for q in oracle["queries"]:
    total += 1
    s, e = q["start"], q["end"]
    want = expand([(r[0], r[1], r[2]) for r in q["runs"]], s, e)
    out = subprocess.run([bamdump, f"{d}/synth.bam", q["chrom"], str(s), str(e)],
                         capture_output=True, text=True)
    got_runs = [(int(p[0]), int(p[1]), int(float(p[2])))
                for p in (ln.split("\t") for ln in out.stdout.splitlines())]
    got = expand(got_runs, s, e)
    if want != got:
        fails += 1
        diffs = [(s + i, want[i], got[i]) for i in range(len(want)) if want[i] != got[i]]
        print(f"FAIL {q['chrom']}:{s}-{e}  ({len(diffs)} base(s) differ)")
        print(f"  first diffs (pos, want, got): {diffs[:6]}")
        if out.stderr.strip():
            print(f"  stderr: {out.stderr.strip()}")

print(f"{total - fails}/{total} queries match" + (" -- ALL PASS" if not fails else f" -- {fails} FAILED"))
sys.exit(1 if fails else 0)
