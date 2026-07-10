#!/usr/bin/env python3
"""Validate the C bamcov junction reader (juncdump) against oracle_junc.json.
Usage: check_junc.py <juncdump-binary> <fixture-dir>   (dir has synth_junc.bam, oracle_junc.json)

juncdump prints one line per junction occurrence (start, end, strand); we tally
identical junctions whose start falls in the query window and compare to the
oracle's aggregated [start, end, strand, count] list."""
import json
import subprocess
import sys

juncdump, d = sys.argv[1], sys.argv[2]
oracle = json.load(open(f"{d}/oracle_junc.json"))
fails = total = 0

for q in oracle["queries"]:
    total += 1
    s, e = q["start"], q["end"]
    out = subprocess.run([juncdump, f"{d}/synth_junc.bam", q["chrom"], str(s), str(e)],
                         capture_output=True, text=True)
    agg = {}
    for ln in out.stdout.splitlines():
        p = ln.split("\t")
        js, je, st = int(p[0]), int(p[1]), p[2]
        if s <= js < e:                       # oracle counts by junction start in window
            agg[(js, je, st)] = agg.get((js, je, st), 0) + 1
    got = sorted([list(k) + [v] for k, v in agg.items()])
    want = sorted([list(x) for x in q["junctions"]])
    if got != want:
        fails += 1
        print(f"FAIL {q['chrom']}:{s}-{e}")
        print(f"  want {want}")
        print(f"  got  {got}")
        if out.stderr.strip():
            print(f"  stderr: {out.stderr.strip()}")

print(f"{total - fails}/{total} junction queries match" + (" -- ALL PASS" if not fails else f" -- {fails} FAILED"))
sys.exit(1 if fails else 0)
