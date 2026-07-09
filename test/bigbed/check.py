#!/usr/bin/env python3
"""Validate the C bigbed reader (bbdump) against oracle.json (from pybigtools).
Usage: check.py <bbdump-binary> <fixture-dir>   (dir has small.bb, many.bb, oracle.json)"""
import json
import subprocess
import sys

bbdump, d = sys.argv[1], sys.argv[2]
oracle = json.load(open(f"{d}/oracle.json"))
fails = total = 0
for fn, spec in oracle["files"].items():
    for q in spec["queries"]:
        total += 1
        out = subprocess.run([bbdump, f"{d}/{fn}", q["chrom"], str(q["start"]), str(q["end"])],
                             capture_output=True, text=True)
        got = []
        for ln in out.stdout.splitlines():
            p = ln.split("\t")
            # bbdump prints: start end rest... ; oracle record: [start, end, name, score, strand]
            got.append([int(p[0]), int(p[1])] + p[2:])
        want = [[r[0], r[1]] + [str(x) for x in r[2:]] for r in q["records"]]
        got.sort(); want.sort()
        if got != want:
            fails += 1
            print(f"FAIL {fn} {q['chrom']}:{q['start']}-{q['end']}")
            print(f"  want {want}")
            print(f"  got  {got}")
            if out.stderr.strip():
                print(f"  stderr: {out.stderr.strip()}")
print(f"{total - fails}/{total} queries match" + (" -- ALL PASS" if not fails else f" -- {fails} FAILED"))
sys.exit(1 if fails else 0)
