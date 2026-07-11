#!/usr/bin/env python3
"""Validate strand-aware coverage (readTxnStrand + filtering) against
oracle_strand.json. Usage: check_strand.py <bamdump-binary> <fixture-dir>

Each query fixes a wanted transcription strand + library mode; bamdump is run
with those, and its coverage runs are compared to the reads the oracle expects
to survive the strand filter."""
import json
import subprocess
import sys

bamdump, d = sys.argv[1], sys.argv[2]
oracle = json.load(open(f"{d}/oracle_strand.json"))
fails = total = 0

for q in oracle["queries"]:
    total += 1
    args = [bamdump, f"{d}/synth_str.bam", q["chrom"], str(q["start"]), str(q["end"]),
            q["want"], q["libmode"]]
    out = subprocess.run(args, capture_output=True, text=True)
    got = sorted([[int(p[0]), int(p[1]), int(float(p[2]))]
                  for p in (ln.split("\t") for ln in out.stdout.splitlines())])
    want = sorted([list(r) for r in q["runs"]])
    if got != want:
        fails += 1
        print(f"FAIL want={q['want']} lib={q['libmode']}")
        print(f"  want {want}")
        print(f"  got  {got}")
        if out.stderr.strip():
            print(f"  stderr: {out.stderr.strip()}")

print(f"{total - fails}/{total} strand queries match" + (" -- ALL PASS" if not fails else f" -- {fails} FAILED"))
sys.exit(1 if fails else 0)
