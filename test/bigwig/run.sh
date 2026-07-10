#!/usr/bin/env bash
#
# Standalone test for the bigWig reader (src/bigwig.c). Builds the bwdump driver
# and checks its range-query output against oracle.json. The .bw fixtures and the
# oracle are committed, so this needs only gcc + zlib + python3 (NOT pybigtools).
# Regenerate the fixtures/oracle with make_fixtures.py in an env with pybigtools.
set -uo pipefail
cd "$(dirname "$0")" || exit 2
ROOT=../..
D=$(mktemp -d)
gcc -I"$ROOT/include" -Wall -O2 "$ROOT/src/bigwig.c" bwdump.c -o "$D/bwdump" -lz || { echo "BUILD FAILED"; exit 2; }
# unit test of emitSection (all three wiggle encodings, incl var/fixedStep)
gcc -I"$ROOT/include" -Wall -O2 section_unit.c -o "$D/section_unit" -lz || { echo "UNIT BUILD FAILED"; exit 2; }
"$D/section_unit" || exit 1
# oracle test of the full reader against pybigtools ground truth
python3 check.py "$D/bwdump" .
