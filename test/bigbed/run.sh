#!/usr/bin/env bash
#
# Standalone test for the bigBed reader (src/bigbed.c). Builds the bbdump driver
# and checks its range-query output against oracle.json. The .bb fixtures and the
# oracle are committed, so this needs only gcc + zlib + python3 (NOT pybigtools).
# Regenerate the fixtures/oracle with make_fixtures.py in an env with pybigtools.
set -uo pipefail
cd "$(dirname "$0")" || exit 2
ROOT=../..
BIN=$(mktemp -d)/bbdump
gcc -I"$ROOT/include" -Wall -O2 "$ROOT/src/bigbed.c" bbdump.c -o "$BIN" -lz || { echo "BUILD FAILED"; exit 2; }
python3 check.py "$BIN" .
