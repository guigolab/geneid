#!/usr/bin/env bash
#
# Standalone test for the BAM coverage reader (src/bamcov.c). Builds the bamdump
# driver against htslib and checks its per-region coverage against oracle.json.
# The synth.bam/.bai + oracle are committed, so this needs a WITH_HTSLIB toolchain
# (htslib headers+lib) + python3 -- NOT samtools/pysam. Regenerate the fixture
# with make_fixture.py (needs samtools).
#
# htslib is opt-in; if its headers aren't found this test SKIPS (exit 0) rather
# than failing, matching the WITH_HTSLIB build being optional.
set -uo pipefail
cd "$(dirname "$0")" || exit 2
ROOT=../..

# Locate an htslib prefix (env override, else common install locations).
PREFIX="${HTSLIB_PREFIX:-}"
if [ -z "$PREFIX" ]; then
  for p in /usr/local /opt/homebrew /usr; do
    [ -f "$p/include/htslib/sam.h" ] && PREFIX="$p" && break
  done
fi
if [ -z "$PREFIX" ] || [ ! -f "$PREFIX/include/htslib/sam.h" ]; then
  echo "SKIP: htslib not found (set HTSLIB_PREFIX); bamcov test needs a WITH_HTSLIB toolchain"
  exit 0
fi

D=$(mktemp -d)
gcc -I"$ROOT/include" -I"$PREFIX/include" -Wall -O2 "$ROOT/src/bamcov.c" bamdump.c \
    -o "$D/bamdump" -L"$PREFIX/lib" -lhts || { echo "BUILD FAILED (bamdump)"; exit 2; }
gcc -I"$ROOT/include" -I"$PREFIX/include" -Wall -O2 "$ROOT/src/bamcov.c" juncdump.c \
    -o "$D/juncdump" -L"$PREFIX/lib" -lhts || { echo "BUILD FAILED (juncdump)"; exit 2; }
python3 check.py "$D/bamdump" . || exit 1
python3 check_junc.py "$D/juncdump" .
