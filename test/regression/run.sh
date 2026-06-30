#!/usr/bin/env bash
#
# geneid regression suite.
#
# Each case runs geneid with a fixed param/sequence/flags and compares the
# GFF output (minus the volatile "# date" header) against a committed golden.
# This is the guardrail for the Phase 2 cleanup: any change that alters the
# output of a covered code path shows up here as a FAIL.
#
#   ./run.sh            build (MEM=medium) and verify every case
#   ./run.sh --no-build verify using the existing bin/geneid
#   ./run.sh --bless    (re)generate the goldens from current output
#   ./run.sh <name>...  run only the named case(s)
#
# Add a case by appending one line to CASES below:  name|param|flags|sequence
# (extra per-case args, e.g. -R evidence or -j/-k bounds, go in the flags field).
#
set -uo pipefail
cd "$(dirname "$0")/../.." || exit 2          # repo root
ROOT=$(pwd)
GOLDEN="$ROOT/test/regression/golden"
BIN="$ROOT/bin/geneid"

# name | param | flags | sequence
# Cases mirror the documented workflows in the Current Protocols geneid
# chapter (Alioto et al.). The two MORC3/chr21 cases focus on a single locus
# via -j/-k so they stay fast (~1s) while keeping native chr21 coordinates.
CASES=(
  # snake: large realistic single-isochore genome, full output paths
  "snake|param/Hemorrhois_hippocrepis.geneid.optimized.param|-3UDTA|samples/rHemHip.H1.SUPER_1.1Mb.fasta"
  # human: 3-isochore selection + protein/cDNA/tDNA on the documented example
  "human|param/human3iso.param|-3UDTA|samples/example1.fa"
  # morc_u12: U12 intron prediction on the MORC3 locus (2 real U12 introns), introns printed
  "morc_u12|param/human3isoU12.param|-3UnDTA -j 36315000 -k 36380000|samples/chr21.fa"
  # rnaseq: RNA-seq evidence -- intron junctions (-R), expression coverage (-S), UTRs (-u)
  "rnaseq|param/human.rnaseq.param|-3U -u -R samples/ENCFF001.1.MORC.introns.gff -S samples/ENCFF001.1.MORC.stranded.expression.shuffled.gff -j 36315000 -k 36380000|samples/chr21.fa"
)

# chr21 cases need the unzipped fasta; derive it from the tracked .gz on demand
# (chr21.fa itself is too big to track -- see gunzip step in the chapter).
ensure_chr21() {
  [ -f samples/chr21.fa ] && return 0
  echo "# unzipping samples/chr21.fa from GRCh38.chr21.fa.gz ..."
  gunzip -c samples/GRCh38.chr21.fa.gz > samples/chr21.fa
}

norm() { grep -v '^# date' "$1"; }            # strip only the volatile date line

BUILD=1 BLESS=0; ONLY=()
for a in "$@"; do case "$a" in
  --no-build) BUILD=0 ;;
  --bless)    BLESS=1 ;;
  -*)         echo "unknown option: $a" >&2; exit 2 ;;
  *)          ONLY+=("$a") ;;
esac; done

if [ "$BUILD" = 1 ]; then
  echo "# building (MEM=medium)..."
  make MEM=medium >/dev/null 2>&1 || { echo "BUILD FAILED"; exit 2; }
fi
[ -x "$BIN" ] || { echo "no binary at $BIN (drop --no-build?)"; exit 2; }

pass=0 fail=0 rc=0
for spec in "${CASES[@]}"; do
  IFS='|' read -r name param flags seq <<<"$spec"
  if [ "${#ONLY[@]}" -gt 0 ] && [[ ! " ${ONLY[*]} " == *" $name "* ]]; then continue; fi
  [ "$seq" = "samples/chr21.fa" ] && ensure_chr21
  printf '%-12s ' "$name"
  act=$(mktemp); g="$GOLDEN/$name.norm"
  # shellcheck disable=SC2086  # flags is intentionally word-split
  if ! "$BIN" -P "$param" $flags "$seq" >"$act" 2>/dev/null; then
    echo "RUN-ERROR (geneid exited nonzero)"; fail=$((fail+1)); rc=1; rm -f "$act"; continue
  fi
  if [ "$BLESS" = 1 ]; then
    norm "$act" >"$g"; echo "blessed ($(wc -l <"$g" | tr -d ' ') lines)"
  elif [ ! -f "$g" ]; then
    echo "NO GOLDEN (run --bless first)"; fail=$((fail+1)); rc=1
  elif diff -q <(norm "$act") "$g" >/dev/null; then
    echo "PASS ($(norm "$act" | wc -l | tr -d ' ') lines)"; pass=$((pass+1))
  else
    echo "FAIL"; diff <(norm "$act") "$g" | head -20; fail=$((fail+1)); rc=1
  fi
  rm -f "$act"
done

echo "# ${pass} passed, ${fail} failed"
exit $rc
