#!/usr/bin/env bash
#
# geneid regression suite.
#
# Each case runs geneid with a fixed param/sequence/flags and compares the
# GFF output (minus the volatile "# date" header) against a committed golden.
# This is the guardrail for the Phase 2 cleanup: any change that alters the
# output of a covered code path shows up here as a FAIL.
#
#   ./run.sh            build and verify every case
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
  # longprot: 153-exon / 23313-aa gene (the longest protein on SUPER_1) -- exercises
  # the growable whole-protein + cDNA/tDNA buffers; window cut from SUPER_1 ~199.14-199.30 Mb
  "longprot|param/Hemorrhois_hippocrepis.geneid.optimized.param|-3UDTA|samples/rHemHip.H1.SUPER_1.longprot.fasta"
  # human: 3-isochore selection + protein/cDNA/tDNA on the documented example
  "human|param/human3iso.param|-3UDTA|samples/example1.fa"
  # morc_u12: U12 intron prediction on the MORC3 locus (2 real U12 introns), introns printed
  "morc_u12|param/human3isoU12.param|-3UnDTA -j 36315000 -k 36380000|samples/chr21.fa"
  # rnaseq: RNA-seq evidence -- intron junctions (-R), expression coverage (-S), UTRs (-u)
  "rnaseq|param/human.rnaseq.param|-3U -u -R samples/ENCFF001.1.MORC.introns.gff -S samples/ENCFF001.1.MORC.stranded.expression.shuffled.gff -j 36315000 -k 36380000|samples/chr21.fa"
  # rnaseq_bw: same RNA-seq case but the -S expression coverage is delivered as two
  # stranded bigWigs (plus.bw,minus.bw) via per-split range queries instead of the
  # text GFF. Byte-identical to the rnaseq golden -> guards bwQuery -> FillCoverage ->
  # HSPScan2 and the 0-based bigWig -> 1-based sr[] coordinate mapping on both strands.
  "rnaseq_bw|param/human.rnaseq.param|-3U -u -R samples/ENCFF001.1.MORC.introns.gff -S samples/ENCFF001.1.MORC.stranded.plus.bw,samples/ENCFF001.1.MORC.stranded.minus.bw -j 36315000 -k 36380000|samples/chr21.fa"
  # rnaseq_bw_nou: bigWig -S coverage WITHOUT -u -- guards that coverage scores
  # exons (via sr[]) with no UTR prediction and no readcount[] allocation (the
  # -u-decoupling): must run (no crash) and predict CDS but emit no UTR lines.
  "rnaseq_bw_nou|param/human.rnaseq.param|-3 -R samples/ENCFF001.1.MORC.introns.gff -S samples/ENCFF001.1.MORC.stranded.plus.bw,samples/ENCFF001.1.MORC.stranded.minus.bw -j 36315000 -k 36380000|samples/chr21.fa"
  # human_intron: soft intron-length penalty on a 250kb Red-masked chr14 slice
  # (GRCh38 chr14:67,900,000-68,150,000, RAD51B locus). The param is chr12-trained
  # and carries Intron_length_model 7.1788 1.51411 with weight 0.2 (feature ON) and
  # a generous 500kb gene-model cap. With the penalty on, two far-band introns
  # (>27kb: 48.6kb + 29.0kb) are penalized away vs weight 0, so this case exercises
  # and pins the convex-hinge penalty + far-band fast-DP path (introns printed via -n).
  "human_intron|param/human.chr12.intron_length.param|-3Un|samples/human.chr14.longintron.fasta"
  # human_intron_multifrag: same soft intron-length penalty on a 550kb slice (chr14
  # 67,700,000-68,250,000) -- >500kb (LENGTHSi) so it is processed in TWO fragments,
  # exercising the near-band deque's cross-fragment maintenance in BackupArrayD (the
  # deque-index rebase/expiry) that the single-fragment human_intron case does not.
  "human_intron_multifrag|param/human.chr12.intron_length.param|-3Un|samples/human.chr14.longintron.multifrag.fasta"
  # human_intron_bb: bigBed -R evidence (per-split range query). The same two Intron
  # junctions as a GFF -R would supply, delivered via bigBed; the routing waives the
  # length penalty so both far-band introns (28.9kb + 48.6kb) are spanned. Exercises
  # bbOpen->bbQuery->AddEvidenceExon end-to-end; byte-identical to the GFF -R result.
  "human_intron_bb|param/human.chr12.intron_length.param|-3Un -R samples/human.chr14.longintron.introns.bb|samples/human.chr14.longintron.fasta"
  # --- assemble-only (-O) cases: feed pre-typed exons, skip ab initio prediction ---
  # human_o: single-locus -O assembly of example1's own 8 exons (fast single-split).
  # Guards the -O path's independent nExons accounting (was an uninitialized read).
  "human_o|param/human3iso.param|-3n -O samples/example1.geneid.gff|samples/example1.fa"
  # human_o_multilocus: -O over a 2-record FASTA with the same 8 exons under each
  # locus -- guards multi-locus assembly (two genes, one per sequence).
  "human_o_multilocus|param/human3iso.param|-3n -O samples/example1.2locus.geneid.gff|samples/example1.2locus.fa"
  # morc_jo_u12: -J -O annotation-scoring + U12 typing over the whole chr21 (multi-
  # split), forcing MORC3's 17 CDS exons. Guards (a) -U allowed under -O, (b) the
  # multi-split -J classify walking the printed GOptim chain (2 U12 introns typed).
  "morc_jo_u12|param/human3isoU12.param|-J -3nU -O samples/MORC.CDS.geneid.gff|samples/chr21.fa"
  # morc_o_utr: -O -u UTR assembly over whole chr21, forcing MORC3's CDS + both UTR
  # halves. Guards -u allowed under -O and the UTR-exon frame/remainder fix (the 3'
  # UTR_Terminal_Half, whose length is not a multiple of 3, must not be dropped).
  "morc_o_utr|param/human.rnaseq.param|-3nUu -O samples/MORC.UTR.geneid.gff|samples/chr21.fa"
  # morc_rc_rnaseq: the MORC3 RNA-seq case reverse-complemented -- the gene is on
  # the '-' strand of a standalone morc_rc sequence, with reversed intron (-R) and
  # stranded coverage (-S) evidence. Exercises the reverse-strand coordinate paths
  # (sr[] genomic<->RSequence flip, '-' evidence) end to end; predicts the mirror
  # 17-CDS/2-UTR gene on '-'. Fixtures regenerated by make_morc_rc.py.
  "morc_rc_rnaseq|param/human.rnaseq.param|-3U -u -R samples/morc_rc.introns.gff -S samples/morc_rc.plus.bw,samples/morc_rc.minus.bw|samples/morc_rc.fa"
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
  echo "# building..."
  # bin/geneid is a tracked binary, so a git checkout can leave it with a
  # newer mtime than the objects and make would skip the relink -- running a
  # stale binary against the goldens. Remove it first to force a fresh link.
  rm -f "$BIN"
  make >/dev/null 2>&1 || { echo "BUILD FAILED"; exit 2; }
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
