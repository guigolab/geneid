
from geneid_train.evaluate import (
    Exon,
    _exact_matches,
    _gene_matches,
    _no_overlap_count,
    _overlap_nt,
    evaluate_files,
    read_annotation_gff,
)

# The synthetic pred/real pair below was cross-checked against Enrique Blanco's
# compiled `evaluation` tool: its #Total line (what the geneid optimiser parses)
# reads SN=0.92 SP=1.00 CC=0.95 SNe=0.80 SPe=0.80 SNSP=0.80 SNg=SPg=SNSPg=0.50.

REAL = (
    "LocusA\tanno\tSequence\t1\t5000\t.\t+\t.\tLocusA\n"
    "LocusA\tanno\tFirst\t401\t500\t.\t+\t.\tgeneA\n"
    "LocusA\tanno\tInternal\t1000\t1100\t.\t+\t.\tgeneA\n"
    "LocusA\tanno\tTerminal\t2000\t2100\t.\t+\t.\tgeneA\n"
    "#$\n"
    "LocusB\tanno\tSequence\t1\t3000\t.\t+\t.\tLocusB\n"
    "LocusB\tanno\tFirst\t401\t600\t.\t+\t.\tgeneB\n"
    "LocusB\tanno\tTerminal\t1500\t1600\t.\t+\t.\tgeneB\n"
)
PRED = (
    "LocusA\tgeneid\tFirst\t401\t500\t5\t+\t0\tgeneA\n"
    "LocusA\tgeneid\tInternal\t1000\t1100\t3\t+\t1\tgeneA\n"
    "LocusA\tgeneid\tTerminal\t2000\t2050\t2\t+\t2\tgeneA\n"
    "#$\n"
    "LocusB\tgeneid\tFirst\t401\t600\t4\t+\t0\tgeneB\n"
    "LocusB\tgeneid\tTerminal\t1500\t1600\t3\t+\t1\tgeneB\n"
)


def _write_pair(tmp_path):
    real = tmp_path / "real.gff"
    pred = tmp_path / "pred.gff"
    real.write_text(REAL)
    pred.write_text(PRED)
    return pred, real


def test_evaluate_files_matches_reference_evaluation_tool(tmp_path):
    pred, real = _write_pair(tmp_path)
    a = evaluate_files(pred, real)
    # accumulated counts (each verified against the C tool's per-locus output)
    t = a.totals
    assert (t.tp, t.cds_real, t.cds_pred) == (553, 603, 553)
    assert (t.tpe, t.exr, t.exp) == (4, 5, 5)
    assert (t.tpg, t.ger, t.gep) == (1, 2, 2)
    # total-level metrics == the tool's #Total line
    assert round(a.sn, 2) == 0.92
    assert round(a.sp, 2) == 1.00
    assert round(a.cc, 2) == 0.95
    assert a.sne == 0.8 and a.spe == 0.8 and a.snsp == 0.8
    assert a.sng == 0.5 and a.spg == 0.5 and a.snspg == 0.5
    assert a.ra_me == 0.0 and a.ra_we == 0.0


def test_annotation_first_line_is_info_not_exon(tmp_path):
    real = tmp_path / "real.gff"
    real.write_text(REAL)
    loci = read_annotation_gff(real)
    assert [locus.name for locus in loci] == ["LocusA", "LocusB"]
    assert loci[0].length == 5000  # field 5 of the consumed info line
    # the Sequence info line is not counted as an exon (3 exons, not 4)
    assert len(loci[0].exons) == 3
    assert loci[1].length == 3000 and len(loci[1].exons) == 2


def test_overlap_nt_partial_and_full():
    a = [Exon(100, 200, "+", "g")]
    b = [Exon(150, 250, "+", "g")]
    assert _overlap_nt(a, b) == 51  # 150..200 inclusive
    assert _overlap_nt(a, a) == 101  # full self-overlap


def test_exact_matches_requires_both_boundaries():
    pred = [Exon(10, 20, "+", "g"), Exon(30, 41, "+", "g")]
    real = [Exon(10, 20, "+", "g"), Exon(30, 40, "+", "g")]
    assert _exact_matches(pred, real) == 1  # only the first is an exact match


def test_no_overlap_count_is_missing_or_wrong():
    real = [Exon(10, 20, "+", "g"), Exon(100, 110, "+", "g")]
    pred = [Exon(12, 18, "+", "g")]  # overlaps first, misses second
    assert _no_overlap_count(real, pred) == 1  # ME: second real exon unmatched
    assert _no_overlap_count(pred, real) == 0  # WE: the prediction overlaps


def test_gene_matches_needs_identical_structure():
    real = [Exon(10, 20, "+", "gA"), Exon(50, 60, "+", "gA")]
    exact = [Exon(10, 20, "+", "x"), Exon(50, 60, "+", "x")]
    off = [Exon(10, 20, "+", "y"), Exon(50, 61, "+", "y")]
    assert _gene_matches(exact, real) == 1
    assert _gene_matches(off, real) == 0


def test_exact_match_is_strand_sensitive():
    # identical coordinates on opposite strands are not the same exon
    pred = [Exon(10, 20, "+", "g")]
    real = [Exon(10, 20, "-", "g")]
    assert _exact_matches(pred, real) == 0
