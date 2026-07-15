"""Truth-prep for the benchmark harness: GENCODE GFF3 -> geneid gp format."""

from geneid_train.benchmark import (
    evaluate_asymmetric,
    fasta_seqlen,
    gencode_cds_gp,
    gencode_mane_gp,
)

# gene A: '+' strand, MANE, protein-coding, 3 CDS -> First/Internal/Terminal.
# gene B: '-' strand, MANE, protein-coding, 2 CDS -> transcription order is the
#         high-coord exon First, low-coord exon Terminal.
# gene C: MANE but NOT protein-coding -> excluded.
# gene D: protein-coding but NOT MANE (an alternative isoform) -> excluded.
GFF3 = """\
##gff-version 3
chr9\tX\tgene\t100\t900\t.\t+\t.\tID=gA
chr9\tX\ttranscript\t100\t900\t.\t+\t.\tID=txA;Parent=gA;gene_type=protein_coding;tag=basic,MANE_Select
chr9\tX\tCDS\t100\t150\t.\t+\t0\tID=c1;Parent=txA
chr9\tX\tCDS\t300\t360\t.\t+\t0\tID=c2;Parent=txA
chr9\tX\tCDS\t800\t900\t.\t+\t0\tID=c3;Parent=txA
chr9\tX\ttranscript\t100\t900\t.\t+\t.\tID=txD;Parent=gA;gene_type=protein_coding;tag=basic
chr9\tX\tCDS\t100\t200\t.\t+\t0\tID=d1;Parent=txD
chr9\tX\ttranscript\t2000\t2600\t.\t-\t.\tID=txB;Parent=gB;gene_type=protein_coding;tag=MANE_Select
chr9\tX\tCDS\t2000\t2100\t.\t-\t0\tID=b1;Parent=txB
chr9\tX\tCDS\t2500\t2600\t.\t-\t0\tID=b2;Parent=txB
chr9\tX\ttranscript\t3000\t3300\t.\t+\t.\tID=txC;Parent=gC;gene_type=lncRNA;tag=MANE_Select
chr9\tX\tCDS\t3000\t3300\t.\t+\t0\tID=cc1;Parent=txC
chr8\tX\ttranscript\t10\t99\t.\t+\t.\tID=txE;Parent=gE;gene_type=protein_coding;tag=MANE_Select
chr8\tX\tCDS\t10\t99\t.\t+\t0\tID=e1;Parent=txE
"""


def _parse(gp: str):
    lines = [ln.split("\t") for ln in gp.strip().splitlines()]
    return lines[0], lines[1:]


def test_gencode_mane_gp_filters_and_types(tmp_path):
    f = tmp_path / "g.gff3"
    f.write_text(GFF3)
    info, exons = _parse(gencode_mane_gp(f, "chr9", 5000))

    # info line: field 5 (index 4) is the sequence length
    assert info[0] == "chr9" and info[4] == "5000"

    # only the two MANE protein-coding chr9 transcripts survive (not C/D/chr8)
    groups = {e[8] for e in exons}
    assert groups == {"txA", "txB"}

    by_tid = {tid: [e for e in exons if e[8] == tid] for tid in groups}
    # '+' gene: types follow genomic order
    assert [e[2] for e in by_tid["txA"]] == ["First", "Internal", "Terminal"]
    # '-' gene: First is the high-coord CDS (transcription order)
    b = by_tid["txB"]
    first = next(e for e in b if e[2] == "First")
    term = next(e for e in b if e[2] == "Terminal")
    assert (first[3], first[4]) == ("2500", "2600")
    assert (term[3], term[4]) == ("2000", "2100")
    assert all(e[6] == "-" for e in b)


def test_single_cds_is_typed_single(tmp_path):
    gff = (
        "chr1\tX\ttranscript\t1\t9\t.\t+\t.\tID=t;Parent=g;"
        "gene_type=protein_coding;tag=MANE_Select\n"
        "chr1\tX\tCDS\t1\t9\t.\t+\t0\tID=c;Parent=t\n"
    )
    f = tmp_path / "s.gff3"
    f.write_text(gff)
    _, exons = _parse(gencode_mane_gp(f, "chr1", 100))
    assert len(exons) == 1 and exons[0][2] == "Single"


def test_all_isoform_truth_includes_non_mane(tmp_path):
    f = tmp_path / "g.gff3"
    f.write_text(GFF3)
    # canonical: only MANE txA/txB; all: also the non-MANE isoform txD (and txC,
    # which carries a CDS) -- every CDS-bearing transcript
    _, mane = _parse(gencode_mane_gp(f, "chr9", 5000))
    _, allx = _parse(gencode_cds_gp(f, "chr9", 5000, canonical_only=False))
    assert {e[8] for e in mane} == {"txA", "txB"}
    assert {"txA", "txB", "txD"} <= {e[8] for e in allx}


def _gp(tmp_path, name, rows):
    p = tmp_path / name
    lines = ["chr1\tx\tinfo\t1\t100000\t.\t.\t.\ti"]
    lines += [f"chr1\tx\t{t}\t{s}\t{e}\t.\t+\t.\t{g}" for t, s, e, g in rows]
    p.write_text("\n".join(lines) + "\n")
    return str(p)


def test_asymmetric_credits_noncanonical_isoform_to_precision(tmp_path):
    # canonical (recall): a 3-exon transcript; all (precision): also an isoform
    # whose single exon is (100,200). A prediction of exactly (100,200) matches
    # NO canonical exon (recall exon SN stays 0) but DOES match the isoform, so
    # precision credits it (exon SP = 1).
    recall = _gp(tmp_path, "r.gp", [
        ("First", 100, 150, "txA"), ("Internal", 300, 360, "txA"),
        ("Terminal", 800, 900, "txA"),
    ])
    precision = _gp(tmp_path, "a.gp", [
        ("First", 100, 150, "txA"), ("Internal", 300, 360, "txA"),
        ("Terminal", 800, 900, "txA"), ("Single", 100, 200, "txD"),
    ])
    pred = tmp_path / "pred.gff"
    pred.write_text("chr1\tgeneid\tSingle\t100\t200\t.\t+\t.\tg1\n")

    a = evaluate_asymmetric(str(pred), recall, precision)
    assert a.sne == 0.0          # exact exon not in the canonical reference
    assert a.spe == 1.0          # but matches an annotated isoform -> not a false call
    assert a.spg == 1.0          # gene structure matches the isoform too
    assert a.sp == 1.0           # every predicted nt lies within an annotated CDS
    assert 0.0 < a.sn < 1.0      # partial nucleotide overlap with the canonical exon


def test_asymmetric_precision_nt_never_exceeds_one_with_overlapping_isoforms(tmp_path):
    # two isoforms whose CDS overlap; a prediction spanning both must not
    # double-count (SP <= 1) thanks to the merged reference.
    recall = _gp(tmp_path, "r.gp", [("Single", 100, 300, "txA")])
    precision = _gp(tmp_path, "a.gp", [
        ("Single", 100, 250, "txA"), ("Single", 200, 300, "txB"),
    ])
    pred = tmp_path / "pred.gff"
    pred.write_text("chr1\tgeneid\tSingle\t100\t300\t.\t+\t.\tg1\n")
    a = evaluate_asymmetric(str(pred), recall, precision)
    assert a.sp == 1.0


def test_fasta_seqlen_picks_named_record(tmp_path):
    fa = tmp_path / "g.fa"
    fa.write_text(">chr1 desc\nACGTAC\nGT\n>chr2\nAAAA\n")
    assert fasta_seqlen(fa, "chr1") == 8
    assert fasta_seqlen(fa, "chr2") == 4
