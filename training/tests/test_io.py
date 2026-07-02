from geneid_train.core.fasta import (
    read_fasta,
    read_tbl,
    write_fasta,
    write_tbl,
)
from geneid_train.core.gff import GffRecord, read_gff, write_gff


def test_fasta_roundtrip_and_header_token(tmp_path):
    src = tmp_path / "s.fa"
    src.write_text(">chr1 description here\nACGTACGT\nAAAA\n>chr2\nTTTT\n")
    recs = read_fasta(src)
    assert recs == {"chr1": "ACGTACGTAAAA", "chr2": "TTTT"}
    out = tmp_path / "out.fa"
    write_fasta(recs, out, width=4)
    assert read_fasta(out) == recs


def test_tbl_roundtrip(tmp_path):
    recs = {"a": "ACGT", "b": "TTTTGGGG"}
    p = tmp_path / "s.tbl"
    write_tbl(recs, p)
    assert read_tbl(p) == recs


def test_gff_roundtrip_and_gene_id(tmp_path):
    p = tmp_path / "s.gff"
    p.write_text(
        "chr1\tsrc\tCDS\t11\t16\t.\t+\t0\tgene1\n"
        "chr1\tsrc\tCDS\t27\t35\t.\t+\t0\tgene1\n"
    )
    recs = read_gff(p)
    assert len(recs) == 2
    assert recs[0].gene_id == "gene1"
    assert recs[0].start == 11 and recs[0].end == 16
    out = tmp_path / "out.gff"
    write_gff(recs, out)
    assert out.read_text() == p.read_text()


def test_gff_gene_id_is_bare_group_token():
    # geneid training GFF2 uses the bare gene id (optionally quoted) as the group.
    assert GffRecord("c", "s", "CDS", 1, 9, ".", "+", "0", "abc").gene_id == "abc"
    assert GffRecord("c", "s", "CDS", 1, 9, ".", "+", "0", '"abc"').gene_id == "abc"
