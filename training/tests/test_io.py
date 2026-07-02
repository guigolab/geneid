from geneid_train.core.fasta import (
    read_fasta,
    read_tbl,
    write_fasta,
    write_tbl,
)
from geneid_train.core.gff import (
    GffRecord,
    parse_attributes,
    read_gff3,
    write_gff3,
)


def test_fasta_roundtrip_and_header_token(tmp_path):
    src = tmp_path / "s.fa"
    src.write_text(">chr1 description here\nACGTACGT\nAAAA\n>chr2\nTTTT\n")
    recs = read_fasta(src)
    assert recs == {"chr1": "ACGTACGTAAAA", "chr2": "TTTT"}
    out = tmp_path / "out.fa"
    write_fasta(recs, out, width=4)
    assert read_fasta(out) == recs


def test_fasta_reads_gzip_transparently(tmp_path):
    import gzip

    p = tmp_path / "s.fa.gz"
    with gzip.open(p, "wt") as fh:
        fh.write(">chr1\nACGTACGT\n>chr2\nTTTT\n")
    assert read_fasta(p) == {"chr1": "ACGTACGT", "chr2": "TTTT"}


def test_tbl_roundtrip(tmp_path):
    recs = {"a": "ACGT", "b": "TTTTGGGG"}
    p = tmp_path / "s.tbl"
    write_tbl(recs, p)
    assert read_tbl(p) == recs


def test_gff3_read_and_attributes(tmp_path):
    p = tmp_path / "s.gff3"
    p.write_text(
        "##gff-version 3\n"
        "chr1\tsrc\tCDS\t11\t16\t.\t+\t0\tID=t1.cds1;Parent=t1\n"
        "chr1\tsrc\tCDS\t27\t35\t.\t+\t0\tID=t1.cds2;Parent=t1\n"
    )
    recs = read_gff3(p)
    assert len(recs) == 2
    assert recs[0].id == "t1.cds1"
    assert recs[0].parents == ["t1"]
    assert recs[0].start == 11 and recs[0].end == 16


def test_gff3_roundtrip(tmp_path):
    recs = [
        GffRecord("chr1", "src", "CDS", 11, 16, ".", "+", "0", {"ID": "t1.cds1", "Parent": "t1"}),
    ]
    out = tmp_path / "out.gff3"
    write_gff3(recs, out)
    reread = read_gff3(out)
    assert reread[0].to_line() == recs[0].to_line()


def test_parse_attributes_multi_parent():
    attrs = parse_attributes("ID=c1;Parent=t1,t2")
    assert attrs == {"ID": "c1", "Parent": "t1,t2"}
    rec = GffRecord("c", "s", "CDS", 1, 9, ".", "+", "0", attrs)
    assert rec.parents == ["t1", "t2"]
