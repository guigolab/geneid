from geneid_train.core.gff import GffRecord
from geneid_train.prepare.base import (
    Exon,
    GeneModel,
    build_models,
    extract_locus,
    filter_complete,
    filter_min_protein,
    filter_non_overlapping,
)

# Synthetic chr1 (1-based): flank[1-10], exon1[11-16]=ATGAAA,
# intron[17-26]=GTAAAAAAAG, exon2[27-35]=AAATTTTAA, flank[36-45].
CHROM = "CCCCCCCCCC" + "ATGAAA" + "GTAAAAAAAG" + "AAATTTTAA" + "GGGGGGGGGG"
GENOME = {"chr1": CHROM}


def _gene1_records():
    return [
        GffRecord("chr1", "s", "CDS", 11, 16, ".", "+", "0", {"ID": "c1", "Parent": "gene1"}),
        GffRecord("chr1", "s", "CDS", 27, 35, ".", "+", "0", {"ID": "c2", "Parent": "gene1"}),
    ]


def test_build_and_cds():
    (m,) = build_models(_gene1_records())
    assert m.gene_id == "gene1"
    assert m.is_multiexonic
    assert m.cds(GENOME) == "ATGAAAAAATTTTAA"
    assert m.protein(GENOME) == "MKKF*"
    assert m.introns() == [(17, 26)]
    assert m.intron_seqs(GENOME) == ["GTAAAAAAAG"]  # donor GT / acceptor AG
    assert m.is_complete(GENOME)


def test_build_rejects_mixed_strand():
    recs = _gene1_records()
    recs[1].strand = "-"
    try:
        build_models(recs)
    except ValueError as e:
        assert "mixed strands" in str(e)
    else:
        raise AssertionError("expected ValueError on mixed strands")


def test_filters_complete_and_min_protein():
    models = build_models(_gene1_records())
    assert len(filter_complete(models, GENOME)) == 1
    assert len(filter_min_protein(models, GENOME, min_aa=4)) == 1  # MKKF
    assert len(filter_min_protein(models, GENOME, min_aa=5)) == 0


def _model(gid, start, end, seqid="chr1", strand="+"):
    return GeneModel(gid, seqid, strand, [Exon(start, end, "0")])


def test_non_overlapping_drops_both_of_a_pair():
    a = _model("a", 10, 100)
    b = _model("b", 90, 200, strand="-")  # overlaps a on opposite strand -> both drop
    c = _model("c", 300, 400)  # clear
    kept = filter_non_overlapping([a, b, c])
    assert [m.gene_id for m in kept] == ["c"]


def test_non_overlapping_keeps_different_seqids():
    a = _model("a", 10, 100, seqid="chr1")
    b = _model("b", 50, 150, seqid="chr2")
    kept = filter_non_overlapping([a, b])
    assert {m.gene_id for m in kept} == {"a", "b"}


def test_extract_locus_shifts_coordinates():
    (m,) = build_models(_gene1_records())
    locus = extract_locus(m, GENOME, flank=5)
    assert locus.seq == CHROM[5:40]  # lo=5 (0-based), hi=40
    assert [(r.start, r.end) for r in locus.records] == [(6, 11), (22, 30)]
    assert all(r.seqid == "gene1" for r in locus.records)
