from geneid_train.prepare.base import Exon, GeneModel
from geneid_train.prepare.classify import classify_report

# chr1: exon[1-6] GT-AG intron[7-16] exon[17-22] GC-AG intron[23-32] exon[33-38]
CHR1 = "AAAAAA" + "GTAAAAAAAG" + "CCCCCC" + "GCAAAAAAAG" + "GGGGGG"
# chr2: exon[1-6] AT-AC intron[7-16] exon[17-22]
CHR2 = "AAAAAA" + "ATAAAAAAAC" + "GGGGGG"
GENOME = {"chr1": CHR1, "chr2": CHR2}

M1 = GeneModel("g1", "chr1", "+", [Exon(1, 6), Exon(17, 22), Exon(33, 38)])
M2 = GeneModel("g2", "chr2", "+", [Exon(1, 6), Exon(17, 22)])


def test_tally_counts_each_class():
    rep = classify_report([M1, M2], GENOME, min_sites=1)
    assert rep.n_introns == 3
    assert rep.pair_counts[("GT", "AG")] == 1
    assert rep.pair_counts[("GC", "AG")] == 1
    assert rep.pair_counts[("AT", "AC")] == 1
    assert rep.donor_counts == {"GT": 1, "GC": 1, "AT": 1}


def test_recommendation_train_vs_transplant():
    by_name = {c.name: c for c in classify_report([M1, M2], GENOME, min_sites=1).classes}
    assert "train de novo" in by_name["U2 GC-AG"].recommendation
    # bulk class is never gated
    assert "bulk" in by_name["U2 GT-AG"].recommendation
    # U12 GT-AG can't be seen by dinucleotide alone
    assert by_name["U12 GT-AG"].count == -1
    assert "bootstrap" in by_name["U12 GT-AG"].recommendation


def test_recommendation_transplant_when_below_threshold():
    by_name = {c.name: c for c in classify_report([M1, M2], GENOME, min_sites=50).classes}
    assert "transplant" in by_name["U2 GC-AG"].recommendation
    assert "transplant" in by_name["U12 AT-AC"].recommendation


def test_n_containing_boundaries_skipped():
    genome = {"c": "AAAAAA" + "GTNAAAANAG" + "GGGGGG"}
    m = GeneModel("g", "c", "+", [Exon(1, 6), Exon(17, 22)])
    rep = classify_report([m], genome, min_sites=1)
    # donor "GT" ok, acceptor "AG" ok here (Ns are internal), so it counts
    assert rep.n_introns == 1
    # now make the acceptor contain N
    genome2 = {"c": "AAAAAA" + "GTAAAAAANN" + "GGGGGG"}
    rep2 = classify_report([m], genome2, min_sites=1)
    assert rep2.n_introns == 0
