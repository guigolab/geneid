from geneid_train.prepare.base import Exon, GeneModel
from geneid_train.prepare.classify import classify_report, detect_u12_gtag

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


def test_u12_gtag_bootstrap_disabled_falls_back_to_punt():
    by_name = {
        c.name: c
        for c in classify_report([M1, M2], GENOME, min_sites=1, bootstrap_u12=False).classes
    }
    assert by_name["U12 GT-AG"].count == -1
    assert "bootstrap" in by_name["U12 GT-AG"].recommendation


def test_u12_gtag_screen_reports_estimate():
    # the lone A-rich GT-AG intron scores far more U2 than U12 -> not a candidate
    rep = classify_report([M1, M2], GENOME, min_sites=1)
    by_name = {c.name: c for c in rep.classes}
    assert by_name["U12 GT-AG"].count == 0
    assert rep.u12_gtag is not None and rep.u12_gtag.n_candidates == 0


def test_recommendation_transplant_when_below_threshold():
    by_name = {c.name: c for c in classify_report([M1, M2], GENOME, min_sites=50).classes}
    assert "transplant" in by_name["U2 GC-AG"].recommendation
    assert "transplant" in by_name["U12 AT-AC"].recommendation


def test_u12_gtag_screen_discriminates_donor_motif():
    # among GT-AG introns, one with the U12 5'SS (GTATCCTT) must score more U12 than
    # U2, while the A-rich ones define the (U2) donor model and score negative
    exon = "AAAAAA"
    u12 = "GTATCCTTAC" + "A" * 30 + "AG"
    plain = ["GTAAGTAAAA" + "A" * 30 + "AG"] * 6  # bulk U2-like donors
    seqs = [u12, *plain]
    # lay the introns out on one contig as consecutive genes
    models, chrom = [], ""
    for i, iseq in enumerate(seqs):
        s = len(chrom)
        chrom += exon + iseq + exon
        e1 = s + 1
        a2 = s + 6 + len(iseq) + 1
        models.append(GeneModel(f"g{i}", "c", "+", [Exon(e1, e1 + 5), Exon(a2, a2 + 5)]))
    est = detect_u12_gtag(models, {"c": chrom}, floor=-1e9)
    assert est.n_scored == len(seqs)
    # the single U12-5'SS intron is the top U12-minus-U2 margin and is positive
    assert est.top_margins[0] > 0
    # with the default floor the A-rich U2 introns are NOT called
    est2 = detect_u12_gtag(models, {"c": chrom})  # bundled floor
    assert est2.n_candidates == 1


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
