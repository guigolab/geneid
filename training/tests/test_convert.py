from geneid_train.core.convert import gff2_to_gff3, gtf_to_gff3
from geneid_train.prepare.base import build_models

GENOME = {"chr1": "CCCCCCCCCC" + "ATGAAA" + "GTAAAAAAAG" + "AAATTTTAA" + "GGGGGGGGGG"}


def test_gff2_to_gff3_builds_hierarchy_and_models():
    lines = [
        "chr1\ts\tCDS\t11\t16\t.\t+\t0\tgene1\n",
        "chr1\ts\tCDS\t27\t35\t.\t+\t0\tgene1\n",
    ]
    recs = gff2_to_gff3(lines)
    types = [r.type for r in recs]
    assert types == ["gene", "mRNA", "CDS", "CDS"]
    # gene id is suffixed to stay unique vs the same-named transcript
    assert recs[0].id == "gene1.gene"  # gene
    assert recs[1].id == "gene1" and recs[1].parents == ["gene1.gene"]  # mRNA -> gene
    assert all(r.parents == ["gene1"] for r in recs if r.type == "CDS")
    # feeds the IR cleanly
    (m,) = build_models(recs)
    assert m.cds(GENOME) == "ATGAAAAAATTTTAA"


def test_gtf_to_gff3_uses_transcript_and_gene_ids():
    lines = [
        'chr1\ts\tCDS\t11\t16\t.\t+\t0\ttranscript_id "t1"; gene_id "g1";\n',
        'chr1\ts\tCDS\t27\t35\t.\t+\t0\ttranscript_id "t1"; gene_id "g1";\n',
    ]
    recs = gtf_to_gff3(lines)
    gene = next(r for r in recs if r.type == "gene")
    mrna = next(r for r in recs if r.type == "mRNA")
    assert gene.id == "g1"
    assert mrna.id == "t1" and mrna.parents == ["g1"]
    assert [r.type for r in recs if r.type == "CDS"] == ["CDS", "CDS"]


def test_gff2_quoted_group_token():
    lines = ['chr1\ts\tCDS\t11\t16\t.\t+\t0\tgene_id "abc"\n']
    recs = gff2_to_gff3(lines)
    assert next(r for r in recs if r.type == "CDS").parents == ["abc"]
