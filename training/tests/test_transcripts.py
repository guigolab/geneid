from geneid_train.core.gff import GffRecord
from geneid_train.prepare.base import build_models, collapse_isoforms
from geneid_train.prepare.transcripts import load


def _rec(type_, start, end, attrs, strand="+"):
    return GffRecord("chr1", "td", type_, start, end, ".", strand, "0" if type_ == "CDS" else ".",
                     attrs)


def _two_isoform_gff():
    # gene g1 with two transcripts: t1 (short, 1 CDS) and t2 (long, 2 CDS)
    return [
        _rec("gene", 100, 500, {"ID": "g1"}),
        _rec("mRNA", 100, 200, {"ID": "t1", "Parent": "g1"}),
        _rec("CDS", 100, 160, {"ID": "t1.c1", "Parent": "t1"}),
        _rec("mRNA", 100, 500, {"ID": "t2", "Parent": "g1"}),
        _rec("CDS", 100, 200, {"ID": "t2.c1", "Parent": "t2"}),
        _rec("CDS", 400, 500, {"ID": "t2.c2", "Parent": "t2"}),
    ]


def test_collapse_picks_longest_cds_per_gene():
    records = _two_isoform_gff()
    models = build_models(records)
    assert len(models) == 2  # two transcripts
    collapsed = collapse_isoforms(models, records)
    assert len(collapsed) == 1  # one gene
    assert collapsed[0].gene_id == "t2"  # longest CDS (61+101 vs 61)
    assert collapsed[0].coding_length == (200 - 100 + 1) + (500 - 400 + 1)


def test_load_transcripts_end_to_end(tmp_path):
    p = tmp_path / "td.gff3"
    from geneid_train.core.gff import write_gff3

    write_gff3(_two_isoform_gff(), p)
    models = load(p)
    assert len(models) == 1
    assert models[0].gene_id == "t2"


def test_model_without_gene_record_is_its_own_gene():
    # bare CDS with Parent but no mRNA/gene record
    recs = [
        GffRecord("chr1", "td", "CDS", 1, 30, ".", "+", "0", {"ID": "c", "Parent": "orphan"}),
    ]
    models = build_models(recs)
    collapsed = collapse_isoforms(models, recs)
    assert len(collapsed) == 1
    assert collapsed[0].gene_id == "orphan"
