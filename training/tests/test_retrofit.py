import pytest

from geneid_train.core.param import Param
from geneid_train.param.gene_model import HUMAN_INTRON_LENGTH_MODEL
from geneid_train.param.retrofit import retrofit_param

# minimal single-isochore CDS-only param: a gene model with the standard
# intragenic rule, and an Exon_weights anchor for the intron-length insertion.
MINI = (
    "number_of_isochores\n1\n"
    "General_Gene_Model\n"
    "First+:Internal+  Internal+:Terminal+  40:9000\n"
    "First+  Intron+  1:1\n"
    "Exon_weights\n1 1 1 1\n"
)


def test_retrofit_utr_adds_rules_and_preserves_range():
    out = retrofit_param(MINI, utr=True)
    p = Param.from_text(out)
    gm = p.data_lines("General_Gene_Model")
    assert any("UTR_First_Half" in ln for ln in gm)
    assert any(ln.split()[-1] == "0:Infinity" for ln in gm)  # UTR intergenic min 0
    intragenic = next(ln for ln in gm if ln.startswith("First+:Internal+"))
    assert intragenic.split()[-1] == "40:9000"  # the param's own intron range kept


def test_retrofit_intron_length_human_inserted_before_exon_weights():
    out = retrofit_param(MINI, intron_length=HUMAN_INTRON_LENGTH_MODEL, intron_length_weight=0.0)
    p = Param.from_text(out)
    assert p.scalar("Intron_length_model") == "7.1788 1.51411"
    assert p.scalar("Intron_length_score_weight") == "0"
    kws = p.keywords()
    assert kws.index("Intron_length_model") < kws.index("Exon_weights")


def test_retrofit_is_idempotent():
    once = retrofit_param(MINI, utr=True, intron_length=HUMAN_INTRON_LENGTH_MODEL)
    twice = retrofit_param(once, utr=True, intron_length=HUMAN_INTRON_LENGTH_MODEL)
    assert once == twice  # second pass detects both are present and injects nothing


def test_retrofit_rejects_multi_isochore():
    multi = MINI.replace("number_of_isochores\n1\n", "number_of_isochores\n3\n")
    with pytest.raises(NotImplementedError):
        retrofit_param(multi, utr=True)


def test_retrofit_roundtrips_other_sections_unchanged():
    # a section we don't touch must be byte-identical afterwards
    src = MINI + "Some_Other_Section\n0.5\n"
    out = retrofit_param(src, utr=True)
    assert "Some_Other_Section\n0.5\n" in out
