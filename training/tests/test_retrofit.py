import pytest

from geneid_train.core.param import Param
from geneid_train.param.gene_model import HUMAN_INTRON_LENGTH_MODEL
from geneid_train.param.retrofit import retrofit_param

# minimal single-isochore CDS-only param: Donor/Acceptor_profile anchors (for the
# U12 insertion), the standard gene-model intragenic rule, and an Exon_weights
# anchor (for the intron-length + U12-threshold insertions).
MINI = (
    "number_of_isochores\n1\n"
    "Donor_profile\n8 1 -7 1 0 1\n1 AA 0\n"
    "Acceptor_profile\n30 28 -7 1 0 1\n1 AA 0\n"
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


def test_retrofit_u12_inserts_trio_and_gates_in_the_right_slots():
    out = retrofit_param(MINI, u12=True)
    p = Param.from_text(out)
    for name in ("U12_Branch_point_profile", "U12gtag_Donor_profile", "U12atac_Donor_profile",
                 "U12gtag_Acceptor_profile", "U12atac_Acceptor_profile"):
        assert p.has(name)
    assert p.scalar("U12_Splice_Score_Threshold") == "9"
    kws = p.keywords()
    # optional profiles precede the required profile they extend
    assert kws.index("U12gtag_Donor_profile") < kws.index("Donor_profile")
    assert kws.index("U12_Branch_point_profile") < kws.index("Acceptor_profile")
    assert kws.index("U12_Splice_Score_Threshold") < kws.index("Exon_weights")


def test_retrofit_is_idempotent():
    kw = dict(utr=True, intron_length=HUMAN_INTRON_LENGTH_MODEL, u12=True)
    once = retrofit_param(MINI, **kw)
    twice = retrofit_param(once, **kw)
    assert once == twice  # second pass detects all present and injects nothing


def test_retrofit_rejects_multi_isochore():
    multi = MINI.replace("number_of_isochores\n1\n", "number_of_isochores\n3\n")
    with pytest.raises(NotImplementedError):
        retrofit_param(multi, utr=True)


def test_retrofit_roundtrips_other_sections_unchanged():
    # a section we don't touch must be byte-identical afterwards
    src = MINI + "Some_Other_Section\n0.5\n"
    out = retrofit_param(src, utr=True)
    assert "Some_Other_Section\n0.5\n" in out
