import pytest

from geneid_train.core.param import Param
from geneid_train.param.gene_model import extract_intron_range, utr_gene_model_lines


def test_utr_lines_use_canonical_names_and_zero_intergenic():
    lines = utr_gene_model_lines("40:10000")
    text = "\n".join(lines)
    # canonical exon-type names (geneid.h sUTR*), not the legacy typo
    assert "UTR_5prime_Internal_Half" in text
    assert "UTR_5Internal_Half" not in text
    assert "40:10000" in text  # the intron range is carried through verbatim
    # the four intergenic rules must use a 0 minimum (UTRs may abut/overlap)
    inter = [ln for ln in lines
             if ln.startswith(("aataaa+:Terminal+", "Promoter-:First-"))]
    assert len(inter) == 4
    assert all(ln.split()[-1] == "0:Infinity" for ln in inter)


def test_utr_lines_preserve_sentinel_for_assembly():
    # when assembling from the template the range is still the sentinel
    assert any("@INTRON_RANGE@" in ln for ln in utr_gene_model_lines("@INTRON_RANGE@"))


def test_extract_intron_range():
    p = Param.from_text(
        "General_Gene_Model\n"
        "First+:Internal+  Internal+:Terminal+  40:23712.339\n"
        "First+  Intron+  1:1\n"
    )
    assert extract_intron_range(p) == "40:23712.339"


def test_extract_intron_range_raises_on_nonstandard():
    p = Param.from_text("General_Gene_Model\nWeird+  Rule+  1:1\n")
    with pytest.raises(ValueError):
        extract_intron_range(p)
