import pytest

from geneid_train.core.param import Param
from geneid_train.param.assemble import assemble_param, load_template
from geneid_train.stats.coding import format_markov_matrix
from geneid_train.stats.genemodel import format_range, intron_range
from geneid_train.stats.sites import SiteWindow, _fmt, format_profile, profile_header

from .conftest import ref_dir


def test_profile_header_order0_and_order1():
    w0 = SiteWindow(start=25, end=32, offset=3, length=8, st=None, nd=None, rd=None)
    assert profile_header(w0, order=0) == ["8", "3", "-7", "0"]
    w1 = SiteWindow(start=29, end=36, offset=1, length=8, st=3, nd=4, rd=5)
    assert profile_header(w1, order=1) == ["8", "1", "-7", "1", "0", "1"]


def test_format_profile_structure():
    matrix = {(1, "A"): 0.5, (1, "C"): -1.25413, (2, "G"): -9999.0}
    lines = format_profile(matrix, ["8", "3", "-7", "0"])
    assert lines[0] == "8 3 -7 0"
    assert lines[1] == "# Transition probabilities at every position"
    assert lines[2] == "1 A 0.5"
    assert lines[3] == "1 C -1.25413"
    assert lines[4] == "2 G -9999"


def test_fmt_integers_bare():
    assert _fmt(-9999.0) == "-9999"
    assert _fmt(0.0) == "0"
    assert _fmt(0.729033) == "0.729033"


def test_replace_block_data_preserves_surroundings():
    text = "Foo\n1\n# a comment\nBar\nx y z\n"
    p = Param.from_text(text)
    p.replace_block_data("Foo", ["9", "9", "9"])
    out = p.to_text()
    # keyword kept, data replaced, the Bar block untouched
    assert out == "Foo\n9\n9\n9\n# a comment\nBar\nx y z\n"


def _tiny_template():
    return (
        "# geneid parameter file: @SPECIES@\n"
        "Start_profile\n@START_PROFILE@\n"
        "Acceptor_profile\n@ACCEPTOR_PROFILE@\n"
        "Donor_profile\n@DONOR_PROFILE@\n"
        "Markov_order\n@MARKOV_ORDER@\n"
        "Markov_Initial_probability_matrix\n@MARKOV_INITIAL@\n"
        "Markov_Transition_probability_matrix\n@MARKOV_TRANSITION@\n"
        "General_Gene_Model\n"
        "A B @INTRON_RANGE@\n"
        "C D @INTERGENIC_RANGE@\n"
    )


def test_assemble_param_injects_all_sections():
    out = assemble_param(
        species="Testus_specius",
        start_profile=["8 3 -7 0", "1 A 0.5"],
        acceptor_profile=["30 28 -7 1 0 1", "1 AA 0.1"],
        donor_profile=["8 1 -7 1 0 1", "1 AA 0.2"],
        markov_order=5,
        markov_initial=["AAAAA 0 0 -0.5"],
        markov_transition=["AAAAAA 0 0 -0.8"],
        intron_range="24.75:25394.023",
        intergenic_range="200:Infinity",
        template=_tiny_template(),
    )
    assert "@" not in out  # every sentinel filled
    assert "# geneid parameter file: Testus_specius" in out
    p = Param.from_text(out)
    assert p.scalar("Markov_order") == "5"
    assert p.profile("Donor_profile").header == [8, 1, -7, 1, 0, 1]
    assert "A B 24.75:25394.023" in out
    assert "C D 200:Infinity" in out


def test_assemble_param_injects_intron_length_model():
    out = assemble_param(
        species="Testus_specius",
        start_profile=["8 3 -7 0", "1 A 0.5"],
        acceptor_profile=["30 28 -7 1 0 1", "1 AA 0.1"],
        donor_profile=["8 1 -7 1 0 1", "1 AA 0.2"],
        markov_order=5,
        markov_initial=["AAAAA 0 0 -0.5"],
        markov_transition=["AAAAAA 0 0 -0.8"],
        intron_range="24.75:25394.023",
        intergenic_range="200:Infinity",
        template=_tiny_template() + "Exon_weights\n1 1 1 1\n",
        intron_length_model=(7.5, 1.25),
    )
    # model + its weight land as an optional block before Exon_weights; the weight
    # defaults to 0.5 (penalty ON)
    assert "Intron_length_model\n7.5 1.25\n" in out
    assert "Intron_length_score_weight\n0.5\n" in out
    assert out.index("Intron_length_model") < out.index("Exon_weights")


def test_assemble_param_intron_length_weight_override():
    kw = dict(
        species="Testus_specius",
        start_profile=["8 3 -7 0", "1 A 0.5"],
        acceptor_profile=["30 28 -7 1 0 1", "1 AA 0.1"],
        donor_profile=["8 1 -7 1 0 1", "1 AA 0.2"],
        markov_order=5,
        markov_initial=["AAAAA 0 0 -0.5"],
        markov_transition=["AAAAAA 0 0 -0.8"],
        intron_range="24.75:25394.023",
        intergenic_range="200:Infinity",
        template=_tiny_template() + "Exon_weights\n1 1 1 1\n",
        intron_length_model=(7.5, 1.25),
    )
    assert "Intron_length_score_weight\n0\n" in assemble_param(**kw, intron_length_weight=0)
    assert "Intron_length_score_weight\n1.5\n" in assemble_param(**kw, intron_length_weight=1.5)


def test_assemble_param_omits_intron_length_model_when_absent():
    out = assemble_param(
        species="Testus_specius",
        start_profile=["8 3 -7 0"],
        acceptor_profile=["30 28 -7 1 0 1"],
        donor_profile=["8 1 -7 1 0 1"],
        markov_order=5,
        markov_initial=["AAAAA 0 0 -0.5"],
        markov_transition=["AAAAAA 0 0 -0.8"],
        intron_range="24.75:25394.023",
        intergenic_range="200:Infinity",
        template=_tiny_template(),
    )
    assert "Intron_length_model" not in out


def test_bundled_template_has_all_sentinels():
    t = load_template()
    for s in (
        "@START_PROFILE@",
        "@ACCEPTOR_PROFILE@",
        "@DONOR_PROFILE@",
        "@MARKOV_ORDER@",
        "@MARKOV_INITIAL@",
        "@MARKOV_TRANSITION@",
        "@INTRON_RANGE@",
        "@INTERGENIC_RANGE@",
        "@SPECIES@",
    ):
        assert s in t, s
    assert t.count("@INTRON_RANGE@") == 2
    assert t.count("@INTERGENIC_RANGE@") == 4


# ---- full assembly against the reference param ------------------------------

REF = ref_dir()
SP = "Xerocrassa_montserratensis"


@pytest.mark.integration
@pytest.mark.skipif(REF is None, reason="set GENEID_TRAIN_REFDIR to the train_geneid dir")
def test_assembled_param_matches_reference():
    refp = Param.read(REF / f"{SP}.geneid.param")

    def prof_lines(name):
        pr = refp.profile(name)
        matrix = {(p, o): v for p, o, v in pr.rows}
        return format_profile(matrix, [_fmt(x) for x in pr.header])

    def markov_lines(name):
        blk = refp._find(name)
        model = {}
        for i in blk.data_line_indices():
            parts = blk.raw[i].split()
            model[(parts[0], int(parts[2]))] = float(parts[3])
        return format_markov_matrix(model)

    lens = [len(ln.split("\t")[1].strip()) for ln in open(REF / f"{SP}.train.intron.tbl")]
    lo, hi = intron_range(lens)

    out = assemble_param(
        species=SP,
        start_profile=prof_lines("Start_profile"),
        acceptor_profile=prof_lines("Acceptor_profile"),
        donor_profile=prof_lines("Donor_profile"),
        markov_order=5,
        markov_initial=markov_lines("Markov_Initial_probability_matrix"),
        markov_transition=markov_lines("Markov_Transition_probability_matrix"),
        intron_range=format_range(lo, hi),
        intergenic_range="200:Infinity",
    )
    assert "@" not in out
    asm = Param.from_text(out)

    # identical section structure and Markov order
    assert asm.keywords() == refp.keywords()
    assert asm.scalar("Markov_order") == "5"

    # every trained profile reproduces the reference exactly
    for name in ("Start_profile", "Acceptor_profile", "Donor_profile"):
        a = {(p, o): v for p, o, v in asm.profile(name).rows}
        b = {(p, o): v for p, o, v in refp.profile(name).rows}
        assert set(a) == set(b)
        assert max(abs(a[k] - b[k]) for k in b) < 1e-6

    # both Markov matrices identical (string-exact)
    for name in ("Markov_Initial_probability_matrix", "Markov_Transition_probability_matrix"):
        ba, bb = asm._find(name), refp._find(name)
        da = {tuple(ba.raw[i].split()[:3]): ba.raw[i].split()[3] for i in ba.data_line_indices()}
        db = {tuple(bb.raw[i].split()[:3]): bb.raw[i].split()[3] for i in bb.data_line_indices()}
        assert da == db

    # gene-model ranges landed
    gm = "".join(asm._find("General_Gene_Model").raw)
    assert format_range(lo, hi) in gm
    assert "200:Infinity" in gm
