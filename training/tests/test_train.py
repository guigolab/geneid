import pytest

from geneid_train.train import _site_order, build_site_profile, train

from .conftest import ref_dir

REF = ref_dir()
SP = "Xerocrassa_montserratensis"
GENOME = REF.parents[2] / "xgXerMont_curated.no_mt.scrubbed.fa.gz" if REF else None


def test_site_order_standard_case():
    # plenty of donor/acceptor sites, few starts -> order (1, 1, 0)
    assert _site_order(9000, 9000, 800) == (1, 1, 0)


def test_site_order_raises_on_sparse_splice_sites():
    with pytest.raises(NotImplementedError, match="order-0 donor/acceptor"):
        _site_order(500, 9000, 800)


def test_site_order_raises_on_order2_start():
    with pytest.raises(NotImplementedError, match="order-2 start"):
        _site_order(9000, 9000, 6000)


# ---- orchestration validated against the reference --------------------------


def _ref_background():
    freq = {}
    for ln in open(REF / f"{SP}_background.info.freq"):
        p = ln.split()
        if len(p) >= 4:
            freq[(int(p[1]), p[0])] = float(p[3])
    from geneid_train.stats.sites import read_matrix

    dimatrix = read_matrix(REF / f"{SP}_background.info.di-matrix")
    return freq, dimatrix


def _ref_site_seqs(name):
    delim = "\t" if name != "start" else None
    path = REF / f"{SP}.canonical.{name}.tbl"
    return [
        (ln.split(delim)[1] if delim else ln.split()[1]).strip() for ln in open(path)
    ]


@pytest.mark.integration
@pytest.mark.skipif(REF is None, reason="set GENEID_TRAIN_REFDIR to the train_geneid dir")
def test_build_site_profile_matches_reference():
    from geneid_train.core.param import Param

    bg_freq, bg_dimatrix = _ref_background()
    refp = Param.read(REF / f"{SP}.geneid.param")

    def rows_of(lines):
        out = {}
        for ln in lines[2:]:  # skip header + comment
            p = ln.split()
            out[(int(p[0]), p[1])] = float(p[2])
        return out

    cases = [
        ("donor", "Donor_profile", "GT", 1, 1e-3),
        ("acceptor", "Acceptor_profile", "AG", 1, 1e-3),
        ("start", "Start_profile", "ATG", 0, 5e-2),
    ]
    for name, kw, anchor, order, tol in cases:
        lines = build_site_profile(
            _ref_site_seqs(name), site=name, anchor=anchor, order=order,
            bg_freq=bg_freq, bg_dimatrix=bg_dimatrix,
        )
        got = rows_of(lines)
        ref = {(p, o): v for p, o, v in refp.profile(kw).rows}
        assert set(got) == set(ref)
        assert max(abs(got[k] - ref[k]) for k in ref) < tol, name


@pytest.mark.integration
@pytest.mark.skipif(
    REF is None or not (GENOME and GENOME.exists()),
    reason="set GENEID_TRAIN_REFDIR and provide the xgXerMont genome",
)
def test_train_end_to_end_produces_valid_param():
    from geneid_train.core.fasta import read_fasta_subset
    from geneid_train.core.gff import read_gff3
    from geneid_train.core.param import Param
    from geneid_train.prepare.base import (
        build_models,
        collapse_isoforms,
        filter_complete,
        filter_non_overlapping,
    )

    gff = REF.parent / "get_candidates" / "good_candidates.1trans.gff3"
    records = read_gff3(gff)
    models = collapse_isoforms(build_models(records), records)
    # scope to a handful of scaffolds (enough to clear the order-1 site threshold)
    want = {f"SUPER_{i}" for i in range(1, 6)}
    models = [m for m in models if m.seqid in want]
    genome = read_fasta_subset(str(GENOME), {m.seqid for m in models})
    models = filter_non_overlapping(filter_complete(models, genome))

    text = train(models, genome, SP, seed=1)
    assert "@" not in text  # every sentinel filled

    p = Param.from_text(text)
    for kw in ("Start_profile", "Acceptor_profile", "Donor_profile"):
        assert p.profile(kw).rows

    # Markov_order must agree with the matrix oligo widths (regression: a small
    # training set picks order 4, and the matrices must be built at that order,
    # else geneid rejects the file)
    order = int(p.scalar("Markov_order"))
    lines = text.splitlines()
    init_i = lines.index("Markov_Initial_probability_matrix")
    trans_i = lines.index("Markov_Transition_probability_matrix")
    assert len(lines[init_i + 1].split()[0]) == order  # initial oligo = order-mer
    assert len(lines[trans_i + 1].split()[0]) == order + 1  # transition = (order+1)-mer
    assert p.has("Markov_Initial_probability_matrix")
    assert p.has("Markov_Transition_probability_matrix")
    gm = "".join(p._find("General_Gene_Model").raw)
    assert "200:Infinity" in gm  # intergenic range injected
    assert ":500000" in gm  # fixed 500 kb max-intron default (not the p99.9 estimate)
    assert "Intron_length_model" in text  # soft intron-length model emitted
    assert "Intron_length_score_weight\n0.5\n" in text  # penalty ON by default (weight 0.5)
