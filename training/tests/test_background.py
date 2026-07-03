import pytest

from geneid_train.stats.background import background_models, from_genome, sample_kmers


def test_sample_kmers_shape_and_alphabet():
    genome = ["ACGT" * 500]  # 2000 bp, all ACGT
    kmers = sample_kmers(genome, k=10, n=20, seed=0)
    assert len(kmers) == 20
    assert all(len(km) == 10 for km in kmers)
    assert all(set(km) <= set("ACGT") for km in kmers)


def test_sample_kmers_is_seeded_reproducible():
    genome = ["ACGTACGTAC" * 200]
    a = sample_kmers(genome, k=8, n=15, seed=42)
    b = sample_kmers(genome, k=8, n=15, seed=42)
    c = sample_kmers(genome, k=8, n=15, seed=43)
    assert a == b
    assert a != c


def test_sample_kmers_rejects_non_acgt():
    # only a short ACGT stretch is long enough for an 8-mer; N regions rejected
    genome = ["N" * 50 + "ACGTACGTACGT" + "N" * 50]
    kmers = sample_kmers(genome, k=8, n=5, seed=0)
    assert all("N" not in km for km in kmers)


def test_sample_kmers_raises_when_too_short():
    with pytest.raises(ValueError):
        sample_kmers(["ACGT"], k=62, n=10)


def test_background_models_shapes():
    genome = ["ACGTACGTAC" * 300]
    freq, dimatrix = from_genome(genome, k=20, n=200, seed=0)
    # order-0 freq keyed (pos, base); order-1 dimatrix keyed (pos, dinuc)
    assert all(len(o) == 1 for _, o in freq)
    assert all(len(o) == 2 for _, o in dimatrix)
    # each position's single-base frequencies sum to ~1
    by_pos: dict[int, float] = {}
    for (pos, _), v in freq.items():
        by_pos[pos] = by_pos.get(pos, 0.0) + v
    assert all(abs(s - 1.0) < 1e-9 for s in by_pos.values())


def test_background_models_direct():
    freq, dimatrix = background_models(["ACGTACGT", "TGCATGCA"])
    assert freq and dimatrix
