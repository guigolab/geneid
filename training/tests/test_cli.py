from geneid_train.cli import build_parser
from geneid_train.stats.genemodel import (
    DEFAULT_INTRON_LENGTH_WEIGHT,
    DEFAULT_MAX_INTRON,
)

BASE = ["train", "--gff", "a.gff3", "--fastas", "g.fa",
        "--species", "sp", "--output", "out.param"]


def test_train_max_intron_defaults_to_fixed_500kb():
    # The gene-model max intron defaults to a genome-independent 500 kb safety bound;
    # the soft intron-length penalty (not this cap) does the real length tuning.
    args = build_parser().parse_args(BASE)
    assert args.max_intron == DEFAULT_MAX_INTRON == 500_000


def test_train_max_intron_override():
    args = build_parser().parse_args([*BASE, "--max-intron", "1000000"])
    assert args.max_intron == 1_000_000


def test_train_intron_length_weight_defaults_to_half():
    # The soft intron-length penalty is ON by default (weight 0.5).
    args = build_parser().parse_args(BASE)
    assert args.intron_length_weight == DEFAULT_INTRON_LENGTH_WEIGHT == 0.5


def test_train_intron_length_weight_override():
    args = build_parser().parse_args([*BASE, "--intron-length-weight", "0"])
    assert args.intron_length_weight == 0.0
