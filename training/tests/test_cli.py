from geneid_train.cli import build_parser
from geneid_train.stats.genemodel import DEFAULT_MAX_INTRON

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
