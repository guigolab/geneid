"""Assemble a complete geneid parameter file from trained components.

Everything geneid needs that is *not* estimated from the training set — the
score factors and weights, the RSS/cutoff defaults, the fixed ``Stop_profile``,
and the GenAmic gene-model rules — lives in a bundled single-isochore template
(``data/template.param``, a real geneid param with the trained sections stubbed
out and the variable ranges marked with sentinels). Assembly swaps in the
trained pieces: the Start/Acceptor/Donor profiles, the two coding Markov
matrices, and the intron/intergenic gene-model ranges.

Byte-identical reproduction of a legacy param is a non-goal (log values differ
in the low decimals); the contract is a valid, correctly structured param whose
trained sections carry our estimates.
"""

from __future__ import annotations

from importlib.resources import files

from ..core.param import Param
from ..stats.genemodel import DEFAULT_INTRON_LENGTH_WEIGHT
from .gene_model import utr_gene_model_lines
from .u12 import U12Sections


def load_template() -> str:
    """Return the bundled template parameter file as text."""
    return files("geneid_train").joinpath("data/template.param").read_text()


def assemble_param(
    *,
    species: str,
    start_profile: list[str],
    acceptor_profile: list[str],
    donor_profile: list[str],
    markov_order: int,
    markov_initial: list[str],
    markov_transition: list[str],
    intron_range: str,
    intergenic_range: str,
    template: str | None = None,
    u12: U12Sections | None = None,
    u12_splice_thresh: float = 9.0,
    u12_exon_thresh: float = 8.0,
    intron_length_model: tuple[float, float] | None = None,
    intron_length_weight: float = DEFAULT_INTRON_LENGTH_WEIGHT,
    utr: bool = False,
) -> str:
    """Build a full geneid param. Profile/Markov args are geneid-format data
    lines (from ``stats.sites.format_profile`` / ``stats.coding.format_markov_matrix``);
    ``intron_range``/``intergenic_range`` are ``min:max`` tokens.

    When ``u12`` is given, the minor-spliceosome profile trio is spliced in before
    the required acceptor/donor profiles so geneid enables U12 splice prediction,
    and the ``U12_Splice_Score_Threshold`` / ``U12_Exon_Score_Threshold`` gates are
    emitted. Those thresholds are essential: without them geneid uses its −1000
    default, which accepts almost any GT-AG intron as U12 (≈6.5% mislabeled on
    xgXerMont vs a plausible <0.5%). The default (9, human value; drosophila uses
    10) is deliberately conservative — U12 is hard to train without a large
    genome-specific set, so we favour precision and let homology/RNA-seq recover
    false negatives later. A SUPER_1 sweep: threshold 8→0.41%, 9≈0.1-0.2%, 10→0.03%.
    """
    p = Param.from_text(template if template is not None else load_template())

    p.replace_block_data("Start_profile", start_profile)
    p.replace_block_data("Acceptor_profile", acceptor_profile)
    p.replace_block_data("Donor_profile", donor_profile)

    p.set_scalar("Markov_order", markov_order)
    p.replace_block_data("Markov_Initial_probability_matrix", markov_initial)
    p.replace_block_data("Markov_Transition_probability_matrix", markov_transition)

    if utr:
        # Swap the CDS-only gene model for a UTR-aware one; the intragenic CDS
        # intron distance keeps the @INTRON_RANGE@ sentinel (filled below), the
        # intergenic minimum drops to 0 (neighbouring UTRs may abut/overlap).
        p.replace_block_data("General_Gene_Model", utr_gene_model_lines("@INTRON_RANGE@"))

    if u12 is not None:
        # order matters only in that each optional profile must precede the fixed
        # profile it extends (geneid reads optional profiles until it hits the
        # required Acceptor_profile / Donor_profile — see ReadProfileSpliceSites)
        p.insert_text_before("Acceptor_profile", u12.acceptor_side)
        p.insert_text_before("Donor_profile", u12.donor_side)
        # the U12 acceptance-score gates (read in the optional-scalar block before
        # Exon_weights) — without them geneid mass-mislabels U2 GT-AG as U12
        p.insert_text_before(
            "Exon_weights",
            f"U12_Splice_Score_Threshold\n{u12_splice_thresh:g}\n"
            f"U12_Exon_Score_Threshold\n{u12_exon_thresh:g}\n",
        )

    if intron_length_model is not None:
        # Soft intron-length penalty: the log-normal (mu, sigma) over ln(length)
        # plus its weight (lambda). The weight defaults to DEFAULT_INTRON_LENGTH_WEIGHT
        # (0.5 = ON); pass 0 to carry the model inert. geneid reads both in the
        # optional-scalar block before Exon_weights.
        mu, sigma = intron_length_model
        p.insert_text_before(
            "Exon_weights",
            f"Intron_length_model\n{mu:.6g} {sigma:.6g}\n"
            f"Intron_length_score_weight\n{intron_length_weight:g}\n",
        )

    text = p.to_text()
    text = text.replace("@INTRON_RANGE@", intron_range)
    text = text.replace("@INTERGENIC_RANGE@", intergenic_range)
    text = text.replace("@SPECIES@", species)
    return text
