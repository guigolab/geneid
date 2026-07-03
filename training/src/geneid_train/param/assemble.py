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
) -> str:
    """Build a full geneid param. Profile/Markov args are geneid-format data
    lines (from ``stats.sites.format_profile`` / ``stats.coding.format_markov_matrix``);
    ``intron_range``/``intergenic_range`` are ``min:max`` tokens."""
    p = Param.from_text(template if template is not None else load_template())

    p.replace_block_data("Start_profile", start_profile)
    p.replace_block_data("Acceptor_profile", acceptor_profile)
    p.replace_block_data("Donor_profile", donor_profile)

    p.set_scalar("Markov_order", markov_order)
    p.replace_block_data("Markov_Initial_probability_matrix", markov_initial)
    p.replace_block_data("Markov_Transition_probability_matrix", markov_transition)

    text = p.to_text()
    text = text.replace("@INTRON_RANGE@", intron_range)
    text = text.replace("@INTERGENIC_RANGE@", intergenic_range)
    text = text.replace("@SPECIES@", species)
    return text
