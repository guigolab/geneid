"""Retrofit an existing geneid parameter file with newer capabilities.

Older params (and anything the trainer emitted before these features) are
CDS-only. ``retrofit`` injects, non-destructively and idempotently, the pieces
needed for the modern RNA-seq workflow onto an existing param:

  * ``--utr``: a UTR-aware gene model (so ``geneid -u`` with ``-S``/``-Y``
    coverage predicts UTRs), reusing the param's own trained intron range and
    dropping the intergenic minimum to 0;
  * ``--intron-length``: the soft intron-length penalty section (log-normal
    mu/sigma + weight), defaulting the weight to 0 (inert) so it never changes
    predictions until deliberately enabled;
  * ``--u12``: the bundled pan-taxon U12 (minor-spliceosome) profile trio plus
    the U12 acceptance-score gates, so ``geneid -U`` predicts U12 introns. The
    profiles carry their own geneid geometry (offset/length, from the reference
    U12 param -- see param/u12.py), so they drop into any param's Donor/Acceptor
    slot; nothing genome-specific is re-registered here.

Only single-isochore params are supported for now (the trainer produces those;
multi-isochore params like human.rnaseq.param are already UTR-capable).
"""

from __future__ import annotations

from ..core.param import Param
from .gene_model import HUMAN_INTRON_LENGTH_MODEL, extract_intron_range, utr_gene_model_lines
from .u12 import load_bundled_u12


def _has_utr(param: Param) -> bool:
    return any("UTR_" in line for line in param.data_lines("General_Gene_Model"))


def retrofit_param(
    text: str,
    *,
    utr: bool = False,
    intron_length: tuple[float, float] | None = None,
    intron_length_weight: float = 0.0,
    u12: bool = False,
    u12_splice_thresh: float = 9.0,
    u12_exon_thresh: float = 8.0,
) -> str:
    """Return the retrofitted param text. ``intron_length`` may be
    ``HUMAN_INTRON_LENGTH_MODEL`` or an explicit ``(mu, sigma)``."""
    p = Param.from_text(text)

    if p.has("number_of_isochores") and p.num_isochores > 1:
        raise NotImplementedError(
            f"retrofit supports single-isochore params only (found "
            f"{p.num_isochores} isochores); multi-isochore support is a follow-up"
        )

    if utr:
        if _has_utr(p):
            print("retrofit: gene model already has UTR rules; leaving it unchanged")
        else:
            intron_range = extract_intron_range(p)
            p.replace_block_data("General_Gene_Model", utr_gene_model_lines(intron_range))

    if intron_length is not None:
        if p.has("Intron_length_model"):
            print("retrofit: param already has an Intron_length_model; leaving it unchanged")
        else:
            mu, sigma = intron_length
            src = "human default" if intron_length == HUMAN_INTRON_LENGTH_MODEL else "supplied"
            p.insert_text_before(
                "Exon_weights",
                "# Soft intron-length penalty: log-normal (mu sigma) over ln(intron length)\n"
                f"# + weight (lambda). mu/sigma below are the {src} values; for another\n"
                "# species either retrain (geneid-train optimize --tune-intron-length) or\n"
                "# keep the weight at 0 to use the hard gene-model max-distance gate instead.\n"
                "Intron_length_model\n"
                f"{mu:.6g} {sigma:.6g}\n"
                "Intron_length_score_weight\n"
                f"{intron_length_weight:g}\n",
            )

    if u12:
        # geneid enables a U12 subtype only when its full profile trio is present
        # (branch + donor + acceptor); the bundle carries GT-AG and AT-AC.
        if p.has("U12gtag_Donor_profile"):
            print("retrofit: param already has U12 profiles; leaving them unchanged")
        else:
            sections = load_bundled_u12()
            # optional profiles must precede the required profile they extend
            p.insert_text_before("Acceptor_profile", sections.acceptor_side)
            p.insert_text_before("Donor_profile", sections.donor_side)
            # the U12 acceptance gates (without them geneid's -1000 default
            # mass-mislabels U2 GT-AG as U12); read before Exon_weights
            if not p.has("U12_Splice_Score_Threshold"):
                p.insert_text_before(
                    "Exon_weights",
                    f"U12_Splice_Score_Threshold\n{u12_splice_thresh:g}\n"
                    f"U12_Exon_Score_Threshold\n{u12_exon_thresh:g}\n",
                )

    return p.to_text()
