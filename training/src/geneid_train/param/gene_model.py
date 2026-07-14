"""GenAmic gene-model rules for UTR-aware prediction.

The bundled ``template.param`` carries a CDS-only gene model. Turning on UTR
prediction (geneid ``-u`` with ``-S``/``-Y`` coverage) needs the gene model to
declare the UTR exon types and their connections, so ``assemble`` (``train
--utr``) and ``retrofit --utr`` swap in the block built here.

Two deliberate choices vs the legacy hand-made human RNA-seq param:
  * the intergenic minimum distance is 0 (not 60): a gene's UTRs extend it, and
    neighbouring genes' UTRs routinely abut or overlap, so 0 is the sane floor;
  * the canonical exon-type name ``UTR_5prime_Internal_Half`` is used throughout
    (the legacy param had a ``UTR_5Internal_Half`` typo that never matched).

The UTR exon types must match geneid's include/geneid.h sUTR* definitions.
"""

from __future__ import annotations

# Human soft intron-length model (log-normal mu, sigma over ln(intron length)),
# from the chr12 training set. A reasonable default to inject when retrofitting a
# param that has no Intron_length_model; species without a fit should keep the
# penalty weight at 0 (see retrofit).
HUMAN_INTRON_LENGTH_MODEL: tuple[float, float] = (7.1788, 1.51411)

UTR_INTERGENIC = "0:Infinity"


def _rule(antecedent: str, consequent: str, drange: str) -> str:
    """One gene-model rule line: antecedent, consequent, distance range. geneid
    tokenizes on whitespace, so the padding is purely for human readability."""
    return f"{antecedent:<56} {consequent:<56} {drange}"


def utr_gene_model_lines(intron_range: str, intergenic: str = UTR_INTERGENIC) -> list[str]:
    """Data lines for a UTR-aware ``General_Gene_Model`` section.

    ``intron_range`` is the ``min:max`` token for the intragenic CDS intron
    distance (a real value when retrofitting, or the ``@INTRON_RANGE@`` sentinel
    when assembling from the template). ``intergenic`` defaults to ``0:Infinity``.
    """
    lines: list[str] = []

    lines.append("# INTRAgenic connections (CDS)")
    lines.append(_rule("First+:Internal+", "Internal+:Terminal+", intron_range))
    lines.append(_rule("Internal-:Terminal-", "First-:Internal-", intron_range))
    for a, b in (
        ("First+", "Intron+"), ("Internal+", "Intron+"),
        ("Intron+", "Internal+"), ("Intron+", "Terminal+"),
        ("Intron-", "First-"), ("Internal-", "Intron-"),
        ("Intron-", "Internal-"), ("Terminal-", "Intron-"),
    ):
        lines.append(_rule(a, b, "1:1"))

    lines.append("# UTR exons: attach to the coding gene and splice among themselves")
    lines.append(_rule("UTR_First_Half+:UTR_5prime_Internal_Half+", "First+:Single+", "1:1"))
    lines.append(_rule("Terminal+:Single+", "UTR_Terminal_Half+:UTR_3prime_Internal_Half+", "1:1"))
    lines.append(_rule("First-:Single-", "UTR_First_Half-:UTR_5prime_Internal_Half-", "1:1"))
    lines.append(_rule("UTR_Terminal_Half-:UTR_3prime_Internal_Half-", "Terminal-:Single-", "1:1"))
    lines.append(_rule("UTR_First+", "UTR_5prime_Intron+:UTR_3prime_Intron+", "1:1"))
    lines.append(_rule("UTR_Internal+", "UTR_5prime_Intron+:UTR_3prime_Intron+", "1:1"))
    lines.append(_rule("UTR_3prime_Internal_Half+", "UTR_3prime_Intron+", "1:1"))
    lines.append(_rule("UTR_5prime_Intron+", "UTR_Internal+:UTR_5prime_Internal_Half+", "1:1"))
    lines.append(_rule("UTR_3prime_Intron+", "UTR_Terminal+", "1:1"))
    lines.append(_rule("UTR_5prime_Intron-:UTR_3prime_Intron-", "UTR_First-", "1:1"))
    lines.append(_rule("UTR_5prime_Intron-:UTR_3prime_Intron-", "UTR_Internal-", "1:1"))
    lines.append(_rule("UTR_3prime_Intron-", "UTR_3prime_Internal_Half-", "1:1"))
    lines.append(_rule("UTR_Internal-:UTR_5prime_Internal_Half-", "UTR_5prime_Intron-", "1:1"))
    lines.append(_rule("UTR_Terminal-", "UTR_3prime_Intron-", "1:1"))

    lines.append("# External features (promoter / poly-A signal)")
    lines.append(_rule("Promoter+", "First+:Single+", "50:4000"))
    lines.append(_rule("First-:Single-", "Promoter-", "50:4000"))
    lines.append(_rule("Terminal+:Single+", "aataaa+", "50:4000"))
    lines.append(_rule("aataaa-", "Terminal-:Single-", "50:4000"))

    lines.append("# INTERgenic connections (UTR-aware; min 0 so neighbouring UTRs may overlap)")
    fwd_end = "aataaa+:Terminal+:Single+:UTR_Terminal+:UTR_Terminal_Half+"
    fwd_start = "Single+:First+:Promoter+:UTR_First+:UTR_First_Half+"
    rev_end = "Promoter-:First-:Single-:UTR_First-:UTR_First_Half-"
    rev_start = "Single-:Terminal-:aataaa-:UTR_Terminal-:UTR_Terminal_Half-"
    lines.append(_rule(fwd_end, fwd_start, intergenic))
    lines.append(_rule(fwd_end, rev_start, intergenic))
    lines.append(_rule(rev_end, fwd_start, intergenic))
    lines.append(_rule(rev_end, rev_start, intergenic))

    lines.append("# BEGINNING and END of prediction")
    lines.append(_rule(
        "Begin+",
        "First+:Internal+:Terminal+:Single+:UTR_First+:UTR_First_Half+:UTR_5prime_Internal_Half+",
        "0:Infinity"))
    lines.append(_rule(
        "Begin-",
        "First-:Internal-:Terminal-:Single-:UTR_Terminal_Half-:UTR_3prime_Internal_Half-:UTR_Terminal-",
        "0:Infinity"))
    lines.append(_rule(
        "First+:Internal+:Terminal+:Single+:UTR_Terminal_Half+:UTR_Terminal+:UTR_3prime_Internal_Half+",
        "End+", "0:Infinity"))
    lines.append(_rule(
        "First-:Internal-:Terminal-:Single-:UTR_First-:UTR_First_Half-:UTR_5prime_Internal_Half-",
        "End-", "0:Infinity"))
    return lines


def extract_intron_range(param) -> str:  # noqa: ANN001 (Param, avoid import cycle)
    """Return the intragenic CDS intron ``min:max`` from a param's gene model
    (the ``First+:Internal+ -> Internal+:Terminal+`` rule). Raises if the gene
    model does not have the expected standard structure."""
    for line in param.data_lines("General_Gene_Model"):
        toks = line.split()
        if len(toks) >= 3 and toks[0] == "First+:Internal+":
            return toks[2]
    raise ValueError(
        "gene model has no 'First+:Internal+' intragenic rule; cannot determine "
        "the intron range to preserve (non-standard gene model?)"
    )
