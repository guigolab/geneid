"""RNA-seq / TransDecoder front-end.

Ingests a TransDecoder (or PASA) genome GFF3 of predicted ORFs and produces the
validated gene-model IR, one representative (longest-CDS) model per gene.

Provenance/quality filtering — removing repeat-overlapping candidates and the
UniProt cross-check (DIAMOND, >=90% coverage) — is done upstream in the annotation
pipeline (e.g. the CNAG `get_candidates` step yields `good_candidates.1trans.gff3`).
This adapter therefore does not re-implement the UniProt check; it consumes the
already-filtered candidate GFF3. Verified against real xgXerMont data: all such
candidates are complete ORFs (ATG start, terminal stop inside the CDS coordinates,
in-frame, no internal stops).
"""

from __future__ import annotations

from pathlib import Path

from ..core.gff import read_gff3
from .base import GeneModel, build_models, collapse_isoforms


def load(gff3_path: str | Path, cds_type: str = "CDS") -> list[GeneModel]:
    """Read a TransDecoder/PASA genome GFF3 -> one GeneModel per gene."""
    records = read_gff3(gff3_path)
    models = build_models(records, cds_type=cds_type)
    return collapse_isoforms(models, records)
