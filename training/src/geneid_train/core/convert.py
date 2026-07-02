"""Converters from other annotation flavors into canonical GFF3.

The pipeline ingests GFF3 only; these converters are the sanctioned on-ramp for
GFF2 and GTF input. Both are reduced to CDS features carrying a transcript id and
a gene id, then re-emitted as a well-formed gene -> mRNA -> CDS GFF3 hierarchy so
downstream code has a single format to reason about.
"""

from __future__ import annotations

import re
from dataclasses import dataclass

from .gff import GffRecord


@dataclass
class _CdsEntry:
    seqid: str
    source: str
    start: int
    end: int
    score: str
    strand: str
    phase: str
    tx_id: str
    gene_id: str


def _assemble_gff3(entries: list[_CdsEntry]) -> list[GffRecord]:
    """Build gene -> mRNA -> CDS records from flat CDS entries."""
    # preserve first-seen order of transcripts and genes
    tx_order: list[str] = []
    tx_cds: dict[str, list[_CdsEntry]] = {}
    tx_gene: dict[str, str] = {}
    for e in entries:
        if e.tx_id not in tx_cds:
            tx_cds[e.tx_id] = []
            tx_order.append(e.tx_id)
            tx_gene[e.tx_id] = e.gene_id
        tx_cds[e.tx_id].append(e)

    gene_txs: dict[str, list[str]] = {}
    for tx in tx_order:
        gene_txs.setdefault(tx_gene[tx], []).append(tx)

    # keep GFF3 IDs unique: when a gene shares its id with one of its own
    # transcripts (the GFF2 bare-group case), give the gene feature a suffix.
    gene_fid = {
        gene: (f"{gene}.gene" if gene in txs else gene) for gene, txs in gene_txs.items()
    }

    out: list[GffRecord] = []
    emitted_gene: set[str] = set()
    for tx in tx_order:
        gene = tx_gene[tx]
        cds = sorted(tx_cds[tx], key=lambda e: e.start)
        first = cds[0]
        if gene not in emitted_gene:
            g_cds = [c for t in gene_txs[gene] for c in tx_cds[t]]
            out.append(
                GffRecord(
                    first.seqid, first.source, "gene",
                    min(c.start for c in g_cds), max(c.end for c in g_cds),
                    ".", first.strand, ".", {"ID": gene_fid[gene]},
                )
            )
            emitted_gene.add(gene)
        out.append(
            GffRecord(
                first.seqid, first.source, "mRNA",
                min(c.start for c in cds), max(c.end for c in cds),
                ".", first.strand, ".", {"ID": tx, "Parent": gene_fid[gene]},
            )
        )
        for i, c in enumerate(cds, start=1):
            out.append(
                GffRecord(
                    c.seqid, c.source, "CDS", c.start, c.end, c.score, c.strand,
                    c.phase, {"ID": f"{tx}.cds{i}", "Parent": tx},
                )
            )
    return out


def _gff2_group_id(col9: str) -> str:
    """Extract the transcript id from a GFF2 group column (bare token, or a
    quoted token, or a ``key "value"`` attribute)."""
    col9 = col9.strip()
    m = re.search(r'"([^"]+)"', col9)
    if m:
        return m.group(1)
    return col9.split()[0].strip('"') if col9.split() else col9


def gff2_to_gff3(lines: list[str], cds_type: str = "CDS") -> list[GffRecord]:
    entries: list[_CdsEntry] = []
    for line in lines:
        if not line.strip() or line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 8 or f[2] != cds_type:
            continue
        tx = _gff2_group_id(f[8] if len(f) > 8 else "")
        entries.append(
            _CdsEntry(f[0], f[1], int(f[3]), int(f[4]), f[5], f[6], f[7], tx, tx)
        )
    return _assemble_gff3(entries)


def _gtf_attr(col9: str, key: str) -> str | None:
    m = re.search(rf'{re.escape(key)}\s+"([^"]+)"', col9)
    return m.group(1) if m else None


def gtf_to_gff3(lines: list[str], cds_type: str = "CDS") -> list[GffRecord]:
    entries: list[_CdsEntry] = []
    for line in lines:
        if not line.strip() or line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != cds_type:
            continue
        tx = _gtf_attr(f[8], "transcript_id")
        gene = _gtf_attr(f[8], "gene_id") or tx
        if tx is None:
            raise ValueError(f"GTF CDS line missing transcript_id: {line!r}")
        entries.append(
            _CdsEntry(f[0], f[1], int(f[3]), int(f[4]), f[5], f[6], f[7], tx, gene)
        )
    return _assemble_gff3(entries)
