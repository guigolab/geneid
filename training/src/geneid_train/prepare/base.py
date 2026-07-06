"""Validated gene-model intermediate representation (IR) and shared filters.

Both training-set front-ends (BUSCO, RNA-seq/TransDecoder) converge on the
``GeneModel`` IR defined here, so all filtering — completeness, minimum protein
length, non-overlap, flank extraction — lives in one place and is exercised by one
test suite.

Coordinates are 1-based inclusive throughout, matching geneid/GFF conventions.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from dataclasses import dataclass

from ..core.gff import GffRecord
from ..core.seq import STANDARD, GeneticCode, revcomp


@dataclass
class Exon:
    start: int  # 1-based inclusive
    end: int  # 1-based inclusive
    phase: str = "."


@dataclass
class GeneModel:
    gene_id: str
    seqid: str
    strand: str
    exons: list[Exon]  # stored sorted by ascending genomic start

    @property
    def start(self) -> int:
        return self.exons[0].start

    @property
    def end(self) -> int:
        return self.exons[-1].end

    @property
    def is_multiexonic(self) -> bool:
        return len(self.exons) > 1

    @property
    def coding_length(self) -> int:
        """Total CDS length in nt (sum of exon lengths) — no genome needed."""
        return sum(e.end - e.start + 1 for e in self.exons)

    def cds(self, genome: Mapping[str, str]) -> str:
        """Spliced CDS in 5'->3' translation orientation."""
        chrom = genome[self.seqid]
        pieces = [chrom[e.start - 1 : e.end] for e in self.exons]
        joined = "".join(pieces)
        return revcomp(joined) if self.strand == "-" else joined

    def protein(self, genome: Mapping[str, str], code: GeneticCode = STANDARD) -> str:
        return code.translate(self.cds(genome))

    def introns(self) -> list[tuple[int, int]]:
        """Genomic (start, end) of each intron, ascending. Empty if single-exon."""
        out = []
        for prev, nxt in zip(self.exons, self.exons[1:], strict=False):
            out.append((prev.end + 1, nxt.start - 1))
        return out

    def intron_seqs(self, genome: Mapping[str, str]) -> list[str]:
        """Intron sequences in transcription orientation (donor at the 5' end)."""
        chrom = genome[self.seqid]
        seqs = [chrom[s - 1 : e] for s, e in self.introns()]
        if self.strand == "-":
            seqs = [revcomp(x) for x in reversed(seqs)]
        return seqs

    def is_complete(self, genome: Mapping[str, str], code: GeneticCode = STANDARD) -> bool:
        from ..core.seq import is_complete_cds

        return is_complete_cds(self.cds(genome), code)


def build_models(records: Iterable[GffRecord], cds_type: str = "CDS") -> list[GeneModel]:
    """Group CDS features into GeneModels by their GFF3 ``Parent`` (transcript id).

    Records whose ``type`` differs from ``cds_type`` are ignored. Raises if a
    transcript has features on mixed seqids or strands, or a CDS lacks a Parent.
    """
    grouped: dict[str, list[GffRecord]] = {}
    for rec in records:
        if rec.type != cds_type:
            continue
        if not rec.parents:
            raise ValueError(f"CDS {rec.id or rec.to_line()!r} has no Parent attribute")
        grouped.setdefault(rec.parents[0], []).append(rec)

    models: list[GeneModel] = []
    for tx_id, recs in grouped.items():
        seqids = {r.seqid for r in recs}
        strands = {r.strand for r in recs}
        if len(seqids) != 1:
            raise ValueError(f"transcript {tx_id!r} spans multiple seqids: {sorted(seqids)}")
        if len(strands) != 1:
            raise ValueError(f"transcript {tx_id!r} has mixed strands: {sorted(strands)}")
        recs_sorted = sorted(recs, key=lambda r: r.start)
        exons = [Exon(r.start, r.end, r.phase) for r in recs_sorted]
        models.append(
            GeneModel(
                gene_id=tx_id,
                seqid=recs_sorted[0].seqid,
                strand=recs_sorted[0].strand,
                exons=exons,
            )
        )
    return models


def gene_of_transcript(records: Iterable[GffRecord]) -> dict[str, str]:
    """Map transcript id -> gene id from mRNA/transcript records' Parent."""
    out: dict[str, str] = {}
    for r in records:
        if r.type in ("mRNA", "transcript") and r.id:
            out[r.id] = r.parents[0] if r.parents else r.id
    return out


def collapse_isoforms(
    models: Iterable[GeneModel], records: Iterable[GffRecord]
) -> list[GeneModel]:
    """Keep one representative model (longest CDS) per gene.

    Must run before ``filter_non_overlapping``: sibling isoforms overlap and would
    otherwise eliminate each other. Genes are resolved via transcript->gene Parent
    links; a model whose transcript has no gene record is treated as its own gene.
    """
    tx2gene = gene_of_transcript(records)
    by_gene: dict[str, list[GeneModel]] = {}
    for m in models:
        gene = tx2gene.get(m.gene_id, m.gene_id)
        by_gene.setdefault(gene, []).append(m)
    return [max(ms, key=lambda x: x.coding_length) for ms in by_gene.values()]


# ---- filters --------------------------------------------------------------


def filter_complete(
    models: Iterable[GeneModel], genome: Mapping[str, str], code: GeneticCode = STANDARD
) -> list[GeneModel]:
    return [m for m in models if m.is_complete(genome, code)]


def filter_min_protein(
    models: Iterable[GeneModel],
    genome: Mapping[str, str],
    min_aa: int,
    code: GeneticCode = STANDARD,
) -> list[GeneModel]:
    kept = []
    for m in models:
        prot = m.protein(genome, code).rstrip("*")
        if len(prot) >= min_aa:
            kept.append(m)
    return kept


def filter_non_overlapping(models: Iterable[GeneModel]) -> list[GeneModel]:
    """Drop every model whose genomic span overlaps any other model on the same
    seqid, regardless of strand (a training-set requirement). Both members of an
    overlapping pair are removed.
    """
    mlist = list(models)
    overlapping: set[int] = set()
    by_seq: dict[str, list[int]] = {}
    for i, m in enumerate(mlist):
        by_seq.setdefault(m.seqid, []).append(i)
    for idxs in by_seq.values():
        idxs.sort(key=lambda i: mlist[i].start)
        for a in range(len(idxs)):
            ma = mlist[idxs[a]]
            for b in range(a + 1, len(idxs)):
                mb = mlist[idxs[b]]
                if mb.start > ma.end:
                    break  # no further overlaps for ma (sorted by start)
                overlapping.add(idxs[a])
                overlapping.add(idxs[b])
    return [m for i, m in enumerate(mlist) if i not in overlapping]


@dataclass
class Locus:
    """A gene placed on a sub-sequence with genomic flanks, coordinates shifted
    so the locus starts at 1. This is the per-gene training unit geneid sees."""

    gene_id: str
    seq: str
    records: list[GffRecord]


def extract_locus(
    model: GeneModel, genome: Mapping[str, str], flank: int, source: str = "geneid_train"
) -> Locus:
    chrom = genome[model.seqid]
    lo = max(0, model.start - 1 - flank)  # 0-based slice start
    hi = min(len(chrom), model.end + flank)  # 0-based slice end (exclusive)
    seq = chrom[lo:hi]
    shift = lo  # subtract from 1-based coords: new = old - lo
    recs = [
        GffRecord(
            seqid=model.gene_id,
            source=source,
            type="CDS",
            start=e.start - shift,
            end=e.end - shift,
            score=".",
            strand=model.strand,
            phase=e.phase,
            attributes={"ID": f"{model.gene_id}.cds{i}", "Parent": model.gene_id},
        )
        for i, e in enumerate(model.exons, start=1)
    ]
    return Locus(gene_id=model.gene_id, seq=seq, records=recs)
