"""Prediction accuracy scorer (replaces the compiled ``evaluation`` tool).

Compares a geneid prediction GFF against an annotation GFF and reports
sensitivity/specificity at the nucleotide, exon and gene levels, matching the
definitions of Enrique Blanco's ``evaluation`` C tool (the metric the parameter
optimiser maximises is the exon-level ``SNSP = (SNe + SPe)/2``).

Both inputs are read per locus in the geneid "gp" convention: the *annotation*
file's first line of each locus is an info line whose 5th field is the sequence
length (and is not itself counted as an exon — this mirrors the C tool's
``ExtractInfo``); every other line is a typed CDS exon (First/Internal/Terminal/
Single) grouped into genes by column 9. Stats are accumulated across loci and the
totals reported (the ``-t`` total line the optimiser reads).
"""

from __future__ import annotations

import math
from collections import defaultdict
from dataclasses import dataclass, field


@dataclass
class Exon:
    start: int
    end: int
    strand: str
    group: str


@dataclass
class Locus:
    name: str
    length: int
    exons: list[Exon]


def _overlap_nt(a: list[Exon], b: list[Exon]) -> int:
    """Total nucleotides covered by both exon sets (same strand). Exon sets are
    internally non-overlapping CDS, so this equals the C tool's computeTP."""
    tp = 0
    for p in a:
        for r in b:
            lo = max(p.start, r.start)
            hi = min(p.end, r.end)
            if hi >= lo:
                tp += hi - lo + 1
    return tp


def _by_strand(exons: list[Exon]) -> dict[str, list[Exon]]:
    out: dict[str, list[Exon]] = defaultdict(list)
    for e in exons:
        out[e.strand].append(e)
    return out


def _exact_matches(pred: list[Exon], real: list[Exon]) -> int:
    """Count predicted exons whose (start, end) exactly match a real exon on the
    same strand (multiset intersection)."""
    real_keys: dict[tuple[int, int, str], int] = defaultdict(int)
    for r in real:
        real_keys[(r.start, r.end, r.strand)] += 1
    tpe = 0
    for p in pred:
        k = (p.start, p.end, p.strand)
        if real_keys.get(k, 0) > 0:
            real_keys[k] -= 1
            tpe += 1
    return tpe


def _no_overlap_count(query: list[Exon], other: list[Exon]) -> int:
    """Number of query exons with no overlapping exon in ``other`` (same strand)."""
    by_strand = _by_strand(other)
    n = 0
    for q in query:
        hits = by_strand.get(q.strand, [])
        if not any(min(q.end, o.end) >= max(q.start, o.start) for o in hits):
            n += 1
    return n


def _genes(exons: list[Exon]) -> dict[tuple[str, str], list[Exon]]:
    """Group exons into genes by (group id, strand), each sorted by start."""
    genes: dict[tuple[str, str], list[Exon]] = defaultdict(list)
    for e in exons:
        genes[(e.group, e.strand)].append(e)
    for g in genes.values():
        g.sort(key=lambda e: e.start)
    return genes


def _gene_matches(pred: list[Exon], real: list[Exon]) -> int:
    """Count predicted genes whose full exon structure (same count, identical
    boundaries) matches a real gene on the same strand."""
    real_structs: dict[tuple, int] = defaultdict(int)
    for (_, strand), exons in _genes(real).items():
        key = (strand, tuple((e.start, e.end) for e in exons))
        real_structs[key] += 1
    tpg = 0
    for (_, strand), exons in _genes(pred).items():
        key = (strand, tuple((e.start, e.end) for e in exons))
        if real_structs.get(key, 0) > 0:
            real_structs[key] -= 1
            tpg += 1
    return tpg


@dataclass
class Totals:
    tp: int = 0
    cds_real: int = 0
    cds_pred: int = 0
    length: int = 0
    tpe: int = 0
    exr: int = 0
    exp: int = 0
    me: int = 0
    we: int = 0
    tpg: int = 0
    ger: int = 0
    gep: int = 0


@dataclass
class Accuracy:
    """Total-level accuracy across all loci."""

    sn: float
    sp: float
    cc: float
    sne: float
    spe: float
    snsp: float
    sng: float
    spg: float
    snspg: float
    ra_me: float
    ra_we: float
    totals: Totals = field(repr=False, default_factory=Totals)


def _safe_div(a: float, b: float) -> float:
    return a / b if b else 0.0


def accumulate(pred_exons: list[Exon], real_exons: list[Exon], length: int, t: Totals) -> None:
    """Add one locus's counts into the running totals."""
    pred_s, real_s = _by_strand(pred_exons), _by_strand(real_exons)
    for strand in ("+", "-"):
        t.tp += _overlap_nt(pred_s.get(strand, []), real_s.get(strand, []))
        t.tpe += _exact_matches(pred_s.get(strand, []), real_s.get(strand, []))
    t.cds_real += sum(e.end - e.start + 1 for e in real_exons)
    t.cds_pred += sum(e.end - e.start + 1 for e in pred_exons)
    t.length += length
    t.exr += len(real_exons)
    t.exp += len(pred_exons)
    t.me += _no_overlap_count(real_exons, pred_exons)
    t.we += _no_overlap_count(pred_exons, real_exons)
    t.tpg += _gene_matches(pred_exons, real_exons)
    t.ger += len(_genes(real_exons))
    t.gep += len(_genes(pred_exons))


def finalize(t: Totals) -> Accuracy:
    """Compute total-level SN/SP/CC and exon/gene metrics from accumulated counts."""
    sn = _safe_div(t.tp, t.cds_real)
    sp = _safe_div(t.tp, t.cds_pred)
    # correlation coefficient (needs the non-coding complement via sequence length)
    rn = t.length - t.cds_real
    pn = t.length - t.cds_pred
    fp = t.cds_pred - t.tp
    fn = t.cds_real - t.tp
    tn = rn - fp
    denom = t.cds_real * rn * t.cds_pred * pn
    cc = ((t.tp * tn) - (fn * fp)) / math.sqrt(denom) if denom > 0 else 0.0
    sne = _safe_div(t.tpe, t.exr)
    spe = _safe_div(t.tpe, t.exp)
    sng = _safe_div(t.tpg, t.ger)
    spg = _safe_div(t.tpg, t.gep)
    return Accuracy(
        sn=sn, sp=sp, cc=cc,
        sne=sne, spe=spe, snsp=(sne + spe) / 2,
        sng=sng, spg=spg, snspg=(sng + spg) / 2,
        ra_me=_safe_div(t.me, t.exr), ra_we=_safe_div(t.we, t.exp),
        totals=t,
    )


_TYPED_EXONS = {"First", "Internal", "Terminal", "Single"}


def read_annotation_gff(path) -> list[Locus]:
    """Parse the annotation (real) GFF in gp convention: the first line of each
    locus is the info line (field 5 = sequence length, not an exon); the rest are
    typed CDS exons grouped by column 9."""
    loci: list[Locus] = []
    current: Locus | None = None
    for raw in open(path):
        line = raw.rstrip("\n")
        if not line or line.startswith("#"):
            continue
        f = line.split("\t")
        locus = f[0]
        if current is None or current.name != locus:
            # first line of a new locus = info line (length in field 5)
            current = Locus(name=locus, length=int(f[4]), exons=[])
            loci.append(current)
            continue
        if f[2] in _TYPED_EXONS:
            current.exons.append(Exon(int(f[3]), int(f[4]), f[6], f[8]))
    return loci


def read_prediction_gff(path) -> dict[str, list[Exon]]:
    """Parse the prediction GFF into ``{locus: [typed CDS exons]}`` (no info line)."""
    out: dict[str, list[Exon]] = defaultdict(list)
    for raw in open(path):
        line = raw.rstrip("\n")
        if not line or line.startswith("#"):
            continue
        f = line.split("\t")
        if len(f) >= 9 and f[2] in _TYPED_EXONS:
            out[f[0]].append(Exon(int(f[3]), int(f[4]), f[6], f[8]))
    return out


def evaluate_files(prediction_gff, annotation_gff) -> Accuracy:
    """Score a prediction GFF against an annotation GFF (total-level accuracy)."""
    real_loci = read_annotation_gff(annotation_gff)
    pred = read_prediction_gff(prediction_gff)
    t = Totals()
    for locus in real_loci:
        accumulate(pred.get(locus.name, []), locus.exons, locus.length, t)
    return finalize(t)
