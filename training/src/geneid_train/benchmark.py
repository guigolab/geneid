"""Benchmark geneid predictions against a GENCODE reference (MANE Select).

Turns "does this RNA-seq mode help?" into SN/SP numbers. Two pieces:

  * ``gencode_mane_gp`` converts a GENCODE GFF3 into geneid's gp-format CDS
    truth for one sequence -- the MANE Select transcript of each protein-coding
    gene, its CDS exons typed First/Internal/Terminal/Single (transcription
    order) and grouped by transcript, with the per-locus info line evaluate.py
    expects (field 5 = sequence length).
  * ``run_benchmark`` runs geneid over the sequence in several evidence modes
    (ab initio, +coverage, +UTR, +BAM, ...) and scores each with
    :func:`geneid_train.evaluate.evaluate_files`, so the deltas show which mode
    (and knob) helps or hurts.

Mode A (refine ab-initio) today; the same harness will score Mode B
(evidence-gated reconstruction) once that lands.
"""

from __future__ import annotations

import gzip
import subprocess
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

from .evaluate import (
    Accuracy,
    Exon,
    Totals,
    _by_strand,
    _exact_matches,
    _gene_matches,
    _no_overlap_count,
    _overlap_nt,
    accumulate,
    read_annotation_gff,
    read_prediction_gff,
)

_TYPED = ("First", "Internal", "Terminal", "Single")


def _merge(exons: list[Exon]) -> dict[str, list[Exon]]:
    """Per-strand union of exon intervals, so overlapping isoforms don't
    double-count in the nucleotide-precision overlap (which assumes a
    non-overlapping reference)."""
    out: dict[str, list[Exon]] = {}
    for strand, es in _by_strand(exons).items():
        merged: list[list[int]] = []
        for s, e in sorted((x.start, x.end) for x in es):
            if merged and s <= merged[-1][1]:
                merged[-1][1] = max(merged[-1][1], e)
            else:
                merged.append([s, e])
        out[strand] = [Exon(s, e, strand, "m") for s, e in merged]
    return out


def evaluate_asymmetric(pred_gff: str, recall_gp: str, precision_gp: str) -> Accuracy:
    """Recall (SN) scored against the canonical (MANE) reference; precision (SP)
    scored against all annotated CDS isoforms. Denominators are the prediction.
    A predicted feature counts toward precision if it matches *any* isoform, so
    predicting a real non-canonical transcript is not penalized as a false call."""
    pred = read_prediction_gff(pred_gff)
    prec = {loc.name: loc for loc in read_annotation_gff(precision_gp)}

    r = Totals()          # recall side (vs canonical): tp/tpe/tpg + pred denominators
    tp_u = tpe_a = tpg_a = we_a = 0   # precision numerators (vs all isoforms)
    for locus in read_annotation_gff(recall_gp):
        p = pred.get(locus.name, [])
        accumulate(p, locus.exons, locus.length, r)
        allx = prec[locus.name].exons if locus.name in prec else []
        pred_s, all_s, all_m = _by_strand(p), _by_strand(allx), _merge(allx)
        for st in ("+", "-"):
            tp_u += _overlap_nt(pred_s.get(st, []), all_m.get(st, []))
            tpe_a += _exact_matches(pred_s.get(st, []), all_s.get(st, []))
        tpg_a += _gene_matches(p, allx)
        we_a += _no_overlap_count(p, allx)

    sn, sp = _sd(r.tp, r.cds_real), _sd(tp_u, r.cds_pred)
    sne, spe = _sd(r.tpe, r.exr), _sd(tpe_a, r.exp)
    sng, spg = _sd(r.tpg, r.ger), _sd(tpg_a, r.gep)
    return Accuracy(
        sn=sn, sp=sp, cc=0.0,
        sne=sne, spe=spe, snsp=(sne + spe) / 2,
        sng=sng, spg=spg, snspg=(sng + spg) / 2,
        ra_me=_sd(r.me, r.exr), ra_we=_sd(we_a, r.exp), totals=r,
    )


def _sd(a: float, b: float) -> float:
    return a / b if b else 0.0


def _open(path: str | Path):
    p = str(path)
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p)


def _attrs(col9: str) -> dict[str, str]:
    return dict(kv.split("=", 1) for kv in col9.split(";") if "=" in kv)


def fasta_seqlen(fasta: str | Path, seqid: str) -> int:
    """Length of the ``seqid`` record in a FASTA (the value geneid sees, so the
    gp info line matches). The header token is the id up to first whitespace."""
    n, active = 0, False
    with _open(fasta) as fh:
        for line in fh:
            if line.startswith(">"):
                if active:
                    break
                active = line[1:].split()[0] == seqid
            elif active:
                n += len(line.strip())
    if not active and n == 0:
        raise ValueError(f"sequence {seqid!r} not found in {fasta}")
    return n


def gencode_cds_gp(gff3_path: str | Path, seqid: str, seqlen: int, *, canonical_only: bool) -> str:
    """gp-format CDS truth for ``seqid``, CDS exons typed and grouped by transcript.

    ``canonical_only=True`` keeps one MANE Select transcript per protein-coding
    gene (the recall reference); ``False`` keeps every CDS-bearing transcript
    (all isoforms -- the precision reference). GENCODE is feature-ordered
    (transcript before its CDS), so a single pass suffices."""
    strand: dict[str, str] = {}
    cds: dict[str, list[tuple[int, int]]] = defaultdict(list)
    with _open(gff3_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[0] != seqid:
                continue
            if f[2] == "transcript":
                a = _attrs(f[8])
                mane = "MANE_Select" in a.get("tag", "") and a.get("gene_type") == "protein_coding"
                if mane or not canonical_only:
                    strand[a["ID"]] = f[6]
            elif f[2] == "CDS":
                parent = _attrs(f[8]).get("Parent", "")
                if parent in strand:
                    cds[parent].append((int(f[3]), int(f[4])))

    # info line first (evaluate reads field 5 as the locus length), then exons
    lines = [f"{seqid}\tgencode\tinfo\t1\t{seqlen}\t.\t.\t.\tinfo"]
    for tid, exons in cds.items():
        exons.sort()
        ordered = exons if strand[tid] == "+" else list(reversed(exons))
        n = len(ordered)
        for i, (s, e) in enumerate(ordered):
            if n == 1:
                typ = "Single"
            else:
                typ = "First" if i == 0 else "Terminal" if i == n - 1 else "Internal"
            lines.append(f"{seqid}\tgencode\t{typ}\t{s}\t{e}\t.\t{strand[tid]}\t.\t{tid}")
    return "\n".join(lines) + "\n"


def gencode_mane_gp(gff3_path: str | Path, seqid: str, seqlen: int) -> str:
    """The canonical (MANE Select) CDS truth -- the recall reference."""
    return gencode_cds_gp(gff3_path, seqid, seqlen, canonical_only=True)


def run_geneid(geneid_bin: str, param: str, fasta: str, flags: tuple[str, ...] = ()) -> str:
    """Run ``geneid -G <flags> -P param fasta`` and keep the typed-CDS-exon
    prediction lines (First/Internal/Terminal/Single)."""
    out = subprocess.run(
        [geneid_bin, "-G", *flags, "-P", param, fasta],
        capture_output=True, text=True, check=True,
    ).stdout
    kept = [ln for ln in out.splitlines()
            if len(f := ln.split("\t")) >= 9 and f[2] in _TYPED]
    return "\n".join(kept) + "\n"


@dataclass
class ModeResult:
    name: str
    accuracy: Accuracy


def run_benchmark(
    geneid_bin: str, param: str, fasta: str, recall_gp: str, precision_gp: str,
    modes: list[tuple[str, tuple[str, ...]]], workdir: str | Path,
) -> list[ModeResult]:
    """Score each ``(name, flags)`` mode: recall vs the canonical truth
    ``recall_gp``, precision vs the all-isoform truth ``precision_gp``."""
    work = Path(workdir)
    work.mkdir(parents=True, exist_ok=True)
    results = []
    for name, flags in modes:
        pred = work / f"pred.{name}.gff"
        pred.write_text(run_geneid(geneid_bin, param, fasta, flags))
        results.append(ModeResult(name, evaluate_asymmetric(str(pred), recall_gp, precision_gp)))
    return results


def format_table(results: list[ModeResult]) -> str:
    """A compact SN/SP table (nucleotide, exon, gene levels)."""
    hdr = f"{'mode':<16}{'nSN':>7}{'nSP':>7}{'eSN':>7}{'eSP':>7}{'eSNSP':>7}{'gSN':>7}{'gSP':>7}"
    rows = [hdr, "-" * len(hdr)]
    for r in results:
        a = r.accuracy
        rows.append(
            f"{r.name:<16}{a.sn:>7.3f}{a.sp:>7.3f}{a.sne:>7.3f}{a.spe:>7.3f}"
            f"{a.snsp:>7.3f}{a.sng:>7.3f}{a.spg:>7.3f}"
        )
    return "\n".join(rows)
