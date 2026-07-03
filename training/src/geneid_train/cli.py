"""``geneid-train`` command-line entry point.

Phase 1 ships a working ``param-info`` inspector and stubs for the pipeline
subcommands so the CLI surface and wiring exist before the statistics land.
"""

from __future__ import annotations

import argparse
import sys

from . import __version__
from .core.convert import gff2_to_gff3, gtf_to_gff3
from .core.fasta import read_fasta, read_fasta_subset, write_fasta
from .core.gff import GffRecord, read_gff3, write_gff3
from .core.param import Param
from .prepare.base import (
    build_models,
    collapse_isoforms,
    filter_complete,
    filter_min_protein,
    filter_non_overlapping,
)
from .prepare.classify import classify_report

_U12_MARKERS = ("U12_Splice_Score_Threshold", "U12_Branch_point_profile")

_STUBS = {
    "evaluate": "score a .param against held-out gene models, U2/U12-aware (phase 5)",
    "jackknife": "leave-group-out cross-validation of a training set (phase 7)",
}


def _cmd_param_info(args: argparse.Namespace) -> int:
    param = Param.read(args.path)
    profiles = param.profile_names()
    is_u12 = any(param.has(m) for m in _U12_MARKERS)
    print(f"file:        {args.path}")
    print(f"isochores:   {param.num_isochores}")
    print(f"U12-aware:   {'yes' if is_u12 else 'no'}")
    print(f"NO_SCORE:    {param.scalar('NO_SCORE')}")
    if param.has("Markov_order"):
        print(f"Markov order:{param.scalar('Markov_order')}")
    print(f"profiles ({len(profiles)}):")
    for name in profiles:
        prof = param.profile(name)
        print(f"  {name:<28} length={prof.length} order={prof.order} rows={len(prof.rows)}")
    return 0


def _cmd_convert(args: argparse.Namespace) -> int:
    lines = open(args.input).readlines()
    if args.from_ == "gff2":
        records = gff2_to_gff3(lines)
    elif args.from_ == "gtf":
        records = gtf_to_gff3(lines)
    else:  # pragma: no cover - argparse restricts choices
        raise ValueError(args.from_)
    write_gff3(records, args.output)
    print(f"wrote {len(records)} GFF3 records to {args.output}")
    return 0


def _cmd_classify(args: argparse.Namespace) -> int:
    records = read_gff3(args.gff)
    models = collapse_isoforms(build_models(records), records)
    if not models:
        sys.stderr.write("no CDS-grouped gene models found in GFF3\n")
        return 1
    genome = read_fasta_subset(args.fastas, {m.seqid for m in models})
    rep = classify_report(models, genome, min_sites=args.min_sites)
    print(f"models={rep.n_models} multi-exonic={rep.n_multiexonic} introns={rep.n_introns}")
    top_d = list(rep.donor_counts.items())[:5]
    top_a = list(rep.acceptor_counts.items())[:5]
    print("top donor dinucs:    " + ", ".join(f"{d}={n}" for d, n in top_d))
    print("top acceptor dinucs: " + ", ".join(f"{a}={n}" for a, n in top_a))
    print("\nsplice classes:")
    for c in rep.classes:
        frac = f"{100 * c.fraction:.2f}%" if c.count >= 0 else "  -  "
        cnt = str(c.count) if c.count >= 0 else "?"
        print(f"  {c.name:12} {cnt:>6} {frac:>7}  -> {c.recommendation}")
        print(f"               profile: {c.profile}")
    return 0


def _cmd_prepare(args: argparse.Namespace) -> int:
    genome = read_fasta(args.fastas)
    records = read_gff3(args.gff)
    models = build_models(records)
    if not models:
        sys.stderr.write("no CDS-grouped gene models found in GFF3\n")
        return 1
    print(f"transcripts:       {len(models)}")
    if not args.no_collapse:
        models = collapse_isoforms(models, records)
        print(f"genes (collapsed): {len(models)}")
    print(f"  multi-exonic:    {sum(m.is_multiexonic for m in models)}")

    kept = filter_complete(models, genome)
    print(f"complete CDS:      {len(kept)}")
    kept = filter_min_protein(kept, genome, args.min_aa)
    print(f">= {args.min_aa} aa:{' ' * max(1, 8 - len(str(args.min_aa)))}{len(kept)}")
    kept = filter_non_overlapping(kept)
    print(f"non-overlapping:   {len(kept)}")

    if args.results:
        recs = [r for m in kept for r in _model_records(m)]
        write_gff3(recs, f"{args.results}.validated.gff3")
        write_fasta({m.gene_id: m.cds(genome) for m in kept}, f"{args.results}.validated.cds.fa")
        write_fasta(
            {m.gene_id: m.protein(genome).rstrip("*") for m in kept},
            f"{args.results}.validated.prot.fa",
        )
        print(f"wrote:             {args.results}.validated.{{gff3,cds.fa,prot.fa}}")
    return 0


def _cmd_train(args: argparse.Namespace) -> int:
    from .train import train

    genome = read_fasta(args.fastas)
    records = read_gff3(args.gff)
    models = collapse_isoforms(build_models(records), records)
    if not models:
        sys.stderr.write("no CDS-grouped gene models found in GFF3\n")
        return 1
    models = filter_non_overlapping(
        filter_min_protein(filter_complete(models, genome), genome, args.min_aa)
    )
    sys.stderr.write(f"training on {len(models)} complete, non-overlapping gene models\n")
    if not models:
        sys.stderr.write("no models survived filtering\n")
        return 1
    try:
        param_text = train(models, genome, args.species, seed=args.seed)
    except NotImplementedError as exc:
        sys.stderr.write(f"geneid-train train: {exc}\n")
        return 2
    with open(args.output, "w") as fh:
        fh.write(param_text)
    print(f"wrote parameter file: {args.output}")
    return 0


def _model_records(model) -> list[GffRecord]:
    return [
        GffRecord(
            model.seqid, "geneid_train", "CDS", e.start, e.end, ".", model.strand, e.phase,
            {"ID": f"{model.gene_id}.cds{i}", "Parent": model.gene_id},
        )
        for i, e in enumerate(model.exons, start=1)
    ]


def _make_stub(name: str, desc: str):
    def _run(_args: argparse.Namespace) -> int:
        sys.stderr.write(f"geneid-train {name}: not implemented yet — {desc}\n")
        return 2

    return _run


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="geneid-train", description=__doc__)
    parser.add_argument("--version", action="version", version=f"geneid-train {__version__}")
    sub = parser.add_subparsers(dest="command", required=True)

    p_info = sub.add_parser("param-info", help="summarize a geneid .param file")
    p_info.add_argument("path", help="path to a .param file")
    p_info.set_defaults(func=_cmd_param_info)

    p_prep = sub.add_parser(
        "prepare", help="build a validated training set from a GFF + genomic FASTA"
    )
    p_prep.add_argument("--gff", required=True, help="GFF3 of CDS features (Parent = transcript)")
    p_prep.add_argument("--fastas", required=True, help="genomic multi-FASTA")
    p_prep.add_argument("--min-aa", type=int, default=100, help="minimum protein length (aa)")
    p_prep.add_argument("--flank", type=int, default=1000, help="flank nt for locus extraction")
    p_prep.add_argument(
        "--no-collapse", action="store_true",
        help="keep every transcript instead of one representative (longest) per gene",
    )
    p_prep.add_argument(
        "--results", help="output path prefix; writes validated.{gff,cds.fa,prot.fa}"
    )
    p_prep.set_defaults(func=_cmd_prepare)

    p_cls = sub.add_parser(
        "classify", help="tally intron splice classes (GT-AG/GC-AG/AT-AC) and recommend profiles"
    )
    p_cls.add_argument("--gff", required=True, help="GFF3 of CDS features (Parent = transcript)")
    p_cls.add_argument("--fastas", required=True, help="genomic multi-FASTA")
    p_cls.add_argument(
        "--min-sites", type=int, default=50, help="min sites to train a rare-class profile de novo"
    )
    p_cls.set_defaults(func=_cmd_classify)

    p_train = sub.add_parser(
        "train", help="estimate site + coding models and assemble a geneid .param file"
    )
    p_train.add_argument("--gff", required=True, help="GFF3 of CDS features (Parent = transcript)")
    p_train.add_argument("--fastas", required=True, help="genomic multi-FASTA")
    p_train.add_argument("--species", required=True, help="species name for the param header")
    p_train.add_argument("--output", required=True, help="output .param path")
    p_train.add_argument("--min-aa", type=int, default=100, help="minimum protein length (aa)")
    p_train.add_argument(
        "--seed", type=int, default=0, help="RNG seed for background sampling (reproducibility)"
    )
    p_train.set_defaults(func=_cmd_train)

    p_conv = sub.add_parser("convert", help="convert GFF2 or GTF annotation to canonical GFF3")
    p_conv.add_argument("--from", dest="from_", required=True, choices=["gff2", "gtf"])
    p_conv.add_argument("input", help="input annotation file (GFF2 or GTF)")
    p_conv.add_argument("output", help="output GFF3 path")
    p_conv.set_defaults(func=_cmd_convert)

    for name, desc in _STUBS.items():
        sp = sub.add_parser(name, help=desc)
        sp.set_defaults(func=_make_stub(name, desc))

    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
