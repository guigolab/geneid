"""``geneid-train`` command-line entry point.

Phase 1 ships a working ``param-info`` inspector and stubs for the pipeline
subcommands so the CLI surface and wiring exist before the statistics land.
"""

from __future__ import annotations

import argparse
import sys

from . import __version__
from .core.fasta import read_fasta, write_fasta
from .core.gff import read_gff, write_gff
from .core.param import Param
from .prepare.base import (
    build_models,
    filter_complete,
    filter_min_protein,
    filter_non_overlapping,
)

_U12_MARKERS = ("U12_Splice_Score_Threshold", "U12_Branch_point_profile")

_STUBS = {
    "train": "estimate site + coding models and assemble a .param file (phases 3-4)",
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


def _cmd_prepare(args: argparse.Namespace) -> int:
    genome = read_fasta(args.fastas)
    models = build_models(read_gff(args.gff))
    print(f"input models:      {len(models)}")
    if not models:
        sys.stderr.write("no CDS-grouped gene models found in GFF\n")
        return 1
    print(f"  multi-exonic:    {sum(m.is_multiexonic for m in models)}")

    kept = filter_complete(models, genome)
    print(f"complete CDS:      {len(kept)}")
    kept = filter_min_protein(kept, genome, args.min_aa)
    print(f">= {args.min_aa} aa:{' ' * max(1, 8 - len(str(args.min_aa)))}{len(kept)}")
    kept = filter_non_overlapping(kept)
    print(f"non-overlapping:   {len(kept)}")

    if args.results:
        recs = [r for m in kept for r in _model_records(m)]
        write_gff(recs, f"{args.results}.validated.gff")
        write_fasta({m.gene_id: m.cds(genome) for m in kept}, f"{args.results}.validated.cds.fa")
        write_fasta(
            {m.gene_id: m.protein(genome).rstrip("*") for m in kept},
            f"{args.results}.validated.prot.fa",
        )
        print(f"wrote:             {args.results}.validated.{{gff,cds.fa,prot.fa}}")
    return 0


def _model_records(model):
    from .core.gff import GffRecord

    return [
        GffRecord(model.seqid, "geneid_train", "CDS", e.start, e.end, ".", model.strand,
                  e.frame, model.gene_id)
        for e in model.exons
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
    p_prep.add_argument("--gff", required=True, help="GFF2 of CDS features grouped by gene id")
    p_prep.add_argument("--fastas", required=True, help="genomic multi-FASTA")
    p_prep.add_argument("--min-aa", type=int, default=100, help="minimum protein length (aa)")
    p_prep.add_argument("--flank", type=int, default=1000, help="flank nt for locus extraction")
    p_prep.add_argument(
        "--results", help="output path prefix; writes validated.{gff,cds.fa,prot.fa}"
    )
    p_prep.set_defaults(func=_cmd_prepare)

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
