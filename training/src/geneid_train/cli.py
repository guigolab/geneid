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

_STUBS: dict[str, str] = {}


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
    rep = classify_report(
        models, genome, min_sites=args.min_sites,
        bootstrap_u12=not args.no_u12_bootstrap, u12_floor=args.u12_floor,
    )
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
    if rep.u12_gtag is not None:
        e = rep.u12_gtag
        top = ", ".join(f"{s:.1f}" for s in e.top_margins)
        print(
            f"\nU12 GT-AG screen (U12-vs-U2 donor, not a calibrated count): "
            f"{e.n_candidates}/{e.n_scored} GT-AG introns score more U12 than U2 "
            f"(margin >= {e.margin:g}, U12 floor >= {e.floor:.2f})"
        )
        print(f"  top U12-minus-U2 donor log-likelihood: {top}")
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


def _cmd_jackknife(args: argparse.Namespace) -> int:
    from .evaluate import read_annotation_gff
    from .jackknife import jackknife

    genome = read_fasta(args.fastas)
    records = read_gff3(args.gff)
    models = collapse_isoforms(build_models(records), records)
    models = filter_non_overlapping(
        filter_min_protein(filter_complete(models, genome), genome, args.min_aa)
    )
    if not models:
        sys.stderr.write("no models survived filtering\n")
        return 1
    locus_fasta = read_fasta(args.eval_fastas)
    annotations = read_annotation_gff(args.eval_gff)
    weights = (args.ewf, args.owf) if args.ewf is not None and args.owf is not None else None
    acc = jackknife(
        models, genome, args.species, locus_fasta, annotations,
        geneid_bin=args.geneid, folds=args.folds, weights=weights,
    )
    print(f"{args.folds}-fold cross-validation ({len(models)} models):")
    print(f"nucleotide  SN={acc.sn:.3f} SP={acc.sp:.3f} CC={acc.cc:.3f}")
    print(f"exon        SNe={acc.sne:.3f} SPe={acc.spe:.3f} SNSP={acc.snsp:.3f}")
    print(f"gene        SNg={acc.sng:.3f} SPg={acc.spg:.3f} SNSPg={acc.snspg:.3f}")
    return 0


def _cmd_evaluate(args: argparse.Namespace) -> int:
    from .evaluate import evaluate_files

    a = evaluate_files(args.predictions, args.annotations)
    print("level        SN     SP     combined")
    print(f"nucleotide  {a.sn:6.3f} {a.sp:6.3f}   CC={a.cc:.3f}")
    print(f"exon        {a.sne:6.3f} {a.spe:6.3f}   SNSP={a.snsp:.3f}")
    print(f"gene        {a.sng:6.3f} {a.spg:6.3f}   SNSPg={a.snspg:.3f}")
    print(f"raME={a.ra_me:.3f} raWE={a.ra_we:.3f}")
    return 0


def _cmd_optimize(args: argparse.Namespace) -> int:
    from .optimize import coordinate_descent, global_optimize, optimize, uniform_point

    base = open(args.param).read()
    types = "First/Internal/Terminal/Single"

    if args.strategy == "global":
        from .optimize import SearchSpace

        # Latin-hypercube global sampling + compass-search refinement (per type)
        keywords = Param.from_text(base).keywords()
        branch_profiles: tuple[str, ...] = ()
        if getattr(args, "tune_branch", False):
            branch = "U12_Branch_point_profile"
            if branch in keywords:
                branch_profiles = (branch,)
            else:
                print(f"note: --tune-branch ignored (no {branch} in {args.param})")
        tune_il = getattr(args, "tune_intron_length", False)
        if tune_il and "Intron_length_score_weight" not in keywords:
            print(f"note: --tune-intron-length ignored (no Intron_length_model in {args.param})")
            tune_il = False
        space = SearchSpace(branch_profiles=branch_profiles, tune_intron_length=tune_il)
        opt_text, res = global_optimize(
            base, args.eval_fastas, args.eval_gff,
            geneid_bin=args.geneid, space=space, n_samples=args.samples,
            workers=args.workers, seed=args.seed,
        )
        print(f"global search: {res.n_evaluations} evals, {len(res.history)} improving moves "
              f"-> SNSP={res.accuracy.snsp:.4f}")
        print(f"  eWF [{types}] = {tuple(round(x, 3) for x in res.point.ewf)}")
        print(f"  oWF [{types}] = {tuple(round(x, 3) for x in res.point.owf)}")
        for name, kn in res.branch:
            print(f"  {name}: acc_context={kn.acc_context} min_dist={kn.min_dist} "
                  f"opt_dist={kn.opt_dist} pen_scale={kn.pen_scale:g}")
        if res.intron_length_weight is not None:
            print(f"  Intron_length_score_weight (lambda) = {res.intron_length_weight:g}")
    else:
        opt_text, results = optimize(
            base, args.eval_fastas, args.eval_gff,
            geneid_bin=args.geneid, workers=args.workers,
        )
        best = results[0]
        print(f"uniform grid: {len(results)} points, best eWF={best.ewf:g} oWF={best.owf:g} "
              f"-> SNSP={best.accuracy.snsp:.4f}")
        if args.strategy == "per-type":
            opt_text, cd, history = coordinate_descent(
                base, args.eval_fastas, args.eval_gff, geneid_bin=args.geneid,
                init=uniform_point(best.ewf, best.owf), workers=args.workers,
            )
            print(f"per-type refine ({len(history)} improving steps) -> "
                  f"SNSP={cd.accuracy.snsp:.4f}")
            print(f"  eWF [{types}] = {cd.point.ewf}")
            print(f"  oWF [{types}] = {cd.point.owf}")

    with open(args.output, "w") as fh:
        fh.write(opt_text)
    print(f"wrote optimized parameter file: {args.output}")
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
        param_text = train(
            models, genome, args.species, seed=args.seed, u12=args.u12,
            u12_splice_thresh=args.u12_splice_thresh, u12_exon_thresh=args.u12_exon_thresh,
            u2_branch=args.u2_branch, branch_weight=args.branch_weight,
        )
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
    p_cls.add_argument(
        "--u12-floor", type=float, default=None,
        help="absolute U12 donor log-likelihood floor for the U12-vs-U2 GT-AG screen "
             "(default: the bundled calibration floor)",
    )
    p_cls.add_argument(
        "--no-u12-bootstrap", action="store_true",
        help="skip the U12 GT-AG branch-model bootstrap scoring",
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
    p_train.add_argument(
        "--u12",
        action="store_true",
        help="include bundled U12 (minor-spliceosome) profiles so geneid -U predicts U12 introns",
    )
    p_train.add_argument(
        "--u12-splice-thresh", type=float, default=9.0,
        help="U12_Splice_Score_Threshold: min combined U12 donor+acceptor score to accept a "
             "U12 join (default 9, conservative; reference U12 params; lower over-calls U12)",
    )
    p_train.add_argument(
        "--u12-exon-thresh", type=float, default=8.0,
        help="U12_Exon_Score_Threshold: min combined exon score for a U12 join (default 8)",
    )
    p_train.add_argument(
        "--u2-branch",
        action="store_true",
        help="discover a U2 Branch_point_profile from the genome's introns (EM) and include it",
    )
    p_train.add_argument(
        "--branch-weight", type=float, default=0.0,
        help="Branch_point_score_weight: contribution of the branch score to the acceptor "
             "score (default 0 = scored and reported via bp_score/bp_pos but not counted)",
    )
    p_train.set_defaults(func=_cmd_train)

    p_opt = sub.add_parser(
        "optimize", help="grid-search exon/site weights against a held-out set (maximise SNSP)"
    )
    p_opt.add_argument("--param", required=True, help="base trained .param file")
    p_opt.add_argument("--eval-fastas", required=True, help="held-out locus FASTA (gp format)")
    p_opt.add_argument("--eval-gff", required=True, help="held-out annotation GFF (gp convention)")
    p_opt.add_argument("--output", required=True, help="output optimized .param path")
    p_opt.add_argument("--geneid", default="geneid", help="path to the geneid binary")
    p_opt.add_argument("--workers", type=int, default=4, help="parallel geneid runs")
    p_opt.add_argument(
        "--strategy", choices=["uniform", "per-type", "global"], default="uniform",
        help="uniform grid (default); per-type coordinate descent; or global "
             "Latin-hypercube sampling + compass refinement over the 8 per-type weights",
    )
    p_opt.add_argument("--samples", type=int, default=32, help="global: LHS sample count")
    p_opt.add_argument("--seed", type=int, default=0, help="global: LHS RNG seed")
    p_opt.add_argument(
        "--tune-branch", action="store_true",
        help="global: also tune the U12 branch-distance knobs (acc_context/min_dist/"
             "opt_dist/pen_scale) if the param has a U12_Branch_point_profile",
    )
    p_opt.add_argument(
        "--tune-intron-length", action="store_true",
        help="global: also tune the soft intron-length penalty weight "
             "(Intron_length_score_weight) if the param has an Intron_length_model",
    )
    p_opt.set_defaults(func=_cmd_optimize)

    p_jk = sub.add_parser(
        "jackknife", help="leave-group-out cross-validation of a training set"
    )
    p_jk.add_argument("--gff", required=True, help="training GFF3 of CDS features")
    p_jk.add_argument("--fastas", required=True, help="genomic multi-FASTA (for training)")
    p_jk.add_argument("--species", required=True, help="species name")
    p_jk.add_argument("--eval-fastas", required=True, help="per-locus prediction FASTA (gp format)")
    p_jk.add_argument("--eval-gff", required=True, help="per-locus annotation GFF (gp convention)")
    p_jk.add_argument("--geneid", default="geneid", help="path to the geneid binary")
    p_jk.add_argument("--folds", type=int, default=10, help="number of cross-validation folds")
    p_jk.add_argument("--min-aa", type=int, default=100, help="minimum protein length (aa)")
    p_jk.add_argument("--ewf", type=float, help="optimized exon weight to apply per fold")
    p_jk.add_argument("--owf", type=float, help="optimized exon factor to apply per fold")
    p_jk.set_defaults(func=_cmd_jackknife)

    p_eval = sub.add_parser(
        "evaluate", help="score a prediction GFF against an annotation GFF (SN/SP)"
    )
    p_eval.add_argument("predictions", help="geneid prediction GFF (typed CDS exons)")
    p_eval.add_argument("annotations", help="annotation GFF in gp convention (per-locus info line)")
    p_eval.set_defaults(func=_cmd_evaluate)

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
