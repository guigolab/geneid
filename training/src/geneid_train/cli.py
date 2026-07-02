"""``geneid-train`` command-line entry point.

Phase 1 ships a working ``param-info`` inspector and stubs for the pipeline
subcommands so the CLI surface and wiring exist before the statistics land.
"""

from __future__ import annotations

import argparse
import sys

from . import __version__
from .core.param import Param

_U12_MARKERS = ("U12_Splice_Score_Threshold", "U12_Branch_point_profile")

_STUBS = {
    "prepare": "build a validated training set from BUSCO or RNA-seq/TransDecoder input (phase 2)",
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
