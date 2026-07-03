#!/usr/bin/env python3
"""Regenerate the bundled pan-taxon U12 profiles from the IAOD intron set.

This is a *build* tool, not part of the runtime: it derives the U12 donor /
acceptor / branch profiles once from the pooled IAOD U12 introns and writes them
into ``src/geneid_train/data/u12_profiles.param`` (shipped in the wheel) plus an
attribution NOTICE. ``geneid-train train --u12`` then loads that bundled file; it
never needs the IAOD source data or a reference param at runtime.

Usage:
    python tools/build_u12_bundle.py \
        --iaod-dir /path/to/dir/of/{genome}_U12.fasta \
        --reference /path/to/human3isoU12.param

The IAOD data (Moyer, Larue, Hershberger, Roy, Padgett, NAR 2020,
doi:10.1093/nar/gkaa464; https://introndb.lerner.ccf.org) is CC-BY; the derived
PWMs are redistributed under that license with attribution.
"""

from __future__ import annotations

import argparse
import glob
from pathlib import Path

from geneid_train.core.param import Param
from geneid_train.param.u12 import build_u12_sections
from geneid_train.prepare.u12 import by_subtype, load_u12_introns

_NOTICE = """\
Bundled U12 (minor-spliceosome) splice profiles
===============================================

The file `u12_profiles.param` contains geneid U12gtag / U12atac donor, acceptor
and branch-point position-weight arrays. Their VALUES are derived from the pooled
U12-type intron set of the Intron Annotation and Orthology Database (IAOD):

    Moyer DC, Larue GE, Hershberger CE, Roy SW, Padgett RA.
    "Comprehensive database and evolutionary dynamics of U12-type introns."
    Nucleic Acids Research 48(13):7066-7078, 2020. doi:10.1093/nar/gkaa464
    https://introndb.lerner.ccf.org

IAOD is distributed under CC-BY; these derived PWMs are redistributed under the
same license with attribution to the authors above.

The profile GEOMETRY (length / offset / cutoff / order / score factors and the
branch-point distance parameters) is inherited from geneid's reference U12
parameter set (human3isoU12.param); the IAOD-derived PWM values are registered
onto that geometry (see geneid_train.param.u12). Positions outside the IAOD
intron-only window are left neutral (0).
"""


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--iaod-dir", required=True, help="dir of {genome}_U12.fasta files")
    ap.add_argument("--reference", required=True, help="reference U12 param for geometry")
    ap.add_argument("--branch-opt-dist", type=int, default=None, help="override branch opt_dist")
    ap.add_argument(
        "--out",
        default=str(Path(__file__).resolve().parents[1] / "src/geneid_train/data"),
        help="output data dir (default: the package data dir)",
    )
    args = ap.parse_args()

    fastas = sorted(glob.glob(f"{args.iaod_dir}/*_U12.fasta"))
    if not fastas:
        ap.error(f"no *_U12.fasta files in {args.iaod_dir}")
    introns = load_u12_introns(fastas)
    groups = by_subtype(introns)
    print(
        f"pooled {len(introns)} U12 introns "
        f"(gtag={len(groups['gtag'])}, atac={len(groups['atac'])}, other={len(groups['other'])}) "
        f"from {len(fastas)} genomes"
    )

    ref = Param.read(args.reference)
    sections = build_u12_sections(introns, ref, branch_opt_dist=args.branch_opt_dist)

    outdir = Path(args.out)
    outdir.mkdir(parents=True, exist_ok=True)
    header = (
        "# Pan-taxon U12 splice profiles derived from the IAOD U12 intron set.\n"
        "# See U12_PROFILES.NOTICE for provenance and license (CC-BY, Roy et al. 2020).\n"
        "# Regenerate with tools/build_u12_bundle.py.\n"
    )
    param_path = outdir / "u12_profiles.param"
    param_path.write_text(header + sections.acceptor_side + sections.donor_side)
    (outdir / "U12_PROFILES.NOTICE").write_text(_NOTICE)
    print(f"wrote {param_path}")
    print(f"wrote {outdir / 'U12_PROFILES.NOTICE'}")


if __name__ == "__main__":
    main()
