"""Emit the U12 (minor-spliceosome) profile sections of a geneid ``.param``.

geneid switches on U12-type splice prediction only when a *complete trio* of
profiles is present for a subtype — a shared ``U12_Branch_point_profile`` plus a
donor and an acceptor profile for that subtype (see ``readparam.c``
``ReadProfileSpliceSites``: the ``u12bp && u12gtagAcc && u12gtagDon`` gate). We
support the two dominant subtypes, GT-AG and AT-AC.

Provenance / design
-------------------
A single genome almost never has enough U12 introns to train these profiles, so
the PWM *values* are derived from the pooled IAOD U12 intron set (Moyer/Larue/
Roy/Padgett 2020; CC-BY — see ``data/U12_PROFILES.NOTICE``), which is intron-only.

The *geometry* (each profile's length / offset / cutoff / order / score factors,
and the branch-point distance knobs) is inherited verbatim from a proven,
geneid-tested reference U12 param (isochore 1 of ``human3isoU12.param``). geneid's
donor/acceptor ``offset`` encodes an exon-context alignment that intron-only data
cannot reproduce from scratch, and getting it wrong shifts predicted splice
coordinates. So instead of guessing the geometry we keep the reference's and
*register* the IAOD PWM onto it: we find the position shift that best aligns the
IAOD log-odds matrix to the reference one (both are log-odds of the same U12
motif, so their column shapes correlate), then overlay the IAOD values on the
positions that shift covers. Positions the intron-only window can't inform (exon
flanks) are left neutral (0), never scored against.

This is the Stage-1 bundle: pan-taxon IAOD signal on a correct geneid geometry.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

from ..core.param import Param, Profile
from ..prepare.u12 import (
    U12Intron,
    acceptor_windows,
    by_subtype,
    donor_windows,
    locate_branches,
    markov_background,
)
from ..stats.sites import Matrix, log_ratio, position_matrix, submatrix

# The reference U12 profiles whose geometry we inherit. Acceptor-side profiles
# are emitted before ``Acceptor_profile``; donor-side before ``Donor_profile``.
ACCEPTOR_SIDE = ("U12_Branch_point_profile", "U12gtag_Acceptor_profile", "U12atac_Acceptor_profile")
DONOR_SIDE = ("U12gtag_Donor_profile", "U12atac_Donor_profile")


@dataclass
class U12Sections:
    """Rendered geneid section text for the two insertion points."""

    acceptor_side: str  # branch + both acceptors, before Acceptor_profile
    donor_side: str  # both donors, before Donor_profile


def _full_matrix(seqs: Sequence[str], order: int, background: Matrix) -> Matrix:
    """Unwindowed order-k log-odds PWM of ``seqs`` against ``background`` — the raw
    material for registration (positions 1..L-order, all oligos)."""
    site = position_matrix(seqs, order=order)
    if not site:
        return {}
    maxpos = max(p for p, _ in site)
    return submatrix(log_ratio(site, background), 1, maxpos)


def register_shift(iaod: Matrix, ref: Matrix) -> int:
    """Find the integer shift ``s`` mapping IAOD position ``j`` to reference
    position ``j + s`` that best aligns the two PWMs.

    Both matrices are log-odds of the same motif, so the correct registration
    maximises their normalised dot product over shared ``(pos, oligo)`` cells.
    We scan every shift that keeps at least three positions overlapping.
    """
    ipos = sorted({p for p, _ in iaod})
    rpos = sorted({p for p, _ in ref})
    if not ipos or not rpos:
        return 0
    best_s, best_score = 0, float("-inf")
    for s in range(rpos[0] - ipos[-1], rpos[-1] - ipos[0] + 1):
        num = 0.0
        overlap = 0
        for (p, oligo), v in iaod.items():
            rv = ref.get((p + s, oligo))
            if rv is not None:
                num += v * rv
                overlap += 1
        # need a few overlapping positions and require positive correlation
        n_overlap_pos = len({p for p, _ in iaod if (p + s) in rpos})
        if n_overlap_pos >= 3 and num > best_score:
            best_s, best_score = s, num
    return best_s


def overlay_profile(ref: Profile, iaod: Matrix, shift: int) -> list[str]:
    """Produce geneid data lines for a profile with the reference *geometry* and
    IAOD-derived *values*. Reference position ``pos`` takes the IAOD value at
    ``pos - shift`` when available, otherwise a neutral ``0``. Every oligo of the
    profile's order is emitted at every position (a complete matrix)."""
    dim, order = ref.length, ref.order
    oligos = _oligos(order)
    header = _fmt_header(ref.header)
    lines = [header, "# Transition probabilities at every position"]
    for pos in range(1, dim + 1):
        for oligo in oligos:
            v = iaod.get((pos - shift, oligo), 0.0)
            lines.append(f"{pos} {oligo} {_fmt(v)}")
    return lines


def _oligos(order: int) -> list[str]:
    from itertools import product

    return ["".join(p) for p in product("ACGT", repeat=order + 1)]


def _fmt(v: float) -> str:
    if v == int(v):
        return str(int(v))
    return f"{v:g}"


def _fmt_header(header: list[float]) -> str:
    return " ".join(_fmt(h) for h in header)


def _section_text(name: str, lines: list[str]) -> str:
    return name + "\n" + "\n".join(lines) + "\n"


def build_u12_sections(
    introns: Sequence[U12Intron],
    reference: Param,
    *,
    branch_opt_dist: int | None = None,
) -> U12Sections:
    """Train the U12 profile trio from ``introns`` and render the geneid sections,
    inheriting geometry from ``reference`` (a parsed U12 param).

    ``branch_opt_dist`` overrides the reference branch-point ``opt_dist`` knob
    (the IAOD data refines it upward of the human default); the other branch
    distance knobs are inherited and are tuned later by the optimiser.
    """
    groups = by_subtype(introns)
    acc_lines: list[str] = []
    don_lines: list[str] = []

    # --- shared branch point: locate it in every intron (seeded by the reference
    #     branch PWM), then retrain against an intronic background ---
    branch_ref = reference.profile("U12_Branch_point_profile")
    branch_pwm = {(p, o): v for p, o, v in branch_ref.rows}
    order = branch_ref.order
    all_seqs = [it.seq for it in introns]
    hits = locate_branches(
        introns, branch_pwm, order=order, offset=int(branch_ref.header[1]),
        acc_context=_hdr(branch_ref, 6, 50), min_dist=_hdr(branch_ref, 7, 7),
    )
    branch_windows = [h.window for h in hits]
    branch_bg = markov_background(all_seqs, order)
    branch_iaod = _full_matrix(branch_windows, order, branch_bg)
    branch_shift = register_shift(branch_iaod, branch_pwm_matrix(branch_ref))
    branch_header = list(branch_ref.header)
    if branch_opt_dist is not None and len(branch_header) >= 9:
        branch_header[8] = float(branch_opt_dist)
    branch_prof = Profile("U12_Branch_point_profile", branch_header, branch_ref.rows)
    acc_lines.append(
        _section_text(
            "U12_Branch_point_profile",
            overlay_profile(branch_prof, branch_iaod, branch_shift),
        )
    )

    # --- donor + acceptor per subtype ---
    for subtype, gname in (("gtag", "U12gtag"), ("atac", "U12atac")):
        subset = groups[subtype]
        # donor: first intronic bases (the U12 5' splice site)
        dref = reference.profile(f"{gname}_Donor_profile")
        d_iaod = _train_side(donor_windows(subset, dref.length + 6), dref.order)
        d_shift = register_shift(d_iaod, _profile_matrix(dref))
        don_lines.append(
            _section_text(f"{gname}_Donor_profile", overlay_profile(dref, d_iaod, d_shift))
        )
        # acceptor: last intronic bases (PPT + acceptor)
        aref = reference.profile(f"{gname}_Acceptor_profile")
        a_iaod = _train_side(acceptor_windows(subset, aref.length + 6), aref.order)
        a_shift = register_shift(a_iaod, _profile_matrix(aref))
        acc_lines.append(
            _section_text(f"{gname}_Acceptor_profile", overlay_profile(aref, a_iaod, a_shift))
        )

    return U12Sections("".join(acc_lines), "".join(don_lines))


def partition_sections(text: str) -> U12Sections:
    """Split a concatenated U12 profile fragment into the two insertion groups by
    profile name: the branch point and both acceptor profiles go before
    ``Acceptor_profile``; both donor profiles go before ``Donor_profile``."""
    p = Param.from_text(text)
    acc: list[str] = []
    don: list[str] = []
    for kw in p.keywords():
        block = p._find(kw)
        raw = "".join(block.raw)
        (don if kw in DONOR_SIDE else acc).append(raw)
    return U12Sections("".join(acc), "".join(don))


def load_bundled_u12() -> U12Sections:
    """Load the pan-taxon U12 profiles bundled with the package (derived from the
    IAOD U12 intron set — see ``data/U12_PROFILES.NOTICE``)."""
    from importlib.resources import files

    text = files("geneid_train").joinpath("data/u12_profiles.param").read_text()
    return partition_sections(text)


def _train_side(windows: Sequence[str], order: int) -> Matrix:
    bg = markov_background(windows, order)
    return _full_matrix(windows, order, bg)


def _profile_matrix(prof: Profile) -> Matrix:
    return {(p, o): v for p, o, v in prof.rows}


# kept as a named helper for readability at the two call sites
branch_pwm_matrix = _profile_matrix


def _hdr(prof: Profile, idx: int, default: int) -> int:
    return int(prof.header[idx]) if len(prof.header) > idx else default
