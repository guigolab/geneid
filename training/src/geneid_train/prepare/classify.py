"""Splice-class tally and per-class training recommendations.

Given a validated gene set, count intron donor/acceptor dinucleotides and map them
to the geneid profile classes, then recommend — per rare class — whether to train
a profile de novo, transplant an existing cross-taxa one, or omit it. This is the
report that drives the optional-profile decision (GC donors, U12) described in
DESIGN.md; it is deliberately lightweight and does not attempt full-fidelity
intron classification.

Note: U12 GT-AG introns are indistinguishable from bulk U2 GT-AG by dinucleotide
alone. To estimate how many are present :func:`detect_u12_gtag` compares each
GT-AG intron's 5' donor under a bundled U12 donor model (IAOD-derived; the U12
5'SS ``GTATCCTT`` is highly distinctive) against a U2 donor model trained from the
genome's own GT-AG introns, and calls U12 only when the U12 score both beats U2
and clears an absolute floor. This is a screen to inform the train-de-novo /
transplant decision, not a per-intron classifier.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from dataclasses import dataclass, field

from .base import GeneModel

_ACGT = frozenset("ACGT")

# (donor, acceptor) -> (human-readable class, geneid profile section it feeds)
_CLASS_MAP = {
    ("GT", "AG"): ("U2 GT-AG", "U2gta_Donor_profile"),
    ("GC", "AG"): ("U2 GC-AG", "U2gcag_Donor_profile"),
    ("AT", "AC"): ("U12 AT-AC", "U12atac_Donor_profile (+acceptor+branch trio)"),
}
# rare classes we make an explicit train/transplant/omit call on
_OPTIONAL = {("GC", "AG"), ("AT", "AC")}


@dataclass
class ClassRecommendation:
    name: str
    profile: str
    count: int
    fraction: float
    recommendation: str


@dataclass
class U12GtagEstimate:
    """Screen for U12-type GT-AG introns by a U12-vs-U2 donor comparison.

    Each GT-AG intron's 5' donor window is scored under a bundled U12 donor PWM
    (IAOD-derived) and a U2 donor PWM trained from *this genome's own* GT-AG
    introns. An intron is a candidate only if BOTH: its U12 score beats its U2
    score by ``margin`` (it looks more U12 than U2), AND the U12 score clears an
    absolute ``floor`` (the U12 motif is actually present) — introns that score
    poorly under both models are cryptic or U2, not U12. ``top_margins`` are the
    largest U12-minus-U2 differences (a ranked shortlist). Still a screen, not a
    calibrated classifier; confirm with a dedicated tool (intronIC / BPP).
    """

    n_scored: int
    n_candidates: int
    margin: float
    floor: float
    top_margins: list[float] = field(default_factory=list)


@dataclass
class SpliceReport:
    n_models: int
    n_multiexonic: int
    n_introns: int
    donor_counts: dict[str, int]
    acceptor_counts: dict[str, int]
    pair_counts: dict[tuple[str, str], int]
    classes: list[ClassRecommendation]
    u12_gtag: U12GtagEstimate | None = None


def detect_u12_gtag(
    models: Iterable[GeneModel],
    genome: Mapping[str, str],
    *,
    margin: float = 0.0,
    floor: float | None = None,
    max_report: int = 8,
) -> U12GtagEstimate:
    """Screen GT-AG introns for likely U12 members by a U12-vs-U2 donor comparison.

    The U12 donor PWM is bundled (IAOD-derived); the U2 donor PWM is trained here
    from the genome's own GT-AG donor windows. An intron is a candidate when its
    U12 donor log-likelihood both beats its U2 log-likelihood by ``margin`` and
    clears an absolute ``floor`` (default: the bundled calibration floor — the low
    percentile of known U12 donor scores). See :class:`U12GtagEstimate`.
    """
    from ..param.u12 import donor_loglik, load_u12_donor_model
    from ..stats.sites import frequency

    u12_freq, length, cal_floor = load_u12_donor_model()
    floor = cal_floor if floor is None else floor

    windows: list[str] = []
    for m in models:
        if m.seqid not in genome:
            continue
        for iseq in m.intron_seqs(genome):
            if len(iseq) < length or iseq[:2] != "GT" or iseq[-2:] != "AG":
                continue
            w = iseq[:length]
            if not set(w) <= _ACGT:
                continue
            windows.append(w)

    u2_freq = frequency(windows)  # the genome's own bulk-GT-AG (mostly U2) donor model
    diffs: list[float] = []
    n_candidates = 0
    for w in windows:
        s12 = donor_loglik(w, u12_freq, length)
        s2 = donor_loglik(w, u2_freq, length)
        diffs.append(s12 - s2)
        if s12 - s2 >= margin and s12 >= floor:
            n_candidates += 1
    diffs.sort(reverse=True)
    return U12GtagEstimate(len(windows), n_candidates, margin, floor, diffs[:max_report])


def _recommend(count: int, min_sites: int) -> str:
    if count == 0:
        return "omit (none found)"
    if count >= min_sites:
        return f"train de novo (compare held-out SN/SP vs transplant); n={count}"
    return f"transplant (too few to train: n={count} < {min_sites})"


def classify_report(
    models: Iterable[GeneModel],
    genome: Mapping[str, str],
    min_sites: int = 50,
    *,
    bootstrap_u12: bool = True,
    u12_floor: float | None = None,
) -> SpliceReport:
    models = list(models)
    donor: dict[str, int] = {}
    acceptor: dict[str, int] = {}
    pair: dict[tuple[str, str], int] = {}
    n_introns = 0
    for m in models:
        if m.seqid not in genome:
            continue
        for iseq in m.intron_seqs(genome):
            if len(iseq) < 4:
                continue
            d, a = iseq[:2], iseq[-2:]
            if "N" in d or "N" in a:
                continue
            n_introns += 1
            donor[d] = donor.get(d, 0) + 1
            acceptor[a] = acceptor.get(a, 0) + 1
            pair[(d, a)] = pair.get((d, a), 0) + 1

    classes: list[ClassRecommendation] = []
    denom = n_introns or 1
    for key, (name, profile) in _CLASS_MAP.items():
        count = pair.get(key, 0)
        if key in _OPTIONAL:
            rec = _recommend(count, min_sites)
        else:
            rec = "bulk class (always trained)"
        classes.append(ClassRecommendation(name, profile, count, count / denom, rec))
    # U12 GT-AG can't be seen by dinucleotide alone — bootstrap-score the branch
    u12_gtag = None
    if bootstrap_u12:
        u12_gtag = detect_u12_gtag(models, genome, floor=u12_floor)
        n = u12_gtag.n_candidates
        screen = "U12-vs-U2 donor screen — confirm with a U12 classifier"
        rec = (
            f"candidates suggest training de novo may be worthwhile ({screen}); n~{n}"
            if n >= min_sites
            else f"transplant the bundled U12 profiles ({screen}); n~{n}"
            if n > 0
            else f"omit or transplant (no U12-like GT-AG donors; {screen})"
        )
        classes.append(
            ClassRecommendation(
                "U12 GT-AG",
                "U12gtag_Donor_profile (+acceptor+branch trio)",
                n,
                n / denom,
                rec,
            )
        )
    else:
        classes.append(
            ClassRecommendation(
                "U12 GT-AG",
                "U12gtag_Donor_profile (+acceptor+branch trio)",
                -1,
                0.0,
                "requires bootstrap scoring of GT-AG introns vs an existing U12 param",
            )
        )
    return SpliceReport(
        n_models=len(models),
        n_multiexonic=sum(m.is_multiexonic for m in models),
        n_introns=n_introns,
        donor_counts=dict(sorted(donor.items(), key=lambda kv: -kv[1])),
        acceptor_counts=dict(sorted(acceptor.items(), key=lambda kv: -kv[1])),
        pair_counts=pair,
        classes=classes,
        u12_gtag=u12_gtag,
    )
