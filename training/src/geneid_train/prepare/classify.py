"""Splice-class tally and per-class training recommendations.

Given a validated gene set, count intron donor/acceptor dinucleotides and map them
to the geneid profile classes, then recommend — per rare class — whether to train
a profile de novo, transplant an existing cross-taxa one, or omit it. This is the
report that drives the optional-profile decision (GC donors, U12) described in
DESIGN.md; it is deliberately lightweight and does not attempt full-fidelity
intron classification.

Note: U12 GT-AG introns are indistinguishable from bulk U2 GT-AG by dinucleotide
alone — separating them requires scoring against an existing U12 branch/donor/
acceptor model (bootstrap), which is a separate step. This report flags that
rather than guessing.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from dataclasses import dataclass

from .base import GeneModel

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
class SpliceReport:
    n_models: int
    n_multiexonic: int
    n_introns: int
    donor_counts: dict[str, int]
    acceptor_counts: dict[str, int]
    pair_counts: dict[tuple[str, str], int]
    classes: list[ClassRecommendation]


def _recommend(count: int, min_sites: int) -> str:
    if count == 0:
        return "omit (none found)"
    if count >= min_sites:
        return f"train de novo (compare held-out SN/SP vs transplant); n={count}"
    return f"transplant (too few to train: n={count} < {min_sites})"


def classify_report(
    models: Iterable[GeneModel], genome: Mapping[str, str], min_sites: int = 50
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
    # U12 GT-AG cannot be seen by dinucleotide alone
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
    )
