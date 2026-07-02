"""Sequence utilities: reverse-complement, translation, CDS completeness.

Genetic-code awareness matters here: some geneid target species use non-standard
codon tables (e.g. ciliates translate TAA/TAG as glutamine), which changes what
counts as a stop and therefore what counts as a *complete* CDS. ``GeneticCode``
captures the stop set so completeness checks stay correct per species.
"""

from __future__ import annotations

from dataclasses import dataclass

_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")

# NCBI translation table 1 (standard)
_STANDARD_CODONS = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}


def revcomp(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


@dataclass(frozen=True)
class GeneticCode:
    """A codon table, specialized by which codons are stops and starts."""

    codons: dict[str, str]
    starts: frozenset[str] = frozenset({"ATG"})

    @property
    def stops(self) -> frozenset[str]:
        return frozenset(c for c, aa in self.codons.items() if aa == "*")

    def translate(self, seq: str) -> str:
        s = seq.upper()
        out = []
        for i in range(0, len(s) - len(s) % 3, 3):
            out.append(self.codons.get(s[i : i + 3], "X"))
        return "".join(out)

    def is_stop(self, codon: str) -> bool:
        return self.codons.get(codon.upper()) == "*"

    def is_start(self, codon: str) -> bool:
        return codon.upper() in self.starts


STANDARD = GeneticCode(codons=dict(_STANDARD_CODONS))


def is_complete_cds(cds: str, code: GeneticCode = STANDARD) -> bool:
    """A complete CDS: length a multiple of 3, canonical start, terminal stop,
    and no in-frame internal stop."""
    s = cds.upper()
    if len(s) < 6 or len(s) % 3 != 0:
        return False
    if not code.is_start(s[:3]):
        return False
    if not code.is_stop(s[-3:]):
        return False
    internal = code.translate(s[:-3])
    return "*" not in internal


def has_internal_stop(cds: str, code: GeneticCode = STANDARD) -> bool:
    """In-frame stop before the final codon (ignores a terminal stop)."""
    s = cds.upper()
    body = s[:-3] if len(s) >= 3 and code.is_stop(s[-3:]) else s
    return "*" in code.translate(body)
