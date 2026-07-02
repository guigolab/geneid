"""GFF3 reader/writer — the tool's one canonical annotation format.

geneid-train ingests GFF3 exclusively (it is the most robust flavor). Other
flavors (GFF2, GTF) are handled by explicit converters in ``core.convert``, never
parsed directly by the pipeline. Coordinates are 1-based inclusive, matching GFF3
and geneid conventions.

CDS features are grouped into transcripts by their ``Parent`` attribute.
"""

from __future__ import annotations

from collections.abc import Iterable, Iterator
from dataclasses import dataclass, field
from pathlib import Path


def parse_attributes(col9: str) -> dict[str, str]:
    """Parse a GFF3 column-9 ``key=value;key=value`` string into a dict.

    Order is preserved (dict insertion order). Percent-decoding is intentionally
    not applied; ids in practice are already safe tokens.
    """
    attrs: dict[str, str] = {}
    for field_ in col9.strip().split(";"):
        field_ = field_.strip()
        if not field_:
            continue
        key, sep, value = field_.partition("=")
        if sep:
            attrs[key.strip()] = value.strip()
    return attrs


def format_attributes(attrs: dict[str, str]) -> str:
    return ";".join(f"{k}={v}" for k, v in attrs.items())


@dataclass
class GffRecord:
    seqid: str
    source: str
    type: str
    start: int  # 1-based inclusive
    end: int  # 1-based inclusive
    score: str  # "." or a float, kept as string
    strand: str  # "+", "-", or "."
    phase: str  # "0", "1", "2", or "."  (GFF3 calls col 8 "phase")
    attributes: dict[str, str] = field(default_factory=dict)

    @property
    def id(self) -> str | None:
        return self.attributes.get("ID")

    @property
    def parents(self) -> list[str]:
        raw = self.attributes.get("Parent")
        return raw.split(",") if raw else []

    def to_line(self) -> str:
        return "\t".join(
            [
                self.seqid,
                self.source,
                self.type,
                str(self.start),
                str(self.end),
                self.score,
                self.strand,
                self.phase,
                format_attributes(self.attributes),
            ]
        )


def iter_gff3(path: str | Path) -> Iterator[GffRecord]:
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 8:
                raise ValueError(f"malformed GFF3 line ({len(f)} cols): {line!r}")
            yield GffRecord(
                seqid=f[0],
                source=f[1],
                type=f[2],
                start=int(f[3]),
                end=int(f[4]),
                score=f[5],
                strand=f[6],
                phase=f[7],
                attributes=parse_attributes(f[8]) if len(f) > 8 else {},
            )


def read_gff3(path: str | Path) -> list[GffRecord]:
    return list(iter_gff3(path))


def write_gff3(records: Iterable[GffRecord], path: str | Path, header: bool = True) -> None:
    with open(path, "w") as fh:
        if header:
            fh.write("##gff-version 3\n")
        for rec in records:
            fh.write(rec.to_line() + "\n")
