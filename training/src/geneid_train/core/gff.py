"""Minimal GFF2 reader/writer using geneid's coordinate conventions.

geneid consumes GFF2: nine tab-separated columns, 1-based inclusive coordinates,
with the ninth column (``group``) naming the gene a feature belongs to. Training
input is typically CDS features grouped by gene id. We keep the group column as a
raw string and expose the primary id (its first token, quotes stripped).
"""

from __future__ import annotations

from collections.abc import Iterable, Iterator
from dataclasses import dataclass
from pathlib import Path


@dataclass
class GffRecord:
    seqid: str
    source: str
    type: str
    start: int  # 1-based inclusive
    end: int  # 1-based inclusive
    score: str  # kept as string; often "." or a float
    strand: str  # "+", "-", or "."
    frame: str  # "0", "1", "2", or "."
    group: str

    @property
    def gene_id(self) -> str:
        tok = self.group.split()[0] if self.group.split() else self.group
        return tok.strip().strip('"').strip("'")

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
                self.frame,
                self.group,
            ]
        )


def iter_gff(path: str | Path) -> Iterator[GffRecord]:
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 8:
                raise ValueError(f"malformed GFF line ({len(f)} cols): {line!r}")
            group = f[8] if len(f) > 8 else ""
            yield GffRecord(
                seqid=f[0],
                source=f[1],
                type=f[2],
                start=int(f[3]),
                end=int(f[4]),
                score=f[5],
                strand=f[6],
                frame=f[7],
                group=group,
            )


def read_gff(path: str | Path) -> list[GffRecord]:
    return list(iter_gff(path))


def write_gff(records: Iterable[GffRecord], path: str | Path) -> None:
    with open(path, "w") as fh:
        for rec in records:
            fh.write(rec.to_line() + "\n")
