"""FASTA reading/writing and the two-column ``tbl`` format.

The legacy pipeline shuttled sequences through a ``tbl`` format (``id<TAB>seq``,
one record per line) via the ``FastaToTbl`` / ``TblToFasta`` shell tools. We keep
that format for interoperability but treat plain dicts as the in-memory currency.
"""

from __future__ import annotations

import gzip
import io
from collections.abc import Iterator, Mapping
from pathlib import Path


def open_text(path: str | Path) -> io.TextIOBase:
    """Open a text file, transparently decompressing gzip (detected by magic
    bytes, so a mislabeled or extensionless ``.gz`` still works)."""
    with open(path, "rb") as fh:
        magic = fh.read(2)
    if magic == b"\x1f\x8b":
        return gzip.open(path, "rt")
    return open(path)


def iter_fasta(path: str | Path) -> Iterator[tuple[str, str]]:
    """Yield ``(id, sequence)`` pairs. The id is the first whitespace-delimited
    token of the header (matching how geneid's tools key sequences).

    Gzip-compressed FASTA (``.fa.gz``) is read transparently."""
    name: str | None = None
    chunks: list[str] = []
    with open_text(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks)
                name = line[1:].split()[0] if line[1:].strip() else ""
                chunks = []
            else:
                chunks.append(line.strip())
    if name is not None:
        yield name, "".join(chunks)


def read_fasta(path: str | Path) -> dict[str, str]:
    out: dict[str, str] = {}
    for name, seq in iter_fasta(path):
        if name in out:
            raise ValueError(f"duplicate sequence id {name!r} in {path}")
        out[name] = seq
    return out


def read_fasta_subset(path: str | Path, ids: set[str], upper: bool = True) -> dict[str, str]:
    """Read only the sequences whose id is in ``ids`` (memory-friendly for large
    genomes). Stops early once every requested id is found."""
    out: dict[str, str] = {}
    for name, seq in iter_fasta(path):
        if name in ids:
            out[name] = seq.upper() if upper else seq
            if len(out) == len(ids):
                break
    return out


def write_fasta(records: Mapping[str, str], path: str | Path, width: int = 60) -> None:
    with open(path, "w") as fh:
        for name, seq in records.items():
            fh.write(f">{name}\n")
            if width and width > 0:
                for i in range(0, len(seq), width):
                    fh.write(seq[i : i + width] + "\n")
            else:
                fh.write(seq + "\n")


def read_tbl(path: str | Path) -> dict[str, str]:
    out: dict[str, str] = {}
    with open(path) as fh:
        for line in fh:
            if not line.strip():
                continue
            name, _, seq = line.rstrip("\n").partition("\t")
            out[name] = seq
    return out


def write_tbl(records: Mapping[str, str], path: str | Path) -> None:
    with open(path, "w") as fh:
        for name, seq in records.items():
            fh.write(f"{name}\t{seq}\n")
