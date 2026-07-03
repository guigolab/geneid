"""Read/write geneid ``.param`` files.

This module is the authority for the geneid parameter-file format. The format is
frozen: the compiled ``geneid`` binary parses these files, so byte-level format
compatibility is mandatory even as the pipeline is free to change how the values
themselves are computed.

Format rule that makes parsing robust: a *section keyword* is always a bare single
token on its own line (``^[A-Za-z][A-Za-z0-9_]*$``, no leading whitespace). Data
lines never match this — they carry spaces, digits, ``+``, ``:`` or ``.`` — so the
keyword pattern alone separates headers from data.

Round-trip guarantee: sections we do not yet model are echoed verbatim, so
``Param.read(p).to_text()`` is byte-identical to the original for every real param.
Typed setters re-render only the single data line they touch.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path

_KEYWORD_RE = re.compile(r"^[A-Za-z][A-Za-z0-9_]*$")


def _is_keyword(line: str) -> bool:
    """True if ``line`` is a section-keyword header line."""
    body = line.rstrip("\n")
    return body == body.lstrip() and bool(_KEYWORD_RE.match(body.strip()))


def _is_data(line: str) -> bool:
    body = line.strip()
    return bool(body) and not body.startswith("#")


@dataclass
class Profile:
    """A geneid site profile (PWA / order-k Markov)."""

    name: str
    header: list[float]  # length offset cutoff order [a b acc_context min_dist opt_dist pen_scale]
    rows: list[tuple[int, str, float]]  # (position, oligo, score)

    @property
    def length(self) -> int:
        return int(self.header[0])

    @property
    def order(self) -> int:
        return int(self.header[3])


@dataclass
class _Block:
    """One section: keyword header plus every raw line up to the next header.

    ``keyword`` is ``None`` for the file preamble. ``raw`` includes the keyword
    line itself and any trailing comment/blank lines, so concatenating every
    block's ``raw`` reproduces the file byte-for-byte.
    """

    keyword: str | None
    raw: list[str] = field(default_factory=list)

    def data_line_indices(self) -> list[int]:
        # skip index 0 (the keyword line) for real sections
        start = 0 if self.keyword is None else 1
        return [i for i in range(start, len(self.raw)) if _is_data(self.raw[i])]

    def first_data_index(self) -> int | None:
        idx = self.data_line_indices()
        return idx[0] if idx else None


class Param:
    """A parsed geneid parameter file."""

    def __init__(self, blocks: list[_Block]):
        self._blocks = blocks

    # ---- construction ---------------------------------------------------
    @classmethod
    def from_text(cls, text: str) -> Param:
        lines = text.splitlines(keepends=True)
        blocks: list[_Block] = []
        current = _Block(keyword=None)
        for line in lines:
            if _is_keyword(line):
                blocks.append(current)
                current = _Block(keyword=line.rstrip("\n").strip(), raw=[line])
            else:
                current.raw.append(line)
        blocks.append(current)
        # drop an empty leading preamble block only if it holds nothing at all
        if blocks and blocks[0].keyword is None and not blocks[0].raw:
            blocks.pop(0)
        return cls(blocks)

    @classmethod
    def read(cls, path: str | Path) -> Param:
        return cls.from_text(Path(path).read_text())

    # ---- serialization --------------------------------------------------
    def to_text(self) -> str:
        return "".join(line for block in self._blocks for line in block.raw)

    def write(self, path: str | Path) -> None:
        Path(path).write_text(self.to_text())

    # ---- block lookup ---------------------------------------------------
    def _find(self, keyword: str, index: int = 0) -> _Block:
        matches = [b for b in self._blocks if b.keyword == keyword]
        if not matches:
            raise KeyError(f"no section {keyword!r} in parameter file")
        try:
            return matches[index]
        except IndexError as exc:
            raise KeyError(f"section {keyword!r} occurrence {index} not found") from exc

    def has(self, keyword: str) -> bool:
        return any(b.keyword == keyword for b in self._blocks)

    def keywords(self) -> list[str]:
        return [b.keyword for b in self._blocks if b.keyword is not None]

    # ---- typed accessors ------------------------------------------------
    def scalar(self, keyword: str, index: int = 0) -> str:
        """Return the first data line of a section as a stripped string."""
        block = self._find(keyword, index)
        di = block.first_data_index()
        if di is None:
            raise ValueError(f"section {keyword!r} has no data line")
        return block.raw[di].strip()

    def vector(self, keyword: str, index: int = 0) -> list[str]:
        return self.scalar(keyword, index).split()

    def set_scalar(self, keyword: str, value: object, index: int = 0) -> None:
        """Replace the single data line of a section, preserving line ending."""
        block = self._find(keyword, index)
        di = block.first_data_index()
        if di is None:
            raise ValueError(f"section {keyword!r} has no data line to set")
        ending = "\n" if block.raw[di].endswith("\n") else ""
        block.raw[di] = f"{value}{ending}"

    def replace_block_data(self, keyword: str, new_lines: list[str], index: int = 0) -> None:
        """Replace a section's data lines with ``new_lines`` (each without a
        trailing newline), keeping the keyword line and any leading/trailing
        comment or blank lines. Used to swap in a freshly trained profile or
        Markov matrix while preserving the file's structure around it.
        """
        block = self._find(keyword, index)
        di = block.data_line_indices()
        if not di:
            raise ValueError(f"section {keyword!r} has no data to replace")
        first, last = di[0], di[-1]
        head = block.raw[:first]
        tail = block.raw[last + 1 :]
        block.raw = head + [ln + "\n" for ln in new_lines] + tail

    @property
    def num_isochores(self) -> int:
        return int(self.scalar("number_of_isochores"))

    def profile_names(self) -> list[str]:
        return [b.keyword for b in self._blocks if b.keyword and b.keyword.endswith("_profile")]

    def profile(self, name: str, index: int = 0) -> Profile:
        block = self._find(name, index)
        di = block.data_line_indices()
        if not di:
            raise ValueError(f"profile {name!r} has no data")
        header = [float(x) for x in block.raw[di[0]].split()]
        rows: list[tuple[int, str, float]] = []
        for i in di[1:]:
            parts = block.raw[i].split()
            rows.append((int(parts[0]), parts[1], float(parts[2])))
        return Profile(name=name, header=header, rows=rows)
