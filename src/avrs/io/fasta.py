"""Minimal FASTA reader/writer with in-memory contig index (demo-scale)."""
from __future__ import annotations

from pathlib import Path


class Fasta:
    """Loads a small FASTA fully into memory. For WGS use the pysam backend."""

    def __init__(self, path: str | Path):
        self.path = Path(path)
        self.contigs: dict[str, str] = {}
        self._load()

    def _load(self):
        name = None
        chunks: list[str] = []
        with open(self.path, "r", encoding="utf-8") as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    if name is not None:
                        self.contigs[name] = "".join(chunks)
                    name = line[1:].split()[0]
                    chunks = []
                else:
                    chunks.append(line.upper())
        if name is not None:
            self.contigs[name] = "".join(chunks)

    @property
    def contig_names(self) -> list[str]:
        return list(self.contigs)

    def contig_len(self, chrom: str) -> int:
        if chrom not in self.contigs:
            raise KeyError(f"Contig {chrom!r} not in {self.path}. Available: {self.contig_names[:10]}...")
        return len(self.contigs[chrom])

    def fetch(self, chrom: str, start: int, end: int) -> str:
        """1-based inclusive range."""
        seq = self.contigs[chrom]
        s = max(0, start - 1)
        e = min(len(seq), end)
        return seq[s:e]


def write_fasta(path: str | Path, contigs: dict[str, str], width: int = 60):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        for name, seq in contigs.items():
            fh.write(f">{name}\n")
            for i in range(0, len(seq), width):
                fh.write(seq[i:i + width] + "\n")


def write_fai(path: str | Path, contigs: dict[str, str], width: int = 60):
    """Write a .fai index for a FASTA written with fixed line width."""
    linesep = width + 1
    offset = 0
    rows = []
    for name, seq in contigs.items():
        rows.append((name, len(seq), offset))
        offset += len(seq) + (len(seq) + width - 1) // width + 1  # seq bytes + newlines + header line
    with open(str(path) + ".fai", "w", encoding="utf-8") as fh:
        for name, ln, off in rows:
            fh.write(f"{name}\t{ln}\t{off}\t{width}\t{linesep}\n")
