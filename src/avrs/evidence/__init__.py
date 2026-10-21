"""Evidence backends: pysam (preferred, Linux/WSL2) and pure-python fallback."""
from __future__ import annotations

from pathlib import Path
from typing import Optional, Protocol

from .core import ReadView, EvidenceParams, compute_evidence
from .metrics import VariantEvidence


class Backend(Protocol):
    name: str

    def evidence(self, chrom: str, pos: int, ref: str, alt: str,
                 variant_type: str, caller: dict | None = None) -> VariantEvidence: ...

    def coverage_at(self, chrom: int_pos) -> dict: ...  # type: ignore[name-defined]


def has_pysam() -> bool:
    try:
        import pysam  # noqa: F401
        return True
    except ImportError:
        return False


class PysamEvidenceBackend:
    """Primary backend on Linux/WSL2 (region queries via BAM index)."""
    name = "pysam"

    def __init__(self, bam_path: str | Path, fasta_path: str | Path,
                 params: EvidenceParams | None = None):
        import pysam
        self.params = params or EvidenceParams()
        self.bam = pysam.AlignmentFile(str(bam_path), "rb")
        self.fasta = pysam.FastaFile(str(fasta_path))
        self.bam_path = str(bam_path)
        ref = self._self_ref()

    def _self_ref(self):
        return _PysamRef(self.fasta)

    def _adapt(self, r) -> ReadView:
        import pysam
        cig = []
        for op, l in r.cigartuples or []:
            from ..io.bam import CIGAR_OPS
            cig.append((l, CIGAR_OPS[op]))
        qual = list(r.query_qualities) if r.query_qualities is not None else []
        return ReadView(
            qname=r.query_name or "", flag=r.flag, pos0=r.reference_start or 0,
            mapq=r.mapping_quality or 0, cigar=cig,
            seq=r.query_sequence or "", qual=qual,
            is_paired=bool(r.flag & 0x1),
        )

    def evidence(self, chrom, pos, ref, alt, variant_type, caller=None,
                 target_region=None) -> VariantEvidence:
        t0 = max(0, pos - 1 - 75)
        t1 = pos + max(len(ref), len(alt)) + 75
        reads = [self._adapt(r) for r in self.bam.fetch(chrom, t0, t1)]
        return compute_evidence(
            reads, chrom, pos, ref, alt, variant_type,
            self._self_ref(), self.params, backend=self.name,
            caller=caller, target_region=target_region,
        )

    def coverage_at(self, chrom: str, pos: int, window: int = 0) -> dict:
        """Base-level usable coverage stats around a position.

        count_coverage returns one array per base; the sum over the four
        arrays at a position IS the depth at that position (no extra divisor).
        """
        import pysam
        col = self.bam.count_coverage(chrom, pos - 1 - window, pos + window,
                                      quality_threshold=self.params.min_baseq)
        total = sum(sum(c) for c in col)
        n = max(1, len(col[0]))
        return {"mean_depth": total / n, "min_baseq": self.params.min_baseq}

    def callable_status(self, chrom: str, pos: int, min_depth: int) -> str:
        depth = self.coverage_at(chrom, pos)["mean_depth"]
        if depth >= min_depth:
            return "callable"
        return "insufficient_coverage"


class _PysamRef:
    def __init__(self, fasta):
        self.fasta = fasta

    def fetch(self, chrom: str, start: int, end: int) -> str:  # 1-based inclusive
        return self.fasta.fetch(chrom, start - 1, end)


class PurePythonEvidenceBackend:
    """Fallback backend (works on native Windows, demo-scale BAMs)."""
    name = "purepython"

    def __init__(self, bam_path: str | Path, fasta_path: str | Path,
                 params: EvidenceParams | None = None):
        from ..io.fasta import Fasta
        from ..io import bam as bamio
        self.params = params or EvidenceParams()
        self.bam_path = Path(bam_path)
        self.fasta = Fasta(fasta_path)
        self._bamio = bamio

    def _adapt(self, rec) -> ReadView:
        return ReadView(
            qname=rec.qname, flag=rec.flag, pos0=rec.pos, mapq=rec.mapq,
            cigar=list(rec.cigar), seq=rec.seq, qual=rec.qual,
            is_paired=bool(rec.flag & 0x1),
        )

    def evidence(self, chrom, pos, ref, alt, variant_type, caller=None,
                 target_region=None) -> VariantEvidence:
        t0 = max(0, pos - 1 - 75)
        t1 = pos + max(len(ref), len(alt)) + 75
        recs = self._bamio.fetch_region(self.bam_path, chrom, t0, t1)
        reads = [self._adapt(r) for r in recs]
        refetcher = self._FastaRef(self.fasta)
        return compute_evidence(
            reads, chrom, pos, ref, alt, variant_type,
            refetcher, self.params, backend=self.name,
            caller=caller, target_region=target_region,
        )

    def coverage_at(self, chrom: str, pos: int, window: int = 0) -> dict:
        t0 = max(0, pos - 1 - window)
        t1 = pos + window
        recs = self._bamio.fetch_region(self.bam_path, chrom, t0, t1)
        pts = range(t0, t1) if window else [pos - 1]
        depths = [sum(1 for r in recs if r.pos <= p < r.pos + r.reference_length()) for p in pts]
        return {"mean_depth": sum(depths) / max(1, len(depths))}

    def callable_status(self, chrom: str, pos: int, min_depth: int) -> str:
        depth = self.coverage_at(chrom, pos)["mean_depth"]
        return "callable" if depth >= min_depth else "insufficient_coverage"

    class _FastaRef:
        def __init__(self, fasta):
            self.fasta = fasta

        def fetch(self, chrom: str, start: int, end: int) -> str:
            return self.fasta.fetch(chrom, start, end)


def make_backend(bam_path: str | Path, fasta_path: str | Path,
                 params: EvidenceParams | None = None,
                 prefer: str = "auto") -> Backend:
    """Select the evidence backend. auto -> pysam when importable."""
    if prefer == "pysam" and not has_pysam():
        raise RuntimeError("pysam requested but not importable")
    if prefer in ("auto", "pysam") and has_pysam():
        return PysamEvidenceBackend(bam_path, fasta_path, params)
    return PurePythonEvidenceBackend(bam_path, fasta_path, params)
