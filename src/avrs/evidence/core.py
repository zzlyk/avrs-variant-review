"""Shared allele-aware evidence computation.

Both evidence backends (pysam on Linux/WSL, pure-python elsewhere) adapt their
reads to `ReadView` and call `compute_evidence`, so the correctness-critical
logic exists exactly once. Indels are handled via CIGAR-aware allele
reconstruction, not single-base pileup (spec 4.0).
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Optional, Protocol

from .metrics import VariantEvidence

SEQ_NT16 = "=ACMGRSVTWYHKDBN"

# SAM flags
FREVERSE = 16
FSECONDARY = 256
FQCFAIL = 512
FDUP = 1024
FSUPPLEMENTARY = 2048
FPROPER_PAIR = 2


@dataclass
class ReadView:
    """Minimal read representation shared by all backends."""
    qname: str
    flag: int
    pos0: int                       # 0-based leftmost aligned position
    mapq: int
    cigar: list[tuple[int, str]]    # [(len, op)]
    seq: str
    qual: list[int]                 # Phred qualities, per base
    is_paired: bool = False

    def ref_end(self) -> int:       # 0-based exclusive
        return self.pos0 + sum(l for l, op in self.cigar if op in "MDN=X")


class RefFetcher(Protocol):
    def fetch(self, chrom: str, start: int, end: int) -> str: ...  # 1-based inclusive


@dataclass
class EvidenceParams:
    min_mapq: int = 20
    min_baseq: int = 13
    min_dist_from_read_end: int = 5
    local_window: int = 50
    homopolymer_scan_window: int = 20
    dup_flag_excludes: bool = True


def _read_allele_over_span(read: ReadView, t0: int, tlen: int) -> Optional[str]:
    """Reconstruct the read's allele across reference span [t0, t0+tlen).

    - M/=/X contribute their base if within the span.
    - D/N within the span contribute nothing (deletion of reference bases).
    - I contributes its inserted bases when the junction sits exactly at
      t0 (the VCF insertion anchor convention: insertion follows the anchor
      base at 0-based t0).
    Returns None if the read does not fully span the region with aligned
    bases (i.e. the read is uninformative for this variant).
    """
    ref = read.pos0
    q = 0
    out: list[str] = []
    covered = 0
    for op_len, op in read.cigar:
        if op in "M=X":
            for k in range(op_len):
                r = ref + k
                if t0 <= r < t0 + tlen:
                    idx = q + k
                    if idx >= len(read.seq):
                        return None
                    out.append(read.seq[idx].upper())
                    covered += 1
            ref += op_len
            q += op_len
        elif op == "I":
            # insertion right after the base at ref-1; junction == ref
            if ref == t0 + 1 and t0 <= ref - 1 < t0 + tlen:
                piece = read.seq[q:q + op_len]
                if len(piece) < op_len:
                    return None
                out.append(piece.upper())
            q += op_len
        elif op in "DN":
            for k in range(op_len):
                r = ref + k
                if t0 <= r < t0 + tlen:
                    if op == "N":
                        return None  # spliced gap: uninformative
                    covered += 1  # deleted base: contributes nothing
            ref += op_len
        elif op in "SH":
            if op == "H":
                continue  # hard clip has no query bases
            q += op_len  # soft clip consumes query only
        elif op == "P":
            continue
    if covered < tlen:
        return None
    return "".join(out)


def _base_qual_at(read: ReadView, t0: int) -> int:
    for op_len, op in read.cigar:
        pass
    # position of base aligned at t0 within query
    ref = read.pos0
    q = 0
    for op_len, op in read.cigar:
        if op in "M=X":
            for k in range(op_len):
                if ref + k == t0:
                    return read.qual[q + k] if q + k < len(read.qual) else 0
            ref += op_len
            q += op_len
        elif op in "DN":
            ref += op_len
        elif op in "IS":
            q += op_len
    return 0


def _alt_base_qual(read: ReadView, t0: int, tlen: int) -> int:
    """Mean quality of the base(s) evidencing the ALT at the locus."""
    quals: list[int] = []
    ref = read.pos0
    q = 0
    for op_len, op in read.cigar:
        if op in "M=X":
            for k in range(op_len):
                if t0 <= ref + k < t0 + tlen and q + k < len(read.qual):
                    quals.append(read.qual[q + k])
            ref += op_len
            q += op_len
        elif op == "I":
            if ref == t0 + 1:
                for k in range(op_len):
                    if q + k < len(read.qual):
                        quals.append(read.qual[q + k])
            q += op_len
        elif op in "DN":
            ref += op_len
        elif op == "S":
            q += op_len
    return sum(quals) // max(1, len(quals)) if quals else 0


def _dist_from_ends(read: ReadView, t0: int) -> int:
    n_aligned = sum(l for l, op in read.cigar if op in "MIS=X")
    left = t0 - read.pos0
    right = n_aligned - 1 - left
    return max(0, min(left, right))


def _has_soft_clip(read: ReadView) -> bool:
    return any(op == "S" and l > 0 for l, op in read.cigar)


def compute_evidence(
    reads: list[ReadView],
    chrom: str,
    pos: int,
    ref: str,
    alt: str,
    variant_type: str,
    refetcher: RefFetcher,
    params: EvidenceParams,
    backend: str,
    caller: Optional[dict] = None,
    target_region: Optional[tuple[int, int]] = None,
) -> VariantEvidence:
    """caller: optional dict with DP / AD / GQ / QUAL / GT from the VCF."""
    t0 = pos - 1            # 0-based anchor
    tlen = len(ref)
    ref_u = ref.upper()
    alt_u = alt.upper()

    ev = VariantEvidence(
        chrom=chrom, pos=pos, ref=ref, alt=alt, variant_type=variant_type, backend=backend
    )
    ev.total_reads = len(reads)

    usable: list[ReadView] = []
    alt_reads: list[ReadView] = []
    alt_quals: list[int] = []
    alt_mapqs: list[int] = []
    alt_strands: list[int] = []   # 1 forward, -1 reverse
    alt_dists: list[int] = []
    start_positions: list[int] = []
    dup_count = 0
    softclip_count = 0
    low_bq_alt_reads = 0
    low_mapq_reads = 0
    usable_alleles: list[str] = []

    for r in reads:
        if r.flag & (FSECONDARY | FSUPPLEMENTARY):
            continue
        if params.dup_flag_excludes and (r.flag & FDUP):
            dup_count += 1
            continue
        allele = _read_allele_over_span(r, t0, tlen)
        if allele is None:
            continue
        # Base quality is REPORTED, not used as a hard filter: a QC system
        # must be able to show "ALT support exists but at low base quality"
        # instead of silently dropping the reads (which would look like
        # "no evidence at all" and hide the artifact signal).
        bq = _alt_base_qual(r, t0, tlen)
        if r.mapq < params.min_mapq:
            low_mapq_reads += 1
            continue
        usable.append(r)
        usable_alleles.append(allele)
        if allele == alt_u:
            alt_reads.append(r)
            alt_quals.append(bq)
            alt_mapqs.append(r.mapq)
            alt_strands.append(-1 if (r.flag & FREVERSE) else 1)
            alt_dists.append(_dist_from_ends(r, t0))
            start_positions.append(r.pos0)
            if bq < params.min_baseq:
                low_bq_alt_reads += 1
        elif allele == ref_u:
            ev.ref_count += 1
        else:
            ev.other_count += 1
        if _has_soft_clip(r):
            softclip_count += 1

    ev.alt_count = len(alt_reads)
    ev.alt_reads_below_baseq = low_bq_alt_reads
    ev.low_mapq_reads = low_mapq_reads
    ev.usable_reads = ev.ref_count + ev.alt_count + ev.other_count
    ev.duplicate_reads = dup_count
    ev.dup_frac = dup_count / (dup_count + ev.usable_reads) if (dup_count + ev.usable_reads) else None
    ev.soft_clipped_reads_near_locus = softclip_count
    ev.softclip_frac = softclip_count / len(reads) if reads else None

    denom = ev.ref_count + ev.alt_count + ev.other_count
    ev.vaf = (ev.alt_count / denom) if denom else None

    if alt_reads:
        ev.strand_forward = sum(1 for s in alt_strands if s == 1)
        ev.strand_reverse = sum(1 for s in alt_strands if s == -1)
        ev.strand_ratio_imbalance = max(ev.strand_forward, ev.strand_reverse) / ev.alt_count
        ev.mean_mapq_alt = sum(alt_mapqs) / len(alt_mapqs)
        ev.min_mapq_alt = min(alt_mapqs)
        ev.mean_baseq_alt = sum(alt_quals) / len(alt_quals)
        ev.min_baseq_alt = min(alt_quals)
        ev.mean_dist_from_read_ends = sum(alt_dists) / len(alt_dists)
        ev.alt_reads_within_end_dist = sum(1 for d in alt_dists if d < params.min_dist_from_read_end)
        if start_positions:
            from collections import Counter
            top = Counter(start_positions).most_common(1)[0][1]
            ev.alt_start_cluster_frac = top / len(start_positions)
        # paired-end consistency: only qnames that actually have two mates
        # are evaluated (single-end or unpaired reads are not penalised)
        read_alleles: dict[str, list[str]] = {}
        for rr, allele in zip(usable, usable_alleles):
            read_alleles.setdefault(rr.qname, []).append(allele)
        for qname, alleles in read_alleles.items():
            if len(alleles) < 2:
                continue
            alt_mates = sum(1 for a in alleles if a == alt_u)
            if alt_mates == len(alleles):
                ev.concordant_pairs += 1
            elif alt_mates > 0:
                ev.discordant_pairs += 1
        pair_total = ev.concordant_pairs + ev.discordant_pairs
        ev.pair_concordance = (ev.concordant_pairs / pair_total) if pair_total else None

    # local depth (usable reads overlapping sampled points in +/- window)
    win = params.local_window
    sample_points = [t0 - win, t0 - win // 2, t0 + tlen // 2, t0 + win // 2, t0 + win]
    depths = []
    for p in sample_points:
        d = sum(1 for r in usable if r.pos0 <= p < r.ref_end())
        depths.append(d)
    ev.local_depth_mean = sum(depths) / len(depths) if depths else None
    if ev.local_depth_mean:
        ev.depth_ratio = (denom / ev.local_depth_mean) if ev.local_depth_mean else None

    # local sequence complexity: longest homopolymer overlapping/adjacent to variant
    scan_start = max(1, pos - params.homopolymer_scan_window)
    scan_end = pos + len(ref) + params.homopolymer_scan_window
    try:
        seq = refetcher.fetch(chrom, scan_start, scan_end).upper()
        best = 1
        run = 1
        for i in range(1, len(seq)):
            if seq[i] == seq[i - 1]:
                run += 1
                best = max(best, run)
            else:
                run = 1
        ev.homopolymer_len_spanning = best
    except Exception:
        ev.warnings.append("reference_unavailable_for_context")

    # cross-check with caller
    if caller:
        ev.caller_dp = caller.get("DP")
        ev.caller_ad_ref = caller.get("AD_REF")
        ev.caller_ad_alt = caller.get("AD_ALT")
        ev.caller_gq = caller.get("GQ")
        ev.caller_qual = caller.get("QUAL")
        ev.genotype = caller.get("GT")
        if ev.caller_ad_ref is not None and ev.caller_ad_alt is not None and denom:
            ad_total = (ev.caller_ad_ref or 0) + (ev.caller_ad_alt or 0)
            if ad_total:
                tol = max(3, 0.25 * ad_total)
                diff = abs(ev.caller_ad_alt - ev.alt_count)
                ev.ad_consistency = diff <= tol
            else:
                ev.ad_consistency = None
        gt = caller.get("GT")
        if gt:
            import re as _re
            alleles = [a for a in _re.split(r"[/|]", gt) if a != "."]
            if alleles:
                n_alt_calls = sum(1 for a in alleles if a != "0")
                if len(alleles) == 2:
                    if n_alt_calls == 2:
                        ev.vaf_expected_for_gt = 1.0
                    elif n_alt_calls == 1:
                        ev.vaf_expected_for_gt = 0.5
                    else:
                        ev.vaf_expected_for_gt = 0.0
                    if ev.vaf is not None:
                        ev.vaf_deviation = abs(ev.vaf - ev.vaf_expected_for_gt)

    if target_region:
        rs, re_ = target_region
        if rs <= pos <= re_:
            ev.in_target_region = True
            ev.dist_to_target_edge = min(abs(pos - rs), abs(pos - re_))
        else:
            ev.in_target_region = False

    if ev.total_reads == 0:
        ev.warnings.append("no_reads_overlap_variant")
    return ev
