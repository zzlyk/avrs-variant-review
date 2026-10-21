"""Evidence metric models. The schema deliberately separates:

- read depth (raw overlapping reads)
- usable depth (reads passing MAPQ/BASEQ/dup/secondary filters)
- VCF AD (caller-reported allele depths)
- pileup-derived ref/alt/other counts

These are different definitions and are never conflated (spec 4.0).
"""
from __future__ import annotations

from typing import Optional

from pydantic import BaseModel, Field


class VariantEvidence(BaseModel):
    chrom: str
    pos: int                      # 1-based, normalized
    ref: str
    alt: str
    variant_type: str

    total_reads: int = 0          # raw reads overlapping the variant span
    usable_reads: int = 0         # after MAPQ/BASEQ/dup/secondary filters
    ref_count: int = 0            # usable reads supporting REF
    alt_count: int = 0            # usable reads supporting ALT
    other_count: int = 0          # usable reads supporting neither
    vaf: Optional[float] = None   # alt / (ref+alt+other) among usable reads

    # strand support (ALT reads)
    strand_forward: int = 0
    strand_reverse: int = 0
    strand_ratio_imbalance: Optional[float] = None  # max(fwd,rev)/alt_count

    # quality
    mean_mapq_alt: Optional[float] = None
    min_mapq_alt: Optional[int] = None
    mean_baseq_alt: Optional[float] = None
    min_baseq_alt: Optional[int] = None
    alt_reads_below_baseq: int = 0    # ALT reads whose base quality is below the configured minimum
    low_mapq_reads: int = 0           # reads filtered out because MAPQ < configured minimum

    # positional
    mean_dist_from_read_ends: Optional[float] = None
    alt_reads_within_end_dist: int = 0   # ALT reads closer than threshold to read end

    # artifact signals
    soft_clipped_reads_near_locus: int = 0
    softclip_frac: Optional[float] = None
    duplicate_reads: int = 0
    dup_frac: Optional[float] = None
    alt_start_cluster_frac: Optional[float] = None  # max frac of ALT reads with identical alignment start
    discordant_pairs: int = 0
    concordant_pairs: int = 0
    pair_concordance: Optional[float] = None

    # local context
    local_depth_mean: Optional[float] = None  # mean usable depth in +/- window
    depth_ratio: Optional[float] = None       # depth(variant)/local_depth_mean
    homopolymer_len_spanning: Optional[int] = None

    # cross-check against caller
    caller_dp: Optional[int] = None
    caller_ad_ref: Optional[int] = None
    caller_ad_alt: Optional[int] = None
    caller_gq: Optional[int] = None
    caller_qual: Optional[float] = None
    ad_consistency: Optional[bool] = None     # VCF AD vs recomputed pileup
    genotype: Optional[str] = None
    vaf_expected_for_gt: Optional[float] = None
    vaf_deviation: Optional[float] = None

    # region
    in_target_region: Optional[bool] = None
    dist_to_target_edge: Optional[int] = None

    warnings: list[str] = Field(default_factory=list)
    backend: str = "unknown"
