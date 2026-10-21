"""Configurable triage rule engine (spec 4.1).

Every hit carries: rule_id, observed_value, configured_threshold, reason,
severity. Rules never declare a variant true/false - they flag risk signals
and feed prioritization + the AI review.
"""
from __future__ import annotations

from typing import Any, Optional

from pydantic import BaseModel, Field

from ..config import TriageProfile
from ..evidence.metrics import VariantEvidence


class RuleHit(BaseModel):
    rule_id: str
    observed_value: Any = None
    configured_threshold: Any = None
    reason: str
    severity: str = "warning"  # info | warning | high


def _gq(ev: VariantEvidence) -> Optional[int]:
    return ev.caller_gq


def _rules_for(ev: VariantEvidence, profile: TriageProfile) -> list:
    p = profile
    R = []

    def add(rule_id, obs, thr, reason, sev="warning", enabled=True):
        if obs is not None and enabled:
            R.append(RuleHit(rule_id=rule_id, observed_value=obs,
                             configured_threshold=thr, reason=reason, severity=sev))

    # 1 low DP
    add("LOW_DP", ev.total_reads if ev.total_reads is not None else None,
        p.min_dp,
        f"Total read depth {ev.total_reads} below configured minimum {p.min_dp}",
        "warning", ev.total_reads is not None and ev.total_reads < p.min_dp)

    # 2 low GQ
    gq = _gq(ev)
    add("LOW_GQ", gq, p.min_gq,
        f"Genotype quality {gq} below configured minimum {p.min_gq}",
        "warning", gq is not None and gq < p.min_gq)

    # 3 low MQ
    add("LOW_MQ", ev.mean_mapq_alt, p.min_mapq,
        f"Mean mapping quality of ALT reads {ev.mean_mapq_alt} below {p.min_mapq}"
        if ev.mean_mapq_alt is not None else None,
        "warning", ev.mean_mapq_alt is not None and ev.mean_mapq_alt < p.min_mapq)

    # 3b every overlapping read below the MAPQ threshold (locus-wide quality issue)
    add("ALL_READS_LOW_MQ", ev.low_mapq_reads, ev.total_reads,
        f"All {ev.total_reads} overlapping reads are below MAPQ {p.min_mapq}"
        if ev.total_reads else None,
        "warning", bool(ev.total_reads) and ev.low_mapq_reads >= ev.total_reads)

    # 4 insufficient ALT reads
    add("LOW_ALT_AD", ev.alt_count, p.min_alt_ad,
        f"Only {ev.alt_count} usable ALT-supporting reads (< {p.min_alt_ad})",
        "high", ev.alt_count is not None and ev.alt_count < p.min_alt_ad)

    # 5 strand bias
    sri = ev.strand_ratio_imbalance
    add("STRAND_BIAS", sri, p.max_strand_ratio_imbalance,
        f"{sri:.2f} of ALT reads on one strand (threshold {p.max_strand_ratio_imbalance})"
        if sri is not None else None,
        "high", sri is not None and sri >= p.max_strand_ratio_imbalance)

    # 6 low base quality of ALT
    add("LOW_ALT_BQ", ev.mean_baseq_alt, p.min_baseq,
        f"Mean ALT base quality {ev.mean_baseq_alt} below {p.min_baseq}"
        if ev.mean_baseq_alt is not None else None,
        "warning", ev.mean_baseq_alt is not None and ev.mean_baseq_alt < p.min_baseq)

    # 6b ALT reads below the base-quality threshold (reported, not filtered)
    add("ALT_READS_LOW_BQ", ev.alt_reads_below_baseq, 0,
        f"{ev.alt_reads_below_baseq} ALT reads below base-quality threshold {p.min_baseq}",
        "warning", bool(ev.alt_reads_below_baseq))

    # 7 ALT near read ends
    add("ALT_NEAR_READ_END", ev.alt_reads_within_end_dist, p.min_dist_from_read_end,
        f"{ev.alt_reads_within_end_dist} ALT reads within {p.min_dist_from_read_end}bp of read end",
        "warning", ev.alt_reads_within_end_dist is not None and ev.alt_reads_within_end_dist > 0
        and (ev.alt_count and ev.alt_reads_within_end_dist / ev.alt_count >= 0.5))

    # 8 soft clipping
    add("SOFT_CLIPPING", ev.softclip_frac, p.max_softclip_frac,
        f"Soft-clipped fraction near locus {ev.softclip_frac:.2f} > {p.max_softclip_frac}"
        if ev.softclip_frac is not None else None,
        "warning", ev.softclip_frac is not None and ev.softclip_frac > p.max_softclip_frac)

    # 9 start-position clustering
    add("START_POS_CLUSTER", ev.alt_start_cluster_frac, p.max_start_pos_cluster_frac,
        f"{ev.alt_start_cluster_frac:.2f} of ALT reads share an alignment start"
        if ev.alt_start_cluster_frac is not None else None,
        "high", ev.alt_start_cluster_frac is not None and ev.alt_start_cluster_frac > p.max_start_pos_cluster_frac)

    # 10 duplicates
    add("HIGH_DUP_RATE", ev.dup_frac, p.max_dup_frac,
        f"Duplicate fraction {ev.dup_frac:.2f} > {p.max_dup_frac}"
        if ev.dup_frac is not None else None,
        "warning", ev.dup_frac is not None and ev.dup_frac > p.max_dup_frac)

    # 11 indel in repeat/homopolymer
    if ev.variant_type in ("INS", "DEL"):
        hp = ev.homopolymer_len_spanning
        add("INDEL_IN_REPEAT", hp, p.homopolymer_len,
            f"Indel within/near homopolymer of length {hp} (>= {p.homopolymer_len})"
            if hp is not None else None,
            "warning", hp is not None and hp >= p.homopolymer_len)

    # 12 low complexity window
    hp = ev.homopolymer_len_spanning
    add("LOW_COMPLEXITY", hp, p.low_complexity_window,
        None, "info", hp is not None and hp >= p.low_complexity_window)
    if hp is not None and hp >= p.low_complexity_window:
        R[-1].reason = f"Local homopolymer run {hp}bp indicates low-complexity context"

    # 13 multi-mapping suspicion: very low min MQ among alt reads
    add("MULTI_MAPPING_RISK", ev.min_mapq_alt, 10,
        f"Minimum ALT mapping quality {ev.min_mapq_alt} suggests multi-mapping reads"
        if ev.min_mapq_alt is not None else None,
        "high", ev.min_mapq_alt is not None and ev.min_mapq_alt < 10)

    # 14 VCF AD vs pileup mismatch
    add("AD_PILEUP_MISMATCH", ev.ad_consistency, True,
        "VCF AD conflicts with recomputed pileup evidence",
        "high", ev.ad_consistency is False)

    # 15 abnormal VAF for genotype
    add("VAF_GT_MISMATCH", ev.vaf_deviation, p.max_vaf_dev,
        f"|VAF {ev.vaf:.2f} - expected for GT {ev.vaf_expected_for_gt:.2f}| = {ev.vaf_deviation:.2f} > {p.max_vaf_dev}"
        if (ev.vaf is not None and ev.vaf_deviation is not None) else None,
        "high", ev.vaf_deviation is not None and ev.vaf_deviation > p.max_vaf_dev)

    # 16 zero depth
    add("NO_COVERAGE", ev.total_reads, 1,
        "No reads overlap the variant position",
        "high", ev.total_reads == 0)

    # 17 local coverage anomaly
    dr = ev.depth_ratio
    add("LOCAL_COVERAGE_ANOMALY", dr, (0.5, 2.0),
        f"Variant depth/local depth ratio {dr:.2f} outside [0.5, 2.0]"
        if dr is not None else None,
        "warning", dr is not None and (dr < 0.5 or dr > 2.0))

    # 18 pair support inconsistency
    add("PAIR_SUPPORT_INCONSISTENT", ev.pair_concordance, 1.0,
        f"Only {ev.concordant_pairs} of {ev.concordant_pairs + ev.discordant_pairs} ALT read-pairs concordant"
        if ev.pair_concordance is not None else None,
        "warning", ev.pair_concordance is not None and ev.concordant_pairs + ev.discordant_pairs >= 2
        and ev.pair_concordance < 1.0)

    # 19 target region edge
    add("TARGET_REGION_EDGE", ev.dist_to_target_edge, p.min_target_edge_dist,
        f"Variant {ev.dist_to_target_edge}bp from target-region edge (< {p.min_target_edge_dist})"
        if ev.dist_to_target_edge is not None else None,
        "info", ev.dist_to_target_edge is not None and ev.dist_to_target_edge < p.min_target_edge_dist)

    # 20 caller QUAL low
    add("LOW_QUAL", ev.caller_qual, 30,
        f"Caller QUAL {ev.caller_qual} < 30"
        if ev.caller_qual is not None else None,
        "info", ev.caller_qual is not None and ev.caller_qual < 30)

    return R


def evaluate_rules(ev: VariantEvidence, profile: TriageProfile) -> list[RuleHit]:
    """Evaluate all rules; returns hits (only rules that fired)."""
    return _rules_for(ev, profile)


SEVERITY_WEIGHT = {"high": 10, "warning": 4, "info": 1}


def suspiciousness_score(hits: list[RuleHit]) -> int:
    """Deterministic risk score (0-100, NOT a validated probability)."""
    raw = sum(SEVERITY_WEIGHT.get(h.severity, 0) for h in hits)
    return min(100, raw)


def prioritize(hits: list[RuleHit], in_high_priority_region: bool = False,
               vcf_filter_failed: bool = False) -> str:
    """HIGH / MEDIUM / LOW review priority.

    Uses: evidence quality anomalies, conflicting data, region importance
    (user-designated), suspicious filter records. A variant's clinical
    importance (e.g. ClinVar) never *reduces* its evidence review.

    Any rule hit at all puts a variant into at least MEDIUM - otherwise a
    single clear artifact signal would never reach a human.
    """
    if not hits:
        return "MEDIUM" if (in_high_priority_region or vcf_filter_failed) else "LOW"
    n_high = sum(1 for h in hits if h.severity == "high")
    n_warn = sum(1 for h in hits if h.severity == "warning")
    if n_high >= 1:
        return "HIGH"
    if n_warn >= 1:
        return "MEDIUM"
    return "MEDIUM" if (in_high_priority_region or vcf_filter_failed) else "LOW"


REVIEW_MODES = ("all", "suspicious_only", "target_genes", "filtered_calls",
                "bed_regions", "selected_variants")


def should_review(mode: str, hits: list[RuleHit], priority: str, gene: Optional[str],
                  target_genes: list[str], vcf_filter_failed: bool,
                  in_bed_region: bool, variant_uid: Optional[str],
                  selected_uids: list[str]) -> bool:
    """Review-mode gating (spec 4.2)."""
    if mode == "all":
        return True
    if mode == "suspicious_only":
        return priority in ("HIGH", "MEDIUM")
    if mode == "target_genes":
        return bool(gene and gene in target_genes)
    if mode == "filtered_calls":
        return vcf_filter_failed
    if mode == "bed_regions":
        return in_bed_region
    if mode == "selected_variants":
        return bool(variant_uid and variant_uid in selected_uids)
    raise ValueError(f"Unknown review mode: {mode}")
