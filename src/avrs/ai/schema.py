"""Strict Pydantic schema for AI review output (spec 6.4).

These labels are AI *recommendations*, not truth-validated TP/FP/FN calls.
The risk score is a deterministic rule-based score (versioned), NOT a
calibrated probability - it must never be reported as such.
"""
from __future__ import annotations

import datetime as _dt
from enum import Enum
from typing import Any, Optional

from pydantic import BaseModel, ConfigDict, Field, field_validator


class AIRecommendation(str, Enum):
    LIKELY_TRUE_VARIANT = "LIKELY_TRUE_VARIANT"
    SUSPECTED_FALSE_POSITIVE = "SUSPECTED_FALSE_POSITIVE"
    INCONCLUSIVE = "INCONCLUSIVE"
    INSUFFICIENT_EVIDENCE = "INSUFFICIENT_EVIDENCE"
    REVIEW_FAILED = "REVIEW_FAILED"
    POSSIBLE_FALSE_NEGATIVE = "POSSIBLE_FALSE_NEGATIVE"  # FN-discovery candidates


class ReviewStatus(str, Enum):
    COMPLETED = "COMPLETED"
    COMPLETED_NO_VISUAL = "COMPLETED_NO_VISUAL"   # quantitative evidence only
    COMPLETED_DEGRADED = "COMPLETED_DEGRADED"     # visual attempted but unusable
    FAILED = "FAILED"
    SKIPPED = "SKIPPED"


class VisualFindings(BaseModel):
    model_config = ConfigDict(extra="forbid")
    performed: bool = False
    alt_support_observable: Optional[str] = None      # yes | no | unclear
    alt_on_independent_reads: Optional[str] = None
    alt_concentrated_at_read_ends: Optional[str] = None
    soft_clip_clustering: Optional[str] = None
    local_alignment_anomaly: Optional[str] = None
    strand_support_appears_one_sided: Optional[str] = None
    nearby_mismatch_cluster: Optional[str] = None
    notes: list[str] = Field(default_factory=list)
    limitations: list[str] = Field(default_factory=list)


class AIReviewResult(BaseModel):
    model_config = ConfigDict(extra="forbid")
    # identity
    sample_id: str
    variant_id: str
    reference_build: str
    chromosome: str
    position: int
    ref: str
    alt: str
    original_vcf_status: Optional[str] = None
    review_status: ReviewStatus = ReviewStatus.COMPLETED
    # recommendation
    ai_recommendation: AIRecommendation
    risk_score: float = Field(ge=0, le=100)
    risk_score_version: str = "rulescore-v1"
    risk_score_is_calibrated_probability: bool = False
    requires_human_review: bool = True
    # evidence
    evidence_summary: str
    bam_metrics: dict[str, Any]
    vcf_metrics: dict[str, Any]
    visual_findings: VisualFindings = Field(default_factory=VisualFindings)
    artifact_signals: list[str] = Field(default_factory=list)
    supporting_evidence: list[str] = Field(default_factory=list)
    conflicting_evidence: list[str] = Field(default_factory=list)
    limitations: list[str] = Field(default_factory=list)
    reason_codes: list[str] = Field(default_factory=list)
    screenshot_paths: list[str] = Field(default_factory=list)
    # provenance
    model_name: str
    model_version: str = ""
    prompt_version: str = ""
    review_timestamp: str = Field(default_factory=lambda: _dt.datetime.now(_dt.timezone.utc).isoformat())
    visual_review_performed: bool = False

    @field_validator("position")
    @classmethod
    def _pos_positive(cls, v: int) -> int:
        if v < 1:
            raise ValueError(f"position must be 1-based, got {v}")
        return v
