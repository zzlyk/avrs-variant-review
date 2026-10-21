"""Decision fusion engine (spec 6.1, three layers):

Layer 1 deterministic evidence engine (BAM/VCF quantitative facts),
Layer 2 IGV visual analysis (only when a vision model actually received the
image), Layer 3 fusion into a structured, schema-validated review.

Quantitative evidence is the foundation; visual observations supplement
interpretation and surface pattern anomalies. A colorful screenshot never
overrides pileup counts.
"""
from __future__ import annotations

from typing import Any, Optional

from ..evidence.metrics import VariantEvidence
from ..triage import RuleHit, suspiciousness_score
from .providers import AIProvider, MockAIProvider
from .prompts import PROMPT_VERSION
from .schema import AIReviewResult, AIRecommendation, ReviewStatus, VisualFindings


def _risk_score(hits: list[RuleHit], recommendation: AIRecommendation) -> tuple[float, str]:
    """Deterministic, versioned rule score (0-100). NOT a calibrated probability."""
    base = suspiciousness_score(hits)
    if recommendation == AIRecommendation.SUSPECTED_FALSE_POSITIVE:
        base = max(base, 60)
    elif recommendation == AIRecommendation.INSUFFICIENT_EVIDENCE:
        base = max(base, 40)
    elif recommendation == AIRecommendation.INCONCLUSIVE:
        base = max(base, 30)
    elif recommendation == AIRecommendation.LIKELY_TRUE_VARIANT:
        base = min(base, 60)
    return float(min(100, base)), "rulescore-v1"


class ReviewEngine:
    """Fuses evidence + rules + provider output into a validated AIReviewResult."""

    def __init__(self, provider: AIProvider, config=None, thresholds: dict | None = None):
        self.provider = provider
        self.config = config
        self.thresholds = thresholds or {}
        self.vision_allowed = provider.supports_vision and (
            getattr(provider, "is_local", False)
            or (config is not None and config.ai.allow_remote_screenshot_transfer)
        )

    def review(self, ev: VariantEvidence, hits: list[RuleHit],
               task_ctx: dict[str, Any], snapshot_paths: list[str]) -> AIReviewResult:
        """task_ctx: sample_id, variant_id/uid, chrom, pos, ref, alt, variant_type,
        vcf_filter, genotype, reference_build."""
        ctx = dict(task_ctx)
        evidence_dict = ev.model_dump()
        if self.thresholds:
            evidence_dict["_thresholds"] = self.thresholds
        rules_json = [h.model_dump() for h in hits]

        images: list[str] = []
        if self.vision_allowed:
            images = [p for p in snapshot_paths]
        if snapshot_paths and not images:
            ctx["_visual_note"] = "screenshots exist but visual transfer is disabled"

        raw = self.provider.review_variant(ctx, evidence_dict, rules_json, images)
        vision_used = bool(raw.pop("_vision_used", False)) and bool(images)

        rec = self._coerce_recommendation(raw.get("ai_recommendation"))
        risk, score_ver = _risk_score(hits, rec)

        vf_raw = raw.get("visual_findings") or {}
        vf = VisualFindings(
            performed=vision_used,
            alt_support_observable=self._tri_state(vf_raw.get("alt_support_observable")),
            alt_on_independent_reads=self._tri_state(vf_raw.get("alt_on_independent_reads")),
            alt_concentrated_at_read_ends=self._tri_state(vf_raw.get("alt_concentrated_at_read_ends")),
            soft_clip_clustering=self._tri_state(vf_raw.get("soft_clip_clustering")),
            local_alignment_anomaly=self._tri_state(vf_raw.get("local_alignment_anomaly")),
            strand_support_appears_one_sided=self._tri_state(vf_raw.get("strand_support_appears_one_sided")),
            nearby_mismatch_cluster=self._tri_state(vf_raw.get("nearby_mismatch_cluster")),
            notes=[str(x) for x in vf_raw.get("notes", [])] if isinstance(vf_raw, dict) else [],
            limitations=[str(x) for x in vf_raw.get("limitations", [])] if isinstance(vf_raw, dict) else [],
        )

        requires_human = rec in (
            AIRecommendation.SUSPECTED_FALSE_POSITIVE,
            AIRecommendation.INCONCLUSIVE,
            AIRecommendation.INSUFFICIENT_EVIDENCE,
            AIRecommendation.POSSIBLE_FALSE_NEGATIVE,
        ) or risk >= 40

        status = (
            ReviewStatus.COMPLETED if vision_used else ReviewStatus.COMPLETED_NO_VISUAL
        )
        limitations = [str(x) for x in raw.get("limitations", [])] or []
        if not vision_used:
            limitations.append(
                "No vision model analyzed the screenshots; review is quantitative-evidence-only."
                if snapshot_paths else "No screenshots available for this variant.")

        return AIReviewResult(
            sample_id=ctx["sample_id"],
            variant_id=ctx["variant_id"],
            reference_build=ctx.get("reference_build", "unknown"),
            chromosome=ctx["chrom"],
            position=int(ctx["pos"]),
            ref=ctx["ref"],
            alt=ctx["alt"],
            original_vcf_status=ctx.get("vcf_filter"),
            review_status=status,
            ai_recommendation=rec,
            risk_score=risk,
            risk_score_version=score_ver,
            risk_score_is_calibrated_probability=False,
            requires_human_review=requires_human,
            evidence_summary=str(raw.get("evidence_summary", "")),
            bam_metrics={k: v for k, v in evidence_dict.items()
                         if k in ("total_reads", "usable_reads", "ref_count", "alt_count",
                                  "other_count", "vaf", "strand_forward", "strand_reverse",
                                  "mean_mapq_alt", "mean_baseq_alt", "dup_frac",
                                  "softclip_frac", "alt_start_cluster_frac", "local_depth_mean")},
            vcf_metrics={k: v for k, v in evidence_dict.items()
                         if k in ("caller_dp", "caller_ad_ref", "caller_ad_alt", "caller_gq",
                                  "caller_qual", "genotype", "vaf_expected_for_gt",
                                  "vaf_deviation", "ad_consistency")},
            visual_findings=vf,
            artifact_signals=[str(x) for x in raw.get("artifact_signals", [])],
            supporting_evidence=[str(x) for x in raw.get("supporting_evidence", [])],
            conflicting_evidence=[str(x) for x in raw.get("conflicting_evidence", [])],
            limitations=limitations,
            reason_codes=[str(x) for x in raw.get("reason_codes", [])]
            + [h.rule_id for h in hits],
            screenshot_paths=list(snapshot_paths),
            model_name=f"{getattr(self.provider, 'name', '?')}:{getattr(self.provider, 'model', '?')}",
            model_version="",
            prompt_version=PROMPT_VERSION,
            visual_review_performed=vision_used,
        )

    @staticmethod
    def _tri_state(v) -> Optional[str]:
        if v is None:
            return None
        s = str(v).strip().lower()
        if s in ("yes", "true"):
            return "yes"
        if s in ("no", "false"):
            return "no"
        return "unclear"

    @staticmethod
    def _coerce_recommendation(v) -> AIRecommendation:
        if isinstance(v, AIRecommendation):
            return v
        try:
            return AIRecommendation(str(v).strip().upper())
        except ValueError:
            return AIRecommendation.INCONCLUSIVE
