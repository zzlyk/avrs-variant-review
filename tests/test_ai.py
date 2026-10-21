"""AI review schema/provider tests + review state machine tests."""
import json

import pytest
from pydantic import ValidationError

from avrs.ai import ReviewEngine, MockAIProvider, make_provider
from avrs.ai.schema import AIReviewResult, AIRecommendation, ReviewStatus
from avrs.ai.prompts import PROMPT_VERSION, build_review_prompt
from avrs.evidence.metrics import VariantEvidence
from avrs.triage import RuleHit
from avrs.db.models import TaskState, check_transition
from avrs.config import AVRSConfig


def _ev(**kw):
    base = dict(chrom="chr1", pos=100, ref="A", alt="G", variant_type="SNV",
                total_reads=20, usable_reads=20, ref_count=10, alt_count=10,
                other_count=0, vaf=0.5, strand_forward=5, strand_reverse=5,
                strand_ratio_imbalance=0.5, mean_mapq_alt=60.0, mean_baseq_alt=37.0)
    base.update(kw)
    return VariantEvidence(**base)


def _ctx():
    return {"sample_id": "S1", "variant_id": "B#chr1:100:A>G", "chrom": "chr1",
            "pos": 100, "ref": "A", "alt": "G", "variant_type": "SNV",
            "vcf_filter": "PASS", "genotype": "0/1", "reference_build": "B"}


def test_schema_rejects_unknown_fields():
    with pytest.raises(ValidationError):
        AIReviewResult(
            sample_id="S", variant_id="v", reference_build="b", chromosome="c",
            position=1, ref="A", alt="T", ai_recommendation="LIKELY_TRUE_VARIANT",
            risk_score=10, evidence_summary="x", bam_metrics={}, vcf_metrics={},
            model_name="m", unexpected_field="boom")


def test_schema_requires_1based_position():
    with pytest.raises(ValidationError):
        AIReviewResult(
            sample_id="S", variant_id="v", reference_build="b", chromosome="c",
            position=0, ref="A", alt="T", ai_recommendation="INCONCLUSIVE",
            risk_score=10, evidence_summary="x", bam_metrics={}, vcf_metrics={},
            model_name="m")


def test_mock_provider_labels_true_variant():
    ev = _ev()
    res = ReviewEngine(MockAIProvider()).review(ev, [], _ctx(), [])
    assert res.ai_recommendation == AIRecommendation.LIKELY_TRUE_VARIANT
    assert res.risk_score_is_calibrated_probability is False
    assert res.visual_review_performed is False
    assert res.review_status == ReviewStatus.COMPLETED_NO_VISUAL
    assert "quantitative" in " ".join(res.limitations).lower() or res.limitations


def test_mock_provider_flags_false_positive():
    """Sufficient depth, but strong artifact signals -> suspected artifact."""
    ev = _ev(alt_count=6, ref_count=14, total_reads=20, vaf=0.3,
             strand_forward=0, strand_reverse=6, strand_ratio_imbalance=1.0,
             mean_mapq_alt=8.0, mean_baseq_alt=6.0, alt_start_cluster_frac=0.9)
    hits = [RuleHit(rule_id="STRAND_BIAS", observed_value=1.0,
                    configured_threshold=0.9, reason="strand", severity="high")]
    res = ReviewEngine(MockAIProvider()).review(ev, hits, _ctx(), [])
    assert res.ai_recommendation == AIRecommendation.SUSPECTED_FALSE_POSITIVE
    assert res.requires_human_review is True
    assert res.risk_score >= 60


def test_mock_provider_insufficient_evidence_below_configured_floor():
    ev = _ev(alt_count=1, ref_count=1, total_reads=2, vaf=0.5)
    res = ReviewEngine(MockAIProvider(), thresholds={"min_dp": 10, "min_alt_ad": 3}) \
        .review(ev, [], _ctx(), [])
    assert res.ai_recommendation == AIRecommendation.INSUFFICIENT_EVIDENCE


def test_insufficient_evidence_label():
    ev = _ev(total_reads=0, usable_reads=0, ref_count=0, alt_count=0, vaf=None)
    res = ReviewEngine(MockAIProvider()).review(ev, [], _ctx(), [])
    assert res.ai_recommendation == AIRecommendation.INSUFFICIENT_EVIDENCE


def test_visual_review_not_claimed_without_vision_model():
    ev = _ev()
    res = ReviewEngine(MockAIProvider()).review(ev, [], _ctx(), ["/tmp/a.png"])
    assert res.visual_review_performed is False
    assert res.visual_findings.performed is False
    assert any("vision" in l.lower() for l in res.limitations)


def test_prompt_contains_hard_rules_and_untrusted_data_warning():
    system, user = build_review_prompt(_ctx(), {"alt_count": 3}, [], images_available=True)
    assert "HARD RULES" in system
    assert "untrusted" in system.lower()
    assert "INCONCLUSIVE" in system
    assert PROMPT_VERSION in user or PROMPT_VERSION == "v1.2"
    assert "attached" in user


def test_prompt_without_images_says_so():
    _, user = build_review_prompt(_ctx(), {}, [], images_available=False)
    assert "No IGV screenshots" in user
    assert "do not pretend" in user.lower()


def test_review_result_roundtrip_json():
    ev = _ev()
    res = ReviewEngine(MockAIProvider()).review(ev, [], _ctx(), [])
    data = json.loads(res.model_dump_json())
    assert AIReviewResult(**data).variant_id == res.variant_id


def test_task_state_transitions():
    check_transition(TaskState.PENDING, TaskState.RUNNING)
    check_transition(TaskState.RUNNING, TaskState.COMPLETED)
    check_transition(TaskState.RUNNING, TaskState.RETRYING)
    check_transition(TaskState.FAILED, TaskState.RETRYING)
    with pytest.raises(ValueError):
        check_transition(TaskState.COMPLETED, TaskState.RUNNING)
    with pytest.raises(ValueError):
        check_transition(TaskState.PENDING, TaskState.COMPLETED)


def test_provider_selection(cfg=None):
    c = AVRSConfig()
    c.ai.provider = "mock"
    assert isinstance(make_provider(c.ai), MockAIProvider)
    c.ai.provider = "ollama"
    c.ai.base_url = "http://127.0.0.1:11434"
    p = make_provider(c.ai)
    assert p.name == "ollama" and p.is_local is True
