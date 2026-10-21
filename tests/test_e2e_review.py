"""End-to-end Mode-B integration test on the demo dataset:
VCF import -> evidence -> triage -> IGV (mock) -> AI review -> DB -> report.
"""
import json
from pathlib import Path

import pytest

from avrs.db.models import Variant, ReviewTask, AIReview, HumanReview, IGVSnapshot
from avrs.io.fasta import Fasta
from avrs.review import (import_sample_bam_vcf, triage_sample, run_reviews,
                         record_human_decision)
from avrs.reports import generate_reports, collect_rows
from avrs.variants import import_vcf


def test_mode_b_end_to_end(cfg, db_session, demo_dataset):
    sample_id = demo_dataset["sample_id"]
    # 1) register + import
    smp = import_sample_bam_vcf(db_session, cfg, sample_id,
                                demo_dataset["bam"], demo_dataset["vcf"])
    fasta = Fasta(demo_dataset["fasta"])
    n, skipped = import_vcf(db_session, demo_dataset["vcf"], smp.id, sample_id,
                            demo_dataset["reference_build"], fasta=fasta)
    db_session.commit()
    assert n == 7, f"expected 7 VCF variants (one scenario is intentionally absent), got {n}"

    # 2) evidence + triage
    tasks = triage_sample(db_session, cfg, sample_id, demo_dataset["fasta"],
                          demo_dataset["bam"])
    assert tasks > 0
    evs = db_session.query(Variant).count()
    assert evs == 7

    # 3) IGV (mock) + AI review
    stats = run_reviews(db_session, cfg, sample_id, demo_dataset["bam"],
                        demo_dataset["vcf"], demo_dataset["fasta"], use_mock_igv=True)
    assert stats["attempted"] == tasks
    assert stats["completed"] == tasks
    assert stats["failed"] == 0
    # mock IGV produced real PNG files bound to the right variant
    snaps = db_session.query(IGVSnapshot).all()
    assert len(snaps) == 2 * tasks
    for s in snaps:
        assert Path(s.path).exists()
        assert s.variant_uid and s.sample_id == smp.id
        assert s.bam_checksum and s.vcf_checksum
        assert s.locus
    # no cross-variant mixing: each variant's snapshot path contains its uid
    for s in snaps:
        assert s.variant_uid.replace(":", "_").replace(">", "_") in s.path or \
            s.variant_uid.split("#")[-1].split(":")[0] in s.path

    # 4) AI reviews persisted with provenance
    ais = db_session.query(AIReview).all()
    assert len(ais) == tasks
    for a in ais:
        assert a.model_name.startswith("mock")
        assert a.prompt_version
        assert a.ai_recommendation in (
            "LIKELY_TRUE_VARIANT", "SUSPECTED_FALSE_POSITIVE", "INCONCLUSIVE",
            "INSUFFICIENT_EVIDENCE", "REVIEW_FAILED", "POSSIBLE_FALSE_NEGATIVE")
        res = a.result
        assert res["risk_score_version"]
        assert res["risk_score_is_calibrated_probability"] is False
        assert res["visual_review_performed"] is False

    # 5) human decision with history
    v = db_session.query(Variant).first()
    r1 = record_human_decision(db_session, v.variant_uid, sample_id, "alice",
                               "reject_variant", reason="artifact", comment="first")
    r2 = record_human_decision(db_session, v.variant_uid, sample_id, "bob",
                               "accept_variant", reason="confirmed")
    db_session.refresh(r1)
    assert r1.superseded_by_id == r2.id
    assert r2.supersedes_id == r1.id
    assert db_session.query(HumanReview).filter_by(variant_id=v.id).count() == 2

    # 6) reports
    out = Path(cfg.paths.output_dir) / "reports"
    files = generate_reports(db_session, out, sample_id=sample_id)
    for key in ("review_results_json", "review_results_csv", "review_summary_json",
                "review_summary_html", "exception_report", "audit_log", "qc_summary"):
        assert Path(files[key]).exists()
    summary = json.loads(Path(files["review_summary_json"]).read_text())
    assert summary["total_variants"] == 7
    assert summary["human_decided"] >= 1
    html = Path(files["review_summary_html"]).read_text()
    assert "Research / validation use only" in html
    assert sample_id in html


def test_review_task_states_progress(cfg, db_session, demo_dataset):
    sample_id = demo_dataset["sample_id"]
    smp = import_sample_bam_vcf(db_session, cfg, sample_id, demo_dataset["bam"],
                                demo_dataset["vcf"])
    triage_sample(db_session, cfg, sample_id, demo_dataset["fasta"], demo_dataset["bam"])
    run_reviews(db_session, cfg, sample_id, demo_dataset["bam"], demo_dataset["vcf"],
                demo_dataset["fasta"], use_mock_igv=True)
    tasks = db_session.query(ReviewTask).all()
    assert all(t.status.value == "COMPLETED" for t in tasks)
    assert all(t.attempts == 1 for t in tasks)


def test_report_rows_separate_ai_and_human(cfg, db_session, demo_dataset):
    sample_id = demo_dataset["sample_id"]
    smp = import_sample_bam_vcf(db_session, cfg, sample_id, demo_dataset["bam"],
                                demo_dataset["vcf"])
    fasta = Fasta(demo_dataset["fasta"])
    import_vcf(db_session, demo_dataset["vcf"], smp.id, sample_id,
               demo_dataset["reference_build"], fasta=fasta)
    db_session.commit()
    triage_sample(db_session, cfg, sample_id, demo_dataset["fasta"], demo_dataset["bam"])
    run_reviews(db_session, cfg, sample_id, demo_dataset["bam"], demo_dataset["vcf"],
                demo_dataset["fasta"], use_mock_igv=True)
    rows = collect_rows(db_session, sample_id)
    assert rows
    for r in rows:
        assert "ai_recommendation" in r and "human_decision" in r
        # AI label must not be presented as a truth label
        assert r["ai_recommendation"] != "TP" and r["ai_recommendation"] != "FP"
