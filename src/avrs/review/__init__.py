"""Review orchestration service: the Mode-B end-to-end loop.

VCF import -> evidence extraction -> triage -> review-task creation ->
(IGV snapshots -> AI review) -> persistence. Both the CLI and the web API
call THIS layer (spec 13: CLI and UI share one business core).

Task lifecycle follows TaskState with validated transitions; one failed
variant never aborts the whole sample batch (spec 15.2).
"""
from __future__ import annotations

import json
import logging
import time
from pathlib import Path
from typing import Optional

from sqlalchemy.orm import Session

from ..db import models as db
from ..db import audit
from ..config import AVRSConfig
from ..evidence import make_backend, EvidenceParams
from ..evidence.metrics import VariantEvidence
from ..triage import evaluate_rules, prioritize, should_review, suspiciousness_score
from ..variants import import_vcf
from ..ai import make_provider, ReviewEngine
from ..igv.client import IGVMCPClient
from ..igv.transport import make_transport
from ..igv.reviewer import IGVReviewer, VariantReviewTask, checksum_if_exists
from ..igv.mock_server import start_mock_igv

log = logging.getLogger("avrs.review")


def import_sample_bam_vcf(session: Session, cfg: AVRSConfig, sample_id: str,
                          bam_path: str | Path, vcf_path: str | Path,
                          sample_type: str = "germline_dna") -> db.Sample:
    """Register a sample with its BAM/VCF (Mode B) or full-pipeline outputs."""
    sample = session.query(db.Sample).filter_by(sample_id=sample_id).one_or_none()
    if sample:
        return sample
    sample = db.Sample(
        sample_id=sample_id, sample_type=sample_type,
        analysis_type=cfg.analysis.data_type,
        reference_build=cfg.reference.build, status="IMPORTED",
    )
    session.add(sample)
    session.flush()
    for kind, path in [("bam", bam_path), ("bam_index", str(bam_path) + ".bai"),
                       ("vcf", vcf_path), ("reference", cfg.reference.fasta)]:
        p = Path(path) if path else None
        rec = db.SampleFile(sample_id=sample.id, kind=kind, path=str(p) if p else "",
                            checksum=checksum_if_exists(p) if p and p.exists() else None)
        session.add(rec)
    session.commit()
    audit(session, actor="system:review.service", action="register_sample",
          entity_type="sample", entity_id=sample_id,
          details={"bam": str(bam_path), "vcf": str(vcf_path)})
    return sample


def triage_sample(session: Session, cfg: AVRSConfig, sample_id: str,
                  fasta_path: str | Path, bam_path: str | Path) -> int:
    """Compute evidence + rules + priority for every imported variant,
    create ReviewTasks per configured review modes."""
    sample = session.query(db.Sample).filter_by(sample_id=sample_id).one()
    profile = cfg.triage.profile_for(cfg.analysis.data_type)
    params = EvidenceParams(min_mapq=profile.min_mapq, min_baseq=profile.min_baseq,
                            min_dist_from_read_end=profile.min_dist_from_read_end)
    backend = make_backend(bam_path, fasta_path, params)
    bed_intervals = _load_bed(cfg.analysis.target_bed) if cfg.analysis.target_bed else []

    def _region_for(chrom: str, pos: int):
        for c, s, e in bed_intervals:
            if c == chrom and s <= pos <= e:
                return (s, e)
        return None

    n_tasks = 0
    for variant in session.query(db.Variant).filter_by(sample_id=sample.id).all():
        ev_rec = session.query(db.VariantEvidence).filter_by(
            variant_id=variant.id).one_or_none()
        if ev_rec is None:
            info = variant.info
            caller = {
                "DP": info.get("DP"), "AD_REF": info.get("AD_REF"),
                "AD_ALT": info.get("AD_ALT"), "GQ": info.get("GQ"),
                "QUAL": variant.qual, "GT": variant.genotype,
            }
            try:
                ev = backend.evidence(
                    variant.chrom, variant.pos, variant.ref, variant.alt,
                    variant.variant_type, caller=caller,
                    target_region=_region_for(variant.chrom, variant.pos))
            except Exception as e:  # noqa: BLE001 - record, don't crash the batch
                log.warning("evidence extraction failed for %s: %s", variant.variant_uid, e)
                ev = VariantEvidence(
                    chrom=variant.chrom, pos=variant.pos, ref=variant.ref,
                    alt=variant.alt, variant_type=variant.variant_type,
                    backend=backend.name)
                ev.warnings.append(f"evidence_extraction_error: {e}")
            hits = evaluate_rules(ev, profile)
            prio = prioritize(
                hits, in_high_priority_region=False,
                vcf_filter_failed=(variant.original_filter not in (None, ".", "PASS")))
            ev_rec = db.VariantEvidence(
                variant_id=variant.id, sample_id=sample.id, backend=backend.name,
                metrics_json=ev.model_dump_json(),
                rule_hits_json=json.dumps([h.model_dump() for h in hits]),
                priority=prio,
            )
            session.add(ev_rec)
        else:
            hits = _hits_from_json(ev_rec.rule_hits)
            prio = ev_rec.priority or "MEDIUM"
        # create the review task if it does not exist yet (idempotent re-runs)
        has_task = session.query(db.ReviewTask).filter_by(
            sample_id=sample.id, variant_id=variant.id).first() is not None
        if not has_task and should_review(
                cfg.review.modes[0], hits, prio, variant.gene, cfg.review.target_genes,
                variant.original_filter not in (None, ".", "PASS"),
                bool(ev_rec.metrics.get("in_target_region")), variant.variant_uid, []):
            task = db.ReviewTask(
                sample_id=sample.id, variant_id=variant.id,
                review_mode=cfg.review.modes[0], priority=prio,
            )
            session.add(task)
            n_tasks += 1
    session.commit()
    audit(session, actor="system:review.service", action="triage_sample",
          entity_type="sample", entity_id=sample_id, details={"review_tasks": n_tasks})
    return n_tasks


def _load_bed(bed_path) -> list[tuple[str, int, int]]:
    """Load BED intervals as (chrom, start1, end). BED is 0-based half-open."""
    out = []
    try:
        for line in Path(bed_path).read_text(encoding="utf-8").splitlines():
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            parts = line.split("\t")
            c, s, e = parts[0], int(parts[1]) + 1, int(parts[2])
            out.append((c, s, e))
    except Exception as e:
        log.warning("failed to parse BED %s: %s", bed_path, e)
    return out


def run_reviews(session: Session, cfg: AVRSConfig, sample_id: str,
                bam_path: str | Path, vcf_path: str | Path,
                fasta_path: str | Path,
                use_mock_igv: bool = False, limit: Optional[int] = None) -> dict:
    """Execute pending review tasks end to end. Returns summary stats."""
    sample = session.query(db.Sample).filter_by(sample_id=sample_id).one()
    tasks = (session.query(db.ReviewTask)
             .filter_by(sample_id=sample.id)
             .filter(db.ReviewTask.status.in_([db.TaskState.PENDING, db.TaskState.RETRYING]))
             .order_by(db.ReviewTask.priority.desc())
             .all())
    if limit:
        tasks = tasks[:limit]
    profile = cfg.triage.profile_for(cfg.analysis.data_type)
    backend = make_backend(bam_path, fasta_path,
                           EvidenceParams(min_mapq=profile.min_mapq,
                                          min_baseq=profile.min_baseq))
    provider = make_provider(cfg.ai)
    engine = ReviewEngine(provider, cfg, thresholds=profile.model_dump())

    mock_server = None
    if use_mock_igv:
        mock_server = start_mock_igv()
        from ..igv.transport import DirectTCPTransport
        transport = DirectTCPTransport("127.0.0.1", mock_server.port, timeout=10.0)
        client = IGVMCPClient(transport)
        reviewer = IGVReviewer(client, cfg, cfg.resolve_output(cfg.paths.snapshot_dir))
        igv_transport_label = "mock"
    else:
        transport = make_transport(cfg.igv)
        client = IGVMCPClient(transport, cfg.igv.wait.load_wait_s, cfg.igv.wait.goto_wait_s)
        reviewer = IGVReviewer(client, cfg, cfg.resolve_output(cfg.paths.snapshot_dir))
        igv_transport_label = cfg.igv.transport

    bam_sum = checksum_if_exists(bam_path)
    vcf_sum = checksum_if_exists(vcf_path)
    genome_id = cfg.reference.fasta or cfg.reference.build

    stats = {"attempted": 0, "completed": 0, "failed": 0, "snapshot_failures": 0,
             "visual_reviews": 0, "durations": [], "igv_transport": igv_transport_label}
    igv_session = db.IGVSession(session_token=client.session_token,
                                transport=igv_transport_label,
                                host=cfg.igv.host, port=cfg.igv.port)
    session.add(igv_session)
    session.flush()

    for task in tasks:
        stats["attempted"] += 1
        t0 = time.time()
        variant = session.get(db.Variant, task.variant_id)
        ev_rec = (session.query(db.VariantEvidence)
                  .filter_by(variant_id=variant.id).one_or_none())
        if not ev_rec:
            task.set_state(db.TaskState.FAILED)
            task.last_error = "no evidence record; run triage first"
            session.commit()
            stats["failed"] += 1
            continue
        task.set_state(db.TaskState.RUNNING)
        task.attempts += 1
        session.commit()
        try:
            igv_task = VariantReviewTask(
                sample_id=sample_id, variant_uid=variant.variant_uid,
                chrom=variant.chrom, pos=variant.pos, ref=variant.ref, alt=variant.alt,
                variant_type=variant.variant_type, genome_build=cfg.reference.build,
                bam_path=str(bam_path),
                bam_index_path=str(bam_path) + ".bai",
                vcf_path=str(vcf_path),
                genome_id=genome_id, bam_checksum=bam_sum, vcf_checksum=vcf_sum,
            )
            igv_result = reviewer.review_variant(igv_task, retries=cfg.review.max_retries - 1)
            good_snaps = [s for s in igv_result.snapshots if s.status == "CAPTURED"]
            bad_snaps = [s for s in igv_result.snapshots if s.status == "FAILED"]
            stats["snapshot_failures"] += len(bad_snaps)

            for snap in igv_result.snapshots:
                rec = db.IGVSnapshot(
                    review_task_id=task.id, session_id=igv_session.id,
                    sample_id=sample.id, variant_id=variant.id,
                    variant_uid=variant.variant_uid, kind=snap.kind, path=snap.path,
                    genome_build=snap.genome_build, locus=snap.locus,
                    bam_checksum=snap.bam_checksum, vcf_checksum=snap.vcf_checksum,
                    status=snap.status,
                )
                session.add(rec)

            ev = VariantEvidence(**ev_rec.metrics)
            hits = _hits_from_json(ev_rec.rule_hits)
            result = engine.review(
                ev, _hits_from_json(ev_rec.rule_hits),
                task_ctx={
                    "sample_id": sample_id, "variant_id": variant.variant_uid,
                    "chrom": variant.chrom, "pos": variant.pos, "ref": variant.ref,
                    "alt": variant.alt, "variant_type": variant.variant_type,
                    "vcf_filter": variant.original_filter, "genotype": variant.genotype,
                    "reference_build": cfg.reference.build,
                },
                snapshot_paths=[s.path for s in good_snaps],
            )
            session.add(db.AIReview(
                review_task_id=task.id, variant_id=variant.id,
                variant_uid=variant.variant_uid,
                result_json=result.model_dump_json(),
                model_name=result.model_name,
                prompt_version=result.prompt_version,
                visual_review_performed=result.visual_review_performed,
                requires_human_review=result.requires_human_review,
                ai_recommendation=result.ai_recommendation.value,
                risk_score=result.risk_score,
            ))
            task.set_state(db.TaskState.COMPLETED)
            task.last_error = igv_result.error
            stats["completed"] += 1
            stats["visual_reviews"] += 1 if result.visual_review_performed else 0
        except Exception as e:  # noqa: BLE001
            log.exception("review failed for %s", variant.variant_uid)
            if task.attempts < cfg.review.max_retries:
                task.set_state(db.TaskState.RETRYING)
            else:
                task.set_state(db.TaskState.FAILED)
            task.last_error = str(e)[:2000]
            stats["failed"] += 1
        stats["durations"].append(round(time.time() - t0, 2))
        session.commit()

    igv_session.ended_at = db.utcnow()
    session.commit()
    audit(session, actor="system:review.service", action="run_reviews",
          entity_type="sample", entity_id=sample_id, details=stats)
    if mock_server:
        mock_server.stop()
    return stats


def _hits_from_json(rule_hits_json: list[dict]):
    from ..triage import RuleHit
    return [RuleHit(**h) for h in rule_hits_json]


def record_human_decision(session: Session, variant_uid: str, sample_id: str,
                          reviewer: str, decision: str, reason: str | None = None,
                          comment: str | None = None) -> db.HumanReview:
    """Human decision with full history - never overwrites a prior decision."""
    if decision not in ("accept_variant", "reject_variant",
                        "needs_further_review", "insufficient_evidence"):
        raise ValueError(f"Unknown decision: {decision}")
    sample = session.query(db.Sample).filter_by(sample_id=sample_id).one()
    variant = (session.query(db.Variant)
               .filter_by(sample_id=sample.id, variant_uid=variant_uid).one())
    prev = (session.query(db.HumanReview)
            .filter_by(variant_id=variant.id)
            .order_by(db.HumanReview.id.desc()).first())
    rec = db.HumanReview(variant_id=variant.id, reviewer=reviewer,
                         decision=decision, reason=reason, comment=comment,
                         supersedes_id=prev.id if prev else None)
    session.add(rec)
    session.flush()
    if prev and prev.superseded_by_id is None:
        prev.superseded_by_id = rec.id
    session.commit()
    audit(session, actor=f"user:{reviewer}", action="human_decision",
          entity_type="variant", entity_id=variant_uid,
          details={"decision": decision, "reason": reason, "comment": comment})
    return rec
