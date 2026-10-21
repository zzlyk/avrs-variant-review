"""FastAPI backend for the human review workbench (spec 8).

The API is a thin layer over avrs.review / avrs.reports - CLI and UI share
the same business logic. All state-changing endpoints write audit records.
"""
from __future__ import annotations

import os
from pathlib import Path
from typing import Any, Optional

from fastapi import FastAPI, HTTPException, Query
from fastapi.responses import FileResponse, HTMLResponse
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel
from sqlalchemy.orm import Session

from ..config import AVRSConfig, load_config
from ..db import models as db
from ..db import open_db, audit
from ..reports import collect_rows, build_summary, generate_reports
from ..review import record_human_decision

FRONTEND_DIR = Path(__file__).resolve().parents[3] / "frontend"


class DecisionIn(BaseModel):
    reviewer: str
    decision: str
    reason: Optional[str] = None
    comment: Optional[str] = None


def create_app(config_path: Optional[str | Path] = None) -> FastAPI:
    cfg: AVRSConfig = load_config(config_path) if config_path else AVRSConfig()
    engine, _ = open_db(cfg.resolve_output(cfg.paths.db))
    app = FastAPI(title="AVRS - AI Variant Review System", version="0.1.0")

    def session() -> Session:
        from ..db import get_session
        return get_session(engine)

    # ----------------------------------------------------------- samples
    @app.get("/api/samples")
    def list_samples():
        with session() as s:
            out = []
            for smp in s.query(db.Sample).order_by(db.Sample.id).all():
                n_var = s.query(db.Variant).filter_by(sample_id=smp.id).count()
                n_task = s.query(db.ReviewTask).filter_by(sample_id=smp.id).count()
                n_pending = s.query(db.ReviewTask).filter_by(
                    sample_id=smp.id, status=db.TaskState.PENDING).count()
                n_failed = s.query(db.ReviewTask).filter_by(
                    sample_id=smp.id, status=db.TaskState.FAILED).count()
                n_ai = s.query(db.AIReview).join(db.Variant).filter(
                    db.Variant.sample_id == smp.id).count()
                n_human = s.query(db.HumanReview).join(db.Variant).filter(
                    db.Variant.sample_id == smp.id).count()
                out.append({
                    "sample_id": smp.sample_id, "status": smp.status,
                    "reference_build": smp.reference_build,
                    "analysis_type": smp.analysis_type,
                    "variants": n_var, "review_tasks": n_task,
                    "pending_reviews": n_pending, "failed_tasks": n_failed,
                    "ai_reviews": n_ai, "human_reviews": n_human,
                    "created_at": smp.created_at.isoformat(),
                })
            return out

    # ----------------------------------------------------------- queue
    @app.get("/api/variants")
    def list_variants(sample_id: Optional[str] = None,
                      priority: Optional[str] = None,
                      ai: Optional[str] = None,
                      decision: Optional[str] = None,
                      q: Optional[str] = None,
                      limit: int = Query(200, le=5000),
                      offset: int = 0):
        rows = collect_rows(session(), sample_id)
        def keep(r):
            if priority and r["triage_priority"] != priority:
                return False
            if ai and r["ai_recommendation"] != ai:
                return False
            if decision == "pending" and r["human_decision"]:
                return False
            if decision and decision != "pending" and r["human_decision"] != decision:
                return False
            if q:
                blob = f"{r['variant_uid']} {r['chrom']}{r['pos']} {r['ref']}{r['alt']} {r['gene']}".lower()
                if q.lower() not in blob:
                    return False
            return True
        rows = [r for r in rows if keep(r)]
        return {"total": len(rows), "rows": rows[offset:offset + limit]}

    @app.get("/api/variants/{variant_uid:path}")
    def variant_detail(variant_uid: str, sample_id: Optional[str] = None):
        with session() as s:
            v = s.query(db.Variant).filter_by(variant_uid=variant_uid).one_or_none()
            if not v:
                raise HTTPException(404, f"variant not found: {variant_uid}")
            ev = s.query(db.VariantEvidence).filter_by(variant_id=v.id).one_or_none()
            ai = s.query(db.AIReview).filter_by(variant_id=v.id).order_by(
                db.AIReview.id.desc()).all()
            humans = s.query(db.HumanReview).filter_by(variant_id=v.id).order_by(
                db.HumanReview.id.desc()).all()
            snaps = s.query(db.IGVSnapshot).filter_by(variant_id=v.id).all()
            tasks = s.query(db.ReviewTask).filter_by(variant_id=v.id).all()
            return {
                "variant": {
                    "variant_uid": v.variant_uid, "sample_id": v.sample.sample_id,
                    "chrom": v.chrom, "pos": v.pos, "ref": v.ref, "alt": v.alt,
                    "variant_type": v.variant_type, "genotype": v.genotype,
                    "vcf_filter": v.original_filter, "qual": v.qual,
                    "original": {"chrom": v.original_chrom, "pos": v.original_pos,
                                 "ref": v.original_ref, "alt": v.original_alt},
                    "info": v.info,
                },
                "evidence": ev.metrics if ev else {},
                "rule_hits": ev.rule_hits if ev else [],
                "priority": ev.priority if ev else None,
                "ai_reviews": [a.result for a in ai],
                "human_reviews": [
                    {"id": h.id, "reviewer": h.reviewer, "decision": h.decision,
                     "reason": h.reason, "comment": h.comment,
                     "created_at": h.created_at.isoformat(),
                     "superseded_by": h.superseded_by_id} for h in humans],
                "snapshots": [
                    {"id": sn.id, "kind": sn.kind, "path": sn.path, "locus": sn.locus,
                     "status": sn.status, "genome_build": sn.genome_build,
                     "bam_checksum": sn.bam_checksum, "vcf_checksum": sn.vcf_checksum,
                     "captured_at": sn.captured_at.isoformat()} for sn in snaps],
                "review_tasks": [{"id": t.id, "status": t.status.value,
                                  "priority": t.priority, "attempts": t.attempts,
                                  "last_error": t.last_error} for t in tasks],
            }

    @app.post("/api/variants/{variant_uid:path}/decision")
    def post_decision(variant_uid: str, payload: DecisionIn):
        with session() as s:
            v = s.query(db.Variant).filter_by(variant_uid=variant_uid).one_or_none()
            if not v:
                raise HTTPException(404, f"variant not found: {variant_uid}")
            try:
                rec = record_human_decision(
                    s, variant_uid, v.sample.sample_id, payload.reviewer,
                    payload.decision, payload.reason, payload.comment)
            except ValueError as e:
                raise HTTPException(400, str(e))
            return {"id": rec.id, "decision": rec.decision,
                    "created_at": rec.created_at.isoformat()}

    @app.get("/api/snapshots/{snapshot_id}/image")
    def snapshot_image(snapshot_id: int):
        with session() as s:
            sn = s.query(db.IGVSnapshot).get(snapshot_id)
            if not sn:
                raise HTTPException(404, "snapshot not found")
            p = Path(sn.path)
            if not p.exists():
                raise HTTPException(404, f"snapshot file missing on disk: {p}")
            return FileResponse(p)

    # ----------------------------------------------------------- metrics
    @app.get("/api/metrics")
    def metrics(sample_id: Optional[str] = None):
        with session() as s:
            rows = collect_rows(s, sample_id)
            summary = build_summary(s, rows, sample_id)
            from ..fn_discovery import summarize as fn_sum
            if sample_id:
                summary["fn_discovery"] = [fn_sum(s, sample_id)]
            return summary

    @app.get("/api/fn-candidates")
    def fn_candidates(sample_id: str):
        with session() as s:
            smp = s.query(db.Sample).filter_by(sample_id=sample_id).one()
            cands = s.query(db.FNDiscoveryCandidate).filter_by(sample_id=smp.id).all()
            return [{
                "candidate_id": c.candidate_id, "locus": c.locus,
                "source": c.candidate_source, "caller_status": c.caller_status,
                "classification": c.classification,
                "possible_missed_variant": c.possible_missed_variant,
                "callable_status": c.callable_status,
                "bam_evidence": c.bam_evidence,
                "human_review_required": c.human_review_required,
            } for c in cands]

    @app.post("/api/reports")
    def make_report(sample_id: Optional[str] = None):
        with session() as s:
            files = generate_reports(s, cfg.resolve_output("reports"),
                                     sample_id=sample_id)
            return files

    @app.get("/api/audit")
    def audit_log(limit: int = 200):
        with session() as s:
            rows = s.query(db.AuditLog).order_by(db.AuditLog.id.desc()).limit(limit).all()
            return [{"ts": a.ts.isoformat(), "actor": a.actor, "action": a.action,
                     "entity_type": a.entity_type, "entity_id": a.entity_id,
                     "details": a.details_json} for a in rows]

    @app.get("/api/config")
    def config_view():
        return {
            "project_name": cfg.project_name,
            "reference_build": cfg.reference.build,
            "analysis": cfg.analysis.model_dump(),
            "ai": {k: v for k, v in cfg.ai.model_dump().items() if k != "api_key"},
            "igv": {"transport": cfg.igv.transport, "host": cfg.igv.host,
                    "port": cfg.igv.port},
            "triage_profiles": sorted(cfg.triage.profiles),
        }

    # ----------------------------------------------------------- UI
    if FRONTEND_DIR.exists():
        app.mount("/static", StaticFiles(directory=str(FRONTEND_DIR)), name="static")

        @app.get("/", response_class=HTMLResponse)
        def index():
            f = FRONTEND_DIR / "index.html"
            if f.exists():
                return HTMLResponse(f.read_text(encoding="utf-8"))
            return HTMLResponse("<h1>AVRS</h1><p>frontend not found</p>")

    return app


app = None  # populated by `avrs serve` (uvicorn factory)
