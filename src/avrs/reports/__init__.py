"""Report generation (spec 16).

Produces review_results.json / .csv, review_summary.html / .json, a QC
summary, an exception report and the audit log. Reports clearly separate:
raw variant call | AI recommendation | human decision | truth-validated
result (never conflated).
"""
from __future__ import annotations

import csv
import datetime as _dt
import html
import json
import logging
from pathlib import Path
from typing import Any, Optional

from sqlalchemy.orm import Session

from ..db import models as db
from ..db import audit
from .. import DISCLAIMER

log = logging.getLogger("avrs.report")

DECISION_LABEL = {
    "accept_variant": "Accept variant",
    "reject_variant": "Reject variant",
    "needs_further_review": "Need further review",
    "insufficient_evidence": "Insufficient evidence",
}


def collect_rows(session: Session, sample_id: Optional[str] = None) -> list[dict[str, Any]]:
    """One row per variant: raw call + evidence + AI + human decision."""
    q = session.query(db.Variant)
    if sample_id:
        sample = session.query(db.Sample).filter_by(sample_id=sample_id).one()
        q = q.filter_by(sample_id=sample.id)
    rows = []
    for v in q.order_by(db.Variant.id).all():
        ev = (session.query(db.VariantEvidence).filter_by(variant_id=v.id).one_or_none())
        ai = (session.query(db.AIReview)
              .filter_by(variant_id=v.id).order_by(db.AIReview.id.desc()).first())
        humans = (session.query(db.HumanReview)
                  .filter_by(variant_id=v.id).order_by(db.HumanReview.id.desc()).all())
        human = humans[0] if humans else None
        snaps = session.query(db.IGVSnapshot).filter_by(variant_id=v.id).all()
        metrics = ev.metrics if ev else {}
        rule_hits = ev.rule_hits if ev else []
        rows.append({
            "sample_id": v.sample.sample_id,
            "variant_uid": v.variant_uid,
            "gene": v.gene or "",
            "chrom": v.chrom,
            "pos": v.pos,
            "ref": v.ref,
            "alt": v.alt,
            "variant_type": v.variant_type,
            "genotype": v.genotype or "",
            "vcf_filter": v.original_filter,
            "vcf_qual": v.qual,
            "caller_dp": metrics.get("caller_dp"),
            "caller_ad_ref": metrics.get("caller_ad_ref"),
            "caller_ad_alt": metrics.get("caller_ad_alt"),
            "caller_gq": metrics.get("caller_gq"),
            "bam_ref_count": metrics.get("ref_count"),
            "bam_alt_count": metrics.get("alt_count"),
            "bam_vaf": metrics.get("vaf"),
            "mean_mapq_alt": metrics.get("mean_mapq_alt"),
            "mean_baseq_alt": metrics.get("mean_baseq_alt"),
            "strand_imbalance": metrics.get("strand_ratio_imbalance"),
            "ad_consistency": metrics.get("ad_consistency"),
            "triage_priority": ev.priority if ev else "",
            "risk_rule_hits": ";".join(h.get("rule_id", "") for h in rule_hits),
            "ai_recommendation": ai.ai_recommendation if ai else "",
            "ai_risk_score": ai.risk_score if ai else None,
            "ai_visual_review_performed": ai.visual_review_performed if ai else None,
            "ai_model": ai.model_name if ai else "",
            "ai_prompt_version": ai.prompt_version if ai else "",
            "ai_reviewed_at": ai.created_at.isoformat() if ai else "",
            "human_decision": human.decision if human else "",
            "human_reviewer": human.reviewer if human else "",
            "human_decided_at": human.created_at.isoformat() if human else "",
            "human_comment": human.comment if human else "",
            "human_decision_count": len(humans),
            "snapshot_paths": ";".join(s.path for s in snaps),
            "snapshot_count": len(snaps),
        })
    return rows


def _software_versions(session: Session) -> dict:
    versions = {}
    for run in session.query(db.PipelineRun).all():
        if run.software_versions_json:
            versions.update(json.loads(run.software_versions_json))
    return versions


def generate_reports(session: Session, out_dir: str | Path, sample_id: Optional[str] = None,
                     batch_label: Optional[str] = None) -> dict[str, str]:
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    rows = collect_rows(session, sample_id)
    stamp = _dt.datetime.now().strftime("%Y%m%d_%H%M%S")
    prefix = f"{sample_id}_{stamp}" if sample_id else f"batch_{stamp}"

    # ---- review_results.json
    rj = out / f"review_results_{prefix}.json"
    rj.write_text(json.dumps({
        "generated_at": _dt.datetime.now().isoformat(),
        "sample_id": sample_id,
        "batch_label": batch_label,
        "disclaimer": DISCLAIMER,
        "software_versions": _software_versions(session),
        "records": rows,
    }, ensure_ascii=False, indent=2, default=str), encoding="utf-8")

    # ---- review_results.csv
    rc = out / f"review_results_{prefix}.csv"
    if rows:
        with open(rc, "w", newline="", encoding="utf-8-sig") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
            w.writeheader()
            w.writerows(rows)
    else:
        rc.write_text("sample_id,variant_uid\n", encoding="utf-8")

    # ---- review_summary.json
    summary = build_summary(session, rows, sample_id)
    sj = out / f"review_summary_{prefix}.json"
    sj.write_text(json.dumps(summary, ensure_ascii=False, indent=2, default=str), encoding="utf-8")

    # ---- review_summary.html
    sh = out / f"review_summary_{prefix}.html"
    sh.write_text(render_html(summary, rows), encoding="utf-8")

    # ---- exception report (failures, retries, undecided)
    exceptions = [r for r in rows
                  if (r["ai_recommendation"] in ("SUSPECTED_FALSE_POSITIVE", "INCONCLUSIVE",
                                                "INSUFFICIENT_EVIDENCE", "")
                      or r["snapshot_count"] == 0)]
    ex = out / f"exception_report_{prefix}.json"
    ex.write_text(json.dumps({"sample_id": sample_id, "exceptions": exceptions},
                             ensure_ascii=False, indent=2, default=str), encoding="utf-8")

    # ---- audit log
    al = out / f"audit_log_{prefix}.json"
    logs = [{"ts": a.ts.isoformat(), "actor": a.actor, "action": a.action,
             "entity_type": a.entity_type, "entity_id": a.entity_id,
             "details": json.loads(a.details_json) if a.details_json else None}
            for a in session.query(db.AuditLog).order_by(db.AuditLog.id).all()]
    al.write_text(json.dumps(logs, ensure_ascii=False, indent=2, default=str), encoding="utf-8")

    # ---- QC summary
    qc = out / f"qc_summary_{prefix}.json"
    qc.write_text(json.dumps(build_qc_summary(session, sample_id), ensure_ascii=False,
                             indent=2, default=str), encoding="utf-8")

    audit(session, actor="system:reports", action="export_report",
          entity_type="sample", entity_id=sample_id,
          details={"files": [str(rj), str(rc), str(sj), str(sh), str(ex), str(al), str(qc)]})
    return {
        "review_results_json": str(rj), "review_results_csv": str(rc),
        "review_summary_json": str(sj), "review_summary_html": str(sh),
        "exception_report": str(ex), "audit_log": str(al), "qc_summary": str(qc),
    }


def build_summary(session: Session, rows: list[dict], sample_id: Optional[str] | None = None) -> dict:
    ai_counts: dict[str, int] = {}
    human_counts: dict[str, int] = {}
    for r in rows:
        ai_counts[r["ai_recommendation"] or "not_reviewed"] = \
            ai_counts.get(r["ai_recommendation"] or "not_reviewed", 0) + 1
        if r["human_decision"]:
            human_counts[r["human_decision"]] = human_counts.get(r["human_decision"], 0) + 1
    reviewed = sum(1 for r in rows if r["ai_recommendation"])
    decided = sum(1 for r in rows if r["human_decision"])
    escalated = sum(1 for r in rows if r["ai_recommendation"] in
                    ("SUSPECTED_FALSE_POSITIVE", "INCONCLUSIVE", "INSUFFICIENT_EVIDENCE"))
    snap_fail = sum(1 for r in rows if r["snapshot_count"] == 0 and r["ai_recommendation"])
    durations = [r for r in rows if r.get("auto_review_seconds")]
    fn = []
    if sample_id:
        from ..fn_discovery import summarize as fn_sum
        fn = [fn_sum(session, sample_id)]
    return {
        "generated_at": _dt.datetime.now().isoformat(),
        "sample_id": sample_id,
        "disclaimer": DISCLAIMER,
        "total_variants": len(rows),
        "ai_reviewed": reviewed,
        "ai_escalated_to_human": escalated,
        "human_decided": decided,
        "pending_human_review": reviewed - decided,
        "suspected_false_positives": ai_counts.get("SUSPECTED_FALSE_POSITIVE", 0),
        "inconclusive": ai_counts.get("INCONCLUSIVE", 0),
        "insufficient_evidence": ai_counts.get("INSUFFICIENT_EVIDENCE", 0),
        "likely_true": ai_counts.get("LIKELY_TRUE_VARIANT", 0),
        "not_reviewed": ai_counts.get("not_reviewed", 0),
        "ai_recommendation_counts": ai_counts,
        "human_decision_counts": human_counts,
        "snapshot_failures": snap_fail,
        "visual_reviews_performed": sum(1 for r in rows if r.get("ai_visual_review_performed")),
        "fn_discovery": fn,
        "time_savings_note": ("Time-saved metrics require a measured human baseline; "
                              "none is claimed without one."),
    }


def build_qc_summary(session: Session, sample_id: Optional[str] = None) -> dict:
    q = session.query(db.QCMetric)
    if sample_id:
        sample = session.query(db.Sample).filter_by(sample_id=sample_id).one()
        q = q.filter_by(sample_id=sample.id)
    out: dict[str, dict] = {}
    for m in q.all():
        out.setdefault(m.source, {})[m.key] = json.loads(m.value_json)
    return {"sample_id": sample_id, "sources": out}


def render_html(summary: dict, rows: list[dict]) -> str:
    def esc(x) -> str:
        return html.escape(str(x if x is not None else ""))

    kpis = [
        ("Total variants", summary["total_variants"]),
        ("AI reviewed", summary["ai_reviewed"]),
        ("Escalated to human", summary["ai_escalated_to_human"]),
        ("Human decided", summary["human_decided"]),
        ("Pending review", summary["pending_human_review"]),
        ("Suspected false positives", summary["suspected_false_positives"]),
        ("Inconclusive", summary["inconclusive"]),
        ("Insufficient evidence", summary["insufficient_evidence"]),
        ("Snapshot failures", summary["snapshot_failures"]),
    ]
    kpi_html = "".join(
        f'<div class="kpi"><div class="kpi-v">{esc(v)}</div>'
        f'<div class="kpi-l">{esc(k)}</div></div>' for k, v in kpis)

    trs = []
    for r in rows:
        ai = r["ai_recommendation"] or "-"
        cls = {"SUSPECTED_FALSE_POSITIVE": "bad", "LIKELY_TRUE_VARIANT": "good",
               "INCONCLUSIVE": "warn", "INSUFFICIENT_EVIDENCE": "warn"}.get(ai, "")
        snaps = [p for p in (r.get("snapshot_paths") or "").split(";") if p]
        snap_html = " ".join(
            f'<a href="file:///{esc(p).lstrip("/")}" target="_blank">'
            f'{esc(p.split("/")[-1])}</a>' for p in snaps) or "&mdash;"
        trs.append(
            f'<tr><td>{esc(r["sample_id"])}</td>'
            f'<td class="mono">{esc(r["variant_uid"])}</td>'
            f'<td>{esc(r["chrom"])}:{esc(r["pos"])}</td>'
            f'<td class="mono">{esc(r["ref"])}&gt;{esc(r["alt"])}</td>'
            f'<td>{esc(r["genotype"])}</td>'
            f'<td>{esc(r["caller_dp"])}</td>'
            f'<td>{esc(r["bam_vaf"])}</td>'
            f'<td>{esc(r["caller_gq"])}</td>'
            f'<td>{esc(r["vcf_filter"])}</td>'
            f'<td class="{cls}">{esc(ai)}</td>'
            f'<td>{esc(r["ai_risk_score"])}</td>'
            f'<td>{esc(r["triage_priority"])}</td>'
            f'<td>{esc(DECISION_LABEL.get(r["human_decision"], r["human_decision"] or "-"))}</td>'
            f'<td>{esc(r["human_reviewer"])}</td>'
            f'<td>{snap_html}</td></tr>')

    return f"""<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8">
<title>AVRS review summary {esc(summary.get('sample_id') or 'batch')}</title>
<style>
 body{{font-family:-apple-system,Segoe UI,Roboto,Helvetica,Arial,sans-serif;margin:0;background:#f7f8fa;color:#1a1d23}}
 header{{background:#fff;border-bottom:1px solid #e4e7ec;padding:18px 28px}}
 h1{{font-size:19px;margin:0 0 4px}} .sub{{color:#667085;font-size:13px}}
 .warn{{background:#fff8e6;border:1px solid #f2d98c;color:#7a5b00;padding:9px 14px;border-radius:8px;font-size:12.5px;margin:14px 28px}}
 .kpis{{display:flex;flex-wrap:wrap;gap:12px;padding:18px 28px}}
 .kpi{{background:#fff;border:1px solid #e4e7ec;border-radius:10px;padding:12px 16px;min-width:132px}}
 .kpi-v{{font-size:22px;font-weight:600}} .kpi-l{{font-size:11.5px;color:#667085;text-transform:uppercase;letter-spacing:.03em}}
 table{{border-collapse:collapse;width:100%;background:#fff;margin:6px 28px 34px;font-size:12.5px}}
 th,td{{border-bottom:1px solid #eceef2;padding:7px 9px;text-align:left;vertical-align:top}}
 th{{background:#fafbfc;font-weight:600;color:#475467;position:sticky;top:0}}
 .mono{{font-family:ui-monospace,SFMono-Regular,Menlo,monospace;font-size:11.5px}}
 .good{{color:#067647;font-weight:600}} .bad{{color:#b42318;font-weight:600}} .warn{{color:#b54708;background:none;padding:0;border:none}}
 a{{color:#175cd3}}
</style></head><body>
<header><h1>AVRS Review Summary {esc(summary.get('sample_id') or '')}</h1>
<div class="sub">Generated {esc(summary['generated_at'])}</div></header>
<div class="warn"><b>Research / validation use only, until appropriately validated and approved.</b>
&nbsp;AI recommendations are suggestions, not truth labels. Human decisions take precedence.</div>
<div class="kpis">{kpi_html}</div>
<table><thead><tr>
<th>Sample</th><th>Variant ID</th><th>Locus</th><th>Alleles</th><th>GT</th>
<th>DP</th><th>VAF</th><th>GQ</th><th>FILTER</th><th>AI recommendation</th>
<th>Risk</th><th>Priority</th><th>Human decision</th><th>Reviewer</th><th>Snapshots</th>
</tr></thead><tbody>{''.join(trs)}</tbody></table>
</body></html>"""
