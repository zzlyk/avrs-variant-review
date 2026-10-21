"""SQLite persistence via SQLAlchemy 2.x. All 13 audit-oriented tables.

Design rules enforced here:
- Raw VCF/BAM files are recorded (checksums) but never mutated by AVRS.
- AI conclusions and human conclusions are stored separately; human decisions
  never silently overwrite prior decisions (decision history is kept via
  superseded_by chaining).
- Every candidate variant carries a stable, allele-level variant_uid that
  includes the reference build.
"""
from __future__ import annotations

import datetime as _dt
import enum
import hashlib
import json
from pathlib import Path
from typing import Any, Optional

from sqlalchemy import (
    DateTime, Enum, Float, ForeignKey, Index, Integer, String, Text, UniqueConstraint,
)
from sqlalchemy.orm import DeclarativeBase, Mapped, mapped_column, relationship


def utcnow() -> _dt.datetime:
    return _dt.datetime.now(_dt.timezone.utc).replace(tzinfo=None)


def file_checksum(path: str | Path, algo: str = "sha256") -> str:
    h = hashlib.new(algo)
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return f"{algo}:{h.hexdigest()}"


class Base(DeclarativeBase):
    pass


class TaskState(str, enum.Enum):
    PENDING = "PENDING"
    RUNNING = "RUNNING"
    COMPLETED = "COMPLETED"
    FAILED = "FAILED"
    RETRYING = "RETRYING"
    NEEDS_MANUAL_ACTION = "NEEDS_MANUAL_ACTION"


_VALID_TRANSITIONS = {
    TaskState.PENDING: {TaskState.RUNNING, TaskState.FAILED},
    TaskState.RUNNING: {TaskState.COMPLETED, TaskState.FAILED, TaskState.RETRYING, TaskState.NEEDS_MANUAL_ACTION},
    TaskState.RETRYING: {TaskState.RUNNING, TaskState.FAILED},
    TaskState.FAILED: {TaskState.RETRYING, TaskState.PENDING},
    TaskState.NEEDS_MANUAL_ACTION: {TaskState.PENDING, TaskState.RUNNING},
    TaskState.COMPLETED: set(),
}


def check_transition(src: TaskState, dst: TaskState):
    if dst not in _VALID_TRANSITIONS[src]:
        raise ValueError(f"Illegal task state transition {src.value} -> {dst.value}")


class Sample(Base):
    __tablename__ = "samples"
    id: Mapped[int] = mapped_column(primary_key=True)
    sample_id: Mapped[str] = mapped_column(String(255), unique=True, index=True)
    sample_type: Mapped[str] = mapped_column(String(64), default="germline_dna")
    analysis_type: Mapped[str] = mapped_column(String(32), default="panel")
    reference_build: Mapped[str] = mapped_column(String(64), default="GRCh38")
    status: Mapped[str] = mapped_column(String(32), default="IMPORTED")
    metadata_json: Mapped[str | None] = mapped_column(Text, nullable=True)
    created_at: Mapped[_dt.datetime] = mapped_column(DateTime, default=utcnow)

    files: Mapped[list["SampleFile"]] = relationship(back_populates="sample", cascade="all, delete-orphan")
    variants: Mapped[list["Variant"]] = relationship(back_populates="sample")

    @property
    def meta(self) -> dict:
        return json.loads(self.metadata_json) if self.metadata_json else {}


class SampleFile(Base):
    __tablename__ = "sample_files"
    id: Mapped[int] = mapped_column(primary_key=True)
    sample_id: Mapped[int] = mapped_column(ForeignKey("samples.id"), index=True)
    kind: Mapped[str] = mapped_column(String(32))  # fastq_r1|fastq_r2|bam|bam_index|vcf|reference|target_bed|gvcf
    path: Mapped[str] = mapped_column(Text)
    checksum: Mapped[str | None] = mapped_column(String(128), nullable=True)
    recorded_at: Mapped[_dt.datetime] = mapped_column(DateTime, default=utcnow)

    sample: Mapped[Sample] = relationship(back_populates="files")


class PipelineRun(Base):
    __tablename__ = "pipeline_runs"
    id: Mapped[int] = mapped_column(primary_key=True)
    sample_id: Mapped[str] = mapped_column(ForeignKey("samples.id"), index=True)
    mode: Mapped[str] = mapped_column(String(32), default="review_only")  # full_pipeline|review_only
    status: Mapped[TaskState] = mapped_column(Enum(TaskState), default=TaskState.PENDING)
    command: Mapped[str | None] = mapped_column(Text, nullable=True)
    log_path: Mapped[str | None] = mapped_column(Text, nullable=True)
    config_snapshot_json: Mapped[str | None] = mapped_column(Text, nullable=True)
    software_versions_json: Mapped[str | None] = mapped_column(Text, nullable=True)
    started_at: Mapped[_dt.datetime | None] = mapped_column(DateTime, nullable=True)
    finished_at: Mapped[_dt.datetime | None] = mapped_column(DateTime, nullable=True)
    created_at: Mapped[_dt.datetime] = mapped_column(DateTime, default=utcnow)


class QCMetric(Base):
    __tablename__ = "qc_metrics"
    id: Mapped[int] = mapped_column(primary_key=True)
    sample_id: Mapped[int] = mapped_column(ForeignKey("samples.id"), index=True)
    source: Mapped[str] = mapped_column(String(64))  # fastp|flagstat|samtools_stats|mosdepth|multiqc|evidence_qc
    key: Mapped[str] = mapped_column(String(128))
    value_json: Mapped[str] = mapped_column(Text)
    recorded_at: Mapped[_dt.datetime] = mapped_column(DateTime, default=utcnow)
    __table_args__ = (Index("ix_qc_sample_source_key", "sample_id", "source", "key"),)


class Variant(Base):
    __tablename__ = "variants"
    id: Mapped[int] = mapped_column(primary_key=True)
    variant_uid: Mapped[str] = mapped_column(String(255), index=True)  # build#chr:pos:ref>alt (normalized, allele-level)
    sample_id: Mapped[int] = mapped_column(ForeignKey("samples.id"), index=True)
    # original VCF representation (immutable record of what the caller reported)
    original_chrom: Mapped[str] = mapped_column(String(64))
    original_pos: Mapped[int] = mapped_column(Integer)
    original_ref: Mapped[str] = mapped_column(Text)
    original_alt: Mapped[str] = mapped_column(Text)
    original_filter: Mapped[str] = mapped_column(String(255), default=".")
    # normalized representation (allele level)
    chrom: Mapped[str] = mapped_column(String(64))
    pos: Mapped[int] = mapped_column(Integer)
    ref: Mapped[str] = mapped_column(Text)
    alt: Mapped[str] = mapped_column(Text)
    variant_type: Mapped[str] = mapped_column(String(16))  # SNV|INS|DEL|COMPLEX|MNV
    genotype: Mapped[str | None] = mapped_column(String(64), nullable=True)
    qual: Mapped[float | None] = mapped_column(Float, nullable=True)
    info_json: Mapped[str | None] = mapped_column(Text, nullable=True)  # DP,GQ,AD,QUAL... caller-reported
    gene: Mapped[str | None] = mapped_column(String(128), nullable=True)  # optional annotation
    annotation_json: Mapped[str | None] = mapped_column(Text, nullable=True)
    created_at: Mapped[_dt.datetime] = mapped_column(DateTime, default=utcnow)

    sample: Mapped[Sample] = relationship(back_populates="variants")
    __table_args__ = (UniqueConstraint("sample_id", "variant_uid", name="uq_sample_variant"),)

    @property
    def info(self) -> dict[str, Any]:
        return json.loads(self.info_json) if self.info_json else {}

    @property
    def locus(self) -> str:
        return f"{self.chrom}:{self.pos}"


class VariantEvidence(Base):
    __tablename__ = "variant_evidence"
    id: Mapped[int] = mapped_column(primary_key=True)
    variant_id: Mapped[int] = mapped_column(ForeignKey("variants.id"), index=True)
    sample_id: Mapped[int] = mapped_column(ForeignKey("samples.id"), index=True)
    backend: Mapped[str] = mapped_column(String(32))  # pysam|purepython
    metrics_json: Mapped[str] = mapped_column(Text)
    rule_hits_json: Mapped[str | None] = mapped_column(Text, nullable=True)
    priority: Mapped[str | None] = mapped_column(String(16), nullable=True)  # HIGH|MEDIUM|LOW
    computed_at: Mapped[_dt.datetime] = mapped_column(DateTime, default=utcnow)

    @property
    def metrics(self) -> dict[str, Any]:
        return json.loads(self.metrics_json)

    @property
    def rule_hits(self) -> list[dict[str, Any]]:
        return json.loads(self.rule_hits_json) if self.rule_hits_json else []


class ReviewTask(Base):
    __tablename__ = "review_tasks"
    id: Mapped[int] = mapped_column(primary_key=True)
    sample_id: Mapped[int] = mapped_column(ForeignKey("samples.id"), index=True)
    variant_id: Mapped[int] = mapped_column(ForeignKey("variants.id"), index=True)
    review_mode: Mapped[str] = mapped_column(String(32), default="suspicious_only")
    priority: Mapped[str] = mapped_column(String(16), default="MEDIUM")  # HIGH|MEDIUM|LOW
    status: Mapped[TaskState] = mapped_column(Enum(TaskState), default=TaskState.PENDING)
    attempts: Mapped[int] = mapped_column(Integer, default=0)
    last_error: Mapped[str | None] = mapped_column(Text, nullable=True)
    created_at: Mapped[_dt.datetime] = mapped_column(DateTime, default=utcnow)
    updated_at: Mapped[_dt.datetime] = mapped_column(DateTime, default=utcnow, onupdate=utcnow)

    snapshots: Mapped[list["IGVSnapshot"]] = relationship(back_populates="review_task")
    ai_reviews: Mapped[list["AIReview"]] = relationship(back_populates="review_task")

    def set_state(self, dst: TaskState):
        check_transition(self.status, dst)
        self.status = dst
        self.updated_at = utcnow()


class IGVSession(Base):
    __tablename__ = "igv_sessions"
    id: Mapped[int] = mapped_column(primary_key=True)
    session_token: Mapped[str] = mapped_column(String(64), unique=True)
    transport: Mapped[str] = mapped_column(String(16), default="direct")
    host: Mapped[str] = mapped_column(String(64), default="127.0.0.1")
    port: Mapped[int] = mapped_column(Integer, default=60151)
    igv_version: Mapped[str | None] = mapped_column(String(64), nullable=True)
    started_at: Mapped[_dt.datetime] = mapped_column(DateTime, default=utcnow)
    ended_at: Mapped[_dt.datetime | None] = mapped_column(DateTime, nullable=True)

    snapshots: Mapped[list["IGVSnapshot"]] = relationship(back_populates="session")


class IGVSnapshot(Base):
    __tablename__ = "igv_snapshots"
    id: Mapped[int] = mapped_column(primary_key=True)
    review_task_id: Mapped[int] = mapped_column(ForeignKey("review_tasks.id"), index=True)
    session_id: Mapped[int | None] = mapped_column(ForeignKey("igv_sessions.id"), nullable=True)
    sample_id: Mapped[int] = mapped_column(ForeignKey("samples.id"), index=True)
    variant_id: Mapped[int] = mapped_column(ForeignKey("variants.id"), index=True)
    variant_uid: Mapped[str] = mapped_column(String(255))
    kind: Mapped[str] = mapped_column(String(16))  # narrow|context
    path: Mapped[str] = mapped_column(Text)
    genome_build: Mapped[str] = mapped_column(String(64))
    locus: Mapped[str] = mapped_column(String(64))
    bam_checksum: Mapped[str | None] = mapped_column(String(128), nullable=True)
    vcf_checksum: Mapped[str | None] = mapped_column(String(128), nullable=True)
    igv_version: Mapped[str | None] = mapped_column(String(64), nullable=True)
    status: Mapped[str] = mapped_column(String(16), default="CAPTURED")  # CAPTURED|RETRIED|FAILED
    captured_at: Mapped[_dt.datetime] = mapped_column(DateTime, default=utcnow)

    review_task: Mapped[ReviewTask] = relationship(back_populates="snapshots")
    session: Mapped[IGVSession | None] = relationship(back_populates="snapshots")


class AIReview(Base):
    __tablename__ = "ai_reviews"
    id: Mapped[int] = mapped_column(primary_key=True)
    review_task_id: Mapped[int] = mapped_column(ForeignKey("review_tasks.id"), index=True)
    variant_id: Mapped[int] = mapped_column(ForeignKey("variants.id"), index=True)
    variant_uid: Mapped[str] = mapped_column(String(255))
    result_json: Mapped[str] = mapped_column(Text)  # full AIReviewResult
    model_name: Mapped[str] = mapped_column(String(128))
    model_version: Mapped[str] = mapped_column(String(64), default="")
    prompt_version: Mapped[str] = mapped_column(String(32), default="")
    visual_review_performed: Mapped[bool] = mapped_column(default=False)
    requires_human_review: Mapped[bool] = mapped_column(default=True)
    ai_recommendation: Mapped[str] = mapped_column(String(48))
    risk_score: Mapped[float | None] = mapped_column(Float, nullable=True)
    created_at: Mapped[_dt.datetime] = mapped_column(DateTime, default=utcnow)

    review_task: Mapped[ReviewTask] = relationship(back_populates="ai_reviews")

    @property
    def result(self) -> dict[str, Any]:
        return json.loads(self.result_json)


class HumanReview(Base):
    __tablename__ = "human_reviews"
    id: Mapped[int] = mapped_column(primary_key=True)
    variant_id: Mapped[int] = mapped_column(ForeignKey("variants.id"), index=True)
    review_task_id: Mapped[int | None] = mapped_column(ForeignKey("review_tasks.id"), nullable=True, index=True)
    reviewer: Mapped[str] = mapped_column(String(128))
    decision: Mapped[str] = mapped_column(String(32))
    # accept_variant | reject_variant | needs_further_review | insufficient_evidence
    reason: Mapped[str | None] = mapped_column(Text, nullable=True)
    comment: Mapped[str | None] = mapped_column(Text, nullable=True)
    supersedes_id: Mapped[int | None] = mapped_column(ForeignKey("human_reviews.id"), nullable=True)
    superseded_by_id: Mapped[int | None] = mapped_column(ForeignKey("human_reviews.id"), nullable=True)
    created_at: Mapped[_dt.datetime] = mapped_column(DateTime, default=utcnow)


class AuditLog(Base):
    __tablename__ = "audit_logs"
    id: Mapped[int] = mapped_column(primary_key=True)
    ts: Mapped[_dt.datetime] = mapped_column(DateTime, default=utcnow)
    actor: Mapped[str] = mapped_column(String(128))  # user:<name> | system:<module> | ai:<model>
    action: Mapped[str] = mapped_column(String(128))  # e.g. import_vcf, capture_snapshot, human_decision, export_report
    entity_type: Mapped[str | None] = mapped_column(String(64), nullable=True)
    entity_id: Mapped[str | None] = mapped_column(String(128), nullable=True)
    details_json: Mapped[str | None] = mapped_column(Text, nullable=True)


class ModelConfiguration(Base):
    __tablename__ = "model_configurations"
    id: Mapped[int] = mapped_column(primary_key=True)
    name: Mapped[str] = mapped_column(String(128), unique=True)
    provider: Mapped[str] = mapped_column(String(64))  # mock|openai_compatible|ollama
    model: Mapped[str] = mapped_column(String(128))
    supports_vision: Mapped[bool] = mapped_column(default=False)
    base_url: Mapped[str | None] = mapped_column(Text, nullable=True)
    config_json: Mapped[str | None] = mapped_column(Text, nullable=True)
    active: Mapped[bool] = mapped_column(default=False)
    created_at: Mapped[_dt.datetime] = mapped_column(DateTime, default=utcnow)


class FNDiscoveryCandidate(Base):
    __tablename__ = "fn_discovery_candidates"
    id: Mapped[int] = mapped_column(primary_key=True)
    sample_id: Mapped[int] = mapped_column(ForeignKey("samples.id"), index=True)
    candidate_id: Mapped[str] = mapped_column(String(255), unique=True)
    candidate_source: Mapped[str] = mapped_column(String(32))  # filtered_recovery|targeted_locus|second_caller
    locus: Mapped[str] = mapped_column(String(128))
    ref: Mapped[str | None] = mapped_column(Text, nullable=True)
    alt: Mapped[str | None] = mapped_column(Text, nullable=True)
    caller_status: Mapped[str] = mapped_column(String(64))  # in_vcf|filtered_in_raw|absent_from_vcf
    bam_evidence_json: Mapped[str | None] = mapped_column(Text, nullable=True)
    callable_status: Mapped[str] = mapped_column(String(48), default="unknown")
    possible_missed_variant: Mapped[bool] = mapped_column(default=False)
    classification: Mapped[str] = mapped_column(String(48), default="undetermined")
    # possible_missed_variant | no_call | insufficient_coverage | no_evidence_of_alt
    reason_codes: Mapped[str | None] = mapped_column(Text, nullable=True)
    human_review_required: Mapped[bool] = mapped_column(default=True)
    created_at: Mapped[_dt.datetime] = mapped_column(DateTime, default=utcnow)

    @property
    def bam_evidence(self) -> dict[str, Any]:
        return json.loads(self.bam_evidence_json) if self.bam_evidence_json else {}
