"""Automated single-variant IGV review executor (spec 5.2).

Per variant: lock -> health check -> (genome/BAM/VCF load once per sample
session) -> goto with type-aware window -> track configuration -> render
wait (file existence + stability + min size, never sleep-only) -> narrow +
context snapshots -> metadata binding -> release.

Snapshot metadata (sample_id, variant_id, genome_build, actual locus, BAM/
VCF checksums, timestamp, session id, path) is attached at capture time so
a screenshot can never be misattributed to the neighbouring variant.
"""
from __future__ import annotations

import re
import time
import uuid
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

from .client import IGVMCPClient, IGVInstanceLock
from ..db.models import file_checksum


def _safe_name(s: str) -> str:
    return re.sub(r"[^A-Za-z0-9._>:>-]", "_", s)[:150]


@dataclass
class SnapshotMeta:
    kind: str          # narrow | context
    path: str
    sample_id: str
    variant_uid: str
    genome_build: str
    locus: str
    bam_checksum: Optional[str]
    vcf_checksum: Optional[str]
    session_id: str
    timestamp: str
    igv_transport: str
    status: str = "CAPTURED"   # CAPTURED | FAILED


@dataclass
class VariantReviewTask:
    sample_id: str
    variant_uid: str
    chrom: str
    pos: int            # 1-based
    ref: str
    alt: str
    variant_type: str
    genome_build: str
    bam_path: str
    bam_index_path: Optional[str] = None
    vcf_path: Optional[str] = None
    genome_id: Optional[str] = None       # IGV genome id or fasta path
    bam_checksum: Optional[str] = None
    vcf_checksum: Optional[str] = None


@dataclass
class ReviewExecutionResult:
    task: VariantReviewTask
    session_id: str
    snapshots: list[SnapshotMeta] = field(default_factory=list)
    error: Optional[str] = None
    attempts: int = 0


class IGVReviewer:
    """Executes IGV automation for one variant at a time."""

    def __init__(self, client: IGVMCPClient, cfg, snapshot_root: str | Path,
                 lock_path: str | Path | None = None):
        self.client = client
        self.cfg = cfg  # AVRSConfig
        self.snapshot_root = Path(snapshot_root)
        self.lock_path = Path(lock_path) if lock_path else self.snapshot_root / ".igv.lock"
        # per-session data-load cache: keyed by (genome, bam, vcf)
        self._loaded_key: Optional[tuple] = None

    # ------------------------------------------------------------------
    def review_variant(self, task: VariantReviewTask, retries: int = 2) -> ReviewExecutionResult:
        rcfg = self.cfg.review
        result = ReviewExecutionResult(task=task, session_id=uuid.uuid4().hex[:16])
        lock = IGVInstanceLock(self.lock_path)
        try:
            lock.acquire(owner=f"avrs-review:{task.sample_id}")
        except Exception as e:
            result.error = f"IGV lock unavailable: {e}"
            return result
        try:
            if not self.client.health_check():
                result.error = "IGV Desktop command port is not reachable"
                return result
            self._ensure_sample_loaded(task)
            self._navigate(task)
            for kind in ("narrow", "context"):
                meta = self._capture(task, kind, result.session_id, retries)
                result.snapshots.append(meta)
            return result
        except Exception as e:  # noqa: BLE001 - convert to result error
            result.error = f"IGV review failed: {e}"
            return result
        finally:
            lock.release()

    # ------------------------------------------------------------------
    def _map_path(self, p: str | Path) -> str:
        """Path as seen by the IGV instance.

        Path mapping is applied only when explicitly enabled, because it is
        only correct when IGV runs on a different filesystem than the
        analysis (Windows IGV + WSL2 pipeline)."""
        if getattr(self.cfg.igv, "path_mapping_enabled", False):
            return self.cfg.map_path_for_igv(p)
        path = Path(p)
        return str(path.resolve()) if path.exists() else str(path)

    def _ensure_sample_loaded(self, task: VariantReviewTask):
        """Load genome/BAM/VCF once per (genome,bam,vcf) within this reviewer."""
        key = (task.genome_id, task.bam_path, task.vcf_path)
        if self._loaded_key == key:
            return
        if task.genome_id:
            self.client.load_genome(self._map_path(task.genome_id))
            time.sleep(self.cfg.igv.wait.load_wait_s)
        else:
            raise ValueError("genome_id/fasta required for IGV loading")
        self.client.load_alignment(self._map_path(task.bam_path),
                                   self._map_path(task.bam_index_path) if task.bam_index_path else None)
        time.sleep(self.cfg.igv.wait.load_wait_s)
        if task.vcf_path:
            self.client.load_vcf(self._map_path(task.vcf_path))
            time.sleep(self.cfg.igv.wait.load_wait_s)
        self.client.set_max_panel_height(self.cfg.igv.max_panel_height)
        self._loaded_key = key

    def _window_for(self, task: VariantReviewTask, kind: str) -> tuple[int, int]:
        rcfg = self.cfg.review
        if kind == "context":
            w = rcfg.window_width_context
        else:
            w = rcfg.window_width_narrow
            if task.variant_type in ("INS", "DEL", "COMPLEX"):
                w = max(w, len(task.ref) + len(task.alt) + 2 * rcfg.window_width_indel_pad)
        half = w // 2
        start = max(1, task.pos - half)
        end = task.pos + half
        return start, end

    def _navigate(self, task: VariantReviewTask):
        rcfg = self.cfg.review
        s, e = self._window_for(task, "narrow")
        # range view gives deterministic width control
        self.client.goto_range(task.chrom, s, e)
        time.sleep(rcfg.snapshot_wait_initial_s * 0.5)
        # read-friendly track configuration
        try:
            self.client.squish()
        except Exception:
            pass
        try:
            self.client.toggle_pair_view(True)
        except Exception:
            pass
        try:
            # sort by base at the variant to surface ALT-supporting reads
            self.client.sort_alignments("BASE", f"{task.chrom}:{task.pos}")
        except Exception:
            pass
        time.sleep(self.cfg.igv.wait.goto_wait_s)

    def _snapshot_path(self, task: VariantReviewTask, kind: str) -> Path:
        vname = _safe_name(task.variant_uid)
        return (self.snapshot_root / _safe_name(task.sample_id) / vname / f"{kind}.png")

    def _capture(self, task: VariantReviewTask, kind: str, session_id: str,
                 retries: int) -> SnapshotMeta:
        rcfg = self.cfg.review
        if kind == "context":
            s, e = self._window_for(task, "context")
            self.client.goto_range(task.chrom, s, e)
            time.sleep(self.cfg.igv.wait.goto_wait_s)
        path = self._snapshot_path(task, kind)
        path.parent.mkdir(parents=True, exist_ok=True)
        snap_dir = self._map_path(path.parent)
        self.client.set_snapshot_directory(snap_dir)
        fname = path.name
        error: Optional[str] = None
        status = "CAPTURED"
        for attempt in range(retries + 1):
            try:
                # clear stale file so verification cannot pass on an old image
                if path.exists():
                    path.unlink()
                self.client.capture_snapshot(fname)
                ok = self._wait_for_file(path, rcfg)
                if ok:
                    return SnapshotMeta(
                        kind=kind, path=str(path), sample_id=task.sample_id,
                        variant_uid=task.variant_uid, genome_build=task.genome_build,
                        locus=f"{task.chrom}:{self._window_for(task, kind)[0]}-{self._window_for(task, kind)[1]}",
                        bam_checksum=task.bam_checksum, vcf_checksum=task.vcf_checksum,
                        session_id=session_id, timestamp=time.strftime("%Y-%m-%dT%H:%M:%S%z"),
                        igv_transport=type(self.client.transport).__name__,
                    )
                error = "snapshot file missing/too small/unstable after wait window"
            except Exception as e:  # noqa: BLE001
                error = f"snapshot attempt {attempt + 1} failed: {e}"
        # failure path - explicit, never silent
        return SnapshotMeta(
            kind=kind, path=str(path), sample_id=task.sample_id,
            variant_uid=task.variant_uid, genome_build=task.genome_build,
            locus=f"{task.chrom}:{task.pos}", bam_checksum=task.bam_checksum,
            vcf_checksum=task.vcf_checksum, session_id=session_id,
            timestamp=time.strftime("%Y-%m-%dT%H:%M:%S%z"),
            igv_transport=type(self.client.transport).__name__,
            status="FAILED",
        )

    def _wait_for_file(self, path: Path, rcfg) -> bool:
        """Render-success detection: existence + min size + two-poll stability.
        Fixed sleep alone is never the success criterion."""
        deadline = time.time() + rcfg.snapshot_timeout_s
        last_size = -1
        stable_polls = 0
        while time.time() < deadline:
            if path.exists():
                size = path.stat().st_size
                if size >= rcfg.require_min_png_bytes:
                    if size == last_size:
                        stable_polls += 1
                        if stable_polls >= 2:
                            return True
                    else:
                        stable_polls = 0
                        last_size = size
            time.sleep(rcfg.snapshot_poll_interval_s)
        return False


def checksum_if_exists(p: str | Path | None) -> Optional[str]:
    if p is None:
        return None
    p = Path(p)
    if not p.exists():
        return None
    try:
        return file_checksum(p)
    except OSError:
        return None
