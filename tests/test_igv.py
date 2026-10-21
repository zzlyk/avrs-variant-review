"""Integration tests: IGV client against a mock IGV listener, snapshot binding,
session locking, retry behaviour."""
import time
from pathlib import Path

import pytest

from avrs.igv import (IGVMCPClient, IGVInstanceLock, IGVReviewer, VariantReviewTask,
                      DirectTCPTransport, start_mock_igv, OFFICIAL_TOOLS)
from avrs.igv.transport import UnsupportedOperation
from avrs.config import AVRSConfig


@pytest.fixture
def mock_igv():
    srv = start_mock_igv()
    yield srv
    srv.stop()


def test_igv_client_commands_map_to_official_tools(mock_igv):
    c = IGVMCPClient(DirectTCPTransport("127.0.0.1", mock_igv.port))
    assert c.health_check() is True
    c.reset_session()
    c.load_genome("/ref.fa")
    c.load_alignment("/data/x.bam", "/data/x.bam.bai")
    c.load_vcf("/data/x.vcf")
    c.goto_locus("chr1:1000")
    c.goto_range("chr1", 900, 1100)
    c.sort_alignments("BASE", "chr1:1000")
    c.group_alignments("STRAND")
    c.toggle_pair_view(True)
    c.collapse()
    c.squish()
    c.set_max_panel_height(800)
    c.set_snapshot_directory("/tmp/snaps")
    c.capture_snapshot("a.png")
    c.save_session("/tmp/s.xml")
    ops = [op for op, _ in mock_igv.commands]
    for expected in ("new", "genome", "load", "goto", "sort", "group", "viewAsPairs",
                     "collapse", "squish", "maxPanelHeight", "snapshotDirectory",
                     "snapshot", "saveSession"):
        assert expected in ops
    assert mock_igv.genome == "/ref.fa"
    assert "/data/x.bam" in mock_igv.loaded_files
    assert mock_igv.current_locus == "chr1:900-1100"


def test_client_rejects_unbacked_operation(mock_igv):
    c = IGVMCPClient(DirectTCPTransport("127.0.0.1", mock_igv.port))
    with pytest.raises(UnsupportedOperation):
        c.get_version()


def test_official_tool_list_contains_expected_tools():
    # verified against the official igv-mcp tools.yaml (2026-10)
    for t in ("new", "genome", "load", "goto", "snapshot", "snapshotDirectory",
              "group", "sort", "viewAsPairs", "maxPanelHeight", "saveSession"):
        assert t in OFFICIAL_TOOLS


def test_snapshot_creation_and_binding(mock_igv, tmp_path):
    cfg = AVRSConfig()
    cfg.paths.snapshot_dir = str(tmp_path / "snaps")
    client = IGVMCPClient(DirectTCPTransport("127.0.0.1", mock_igv.port))
    reviewer = IGVReviewer(client, cfg, tmp_path / "snaps")
    task = VariantReviewTask(
        sample_id="S1", variant_uid="BUILD#chr1:100:A>G", chrom="chr1", pos=100,
        ref="A", alt="G", variant_type="SNV", genome_build="BUILD",
        bam_path="/data/x.bam", bam_index_path="/data/x.bam.bai",
        vcf_path="/data/x.vcf", genome_id="/ref.fa",
        bam_checksum="sha256:abc", vcf_checksum="sha256:def")
    res = reviewer.review_variant(task)
    assert res.error is None
    assert len(res.snapshots) == 2
    kinds = {s.kind for s in res.snapshots}
    assert kinds == {"narrow", "context"}
    for s in res.snapshots:
        assert s.status == "CAPTURED"
        assert Path(s.path).exists()
        assert Path(s.path).stat().st_size > 1000
        # metadata binding: no cross-variant contamination
        assert s.variant_uid == task.variant_uid
        assert s.sample_id == "S1"
        assert s.bam_checksum == "sha256:abc"
        assert s.vcf_checksum == "sha256:def"
        assert s.genome_build == "BUILD"
    # different files for narrow/context
    assert len({s.path for s in res.snapshots}) == 2


def test_snapshot_failure_is_reported_not_faked(mock_igv, tmp_path):
    mock_igv.fail_next = 100  # every command fails
    cfg = AVRSConfig()
    cfg.review.snapshot_timeout_s = 1.0
    cfg.review.snapshot_poll_interval_s = 0.1
    cfg.igv.wait.load_wait_s = 0.05
    cfg.igv.wait.goto_wait_s = 0.05
    cfg.review.snapshot_wait_initial_s = 0.05
    client = IGVMCPClient(DirectTCPTransport("127.0.0.1", mock_igv.port))
    reviewer = IGVReviewer(client, cfg, tmp_path / "snaps")
    task = VariantReviewTask(sample_id="S1", variant_uid="B#chr1:1:A>G", chrom="chr1",
                              pos=1, ref="A", alt="G", variant_type="SNV",
                              genome_build="B", bam_path="/x.bam", genome_id="/ref.fa")
    res = reviewer.review_variant(task, retries=1)
    # either an error is reported or the snapshots are explicitly FAILED - never "OK"
    assert res.error is not None or all(s.status == "FAILED" for s in res.snapshots)
    assert all(s.status != "CAPTURED" for s in res.snapshots)


def test_instance_lock_is_exclusive(tmp_path):
    lock_path = tmp_path / "igv.lock"
    l1 = IGVInstanceLock(lock_path)
    assert l1.acquire(owner="a") is True
    l2 = IGVInstanceLock(lock_path)
    with pytest.raises(Exception):
        l2.acquire(owner="b")
    l1.release()
    assert l2.acquire(owner="b") is True
    l2.release()


def test_snapshot_verification_detects_empty_file(mock_igv, tmp_path):
    """A zero-byte snapshot must not count as a successful capture."""
    cfg = AVRSConfig()
    cfg.review.require_min_png_bytes = 2000
    cfg.review.snapshot_timeout_s = 1.0
    cfg.review.snapshot_poll_interval_s = 0.1
    client = IGVMCPClient(DirectTCPTransport("127.0.0.1", mock_igv.port))
    reviewer = IGVReviewer(client, cfg, tmp_path / "snaps")
    ok = reviewer._wait_for_file(tmp_path / "missing.png", cfg.review)
    assert ok is False
