"""Application-level IGV MCP client.

App-level method names follow the spec (connect/health_check/load_genome/
load_alignment/goto_locus/...); each maps onto an *actual* tool from the
official igv-mcp tool set. Operations without a backing tool raise
UnsupportedOperation - never silently faked.
"""
from __future__ import annotations

import uuid
from pathlib import Path
from typing import Optional

from .transport import (
    DirectTCPTransport, MCPStdioTransport, Transport, UnsupportedOperation,
    IGVConnectionError, IGVCommandError, OFFICIAL_TOOLS,
)


class IGVMCPClient:
    def __init__(self, transport: Transport, wait_load_s: float = 3.0,
                 wait_goto_s: float = 2.0):
        self.transport = transport
        self.wait_load_s = wait_load_s
        self.wait_goto_s = wait_goto_s
        self.session_token = uuid.uuid4().hex[:16]

    # ---- connection / session -------------------------------------------
    def connect(self):
        if isinstance(self.transport, MCPStdioTransport):
            tools = self.transport.list_tools()
            names = [t["name"] for t in tools]
            missing = OFFICIAL_TOOLS - set(names)
            if missing:
                raise IGVConnectionError(
                    f"igv-mcp tool list missing expected tools: {sorted(missing)}")
            return tools
        ok = self.health_check()
        if not ok:
            raise IGVConnectionError("IGV Desktop command port is not reachable")
        return None

    def health_check(self) -> bool:
        """Liveness: TCP port listener accepting connections (direct) or a
        successful MCP initialize+tools/list (mcp). A positive 'version'
        probe is not part of the official tool set, so it is never faked."""
        if isinstance(self.transport, DirectTCPTransport):
            return self.transport.connectivity_check()
        try:
            self.transport.list_tools()
            return True
        except IGVError:
            return False

    def reset_session(self):
        """`new` - unload all user tracks (keeps genome annotations)."""
        return self.transport.send_command("new")

    def save_session(self, filename: str):
        return self.transport.send_command(f"saveSession {filename}")

    # ---- data loading ----------------------------------------------------
    def load_genome(self, id_or_path: str):
        return self.transport.send_command(f"genome {id_or_path}")

    def load_alignment(self, bam_path: str, index_path: Optional[str] = None):
        if index_path:
            return self.transport.send_command(f"load {bam_path} {index_path}")
        return self.transport.send_command(f"load {bam_path}")

    def load_vcf(self, vcf_path: str):
        return self.transport.send_command(f"load {vcf_path}")

    def load_bed(self, bed_path: str):
        return self.transport.send_command(f"load {bed_path}")

    # ---- navigation / view -----------------------------------------------
    def goto_locus(self, locus: str):
        return self.transport.send_command(f"goto {locus}")

    def goto_position(self, chrom: str, pos: int):
        # IGV centers on a single position with a ~40bp window; for wider
        # context use goto_range.
        return self.goto_locus(f"{chrom}:{pos}")

    def goto_range(self, chrom: str, start: int, end: int):
        return self.goto_locus(f"{chrom}:{start}-{end}")

    def collapse(self, track_name: Optional[str] = None):
        return self.transport.send_command(f"collapse {track_name}".strip())

    def squish(self, track_name: Optional[str] = None):
        return self.transport.send_command(f"squish {track_name}".strip())

    def expand(self, track_name: Optional[str] = None):
        return self.transport.send_command(f"expand {track_name}".strip())

    def set_color(self, color: str, track_name: Optional[str] = None):
        cmd = f"setColor {color} {track_name}".strip() if track_name else f"setColor {color}"
        return self.transport.send_command(cmd)

    def sort_alignments(self, attribute: str, locus: Optional[str] = None,
                        direction: Optional[str] = None):
        """attribute: BASE|STRAND|POSITION|QUALITY|MATECHR|INSERTSIZE|..."""
        parts = ["sort", attribute]
        if direction and not locus:
            raise ValueError("direction requires locus (official tool constraint)")
        if locus:
            parts.append(locus)
        if direction:
            parts.append(direction)
        return self.transport.send_command(" ".join(parts))

    def group_alignments(self, attribute: str, tag: Optional[str] = None):
        """attribute: STRAND|BASE|READ_GROUP|SAMPLE|HAPLOTYPE|..."""
        return self.transport.send_command(f"group {attribute} {tag}".strip() if tag else f"group {attribute}")

    def toggle_pair_view(self, enabled: bool):
        return self.transport.send_command(f"viewAsPairs {'true' if enabled else 'false'}")

    def zoom_in(self):
        return self.transport.send_command("zoomin")

    def zoom_out(self):
        return self.transport.send_command("zoomout")

    def region_of_interest(self, chrom: str, start: int, end: int):
        return self.transport.send_command(f"region {chrom} {start} {end}")

    # ---- snapshots ---------------------------------------------------------
    def set_snapshot_directory(self, directory: str):
        return self.transport.send_command(f"snapshotDirectory {directory}")

    def set_max_panel_height(self, height: int):
        return self.transport.send_command(f"maxPanelHeight {height}")

    def capture_snapshot(self, filename: Optional[str] = None):
        cmd = f"snapshot {filename}".strip() if filename else "snapshot"
        return self.transport.send_command(cmd)

    # ---- explicitly unsupported (no official tool exists) -----------------
    def get_version(self) -> str:
        """No 'version' tool exists in official igv-mcp - report honestly."""
        raise UnsupportedOperation(
            "IGV version query is not part of the official igv-mcp tool set; "
            "read it manually from IGV Help > About.")

    def disconnect(self):
        self.transport.close()


class IGVInstanceLock:
    """Serialize access to the single IGV GUI instance (file lock).

    One IGV Desktop instance must never be driven by two concurrent tasks
    (spec 15.3). Lock file carries owner PID + timestamp; stale locks older
    than max_age_s can be reclaimed.
    """
    def __init__(self, lock_path: str | Path, max_age_s: float = 3600.0):
        self.lock_path = Path(lock_path)
        self.max_age_s = max_age_s
        self._acquired = False

    def acquire(self, owner: str = "avrs") -> bool:
        self.lock_path.parent.mkdir(parents=True, exist_ok=True)
        if self.lock_path.exists():
            import time
            age = time.time() - self.lock_path.stat().st_mtime
            content = self.lock_path.read_text(encoding="utf-8", errors="replace")
            if age < self.max_age_s:
                raise IGVError(
                    f"IGV instance lock held by {content!r} ({age:.0f}s old) at {self.lock_path}")
            # stale
            try:
                self.lock_path.unlink()
            except OSError:
                pass
        import os
        import time
        try:
            fd = os.open(self.lock_path, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
        except FileExistsError:
            raise IGVError("IGV instance lock contention (concurrent acquire)")
        with os.fdopen(fd, "w", encoding="utf-8") as fh:
            fh.write(f"{owner} pid={os.getpid()} ts={time.time()}")
        self._acquired = True
        return True

    def release(self):
        if self._acquired:
            try:
                self.lock_path.unlink(missing_ok=True)
            except OSError:
                pass
            self._acquired = False

    def __enter__(self):
        self.acquire()
        return self

    def __exit__(self, *exc):
        self.release()
