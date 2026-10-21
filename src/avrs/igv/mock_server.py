"""Mock IGV Desktop command-listener server for tests and demo fallback.

Implements the same wire protocol as the real listener (verified against
official igv-mcp): accepts `command arg...\\n`, replies one line ("OK" or
"Error: ..."), and `snapshot <file>` writes a real (tiny but valid) PNG so
the snapshot verification pipeline (existence, min size, stability) runs
exactly as against real IGV. The PNG is a generated placeholder image with
labeled bands - it is NEVER presented as a real IGV rendering.
"""
from __future__ import annotations

import socket
import socketserver
import struct
import threading
import zlib
from pathlib import Path


def make_placeholder_png(width: int = 640, height: int = 360, label: str = "MOCK-IGV") -> bytes:
    """Render a minimal valid PNG with deterministic noise.

    The noise keeps the compressed size well above the configured minimum
    PNG size used for snapshot verification, so the same validation path runs
    against mock and real IGV.
    """
    raw = bytearray()
    seed = sum(ord(c) for c in label) or 1
    for y in range(height):
        raw.append(0)  # filter type 0
        for x in range(width):
            band = (y // 40) % 3
            if band == 0:
                v = 235
            elif band == 1:
                v = 205
            else:
                v = 245
            # deterministic pseudo-random speckle (read-alignment-like texture)
            h = (x * 1103515245 + y * 12345 + seed * 2654435761) & 0xFFFFFFFF
            if (h >> 16) % 97 < 18:
                v = 40 + (h >> 8) % 120
            if abs(y - 320) < 6 and (x // 8) % 2 == 0:  # "text stripe"
                v = 30
            raw.extend([v, v, v])
    raw = bytes(raw)

    def chunk(tag: bytes, data: bytes) -> bytes:
        c = struct.pack(">I", len(data)) + tag + data
        return c + struct.pack(">I", zlib.crc32(tag + data) & 0xFFFFFFFF)

    ihdr = struct.pack(">IIBBBBB", width, height, 8, 2, 0, 0, 0)
    return (b"\x89PNG\r\n\x1a\n"
            + chunk(b"IHDR", ihdr)
            + chunk(b"IDAT", zlib.compress(raw, 6))
            + chunk(b"IEND", b""))


class _Handler(socketserver.BaseRequestHandler):
    def handle(self):
        f = self.request.makefile("r")
        for line in f:
            cmd = line.strip()
            if not cmd:
                continue
            parts = cmd.split(" ", 1)
            op = parts[0]
            arg = parts[1] if len(parts) > 1 else ""
            try:
                resp = self.server.dispatch(op, arg)  # type: ignore[attr-defined]
            except Exception as e:  # noqa: BLE001 - mock must never die
                resp = f"Error: {e}"
            self.request.sendall((resp + "\n").encode("utf-8"))
        f.close()


class MockIGVServer(socketserver.ThreadingTCPServer):
    """A TCP server that mimics the IGV Desktop command listener."""

    allow_reuse_address = True
    daemon_threads = True

    def __init__(self, host: str = "127.0.0.1", port: int = 0):
        super().__init__((host, port), _Handler)
        self.loaded_files: list[str] = []
        self.genome: str | None = None
        self.current_locus: str | None = None
        self.commands: list[tuple[str, str]] = []
        self.snapshot_dir: str | None = None
        self.fail_next = 0
        self._lock = threading.Lock()
        self.port = self.server_address[1]

    def dispatch(self, op: str, arg: str) -> str:
        with self._lock:
            self.commands.append((op, arg))
        if self.fail_next > 0:
            self.fail_next -= 1
            return "Error: simulated failure"
        if op == "new":
            self.loaded_files.clear()
            return "OK"
        if op == "genome":
            self.genome = arg
            return "OK"
        if op == "load":
            self.loaded_files.extend(arg.split())
            return "OK"
        if op == "goto":
            self.current_locus = arg
            return "OK"
        if op in ("collapse", "squish", "expand", "zoomin", "zoomout"):
            return "OK"
        if op in ("sort", "group", "viewAsPairs", "setColor", "region"):
            return "OK"
        if op == "maxPanelHeight":
            return "OK"
        if op == "snapshotDirectory":
            self.snapshot_dir = arg
            return "OK"
        if op == "snapshot":
            target = arg
            if not target and self.snapshot_dir:
                target = str(Path(self.snapshot_dir) / "igv_snapshot.png")
            if not target:
                return "Error: no snapshot target and no snapshotDirectory set"
            p = Path(target)
            if not p.is_absolute() and self.snapshot_dir:
                # IGV resolves relative snapshot names inside the snapshot directory
                p = Path(self.snapshot_dir) / p
            p.parent.mkdir(parents=True, exist_ok=True)
            label = self.current_locus or "unknown-locus"
            p.write_bytes(make_placeholder_png(label=label))
            return "OK"
        if op == "saveSession":
            p = Path(arg)
            p.parent.mkdir(parents=True, exist_ok=True)
            p.write_text("<session>mock</session>", encoding="utf-8")
            return "OK"
        return f"Error: unknown command '{op}'"

    def start_background(self):
        t = threading.Thread(target=self.serve_forever, daemon=True)
        t.start()
        return t

    def stop(self):
        self.shutdown()
        self.server_close()


def start_mock_igv(port: int = 0) -> MockIGVServer:
    srv = MockIGVServer(port=port)
    srv.start_background()
    return srv
