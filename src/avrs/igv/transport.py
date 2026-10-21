"""IGV transports, faithfully implementing the official igv-mcp protocol.

Verified against https://github.com/igvteam/igv-mcp (src/mcpServer.js,
src/igvConnection.js, src/tools.yaml.js):

- Direct transport: TCP socket to IGV Desktop command listener (default
  60151); send `<command> <arg...>\\n`, read single-line response
  (empty response == "OK").
- MCP transport: spawn the official igv-mcp Node server over stdio and
  speak JSON-RPC 2.0 (initialize -> tools/list -> tools/call).

Only tools that actually exist in the official tool list are exposed.
The application-level client never invents parameters.
"""
from __future__ import annotations

import json
import subprocess
import socket
import time
from pathlib import Path
from typing import Optional, Protocol

# Tool names as defined by the official igv-mcp tools.yaml (verified 2026-10).
OFFICIAL_TOOLS = {
    "new", "saveSession", "genome", "load", "goto", "collapse", "squish",
    "expand", "region", "setColor", "group", "sort", "viewAsPairs",
    "snapshot", "snapshotDirectory", "maxPanelHeight", "setSequenceStrand",
    "setSequenceShowTranslation", "overlay", "separate", "zoomin", "zoomout",
}


class IGVError(Exception):
    pass


class IGVConnectionError(IGVError):
    pass


class IGVCommandError(IGVError):
    pass


class UnsupportedOperation(IGVError):
    """The app-level operation has no backing official MCP tool."""


class Transport(Protocol):
    def send_command(self, command: str) -> str: ...
    def close(self) -> None: ...


class DirectTCPTransport:
    """Same wire protocol the official igv-mcp uses under the hood."""

    def __init__(self, host: str = "127.0.0.1", port: int = 60151,
                 timeout: float = 30.0):
        self.host = host
        self.port = port
        self.timeout = timeout

    def connectivity_check(self, timeout: float = 3.0) -> bool:
        try:
            with socket.create_connection((self.host, self.port), timeout=timeout):
                return True
        except OSError:
            return False

    def send_command(self, command: str) -> str:
        try:
            with socket.create_connection((self.host, self.port), timeout=self.timeout) as sock:
                sock.settimeout(self.timeout)
                sock.sendall(f"{command}\n".encode("utf-8"))
                response = b""
                while not response.endswith(b"\n"):
                    chunk = sock.recv(4096)
                    if not chunk:
                        break
                    response += chunk
                text = response.decode("utf-8", errors="replace").strip()
                text = text or "OK"
                if text.lower().startswith("error"):
                    # official igv-mcp treats any non-"OK" reply as an error
                    raise IGVCommandError(f"IGV rejected command {command!r}: {text}")
                return text
        except OSError as e:
            raise IGVConnectionError(
                f"Cannot reach IGV at {self.host}:{self.port}: {e}. "
                "Start IGV Desktop and enable View > Preferences > Advanced > "
                "port command listener (default 60151)."
            ) from e

    def close(self) -> None:
        pass  # one connection per command (same as official igv-mcp)


class MCPStdioTransport:
    """Drive the official igv-mcp Node server over stdio JSON-RPC."""

    def __init__(self, command: str, args: list[str], timeout: float = 60.0):
        self.command = command
        self.args = args
        self.timeout = timeout
        self.proc: Optional[subprocess.Popen] = None
        self._id = 0
        self._tools: list[dict] = []

    def start(self):
        if self.proc and self.proc.poll() is None:
            return
        try:
            self.proc = subprocess.Popen(
                [self.command, *self.args],
                stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                stderr=subprocess.PIPE, text=True, encoding="utf-8",
            )
        except OSError as e:
            raise IGVConnectionError(f"Cannot start igv-mcp server {self.command} {self.args}: {e}") from e
        self._handshake()

    def _handshake(self):
        result = self._request("initialize", {
            "protocolVersion": "2024-11-05",
            "capabilities": {},
            "clientInfo": {"name": "avrs", "version": "0.1.0"},
        })
        server_info = result.get("serverInfo", {})
        self.server_name = server_info.get("name", "unknown")
        self.server_version = server_info.get("version", "unknown")
        self._notify("notifications/initialized")

    def _request(self, method: str, params: dict) -> dict:
        self._id += 1
        msg = {"jsonrpc": "2.0", "id": self._id, "method": method, "params": params}
        assert self.proc and self.proc.stdin and self.proc.stdout
        self.proc.stdin.write(json.dumps(msg) + "\n")
        self.proc.stdin.flush()
        deadline = time.time() + self.timeout
        while time.time() < deadline:
            line = self.proc.stdout.readline()
            if not line:
                raise IGVConnectionError("igv-mcp server closed stdio")
            data = json.loads(line)
            if data.get("id") == self._id:
                if "error" in data:
                    raise IGVCommandError(f"igv-mcp error: {data['error']}")
                return data.get("result", {})
        raise IGVConnectionError(f"Timeout waiting for {method} response")

    def _notify(self, method: str, params: dict | None = None):
        msg = {"jsonrpc": "2.0", "method": method}
        if params:
            msg["params"] = params
        assert self.proc and self.proc.stdin
        self.proc.stdin.write(json.dumps(msg) + "\n")
        self.proc.stdin.flush()

    def list_tools(self) -> list[dict]:
        result = self._request("tools/list", {})
        self._tools = result.get("tools", [])
        return self._tools

    def send_command(self, command: str) -> str:
        """Translate an app-level command string into a tools/call invocation.

        command must be `<tool> <arg>...` with a tool from OFFICIAL_TOOLS.
        """
        parts = command.split(" ", 1)
        tool = parts[0]
        if tool not in OFFICIAL_TOOLS:
            raise UnsupportedOperation(f"Tool {tool!r} is not part of the official igv-mcp tool set")
        if not self._tools and self.proc:
            self.list_tools()
        tool_def = next((t for t in self._tools if t["name"] == tool), None)
        args: dict[str, str] = {}
        if tool_def:
            prop_names = list(tool_def.get("inputSchema", {}).get("properties", {}).keys())
            if len(parts) > 1:
                vals = parts[1].split(" ")
                for name, val in zip(prop_names, vals):
                    args[name] = val
        result = self._request("tools/call", {"name": tool, "arguments": args})
        content = result.get("content", [])
        text = content[0].get("text", "") if content else ""
        if result.get("isError"):
            raise IGVCommandError(f"IGV tool {tool} failed: {text}")
        return text or "OK"

    def close(self) -> None:
        if self.proc and self.proc.poll() is None:
            try:
                self.proc.terminate()
                self.proc.wait(timeout=5)
            except Exception:
                self.proc.kill()
        self.proc = None


def make_transport(igv_cfg) -> Transport:
    from ..config import IGVConfig
    if isinstance(igv_cfg, IGVConfig):
        if igv_cfg.transport == "mcp":
            t = MCPStdioTransport(igv_cfg.mcp_command, igv_cfg.mcp_args,
                                  timeout=igv_cfg.wait.command_timeout_s)
            t.start()
            return t
        return DirectTCPTransport(igv_cfg.host, igv_cfg.port,
                                  timeout=igv_cfg.wait.command_timeout_s)
    raise TypeError("igv_cfg must be an IGVConfig")
