"""Normalise line endings to LF for all text files in the repository.

Windows-side edits can leave CRLF behind, which breaks the Snakemake parser
(and creates noisy diffs). Run from the repo root: python scripts/normalize_eol.py
"""
from __future__ import annotations

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SKIP_DIRS = {".git", "outputs", "smoke_run", "__pycache__", ".workbuddy", "node_modules"}
SUFFIXES = {".py", ".md", ".yaml", ".yml", ".sh", ".html", ".js", ".css", ".toml", ".txt", ".cfg"}


def main() -> int:
    changed = []
    for path in ROOT.rglob("*"):
        if not path.is_file():
            continue
        if any(part in SKIP_DIRS for part in path.parts):
            continue
        if path.suffix not in SUFFIXES and path.name not in ("Snakefile", "Makefile"):
            continue
        data = path.read_bytes()
        if b"\r\n" in data:
            path.write_bytes(data.replace(b"\r\n", b"\n"))
            changed.append(str(path.relative_to(ROOT)))
    for c in changed:
        print("LF normalised:", c)
    print(f"{len(changed)} file(s) updated")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
