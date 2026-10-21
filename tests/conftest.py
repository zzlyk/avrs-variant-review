import os
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from avrs.config import AVRSConfig, TriageProfile  # noqa: E402
from avrs.db import open_db  # noqa: E402


@pytest.fixture
def cfg(tmp_path) -> AVRSConfig:
    c = AVRSConfig()
    c.paths.output_dir = str(tmp_path)
    c.paths.snapshot_dir = str(tmp_path / "snapshots")
    c.paths.db = str(tmp_path / "avrs.db")
    c.reference.build = "TESTBUILD"
    return c


@pytest.fixture
def db_session(cfg):
    engine, session = open_db(cfg.paths.db)
    yield session
    session.close()
    engine.dispose()


@pytest.fixture
def demo_dataset(tmp_path):
    from avrs.demo.generator import generate_demo
    manifest = generate_demo(tmp_path / "demo")
    return manifest
