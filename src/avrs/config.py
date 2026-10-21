"""Typed configuration system (YAML + Pydantic) with friendly validation errors."""
from __future__ import annotations

import os
from pathlib import Path
from typing import Any, Literal, Optional

import yaml
from pydantic import BaseModel, Field, field_validator, model_validator


class ConfigError(Exception):
    """Raised when configuration is invalid."""


class ReferenceConfig(BaseModel):
    build: str = "GRCh38"  # e.g. GRCh38 / hg38 - part of variant identity
    fasta: Optional[str] = None
    chr_prefix: Optional[str] = None  # None = auto-detect from FASTA/VCF; never blind-rewrite
    known_sites: list[str] = Field(default_factory=list)  # BQSR only if provided

    @property
    def fasta_path(self) -> Optional[Path]:
        return Path(self.fasta).expanduser() if self.fasta else None


class PathConfig(BaseModel):
    output_dir: str = "outputs"
    snapshot_dir: str = "outputs/snapshots"
    db: str = "outputs/avrs.db"
    path_mapping: list[dict[str, str]] = Field(default_factory=list)
    """e.g. [{from: '/mnt/c/data', to: 'C:/data'}] for WSL2<->Windows translation.
    Applied to paths sent to IGV (which may run on Windows while the analysis
    pipeline runs in WSL2). Configured, never hardcoded."""


class AnalysisConfig(BaseModel):
    mode: Literal["review_only", "full_pipeline"] = "review_only"
    sequencing: Literal["germline", "somatic"] = "germline"
    data_type: Literal["wgs", "wes", "panel"] = "panel"
    target_bed: Optional[str] = None
    platform: str = "ILLUMINA"

    @field_validator("target_bed")
    @classmethod
    def _bed_exists(cls, v):
        if v and not Path(v).exists():
            raise ConfigError(f"target_bed does not exist: {v}")
        return v


class TriageRuleConfig(BaseModel):
    enabled: bool = True
    severity: Literal["info", "warning", "high"] = "warning"
    thresholds: dict[str, Any] = Field(default_factory=dict)


class TriageProfile(BaseModel):
    min_dp: int = 10
    min_gq: int = 20
    min_alt_ad: int = 3
    min_mapq: int = 20
    min_baseq: int = 13
    max_strand_ratio_imbalance: float = 0.95  # fraction of ALT on one strand
    max_softclip_frac: float = 0.3
    max_dup_frac: float = 0.2
    min_dist_from_read_end: int = 5
    max_start_pos_cluster_frac: float = 0.5
    max_vaf_dev: float = 0.35  # |VAF - expected for GT| tolerance
    homopolymer_len: int = 6
    low_complexity_window: int = 30
    min_target_edge_dist: int = 10
    rules: dict[str, TriageRuleConfig] = Field(default_factory=dict)


def _default_profiles() -> dict[str, TriageProfile]:
    """Per-data-type starting points. These are DEFAULTS, not universal
    standards - every value is meant to be overridden per laboratory."""
    wgs = TriageProfile(min_dp=20, min_gq=30, min_alt_ad=5, min_mapq=30,
                        min_baseq=20, max_strand_ratio_imbalance=0.90,
                        max_softclip_frac=0.25, max_dup_frac=0.15,
                        min_dist_from_read_end=5, max_start_pos_cluster_frac=0.40,
                        max_vaf_dev=0.30)
    wes = TriageProfile(min_dp=15, min_gq=25, min_alt_ad=4, min_mapq=25,
                        min_baseq=15, max_strand_ratio_imbalance=0.90,
                        max_softclip_frac=0.25, max_dup_frac=0.20,
                        min_dist_from_read_end=5, max_start_pos_cluster_frac=0.45,
                        max_vaf_dev=0.35)
    panel = TriageProfile(min_dp=30, min_gq=40, min_alt_ad=8, min_mapq=30,
                          min_baseq=20, max_strand_ratio_imbalance=0.85,
                          max_softclip_frac=0.20, max_dup_frac=0.10,
                          min_dist_from_read_end=8, max_start_pos_cluster_frac=0.35,
                          max_vaf_dev=0.25, homopolymer_len=5,
                          low_complexity_window=25, min_target_edge_dist=15)
    return {"wgs": wgs, "wes": wes, "panel": panel}


class TriageConfig(BaseModel):
    profiles: dict[str, TriageProfile] = Field(default_factory=_default_profiles)

    def profile_for(self, data_type: str) -> TriageProfile:
        """Profile for the analysis type; falls back to any configured profile."""
        if data_type in self.profiles:
            return self.profiles[data_type]
        if "wes" in self.profiles:
            return self.profiles["wes"]
        if self.profiles:
            return next(iter(self.profiles.values()))
        return TriageProfile()


class ReviewConfig(BaseModel):
    modes: list[Literal[
        "all", "suspicious_only", "target_genes", "filtered_calls",
        "bed_regions", "selected_variants"]] = Field(default_factory=lambda: ["suspicious_only"])
    target_genes: list[str] = Field(default_factory=list)
    bed_regions: Optional[str] = None
    batch_size: int = 50
    max_retries: int = 3
    # render-wait: never rely on a fixed sleep as the only success criterion
    snapshot_wait_initial_s: float = 1.5
    snapshot_poll_interval_s: float = 0.5
    snapshot_timeout_s: float = 30.0
    require_file_stable_bytes: int = 1024
    require_min_png_bytes: int = 2000
    window_width_narrow: int = 60  # bp, dynamically adjusted by variant type
    window_width_context: int = 500
    window_width_indel_pad: int = 30


class IGVWaitConfig(BaseModel):
    load_wait_s: float = 3.0
    goto_wait_s: float = 2.0
    command_timeout_s: float = 30.0


class IGVConfig(BaseModel):
    transport: Literal["direct", "mcp"] = "direct"
    host: str = "127.0.0.1"
    port: int = 60151
    # For transport=mcp: the official igv-mcp server (https://github.com/igvteam/igv-mcp)
    mcp_command: Optional[str] = None   # e.g. "node"
    mcp_args: list[str] = Field(default_factory=list)  # e.g. ["/path/igv-mcp/src/main.js", "port", "60151"]
    # Enable ONLY when the IGV instance runs on a different filesystem than
    # the analysis (e.g. IGV Desktop on Windows, pipeline in WSL2). When the
    # analysis and IGV share the same filesystem (e.g. a mock/local IGV inside
    # WSL2) keep this false so native paths are used unchanged.
    path_mapping_enabled: bool = False
    wait: IGVWaitConfig = Field(default_factory=IGVWaitConfig)
    max_panel_height: int = 1000

    @model_validator(mode="after")
    def _mcp_requires_command(self):
        if self.transport == "mcp" and not self.mcp_command:
            raise ConfigError("transport=mcp requires mcp_command (e.g. 'node')")
        return self


class AIProviderConfig(BaseModel):
    provider: Literal["mock", "openai_compatible", "ollama"] = "mock"
    model: str = "mock-reviewer-v1"
    base_url: Optional[str] = None
    api_key_env: Optional[str] = None  # name of env var holding the key; key never logged/stored
    supports_vision: bool = False
    # Security: sending screenshots to a remote provider is opt-in and explicit.
    allow_remote_screenshot_transfer: bool = False
    timeout_s: float = 120.0
    retries: int = 2
    temperature: float = 0.0
    max_output_tokens: int = 2000


class SecurityConfig(BaseModel):
    local_first: bool = True
    redact_names_in_requests: bool = True
    max_llm_context_variants: int = 20  # WGS can have millions of calls - never bulk-load


class LoggingConfig(BaseModel):
    level: Literal["DEBUG", "INFO", "WARNING", "ERROR"] = "INFO"
    file: Optional[str] = None


class AVRSConfig(BaseModel):
    project_name: str = "avrs-project"
    reference: ReferenceConfig = Field(default_factory=ReferenceConfig)
    paths: PathConfig = Field(default_factory=PathConfig)
    analysis: AnalysisConfig = Field(default_factory=AnalysisConfig)
    triage: TriageConfig = Field(default_factory=TriageConfig)
    review: ReviewConfig = Field(default_factory=ReviewConfig)
    igv: IGVConfig = Field(default_factory=IGVConfig)
    ai: AIProviderConfig = Field(default_factory=AIProviderConfig)
    security: SecurityConfig = Field(default_factory=SecurityConfig)
    logging: LoggingConfig = Field(default_factory=LoggingConfig)
    samples: Optional[str] = None  # manifest path

    def resolve_output(self, p: str | Path) -> Path:
        """Resolve an output-relative path (idempotent: passing the output
        dir itself does not duplicate it)."""
        base = Path(self.paths.output_dir).expanduser()
        path = Path(p).expanduser()
        if path.is_absolute():
            return path
        try:
            if path.parts and base.parts and path.parts[0] == base.parts[0]:
                return path
        except IndexError:
            pass
        return base / path

    def map_path_for_igv(self, p: str | Path) -> str:
        """Translate an analysis-side path to the IGV-side filesystem (WSL2<->Windows etc.)."""
        s = str(Path(p).expanduser().resolve()) if Path(p).exists() else str(p)
        for m in self.paths.path_mapping:
            src, dst = m.get("from", ""), m.get("to", "")
            if src and s.startswith(src):
                s = dst + s[len(src):]
                s = s.replace("\\", "/") if dst and (":" in dst) else s
                break
        return s


def load_config(path: str | Path) -> AVRSConfig:
    path = Path(path)
    if not path.exists():
        raise ConfigError(f"Config file not found: {path}")
    try:
        raw = yaml.safe_load(path.read_text(encoding="utf-8"))
    except yaml.YAMLError as e:
        raise ConfigError(f"Invalid YAML in {path}: {e}") from e
    if not isinstance(raw, dict):
        raise ConfigError(f"Config root must be a mapping, got {type(raw).__name__}")
    try:
        return AVRSConfig.model_validate(raw)
    except Exception as e:  # pydantic ValidationError
        raise ConfigError(f"Invalid configuration in {path}:\n{e}") from e


def setup_logging(cfg: LoggingConfig):
    import logging
    fmt = "%(asctime)s %(levelname)s [%(name)s] %(message)s"
    handlers = [logging.StreamHandler()]
    if cfg.file:
        Path(cfg.file).parent.mkdir(parents=True, exist_ok=True)
        handlers.append(logging.FileHandler(cfg.file, encoding="utf-8"))
    logging.basicConfig(level=cfg.level, format=fmt, handlers=handlers, force=True)
    logging.getLogger("httpx").setLevel("WARNING")
    return logging.getLogger("avrs")
