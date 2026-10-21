"""FASTQ -> VCF pipeline support: input validation and Snakemake runner.

Validation is strict by design: reference/VCF naming mismatches are reported,
never silently rewritten (spec 3.1).
"""
from __future__ import annotations

import gzip
import hashlib
import json
import logging
import shutil
import subprocess
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

log = logging.getLogger("avrs.pipeline")


@dataclass
class ValidationIssue:
    level: str      # error | warning | info
    code: str
    message: str


@dataclass
class ValidationResult:
    issues: list[ValidationIssue] = field(default_factory=list)

    @property
    def errors(self):
        return [i for i in self.issues if i.level == "error"]

    @property
    def ok(self) -> bool:
        return not self.errors

    def add(self, level, code, message):
        self.issues.append(ValidationIssue(level, code, message))

    def report(self) -> str:
        if not self.issues:
            return "All input checks passed."
        lines = []
        for i in self.issues:
            lines.append(f"[{i.level.upper():7}] {i.code}: {i.message}")
        return "\n".join(lines)


def _peek_fastq(path: Path, n: int = 200) -> list[str]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", errors="replace") as fh:
        return [next(fh, "") for _ in range(n)]


def validate_fastq_pair(r1: Path, r2: Path, result: ValidationResult):
    for p, tag in ((r1, "R1"), (r2, "R2")):
        if not p.exists():
            result.add("error", "FASTQ_MISSING", f"{tag} not found: {p}")
            continue
        if p.stat().st_size == 0:
            result.add("error", "FASTQ_EMPTY", f"{tag} is empty: {p}")
            continue
        lines = _peek_fastq(p)
        if len(lines) < 4:
            result.add("error", "FASTQ_TRUNCATED", f"{tag} has fewer than 4 lines: {p}")
            continue
        if not lines[0].startswith("@"):
            result.add("error", "FASTQ_FORMAT", f"{tag} does not start with '@': {p}")
        if not lines[2].startswith("+"):
            result.add("error", "FASTQ_FORMAT", f"{tag} line 3 is not '+': {p}")
        seq_len = len(lines[1].strip())
        qual_len = len(lines[3].strip())
        if seq_len != qual_len:
            result.add("error", "FASTQ_LENGTH_MISMATCH",
                       f"{tag} sequence length {seq_len} != quality length {qual_len}")
    if r1.exists() and r2.exists():
        n1 = _count_reads(r1)
        n2 = _count_reads(r2)
        if n1 != n2:
            result.add("error", "FASTQ_PAIR_COUNT",
                       f"R1 has {n1} reads, R2 has {n2} reads - files are not a matched pair")
        else:
            result.add("info", "FASTQ_PAIR_OK", f"{n1} read pairs")


def _count_reads(path: Path) -> int:
    opener = gzip.open if path.suffix == ".gz" else open
    n = 0
    with opener(path, "rb") as fh:
        for i, _ in enumerate(fh):
            if i % 4 == 0:
                n += 1
    return n


def validate_reference(fasta: Optional[Path], vcf: Optional[Path],
                        result: ValidationResult):
    if not fasta:
        result.add("warning", "NO_REFERENCE", "No reference FASTA configured")
        return
    if not fasta.exists():
        result.add("error", "REFERENCE_MISSING", f"Reference FASTA not found: {fasta}")
        return
    if not Path(str(fasta) + ".fai").exists():
        result.add("warning", "REFERENCE_NO_INDEX",
                   f"Missing FASTA index {fasta}.fai (run samtools faidx)")
    if vcf and Path(vcf).exists():
        from ..io.vcf import VcfReader
        try:
            reader = VcfReader(vcf)
            vcf_contigs = set(reader.header.contigs)
        except Exception as e:
            result.add("error", "VCF_HEADER", f"Cannot read VCF header: {e}")
            return
        try:
            from ..io.fasta import Fasta
            fa = Fasta(fasta)
            ref_contigs = set(fa.contig_names)
        except Exception:
            ref_contigs = set()
        if ref_contigs and vcf_contigs:
            missing = vcf_contigs - ref_contigs
            if missing:
                style = "chr-prefix" if all(c.startswith("chr") for c in missing) else "other"
                result.add("error", "CONTIG_NAMING_MISMATCH",
                           f"VCF contigs absent from reference ({style}): "
                           f"{sorted(missing)[:10]}. Fix upstream (bcftools reheader/annotate) - "
                           f"AVRS never rewrites contig names automatically.")


def validate_known_sites(known: list[str], result: ValidationResult):
    for k in known:
        if not Path(k).exists():
            result.add("error", "KNOWN_SITES_MISSING",
                       f"known-sites resource missing (needed for BQSR): {k}")


def checksum(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return f"sha256:{h.hexdigest()}"


def run_pipeline(pipeline_config: str, sample: Optional[str] = None,
                 dry_run: bool = False, cores: int = 4) -> dict:
    """Run the Snakemake workflow (workflow/Snakefile) with the given config."""
    import yaml
    cfg_path = Path(pipeline_config)
    if not cfg_path.exists():
        raise FileNotFoundError(cfg_path)
    raw = yaml.safe_load(cfg_path.read_text(encoding="utf-8")) or {}
    root = Path(raw.get("workflow_dir", Path(__file__).resolve().parents[2] / "workflow"))
    snakefile = root / "Snakefile"
    if not snakefile.exists():
        raise FileNotFoundError(f"Snakefile not found: {snakefile}")
    if not shutil.which("snakemake"):
        raise RuntimeError("snakemake is not on PATH (install in WSL2/Linux env)")
    workdir = cfg_path.parent / "snakemake_work"
    workdir.mkdir(parents=True, exist_ok=True)
    cmd = [
        "snakemake", "--snakefile", str(snakefile), "--configfile", str(cfg_path),
        "--cores", str(cores), "--directory", str(workdir),
        "--rerun-incomplete", "--printshellcmds",
    ]
    if sample:
        cmd += ["--allowed-rules"]  # placeholder to keep snakemake invocation explicit
    if dry_run:
        cmd.append("--dry-run")
    console_msg = " ".join(cmd)
    log.info("running: %s", console_msg)
    t0 = time.time()
    proc = subprocess.run(cmd, capture_output=True, text=True)
    duration = round(time.time() - t0, 2)
    log_path = workdir / "snakemake.log"
    log_path.write_text((proc.stdout or "") + "\n" + (proc.stderr or ""), encoding="utf-8")
    return {
        "command": console_msg, "returncode": proc.returncode,
        "duration_s": duration, "log": str(log_path),
        "stdout_tail": (proc.stdout or "")[-4000:], "stderr_tail": (proc.stderr or "")[-4000:],
    }
