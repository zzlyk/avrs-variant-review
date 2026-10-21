"""AVRS command line interface.

CLI and the web UI call the same service layer (avrs.review / avrs.reports),
so behaviour cannot drift between them.
"""
from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path
from typing import Optional

import typer
from rich.console import Console
from rich.table import Table

from . import __version__, DISCLAIMER
from .config import AVRSConfig, ConfigError, load_config, setup_logging

app = typer.Typer(add_completion=False, help="AVRS - AI Variant Review System")
variants_app = typer.Typer(help="Variant import and triage")
pipeline_app = typer.Typer(help="FASTQ -> VCF pipeline")
app.add_typer(variants_app, name="variants")
app.add_typer(pipeline_app, name="pipeline")

console = Console()
CONFIG_ENV = "AVRS_CONFIG"


def _config(path: Optional[str]) -> AVRSConfig:
    cfg_path = path or os.environ.get(CONFIG_ENV)
    if not cfg_path:
        raise typer.BadParameter(
            "No config file. Use --config /path/config.yaml or set AVRS_CONFIG.")
    try:
        cfg = load_config(cfg_path)
    except ConfigError as e:
        console.print(f"[red]Configuration error:[/red] {e}")
        raise typer.Exit(2)
    setup_logging(cfg.logging)
    return cfg


def _db(cfg: AVRSConfig):
    from .db import open_db
    return open_db(cfg.resolve_output(cfg.paths.db))


# --------------------------------------------------------------------- init
@app.command()
def init(config: str = typer.Option("config/config.yaml", help="Path of the config to create"),
         force: bool = typer.Option(False, help="Overwrite existing file")):
    """Create a starter configuration file and output directory layout."""
    p = Path(config)
    if p.exists() and not force:
        console.print(f"[yellow]{p} already exists (use --force to overwrite)[/yellow]")
        raise typer.Exit(1)
    p.parent.mkdir(parents=True, exist_ok=True)
    example = Path(__file__).resolve().parents[2] / "config" / "config.example.yaml"
    if example.exists():
        shutil.copyfile(example, p)
    else:
        p.write_text("project_name: avrs-project\n", encoding="utf-8")
    for sub in ("outputs", "outputs/snapshots", "outputs/reports"):
        Path(sub).mkdir(parents=True, exist_ok=True)
    console.print(f"[green]Created[/green] {p} and outputs/ layout")
    console.print("Next: edit the config, then run [bold]avrs doctor --config " + str(p) + "[/bold]")


# ------------------------------------------------------------------- doctor
@app.command()
def doctor(config: Optional[str] = typer.Option(None, "--config", "-c")):
    """Environment and dependency check (IGV, MCP server, bio tools, AI provider)."""
    cfg = _config(config) if config or os.environ.get(CONFIG_ENV) else AVRSConfig()
    t = Table(title="AVRS doctor", show_lines=False)
    t.add_column("Check")
    t.add_column("Status")
    t.add_column("Detail", overflow="fold")

    def row(name, ok, detail):
        t.add_row(name, "[green]OK[/green]" if ok else "[red]MISSING[/red]", detail)

    row("AVRS version", True, __version__)
    try:
        import pysam
        row("pysam", True, f"{pysam.__version__} (primary BAM/VCF backend)")
    except ImportError:
        row("pysam", False, "not importable -> pure-python evidence backend "
                            "(recommended: install in WSL2/Linux)")

    for tool, args in [("samtools", ["--version"]), ("bcftools", ["--version"]),
                       ("bwa-mem2", ["version"]), ("fastp", ["--version"]),
                       ("gatk", ["--version"]), ("mosdepth", ["--version"]),
                       ("snakemake", ["--version"]), ("node", ["--version"])]:
        exe = shutil.which(tool)
        if exe:
            try:
                out = subprocess.run([exe, *args], capture_output=True, text=True, timeout=30)
                first = (out.stdout or out.stderr).strip().splitlines()
                row(tool, True, f"{exe} - {first[0] if first else ''}")
            except Exception as e:
                row(tool, True, f"{exe} (version probe failed: {e})")
        else:
            row(tool, False, "not on PATH")

    # IGV
    from .igv.transport import DirectTCPTransport
    tr = DirectTCPTransport(cfg.igv.host, cfg.igv.port)
    if tr.connectivity_check():
        row("IGV Desktop port", True, f"{cfg.igv.host}:{cfg.igv.port} reachable")
    else:
        row("IGV Desktop port", False,
            f"{cfg.igv.host}:{cfg.igv.port} not reachable - start IGV Desktop and enable "
            "View > Preferences > Advanced > port command listener")
    if cfg.igv.transport == "mcp":
        if cfg.igv.mcp_command and shutil.which(cfg.igv.mcp_command):
            row("igv-mcp server", True, f"{cfg.igv.mcp_command} {cfg.igv.mcp_args}")
        else:
            row("igv-mcp server", False,
                "mcp_command missing/not on PATH - see docs/DEPLOYMENT.md for the official "
                "igv-mcp setup (https://github.com/igvteam/igv-mcp)")
    else:
        row("igv-mcp server", True, "not required (transport=direct uses the same TCP protocol)")

    # AI provider
    if cfg.ai.provider == "mock":
        row("AI provider", True, "mock (deterministic rule-based review)")
    elif cfg.ai.provider == "ollama":
        row("AI provider", True, f"ollama @ {cfg.ai.base_url or 'http://127.0.0.1:11434'} "
                                 f"model={cfg.ai.model} vision={cfg.ai.supports_vision}")
    else:
        key_ok = bool(cfg.ai.api_key_env and os.environ.get(cfg.ai.api_key_env))
        row("AI provider", key_ok, f"{cfg.ai.provider} @ {cfg.ai.base_url} "
                                   f"key_env={cfg.ai.api_key_env} "
                                   f"{'set' if key_ok else 'NOT SET'}")
    row("Remote screenshot transfer", cfg.ai.allow_remote_screenshot_transfer,
        "enabled" if cfg.ai.allow_remote_screenshot_transfer
        else "disabled (default: screenshots stay local)")
    console.print(t)
    console.print(f"\n[dim]{DISCLAIMER}[/dim]")


# ----------------------------------------------------------------- variants
@variants_app.command("import")
def variants_import(vcf: str = typer.Option(..., help="VCF path (.vcf or .vcf.gz)"),
                    sample: str = typer.Option(..., help="Sample id"),
                    config: Optional[str] = typer.Option(None, "--config", "-c"),
                    build: Optional[str] = typer.Option(None, help="Reference build override")):
    """Import a VCF into the database (allele-level, normalized)."""
    from .variants import import_vcf
    from .io.fasta import Fasta
    cfg = _config(config)
    engine, session = _db(cfg)
    fasta = Fasta(cfg.reference.fasta) if cfg.reference.fasta and \
        Path(cfg.reference.fasta).exists() else None
    smp = session.query(__import__("avrs.db.models", fromlist=["Sample"]).Sample) \
        .filter_by(sample_id=sample).one_or_none()
    if smp is None:
        from .db.models import Sample
        smp = Sample(sample_id=sample, reference_build=build or cfg.reference.build,
                     analysis_type=cfg.analysis.data_type)
        session.add(smp)
        session.commit()
    n, skipped = import_vcf(session, vcf, smp.id, sample, build or cfg.reference.build,
                            fasta=fasta)
    console.print(f"[green]Imported {n} variants[/green] (skipped {skipped}) for {sample}")


@variants_app.command("list")
def variants_list(sample: Optional[str] = typer.Option(None),
                  config: Optional[str] = typer.Option(None, "--config", "-c"),
                  limit: int = 50):
    """List imported variants."""
    from .db.models import Variant
    cfg = _config(config)
    engine, session = _db(cfg)
    q = session.query(Variant)
    if sample:
        from .db.models import Sample
        smp = session.query(Sample).filter_by(sample_id=sample).one()
        q = q.filter_by(sample_id=smp.id)
    t = Table(title="Variants")
    for c in ("variant_uid", "chrom", "pos", "ref", "alt", "type", "GT", "FILTER"):
        t.add_column(c)
    for v in q.limit(limit).all():
        t.add_row(v.variant_uid, v.chrom, str(v.pos), v.ref, v.alt, v.variant_type,
                  v.genotype or "-", v.original_filter)
    console.print(t)


# ------------------------------------------------------------------- triage
@app.command()
def triage(sample: str = typer.Option(...), bam: str = typer.Option(...),
           fasta: Optional[str] = typer.Option(None),
           config: Optional[str] = typer.Option(None, "--config", "-c")):
    """Compute BAM evidence + risk rules and create review tasks."""
    from .review import triage_sample
    cfg = _config(config)
    engine, session = _db(cfg)
    fasta_path = fasta or cfg.reference.fasta
    if not fasta_path:
        console.print("[red]--fasta or reference.fasta in config is required[/red]")
        raise typer.Exit(2)
    n = triage_sample(session, cfg, sample, fasta_path, bam)
    console.print(f"[green]Created {n} review tasks[/green] for {sample}")


# ---------------------------------------------------------------------- igv
@app.command("igv")
def igv_cmd(check: bool = typer.Option(False, help="Only run the connectivity check"),
            config: Optional[str] = typer.Option(None, "--config", "-c")):
    """IGV / igv-mcp connectivity check."""
    cfg = _config(config)
    from .igv.transport import DirectTCPTransport, MCPStdioTransport
    t = Table(title="IGV status")
    t.add_column("Item")
    t.add_column("Status")
    t.add_column("Detail", overflow="fold")
    tr = DirectTCPTransport(cfg.igv.host, cfg.igv.port)
    ok = tr.connectivity_check()
    t.add_row("command listener", "[green]OK[/green]" if ok else "[red]DOWN[/red]",
              f"{cfg.igv.host}:{cfg.igv.port}")
    if cfg.igv.transport == "mcp" and cfg.igv.mcp_command:
        try:
            m = MCPStdioTransport(cfg.igv.mcp_command, cfg.igv.mcp_args)
            m.start()
            tools = m.list_tools()
            t.add_row("igv-mcp tools", "[green]OK[/green]", f"{len(tools)} tools: "
                                                        + ", ".join(x['name'] for x in tools[:8]) + " ...")
            m.close()
        except Exception as e:
            t.add_row("igv-mcp tools", "[red]ERROR[/red]", str(e))
    else:
        t.add_row("igv-mcp tools", "[dim]n/a[/dim]", "transport=direct")
    console.print(t)
    if not ok:
        console.print("[yellow]Start IGV Desktop and enable the port command listener "
                      "(View > Preferences > Advanced > Port enabled, default 60151).[/yellow]")


# ------------------------------------------------------------------- review
@app.command()
def review(sample: Optional[str] = typer.Option(None, help="Sample id"),
           variant: Optional[str] = typer.Option(None, help="Single variant uid"),
           bam: Optional[str] = typer.Option(None), vcf: Optional[str] = typer.Option(None),
           fasta: Optional[str] = typer.Option(None),
           config: Optional[str] = typer.Option(None, "--config", "-c"),
           mock_igv: bool = typer.Option(False, help="Use the built-in mock IGV server"),
           limit: Optional[int] = typer.Option(None, help="Max review tasks")):
    """Run the automated IGV + AI review for a sample (or a single variant)."""
    from .review import run_reviews, triage_sample, import_sample_bam_vcf
    cfg = _config(config)
    engine, session = _db(cfg)
    if not sample:
        console.print("[red]--sample is required[/red]")
        raise typer.Exit(2)
    bam_path = bam or _find_file(session, sample, "bam")
    vcf_path = vcf or _find_file(session, sample, "vcf")
    fasta_path = fasta or cfg.reference.fasta
    if not (bam_path and vcf_path and fasta_path):
        console.print(f"[red]Missing inputs (bam={bam_path}, vcf={vcf_path}, "
                      f"fasta={fasta_path})[/red]")
        raise typer.Exit(2)
    import_sample_bam_vcf(session, cfg, sample, bam_path, vcf_path)
    if variant:
        # single variant: create/refresh its task then review just that one
        from .db.models import Variant, ReviewTask
        v = session.query(Variant).filter_by(variant_uid=variant).one()
        t = session.query(ReviewTask).filter_by(variant_id=v.id).one_or_none()
        if not t:
            t = ReviewTask(sample_id=v.sample_id, variant_id=v.id)
            session.add(t)
            session.commit()
        t.status = __import__("avrs.db.models", fromlist=["TaskState"]).TaskState.PENDING
        session.commit()
    else:
        triage_sample(session, cfg, sample, fasta_path, bam_path)
    stats = run_reviews(session, cfg, sample, bam_path, vcf_path, fasta_path,
                        use_mock_igv=mock_igv, limit=limit)
    console.print(json.dumps(stats, indent=2, default=str))


def _find_file(session, sample_id: str, kind: str) -> Optional[str]:
    from .db.models import Sample, SampleFile
    smp = session.query(Sample).filter_by(sample_id=sample_id).one_or_none()
    if not smp:
        return None
    f = session.query(SampleFile).filter_by(sample_id=smp.id, kind=kind).first()
    return f.path if f else None


# -------------------------------------------------------------- fn discovery
@app.command("fn-discovery")
def fn_discovery(sample: str = typer.Option(...), bam: Optional[str] = typer.Option(None),
                 vcf: Optional[str] = typer.Option(None),
                 raw_vcf: Optional[str] = typer.Option(None, help="Unfiltered/raw VCF"),
                 hotspots: Optional[str] = typer.Option(None, help="hotspot TSV/BED"),
                 fasta: Optional[str] = typer.Option(None),
                 config: Optional[str] = typer.Option(None, "--config", "-c")):
    """Discover false-negative CANDIDATES (filtered recovery + targeted loci)."""
    from .fn_discovery import recover_filtered_candidates, targeted_locus_review, \
        load_hotspot_file, summarize
    cfg = _config(config)
    engine, session = _db(cfg)
    bam_path = bam or _find_file(session, sample, "bam")
    vcf_path = vcf or _find_file(session, sample, "vcf")
    fasta_path = fasta or cfg.reference.fasta
    if not (bam_path and fasta_path):
        console.print("[red]--bam and fasta are required[/red]")
        raise typer.Exit(2)
    all_cands = []
    if raw_vcf:
        all_cands += recover_filtered_candidates(session, sample, cfg, raw_vcf,
                                                 bam_path, fasta_path)
    if hotspots:
        loci = load_hotspot_file(hotspots)
        all_cands += targeted_locus_review(session, sample, cfg, bam_path, fasta_path, loci)
    console.print(json.dumps(summarize(session, sample), indent=2, ensure_ascii=False))


# -------------------------------------------------------------------- report
@app.command()
def report(sample: Optional[str] = typer.Option(None),
           config: Optional[str] = typer.Option(None, "--config", "-c"),
           out: Optional[str] = typer.Option(None, help="Output directory")):
    """Export JSON/CSV/HTML reports + QC summary + audit log."""
    from .reports import generate_reports
    cfg = _config(config)
    engine, session = _db(cfg)
    out_dir = out or str(cfg.resolve_output("reports"))
    files = generate_reports(session, out_dir, sample_id=sample)
    for k, v in files.items():
        console.print(f"[green]{k}[/green]: {v}")


# --------------------------------------------------------------------- qc
@app.command()
def qc(sample: str = typer.Option(...),
       config: Optional[str] = typer.Option(None, "--config", "-c")):
    """Show QC metrics recorded for a sample."""
    from .db.models import QCMetric, Sample
    cfg = _config(config)
    engine, session = _db(cfg)
    smp = session.query(Sample).filter_by(sample_id=sample).one()
    rows = session.query(QCMetric).filter_by(sample_id=smp.id).all()
    t = Table(title=f"QC metrics - {sample}")
    t.add_column("source")
    t.add_column("key")
    t.add_column("value", overflow="fold")
    for m in rows:
        t.add_row(m.source, m.key, m.value_json)
    console.print(t)


# ------------------------------------------------------------------ pipeline
@pipeline_app.command("run")
def pipeline_run(config: str = typer.Option(..., help="Pipeline config yaml"),
                 sample: Optional[str] = typer.Option(None),
                 dry_run: bool = typer.Option(False),
                 cores: int = 4):
    """Run the Snakemake FASTQ->VCF workflow (requires snakemake + bio tools)."""
    from .pipeline.runner import run_pipeline
    run_pipeline(config, sample=sample, dry_run=dry_run, cores=cores)


# --------------------------------------------------------------------- serve
@app.command()
def serve(config: Optional[str] = typer.Option(None, "--config", "-c"),
          host: str = "127.0.0.1", port: int = 8000, reload: bool = False):
    """Start the FastAPI review workbench."""
    cfg = _config(config)
    import uvicorn
    from .api.app import create_app
    console.print(f"[green]AVRS workbench:[/green] http://{host}:{port}")
    uvicorn.run(create_app(config or os.environ.get(CONFIG_ENV)),
                host=host, port=port, reload=reload, factory=False)


# ---------------------------------------------------------------------- demo
@app.command()
def demo(config: Optional[str] = typer.Option(None, "--config", "-c"),
         out: Optional[str] = typer.Option(None, help="Demo output dir"),
         mock_igv: bool = typer.Option(True, help="Use the mock IGV server"),
         with_fake_igv: bool = typer.Option(False, "--real-igv",
                                            help="Use the real IGV Desktop instead of the mock")):
    """Generate the demo dataset and run the full Mode-B loop end to end."""
    from .demo.generator import generate_demo
    from .review import (import_sample_bam_vcf, triage_sample, run_reviews)
    from .reports import generate_reports
    from .variants import import_vcf
    from .io.fasta import Fasta
    cfg = _config(config)
    out_dir = Path(out or cfg.resolve_output("demo"))
    console.print(f"[bold]1/5[/bold] generating demo dataset in {out_dir}")
    manifest = generate_demo(out_dir)
    engine, session = _db(cfg)
    sample_id = manifest["sample_id"]
    console.print(f"[bold]2/5[/bold] importing VCF")
    smp = import_sample_bam_vcf(session, cfg, sample_id, manifest["bam"], manifest["vcf"])
    fasta = Fasta(manifest["fasta"])
    import_vcf(session, manifest["vcf"], smp.id, sample_id,
               manifest["reference_build"], fasta=fasta)
    session.commit()
    console.print(f"[bold]3/5[/bold] extracting evidence + triage")
    n_tasks = triage_sample(session, cfg, sample_id, manifest["fasta"], manifest["bam"])
    console.print(f"      {n_tasks} review tasks")
    console.print(f"[bold]4/5[/bold] running IGV + AI review "
                  f"({'mock IGV' if mock_igv and not with_fake_igv else 'real IGV'})")
    stats = run_reviews(session, cfg, sample_id, manifest["bam"], manifest["vcf"],
                        manifest["fasta"], use_mock_igv=(mock_igv and not with_fake_igv))
    console.print(f"      {json.dumps(stats, default=str)}")
    console.print(f"[bold]5/5[/bold] exporting reports")
    files = generate_reports(session, cfg.resolve_output("reports"),
                             sample_id=sample_id)
    console.print(json.dumps(files, indent=2))
    console.print(f"\n[green]Demo complete.[/green] Start the workbench with: "
                  f"avrs serve --config {config or os.environ.get(CONFIG_ENV)}")


@app.command()
def version():
    """Print version and disclaimer."""
    console.print(f"AVRS {__version__}")
    console.print(f"[dim]{DISCLAIMER}[/dim]")


if __name__ == "__main__":
    app()
