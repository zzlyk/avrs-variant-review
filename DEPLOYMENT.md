# DEPLOYMENT

Supported runtime: **WSL2 Ubuntu** (analysis) + **Windows** (IGV Desktop).
A pure Linux deployment works identically.

---

## 1. Architecture of a deployment

```
┌──────────────────────────── Windows host ─────────────────────────────┐
│  IGV Desktop (port listener 60151)        AVRS workbench (uvicorn)    │
│        ▲                                        ▲                    │
│        │ TCP (localhost only)                  │ HTTP (localhost)     │
└────────┼────────────────────────────────────────┼────────────────────┘
         │                                        │
┌────────▼──────── WSL2 Ubuntu ──────────────────▼────────────────────┐
│  AVRS CLI / FastAPI / review engine                                  │
│  Snakemake workflow (fastp · bwa-mem2 · samtools · GATK · bcftools)  │
│  SQLite DB + snapshots on /mnt/f (or a native ext4 volume)           │
│  Optional: Ollama (local vision model)                              │
└──────────────────────────────────────────────────────────────────────┘
```

* The analysis never needs the Windows filesystem, but sharing `/mnt/f` keeps
  BAM/VCF and IGV on the same volume.
* For speed, prefer a native WSL2 ext4 volume (`\\wsl$\...` or a VHD) for
  the working directory; `/mnt/f` is fine for reference data.

---

## 2. WSL2 installation (no root required)

```bash
# from WSL
cd /mnt/f/codex/mcp-igv
bash scripts/wsl_setup.sh            # micromamba + bioconda toolchain
source ~/micromamba/etc/profile.d/micromamba.sh
micromamba activate avrs
pip install -e .
avrs doctor                          # environment check
```

`scripts/wsl_setup.sh` installs into `~/micromamba/envs/avrs`:

| Tool | Version installed | Purpose |
|---|---|---|
| Python | 3.11 | runtime |
| pysam | 0.24.1 | primary BAM/VCF evidence backend |
| samtools | 1.24 | BAM QC, flagstat/stats, indexing |
| bcftools | 1.24 | VCF normalisation/filtering |
| bwa-mem2 | 2.3 | alignment |
| fastp | 1.4.0 | FASTQ QC + trimming |
| mosdepth | 0.3.14 | coverage |
| gatk4 | 4.6.2.0 (OpenJDK 17) | HaplotypeCaller / GenotypeGVCFs |
| snakemake | 9.27 | workflow engine |

With root access you may instead use `apt`/`mamba`; the only requirement is
that the tools above are on `PATH`.

### WSL2 specifics

* Enable systemd in WSL if you want services; otherwise run the workbench
  under `nohup`/`tmux`.
* `/mnt/f` is slower than the native filesystem — put the SQLite DB and the
  snapshot directory on a native volume if throughput matters.
* Memory: GATK needs `--java-options '-Xmx4g'` (already set in the workflow).

---

## 3. IGV Desktop (Windows)

1. Install IGV Desktop (https://igv.org/download.html).
2. Enable the command listener: `View ▸ Preferences ▸ Advanced ▸ Port
   enabled` (default port **60151**).
3. Verify from WSL:

```bash
python - <<'PY'
import socket
s = socket.create_connection(("127.0.0.1", 60151), timeout=3)
print("IGV port reachable from WSL")
s.close()
PY
```

Windows Firewall may block WSL→Windows connections; if the check fails,
allow the port for the local profile or use
`netsh interface portproxy`.

### Official MCP server (optional)

AVRS can drive IGV either directly (same TCP protocol) or through the
official `igv-mcp` server:

```bash
git clone https://github.com/igvteam/igv-mcp /opt/igv-mcp
cd /opt/igv-mcp && npm install
node src/main.js --host 127.0.0.1:60151      # stdio MCP server
```

```yaml
igv:
  transport: mcp
  mcp_command: node
  mcp_args: ["/opt/igv-mcp/src/main.js", "port", "60151"]
```

`avrs igv check` lists the tools the server actually exposes.

---

## 4. Linux (native) deployment

```bash
sudo apt install -y python3.11 python3-pip python3-venv samtools bcftools \
                    fastp openjdk-17-jre-headless
python3.11 -m venv ~/venvs/avrs && source ~/venvs/avrs/bin/activate
pip install -e .
# bwa-mem2, mosdepth, gatk4 and snakemake: install via bioconda/conda
avrs doctor
```

Everything else is identical; set `igv.path_mapping_enabled: false` unless
IGV runs on a different filesystem.

---

## 5. Configuration

```yaml
project_name: lab-demo
reference:
  build: GRCh38
  fasta: /data/reference/GRCh38.fa
paths:
  output_dir: /data/avrs/out
  snapshot_dir: /data/avrs/out/snapshots
  db: /data/avrs/out/avrs.db
  path_mapping:                       # only when IGV is on another filesystem
    - from: /mnt/f/lab
      to: F:/lab
analysis:
  mode: review_only
  sequencing: germline
  data_type: panel
  target_bed: /data/targets/panel.bed
triage:
  profiles:
    panel: { min_dp: 30, min_gq: 40, min_alt_ad: 8, ... }   # lab-specific
review:
  modes: [suspicious_only]
  max_retries: 3
igv:
  transport: direct
  host: 127.0.0.1
  port: 60151
  path_mapping_enabled: true
ai:
  provider: ollama
  model: qwen2.5vl:7b
  base_url: http://127.0.0.1:11434
  supports_vision: true
security:
  local_first: true
```

`avrs init` writes a starter config; `avrs doctor` validates it and the
environment.

---

## 6. Running as a service

### Workbench (FastAPI/uvicorn)

```bash
cat > ~/avrs.service <<'EOF'
[Unit]
Description=AVRS review workbench
After=network.target

[Service]
User=%u
WorkingDirectory=/mnt/f/lab/mcp-igv
Environment=PATH=%h/micromamba/envs/avrs/bin:/usr/bin:/bin
Environment=AVRS_CONFIG=/mnt/f/lab/mcp-igv/config/lab.yaml
ExecStart=%h/micromamba/envs/avrs/bin/uvicorn avrs.api.app:app --host 127.0.0.1 --port 8000
Restart=on-failure

[Install]
WantedBy=default.target
EOF
```

> `avrs.api.app:app` is a placeholder; use the CLI instead:
> `ExecStart=%h/micromamba/envs/avrs/bin/avrs serve --config /mnt/f/lab/mcp-igv/config/lab.yaml`

### Scheduled reviews

Use cron or the Snakemake scheduler to run `avrs review` per sample; the
task state machine makes re-runs safe (completed tasks are skipped, failed
tasks are retried up to `review.max_retries`).

---

## 7. Operational runbook

| Situation | Action |
|---|---|
| IGV not reachable | `avrs igv check` → start IGV, enable the port listener, check the firewall |
| A sample's reviews all failed | `avrs review --sample SAMPLE --config …` re-runs failed tasks; inspect `review_tasks.last_error` in the DB |
| Stale IGV lock | delete `outputs/snapshots/.igv.lock` after confirming no review is running |
| Model unavailable | the system still works with `ai.provider: mock`; switch back when the endpoint is healthy |
| Disk filling up | snapshots dominate; archive `outputs/snapshots/<sample>` after the review is signed off |
| New reference build | change `reference.build`; variant IDs embed the build, so IDs never collide across builds |
| Want to prove nothing changed | `variant_evidence.bam_checksum` / `igv_snapshots.bam_checksum` record the exact inputs |

---

## 8. Backup & retention

* **Database** (`avrs.db`) — the review record. Back up with the snapshot
  directory; the two must be restored together.
* **Snapshots** — visual evidence; keep as long as the review record is
  legally required.
* **Raw data** — never modified by AVRS; managed by the LIMS/backup policy.
* **Reports** — regenerate at any time with `avrs report`.

---

## 9. Upgrading

```bash
git pull
source ~/micromamba/etc/profile.d/micromamba.sh && micromamba activate avrs
pip install -e .
python -m pytest tests/ -q
avrs doctor --config config/lab.yaml
```

Schema changes are additive; `avrs.db` is created with
`create_all()` and existing tables are left untouched. For destructive
schema changes, export the reports first and recreate the database.
