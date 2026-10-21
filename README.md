# AVRS — AI Variant Review System

> **Research / validation use only, until appropriately validated and approved.**
> Not for clinical diagnostic use. See [SECURITY.md](SECURITY.md) and [VALIDATION.md](VALIDATION.md).

[![CI](https://github.com/zzlyk/avrs-variant-review/actions/workflows/ci.yml/badge.svg)](https://github.com/zzlyk/avrs-variant-review/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

AVRS automates the repetitive part of germline variant review: it extracts
quantitative evidence from BAM/VCF, applies configurable risk rules, drives
**IGV Desktop** through the official `igv-mcp` interface, captures
screenshots, asks a (vision-capable) model for a **structured recommendation**,
and presents everything in a human review workbench with a full audit trail.

The guiding principle:

> **Not "let an AI look at IGV screenshots" — but "sequencing evidence +
> IGV visualisation + AI-assisted interpretation + human final sign-off".**

---

## 1. What the system does

| Stage | Module | Output |
|---|---|---|
| FASTQ → VCF (Mode A) | `workflow/Snakefile`, `avrs.pipeline` | coordinate-sorted BAM, bgzipped/indexed VCF, QC |
| VCF import (Mode B) | `avrs.variants` | allele-level, normalised variants with stable build-scoped IDs |
| Evidence extraction | `avrs.evidence` | per-variant BAM metrics (SNV **and** indel-aware) |
| Risk triage | `avrs.triage` | 20+ configurable rules → HIGH/MEDIUM/LOW |
| IGV automation | `avrs.igv` | narrow + context screenshots, bound to the exact variant |
| AI review | `avrs.ai` | schema-validated recommendation (mock / OpenAI-compatible / Ollama) |
| Human review | `avrs.api` + `frontend/` | workbench with accept/reject/escalate and decision history |
| FN discovery | `avrs.fn_discovery` | filtered-recovery + targeted-locus **candidates** |
| Reporting | `avrs.reports` | JSON / CSV / HTML + QC summary + exceptions + audit log |

Two run modes:

* **Mode A — Full pipeline**: FASTQ → QC → alignment → BAM → VCF → review.
* **Mode B — Review only**: existing BAM + BAI + VCF + reference → review.
  *(This is the fastest path to value and is fully implemented.)*

---

## 2. Architecture

```
                     ┌──────────────────────────────────────────┐
  FASTQ ──────────▶ │ Snakemake workflow (Mode A)              │
  (Mode A)          │ fastp → bwa-mem2 → samtools → GATK →     │
                     │ bcftools norm/filter                     │
                     └───────────────┬──────────────────────────┘
                                     │ BAM + BAI + VCF
                     ┌───────────────▼──────────────────────────┐
                     │ avrs.variants   import / normalise / ID │
                     │ avrs.evidence   allele-aware BAM metrics│
                     │ avrs.triage     configurable risk rules │
                     └───────────────┬──────────────────────────┘
                                     │ review tasks
             ┌───────────────────────▼──────────────────────────┐
             │ avrs.igv   IGV MCP (official igv-mcp protocol)  │
             │  lock → load → goto → render-wait → snapshot×2  │
             └───────────────────────┬──────────────────────────┘
                                     │ PNG + metadata binding
             ┌───────────────────────▼──────────────────────────┐
             │ avrs.ai  L1 evidence → L2 vision → L3 fusion     │
             │          strict Pydantic AIReviewResult         │
             └───────────────────────┬──────────────────────────┘
                                     │
             ┌───────────────────────▼──────────────────────────┐
             │ avrs.api + frontend/  human review workbench    │
             │ SQLite (13 tables) + audit log + reports         │
             └──────────────────────────────────────────────────┘
```

Details: [ARCHITECTURE.md](ARCHITECTURE.md).

---

## 3. Requirements

**Analysis side (WSL2 Ubuntu / Linux — the supported runtime):**

```bash
# one-shot, no root required (micromamba + bioconda)
bash scripts/wsl_setup.sh
```

This installs Python 3.11, `pysam`, `samtools`, `bcftools`, `bwa-mem2`, `fastp`,
`mosdepth`, `gatk4` (OpenJDK 17) and `snakemake` into `~/micromamba/envs/avrs`.

**IGV side (Windows desktop):**

* IGV Desktop with the *port command listener* enabled
  (`View ▸ Preferences ▸ Advanced ▸ Port enabled`, default `60151`).
* Optional: the official MCP server <https://github.com/igvteam/igv-mcp>
  (AVRS can talk to IGV either through that MCP server **or** directly over
  the identical TCP command protocol).

`avrs doctor` checks every one of these and prints a table.

---

## 4. Install

```bash
git clone <this repo> && cd mcp-igv
cd /mnt/f/codex/mcp-igv                       # inside WSL2
bash scripts/wsl_setup.sh                    # toolchain (once)
source ~/micromamba/etc/profile.d/micromamba.sh
micromamba activate avrs
pip install -e .                             # AVRS itself
```

---

## 5. Quick start — the demo (works fully offline)

```bash
avrs demo --config config/demo.yaml
```

The demo generates a synthetic 2-contig genome with **7 review scenarios**,
then runs the complete Mode-B loop (import → evidence → triage → IGV → AI →
report) against a built-in mock IGV listener, and writes everything to
`outputs/`.

Expected scenario behaviour (mock IGV, mock reviewer):

| Scenario | Designed signal | Typical outcome |
|---|---|---|
| clean SNV | balanced ALT support, high BQ/MQ | no rule hit → not queued |
| low-quality SNV | all reads below MAPQ floor | `INSUFFICIENT_EVIDENCE` |
| strand-biased SNV | all ALT reads on one strand | `SUSPECTED_FALSE_POSITIVE` |
| homopolymer indel | insertion inside an 11 bp poly-A | `INCONCLUSIVE` (representation risk) |
| low-coverage locus | DP below floor | `INSUFFICIENT_EVIDENCE` |
| caller miss | ALT support but absent from VCF | FN **candidate** (`possible_missed_variant`) |
| equivocal | caller AD ≠ pileup AD | `INCONCLUSIVE` |

> Demo labels are *scenario* labels for pipeline testing. They are **not**
> clinical truth and must never be used as a truth set.

---

## 6. Running against your own data (Mode B)

```bash
# 1) configure
cp config/config.example.yaml config/my_lab.yaml
$EDITOR config/my_lab.yaml          # reference, thresholds, IGV host, AI provider

# 2) check the environment
avrs doctor --config config/my_lab.yaml

# 3) import the VCF (allele-level, normalised, stable IDs)
avrs variants import --vcf data/SAMPLE001.vcf.gz --sample SAMPLE001 \
    --config config/my_lab.yaml

# 4) evidence + triage
avrs triage --sample SAMPLE001 --bam data/SAMPLE001.bam \
    --fasta reference/GRCh38.fa --config config/my_lab.yaml

# 5) automated IGV + AI review
avrs review --sample SAMPLE001 --config config/my_lab.yaml

# 6) false-negative candidate discovery
avrs fn-discovery --sample SAMPLE001 --hotspots data/hotspots.tsv \
    --config config/my_lab.yaml

# 7) reports
avrs report --sample SAMPLE001 --config config/my_lab.yaml

# 8) workbench
avrs serve --config config/my_lab.yaml     # http://127.0.0.1:8000
```

### Mode A (FASTQ → VCF)

```bash
cp config/pipeline.example.yaml config/pipeline.yaml
$EDITOR config/pipeline.yaml            # fastq paths, reference, filters
avrs pipeline run --config config/pipeline.yaml
```

The workflow runs fastp → bwa-mem2 → samtools sort/markdup/index → GATK
HaplotypeCaller (gVCF) → GenotypeGVCFs → configurable hard filters (VQSR only
when truth resources are configured) → `bcftools norm` → tabix, plus
flagstat/stats/mosdepth/MultiQC. The original FASTQ is never modified.

---

## 7. WSL2 + Windows IGV (path mapping)

The analysis runs in WSL2 (`/mnt/f/...`), IGV Desktop runs on Windows
(`F:\...`). Configure the translation explicitly:

```yaml
paths:
  path_mapping:
    - from: /mnt/f/codex/mcp-igv
      to: F:/codex/mcp-igv
igv:
  path_mapping_enabled: true      # only true when IGV is on another filesystem
  host: 127.0.0.1
  port: 60151
```

When IGV and the analysis share a filesystem (e.g. an IGV inside WSL, or the
built-in mock), leave `path_mapping_enabled: false` so native paths are used.

---

## 8. AI providers

```yaml
ai:
  provider: mock                 # deterministic, offline, no key required
  model: mock-reviewer-v1
  # provider: ollama             # local models (llava, qwen2-vl, gpt-oss ...)
  # base_url: http://127.0.0.1:11434
  # supports_vision: true        # screenshots stay local with a local provider
  # provider: openai_compatible  # any OpenAI-compatible endpoint
  # base_url: https://api.example.com/v1
  # api_key_env: OPENAI_API_KEY  # the key is read from the environment only
  allow_remote_screenshot_transfer: false   # explicit opt-in for cloud vision
```

Guarantees enforced in code:

* Without a key/model service the system still runs the deterministic
  rule-based review (`provider: mock`).
* If the model has no vision support, **no image is sent** and the result is
  marked `visual_review_performed=false` — visual review is never faked.
* Screenshot transfer to a remote provider requires
  `allow_remote_screenshot_transfer: true`.
* VCF/annotation text is passed as *untrusted data*; the model never gets
  shell or file-system access.

---

## 9. Outputs

```
outputs/
├── avrs.db                     # SQLite (13 tables)
├── demo/                       # generated demo dataset
├── snapshots/<sample>/<variant_uid>/{narrow,context}.png
└── reports/
    ├── review_results_<sample>_<ts>.json / .csv
    ├── review_summary_<sample>_<ts>.json / .html
    ├── exception_report_<sample>_<ts>.json
    ├── qc_summary_<sample>_<ts>.json
    └── audit_log_<sample>_<ts>.json
```

Every snapshot is bound to `sample_id`, `variant_uid`, `genome_build`, the
actual locus, the BAM/VCF checksums, the IGV session and a timestamp — a
screenshot can never be attributed to the neighbouring variant.

---

## 10. Tests & validation

```bash
python -m pytest tests/ -q            # unit + integration (42+ tests)
bash scripts/wsl_pipeline_smoke.sh    # real FASTQ -> VCF -> review smoke test
```

The suite covers FASTQ/manifest validation, reference consistency, VCF parsing,
allele-level normalisation, SNV/indel evidence extraction, CIGAR edge cases,
empty coverage, variant-ID stability, the rule engine, AI JSON validation,
review state transitions, mock-IGV load/goto/snapshot/retry/locking, snapshot
locus association and the full Mode-B end-to-end loop.

See [VALIDATION.md](VALIDATION.md) for what is **and is not** validated.

---

## 11. Documentation

| File | Content |
|---|---|
| [ARCHITECTURE.md](ARCHITECTURE.md) | modules, data flow, data model, extension points |
| [VALIDATION.md](VALIDATION.md) | test matrix, benchmark adapters (GIAB/hap.py), acceptance criteria |
| [SECURITY.md](SECURITY.md) | data protection, model access, threat model, clinical limits |
| [DEPLOYMENT.md](DEPLOYMENT.md) | WSL2/Linux/Windows setup, IGV MCP, ops runbook |

---

## 12. Limitations (read before use)

* **Not validated for clinical use.** No HIPAA/CLIA/CAP compliance is claimed.
* AI labels are **recommendations**, not truth-validated TP/FP/FN calls.
* The `risk_score` is a **deterministic rule score**, not a calibrated
  probability; the schema carries `risk_score_is_calibrated_probability=false`.
* “Time saved” is only reported when a measured human baseline exists.
* Second-caller comparison is a pluggable interface — no DeepVariant run is
  wired in yet.
* Somatic/trio/CNV/SV analysis types are **not** implemented; the config
  separates germline from somatic, but only germline logic exists today.

---

## 13. License & citation

MIT — see [LICENSE](LICENSE).

If you use AVRS in research, please cite the software and state clearly that
it was **not clinically validated**:

```bibtex
@software{avrs2026,
  title  = {AVRS: AI-Assisted Variant Review System},
  author = {Yanke (zzlyk)},
  year   = {2026},
  url    = {https://github.com/zzlyk/avrs-variant-review},
  note   = {Research / validation use only; not validated for clinical diagnostic use.}
}
```

## 14. Acknowledgements

* The IGV integration follows the official
  [igv-mcp](https://github.com/igvteam/igv-mcp) server and its documented
  command protocol.
* The demo dataset is fully synthetic and contains no patient data.
