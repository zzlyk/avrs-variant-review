# ARCHITECTURE

## 1. Design goals

1. **Evidence first.** Quantitative BAM/VCF evidence is the foundation; the
   visual layer only adds interpretation and anomaly detection.
2. **Configurable, not hardcoded.** Every threshold, path, tool and model is
   configuration. The code contains no universal DP/MQ/VAF "standard".
3. **Separation of concerns.** Variant *authenticity* QC is strictly
   separated from variant *pathogenicity* annotation. AVRS never states that
   a variant is pathogenic.
4. **Auditability.** Raw calls, AI conclusions and human conclusions live in
   separate tables; nothing is overwritten; every decision is versioned.
5. **Fail loudly.** A missing screenshot, a failed model call or a failed
   task is recorded as such — never silently faked or downgraded.

## 2. Module map

```
src/avrs/
├── config.py            Pydantic configuration (typed + friendly errors)
├── cli.py               Typer CLI (same service layer as the web UI)
├── io/                  Dependency-light readers/writers
│   ├── bgzf.py          BGZF blocks + virtual offsets
│   ├── bam.py           BAM/BAI read+write, CIGAR utilities
│   ├── vcf.py           VCF parsing, multiallelic split, normalisation
│   └── fasta.py         FASTA reader/writer (+ .fai)
├── variants/            VCF import → allele-level normalised variants
├── evidence/            Backend abstraction + allele-aware evidence engine
│   ├── core.py          ReadView adapter + shared computation core
│   ├── metrics.py       VariantEvidence schema
│   └── __init__.py      pysam backend (primary) / pure-python (fallback)
├── triage/              20+ configurable rules, prioritisation, review modes
├── igv/                 IGV automation
│   ├── transport.py     DirectTCP + MCP-stdio transports (official igv-mcp)
│   ├── client.py        IGVMCPClient (app-level ops → real tools)
│   ├── reviewer.py      Per-variant review executor + snapshot binding
│   ├── mock_server.py   Mock IGV listener for tests / offline demo
│   └── client lock      IGVInstanceLock (single-instance serialisation)
├── ai/                  Review engine
│   ├── schema.py        Strict AIReviewResult (extra="forbid")
│   ├── prompts.py       Versioned prompt (v1.2) with hard rules
│   ├── providers.py     mock / OpenAI-compatible / Ollama
│   └── engine.py        L1 evidence → L2 vision → L3 fusion
├── fn_discovery/        Filtered recovery + targeted locus review
├── review/              Orchestration (triage, run, human decisions)
├── db/                  SQLAlchemy models (13 tables) + sessions
├── api/                 FastAPI backend
├── reports/             JSON/CSV/HTML + QC + exceptions + audit
├── pipeline/            FASTQ validation + Snakemake runner
└── demo/                Reproducible demo dataset generator
```

## 3. Data flow (Mode B)

```
VCF ──► variants.import_vcf
          │  split multiallelic → normalise (left-align) → stable variant_uid
          │  {build}#{chrom}:{pos}:{ref}>{alt}
          ▼
       variants table ──► evidence backend (pysam | pure-python)
          │                  • allele-aware read reconstruction (CIGAR)
          │                  • ref/alt/other counts, strand, MAPQ/BQ,
          │                    duplication, soft-clip, start clustering,
          │                    pair concordance, local depth, homopolymer
          ▼
       variant_evidence + rule hits + priority
          │
          ▼
       review_tasks ──► IGVReviewer (lock → load → goto → render-wait → ×2 PNG)
          │                         │
          │                         ▼
          │                   igv_snapshots (bound to variant + checksums)
          ▼
       AI ReviewEngine ──► ai_reviews (AIReviewResult JSON + provenance)
          │
          ▼
       human_reviews (append-only, supersede chain) + audit_logs
          │
          ▼
       reports (JSON/CSV/HTML) + workbench (FastAPI + static UI)
```

## 4. Evidence engine — why it is structured this way

`compute_evidence()` is the only implementation of the metric logic. Both
backends convert their native read objects into a `ReadView`
(name/flag/pos/mapq/cigar/seq/qual) and call the same function. This keeps
the correctness-critical code testable without a BAM file and identical on
Windows and Linux.

Key properties:

* **Allele reconstruction instead of base pileup.** For each read the engine
  walks the CIGAR and rebuilds the read's allele over the variant's
  reference span (insertions contribute their inserted bases, deletions
  contribute nothing, `N` makes the read uninformative). A read supports
  ALT only if its reconstructed allele equals the ALT allele — which is how
  indels are handled correctly instead of by single-base pileup.
* **Definitions are never conflated**: `total_reads` (raw overlap),
  `usable_reads` (after MAPQ/dup/secondary filters), `ref_count` /
  `alt_count` / `other_count` (allele classes), and `caller_dp` /
  `caller_ad_*` (what the caller reported) are separate fields.
* **Quality is reported, not hidden.** Low base quality is *counted*
  (`alt_reads_below_baseq`) instead of silently dropping the reads, so the
  "low-quality artifact" signal survives into the rule engine.
* **VCF-vs-pileup cross-check.** `ad_consistency` compares the caller's AD
  with the recomputed pileup — a mismatch is a high-severity signal.

## 5. Triage rules

Rules live in `avrs/triage/__init__.py` and are driven by a `TriageProfile`
(per data type: WGS / WES / Panel). Each hit carries:

```python
RuleHit(rule_id, observed_value, configured_threshold, reason, severity)
```

Rule ids cover: `LOW_DP`, `LOW_GQ`, `LOW_MQ`, `ALL_READS_LOW_MQ`, `LOW_ALT_AD`,
`STRAND_BIAS`, `LOW_ALT_BQ`, `ALT_READS_LOW_BQ`, `ALT_NEAR_READ_END`,
`SOFT_CLIPPING`, `START_POS_CLUSTER`, `HIGH_DUP_RATE`, `INDEL_IN_REPEAT`,
`LOW_COMPLEXITY`, `MULTI_MAPPING_RISK`, `AD_PILEUP_MISMATCH`,
`VAF_GT_MISMATCH`, `NO_COVERAGE`, `LOCAL_COVERAGE_ANOMALY`,
`PAIR_SUPPORT_INCONSISTENT`, `TARGET_REGION_EDGE`, `LOW_QUAL`.

Rules never assert truth — they flag signals. Priority is derived from rule
severity (any hit ⇒ at least MEDIUM; a `high` severity hit ⇒ HIGH). A
variant's clinical importance (e.g. ClinVar) never *lowers* its evidence
review requirements.

## 6. IGV integration

The official `igv-mcp` server was inspected directly (source of
`src/tools.yaml.js`, `src/mcpServer.js`, `src/igvConnection.js`) to derive
the real tool list and wire protocol:

* transport: TCP to IGV Desktop's command listener (default `60151`),
  one command per connection, `command arg…\n`, single-line reply;
* the MCP server wraps exactly that protocol over stdio JSON-RPC.

AVRS therefore ships two transports behind one interface:

| Transport | Use case |
|---|---|
| `direct` | IGV reachable directly; identical wire protocol, no Node dependency |
| `mcp` | spawn the official `igv-mcp` server and speak JSON-RPC to it |

Only tools that exist in the official tool set are used
(`new, saveSession, genome, load, goto, collapse, squish, expand, region,
setColor, group, sort, viewAsPairs, snapshot, snapshotDirectory, maxPanelHeight,
setSequenceStrand, setSequenceShowTranslation, overlay, separate, zoomin,
zoomout`). Operations with no backing tool (e.g. querying the IGV version)
raise `UnsupportedOperation` instead of pretending to work.

Per-variant review sequence (spec §5.2): acquire instance lock → health
check → load genome/BAM/VCF (once per sample) → goto with type-aware window
→ track configuration (squish, pair view, sort by base) → **render wait**
(file existence + minimum size + two-poll stability, never a fixed sleep as
the only criterion) → narrow + context snapshots → metadata binding →
release lock.

## 7. AI review engine

Three layers, fused:

1. **L1 deterministic evidence** — BAM/VCF metrics (authoritative).
2. **L2 IGV visual analysis** — only when a vision-capable model actually
   received the image. If the model has no vision support, no image is sent
   and `visual_review_performed=false` is recorded.
3. **L3 fusion** — the model output is validated against the strict
   `AIReviewResult` schema (`extra="forbid"`), the deterministic rule score
   is computed (`rulescore-v1`, explicitly *not* a calibrated probability),
   and the recommendation is combined with rule hits.

Provider abstraction: `mock` (deterministic, offline), `openai_compatible`
(any OpenAI-compatible endpoint), `ollama` (local). The prompt is versioned
(`v1.2`) and constrains the model: only use supplied evidence, separate
observation from inference, never invent read counts, never claim pathogenicity,
return `INCONCLUSIVE` when unsure. VCF/annotation text is passed as untrusted
data; the model gets no shell or filesystem access.

## 8. Data model (SQLite, 13 tables)

| Table | Purpose |
|---|---|
| `samples` | sample identity, build, analysis type, status |
| `sample_files` | file inventory + checksums (BAM/VCF/FASTQ/reference) |
| `pipeline_runs` | Snakemake/AVRS run state, logs, software versions |
| `qc_metrics` | fastp / flagstat / mosdepth / evidence QC |
| `variants` | raw + normalised representation, stable `variant_uid` |
| `variant_evidence` | computed metrics + rule hits + priority |
| `review_tasks` | per-variant work item with validated state machine |
| `igv_sessions` | IGV session bookkeeping |
| `igv_snapshots` | PNG path + variant/build/locus/checksum binding |
| `ai_reviews` | AI result JSON + model/prompt provenance |
| `human_reviews` | human decisions (append-only, supersede chain) |
| `audit_logs` | every state-changing action |
| `model_configurations` | model/provider registry |
| `fn_discovery_candidates` | FN candidates with explicit classification |

Task state machine: `PENDING → RUNNING → COMPLETED | FAILED | RETRYING |
NEEDS_MANUAL_ACTION`, transitions validated by `check_transition()`.

## 9. Scale and performance

* WGS produces millions of calls: variants are read from an indexed VCF,
  BAM access is region-scoped through the BAM index, and the database is
  written incrementally.
* The LLM never sees the whole VCF: only one variant plus its filtered,
  de-identified evidence block per request (`max_llm_context_variants`).
* Reviews run in bounded batches (`review.batch_size`); a failing variant is
  retried or marked `FAILED` without aborting the sample.
* IGV is a single GUI instance — access is serialised by a file lock
  (`IGVInstanceLock`) so two tasks can never interleave screenshots.

## 10. Extension points

| Want to add | How |
|---|---|
| A new AI provider | implement `AIProvider` in `avrs/ai/providers.py`, register in `make_provider` |
| A new evidence backend | implement the `Backend` protocol, return `ReadView`s, call `compute_evidence` |
| A new triage rule | add an `add(...)` call in `avrs/triage/__init__.py`, expose thresholds in `TriageProfile` |
| A second variant caller | subclass `SecondCallerAdapter` (fn_discovery) and register a caller name |
| A new report | extend `avrs/reports/__init__.py` |
| A new API consumer | the FastAPI app is a thin layer over `avrs.review` / `avrs.reports` |
