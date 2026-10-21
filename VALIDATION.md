# VALIDATION

This document states exactly what **is** and **is not** validated in AVRS.
It is deliberately explicit: an automated review assistant that overstates its
own validation status is more dangerous than one that admits its limits.

> **Research / validation use only, until appropriately validated and approved.**

---

## 1. Validation levels used here

| Level | Meaning |
|---|---|
| ✅ **Implemented + tested** | Code exists, is exercised by the automated test suite, and the test passes in the target environment. |
| 🟡 **Implemented, not validated against external truth** | Code works, but no independent truth set has been applied. |
| 🔵 **Interface only** | Pluggable contract exists; no concrete implementation shipped. |
| ⬜ **Not implemented** | Explicitly out of scope for this release. |

---

## 2. Test matrix (as executed)

Environment: WSL2 Ubuntu 22.04, Python 3.11.17, pysam 0.24.1, samtools 1.24,
bcftools 1.24, bwa-mem2 2.3, fastp 1.4.0, mosdepth 0.3.14, GATK 4.6.2.0,
snakemake 9.27.0.

```bash
python -m pytest tests/ -q          # 43 tests, all passing
bash scripts/wsl_pipeline_smoke.sh  # real FASTQ -> VCF -> AVRS review
```

| Test file | Level | Covers |
|---|---|---|
| `tests/test_vcf.py` | ✅ | VCF header/record parsing, multiallelic splitting + AD subsetting, missing-ALT handling, indel left-alignment, SNV normalisation, build-scoped stable variant IDs, variant type classification |
| `tests/test_evidence.py` | ✅ | BAM write/read round-trip, CIGAR parsing, BAI index vs sequential-scan equivalence (200 reads, 4 windows), SNV ref/alt counting, **indel support via CIGAR `D` ops (not base pileup)**, empty coverage, strand-bias detection, rule engine hits, review-mode gating |
| `tests/test_igv.py` | ✅ | app-level client → official IGV tool mapping, unsupported operation reporting, snapshot capture + metadata binding, snapshot failure handling, instance-lock exclusivity, render-verification rejection of missing files |
| `tests/test_ai.py` | ✅ | strict schema (extra fields rejected, 1-based position), recommendation labels, insufficient-evidence path, "no visual review claimed without a vision model", prompt hard-rules, JSON round-trip, task state transitions, provider selection |
| `tests/test_e2e_review.py` | ✅ | full Mode-B loop on the demo dataset: import → evidence → triage → mock-IGV review → AI review → DB → human decision history → JSON/CSV/HTML/QC/exception/audit reports |
| `tests/test_fn_discovery.py` | ✅ | hotspot parsing, targeted-locus classification (`in_vcf` / `possible_missed_variant` / `insufficient_coverage` / `no_evidence_of_alt`), filtered-recovery candidates, disclaimer text |

---

## 3. Component validation status

| Component | Level | Notes |
|---|---|---|
| Pure-Python BGZF/BAM/BAI | 🟡 | Round-trip and index-vs-scan equivalence are tested; **cross-validation against `samtools index`/`samtools view` has not been performed in this environment** (the production path uses pysam, which is samtools-backed). |
| pysam evidence backend | ✅ | Exercised by the demo and e2e tests in WSL (pysam 0.24.1). |
| Pure-python evidence backend | ✅ | Exercised by unit tests; it is the fallback for platforms where pysam cannot be built. |
| VCF import + normalisation | 🟡 | Logic tested; the production path additionally runs `bcftools norm` in the workflow. Agreement between the two normalisers has **not** been measured on a large panel. |
| Triage rules | 🟡 | Unit-tested for firing behaviour; thresholds are defaults, **not** validated against a truth set. |
| IGV integration (wire protocol) | ✅ | Command mapping and failure handling tested against a mock listener that implements the real protocol. |
| IGV integration (**real IGV Desktop**) | 🟡 | The client was written against the official `igv-mcp` source; the *real* IGV Desktop was **not running in this environment**, so end-to-end capture with a genuine IGV rendering is **not yet verified**. See §6. |
| Vision review | 🔵 | The provider interfaces (Ollama / OpenAI-compatible) are implemented and unit-tested with a mock provider; **no real vision model was invoked in this environment**, therefore no visual finding has been validated. |
| AI JSON contract | ✅ | Strict schema with `extra="forbid"`; provider output is validated and a failure is recorded as `REVIEW_FAILED`/`INCONCLUSIVE`, never silently accepted. |
| Human review workbench | ✅ | API + UI exercised through the e2e test (decision history, supersede chain, audit rows). Browser rendering verified manually during development. |
| FN discovery | 🟡 | Candidate generation and classification are tested; **no confirmed false-negative rate can be claimed** without external truth. |
| Report generation | ✅ | All seven artefacts produced and content-checked in the e2e test. |
| FASTQ → VCF workflow | ✅ | The Snakemake workflow was executed end to end on a synthetic 120 kb reference with 18 000 read pairs (fastp → bwa-mem2 → samtools → GATK HaplotypeCaller → bcftools norm/filter → AVRS review). |
| Somatic / trio / CNV / SV / mosaic / long-read | ⬜ | Not implemented. The configuration separates germline from somatic, but only germline logic exists. |
| DeepVariant / second-caller comparison | 🔵 | `SecondCallerAdapter` interface + normalised set comparison exist; no concrete second caller is wired up. |
| GIAB / hap.py benchmark harness | 🔵 | Designed (see §5) but not executed here. |
| Authentication / RBAC / encryption at rest | ⬜ | Interface-level only (see [SECURITY.md](SECURITY.md)). |

---

## 4. Demo dataset (not a truth set)

`avrs demo` generates a synthetic 2-contig genome (30 kb + 15 kb) with eight
simulated scenarios and no patient data. The scenario labels exist to prove
that the pipeline *reaches* the expected code paths:

| Scenario | Signal | Path exercised |
|---|---|---|
| clean SNV | balanced strands, high BQ/MQ | no rule hit → not queued |
| low-quality SNV | all reads below the MAPQ floor | `ALL_READS_LOW_MQ`, `LOW_ALT_AD` → `INSUFFICIENT_EVIDENCE` |
| strand-biased SNV | all ALT reads on one strand | `STRAND_BIAS` → `SUSPECTED_FALSE_POSITIVE` |
| homopolymer indel | insertion in an 11 bp poly-A | `INDEL_IN_REPEAT` → `INCONCLUSIVE` (representation risk) |
| low-coverage locus | DP below floor | `LOW_DP` → `INSUFFICIENT_EVIDENCE` |
| caller miss | ALT support, absent from VCF | FN candidate (`possible_missed_variant`) |
| equivocal | caller AD ≠ pileup AD | `AD_PILEUP_MISMATCH` → `INCONCLUSIVE` |
| second contig | clean SNV on chr2 | contig handling |

**These labels are not clinical truth and must never be used as a truth set.**

---

## 5. Truth-set benchmark (design; not yet executed)

Intended harness for GIAB HG002 (or an equivalent truth set):

```
truth VCF + high-confidence BED  ─┐
sample FASTQ/BAM (same reference) ─┼─► hap.py (haplotype-aware) ─► TP/FP/FN
                                  ─┘
```

Requirements for the implementation (not yet in the codebase):

* metrics only inside the **defined confident regions**;
* SNV and indel reported separately;
* stratification by genomic region (e.g. low-complexity / homopolymer);
* evaluation **before and after** the AI review stage, so the assistant's
  effect on the workload is measurable;
* candidate detection rate, human workload and misjudgement risk evaluated
  separately;
* AI prediction labels are **never** used as truth.

GIAB resources used in such a benchmark:

| Resource | Use |
|---|---|
| `HG002_GRCh38_1_22_v4.2.1_benchmark_vcf.gz` | truth calls |
| `HG002_GRCh38_1_22_v4.2.1_benchmark_vcf.bcf` | truth (BCF) |
| `HG002_GRCh38_1_22_v4.2.1_highconf.bed` | confident regions |
| `HG002_GRCh38_1_22_v4.2.1.1_alignment.bam` | pre-aligned BAM |
| `HG002_GRCh38_1_22_v4.2.1_fastq/` | raw FASTQ |

---

## 6. Acceptance criteria for production use

| # | Criterion | Status |
|---|---|---|
| 1 | Unit + integration suite green in the target environment | ✅ met (43 tests) |
| 2 | Real FASTQ → VCF → review smoke test reproducible | ✅ met (`scripts/wsl_pipeline_smoke.sh`) |
| 3 | Snapshot ↔ variant binding verified against **real IGV Desktop** | ⬜ open — requires an IGV instance with the port listener enabled |
| 4 | Vision review executed by a real vision model and its findings audited | ⬜ open — needs an Ollama/OpenAI-compatible endpoint with vision support |
| 5 | Rule thresholds calibrated on a local truth set | ⬜ open — lab-specific work |
| 6 | Truth-set benchmark (hap.py) executed; TP/FP/FN/precision/recall reported per variant class | ⬜ open |
| 7 | Inter-reviewer agreement study (human vs AI) | ⬜ open |
| 8 | Workload-reduction baseline measured (same locus set, timed) | ⬜ open |
| 9 | Security review of the deployment environment (auth, encryption, key management) | ⬜ open |
| 10 | Documented risk analysis + sign-off for the intended use | ⬜ open |

**Acceptance principle:** reducing human review time must not come at the
cost of an unacceptable false-negative risk for real variants. Until criteria
3–10 are met, AVRS must be used as an *assistive* tool with human sign-off on
every reported variant — never as an autonomous reporter.

---

## 7. Reproducing the checks

```bash
# unit + integration
python -m pytest tests/ -q

# demo end-to-end (offline, mock IGV)
avrs demo --config config/demo.yaml

# real pipeline smoke test (needs the WSL toolchain)
bash scripts/wsl_pipeline_smoke.sh

# environment check
avrs doctor --config config/demo.yaml
```
