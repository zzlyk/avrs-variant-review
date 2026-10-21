#!/bin/bash
# Real FASTQ -> VCF smoke test using the actual toolchain (WSL2/Linux).
# Usage: wsl bash /mnt/f/codex/mcp-igv/scripts/wsl_pipeline_smoke.sh
set -e
export MAMBA_ROOT_PREFIX="$HOME/micromamba"
export PATH="$HOME/micromamba/envs/avrs/bin:$PATH"
cd /mnt/f/codex/mcp-igv
WORK="$(pwd)/smoke_run"
rm -rf "$WORK"
mkdir -p "$WORK"

echo "== 1/4 generating FASTQ + reference =="
python scripts/make_smoke_fastq.py --out "$WORK/data"

# WSL2 note: bwa-mem2 cannot write its index on DrvFs (/mnt/...), so the
# reference is staged on a native filesystem and indexed there.
NATIVE_REF_DIR="$HOME/avrs_smoke_ref"
rm -rf "$NATIVE_REF_DIR" && mkdir -p "$NATIVE_REF_DIR"
cp "$WORK/data/smoke_reference.fa" "$NATIVE_REF_DIR/reference.fa"
cp "$WORK/data/smoke_reference.fa.fai" "$NATIVE_REF_DIR/reference.fa.fai"
cp "$WORK/data/smoke_reference.fa.dict" "$NATIVE_REF_DIR/reference.fa.dict"
cp "$WORK/data/smoke_reference.dict" "$NATIVE_REF_DIR/reference.dict"
bwa-mem2 index "$NATIVE_REF_DIR/reference.fa" > "$WORK/bwa_index.log" 2>&1
echo "reference staged + indexed on native fs: $NATIVE_REF_DIR"

cat > "$WORK/pipeline.yaml" <<YAML
sample_id: SMOKE001
sample_type: germline_dna
threads: 4
output_dir: $WORK/out
fastq:
  r1: $WORK/data/SMOKE_R1.fastq.gz
  r2: $WORK/data/SMOKE_R2.fastq.gz
reference:
  fasta: $NATIVE_REF_DIR/reference.fa
  build: SMOKE1
  known_sites: []
target_bed: ""
caller:
  extra: ""
filtering:
  mode: hard
  snp_min_gq: 20
  snp_min_dp: 10
  snp_min_ad: 3
  snp_min_qs: 20
  ind_min_gq: 20
  ind_min_dp: 10
  ind_min_ad: 3
  ind_min_qs: 20
YAML

echo "== 2/4 running snakemake (fastp -> bwa-mem2 -> samtools -> GATK -> bcftools) =="
snakemake -s workflow/Snakefile --configfile "$WORK/pipeline.yaml" \
  --directory "$WORK" --cores 4 --rerun-incomplete --printshellcmds \
  > "$WORK/snakemake.log" 2>&1 || {
    echo "SNAKEMAKE FAILED - tail of log:"; tail -40 "$WORK/snakemake.log"; exit 1; }

echo "== 3/4 verifying outputs =="
VCF="$WORK/out/variants/final.filtered.norm.vcf.gz"
ls -la "$VCF" "$VCF.tbi" "$WORK/out/alignment/flagstat.txt" "$WORK/out/qc/fastp.json" || true
echo "--- variant calls ---"
bcftools view -H "$VCF" | head -10
echo "--- call count ---"
bcftools view -H "$VCF" | wc -l
echo "--- flagstat ---"
cat "$WORK/out/alignment/flagstat.txt"
echo "--- fastp summary ---"
python - <<PY
import json
d = json.load(open("$WORK/out/qc/fastp.json"))
print({k: d["summary"][k] for k in ("before_filtering", "after_filtering")})
PY

echo "== 4/4 feeding the result into AVRS review-only mode =="
python - <<PY
import sys
sys.path.insert(0, "src")
from avrs.config import AVRSConfig
from avrs.db import open_db
from avrs.review import import_sample_bam_vcf, triage_sample, run_reviews
from avrs.variants import import_vcf
from avrs.io.fasta import Fasta

cfg = AVRSConfig()
cfg.paths.output_dir = "$WORK/avrs"
cfg.paths.snapshot_dir = "$WORK/avrs/snapshots"
cfg.paths.db = "$WORK/avrs/avrs.db"
cfg.reference.build = "SMOKE1"
cfg.reference.fasta = "$NATIVE_REF_DIR/reference.fa"
cfg.analysis.data_type = "panel"
cfg.review.snapshot_wait_initial_s = 0.3
cfg.review.snapshot_poll_interval_s = 0.2
cfg.review.snapshot_timeout_s = 15
cfg.igv.wait.load_wait_s = 0.3
cfg.igv.wait.goto_wait_s = 0.2

engine, session = open_db(cfg.paths.db)
smp = import_sample_bam_vcf(session, cfg, "SMOKE001",
                            "$WORK/out/alignment/dedup.bam",
                            "$WORK/out/variants/final.filtered.norm.vcf.gz")
fasta = Fasta(cfg.reference.fasta)
n, skipped = import_vcf(session, "$WORK/out/variants/final.filtered.norm.vcf.gz",
                        smp.id, "SMOKE001", "SMOKE1", fasta=fasta)
session.commit()
print("imported variants:", n, "skipped:", skipped)
tasks = triage_sample(session, cfg, "SMOKE001", cfg.reference.fasta,
                       "$WORK/out/alignment/dedup.bam")
print("review tasks:", tasks)
stats = run_reviews(session, cfg, "SMOKE001", "$WORK/out/alignment/dedup.bam",
                    "$WORK/out/variants/final.filtered.norm.vcf.gz",
                    cfg.reference.fasta, use_mock_igv=True)
print("review stats:", stats)
PY

echo "SMOKE TEST OK"
