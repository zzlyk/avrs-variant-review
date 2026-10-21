#!/bin/bash
# WSL environment check for AVRS
for t in samtools bcftools bwa-mem2 bwa snakemake fastp java node conda mosdepth pip3 pip tabix wget curl unzip python3.10 python3.11 python3.12 virtualenv; do
  p=$(command -v "$t" 2>/dev/null)
  if [ -n "$p" ]; then echo "HAVE $t -> $p"; else echo "MISS $t"; fi
done
echo "python3: $(python3 --version 2>&1)"
java -version 2>&1 | head -1 || true
apt list --installed 2>/dev/null | grep -E "^(samtools|bcftools|bwa)" | head -5
ls -la /mnt/f/codex/mcp-igv/ | head -8
