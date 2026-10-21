#!/bin/bash
# AVRS WSL2 environment setup - no root required (micromamba + bioconda)
set -e
export MAMBA_ROOT_PREFIX="$HOME/micromamba"
MAMBA="$HOME/bin/micromamba"

if [ ! -x "$MAMBA" ]; then
  mkdir -p "$HOME/bin"
  curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xj -C "$HOME/bin" --strip-components=1 bin/micromamba
fi
"$MAMBA" --version

echo "[avrs-setup] creating conda env: avrs (python3.11 + bio tools + pysam) ..."
"$MAMBA" create -y -n avrs -c conda-forge -c bioconda \
  python=3.11 pip pysam samtools bcftools bwa-mem2 fastp mosdepth gatk4

echo "[avrs-setup] installing python project deps ..."
"$MAMBA" run -n avrs pip install --quiet pydantic pyyaml sqlalchemy typer jinja2 httpx pytest fastapi uvicorn snakemake

echo "[avrs-setup] verifying ..."
"$MAMBA" run -n avrs bash -c '
python --version && python -c "import pysam; print(\"pysam\", pysam.__version__)" &&
samtools --version | head -1 && bcftools --version | head -1 &&
bwa-mem2 version 2>&1 | head -1 && fastp --version 2>&1 | head -1 &&
mosdepth --version && (gatk --version 2>&1 | grep -i gatk | head -1 || true) &&
snakemake --version'
echo "[avrs-setup] DONE"
