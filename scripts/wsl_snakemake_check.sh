#!/bin/bash
export MAMBA_ROOT_PREFIX="$HOME/micromamba"
export PATH="$HOME/micromamba/envs/avrs/bin:$PATH"
cd /mnt/f/codex/mcp-igv
snakemake -s workflow/Snakefile --configfile "${1:-config/pipeline.example.yaml}" --list "$@" 2>&1 | tail -30
