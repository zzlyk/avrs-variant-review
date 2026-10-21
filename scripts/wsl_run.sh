#!/bin/bash
# Run any command inside the AVRS micromamba env in WSL2.
# Usage (from Windows): wsl bash /mnt/f/codex/mcp-igv/scripts/wsl_run.sh <args...>
export MAMBA_ROOT_PREFIX="$HOME/micromamba"
export PATH="$HOME/micromamba/envs/avrs/bin:$PATH"
export PYTHONPATH="/mnt/f/codex/mcp-igv/src${PYTHONPATH:+:$PYTHONPATH}"
export PYTHONUNBUFFERED=1
cd /mnt/f/codex/mcp-igv || exit 1
exec "$@"
