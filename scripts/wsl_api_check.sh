#!/bin/bash
export MAMBA_ROOT_PREFIX="$HOME/micromamba"
export PATH="$HOME/micromamba/envs/avrs/bin:$PATH"
export PYTHONPATH="/mnt/f/codex/mcp-igv/src"
cd /mnt/f/codex/mcp-igv
python -m uvicorn avrs.api.app:create_app --factory --host 127.0.0.1 --port 8077 --log-level warning &
SRV=$!
sleep 6
echo "=== /api/samples ==="
curl -s http://127.0.0.1:8077/api/samples | head -c 600; echo
echo "=== /api/variants ==="
curl -s "http://127.0.0.1:8077/api/variants?limit=3" | head -c 700; echo
echo "=== /api/metrics ==="
curl -s http://127.0.0.1:8077/api/metrics | head -c 500; echo
echo "=== / (UI) ==="
curl -s -o /tmp/ui.html -w "http_status=%{http_code} bytes=%{size_download}\n" http://127.0.0.1:8077/
grep -c "AVRS" /tmp/ui.html
echo "=== /api/fn-candidates ==="
curl -s "http://127.0.0.1:8077/api/fn-candidates?sample_id=DEMO_SAMPLE_001" | head -c 300; echo
kill $SRV 2>/dev/null
echo "API CHECK DONE"
