#!/usr/bin/env bash
# Launch NanoAssembler from the source tree without building the .app.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PY="$ROOT/vendor/python/bin/python3"
[ -x "$PY" ] || { echo "run scripts/fetch_tools.sh first" >&2; exit 1; }
exec env PYTHONPATH="$ROOT/src" "$PY" -m nanoassembler "$@"
