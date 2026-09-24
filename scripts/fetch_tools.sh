#!/usr/bin/env bash
# Build the bundled arm64 toolchain into vendor/bin.
#
# Everything here is compiled from upstream release tarballs against the macOS
# SDK, so the resulting binaries are native Apple Silicon and carry no conda or
# Homebrew dependency. Python-based tools go into the app's own .venv.
set -uo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BIN="$ROOT/vendor/bin"
BUILD="${NANOASSEMBLER_BUILD_DIR:-$ROOT/vendor/build}"
PYROOT="$ROOT/vendor/python"
PY="$PYROOT/bin/python3"
LOG="$BUILD/build.log"

MINIMAP2_VER=2.28
SAMTOOLS_VER=1.21
BCFTOOLS_VER=1.21
FLYE_VER=2.9.5
RACON_VER=1.5.0
# A relocatable CPython: a normal venv records an absolute path to the base
# interpreter, so it cannot be moved into an .app bundle or onto another Mac.
PYSTANDALONE_TAG=20260901
PYSTANDALONE_VER=3.12.14

mkdir -p "$BIN" "$BUILD"
: > "$LOG"
JOBS=$(sysctl -n hw.ncpu 2>/dev/null || echo 4)

say()  { printf '\n=== %s ===\n' "$*" | tee -a "$LOG"; }
ok()   { printf '  [ok]   %s\n' "$*"; }
warn() { printf '  [warn] %s\n' "$*"; }

if ! xcode-select -p >/dev/null 2>&1; then
  echo "Xcode command line tools are required: xcode-select --install" >&2
  exit 1
fi

fetch() { # url outfile
  [ -f "$2" ] || curl -fsSL "$1" -o "$2"
}

# --- minimap2 -------------------------------------------------------------
build_minimap2() {
  say "minimap2 $MINIMAP2_VER"
  local d="$BUILD/minimap2-$MINIMAP2_VER"
  fetch "https://github.com/lh3/minimap2/archive/refs/tags/v$MINIMAP2_VER.tar.gz" \
        "$BUILD/minimap2.tar.gz" || return 1
  [ -d "$d" ] || tar -xzf "$BUILD/minimap2.tar.gz" -C "$BUILD"
  ( cd "$d" && make -j"$JOBS" arm_neon=1 aarch64=1 ) >>"$LOG" 2>&1 || return 1
  cp "$d/minimap2" "$BIN/" && ok "minimap2"
}

# --- htslib-based tools ---------------------------------------------------
build_htslib_tool() { # name version
  local name=$1 ver=$2
  say "$name $ver"
  local d="$BUILD/$name-$ver"
  fetch "https://github.com/samtools/$name/releases/download/$ver/$name-$ver.tar.bz2" \
        "$BUILD/$name.tar.bz2" || return 1
  [ -d "$d" ] || tar -xjf "$BUILD/$name.tar.bz2" -C "$BUILD"
  ( cd "$d" \
    && ./configure --disable-libcurl --without-curses --disable-plugins \
         --prefix="$BUILD/inst" >>"$LOG" 2>&1 \
    && make -j"$JOBS" >>"$LOG" 2>&1 ) || return 1
  cp "$d/$name" "$BIN/" && ok "$name"
}

# --- racon (optional polisher) -------------------------------------------
build_racon() {
  say "racon $RACON_VER (optional)"
  command -v cmake >/dev/null 2>&1 || { warn "cmake not found; skipping racon"; return 1; }
  local d="$BUILD/racon-v$RACON_VER"
  fetch "https://github.com/lbcb-sci/racon/releases/download/$RACON_VER/racon-v$RACON_VER.tar.gz" \
        "$BUILD/racon.tar.gz" || return 1
  [ -d "$d" ] || tar -xzf "$BUILD/racon.tar.gz" -C "$BUILD"
  ( cd "$d" && mkdir -p build && cd build \
    && cmake -DCMAKE_BUILD_TYPE=Release -Dracon_build_wrapper=OFF .. >>"$LOG" 2>&1 \
    && make -j"$JOBS" >>"$LOG" 2>&1 ) || return 1
  cp "$d/build/bin/racon" "$BIN/" && ok "racon"
}

# --- relocatable python runtime -------------------------------------------
build_portable_python() {
  say "portable CPython $PYSTANDALONE_VER"
  if [ -x "$PY" ]; then ok "already present"; return 0; fi
  local url="https://github.com/astral-sh/python-build-standalone/releases/download"
  url="$url/$PYSTANDALONE_TAG/cpython-${PYSTANDALONE_VER}%2B${PYSTANDALONE_TAG}-aarch64-apple-darwin-install_only.tar.gz"
  fetch "$url" "$BUILD/python.tar.gz" || return 1
  rm -rf "$BUILD/python"
  tar -xzf "$BUILD/python.tar.gz" -C "$BUILD" || return 1
  mkdir -p "$(dirname "$PYROOT")"
  rm -rf "$PYROOT"
  mv "$BUILD/python" "$PYROOT" || return 1
  "$PY" -m pip -q install --upgrade pip setuptools wheel >>"$LOG" 2>&1
  ok "python $("$PY" -c 'import platform;print(platform.python_version())')"
}

# --- python-based tools ---------------------------------------------------
build_python_tools() {
  say "python tools (flye, ragtag, medaka, GUI)"
  [ -x "$PY" ] || { warn "no portable python; skipping"; return 1; }

  "$PY" -m pip -q install PySide6 >>"$LOG" 2>&1 && ok "PySide6" \
    || warn "PySide6 install failed - the GUI will not start"

  if "$PY" -m pip install \
       "https://github.com/mikolmogorov/Flye/archive/refs/tags/$FLYE_VER.tar.gz" \
       >>"$LOG" 2>&1; then ok "flye"; else warn "flye build failed - see $LOG"; fi

  "$PY" -m pip -q install RagTag >>"$LOG" 2>&1 && ok "ragtag" || warn "ragtag failed"

  # Medaka has no guaranteed macOS/arm64 wheel. If it does not install the app
  # still runs: samtools consensus and bcftools take over, and the detected
  # basecalling model is recorded but not used for neural polishing.
  if "$PY" -m pip install medaka pyabpoa >>"$LOG" 2>&1; then
    ok "medaka"
  else
    warn "medaka unavailable - falling back to samtools/bcftools"
  fi

  # pip writes an unquoted absolute interpreter path into every console
  # script's shebang, which the kernel cannot parse once the install directory
  # contains a space. RagTag also shells out to its own helper scripts by name,
  # so all of them have to be fixed, not only the entry points wrapped below.
  "$PY" "$ROOT/scripts/fix_shebangs.py" "$PYROOT" 2>&1 | tee -a "$LOG" | sed 's/^/  /'

  # Wrappers resolve the interpreter relative to their own location, so
  # vendor/ can be copied into an .app bundle or another Mac unchanged.
  for name in flye ragtag.py medaka; do
    [ -f "$PYROOT/bin/$name" ] || continue
    cat > "$BIN/$name" <<'WRAPPER'
#!/bin/sh
DIR=$(cd "$(dirname "$0")" && pwd)
exec "$DIR/../python/bin/python3" "$DIR/../python/bin/__NAME__" "$@"
WRAPPER
    sed -i '' "s|__NAME__|$name|" "$BIN/$name"
    chmod +x "$BIN/$name"
  done
}

FAILED=()
build_minimap2                        || FAILED+=(minimap2)
build_htslib_tool samtools "$SAMTOOLS_VER" || FAILED+=(samtools)
build_htslib_tool bcftools "$BCFTOOLS_VER" || FAILED+=(bcftools)
build_racon                           || true
build_portable_python                 || FAILED+=(python)
build_python_tools

say "installed in $BIN"
ls -1 "$BIN" 2>/dev/null | sed 's/^/  /'
if [ ${#FAILED[@]} -gt 0 ]; then
  echo "FAILED (required): ${FAILED[*]}  -- see $LOG" >&2
  exit 1
fi
echo "toolchain ready"
