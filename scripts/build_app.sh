#!/usr/bin/env bash
# Package NanoAssembler.app for Apple Silicon.
#
# The bundle is assembled by hand rather than with PyInstaller: vendor/python is
# already a relocatable interpreter carrying PySide6, Flye, Medaka and RagTag,
# so freezing would only duplicate PySide6 and torch inside the same .app.
# Everything here is a copy, which keeps the result inspectable and reproducible.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DIST="${1:-$ROOT/dist}"
APP="$DIST/NanoAssembler.app"
RES="$APP/Contents/Resources"
VERSION=$(sed -n 's/^__version__ = "\(.*\)"/\1/p' "$ROOT/src/nanoassembler/__init__.py")

for required in "$ROOT/vendor/python/bin/python3" "$ROOT/vendor/bin/minimap2" \
                "$ROOT/vendor/bin/samtools"; do
  [ -e "$required" ] || { echo "missing $required - run scripts/fetch_tools.sh first" >&2; exit 1; }
done

rm -rf "$APP"
mkdir -p "$APP/Contents/MacOS" "$RES"

echo "copying interpreter and tools…"
cp -R "$ROOT/vendor/python" "$RES/python"
cp -R "$ROOT/vendor/bin" "$RES/bin"
mkdir -p "$RES/src"
cp -R "$ROOT/src/nanoassembler" "$RES/src/nanoassembler"
find "$RES/src" -name '__pycache__' -type d -exec rm -rf {} + 2>/dev/null || true

cat > "$APP/Contents/MacOS/NanoAssembler" <<'LAUNCH'
#!/bin/sh
# Resolve everything relative to the bundle so it runs from any location.
RES=$(cd "$(dirname "$0")/../Resources" && pwd)
export PYTHONPATH="$RES/src${PYTHONPATH:+:$PYTHONPATH}"
export PYTHONDONTWRITEBYTECODE=1
exec "$RES/python/bin/python3" -m nanoassembler "$@"
LAUNCH
chmod +x "$APP/Contents/MacOS/NanoAssembler"

cat > "$APP/Contents/Info.plist" <<PLIST
<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE plist PUBLIC "-//Apple//DTD PLIST 1.0//EN"
  "http://www.apple.com/DTDs/PropertyList-1.0.dtd">
<plist version="1.0">
<dict>
  <key>CFBundleName</key><string>NanoAssembler</string>
  <key>CFBundleDisplayName</key><string>NanoAssembler</string>
  <key>CFBundleExecutable</key><string>NanoAssembler</string>
  <key>CFBundleIdentifier</key><string>org.local.nanoassembler</string>
  <key>CFBundleVersion</key><string>$VERSION</string>
  <key>CFBundleShortVersionString</key><string>$VERSION</string>
  <key>CFBundlePackageType</key><string>APPL</string>
  <key>LSMinimumSystemVersion</key><string>12.0</string>
  <key>NSHighResolutionCapable</key><true/>
  <key>LSApplicationCategoryType</key><string>public.app-category.education</string>
</dict>
</plist>
PLIST

# Ad-hoc signature: without one, Gatekeeper kills a bundle carrying unsigned
# dylibs on first launch.
echo "signing…"
codesign --force --deep --sign - "$APP" >/dev/null 2>&1 \
  && echo "  ad-hoc signed" \
  || echo "  codesign failed - the app still opens via right-click > Open"

echo
echo "built: $APP  ($(du -sh "$APP" | cut -f1))"
echo "test it with:  open '$APP'"
