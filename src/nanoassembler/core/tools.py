"""Locate the external binaries the pipeline shells out to.

Search order is bundle-first so a packaged .app never picks up a stray copy
from the user's PATH:
    1. <App>.app/Contents/Resources/bin      (the shipped bundle)
    2. <repo>/vendor/bin                     (developer checkout)
    3. $NANOASSEMBLER_BIN
    4. PATH
"""
from __future__ import annotations

import os
import re
import shutil
import subprocess
import sys
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Optional


@dataclass(frozen=True)
class ToolSpec:
    name: str
    required: bool
    purpose: str
    version_args: tuple[str, ...] = ("--version",)
    version_re: str = r"([0-9]+\.[0-9]+(?:\.[0-9]+)?)"


TOOLS: tuple[ToolSpec, ...] = (
    ToolSpec("minimap2", True, "read mapping"),
    ToolSpec("samtools", True, "BAM handling, depth, consensus"),
    ToolSpec("bcftools", True, "variant calling fallback"),
    ToolSpec("flye", True, "de novo assembly", ("--version",)),
    ToolSpec("racon", False, "assembly polishing when Medaka is unavailable"),
    ToolSpec("medaka", False, "model-aware consensus and variant calling"),
    ToolSpec("ragtag.py", False, "ordering contigs against a reference"),
)


def bundle_bin_dirs() -> list[Path]:
    dirs: list[Path] = []
    # Inside a PyInstaller/py2app bundle sys.executable is .../Contents/MacOS/x
    exe_dir = Path(sys.executable).resolve().parent
    for parent in (exe_dir, *exe_dir.parents):
        if parent.name == "Contents":
            dirs.append(parent / "Resources" / "bin")
            break
    meipass = getattr(sys, "_MEIPASS", None)
    if meipass:
        dirs.append(Path(meipass) / "bin")
    # Developer checkout: src/nanoassembler/core/tools.py -> repo root
    repo_root = Path(__file__).resolve().parents[3]
    dirs.append(repo_root / "vendor" / "bin")
    env = os.environ.get("NANOASSEMBLER_BIN")
    if env:
        dirs.extend(Path(p) for p in env.split(os.pathsep) if p)
    # Flye ships its own flye-minimap2 / flye-samtools next to the interpreter
    # that installed it, and calls them by name, so that directory has to be on
    # PATH as well.
    dirs.append(Path(sys.executable).parent)   # venv bin, symlink not resolved
    dirs.append(Path(sys.prefix) / "bin")
    seen: set[Path] = set()
    unique = []
    for d in dirs:
        if d.is_dir() and d not in seen:
            seen.add(d)
            unique.append(d)
    return unique


def find_tool(name: str) -> Optional[Path]:
    for d in bundle_bin_dirs():
        cand = d / name
        if cand.is_file() and os.access(cand, os.X_OK):
            return cand
    found = shutil.which(name)
    return Path(found) if found else None


def tool_env(threads: int = 1) -> dict[str, str]:
    """Environment for subprocesses: bundled bin dirs first on PATH.

    Medaka, Flye and RagTag call minimap2/samtools themselves, so they must see
    the same binaries we do. ``threads`` caps the BLAS/OpenMP thread pools that
    Medaka's torch backend would otherwise size to the whole machine, which
    matters when several samples run at once.
    """
    env = dict(os.environ)
    dirs = [str(d) for d in bundle_bin_dirs()]
    if dirs:
        env["PATH"] = os.pathsep.join(dirs + [env.get("PATH", "")])
    for var in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS",
                "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
        env[var] = str(max(1, threads))
    return env


@dataclass
class ToolStatus:
    spec: ToolSpec
    path: Optional[Path]
    version: str = ""
    error: str = ""

    @property
    def ok(self) -> bool:
        return self.path is not None and not self.error

    @property
    def blocking(self) -> bool:
        return self.spec.required and not self.ok


def probe(spec: ToolSpec) -> ToolStatus:
    path = find_tool(spec.name)
    if path is None:
        return ToolStatus(spec, None, error="not found")
    try:
        proc = subprocess.run(
            [str(path), *spec.version_args],
            capture_output=True,
            text=True,
            timeout=60,
            env=tool_env(),
        )
    except (OSError, subprocess.SubprocessError) as exc:
        return ToolStatus(spec, path, error=f"failed to run: {exc}")
    text = (proc.stdout + proc.stderr).strip()
    m = re.search(spec.version_re, text)
    return ToolStatus(spec, path, version=m.group(1) if m else text.splitlines()[0][:40])


@lru_cache(maxsize=1)
def probe_all_cached() -> tuple[ToolStatus, ...]:
    return tuple(probe(spec) for spec in TOOLS)


def probe_all(refresh: bool = False) -> list[ToolStatus]:
    if refresh:
        probe_all_cached.cache_clear()
    return list(probe_all_cached())


def status_map(refresh: bool = False) -> dict[str, ToolStatus]:
    return {s.spec.name: s for s in probe_all(refresh)}


def missing_required(refresh: bool = False) -> list[str]:
    return [s.spec.name for s in probe_all(refresh) if s.blocking]


def have(name: str) -> bool:
    st = status_map().get(name)
    return bool(st and st.ok)


def require(name: str) -> Path:
    st = status_map().get(name)
    if not st or not st.ok or st.path is None:
        raise FileNotFoundError(
            f"{name} is not available. Run scripts/fetch_tools.sh to build the "
            f"bundled toolchain, or install it and restart the app."
        )
    return st.path
