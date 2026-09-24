#!/usr/bin/env python3
"""Make a venv's console scripts safe to run from a path containing spaces.

pip writes each console script with an unquoted absolute interpreter path in
its shebang, which the kernel cannot parse once the path has a space in it
(~/Library/Application Support, Google Drive, ...). Rewriting them into the
/bin/sh + exec form fixes every script, including helpers that tools like
RagTag invoke by name.
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

SHIM = "#!/bin/sh\n'''exec' \"{py}\" \"$0\" \"$@\"\n' '''\n"


def fix(venv: Path) -> list[str]:
    py = venv / "bin" / "python"
    fixed: list[str] = []
    for f in sorted((venv / "bin").iterdir()):
        if not f.is_file() or f.is_symlink():
            continue
        try:
            head = f.open("rb").readline()
        except OSError:
            continue
        if not head.startswith(b"#!") or b" " not in head[2:].strip():
            continue
        body = f.read_bytes().split(b"\n", 1)[1]
        f.write_bytes(SHIM.format(py=py).encode() + body)
        os.chmod(f, 0o755)
        fixed.append(f.name)
    return fixed


if __name__ == "__main__":
    target = Path(sys.argv[1] if len(sys.argv) > 1
                  else Path(__file__).resolve().parents[1] / ".venv")
    names = fix(target)
    print(f"rewrote {len(names)} shebang(s): {', '.join(names) if names else 'none needed'}")
