"""Regression test: sample names and output folders containing spaces.

Flye refuses any path with a space in it, and FASTQ files named after samples
("SGGFLOX A.fastq") are common, so both the sample-folder naming and the Flye
staging fallback are exercised here with de novo assembly.
"""
from __future__ import annotations

import shutil
import sys
import tempfile
import threading
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from nanoassembler.core.models import Mode, RunConfig, Sample  # noqa: E402
from nanoassembler.core.pipeline import run_sample, safe_dirname, unique_dirnames  # noqa: E402

data = Path(sys.argv[1])

assert safe_dirname("112014079811_SGGFLOX A") == "112014079811_SGGFLOX_A"
assert safe_dirname("a/b:c") == "a_b_c"
names = unique_dirnames(["A B", "A_B", "a b"])
assert len({n.lower() for n in names.values()}) == 3, names
print("folder naming OK")

root = Path(tempfile.mkdtemp(prefix="na spaces "))  # deliberately contains a space
try:
    # A spaced FASTQ path as well, like the files that triggered the bug
    fq = root / "SGGFLOX A.fastq.gz"
    shutil.copy(data / "sample01.fastq.gz", fq)
    s = Sample(name="SGGFLOX A", fastqs=[fq], mode=Mode.DENOVO)
    cfg = RunConfig(out_dir=root / "my runs", threads_per_job=4,
                    min_read_len=500, min_read_q=7.0)
    r = run_sample(s, {}, cfg, root / "my runs" / "run", lambda m: None, threading.Event())
    assert r.ok, r.message
    assert r.out_dir.name == "SGGFLOX_A", r.out_dir
    assert (r.out_dir / "assembly.fasta").stat().st_size > 0
    assert (r.out_dir / "work" / "flye" / "flye.log").exists(), "flye.log not copied back"
    leftovers = list(Path(tempfile.gettempdir()).glob("nanoassembler_flye_*"))
    assert not leftovers, f"staging directories left behind: {leftovers}"
    print(f"de novo through spaced paths OK: {r.metrics['assembly_bp']:,} bp")
finally:
    shutil.rmtree(root, ignore_errors=True)
print("SPACES TEST PASSED")
