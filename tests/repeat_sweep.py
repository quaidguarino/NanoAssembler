"""Map out when de novo assembly of a repeat-containing construct breaks.

Builds constructs with N identical cassettes of a given size, simulates ONT
reads of a given mean length, assembles, and reports the assembled size and the
number of repeat copies actually recovered.
"""
from __future__ import annotations

import random
import subprocess
import sys
import threading
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "tests"))

from make_repeat_test import simulate  # noqa: E402
from make_test_data import random_seq, write_fasta  # noqa: E402
from nanoassembler.core.models import Mode, RunConfig, Sample  # noqa: E402
from nanoassembler.core.pipeline import run_sample  # noqa: E402

MM2 = ROOT / "vendor" / "bin" / "minimap2"


def build(rng, repeat_len, copies):
    repeat = random_seq(repeat_len, rng)
    spacers = [random_seq(rng.choice([1000, 1400, 1600]), rng) for _ in range(copies + 1)]
    seq = "".join(spacers[i] + repeat for i in range(copies)) + spacers[copies]
    return seq, repeat


def count_copies(asm: Path, repeat_fa: Path, repeat_len: int) -> int:
    out = subprocess.run(
        [str(MM2), "-c", "-x", "asm10", "-N", "50", "-p", "0.3", str(asm), str(repeat_fa)],
        capture_output=True, text=True,
    )
    return sum(1 for line in out.stdout.splitlines()
               if line and int(line.split("\t")[10]) > repeat_len * 0.8)


def main(work: Path):
    work.mkdir(parents=True, exist_ok=True)
    print(f"{'repeat':>7} {'copies':>6} {'construct':>10} {'read mean':>10} "
          f"{'read N50':>9} {'assembled':>10} {'copies found':>13}  verdict")
    print("-" * 92)
    for repeat_len in (2500, 4000):
        rng = random.Random(100 + repeat_len)
        copies = 3
        construct, repeat = build(rng, repeat_len, copies)
        d = work / f"rep{repeat_len}"
        d.mkdir(exist_ok=True)
        write_fasta(d / "truth.fasta", "construct", construct)
        write_fasta(d / "repeat.fasta", "repeat", repeat)
        for mean_len in (1500, 3000, 6000, 12000):
            if mean_len > len(construct):
                continue
            name = f"r{mean_len}"
            import io, contextlib
            with contextlib.redirect_stdout(io.StringIO()):
                fq = simulate(d, name, mean_len=mean_len, depth=200,
                              seed=mean_len, construct=construct)
            cfg = RunConfig(out_dir=d / "out", threads_per_job=4, min_read_len=500,
                            min_read_q=7.0, target_coverage=200)
            s = Sample(name=f"{repeat_len}_{mean_len}", fastqs=[fq], mode=Mode.DENOVO,
                       expected_size=len(construct))
            r = run_sample(s, {}, cfg, d / "out", lambda m: None, threading.Event())
            if not r.ok:
                print(f"{repeat_len:>7} {copies:>6} {len(construct):>10} {mean_len:>10} "
                      f"{'-':>9} {'FAILED':>10} {'-':>13}  {r.message[:30]}")
                continue
            asm = Path(r.outputs["assembly"])
            found = count_copies(asm, d / "repeat.fasta", repeat_len)
            size = r.metrics["assembly_bp"]
            ratio = size / len(construct)
            if found == copies and 0.97 < ratio < 1.03:
                verdict = "correct"
            elif found < copies or ratio < 0.97:
                verdict = f"TRUNCATED (lost {copies - found} copies)"
            else:
                verdict = f"over-assembled ({ratio:.1f}x, circular duplication)"
            print(f"{repeat_len:>7} {copies:>6} {len(construct):>10} {mean_len:>10} "
                  f"{r.metrics['read_n50']:>9} {size:>10} {found:>13}  {verdict}")


if __name__ == "__main__":
    main(Path(sys.argv[1]))
