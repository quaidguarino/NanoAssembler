"""End-to-end check: run every mode against the synthetic dataset."""
from __future__ import annotations

import sys
import threading
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from nanoassembler.core.models import Mode, Reference, RunConfig, Sample
from nanoassembler.core.pipeline import run_sample
from nanoassembler.core.report import write_report

data = Path(sys.argv[1])
out = Path(sys.argv[2])
only = sys.argv[3:] or None

ref = Reference(path=data / "sample01_reference.fasta", name="plasmid_ref")
from nanoassembler.core.fastq import read_fasta_lengths
ref.length = sum(read_fasta_lengths(ref.path).values())
refs = {ref.name: ref}

cfg = RunConfig(out_dir=out, threads_per_job=4, parallel_jobs=1,
                min_read_len=500, min_read_q=7.0, target_coverage=150,
                min_consensus_depth=10)

modes = [Mode.DENOVO, Mode.REF_CONSENSUS, Mode.REF_VARIANTS, Mode.DENOVO_SCAFFOLD]
if only:
    modes = [Mode(m) for m in only]

results = []
for mode in modes:
    s = Sample(name=f"test_{mode.value}", fastqs=[data / "sample01.fastq.gz"],
               mode=mode, reference=ref.name if mode.needs_reference else None)
    r = run_sample(s, refs, cfg, out, lambda m: print(m, flush=True), threading.Event())
    results.append(r)
    print(f"### {mode.value}: {'OK' if r.ok else 'FAILED - ' + r.message}", flush=True)
    print(f"### metrics: {r.metrics}\n", flush=True)

print("report:", write_report(out, results, cfg))
print("SUMMARY:", [(r.sample, r.ok) for r in results])
