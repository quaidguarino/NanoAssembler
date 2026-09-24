"""Simulate a plasmid carrying several identical cassettes (e.g. repeated EF1A
promoters) to measure how de novo assembly handles the repeats.

Layout (9,000 bp, circular):
    uniqueA 1000 | REPEAT 1200 | uniqueB 1600 | REPEAT 1200 |
    uniqueC 1600 | REPEAT 1200 | uniqueD 1200
The three REPEAT copies are byte-identical, which is the worst case for a
repeat-graph assembler.
"""
from __future__ import annotations

import gzip
import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from make_test_data import add_errors, random_seq, write_fasta

MODEL = "dna_r10.4.1_e8.2_400bps_sup@v5.0.0"


def build_construct(rng: random.Random, repeat_len: int = 1200, copies: int = 3):
    repeat = random_seq(repeat_len, rng)
    uniques = [random_seq(n, rng) for n in (1000, 1600, 1600, 1200)]
    seq = ""
    for i in range(copies):
        seq += uniques[i] + repeat
    seq += uniques[copies]
    return seq, repeat


def simulate(out_dir: Path, name: str, mean_len: int, depth: int, seed: int,
             construct: str) -> Path:
    """Circular plasmid: reads start anywhere and wrap around the origin."""
    rng = random.Random(seed)
    doubled = construct * 2
    genome = len(construct)
    fq = out_dir / f"{name}.fastq.gz"
    written, n = 0, 0
    comp = {"A": "T", "C": "G", "G": "C", "T": "A"}
    with gzip.open(fq, "wt") as fh:
        while written < genome * depth:
            length = int(rng.gauss(mean_len, mean_len * 0.45))
            length = max(500, min(length, genome))  # no concatemers
            start = rng.randrange(genome)
            frag = doubled[start:start + length]
            if rng.random() < 0.5:
                frag = "".join(comp[b] for b in reversed(frag))
            read = add_errors(frag, rng)
            qual = "".join(chr(33 + max(5, min(40, int(rng.gauss(18, 4))))) for _ in read)
            rid = "".join(rng.choice("0123456789abcdef") for _ in range(8))
            fh.write(f"@{rid}-{n} runid=repeattest read={n} ch={rng.randrange(1,512)} "
                     f"basecall_model_version_id={MODEL}\n{read}\n+\n{qual}\n")
            written += len(read)
            n += 1
    print(f"  {name}: {n} reads, mean {written//n} bp, {written/genome:.0f}x")
    return fq


if __name__ == "__main__":
    out = Path(sys.argv[1])
    out.mkdir(parents=True, exist_ok=True)
    rng = random.Random(42)
    construct, repeat = build_construct(rng)
    write_fasta(out / "construct_truth.fasta", "construct_9kb_3xEF1A", construct)
    write_fasta(out / "construct_reference.fasta", "construct_9kb_3xEF1A", construct)
    write_fasta(out / "repeat_unit.fasta", "EF1A_like_repeat", repeat)
    print(f"construct: {len(construct)} bp with 3 identical {len(repeat)} bp cassettes")
    simulate(out, "longreads", mean_len=6000, depth=200, seed=1, construct=construct)
    simulate(out, "midreads", mean_len=3000, depth=200, seed=2, construct=construct)
    simulate(out, "shortreads", mean_len=1500, depth=200, seed=3, construct=construct)
