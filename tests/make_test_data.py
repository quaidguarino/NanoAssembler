"""Generate a small synthetic ONT dataset for end-to-end testing.

Writes a reference FASTA, a mutated 'true' sequence, and FASTQ files whose
headers carry Dorado-style basecall_model_version_id tags.
"""
from __future__ import annotations

import gzip
import random
import sys
from pathlib import Path

BASES = "ACGT"


def random_seq(n: int, rng: random.Random) -> str:
    return "".join(rng.choice(BASES) for _ in range(n))


def mutate(seq: str, n_snps: int, rng: random.Random) -> tuple[str, list[int]]:
    s = list(seq)
    positions = rng.sample(range(100, len(seq) - 100), n_snps)
    for p in positions:
        s[p] = rng.choice([b for b in BASES if b != s[p]])
    return "".join(s), sorted(positions)


def add_errors(seq: str, rng: random.Random, sub=0.015, ins=0.008, dele=0.012) -> str:
    out = []
    for base in seq:
        r = rng.random()
        if r < dele:
            continue
        if r < dele + ins:
            out.append(base)
            out.append(rng.choice(BASES))
            continue
        if r < dele + ins + sub:
            out.append(rng.choice([b for b in BASES if b != base]))
            continue
        out.append(base)
    return "".join(out)


def write_fasta(path: Path, name: str, seq: str, width: int = 60) -> None:
    with open(path, "w") as fh:
        fh.write(f">{name}\n")
        for i in range(0, len(seq), width):
            fh.write(seq[i : i + width] + "\n")


def simulate(out_dir: Path, seed: int = 7, genome: int = 20000, depth: int = 120,
             model: str = "dna_r10.4.1_e8.2_400bps_sup@v5.0.0",
             sample: str = "sample01", n_snps: int = 12) -> dict:
    rng = random.Random(seed)
    out_dir.mkdir(parents=True, exist_ok=True)
    ref = random_seq(genome, rng)
    truth, snp_pos = mutate(ref, n_snps, rng)
    ref_path = out_dir / f"{sample}_reference.fasta"
    write_fasta(ref_path, "plasmid_ref", ref)
    write_fasta(out_dir / f"{sample}_truth.fasta", "plasmid_truth", truth)

    fq = out_dir / f"{sample}.fastq.gz"
    target = genome * depth
    written = 0
    n = 0
    runid = "".join(rng.choice("0123456789abcdef") for _ in range(20))
    with gzip.open(fq, "wt") as fh:
        while written < target:
            length = min(genome, max(1000, int(rng.gauss(6000, 2500))))
            start = rng.randrange(0, max(1, genome - length))
            frag = truth[start : start + length]
            if rng.random() < 0.5:  # reverse strand
                comp = {"A": "T", "C": "G", "G": "C", "T": "A"}
                frag = "".join(comp[b] for b in reversed(frag))
            read = add_errors(frag, rng)
            if len(read) < 500:
                continue
            qual = "".join(chr(33 + max(5, min(40, int(rng.gauss(18, 4))))) for _ in read)
            rid = "".join(rng.choice("0123456789abcdef") for _ in range(8))
            fh.write(
                f"@{rid}-{n} runid={runid} read={n} ch={rng.randrange(1,512)} "
                f"start_time=2026-01-01T00:00:00Z "
                f"basecall_model_version_id={model}\n{read}\n+\n{qual}\n"
            )
            written += len(read)
            n += 1
    return {"reference": ref_path, "fastq": fq, "reads": n, "snps": snp_pos,
            "genome": genome, "model": model, "sample": sample}


if __name__ == "__main__":
    out = Path(sys.argv[1] if len(sys.argv) > 1 else "testdata")
    meta = simulate(out)
    print(f"{meta['reads']} reads -> {meta['fastq']}")
    print(f"reference: {meta['reference']}  ({meta['genome']} bp, {len(meta['snps'])} SNPs)")
