"""FASTQ reading: basecall-model detection, QC filtering, subsampling."""
from __future__ import annotations

import gzip
import math
import re
from pathlib import Path
from typing import Callable, Iterator, Optional

from .models import BasecallInfo, ReadStats

# Dorado / Guppy write the basecalling model into the FASTQ header comment as a
# key=value tag; when FASTQ is produced from BAM with `samtools fastq -T RG`
# it instead shows up inside the read-group string.
_MODEL_TAGS = (
    re.compile(r"\bbasecall_model_version_id=(\S+)"),
    re.compile(r"\bmodel_version_id=(\S+)"),
    re.compile(r"\bmodel=(\S+)"),
)
# Last resort: find a model-shaped token anywhere in the header (this is what
# lets us recover the model from RG:Z:<runid>_<model>_<barcode> strings).
_MODEL_SHAPE = re.compile(
    r"((?:dna|rna)[_a-z0-9.]*?_?"
    r"(?:\d+bps)?_?(?:fast|hac|sup)"
    r"(?:@v[\d.]+)?)",
    re.IGNORECASE,
)
_RUNID_TAGS = (
    re.compile(r"\brunid=(\S+)"),
    re.compile(r"\bRD:Z:(\S+)"),
)


def open_text(path: Path):
    """Open a plain or gzipped FASTQ for text reading."""
    path = Path(path)
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "rt", encoding="utf-8", errors="replace")


def iter_records(path: Path) -> Iterator[tuple[str, str, str]]:
    """Yield (header_without_@, sequence, qualities) from a FASTQ file."""
    with open_text(path) as fh:
        while True:
            header = fh.readline()
            if not header:
                return
            header = header.rstrip("\n")
            seq = fh.readline().rstrip("\n")
            fh.readline()  # '+'
            qual = fh.readline().rstrip("\n")
            if not qual:
                return
            if not header.startswith("@"):
                raise ValueError(f"{path.name}: malformed FASTQ near {header[:60]!r}")
            yield header[1:], seq, qual


def extract_model(header: str) -> Optional[str]:
    """Pull the basecalling model id out of a single FASTQ header line."""
    for pat in _MODEL_TAGS:
        m = pat.search(header)
        if m:
            return m.group(1)
    m = _MODEL_SHAPE.search(header)
    if m:
        return m.group(1)
    return None


def extract_runid(header: str) -> Optional[str]:
    for pat in _RUNID_TAGS:
        m = pat.search(header)
        if m:
            return m.group(1)
    return None


def scan_basecall_info(paths: list[Path], max_reads: int = 4000) -> BasecallInfo:
    """Detect the basecalling model(s) used, by sampling read headers.

    Every FASTQ belonging to the sample is sampled so that a merged file set
    basecalled with different models is reported as mixed rather than silently
    taking the first one.
    """
    info = BasecallInfo()
    if not paths:
        return info
    per_file = max(1, max_reads // len(paths))
    for path in paths:
        n = 0
        try:
            for header, _seq, _qual in iter_records(Path(path)):
                model = extract_model(header)
                if model:
                    info.models[model] = info.models.get(model, 0) + 1
                run = extract_runid(header)
                if run:
                    info.run_ids.add(run)
                n += 1
                info.reads_scanned += 1
                if n >= per_file:
                    break
        except (OSError, ValueError):
            continue
    return info


def mean_quality(qual: str) -> float:
    """Mean read quality as Phred of the mean error probability (ONT convention)."""
    if not qual:
        return 0.0
    err = 0.0
    for ch in qual:
        err += 10.0 ** (-(ord(ch) - 33) / 10.0)
    err /= len(qual)
    if err <= 0:
        return 60.0
    return min(60.0, -10.0 * math.log10(err))


def _n50(lengths: list[int]) -> int:
    if not lengths:
        return 0
    ordered = sorted(lengths, reverse=True)
    half = sum(ordered) / 2
    acc = 0
    for length in ordered:
        acc += length
        if acc >= half:
            return length
    return ordered[-1]


def prepare_reads(
    paths: list[Path],
    out_path: Path,
    min_len: int = 0,
    max_len: int = 0,
    min_q: float = 0.0,
    target_bases: int = 0,
    progress: Optional[Callable[[str], None]] = None,
) -> ReadStats:
    """Merge, filter and (optionally) subsample reads into ``out_path``.

    Subsampling keeps the longest reads until ``target_bases`` is reached, which
    is what you want for assembly and harmless for consensus calling.
    Returns statistics for the reads that were actually written.
    """
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    tmp = out_path.with_suffix(out_path.suffix + ".pass1")

    lengths: list[int] = []
    quals: list[float] = []
    dropped = 0
    with gzip.open(tmp, "wt", compresslevel=1) as out:
        for path in paths:
            for header, seq, qual in iter_records(Path(path)):
                length = len(seq)
                if length < min_len or (max_len and length > max_len):
                    dropped += 1
                    continue
                q = mean_quality(qual)
                if q < min_q:
                    dropped += 1
                    continue
                lengths.append(length)
                quals.append(q)
                out.write(f"@{header}\n{seq}\n+\n{qual}\n")
    if progress:
        progress(f"filter: kept {len(lengths):,} reads, dropped {dropped:,}")

    total = sum(lengths)
    if target_bases and total > target_bases:
        cutoff = _pick_length_cutoff(lengths, target_bases)
        kept_lengths: list[int] = []
        kept_quals: list[float] = []
        written = 0
        with gzip.open(tmp, "rt") as src, gzip.open(out_path, "wt", compresslevel=1) as out:
            while True:
                header = src.readline()
                if not header:
                    break
                seq = src.readline().rstrip("\n")
                src.readline()
                qual = src.readline().rstrip("\n")
                if len(seq) < cutoff or written >= target_bases:
                    continue
                written += len(seq)
                kept_lengths.append(len(seq))
                kept_quals.append(mean_quality(qual))
                out.write(f"{header.rstrip()}\n{seq}\n+\n{qual}\n")
        tmp.unlink(missing_ok=True)
        lengths, quals = kept_lengths, kept_quals
        if progress:
            progress(
                f"subsample: kept {len(lengths):,} reads >= {cutoff:,} bp "
                f"({sum(lengths)/1e6:.1f} Mb)"
            )
    else:
        tmp.replace(out_path)

    if not lengths:
        return ReadStats()
    ordered = sorted(lengths)
    return ReadStats(
        n_reads=len(lengths),
        total_bases=sum(lengths),
        n50=_n50(lengths),
        median_len=ordered[len(ordered) // 2],
        max_len=ordered[-1],
        mean_q=sum(quals) / len(quals),
    )


def _pick_length_cutoff(lengths: list[int], target_bases: int) -> int:
    """Smallest read length such that keeping reads >= it stays near the target."""
    acc = 0
    cutoff = 0
    for length in sorted(lengths, reverse=True):
        acc += length
        cutoff = length
        if acc >= target_bases:
            break
    return cutoff


def read_fasta_lengths(path: Path) -> dict[str, int]:
    """Sequence name -> length for a (optionally gzipped) FASTA."""
    lengths: dict[str, int] = {}
    name = None
    count = 0
    with open_text(Path(path)) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    lengths[name] = count
                name = line[1:].split()[0] if len(line) > 1 else "unnamed"
                count = 0
            elif name is not None:
                count += len(line.strip())
    if name is not None:
        lengths[name] = count
    return lengths
