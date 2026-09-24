"""Structural QC for assembled constructs.

Repeat-graph assemblers collapse repeats they cannot span, which shows up as a
short assembly even though the molecule is intact. The reverse also happens:
when reads approach the length of a circular molecule, Flye can emit the circle
twice. Neither failure announces itself, so both are detected here from
evidence that does not depend on the assembler:

  * read depth mapped back onto the assembly - a cassette present three times
    but assembled once carries roughly three times the median depth;
  * self-alignment of the contig - a doubled circle aligns to its own other half;
  * the assembled size against the size you expect;
  * how many single reads span the whole construct, which settles whether a
    truncation is real or an artefact.
"""
from __future__ import annotations

import statistics
import subprocess
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

from . import tools
from .fastq import iter_records, read_fasta_lengths

WINDOW = 200          # bp, depth is averaged over windows this size
HIGH_RATIO = 1.7      # depth/median above which a region looks multi-copy
LOW_RATIO = 0.35      # depth/median below which a region looks unsupported
MIN_REGION = 300      # bp, ignore shorter anomalies


@dataclass
class Region:
    contig: str
    start: int
    end: int
    ratio: float

    @property
    def length(self) -> int:
        return self.end - self.start

    @property
    def est_copies(self) -> int:
        return max(2, round(self.ratio))

    def __str__(self) -> str:
        return (f"{self.contig}:{self.start:,}-{self.end:,} "
                f"({self.length:,} bp at {self.ratio:.1f}x median)")


@dataclass
class ConstructQC:
    warnings: list[str] = field(default_factory=list)
    metrics: dict[str, object] = field(default_factory=dict)
    high_regions: list[Region] = field(default_factory=list)
    low_regions: list[Region] = field(default_factory=list)
    dedup_path: Optional[Path] = None


def read_depths(depth_tsv: Path) -> dict[str, list[int]]:
    per_contig: dict[str, list[int]] = {}
    with open(depth_tsv) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 3:
                try:
                    per_contig.setdefault(parts[0], []).append(int(parts[2]))
                except ValueError:
                    continue
    return per_contig


def find_anomalies(
    depth_tsv: Path, edge_margin: int = 0
) -> tuple[list[Region], list[Region], float]:
    """Windows whose depth departs from the assembly-wide median.

    ``edge_margin`` suppresses low-coverage calls within that distance of a
    contig end. Depth always tapers over roughly one read length at the ends of
    a linear sequence because no read can extend past it, so without this every
    sample reports two spurious "deletions". Pass the read N50. High-coverage
    calls are never suppressed - a spike at an end is real.
    """
    per_contig = read_depths(depth_tsv)
    all_depths = [d for v in per_contig.values() for d in v]
    if not all_depths:
        return [], [], 0.0
    median = statistics.median(all_depths)
    if median <= 0:
        return [], [], 0.0

    high: list[Region] = []
    low: list[Region] = []
    for contig, depths in per_contig.items():
        for bucket, target in ((high, "high"), (low, "low")):
            run_start = None
            run_vals: list[float] = []
            for i in range(0, len(depths), WINDOW):
                window = depths[i:i + WINDOW]
                ratio = (sum(window) / len(window)) / median
                hit = ratio >= HIGH_RATIO if target == "high" else ratio <= LOW_RATIO
                if hit:
                    if run_start is None:
                        run_start = i
                        run_vals = []
                    run_vals.append(ratio)
                elif run_start is not None:
                    _emit(bucket, contig, run_start, i, run_vals, len(depths),
                          edge_margin if target == "low" else 0)
                    run_start = None
            if run_start is not None:
                _emit(bucket, contig, run_start, len(depths), run_vals, len(depths),
                      edge_margin if target == "low" else 0)
    return high, low, median


def _emit(bucket: list[Region], contig: str, start: int, end: int,
          ratios: list[float], contig_len: int, edge_margin: int) -> None:
    if end - start < MIN_REGION or not ratios:
        return
    if edge_margin:
        margin = min(edge_margin, contig_len // 4)
        if start < margin or end > contig_len - margin:
            return
    bucket.append(Region(contig, start, end, sum(ratios) / len(ratios)))


def detect_self_duplication(assembly: Path, workdir: Path) -> dict[str, float]:
    """Fraction of each contig's first half that aligns to its second half.

    A circular molecule assembled twice end-to-end aligns to itself almost
    perfectly at exactly half its length.
    """
    minimap2 = tools.find_tool("minimap2")
    if minimap2 is None:
        return {}
    seqs: dict[str, list[str]] = {}
    name = None
    with open(assembly) as fh:
        for line in fh:
            if line.startswith(">"):
                name = line[1:].split()[0]
                seqs[name] = []
            elif name:
                seqs[name].append(line.strip())
    result: dict[str, float] = {}
    for contig, chunks in seqs.items():
        seq = "".join(chunks)
        half = len(seq) // 2
        if half < 500:
            continue
        a = workdir / "selfdup_a.fasta"
        b = workdir / "selfdup_b.fasta"
        a.write_text(f">a\n{seq[:half]}\n")
        b.write_text(f">b\n{seq[half:]}\n")
        proc = subprocess.run(
            [str(minimap2), "-c", "-x", "asm10", str(a), str(b)],
            capture_output=True, text=True, env=tools.tool_env(),
        )
        covered = 0
        for line in proc.stdout.splitlines():
            cols = line.split("\t")
            if len(cols) > 10:
                try:
                    covered += int(cols[10])
                except ValueError:
                    continue
        a.unlink(missing_ok=True)
        b.unlink(missing_ok=True)
        result[contig] = covered / half if half else 0.0
    return result


def write_deduplicated(assembly: Path, dest: Path, duplicated: set[str]) -> None:
    """Write a copy with each doubled contig cut back to a single circle."""
    seqs: list[tuple[str, str]] = []
    name, chunks = None, []
    with open(assembly) as fh:
        for line in fh:
            if line.startswith(">"):
                if name:
                    seqs.append((name, "".join(chunks)))
                name, chunks = line[1:].rstrip("\n"), []
            elif name:
                chunks.append(line.strip())
    if name:
        seqs.append((name, "".join(chunks)))
    with open(dest, "w") as out:
        for header, seq in seqs:
            short = header.split()[0]
            if short in duplicated:
                seq = seq[: len(seq) // 2]
                header = f"{header} deduplicated_from_doubled_circle"
            out.write(f">{header}\n")
            for i in range(0, len(seq), 60):
                out.write(seq[i:i + 60] + "\n")


def count_spanning_reads(fastq: Path, expected_size: int, fraction: float = 0.95) -> int:
    """Reads long enough to cover the whole construct on their own."""
    if not expected_size:
        return 0
    threshold = int(expected_size * fraction)
    return sum(1 for _h, seq, _q in iter_records(fastq) if len(seq) >= threshold)


def analyse(
    assembly: Path,
    depth_tsv: Optional[Path],
    reads: Path,
    expected_size: int,
    workdir: Path,
    flye_info: Optional[Path] = None,
    edge_margin: int = 0,
) -> ConstructQC:
    qc = ConstructQC()
    lengths = read_fasta_lengths(assembly)
    total = sum(lengths.values())
    qc.metrics["assembly_bp"] = total

    # --- size against expectation -------------------------------------
    if expected_size:
        ratio = total / expected_size
        qc.metrics["size_vs_expected"] = round(ratio, 3)
        if ratio < 0.95:
            qc.warnings.append(
                f"Assembly is {total:,} bp, {100 * (1 - ratio):.0f}% SHORT of the "
                f"expected {expected_size:,} bp. Repeated features that the reads "
                f"cannot span are the usual cause."
            )
        elif ratio > 1.5:
            qc.warnings.append(
                f"Assembly is {total:,} bp, {ratio:.1f}x the expected "
                f"{expected_size:,} bp - check for a doubled circular contig."
            )

    # --- single reads spanning the whole construct ---------------------
    # Recorded for every sample, but only worth warning about when something
    # else already looks wrong: a large construct legitimately has no
    # full-length read and still assembles correctly.
    spanning = count_spanning_reads(reads, expected_size)
    qc.metrics["full_length_reads"] = spanning
    short_assembly = bool(expected_size) and total / expected_size < 0.95

    # --- doubled circle -------------------------------------------------
    dup = detect_self_duplication(assembly, workdir)
    duplicated = {c for c, frac in dup.items() if frac > 0.9}
    qc.metrics["duplicated_contigs"] = len(duplicated)
    if duplicated:
        qc.warnings.append(
            f"Contig(s) {', '.join(sorted(duplicated))} align to their own second "
            f"half: the circular molecule has been assembled twice end to end. A "
            f"corrected copy is written as assembly.deduplicated.fasta."
        )
        dest = assembly.parent / "assembly.deduplicated.fasta"
        write_deduplicated(assembly, dest, duplicated)
        qc.dedup_path = dest

    # --- collapsed repeats from read depth ------------------------------
    if depth_tsv and depth_tsv.exists():
        high, low, median = find_anomalies(depth_tsv, edge_margin)
        qc.high_regions, qc.low_regions = high, low
        qc.metrics["median_depth"] = round(median, 1)
        if high:
            copies = sum(r.est_copies - 1 for r in high)
            qc.metrics["collapsed_repeat_regions"] = len(high)
            qc.metrics["est_missing_copies"] = copies
            detail = "; ".join(str(r) for r in high[:4])
            qc.warnings.append(
                f"{len(high)} region(s) carry well above the median read depth, "
                f"which is what a collapsed repeat looks like - roughly {copies} "
                f"extra copy(ies) present in the reads but not in the assembly: "
                f"{detail}"
            )
        if low:
            qc.metrics["low_support_regions"] = len(low)
            qc.warnings.append(
                f"{len(low)} region(s) have very low read support and may be "
                f"misassembled: {'; '.join(str(r) for r in low[:3])}"
            )

    # --- can a read settle it? -------------------------------------------
    suspect = short_assembly or bool(qc.high_regions)
    if expected_size and suspect:
        if spanning:
            qc.warnings.append(
                f"{spanning:,} read(s) are long enough to cover the expected "
                f"{expected_size:,} bp on their own - inspect those reads directly "
                f"to decide whether the construct is genuinely truncated."
            )
        else:
            qc.warnings.append(
                f"No single read reaches the expected {expected_size:,} bp, so "
                f"nothing here can distinguish a collapsed repeat from a real "
                f"deletion. A longer-fragment prep would settle it."
            )

    # --- what Flye itself thought ---------------------------------------
    if flye_info and flye_info.exists():
        multi = []
        for line in flye_info.read_text().splitlines()[1:]:
            cols = line.split("\t")
            if len(cols) >= 6:
                try:
                    mult = int(cols[5])
                except ValueError:
                    continue
                if mult > 1 or cols[4].strip().upper() == "Y":
                    multi.append(f"{cols[0]} (multiplicity {mult})")
        if multi:
            qc.metrics["flye_repeat_contigs"] = len(multi)
            qc.warnings.append(
                "Flye marked these contigs as repetitive: " + ", ".join(multi[:5])
            )
    return qc
