"""Data model for samples, references and run configuration."""
from __future__ import annotations

import json
from dataclasses import dataclass, field, asdict
from enum import Enum
from pathlib import Path
from typing import Optional


class Mode(str, Enum):
    """What to do with a sample."""

    DENOVO = "denovo"
    REF_CONSENSUS = "ref_consensus"
    REF_VARIANTS = "ref_variants"
    REF_CONSENSUS_VARIANTS = "ref_consensus_variants"
    DENOVO_SCAFFOLD = "denovo_scaffold"

    @property
    def label(self) -> str:
        return {
            Mode.DENOVO: "De novo assembly",
            Mode.REF_CONSENSUS: "Reference: consensus",
            Mode.REF_VARIANTS: "Reference: variants",
            Mode.REF_CONSENSUS_VARIANTS: "Reference: consensus + variants",
            Mode.DENOVO_SCAFFOLD: "De novo + scaffold to reference",
        }[self]

    @property
    def needs_reference(self) -> bool:
        return self is not Mode.DENOVO

    @property
    def does_assembly(self) -> bool:
        return self in (Mode.DENOVO, Mode.DENOVO_SCAFFOLD)

    @property
    def does_mapping(self) -> bool:
        return self in (
            Mode.REF_CONSENSUS,
            Mode.REF_VARIANTS,
            Mode.REF_CONSENSUS_VARIANTS,
        )

    @property
    def does_consensus(self) -> bool:
        return self in (Mode.REF_CONSENSUS, Mode.REF_CONSENSUS_VARIANTS)

    @property
    def does_variants(self) -> bool:
        return self in (Mode.REF_VARIANTS, Mode.REF_CONSENSUS_VARIANTS)

    @classmethod
    def from_label(cls, label: str) -> "Mode":
        for m in cls:
            if m.label == label:
                return m
        raise ValueError(f"unknown mode label: {label!r}")


@dataclass
class Reference:
    """A reference FASTA registered in the run."""

    path: Path
    name: str = ""
    length: int = 0
    n_seqs: int = 0

    def __post_init__(self) -> None:
        self.path = Path(self.path)
        if not self.name:
            self.name = self.path.stem


@dataclass
class BasecallInfo:
    """Basecalling model(s) detected in a FASTQ's read headers."""

    models: dict[str, int] = field(default_factory=dict)  # model id -> read count
    reads_scanned: int = 0
    run_ids: set[str] = field(default_factory=set)

    @property
    def primary(self) -> Optional[str]:
        if not self.models:
            return None
        return max(self.models.items(), key=lambda kv: kv[1])[0]

    @property
    def is_mixed(self) -> bool:
        return len(self.models) > 1

    @property
    def display(self) -> str:
        """Short form for the table; ``primary`` keeps the full id."""
        if not self.models:
            return "not in headers"
        short = (self.primary or "").removeprefix("dna_")
        if self.is_mixed:
            return f"{short}  (+{len(self.models) - 1} other)"
        return short


@dataclass
class ReadStats:
    """Summary statistics for a sample's reads."""

    n_reads: int = 0
    total_bases: int = 0
    n50: int = 0
    median_len: int = 0
    max_len: int = 0
    mean_q: float = 0.0

    @property
    def display(self) -> str:
        if not self.n_reads:
            return ""
        mb = self.total_bases / 1e6
        return f"{self.n_reads:,} reads / {mb:.1f} Mb / N50 {self.n50:,} / Q{self.mean_q:.1f}"


@dataclass
class Sample:
    """One sample: one or more FASTQ files, a mode, and optionally a reference."""

    name: str
    fastqs: list[Path] = field(default_factory=list)
    mode: Mode = Mode.DENOVO
    reference: Optional[str] = None  # Reference.name
    expected_size: int = 0  # bp; 0 = infer from reference or auto
    basecall: BasecallInfo = field(default_factory=BasecallInfo)
    stats: ReadStats = field(default_factory=ReadStats)
    medaka_model: str = ""  # resolved; "" = none/auto-fallback
    medaka_model_source: str = ""  # "auto" | "manual" | "fallback"
    enabled: bool = True

    def __post_init__(self) -> None:
        self.fastqs = [Path(p) for p in self.fastqs]

    def validate(self, ref_names: set[str]) -> list[str]:
        problems: list[str] = []
        if not self.fastqs:
            problems.append(f"{self.name}: no FASTQ files")
        for f in self.fastqs:
            if not f.exists():
                problems.append(f"{self.name}: missing file {f.name}")
        if self.mode.needs_reference:
            if not self.reference:
                problems.append(f"{self.name}: mode '{self.mode.label}' needs a reference")
            elif self.reference not in ref_names:
                problems.append(f"{self.name}: reference '{self.reference}' is not loaded")
        return problems


@dataclass
class RunConfig:
    """Global settings for a run."""

    out_dir: Path = Path.home() / "NanoAssembler_runs"
    threads_per_job: int = 4
    parallel_jobs: int = 2
    min_read_len: int = 200
    max_read_len: int = 0  # 0 = no cap
    min_read_q: float = 9.0
    target_coverage: int = 200  # 0 = no subsampling
    keep_intermediates: bool = True
    use_medaka: bool = True
    min_consensus_depth: int = 10

    def __post_init__(self) -> None:
        self.out_dir = Path(self.out_dir)

    def to_json(self) -> str:
        d = asdict(self)
        d["out_dir"] = str(self.out_dir)
        return json.dumps(d, indent=2)
