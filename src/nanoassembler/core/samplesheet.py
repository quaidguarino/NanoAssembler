"""Import/export the sample table as CSV so runs can be reproduced."""
from __future__ import annotations

import csv
from pathlib import Path

from .models import Mode, Reference, Sample

HEADER = ["sample", "fastq", "mode", "reference", "expected_size"]


def write(path: Path, samples: list[Sample], references: dict[str, Reference]) -> None:
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(HEADER)
        for s in samples:
            ref = references.get(s.reference or "")
            w.writerow([
                s.name,
                ";".join(str(p) for p in s.fastqs),
                s.mode.value,
                str(ref.path) if ref else "",
                s.expected_size or "",
            ])


def read(path: Path) -> tuple[list[Sample], list[Reference], list[str]]:
    """Load a sample sheet. Returns (samples, references, warnings)."""
    samples: list[Sample] = []
    refs: dict[str, Reference] = {}
    warnings: list[str] = []
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh)
        missing = [c for c in ("sample", "fastq", "mode") if c not in (reader.fieldnames or [])]
        if missing:
            raise ValueError(f"sample sheet is missing column(s): {', '.join(missing)}")
        for i, row in enumerate(reader, start=2):
            name = (row.get("sample") or "").strip()
            if not name:
                continue
            fastqs = [Path(p.strip()) for p in (row.get("fastq") or "").split(";") if p.strip()]
            for f in fastqs:
                if not f.exists():
                    warnings.append(f"line {i}: {f} does not exist")
            try:
                mode = Mode(row.get("mode", "").strip() or "denovo")
            except ValueError:
                warnings.append(f"line {i}: unknown mode {row.get('mode')!r}, using de novo")
                mode = Mode.DENOVO
            ref_name = None
            ref_path = (row.get("reference") or "").strip()
            if ref_path:
                rp = Path(ref_path)
                if rp.exists():
                    ref = refs.get(rp.stem) or Reference(path=rp)
                    refs[ref.name] = ref
                    ref_name = ref.name
                else:
                    warnings.append(f"line {i}: reference {rp} does not exist")
            size = 0
            raw_size = (row.get("expected_size") or "").strip()
            if raw_size:
                try:
                    size = int(float(raw_size))
                except ValueError:
                    warnings.append(f"line {i}: bad expected_size {raw_size!r}")
            samples.append(
                Sample(name=name, fastqs=fastqs, mode=mode, reference=ref_name,
                       expected_size=size)
            )
    return samples, list(refs.values()), warnings
