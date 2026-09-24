"""Per-sample analysis pipeline.

Each sample runs independently: QC -> (assembly | mapping) -> consensus /
variants -> summary. Every external command is logged verbatim into the
sample's run_info.json so a run can be reproduced or audited later.
"""
from __future__ import annotations

import json
import os
import re
import shlex
import shutil
import statistics
import subprocess
import tempfile
import threading
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Optional

from . import constructqc, medaka_models, tools
from .fastq import prepare_reads, read_fasta_lengths, scan_basecall_info
from .models import Mode, Reference, RunConfig, Sample

LogFn = Callable[[str], None]


class Cancelled(Exception):
    """Raised when the user stops a run."""


@dataclass
class StepResult:
    name: str
    ok: bool
    detail: str = ""


@dataclass
class SampleResult:
    sample: str
    mode: str
    ok: bool = False
    out_dir: Optional[Path] = None
    message: str = ""
    steps: list[StepResult] = field(default_factory=list)
    outputs: dict[str, str] = field(default_factory=dict)
    metrics: dict[str, object] = field(default_factory=dict)
    warnings: list[str] = field(default_factory=list)
    duration_s: float = 0.0


class Runner:
    """Runs external commands with logging, cancellation and timing."""

    def __init__(self, log: LogFn, cancel: threading.Event, workdir: Path,
                 threads: int = 1):
        self.log = log
        self.cancel = cancel
        self.workdir = workdir
        self.threads = threads
        self.commands: list[str] = []
        self._proc: Optional[subprocess.Popen] = None

    def check(self) -> None:
        if self.cancel.is_set():
            raise Cancelled()

    def run(
        self,
        argv: list[str],
        stdout_path: Optional[Path] = None,
        allow_fail: bool = False,
        label: str = "",
    ) -> subprocess.CompletedProcess:
        self.check()
        pretty = " ".join(shlex.quote(str(a)) for a in argv)
        if stdout_path:
            pretty += f" > {shlex.quote(str(stdout_path))}"
        self.commands.append(pretty)
        self.log(f"$ {pretty}")
        out = open(stdout_path, "wb") if stdout_path else subprocess.PIPE
        try:
            proc = subprocess.Popen(
                [str(a) for a in argv],
                stdout=out,
                stderr=subprocess.PIPE,
                cwd=str(self.workdir),
                env=tools.tool_env(self.threads),
            )
            self._proc = proc
            _stdout, stderr = proc.communicate()
        finally:
            if stdout_path:
                out.close()
            self._proc = None
        if self.cancel.is_set():
            raise Cancelled()
        err = (stderr or b"").decode("utf-8", "replace").strip()
        if proc.returncode != 0:
            tail = "\n".join(err.splitlines()[-12:])
            msg = f"{label or argv[0]} failed (exit {proc.returncode})\n{tail}"
            if allow_fail:
                self.log(f"  ! {msg}")
            else:
                raise RuntimeError(msg)
        elif err:
            for line in err.splitlines()[-4:]:
                self.log(f"  {line}")
        return subprocess.CompletedProcess(argv, proc.returncode, b"", stderr or b"")

    def pipe_to_sorted_bam(self, argv: list[str], bam: Path, threads: int) -> None:
        """minimap2 | samtools sort -o bam, without a giant intermediate SAM."""
        self.check()
        samtools = tools.require("samtools")
        sort_argv = [str(samtools), "sort", "-@", str(max(1, threads - 1)), "-o", str(bam), "-"]
        pretty = (
            " ".join(shlex.quote(str(a)) for a in argv)
            + " | "
            + " ".join(shlex.quote(str(a)) for a in sort_argv)
        )
        self.commands.append(pretty)
        self.log(f"$ {pretty}")
        env = tools.tool_env(self.threads)
        p1 = subprocess.Popen(
            [str(a) for a in argv],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            cwd=str(self.workdir),
            env=env,
        )
        p2 = subprocess.Popen(
            sort_argv,
            stdin=p1.stdout,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.PIPE,
            cwd=str(self.workdir),
            env=env,
        )
        if p1.stdout:
            p1.stdout.close()
        err2 = p2.communicate()[1]
        err1 = p1.stderr.read() if p1.stderr else b""
        p1.wait()
        if self.cancel.is_set():
            raise Cancelled()
        if p1.returncode != 0:
            raise RuntimeError(f"minimap2 failed: {err1.decode('utf-8','replace')[-500:]}")
        if p2.returncode != 0:
            raise RuntimeError(f"samtools sort failed: {err2.decode('utf-8','replace')[-500:]}")


# --------------------------------------------------------------------------- #
# helpers
# --------------------------------------------------------------------------- #

def safe_dirname(name: str) -> str:
    """Folder name for a sample: whitespace and shell-hostile characters become _.

    The sample keeps its real name everywhere it is displayed; only the folder
    is renamed. Spaces in paths break Flye outright and trip up most shell
    tooling people run on the results afterwards.
    """
    cleaned = re.sub(r"[^\w.+-]+", "_", name.strip()).strip("._")
    return cleaned or "sample"


def unique_dirnames(names: list[str]) -> dict[str, str]:
    """Map sample names to distinct safe folder names ("A B" and "A_B" collide)."""
    taken: set[str] = set()
    result: dict[str, str] = {}
    for name in names:
        base = safe_dirname(name)
        candidate, i = base, 2
        while candidate.lower() in taken:  # macOS filesystems are case-insensitive
            candidate = f"{base}_{i}"
            i += 1
        taken.add(candidate.lower())
        result[name] = candidate
    return result


def choose_flye_mode(basecall_model: str) -> str:
    """Pick Flye's read-type flag from the basecalling model.

    Dorado sup/hac models from v4 onward are accurate enough for --nano-hq;
    anything older or 'fast' gets the tolerant --nano-raw error model.
    """
    m = (basecall_model or "").lower()
    if "sup" in m:
        return "--nano-hq"
    if "hac" in m:
        version = 0.0
        if "@v" in m:
            try:
                version = float(".".join(m.split("@v")[1].split(".")[:2]))
            except ValueError:
                version = 0.0
        return "--nano-hq" if version >= 4.0 else "--nano-raw"
    return "--nano-raw"


def fasta_stats(path: Path) -> dict[str, object]:
    lengths = read_fasta_lengths(path)
    vals = sorted(lengths.values(), reverse=True)
    if not vals:
        return {"contigs": 0, "total_bp": 0, "n50": 0, "largest_bp": 0}
    half = sum(vals) / 2
    acc = 0
    n50 = vals[-1]
    for v in vals:
        acc += v
        if acc >= half:
            n50 = v
            break
    return {
        "contigs": len(vals),
        "total_bp": sum(vals),
        "n50": n50,
        "largest_bp": vals[0],
        "names": list(lengths)[:20],
    }


def depth_summary(depth_tsv: Path, min_depth: int) -> dict[str, object]:
    """Mean/median depth and the fraction of the reference covered."""
    depths: list[int] = []
    with open(depth_tsv) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 3:
                try:
                    depths.append(int(parts[2]))
                except ValueError:
                    continue
    if not depths:
        return {"mean_depth": 0.0, "median_depth": 0, "breadth_1x": 0.0, "breadth_min": 0.0,
                "ref_bp": 0}
    n = len(depths)
    return {
        "ref_bp": n,
        "mean_depth": round(sum(depths) / n, 1),
        "median_depth": int(statistics.median(depths)),
        "breadth_1x": round(100.0 * sum(1 for d in depths if d >= 1) / n, 2),
        "breadth_min": round(100.0 * sum(1 for d in depths if d >= min_depth) / n, 2),
        "min_depth_threshold": min_depth,
    }


def count_variants(vcf: Path) -> int:
    n = 0
    with open(vcf) as fh:
        for line in fh:
            if line and not line.startswith("#"):
                n += 1
    return n


def medaka_ready(cfg) -> bool:
    """Medaka is usable when it is installed and the user has not disabled it."""
    return bool(cfg.use_medaka and tools.have("medaka"))


def medaka_infer_and_stitch(
    runner: Runner,
    bam: Path,
    draft: Path,
    out_path: Path,
    model: str,
    threads: int,
    kind: str = "consensus",
    min_depth: int = 0,
    fill_char: str = "",
) -> bool:
    """Run Medaka's inference then stitch the result into FASTA or VCF.

    Medaka's medaka_consensus/medaka_variant shell wrappers do not quote their
    arguments, so they fail on any path containing a space. Calling the Python
    subcommands directly avoids that, and reuses the BAM we already made
    instead of realigning the reads.
    """
    medaka = tools.require("medaka")
    hdf = runner.workdir / f"medaka_{kind}.hdf"
    hdf.unlink(missing_ok=True)
    proc = runner.run(
        [medaka, "inference", str(bam), str(hdf), "--model", model,
         "--bam_workers", str(max(1, threads // 2))],
        allow_fail=True, label="medaka inference",
    )
    if proc.returncode != 0 or not hdf.exists():
        return False
    if kind == "variant":
        argv = [medaka, "vcf", str(hdf), str(draft), str(out_path)]
    else:
        argv = [medaka, "sequence", str(hdf), str(draft), str(out_path),
                "--threads", str(threads)]
        if min_depth:
            argv += ["--min_depth", str(min_depth)]
        if fill_char:
            argv += ["--fill_char", fill_char]
    proc = runner.run(argv, allow_fail=True, label=f"medaka {kind}")
    return proc.returncode == 0 and out_path.exists() and out_path.stat().st_size > 0


def bcftools_ont_profile(basecall_model: str, available: bool) -> Optional[str]:
    """Pick the bcftools error profile that matches the basecalling model."""
    if not available:
        return None
    return "ont-sup" if "sup" in (basecall_model or "").lower() else "ont"


def bcftools_has_profiles(bcftools: Path) -> bool:
    """Older bcftools builds have no -X platform profiles."""
    try:
        proc = subprocess.run([str(bcftools), "mpileup", "-X", "list"],
                              capture_output=True, timeout=30,
                              env=tools.tool_env())
    except (OSError, subprocess.SubprocessError):
        return False
    return b"ont" in proc.stdout + proc.stderr


# --------------------------------------------------------------------------- #
# main entry point
# --------------------------------------------------------------------------- #

def run_sample(
    sample: Sample,
    references: dict[str, Reference],
    cfg: RunConfig,
    run_dir: Path,
    log: LogFn,
    cancel: threading.Event,
    dir_name: str = "",
) -> SampleResult:
    started = time.time()
    result = SampleResult(sample=sample.name, mode=sample.mode.label)
    out_dir = Path(run_dir) / (dir_name or safe_dirname(sample.name))
    out_dir.mkdir(parents=True, exist_ok=True)
    result.out_dir = out_dir
    work = out_dir / "work"
    work.mkdir(exist_ok=True)
    runner = Runner(log, cancel, work, cfg.threads_per_job)
    info: dict[str, object] = {
        "sample": sample.name,
        "mode": sample.mode.value,
        "inputs": [str(p) for p in sample.fastqs],
        "started": time.strftime("%Y-%m-%d %H:%M:%S"),
    }

    try:
        ref: Optional[Reference] = None
        if sample.mode.needs_reference:
            ref = references.get(sample.reference or "")
            if ref is None:
                raise RuntimeError(f"reference {sample.reference!r} is not loaded")
            info["reference"] = str(ref.path)

        # -- 1. basecalling model ------------------------------------------
        log(f"[{sample.name}] detecting basecalling model from read headers")
        bc = sample.basecall
        if not bc.models:
            bc = scan_basecall_info(sample.fastqs)
            sample.basecall = bc
        if bc.primary:
            log(f"[{sample.name}] basecall model: {bc.display}")
            if bc.is_mixed:
                log(
                    f"[{sample.name}] WARNING: reads from more than one basecalling "
                    f"model ({', '.join(bc.models)}); using the most common one"
                )
        else:
            log(f"[{sample.name}] WARNING: no basecalling model in read headers")
        info["basecall_model"] = bc.primary or ""
        info["basecall_models_seen"] = bc.models
        result.metrics["basecall_model"] = bc.primary or "unknown"
        result.metrics["mixed_models"] = bc.is_mixed

        medaka_bin = tools.status_map().get("medaka")
        medaka_path = str(medaka_bin.path) if medaka_bin and medaka_bin.ok else None
        if sample.medaka_model and sample.medaka_model_source == "manual":
            medaka_model, model_src, model_note = sample.medaka_model, "manual", ""
        else:
            medaka_model, model_src, model_note = medaka_models.resolve(
                bc.primary or "", medaka_path
            )
            sample.medaka_model = medaka_model or ""
            sample.medaka_model_source = model_src
        if model_note:
            log(f"[{sample.name}] {model_note}")
        if medaka_model:
            log(f"[{sample.name}] medaka model: {medaka_model} ({model_src})")
        info["medaka_model"] = medaka_model or ""
        info["medaka_model_source"] = model_src
        result.metrics["medaka_model"] = medaka_model or "none"
        result.steps.append(
            StepResult("basecall model", True, bc.primary or "not found in headers")
        )

        # -- 2. read QC ------------------------------------------------------
        genome_size = sample.expected_size or (ref.length if ref else 0)
        target_bases = cfg.target_coverage * genome_size if (cfg.target_coverage and genome_size) else 0
        reads = work / "reads.fastq.gz"
        log(f"[{sample.name}] filtering reads (len>={cfg.min_read_len}, Q>={cfg.min_read_q})")
        stats = prepare_reads(
            sample.fastqs,
            reads,
            min_len=cfg.min_read_len,
            max_len=cfg.max_read_len,
            min_q=cfg.min_read_q,
            target_bases=target_bases,
            progress=lambda m: log(f"[{sample.name}] {m}"),
        )
        sample.stats = stats
        if stats.n_reads == 0:
            raise RuntimeError(
                "no reads survived filtering - lower the minimum length/quality"
            )
        log(f"[{sample.name}] {stats.display}")
        info["reads"] = {
            "n_reads": stats.n_reads,
            "total_bases": stats.total_bases,
            "n50": stats.n50,
            "mean_q": round(stats.mean_q, 2),
        }
        result.metrics.update(
            reads=stats.n_reads, yield_mb=round(stats.total_bases / 1e6, 2),
            read_n50=stats.n50, mean_q=round(stats.mean_q, 1),
        )
        if genome_size:
            cov = stats.total_bases / genome_size
            result.metrics["est_coverage"] = round(cov, 1)
            log(f"[{sample.name}] estimated coverage: {cov:.0f}x of {genome_size:,} bp")
        result.steps.append(StepResult("read QC", True, stats.display))

        # -- 3. mode-specific work -------------------------------------------
        if sample.mode.does_assembly:
            _assemble(sample, cfg, runner, work, out_dir, reads, genome_size, bc.primary or "",
                      medaka_model, result, log)
        if sample.mode is Mode.DENOVO_SCAFFOLD and ref is not None:
            _scaffold(cfg, runner, work, out_dir, ref, result, log, sample.name)
        if sample.mode.does_mapping and ref is not None:
            _map_and_call(sample, cfg, runner, work, out_dir, reads, ref, medaka_model,
                          result, log)

        result.ok = True
        result.message = "completed"
    except Cancelled:
        result.ok = False
        result.message = "cancelled"
        log(f"[{sample.name}] cancelled")
    except Exception as exc:  # surfaced in the UI, not swallowed
        result.ok = False
        result.message = str(exc)
        log(f"[{sample.name}] ERROR: {exc}")
        result.steps.append(StepResult("failed", False, str(exc)))
    finally:
        result.duration_s = time.time() - started
        info["commands"] = runner.commands
        info["metrics"] = result.metrics
        info["outputs"] = result.outputs
        info["warnings"] = result.warnings
        info["ok"] = result.ok
        info["message"] = result.message
        info["duration_s"] = round(result.duration_s, 1)
        info["tool_versions"] = {
            s.spec.name: s.version for s in tools.probe_all() if s.ok
        }
        (out_dir / "run_info.json").write_text(json.dumps(info, indent=2, default=str))
        if not cfg.keep_intermediates and work.exists():
            shutil.rmtree(work, ignore_errors=True)
    return result


def _assemble(sample, cfg, runner, work, out_dir, reads, genome_size, basecall_model,
              medaka_model, result, log) -> None:
    flye = tools.require("flye")
    asm_dir = work / "flye"
    read_flag = choose_flye_mode(basecall_model)
    log(f"[{sample.name}] de novo assembly with Flye ({read_flag})")
    with _space_free(reads, asm_dir, log, sample.name) as (flye_reads, flye_out):
        argv = [flye, read_flag, str(flye_reads), "--out-dir", str(flye_out),
                "--threads", str(cfg.threads_per_job)]
        if genome_size:
            argv += ["--genome-size", str(genome_size)]
        runner.run(argv, label="flye")
    assembly = asm_dir / "assembly.fasta"
    if not assembly.exists() or assembly.stat().st_size == 0:
        raise RuntimeError("Flye produced no contigs - check coverage and read quality")
    st = fasta_stats(assembly)
    log(f"[{sample.name}] assembly: {st['contigs']} contigs, {st['total_bp']:,} bp, "
        f"N50 {st['n50']:,}")

    info_file = asm_dir / "assembly_info.txt"
    circular = 0
    if info_file.exists():
        for line in info_file.read_text().splitlines()[1:]:
            cols = line.split("\t")
            if len(cols) > 3 and cols[3].strip().lower() in ("y", "yes", "+"):
                circular += 1
        shutil.copy(info_file, out_dir / "assembly_info.txt")
    result.metrics.update(contigs=st["contigs"], assembly_bp=st["total_bp"],
                          assembly_n50=st["n50"], circular_contigs=circular)

    polished = _polish(sample, cfg, runner, work, reads, assembly, medaka_model, log)
    final = out_dir / "assembly.fasta"
    shutil.copy(polished, final)
    result.outputs["assembly"] = str(final)
    result.steps.append(
        StepResult("de novo assembly", True,
                   f"{st['contigs']} contigs, {st['total_bp']:,} bp, {circular} circular")
    )
    _construct_qc(sample, cfg, runner, work, out_dir, reads, final, genome_size,
                  result, log)


def _construct_qc(sample, cfg, runner, work, out_dir, reads, assembly, expected_size,
                  result, log) -> None:
    """Check an assembly for collapsed repeats and doubled circles.

    Reads are mapped back onto the finished assembly: a cassette that is present
    several times but assembled once shows up as a depth spike, which is the one
    piece of evidence a repeat-graph assembler cannot hide.
    """
    log(f"[{sample.name}] checking assembly structure (repeats, size, coverage)")
    bam = work / "assembly_depth.bam"
    depth_tsv = work / "assembly_depth.tsv"
    try:
        runner.pipe_to_sorted_bam(
            [tools.require("minimap2"), "-ax", "map-ont", "-t",
             str(cfg.threads_per_job), str(assembly), str(reads)],
            bam, cfg.threads_per_job,
        )
        runner.run([tools.require("samtools"), "index", str(bam)], label="samtools index")
        runner.run([tools.require("samtools"), "depth", "-a", "-J", str(bam)],
                   stdout_path=depth_tsv, label="samtools depth")
    except RuntimeError as exc:
        log(f"[{sample.name}] could not map reads back for structural QC: {exc}")
        depth_tsv = None

    qc = constructqc.analyse(
        assembly, depth_tsv, reads, expected_size, work,
        flye_info=out_dir / "assembly_info.txt",
        edge_margin=sample.stats.n50,
    )
    result.metrics.update(qc.metrics)
    result.warnings.extend(qc.warnings)
    if depth_tsv and depth_tsv.exists():
        shutil.copy(depth_tsv, out_dir / "assembly_depth.tsv")
        result.outputs["assembly_depth"] = str(out_dir / "assembly_depth.tsv")
    if qc.dedup_path:
        result.outputs["deduplicated_assembly"] = str(qc.dedup_path)
    for warning in qc.warnings:
        log(f"[{sample.name}] WARNING: {warning}")
    if qc.warnings:
        result.steps.append(
            StepResult("structural QC", False,
                       f"{len(qc.warnings)} issue(s) - see warnings")
        )
    else:
        detail = "size and coverage consistent"
        if qc.metrics.get("full_length_reads"):
            detail += f", {qc.metrics['full_length_reads']:,} full-length reads"
        log(f"[{sample.name}] structural QC clean: {detail}")
        result.steps.append(StepResult("structural QC", True, detail))


class _space_free:
    """Give Flye paths without spaces, which it refuses outright.

    When the reads or output folder contain a space, the reads are symlinked
    into a fresh temporary directory, Flye runs there, and its whole output
    folder is copied back afterwards - including on failure, so flye.log is
    still there to read. Otherwise the real paths are used directly.
    """

    def __init__(self, reads: Path, out_dir: Path, log: LogFn, name: str):
        self.reads, self.out_dir, self.log, self.name = reads, out_dir, log, name
        self.stage: Optional[Path] = None

    def __enter__(self) -> tuple[Path, Path]:
        if " " not in str(self.reads) and " " not in str(self.out_dir):
            return self.reads, self.out_dir
        self.stage = Path(tempfile.mkdtemp(prefix="nanoassembler_flye_"))
        if " " in str(self.stage):
            raise RuntimeError(
                f"Flye cannot run: the output path contains spaces and so does the "
                f"temporary directory ({self.stage}). Choose an output folder "
                f"without spaces."
            )
        staged = self.stage / "reads.fastq.gz"
        os.symlink(self.reads.resolve(), staged)
        self.log(f"[{self.name}] path contains spaces, which Flye rejects - "
                 f"running it in {self.stage}")
        return staged, self.stage / "flye"

    def __exit__(self, *_exc) -> None:
        if self.stage is None:
            return
        try:
            produced = self.stage / "flye"
            if produced.exists():
                shutil.copytree(produced, self.out_dir, dirs_exist_ok=True)
        finally:
            shutil.rmtree(self.stage, ignore_errors=True)


def _polish(sample, cfg, runner, work, reads, draft: Path, medaka_model, log) -> Path:
    """Polish a draft assembly with Medaka if we can, else Racon, else leave it."""
    if medaka_ready(cfg) and medaka_model:
        log(f"[{sample.name}] polishing with Medaka ({medaka_model})")
        bam = work / "polish.bam"
        runner.pipe_to_sorted_bam(
            [tools.require("minimap2"), "-ax", "map-ont", "-t",
             str(cfg.threads_per_job), str(draft), str(reads)],
            bam, cfg.threads_per_job,
        )
        runner.run([tools.require("samtools"), "index", str(bam)],
                   label="samtools index")
        out = work / "medaka_polished.fasta"
        if medaka_infer_and_stitch(runner, bam, draft, out, medaka_model,
                                   cfg.threads_per_job, kind="consensus"):
            return out
        log(f"[{sample.name}] Medaka polishing failed; falling back")
    if tools.have("racon"):
        log(f"[{sample.name}] polishing with Racon")
        paf = work / "polish.paf"
        runner.run([tools.require("minimap2"), "-x", "map-ont", "-t",
                    str(cfg.threads_per_job), str(draft), str(reads)],
                   stdout_path=paf, label="minimap2")
        out = work / "racon.fasta"
        proc = runner.run([tools.require("racon"), "-t", str(cfg.threads_per_job),
                           str(reads), str(paf), str(draft)],
                          stdout_path=out, allow_fail=True, label="racon")
        if proc.returncode == 0 and out.exists() and out.stat().st_size:
            return out
    log(f"[{sample.name}] no polisher available - returning the unpolished draft")
    return draft


def _scaffold(cfg, runner, work, out_dir, ref, result, log, sample_name) -> None:
    if not tools.have("ragtag.py"):
        log(f"[{sample_name}] RagTag not installed - skipping reference ordering")
        result.steps.append(StepResult("scaffold to reference", False, "RagTag not installed"))
        return
    asm = out_dir / "assembly.fasta"
    sc_dir = work / "ragtag"
    log(f"[{sample_name}] ordering contigs against {ref.name}")
    runner.run([tools.require("ragtag.py"), "scaffold", str(ref.path), str(asm),
                "-o", str(sc_dir), "-t", str(cfg.threads_per_job)], label="ragtag")
    scaffolded = sc_dir / "ragtag.scaffold.fasta"
    if scaffolded.exists():
        dest = out_dir / "assembly.scaffolded.fasta"
        shutil.copy(scaffolded, dest)
        result.outputs["scaffolded_assembly"] = str(dest)
        st = fasta_stats(dest)
        result.metrics["scaffolds"] = st["contigs"]
        result.steps.append(
            StepResult("scaffold to reference", True, f"{st['contigs']} scaffolds")
        )


def _map_and_call(sample, cfg, runner, work, out_dir, reads, ref, medaka_model,
                  result, log) -> None:
    minimap2 = tools.require("minimap2")
    samtools = tools.require("samtools")
    bam = out_dir / "aligned.bam"
    log(f"[{sample.name}] mapping to {ref.name}")
    runner.pipe_to_sorted_bam(
        [minimap2, "-ax", "map-ont", "--MD", "-t", str(cfg.threads_per_job),
         str(ref.path), str(reads)],
        bam, cfg.threads_per_job,
    )
    runner.run([samtools, "index", str(bam)], label="samtools index")
    result.outputs["alignment"] = str(bam)

    flagstat = out_dir / "mapping_stats.txt"
    runner.run([samtools, "flagstat", str(bam)], stdout_path=flagstat, label="samtools flagstat")
    mapped_pct = _parse_flagstat_mapped(flagstat)
    if mapped_pct is not None:
        result.metrics["mapped_pct"] = mapped_pct
        log(f"[{sample.name}] {mapped_pct:.1f}% of reads mapped")

    depth_tsv = work / "depth.tsv"
    runner.run([samtools, "depth", "-a", "-J", str(bam)], stdout_path=depth_tsv,
               label="samtools depth")
    cov = depth_summary(depth_tsv, cfg.min_consensus_depth)
    result.metrics.update(cov)
    shutil.copy(depth_tsv, out_dir / "depth.tsv")
    result.outputs["depth"] = str(out_dir / "depth.tsv")
    log(f"[{sample.name}] mean depth {cov['mean_depth']}x, "
        f"{cov['breadth_min']}% of reference at >={cfg.min_consensus_depth}x")
    result.steps.append(
        StepResult("mapping", True,
                   f"{cov['mean_depth']}x mean depth, {cov['breadth_1x']}% covered")
    )

    # Depth against the expected map: a region at several times the median
    # carries more copies than the reference says, a dropout is a deletion.
    high, low, median = constructqc.find_anomalies(depth_tsv, sample.stats.n50)
    if high:
        result.metrics["extra_copy_regions"] = len(high)
        detail = "; ".join(str(r) for r in high[:4])
        result.warnings.append(
            f"{len(high)} region(s) of {ref.name} are covered far above the median "
            f"({median:.0f}x) - the sample appears to carry extra copies there: {detail}"
        )
    if low:
        result.metrics["dropout_regions"] = len(low)
        result.warnings.append(
            f"{len(low)} region(s) of {ref.name} have almost no coverage - possible "
            f"deletion relative to the reference: {'; '.join(str(r) for r in low[:3])}"
        )
    for warning in result.warnings[-2:]:
        if warning:
            log(f"[{sample.name}] WARNING: {warning}")

    if sample.mode.does_consensus:
        _consensus(sample, cfg, runner, work, out_dir, reads, ref, bam, medaka_model,
                   result, log)
    if sample.mode.does_variants:
        _variants(sample, cfg, runner, out_dir, ref, bam, result, log)


def _consensus(sample, cfg, runner, work, out_dir, reads, ref, bam, medaka_model,
               result, log) -> None:
    """Reference-guided consensus: Medaka when available, samtools otherwise."""
    dest = out_dir / "consensus.fasta"
    caller = ""
    if medaka_ready(cfg) and medaka_model:
        log(f"[{sample.name}] reference consensus with Medaka ({medaka_model})")
        if medaka_infer_and_stitch(
            runner, bam, ref.path, dest, medaka_model, cfg.threads_per_job,
            kind="consensus", min_depth=cfg.min_consensus_depth, fill_char="N",
        ):
            caller = f"medaka {medaka_model}"
        else:
            log(f"[{sample.name}] Medaka consensus failed; using samtools consensus")
    if not caller:
        samtools = tools.require("samtools")
        log(f"[{sample.name}] reference consensus with samtools "
            f"(min depth {cfg.min_consensus_depth}, low-coverage bases masked as N)")
        runner.run(
            [samtools, "consensus", "-a", "--show-del", "yes", "--mode", "simple",
             "--min-depth", str(cfg.min_consensus_depth), "--call-fract", "0.5",
             "-o", str(dest), str(bam)],
            label="samtools consensus",
        )
        caller = "samtools consensus"

    st = fasta_stats(dest)
    ns = _count_ns(dest)
    result.metrics["consensus_bp"] = st["total_bp"]
    result.metrics["consensus_n_pct"] = (
        round(100.0 * ns / st["total_bp"], 2) if st["total_bp"] else 0.0
    )
    result.metrics["consensus_caller"] = caller
    result.outputs["consensus"] = str(dest)
    log(f"[{sample.name}] consensus: {st['total_bp']:,} bp, "
        f"{result.metrics['consensus_n_pct']}% N ({caller})")
    result.steps.append(
        StepResult("consensus", True,
                   f"{st['total_bp']:,} bp, {result.metrics['consensus_n_pct']}% N, {caller}")
    )


def _variants(sample, cfg, runner, out_dir, ref, bam, result, log) -> None:
    """Call variants with Medaka's variant network when we have one, else bcftools."""
    vcf = out_dir / "variants.vcf"
    if medaka_ready(cfg) and sample.basecall.primary:
        var_model, _src, note = medaka_models.resolve(
            sample.basecall.primary, str(tools.require("medaka")), kind="variant"
        )
        if var_model:
            if note:
                log(f"[{sample.name}] {note}")
            log(f"[{sample.name}] calling variants with Medaka ({var_model})")
            if medaka_infer_and_stitch(runner, bam, ref.path, vcf, var_model,
                                       cfg.threads_per_job, kind="variant"):
                n = count_variants(vcf)
                result.metrics.update(variants=n, variant_caller=f"medaka {var_model}")
                result.outputs["variants"] = str(vcf)
                log(f"[{sample.name}] {n} variants (medaka {var_model})")
                result.steps.append(
                    StepResult("variant calling", True, f"{n} variants, medaka")
                )
                return
            log(f"[{sample.name}] Medaka variant calling failed; using bcftools")

    bcftools = tools.require("bcftools")
    log(f"[{sample.name}] calling variants against {ref.name} with bcftools")
    raw = out_dir / "variants.raw.bcf"
    mpileup = [bcftools, "mpileup", "-f", str(ref.path), "-d", "1000",
               "-q", "10", "-Q", "5", "-a", "AD,DP", "-Ou", "-o", str(raw)]
    profile = bcftools_ont_profile(sample.basecall.primary or "",
                                   bcftools_has_profiles(bcftools))
    if profile:
        mpileup[2:2] = ["-X", profile]
        log(f"[{sample.name}] using bcftools error profile '{profile}'")
    mpileup.append(str(bam))
    runner.run(mpileup, label="bcftools mpileup")
    called = out_dir / "variants.called.bcf"
    normed = out_dir / "variants.norm.bcf"
    runner.run([bcftools, "call", "-mv", "--ploidy", "1", "-Ou", "-o", str(called),
                str(raw)], label="bcftools call")
    runner.run([bcftools, "norm", "-f", str(ref.path), "-Ou", "-o", str(normed),
                str(called)], label="bcftools norm")
    runner.run([bcftools, "filter",
                "-i", f"QUAL>=20 && INFO/DP>={cfg.min_consensus_depth}",
                "-Ov", "-o", str(vcf), str(normed)], label="bcftools filter")
    for tmp in (raw, called, normed):
        tmp.unlink(missing_ok=True)
    n = count_variants(vcf)
    result.metrics.update(variants=n, variant_caller="bcftools")
    result.outputs["variants"] = str(vcf)
    log(f"[{sample.name}] {n} variants passing QUAL>=20 and DP>={cfg.min_consensus_depth}")
    result.steps.append(StepResult("variant calling", True, f"{n} variants, bcftools"))


def _parse_flagstat_mapped(path: Path) -> Optional[float]:
    try:
        for line in path.read_text().splitlines():
            if " mapped (" in line and "primary" not in line:
                pct = line.split("(")[1].split("%")[0]
                return float(pct)
    except (OSError, IndexError, ValueError):
        return None
    return None


def _count_ns(fasta: Path) -> int:
    n = 0
    with open(fasta) as fh:
        for line in fh:
            if not line.startswith(">"):
                n += line.upper().count("N")
    return n
