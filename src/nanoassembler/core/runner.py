"""Background orchestration: run many samples concurrently off the GUI thread."""
from __future__ import annotations

import threading
import time
from concurrent.futures import Future, ThreadPoolExecutor
from pathlib import Path

from PySide6.QtCore import QObject, QRunnable, QThreadPool, Signal

from .models import Reference, RunConfig, Sample
from .pipeline import SampleResult, run_sample, unique_dirnames


class RunSignals(QObject):
    log = Signal(str)
    sample_started = Signal(str)
    sample_finished = Signal(object)  # SampleResult
    progress = Signal(int, int)  # done, total
    finished = Signal(object)  # list[SampleResult]


class RunJob(QRunnable):
    """Executes a whole run; individual samples go to a small thread pool."""

    def __init__(
        self,
        samples: list[Sample],
        references: dict[str, Reference],
        cfg: RunConfig,
        run_dir: Path,
    ):
        super().__init__()
        self.samples = samples
        self.references = references
        self.cfg = cfg
        self.run_dir = Path(run_dir)
        self.signals = RunSignals()
        self.cancel = threading.Event()
        self._lock = threading.Lock()
        self._done = 0
        self.dir_names = unique_dirnames([s.name for s in samples])

    def stop(self) -> None:
        self.cancel.set()

    def _log(self, message: str) -> None:
        self.signals.log.emit(message)

    def _one(self, sample: Sample) -> SampleResult:
        self.signals.sample_started.emit(sample.name)
        self._log(f"--- {sample.name}: {sample.mode.label} ---")
        result = run_sample(
            sample, self.references, self.cfg, self.run_dir, self._log, self.cancel,
            dir_name=self.dir_names[sample.name],
        )
        with self._lock:
            self._done += 1
            done = self._done
        self.signals.sample_finished.emit(result)
        self.signals.progress.emit(done, len(self.samples))
        status = "OK" if result.ok else f"FAILED: {result.message}"
        self._log(f"--- {sample.name}: {status} ({result.duration_s:.0f}s) ---")
        return result

    def run(self) -> None:
        started = time.time()
        self.run_dir.mkdir(parents=True, exist_ok=True)
        (self.run_dir / "run_config.json").write_text(self.cfg.to_json())
        self._log(f"Run directory: {self.run_dir}")
        self._log(
            f"{len(self.samples)} sample(s), {self.cfg.parallel_jobs} at a time, "
            f"{self.cfg.threads_per_job} threads each"
        )
        results: list[SampleResult] = []
        workers = max(1, min(self.cfg.parallel_jobs, len(self.samples)))
        with ThreadPoolExecutor(max_workers=workers) as pool:
            futures: list[Future] = [pool.submit(self._one, s) for s in self.samples]
            for fut in futures:
                try:
                    results.append(fut.result())
                except Exception as exc:  # pragma: no cover - defensive
                    self._log(f"internal error: {exc}")
        elapsed = time.time() - started
        n_ok = sum(1 for r in results if r.ok)
        self._log(
            f"Run finished in {elapsed/60:.1f} min - {n_ok}/{len(results)} succeeded"
        )
        self.signals.finished.emit(results)


def start_run(job: RunJob) -> None:
    QThreadPool.globalInstance().start(job)
