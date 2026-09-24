"""Main application window: sample table, references, settings, run control."""
from __future__ import annotations

import os
import subprocess
import time
from pathlib import Path
from typing import Optional

from PySide6.QtCore import QObject, QRunnable, Qt, QThreadPool, Signal, Slot
from PySide6.QtGui import QAction, QColor, QFont, QKeySequence
from PySide6.QtWidgets import (
    QAbstractItemView, QCheckBox, QComboBox, QDoubleSpinBox, QFileDialog,
    QFormLayout, QGroupBox, QHBoxLayout, QHeaderView, QLabel, QLineEdit,
    QListWidget, QListWidgetItem, QMainWindow, QMessageBox, QPlainTextEdit,
    QProgressBar, QPushButton, QSpinBox, QSplitter, QStatusBar, QTabWidget,
    QTableWidget, QTableWidgetItem, QVBoxLayout, QWidget,
)

from .. import APP_NAME, __version__
from ..core import medaka_models, report, tools
from ..core.fastq import read_fasta_lengths, scan_basecall_info
from ..core.models import Mode, Reference, RunConfig, Sample
from ..core.runner import RunJob
from ..core.samplesheet import read as read_sheet, write as write_sheet

FASTQ_SUFFIXES = (".fastq", ".fq", ".fastq.gz", ".fq.gz")
FASTA_SUFFIXES = (".fa", ".fasta", ".fna", ".fa.gz", ".fasta.gz", ".fna.gz")

COL_ON, COL_NAME, COL_FILES, COL_MODE, COL_REF, COL_SIZE, COL_MODEL, \
    COL_MEDAKA, COL_READS, COL_STATUS = range(10)
HEADERS = ["", "Sample", "FASTQ", "Mode", "Reference", "Expected size",
           "Basecall model", "Medaka model", "Reads", "Status"]


def sample_name_from(path: Path) -> str:
    name = path.name
    for suf in (".gz", ".fastq", ".fq"):
        if name.endswith(suf):
            name = name[: -len(suf)]
    return name or path.stem


def is_fastq(path: Path) -> bool:
    return path.name.lower().endswith(FASTQ_SUFFIXES)


def is_fasta(path: Path) -> bool:
    return path.name.lower().endswith(FASTA_SUFFIXES)


class ScanSignals(QObject):
    done = Signal(str, object, object)  # sample name, BasecallInfo, medaka resolution


class ScanJob(QRunnable):
    """Read the head of a sample's FASTQs to detect the basecalling model."""

    def __init__(self, name: str, paths: list[Path], medaka_bin: Optional[str]):
        super().__init__()
        self.name = name
        self.paths = paths
        self.medaka_bin = medaka_bin
        self.signals = ScanSignals()

    def run(self) -> None:
        info = scan_basecall_info(self.paths)
        resolution = medaka_models.resolve(info.primary or "", self.medaka_bin)
        self.signals.done.emit(self.name, info, resolution)


class MainWindow(QMainWindow):
    def __init__(self) -> None:
        super().__init__()
        self.setWindowTitle(f"{APP_NAME} {__version__}")
        self.resize(1440, 860)
        self.setAcceptDrops(True)

        self.samples: list[Sample] = []
        self.references: dict[str, Reference] = {}
        self.cfg = RunConfig()
        self.job: Optional[RunJob] = None
        self.results: list = []
        self.run_dir: Optional[Path] = None
        self._medaka_choices: list[str] = []

        self._build_ui()
        self._build_menu()
        self._refresh_tool_status()

    # ------------------------------------------------------------------ UI
    def _build_ui(self) -> None:
        splitter = QSplitter(Qt.Orientation.Vertical)

        top = QSplitter(Qt.Orientation.Horizontal)
        top.addWidget(self._build_sample_panel())
        top.addWidget(self._build_side_panel())
        top.setStretchFactor(0, 4)
        top.setStretchFactor(1, 1)
        top.setSizes([1100, 340])

        splitter.addWidget(top)
        splitter.addWidget(self._build_bottom_tabs())
        splitter.setStretchFactor(0, 3)
        splitter.setStretchFactor(1, 2)
        splitter.setSizes([520, 300])

        container = QWidget()
        layout = QVBoxLayout(container)
        layout.setContentsMargins(10, 10, 10, 8)
        layout.addWidget(splitter)
        layout.addLayout(self._build_run_bar())
        self.setCentralWidget(container)

        self.status = QStatusBar()
        self.setStatusBar(self.status)

    def _build_sample_panel(self) -> QWidget:
        box = QWidget()
        layout = QVBoxLayout(box)
        layout.setContentsMargins(0, 0, 0, 0)

        bar = QHBoxLayout()
        for label, slot in (
            ("Add FASTQ files…", self.add_fastq_files),
            ("Add folder…", self.add_fastq_folder),
            ("Merge selected", self.merge_selected),
            ("Remove selected", self.remove_selected),
        ):
            btn = QPushButton(label)
            btn.clicked.connect(slot)
            bar.addWidget(btn)
        bar.addStretch(1)
        self.apply_all_btn = QPushButton("Apply mode/reference to all…")
        self.apply_all_btn.clicked.connect(self.apply_to_all)
        bar.addWidget(self.apply_all_btn)
        layout.addLayout(bar)

        self.table = QTableWidget(0, len(HEADERS))
        self.table.setHorizontalHeaderLabels(HEADERS)
        self.table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.table.setAlternatingRowColors(True)
        self.table.verticalHeader().setVisible(False)
        header = self.table.horizontalHeader()
        header.setSectionResizeMode(QHeaderView.ResizeMode.Interactive)
        header.setSectionResizeMode(COL_ON, QHeaderView.ResizeMode.Fixed)
        self.table.setColumnWidth(COL_ON, 28)
        header.setStretchLastSection(True)
        for col, width in ((COL_NAME, 125), (COL_FILES, 65), (COL_MODE, 195),
                           (COL_REF, 125), (COL_SIZE, 95), (COL_MODEL, 195),
                           (COL_MEDAKA, 175), (COL_READS, 165)):
            self.table.setColumnWidth(col, width)
        self.table.itemChanged.connect(self._on_item_changed)
        layout.addWidget(self.table)

        hint = QLabel(
            "Drag FASTQ files here to add samples, or FASTA files to add references. "
            "One FASTQ = one sample; select several rows and use Merge to combine them."
        )
        hint.setStyleSheet("color: palette(mid); padding: 2px;")
        layout.addWidget(hint)
        return box

    def _build_side_panel(self) -> QWidget:
        panel = QWidget()
        layout = QVBoxLayout(panel)
        layout.setContentsMargins(0, 0, 0, 0)

        ref_box = QGroupBox("References")
        ref_layout = QVBoxLayout(ref_box)
        self.ref_list = QListWidget()
        self.ref_list.setToolTip("FASTA references available to the sample table")
        ref_layout.addWidget(self.ref_list)
        row = QHBoxLayout()
        add_ref = QPushButton("Add FASTA…")
        add_ref.clicked.connect(self.add_references)
        drop_ref = QPushButton("Remove")
        drop_ref.clicked.connect(self.remove_reference)
        row.addWidget(add_ref)
        row.addWidget(drop_ref)
        ref_layout.addLayout(row)
        layout.addWidget(ref_box, 2)

        layout.addWidget(self._build_settings_box(), 3)

        self.tools_box = QGroupBox("Toolchain")
        tl = QVBoxLayout(self.tools_box)
        self.tools_label = QLabel("checking…")
        self.tools_label.setWordWrap(True)
        self.tools_label.setTextFormat(Qt.TextFormat.RichText)
        tl.addWidget(self.tools_label)
        recheck = QPushButton("Re-check")
        recheck.clicked.connect(lambda: self._refresh_tool_status(refresh=True))
        tl.addWidget(recheck)
        layout.addWidget(self.tools_box, 1)
        return panel

    def _build_settings_box(self) -> QGroupBox:
        box = QGroupBox("Run settings")
        form = QFormLayout(box)
        form.setLabelAlignment(Qt.AlignmentFlag.AlignRight)

        out_row = QHBoxLayout()
        self.out_edit = QLineEdit(str(self.cfg.out_dir))
        browse = QPushButton("…")
        browse.setFixedWidth(30)
        browse.clicked.connect(self.pick_out_dir)
        out_row.addWidget(self.out_edit)
        out_row.addWidget(browse)
        form.addRow("Output folder", out_row)

        self.parallel_spin = QSpinBox()
        self.parallel_spin.setRange(1, 16)
        self.parallel_spin.setValue(self.cfg.parallel_jobs)
        self.parallel_spin.setToolTip("How many samples to process at the same time")
        form.addRow("Samples in parallel", self.parallel_spin)

        self.threads_spin = QSpinBox()
        self.threads_spin.setRange(1, max(1, os.cpu_count() or 4))
        self.threads_spin.setValue(self.cfg.threads_per_job)
        self.threads_spin.setToolTip("Threads given to each sample's tools")
        form.addRow("Threads per sample", self.threads_spin)

        self.minlen_spin = QSpinBox()
        self.minlen_spin.setRange(0, 100000)
        self.minlen_spin.setSingleStep(50)
        self.minlen_spin.setValue(self.cfg.min_read_len)
        form.addRow("Min read length (bp)", self.minlen_spin)

        self.maxlen_spin = QSpinBox()
        self.maxlen_spin.setRange(0, 10_000_000)
        self.maxlen_spin.setSingleStep(1000)
        self.maxlen_spin.setValue(self.cfg.max_read_len)
        self.maxlen_spin.setSpecialValueText("no limit")
        form.addRow("Max read length (bp)", self.maxlen_spin)

        self.minq_spin = QDoubleSpinBox()
        self.minq_spin.setRange(0, 40)
        self.minq_spin.setSingleStep(0.5)
        self.minq_spin.setValue(self.cfg.min_read_q)
        form.addRow("Min mean read Q", self.minq_spin)

        self.cov_spin = QSpinBox()
        self.cov_spin.setRange(0, 5000)
        self.cov_spin.setSingleStep(25)
        self.cov_spin.setValue(self.cfg.target_coverage)
        self.cov_spin.setSpecialValueText("off")
        self.cov_spin.setToolTip(
            "Subsample the longest reads down to this depth. Amplicons are often "
            "sequenced to thousands of x, which slows assembly without helping it."
        )
        form.addRow("Target coverage", self.cov_spin)

        self.depth_spin = QSpinBox()
        self.depth_spin.setRange(1, 1000)
        self.depth_spin.setValue(self.cfg.min_consensus_depth)
        self.depth_spin.setToolTip("Consensus positions below this depth are masked as N")
        form.addRow("Min consensus depth", self.depth_spin)

        self.medaka_check = QCheckBox("Use Medaka (model from read headers)")
        self.medaka_check.setChecked(self.cfg.use_medaka)
        form.addRow("", self.medaka_check)

        self.keep_check = QCheckBox("Keep intermediate files")
        self.keep_check.setChecked(self.cfg.keep_intermediates)
        form.addRow("", self.keep_check)
        return box

    def _build_bottom_tabs(self) -> QWidget:
        self.tabs = QTabWidget()
        self.log_view = QPlainTextEdit()
        self.log_view.setReadOnly(True)
        self.log_view.setMaximumBlockCount(20000)
        self.log_view.setFont(QFont("Menlo", 11))
        self.tabs.addTab(self.log_view, "Log")

        self.results_table = QTableWidget(0, 6)
        self.results_table.setHorizontalHeaderLabels(
            ["Sample", "Status", "Basecall model", "Medaka model", "Key metrics", "Output"]
        )
        self.results_table.horizontalHeader().setStretchLastSection(True)
        self.results_table.setColumnWidth(0, 160)
        self.results_table.setColumnWidth(2, 240)
        self.results_table.setColumnWidth(3, 230)
        self.results_table.setColumnWidth(4, 330)
        self.results_table.cellDoubleClicked.connect(self._open_result_row)
        self.tabs.addTab(self.results_table, "Results")
        return self.tabs

    def _build_run_bar(self) -> QHBoxLayout:
        bar = QHBoxLayout()
        self.progress = QProgressBar()
        self.progress.setTextVisible(True)
        self.progress.setFormat("%v / %m samples")
        bar.addWidget(self.progress, 1)
        self.report_btn = QPushButton("Open report")
        self.report_btn.setEnabled(False)
        self.report_btn.clicked.connect(self.open_report)
        bar.addWidget(self.report_btn)
        self.stop_btn = QPushButton("Stop")
        self.stop_btn.setEnabled(False)
        self.stop_btn.clicked.connect(self.stop_run)
        bar.addWidget(self.stop_btn)
        self.run_btn = QPushButton("Run")
        self.run_btn.setDefault(True)
        self.run_btn.setMinimumWidth(120)
        self.run_btn.clicked.connect(self.start_run)
        bar.addWidget(self.run_btn)
        return bar

    def _build_menu(self) -> None:
        m = self.menuBar()
        file_menu = m.addMenu("&File")
        for label, slot, shortcut in (
            ("Add FASTQ files…", self.add_fastq_files, QKeySequence.StandardKey.Open),
            ("Add folder of FASTQ…", self.add_fastq_folder, None),
            ("Add reference FASTA…", self.add_references, None),
        ):
            act = QAction(label, self)
            act.triggered.connect(slot)
            if shortcut:
                act.setShortcut(shortcut)
            file_menu.addAction(act)
        file_menu.addSeparator()
        for label, slot in (("Load sample sheet…", self.load_sheet),
                            ("Save sample sheet…", self.save_sheet),
                            ("Open output folder", self.open_out_dir)):
            act = QAction(label, self)
            act.triggered.connect(slot)
            file_menu.addAction(act)

        run_menu = m.addMenu("&Run")
        start = QAction("Start run", self)
        start.setShortcut("Ctrl+R")
        start.triggered.connect(self.start_run)
        run_menu.addAction(start)
        stop = QAction("Stop run", self)
        stop.triggered.connect(self.stop_run)
        run_menu.addAction(stop)

        help_menu = m.addMenu("&Help")
        about = QAction("About", self)
        about.triggered.connect(self.show_about)
        help_menu.addAction(about)

    # ------------------------------------------------------- tool status
    def _refresh_tool_status(self, refresh: bool = False) -> None:
        statuses = tools.probe_all(refresh=refresh)
        lines = []
        for st in statuses:
            if st.ok:
                lines.append(f"<span style='color:#1a7f52'>&#10003;</span> "
                             f"{st.spec.name} {st.version}")
            else:
                colour = "#b3261e" if st.spec.required else "#8a6100"
                lines.append(f"<span style='color:{colour}'>&#10007;</span> "
                             f"{st.spec.name} &mdash; {st.error}")
        self.tools_label.setText("<br>".join(lines))
        missing = tools.missing_required()
        if missing:
            self.status.showMessage(
                f"Missing required tools: {', '.join(missing)} - run scripts/fetch_tools.sh"
            )
        medaka_st = tools.status_map().get("medaka")
        if medaka_st and medaka_st.ok and medaka_st.path:
            self._medaka_choices = list(medaka_models.available_models(str(medaka_st.path)))
        self.medaka_check.setEnabled(bool(medaka_st and medaka_st.ok))
        if not (medaka_st and medaka_st.ok):
            self.medaka_check.setChecked(False)
            self.medaka_check.setToolTip(
                "Medaka is not installed; consensus falls back to samtools and "
                "variants to bcftools."
            )

    def _medaka_bin(self) -> Optional[str]:
        st = tools.status_map().get("medaka")
        return str(st.path) if st and st.ok and st.path else None

    # ---------------------------------------------------------- sample IO
    def add_fastq_files(self) -> None:
        paths, _ = QFileDialog.getOpenFileNames(
            self, "Add FASTQ files", str(Path.home()),
            "FASTQ (*.fastq *.fq *.fastq.gz *.fq.gz);;All files (*)",
        )
        self._add_paths([Path(p) for p in paths])

    def add_fastq_folder(self) -> None:
        folder = QFileDialog.getExistingDirectory(self, "Add folder of FASTQ files",
                                                  str(Path.home()))
        if not folder:
            return
        found = sorted(p for p in Path(folder).rglob("*") if p.is_file() and is_fastq(p))
        if not found:
            QMessageBox.information(self, APP_NAME, "No FASTQ files found in that folder.")
            return
        self._add_paths(found)

    def add_references(self) -> None:
        paths, _ = QFileDialog.getOpenFileNames(
            self, "Add reference FASTA", str(Path.home()),
            "FASTA (*.fa *.fasta *.fna *.fa.gz *.fasta.gz *.fna.gz);;All files (*)",
        )
        self._add_paths([Path(p) for p in paths])

    def _add_paths(self, paths: list[Path]) -> None:
        added_samples = 0
        for path in paths:
            if is_fastq(path):
                name = self._unique_name(sample_name_from(path))
                self.samples.append(Sample(name=name, fastqs=[path]))
                added_samples += 1
            elif is_fasta(path):
                self._add_reference(path)
        if added_samples:
            self._rebuild_table()
            self._scan_new_samples()
        self._sync_reference_list()
        self.status.showMessage(
            f"{len(self.samples)} sample(s), {len(self.references)} reference(s)"
        )

    def _add_reference(self, path: Path) -> None:
        lengths = read_fasta_lengths(path)
        name = path.stem
        base = name
        i = 2
        while name in self.references and self.references[name].path != path:
            name = f"{base}_{i}"
            i += 1
        self.references[name] = Reference(
            path=path, name=name, length=sum(lengths.values()), n_seqs=len(lengths)
        )
        self._refresh_ref_combos()

    def _unique_name(self, name: str) -> str:
        existing = {s.name for s in self.samples}
        if name not in existing:
            return name
        i = 2
        while f"{name}_{i}" in existing:
            i += 1
        return f"{name}_{i}"

    def remove_selected(self) -> None:
        rows = sorted({i.row() for i in self.table.selectedIndexes()}, reverse=True)
        for row in rows:
            if 0 <= row < len(self.samples):
                del self.samples[row]
        self._rebuild_table()

    def merge_selected(self) -> None:
        rows = sorted({i.row() for i in self.table.selectedIndexes()})
        if len(rows) < 2:
            QMessageBox.information(
                self, APP_NAME, "Select two or more rows to merge into one sample."
            )
            return
        first = self.samples[rows[0]]
        for row in rows[1:]:
            first.fastqs.extend(self.samples[row].fastqs)
        for row in reversed(rows[1:]):
            del self.samples[row]
        first.basecall.models.clear()
        self._rebuild_table()
        self._scan_new_samples()

    def remove_reference(self) -> None:
        item = self.ref_list.currentItem()
        if not item:
            return
        name = item.data(Qt.ItemDataRole.UserRole)
        self.references.pop(name, None)
        for s in self.samples:
            if s.reference == name:
                s.reference = None
        self._sync_reference_list()
        self._refresh_ref_combos()

    def apply_to_all(self) -> None:
        rows = sorted({i.row() for i in self.table.selectedIndexes()})
        if not rows:
            QMessageBox.information(
                self, APP_NAME,
                "Select the row whose mode and reference you want copied to every sample."
            )
            return
        src = self.samples[rows[0]]
        for s in self.samples:
            s.mode = src.mode
            s.reference = src.reference
        self._rebuild_table()

    # ------------------------------------------------------------- table
    def _rebuild_table(self) -> None:
        self.table.blockSignals(True)
        self.table.setRowCount(len(self.samples))
        for row, sample in enumerate(self.samples):
            self._fill_row(row, sample)
        self.table.blockSignals(False)

    def _fill_row(self, row: int, sample: Sample) -> None:
        on = QTableWidgetItem()
        on.setFlags(Qt.ItemFlag.ItemIsUserCheckable | Qt.ItemFlag.ItemIsEnabled)
        on.setCheckState(Qt.CheckState.Checked if sample.enabled else Qt.CheckState.Unchecked)
        self.table.setItem(row, COL_ON, on)

        name = QTableWidgetItem(sample.name)
        name.setToolTip("Double-click to rename; this names the output folder")
        self.table.setItem(row, COL_NAME, name)

        files = QTableWidgetItem(
            f"{len(sample.fastqs)} file" + ("s" if len(sample.fastqs) != 1 else "")
        )
        files.setFlags(Qt.ItemFlag.ItemIsEnabled | Qt.ItemFlag.ItemIsSelectable)
        files.setToolTip("\n".join(str(p) for p in sample.fastqs))
        self.table.setItem(row, COL_FILES, files)

        mode_combo = QComboBox()
        for m in Mode:
            # Store the value, not the enum: Qt round-trips item data through
            # QVariant, which turns a str-based enum back into a plain str.
            mode_combo.addItem(m.label, m.value)
        mode_combo.setCurrentIndex(list(Mode).index(sample.mode))
        mode_combo.currentIndexChanged.connect(
            lambda _i, r=row: self._on_mode_changed(r)
        )
        self.table.setCellWidget(row, COL_MODE, mode_combo)

        ref_combo = QComboBox()
        self._populate_ref_combo(ref_combo, sample)
        ref_combo.currentIndexChanged.connect(lambda _i, r=row: self._on_ref_changed(r))
        self.table.setCellWidget(row, COL_REF, ref_combo)

        size = QTableWidgetItem(f"{sample.expected_size:,}" if sample.expected_size else "")
        size.setToolTip(
            "Expected target size in bp. Used to estimate coverage and to decide "
            "how far to subsample. Taken from the reference when one is set; "
            "type a value here for de novo samples."
        )
        self.table.setItem(row, COL_SIZE, size)

        model = QTableWidgetItem(sample.basecall.display or "scanning…")
        model.setFlags(Qt.ItemFlag.ItemIsEnabled | Qt.ItemFlag.ItemIsSelectable)
        if sample.basecall.is_mixed:
            model.setForeground(QColor("#8a6100"))
            model.setToolTip(
                "More than one basecalling model in this sample:\n"
                + "\n".join(f"{k}: {v} reads" for k, v in sample.basecall.models.items())
            )
        self.table.setItem(row, COL_MODEL, model)

        medaka_combo = QComboBox()
        medaka_combo.setEditable(True)
        self._populate_medaka_combo(medaka_combo, sample)
        medaka_combo.currentTextChanged.connect(
            lambda text, r=row: self._on_medaka_changed(r, text)
        )
        self.table.setCellWidget(row, COL_MEDAKA, medaka_combo)

        reads = QTableWidgetItem(sample.stats.display)
        reads.setFlags(Qt.ItemFlag.ItemIsEnabled | Qt.ItemFlag.ItemIsSelectable)
        self.table.setItem(row, COL_READS, reads)

        status = QTableWidgetItem("ready")
        status.setFlags(Qt.ItemFlag.ItemIsEnabled | Qt.ItemFlag.ItemIsSelectable)
        self.table.setItem(row, COL_STATUS, status)

    def _populate_ref_combo(self, combo: QComboBox, sample: Sample) -> None:
        combo.blockSignals(True)
        combo.clear()
        combo.addItem("— none —", None)
        for name, ref in sorted(self.references.items()):
            combo.addItem(f"{name} ({ref.length:,} bp)", name)
        idx = combo.findData(sample.reference)
        combo.setCurrentIndex(idx if idx >= 0 else 0)
        combo.setEnabled(sample.mode.needs_reference)
        combo.blockSignals(False)

    def _populate_medaka_combo(self, combo: QComboBox, sample: Sample) -> None:
        combo.blockSignals(True)
        combo.clear()
        auto = sample.medaka_model or ""
        combo.addItem(f"Auto: {auto}" if auto else "Auto (none found)", "")
        for m in self._medaka_choices:
            combo.addItem(m, m)
        if sample.medaka_model_source == "manual" and sample.medaka_model:
            idx = combo.findData(sample.medaka_model)
            combo.setCurrentIndex(idx if idx >= 0 else 0)
        else:
            combo.setCurrentIndex(0)
        combo.setToolTip(
            "Chosen automatically from the basecalling model in the read headers. "
            "Pick a specific model to override."
        )
        combo.blockSignals(False)

    def _refresh_ref_combos(self) -> None:
        for row, sample in enumerate(self.samples):
            combo = self.table.cellWidget(row, COL_REF)
            if isinstance(combo, QComboBox):
                self._populate_ref_combo(combo, sample)

    def _sync_reference_list(self) -> None:
        self.ref_list.clear()
        for name, ref in sorted(self.references.items()):
            item = QListWidgetItem(
                f"{name}  ·  {ref.length:,} bp  ·  {ref.n_seqs} seq"
            )
            item.setData(Qt.ItemDataRole.UserRole, name)
            item.setToolTip(str(ref.path))
            self.ref_list.addItem(item)

    def _on_item_changed(self, item: QTableWidgetItem) -> None:
        row = item.row()
        if row >= len(self.samples):
            return
        sample = self.samples[row]
        if item.column() == COL_ON:
            sample.enabled = item.checkState() == Qt.CheckState.Checked
        elif item.column() == COL_SIZE:
            text = item.text().replace(",", "").replace(" ", "").strip()
            if not text:
                sample.expected_size = 0
            else:
                try:
                    sample.expected_size = max(0, int(float(text)))
                except ValueError:
                    pass
            item.setText(f"{sample.expected_size:,}" if sample.expected_size else "")
        elif item.column() == COL_NAME:
            new = item.text().strip()
            if new and new != sample.name:
                sample.name = self._unique_name(new) if any(
                    s is not sample and s.name == new for s in self.samples
                ) else new
                item.setText(sample.name)

    def _on_mode_changed(self, row: int) -> None:
        combo = self.table.cellWidget(row, COL_MODE)
        if isinstance(combo, QComboBox) and row < len(self.samples):
            sample = self.samples[row]
            sample.mode = Mode(combo.currentData())
            ref_combo = self.table.cellWidget(row, COL_REF)
            if isinstance(ref_combo, QComboBox):
                self._populate_ref_combo(ref_combo, sample)

    def _on_ref_changed(self, row: int) -> None:
        combo = self.table.cellWidget(row, COL_REF)
        if isinstance(combo, QComboBox) and row < len(self.samples):
            self.samples[row].reference = combo.currentData()

    def _on_medaka_changed(self, row: int, text: str) -> None:
        if row >= len(self.samples):
            return
        sample = self.samples[row]
        if text.startswith("Auto"):
            sample.medaka_model_source = "auto"
        elif text.strip():
            sample.medaka_model = text.strip()
            sample.medaka_model_source = "manual"

    # ---------------------------------------------------- background scan
    def _scan_new_samples(self) -> None:
        pool = QThreadPool.globalInstance()
        for sample in self.samples:
            if sample.basecall.models or not sample.fastqs:
                continue
            job = ScanJob(sample.name, list(sample.fastqs), self._medaka_bin())
            job.signals.done.connect(self._on_scan_done)
            pool.start(job)

    @Slot(str, object, object)
    def _on_scan_done(self, name: str, info, resolution) -> None:
        model, source, note = resolution
        for row, sample in enumerate(self.samples):
            if sample.name != name:
                continue
            sample.basecall = info
            if sample.medaka_model_source != "manual":
                sample.medaka_model = model or ""
                sample.medaka_model_source = source
            item = self.table.item(row, COL_MODEL)
            if item:
                item.setText(info.display or "not in headers")
                if info.is_mixed:
                    item.setForeground(QColor("#8a6100"))
                    item.setToolTip(
                        "More than one basecalling model in this sample:\n"
                        + "\n".join(f"{k}: {v} reads" for k, v in info.models.items())
                    )
                elif not info.models:
                    item.setForeground(QColor("#b3261e"))
                    item.setToolTip(
                        "No basecalling model tag found. Medaka model selection will "
                        "fall back to samtools/bcftools unless you set one manually."
                    )
            combo = self.table.cellWidget(row, COL_MEDAKA)
            if isinstance(combo, QComboBox) and sample.medaka_model_source != "manual":
                self._populate_medaka_combo(combo, sample)
                if note:
                    combo.setToolTip(note)
            break

    # ------------------------------------------------------- sample sheet
    def load_sheet(self) -> None:
        path, _ = QFileDialog.getOpenFileName(self, "Load sample sheet", str(Path.home()),
                                              "CSV (*.csv);;All files (*)")
        if not path:
            return
        try:
            samples, refs, warnings = read_sheet(Path(path))
        except (OSError, ValueError) as exc:
            QMessageBox.critical(self, APP_NAME, f"Could not read sample sheet:\n{exc}")
            return
        for ref in refs:
            self._add_reference(ref.path)
        for s in samples:
            s.name = self._unique_name(s.name)
            self.samples.append(s)
        self._rebuild_table()
        self._sync_reference_list()
        self._scan_new_samples()
        if warnings:
            QMessageBox.warning(self, APP_NAME, "Sample sheet warnings:\n\n" +
                                "\n".join(warnings[:20]))

    def save_sheet(self) -> None:
        path, _ = QFileDialog.getSaveFileName(self, "Save sample sheet",
                                              str(Path.home() / "samples.csv"),
                                              "CSV (*.csv)")
        if path:
            write_sheet(Path(path), self.samples, self.references)
            self.status.showMessage(f"Saved sample sheet to {path}")

    def pick_out_dir(self) -> None:
        folder = QFileDialog.getExistingDirectory(self, "Output folder", self.out_edit.text())
        if folder:
            self.out_edit.setText(folder)

    def open_out_dir(self) -> None:
        target = self.run_dir or Path(self.out_edit.text())
        if target.exists():
            subprocess.run(["open", str(target)], check=False)
        else:
            QMessageBox.information(self, APP_NAME, "That folder does not exist yet.")

    def open_report(self) -> None:
        if self.run_dir and (self.run_dir / "report.html").exists():
            subprocess.run(["open", str(self.run_dir / "report.html")], check=False)

    def _open_result_row(self, row: int, _col: int) -> None:
        item = self.results_table.item(row, 5)
        if item and item.text():
            subprocess.run(["open", item.text()], check=False)

    # ---------------------------------------------------------------- run
    def _collect_config(self) -> RunConfig:
        return RunConfig(
            out_dir=Path(self.out_edit.text()).expanduser(),
            threads_per_job=self.threads_spin.value(),
            parallel_jobs=self.parallel_spin.value(),
            min_read_len=self.minlen_spin.value(),
            max_read_len=self.maxlen_spin.value(),
            min_read_q=self.minq_spin.value(),
            target_coverage=self.cov_spin.value(),
            keep_intermediates=self.keep_check.isChecked(),
            use_medaka=self.medaka_check.isChecked(),
            min_consensus_depth=self.depth_spin.value(),
        )

    def start_run(self) -> None:
        if self.job is not None:
            QMessageBox.information(self, APP_NAME, "A run is already in progress.")
            return
        missing = tools.missing_required(refresh=True)
        if missing:
            QMessageBox.critical(
                self, APP_NAME,
                "These required tools are missing:\n\n  " + "\n  ".join(missing) +
                "\n\nRun scripts/fetch_tools.sh to build the bundled toolchain."
            )
            return
        active = [s for s in self.samples if s.enabled]
        if not active:
            QMessageBox.information(self, APP_NAME, "Add at least one sample first.")
            return
        problems: list[str] = []
        for s in active:
            problems.extend(s.validate(set(self.references)))
        if problems:
            QMessageBox.critical(self, APP_NAME,
                                 "Fix these before running:\n\n" + "\n".join(problems[:20]))
            return

        total_threads = self.parallel_spin.value() * self.threads_spin.value()
        cpus = os.cpu_count() or 4
        if total_threads > cpus * 1.5:
            reply = QMessageBox.question(
                self, APP_NAME,
                f"{self.parallel_spin.value()} samples x {self.threads_spin.value()} "
                f"threads = {total_threads} threads on {cpus} cores. This will "
                f"oversubscribe the machine and may run slower. Continue?",
            )
            if reply != QMessageBox.StandardButton.Yes:
                return

        self.cfg = self._collect_config()
        self.run_dir = self.cfg.out_dir / time.strftime("run_%Y%m%d_%H%M%S")
        self.results = []
        self.results_table.setRowCount(0)
        self.log_view.clear()
        self.progress.setRange(0, len(active))
        self.progress.setValue(0)
        self.run_btn.setEnabled(False)
        self.stop_btn.setEnabled(True)
        self.report_btn.setEnabled(False)
        for row in range(self.table.rowCount()):
            item = self.table.item(row, COL_STATUS)
            if item:
                item.setText("queued" if self.samples[row].enabled else "skipped")

        self.job = RunJob(active, dict(self.references), self.cfg, self.run_dir)
        self.job.signals.log.connect(self.append_log)
        self.job.signals.sample_started.connect(self._on_sample_started)
        self.job.signals.sample_finished.connect(self._on_sample_finished)
        self.job.signals.progress.connect(lambda done, _t: self.progress.setValue(done))
        self.job.signals.finished.connect(self._on_run_finished)
        QThreadPool.globalInstance().start(self.job)
        self.tabs.setCurrentIndex(0)
        self.status.showMessage(f"Running - output in {self.run_dir}")

    def stop_run(self) -> None:
        if self.job:
            self.job.stop()
            self.append_log("Stop requested - finishing the current step, then halting.")
            self.stop_btn.setEnabled(False)

    @Slot(str)
    def append_log(self, message: str) -> None:
        self.log_view.appendPlainText(message)

    @Slot(str)
    def _on_sample_started(self, name: str) -> None:
        self._set_status(name, "running…")

    @Slot(object)
    def _on_sample_finished(self, result) -> None:
        if not result.ok:
            status = "failed"
        elif result.warnings:
            status = f"done ({len(result.warnings)} warnings)"
        else:
            status = "done"
        self._set_status(result.sample, status)
        for row, sample in enumerate(self.samples):
            if sample.name == result.sample:
                item = self.table.item(row, COL_READS)
                if item:
                    item.setText(sample.stats.display)
                break
        self.results.append(result)
        self._append_result_row(result)

    def _set_status(self, name: str, text: str) -> None:
        for row, sample in enumerate(self.samples):
            if sample.name == name:
                item = self.table.item(row, COL_STATUS)
                if item:
                    item.setText(text)
                    if text == "failed":
                        item.setForeground(QColor("#b3261e"))
                    elif text.startswith("done ("):
                        item.setForeground(QColor("#8a6100"))
                    elif text == "done":
                        item.setForeground(QColor("#1a7f52"))
                break

    def _append_result_row(self, result) -> None:
        row = self.results_table.rowCount()
        self.results_table.insertRow(row)
        metrics = result.metrics
        interesting = [
            ("reads", "reads"), ("est_coverage", "x est. cov"),
            ("mean_depth", "x depth"), ("breadth_min", "% breadth"),
            ("contigs", "contigs"), ("assembly_bp", "bp assembly"),
            ("consensus_n_pct", "% N"), ("variants", "variants"),
        ]
        summary = ", ".join(
            f"{metrics[k]:,}{'' if u.startswith('%') or u.startswith('x') else ' '}{u}"
            if isinstance(metrics[k], int) else f"{metrics[k]}{u}"
            for k, u in interesting if k in metrics
        )
        if not result.ok:
            status_text = f"FAILED - {result.message}"
        elif result.warnings:
            status_text = f"OK - {len(result.warnings)} warnings"
        else:
            status_text = "OK"
        values = [
            result.sample,
            status_text,
            str(metrics.get("basecall_model", "")),
            str(metrics.get("medaka_model", "")),
            summary,
            str(result.out_dir or ""),
        ]
        for col, value in enumerate(values):
            item = QTableWidgetItem(value)
            item.setFlags(Qt.ItemFlag.ItemIsEnabled | Qt.ItemFlag.ItemIsSelectable)
            if col == 1:
                if not result.ok:
                    item.setForeground(QColor("#b3261e"))
                elif result.warnings:
                    item.setForeground(QColor("#8a6100"))
                    item.setToolTip("\n\n".join(result.warnings))
                else:
                    item.setForeground(QColor("#1a7f52"))
            self.results_table.setItem(row, col, item)

    @Slot(object)
    def _on_run_finished(self, results) -> None:
        self.job = None
        self.run_btn.setEnabled(True)
        self.stop_btn.setEnabled(False)
        if self.run_dir:
            try:
                path = report.write_report(self.run_dir, list(results), self.cfg)
                self.append_log(f"Report: {path}")
                self.report_btn.setEnabled(True)
            except OSError as exc:
                self.append_log(f"Could not write report: {exc}")
        n_ok = sum(1 for r in results if r.ok)
        self.status.showMessage(f"Finished: {n_ok}/{len(results)} samples succeeded")
        self.tabs.setCurrentIndex(1)

    # -------------------------------------------------------- drag & drop
    def dragEnterEvent(self, event) -> None:
        if event.mimeData().hasUrls():
            event.acceptProposedAction()

    def dropEvent(self, event) -> None:
        paths: list[Path] = []
        for url in event.mimeData().urls():
            p = Path(url.toLocalFile())
            if p.is_dir():
                paths.extend(sorted(f for f in p.rglob("*")
                                    if f.is_file() and (is_fastq(f) or is_fasta(f))))
            elif p.is_file():
                paths.append(p)
        self._add_paths(paths)
        event.acceptProposedAction()

    def show_about(self) -> None:
        QMessageBox.about(
            self, f"About {APP_NAME}",
            f"<b>{APP_NAME} {__version__}</b><br><br>"
            "Nanopore assembly and reference analysis for small genomes "
            "(amplicons, plasmids, viral genomes) on Apple Silicon.<br><br>"
            "The basecalling model is read from each FASTQ's headers and used to "
            "pick the matching Medaka model and Flye read profile.<br><br>"
            "Uses minimap2, samtools, bcftools, Flye, Medaka and RagTag."
        )

    def closeEvent(self, event) -> None:
        if self.job is not None:
            reply = QMessageBox.question(self, APP_NAME,
                                         "A run is in progress. Stop it and quit?")
            if reply != QMessageBox.StandardButton.Yes:
                event.ignore()
                return
            self.job.stop()
        event.accept()
