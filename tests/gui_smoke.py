"""Headless GUI check: add files, confirm the table fills in from the headers."""
import sys, time
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from PySide6.QtCore import QTimer
from PySide6.QtWidgets import QApplication, QComboBox
from nanoassembler.ui.main_window import (MainWindow, COL_NAME, COL_FILES, COL_MODE,
                                          COL_REF, COL_MODEL, COL_MEDAKA, COL_STATUS)
from nanoassembler.core.models import Mode
from nanoassembler.core.samplesheet import read as read_sheet, write as write_sheet

data = Path(sys.argv[1])
app = QApplication([])
w = MainWindow()

# Simulate a drag-and-drop of two FASTQs and one reference
fq = data / "sample01.fastq.gz"
fq2 = data / "sample02.fastq.gz"
import shutil; shutil.copy(fq, fq2)
w._add_paths([fq, fq2, data / "sample01_reference.fasta"])
assert len(w.samples) == 2, w.samples
assert len(w.references) == 1, w.references
print(f"added {len(w.samples)} samples, {len(w.references)} reference")

# Give the background scan threads a moment
deadline = time.time() + 30
while time.time() < deadline and not all(s.basecall.models for s in w.samples):
    app.processEvents(); time.sleep(0.1)

for row, s in enumerate(w.samples):
    print(f"  row {row}: name={w.table.item(row,COL_NAME).text()!r} "
          f"files={w.table.item(row,COL_FILES).text()} "
          f"model={w.table.item(row,COL_MODEL).text()!r} "
          f"medaka={w.table.cellWidget(row,COL_MEDAKA).currentText()!r}")
    assert s.basecall.primary == "dna_r10.4.1_e8.2_400bps_sup@v5.0.0", s.basecall.models
    assert s.medaka_model == "r1041_e82_400bps_sup_v5.0.0", s.medaka_model

# Switch mode to a reference mode and assign the reference
combo = w.table.cellWidget(0, COL_MODE)
combo.setCurrentIndex(list(Mode).index(Mode.REF_CONSENSUS_VARIANTS))
assert w.samples[0].mode is Mode.REF_CONSENSUS_VARIANTS
ref_combo = w.table.cellWidget(0, COL_REF)
assert ref_combo.isEnabled(), "reference combo should enable for reference modes"
ref_combo.setCurrentIndex(1)
assert w.samples[0].reference == "sample01_reference", w.samples[0].reference
print(f"  mode/reference wiring OK -> {w.samples[0].mode.label} / {w.samples[0].reference}")

# De novo mode must not require a reference
assert w.samples[1].validate(set(w.references)) == []
# A reference mode without a reference must be caught
w.samples[1].mode = Mode.REF_CONSENSUS
problems = w.samples[1].validate(set(w.references))
assert problems and "needs a reference" in problems[0], problems
print(f"  validation OK -> {problems[0]}")
w.samples[1].mode = Mode.DENOVO

# Sample sheet round trip
sheet = data / "sheet.csv"
write_sheet(sheet, w.samples, w.references)
loaded, refs, warn = read_sheet(sheet)
assert len(loaded) == 2 and not warn, (loaded, warn)
assert loaded[0].mode is Mode.REF_CONSENSUS_VARIANTS and loaded[0].reference
print(f"  sample sheet round-trip OK ({len(loaded)} samples, {len(refs)} refs)")

# Manual medaka override
mc = w.table.cellWidget(0, COL_MEDAKA)
mc.setCurrentText("r1041_e82_400bps_hac_v4.2.0")
assert w.samples[0].medaka_model == "r1041_e82_400bps_hac_v4.2.0"
assert w.samples[0].medaka_model_source == "manual"
print("  manual medaka override OK")

print("GUI SMOKE TEST PASSED")
