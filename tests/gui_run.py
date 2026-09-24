"""Drive a real two-sample run through the GUI's threading machinery."""
import sys, time
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from PySide6.QtWidgets import QApplication, QComboBox
from nanoassembler.ui.main_window import MainWindow, COL_MODE, COL_REF, COL_STATUS
from nanoassembler.core.models import Mode

data, out = Path(sys.argv[1]), Path(sys.argv[2])
app = QApplication([])
w = MainWindow()
w.out_edit.setText(str(out))
w.parallel_spin.setValue(2)
w.threads_spin.setValue(4)
w.minlen_spin.setValue(500)
w.minq_spin.setValue(7.0)

w._add_paths([data / "sample01.fastq.gz", data / "sample02.fastq.gz",
              data / "sample01_reference.fasta"])
deadline = time.time() + 30
while time.time() < deadline and not all(s.basecall.models for s in w.samples):
    app.processEvents(); time.sleep(0.1)

# sample01 -> reference consensus + variants; sample02 -> de novo
w.table.cellWidget(0, COL_MODE).setCurrentIndex(list(Mode).index(Mode.REF_CONSENSUS_VARIANTS))
w.table.cellWidget(0, COL_REF).setCurrentIndex(1)
w.table.cellWidget(1, COL_MODE).setCurrentIndex(list(Mode).index(Mode.DENOVO))
print(f"sample01 -> {w.samples[0].mode.label} / {w.samples[0].reference}")
print(f"sample02 -> {w.samples[1].mode.label}")

w.start_run()
deadline = time.time() + 1500
while time.time() < deadline and w.job is not None:
    app.processEvents(); time.sleep(0.2)

print(f"\nrun dir: {w.run_dir}")
for row in range(w.table.rowCount()):
    print(f"  table status[{row}] = {w.table.item(row, COL_STATUS).text()}")
for r in w.results:
    print(f"  {r.sample}: {'OK' if r.ok else 'FAILED ' + r.message} "
          f"({r.duration_s:.0f}s) outputs={sorted(r.outputs)}")
assert w.results and all(r.ok for r in w.results), "a sample failed"
assert (w.run_dir / "report.html").exists(), "no report written"
assert w.report_btn.isEnabled()
print(f"  results tab rows: {w.results_table.rowCount()}")
print("GUI RUN TEST PASSED")

# capture the window for a screenshot without needing a display
w.resize(1440, 880)
w.grab().save(str(out / "nanoassembler_ui.png"))
print("screenshot:", out / "nanoassembler_ui.png")
