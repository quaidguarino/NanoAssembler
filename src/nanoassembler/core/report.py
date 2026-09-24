"""Self-contained HTML summary for a finished run."""
from __future__ import annotations

import html
import json
import time
from pathlib import Path

from .models import RunConfig
from .pipeline import SampleResult

_CSS = """
:root { color-scheme: light dark;
  --bg:#fbfbfa; --card:#fff; --ink:#1a1a19; --muted:#6b6b68; --line:#e4e4e1;
  --ok:#1a7f52; --bad:#b3261e; --warn:#8a6100; --accent:#2b5f9e; }
@media (prefers-color-scheme: dark) { :root {
  --bg:#16171a; --card:#1e2024; --ink:#e9e9e7; --muted:#9a9a96; --line:#31343a;
  --ok:#4cc38a; --bad:#f2726a; --warn:#e0a72e; --accent:#7aa9e0; } }
* { box-sizing:border-box }
body { margin:0; padding:32px; background:var(--bg); color:var(--ink);
  font:14px/1.55 -apple-system,BlinkMacSystemFont,"SF Pro Text",Helvetica,sans-serif; }
h1 { font-size:22px; margin:0 0 4px }
h2 { font-size:15px; margin:28px 0 10px; letter-spacing:.02em; text-transform:uppercase;
     color:var(--muted); font-weight:600 }
.sub { color:var(--muted); margin-bottom:24px }
table { border-collapse:collapse; width:100%; background:var(--card);
  border:1px solid var(--line); border-radius:10px; overflow:hidden }
th,td { padding:9px 12px; text-align:left; border-bottom:1px solid var(--line);
  font-variant-numeric:tabular-nums; white-space:nowrap }
th { background:color-mix(in srgb, var(--card) 80%, var(--line)); font-weight:600;
  font-size:12px; text-transform:uppercase; letter-spacing:.03em; color:var(--muted) }
tr:last-child td { border-bottom:none }
.wrap { overflow-x:auto }
.ok { color:var(--ok); font-weight:600 } .bad { color:var(--bad); font-weight:600 }
.warn { color:var(--warn) }
.card { background:var(--card); border:1px solid var(--line); border-radius:10px;
  padding:16px 18px; margin-bottom:14px }
.card h3 { margin:0 0 4px; font-size:16px }
.meta { color:var(--muted); font-size:13px; margin-bottom:10px }
.kv { display:grid; grid-template-columns:repeat(auto-fill,minmax(190px,1fr)); gap:8px 18px }
.kv div { font-size:13px } .kv b { color:var(--muted); font-weight:500 }
code { font:12px ui-monospace,SFMono-Regular,Menlo,monospace; color:var(--accent);
  word-break:break-all }
.files a { color:var(--accent); text-decoration:none; font:12px ui-monospace,Menlo,monospace }
.files a:hover { text-decoration:underline }
"""


def _fmt(value: object) -> str:
    if isinstance(value, float):
        return f"{value:,.2f}".rstrip("0").rstrip(".")
    if isinstance(value, int):
        return f"{value:,}"
    return html.escape(str(value))


_SUMMARY_COLS = [
    ("sample", "Sample"),
    ("mode", "Mode"),
    ("status", "Status"),
    ("basecall_model", "Basecall model"),
    ("medaka_model", "Medaka model"),
    ("reads", "Reads"),
    ("yield_mb", "Yield (Mb)"),
    ("read_n50", "Read N50"),
    ("mean_q", "Mean Q"),
    ("est_coverage", "Est. cov"),
    ("mean_depth", "Mean depth"),
    ("breadth_min", "Breadth %"),
    ("contigs", "Contigs"),
    ("assembly_bp", "Assembly bp"),
    ("consensus_n_pct", "% N"),
    ("variants", "Variants"),
]


def write_report(run_dir: Path, results: list[SampleResult], cfg: RunConfig) -> Path:
    run_dir = Path(run_dir)
    rows = []
    used = set()
    for r in results:
        m = dict(r.metrics)
        m["sample"] = r.sample
        m["mode"] = r.mode
        m["status"] = "OK" if r.ok else "FAILED"
        m.pop("mixed_models", None)  # surfaced as a banner, not a column
        rows.append(m)
        used.update(k for k, v in m.items() if v not in (None, "", 0))

    cols = [(k, label) for k, label in _SUMMARY_COLS if k in used or k in
            ("sample", "mode", "status", "basecall_model")]

    head = "".join(f"<th>{html.escape(label)}</th>" for _k, label in cols)
    body = []
    for m in rows:
        cells = []
        for k, _label in cols:
            v = m.get(k, "")
            if k == "status":
                cls = "ok" if v == "OK" else "bad"
                cells.append(f'<td class="{cls}">{v}</td>')
            else:
                cells.append(f"<td>{_fmt(v) if v not in (None, '') else '&ndash;'}</td>")
        body.append("<tr>" + "".join(cells) + "</tr>")

    cards = []
    for r in results:
        steps = "".join(
            f'<div><b>{html.escape(s.name)}</b> '
            f'<span class="{"ok" if s.ok else "bad"}">{"OK" if s.ok else "failed"}</span> '
            f"{html.escape(s.detail)}</div>"
            for s in r.steps
        )
        files = "".join(
            f'<div class="files">{html.escape(k)}: '
            f'<a href="file://{html.escape(v)}">{html.escape(Path(v).name)}</a></div>'
            for k, v in sorted(r.outputs.items())
        )
        note = "" if r.ok else f'<div class="bad">{html.escape(r.message)}</div>'
        if r.warnings:
            items = "".join(f"<li>{html.escape(w)}</li>" for w in r.warnings)
            note += (f'<div class="warns"><b>Structural warnings '
                     f'({len(r.warnings)})</b><ul>{items}</ul></div>')
        cards.append(
            f'<div class="card"><h3>{html.escape(r.sample)}</h3>'
            f'<div class="meta">{html.escape(r.mode)} &middot; '
            f"{r.duration_s/60:.1f} min &middot; "
            f'<code>{html.escape(str(r.out_dir or ""))}</code></div>'
            f"{note}<div class='kv'>{steps}</div>{files}</div>"
        )

    mixed = [r.sample for r in results if r.metrics.get("mixed_models")]
    warn = ""
    if mixed:
        warn = (f'<p class="warn">Mixed basecalling models detected in: '
                f'{html.escape(", ".join(mixed))}</p>')

    flagged = [r.sample for r in results if r.warnings]
    if flagged:
        warn += (f'<p class="warn">Structural warnings raised for: '
                 f'{html.escape(", ".join(flagged))} &mdash; see the per-sample '
                 f'sections below.</p>')

    n_ok = sum(1 for r in results if r.ok)
    doc = f"""<!doctype html><html><head><meta charset="utf-8">
<title>NanoAssembler run - {html.escape(run_dir.name)}</title>
<style>{_CSS}</style></head><body>
<h1>NanoAssembler run report</h1>
<div class="sub">{html.escape(run_dir.name)} &middot;
{time.strftime('%Y-%m-%d %H:%M')} &middot;
{n_ok}/{len(results)} samples completed</div>
{warn}
<h2>Summary</h2>
<div class="wrap"><table><thead><tr>{head}</tr></thead><tbody>{''.join(body)}</tbody></table></div>
<h2>Samples</h2>
{''.join(cards)}
<h2>Run settings</h2>
<div class="card"><pre style="margin:0;white-space:pre-wrap"><code>{html.escape(cfg.to_json())}</code></pre></div>
</body></html>"""
    out = run_dir / "report.html"
    out.write_text(doc, encoding="utf-8")

    (run_dir / "summary.json").write_text(json.dumps(rows, indent=2, default=str))
    with open(run_dir / "summary.tsv", "w") as fh:
        fh.write("\t".join(label for _k, label in cols) + "\n")
        for m in rows:
            fh.write("\t".join(str(m.get(k, "")) for k, _l in cols) + "\n")
    return out
