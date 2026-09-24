# NanoAssembler

A macOS (Apple Silicon) GUI for Oxford Nanopore data: de novo assembly, reference
consensus and variant calling for small genomes — amplicons, plasmids and viral
genomes. Many samples and many references in one run, each sample paired with
whichever reference you choose.

The basecalling model is read out of each FASTQ's headers and used automatically
to pick the matching Medaka model, the Medaka variant model, Flye's read profile
and the bcftools error profile.

## Setup

Two commands, both one-time:

```bash
./scripts/fetch_tools.sh
```

Builds the whole toolchain into `vendor/`. minimap2, samtools and bcftools are
compiled from upstream release tarballs as native arm64 binaries; a relocatable
CPython goes into `vendor/python` and receives PySide6, Flye, Medaka, pyabpoa
and RagTag. No conda, no Homebrew, and nothing installed system-wide. Takes a
few minutes and about 3 GB.

A relocatable interpreter is used deliberately: a normal `venv` records an
absolute path to the Python it was created from, so it cannot be moved into an
.app bundle or copied to another Mac.

```bash
./scripts/build_app.sh
```

Assembles `dist/NanoAssembler.app` by copying `vendor/` and the source into a
bundle and ad-hoc signing it. The result runs on any Apple Silicon Mac with no
Python, conda or Xcode installed — copy it to /Applications or hand it to a
colleague. To run from the source tree instead:

```bash
./scripts/run.sh
```

## Using it

1. **Add samples.** Drag FASTQ files onto the table, or use *Add FASTQ files* /
   *Add folder*. One FASTQ becomes one sample, named from the filename. To pool
   several files into one sample, select the rows and press *Merge selected*.
2. **Add references.** Drag FASTA files in, or use *Add FASTA* in the References
   panel. Load as many as you like.
3. **Set the mode per sample** in the Mode column, and pick that sample's
   reference in the Reference column. Different samples can use different modes
   and different references in the same run.
4. **Press Run.** Samples are processed in parallel; the Log tab streams every
   command, and the Results tab fills in as samples finish.

| Mode | What it does | Output |
|---|---|---|
| De novo assembly | Flye, then Medaka polishing | `assembly.fasta` |
| Reference: consensus | minimap2 → Medaka consensus (samtools consensus if Medaka is off) | `consensus.fasta`, `aligned.bam`, `depth.tsv` |
| Reference: variants | minimap2 → Medaka variant network (bcftools if unavailable) | `variants.vcf`, `aligned.bam` |
| Reference: consensus + variants | both of the above from one alignment | both |
| De novo + scaffold to reference | Flye, polish, then RagTag ordering | `assembly.fasta`, `assembly.scaffolded.fasta` |

Every run writes a timestamped folder containing per-sample subfolders, a
`report.html` summary, `summary.tsv`, `summary.json`, and a `run_info.json` per
sample recording the basecalling model, the models chosen from it, tool versions
and the exact command line of every step.

## How the basecalling model is detected

Dorado and Guppy stamp the model into each read header:

```
@0a1b… runid=… ch=105 basecall_model_version_id=dna_r10.4.1_e8.2_400bps_sup@v5.0.0
```

The app samples ~4000 headers per sample and reads that tag; if it is absent it
falls back to the `RG:Z:` read-group string that `samtools fastq -T RG` produces.
The id is then rewritten into Medaka's convention
(`r1041_e82_400bps_sup_v5.0.0`), checked against the models your Medaka actually
ships, and used to select:

- the Medaka consensus model,
- the Medaka **variant** model (a separate network, `…_sup_variant_v5.0.0`),
- Flye's read profile (`--nano-hq` for sup and recent hac, `--nano-raw` otherwise),
- the bcftools platform profile (`-X ont-sup` for sup, `-X ont` otherwise).

If the exact model is not installed, the closest one of the same chemistry and
accuracy is used and the substitution is logged. Models needing signal-level
data absent from FASTQ (`…_rl_lstm384_dwells`) are never chosen automatically.
Any row's model can be overridden by hand in the Medaka model column, and a
sample containing reads from more than one model is flagged in amber.

## Repeated features (multiple promoters, repeated cassettes)

Constructs carrying several identical cassettes are the classic way for de novo
assembly to go quietly wrong. A repeat-graph assembler can only place a repeat
correctly if reads span the whole repeat *plus* unique sequence on both sides;
otherwise it collapses the copies and returns a short assembly, with no error.
That is why a restriction digest can show the right size while the assembly does
not.

Measured on simulated constructs with three identical cassettes:

| Cassette | Read N50 | Assembled | Copies recovered | Result |
|---|---|---|---|---|
| 2.5 kb | 1.8 kb | 12,920 / 13,900 bp | 2 of 3 | truncated |
| 2.5 kb | 3.6 kb | 13,900 bp | 3 of 3 | correct |
| 4 kb | 1.8 kb | 9,642 / 17,200 bp | 1 of 3 | truncated 44% |
| 4 kb | 3.5 kb | 17,200 bp | 3 of 3 | correct |
| 4 kb | 14.2 kb | 34,413 bp | 6 | circle assembled twice |

The pattern: reads must be comfortably longer than one cassette, and once reads
approach the length of the whole plasmid, Flye starts emitting doubled circles
instead.

**The app checks for both failures rather than letting them pass.** After
assembly it maps the reads back onto their own assembly and compares:

- **Assembled size against the expected size** you type in the Expected size
  column. Anything under 95% is flagged.
- **Read depth across the assembly.** Three identical cassettes assembled into
  one carry roughly three times the median depth, so a collapsed repeat shows up
  as a depth spike with an estimated number of missing copies. This is the
  decisive evidence, because it does not depend on the assembler's own opinion.
- **Self-alignment of each contig**, catching a circular molecule assembled
  twice; a corrected `assembly.deduplicated.fasta` is written alongside.
- **Flye's own multiplicity calls** from `assembly_info.txt`.
- **How many single reads span the full expected length.** When a truncation is
  suspected this is what settles whether it is real: a read that covers the
  whole construct cannot have collapsed anything.

Warnings appear in the log, amber in the Status column, and in a highlighted
block per sample in `report.html`. Coverage taper at contig ends is excluded, so
a clean sample raises nothing.

In reference mode the same depth analysis runs against your expected map:
regions covered far above the median mean the sample carries extra copies, and
internal dropouts mean a deletion relative to the map.

**Practical advice for your constructs**

1. Prefer **Reference: consensus** when you have the expected plasmid map.
   Mapping never collapses repeats, so it returns the full-length sequence and
   the depth profile tells you where reality departs from the design.
2. For de novo, always fill in **Expected size** — without it the size check
   cannot run.
3. Aim for a prep whose read N50 is at least ~1.5x your largest repeated
   cassette. Skipping harsh shearing matters more than adding coverage.
4. If reads span the whole construct, that is the strongest evidence available;
   the app counts them for you.

## Settings that matter

- **Target coverage** — amplicons are often sequenced to thousands of ×, which
  slows assembly without improving it. The default subsamples the longest reads
  down to 200× of the reference length (or your stated expected size).
- **Min consensus depth** — positions below it are masked `N` rather than
  reported as reference. Raise it if you care about not over-calling.
- **Samples in parallel × threads per sample** should stay near your core count;
  the app warns if you oversubscribe.

## Known limits

- Sized for targets under a few hundred kb. A bacterial isolate will work; a
  eukaryotic genome will not fit in 16 GB.
- Racon (the fallback polisher, used only when Medaka is unavailable) needs
  `cmake` at build time and is skipped if it is missing. Medaka installs
  natively on Apple Silicon, so this rarely matters.
- The bundle is large (roughly 2–3 GB) because Medaka pulls in PyTorch.
- Medaka's own `medaka_consensus`/`medaka_variant` shell wrappers do not quote
  their arguments and break on paths containing spaces. This app calls Medaka's
  Python subcommands directly and never uses those wrappers.
- R9.4.1 and R10.3 basecalling models are mapped to the MinION/GridION Medaka
  models. On PromethION data, override to the `r941_prom_*` model.

## Testing

```bash
PY=vendor/python/bin/python3
$PY tests/make_test_data.py testdata        # simulated ONT reads with model tags
$PY tests/e2e.py testdata out               # every mode, end to end
QT_QPA_PLATFORM=offscreen PYTHONPATH=src $PY tests/gui_smoke.py testdata
QT_QPA_PLATFORM=offscreen PYTHONPATH=src $PY tests/gui_run.py testdata out
$PY tests/spaces.py testdata                  # sample names / folders with spaces
```

`tests/e2e.py` checks that each mode completes; on the simulated data the
reference consensus comes back 100% identical to the known truth sequence over
the region covered at the minimum depth.
