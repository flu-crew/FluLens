![FluLens — influenza variant visualizer](docs/img/banner.png)

FluLens is a visualizer for inspecting and filtering influenza variant calls. **samples × codons**, one cell per call and
coloured by allele frequency or viewable as a consensus.

It reads the output of [Flumina](https://github.com/flu-crew/Flumina) directly and
answers two questions faster than a spreadsheet does: *is this variant real?* and
*is this sample worth keeping?*

**[▶ Try it in your browser](https://flu-crew.github.io/FluLens/?run=example_run)** —
no installation needed and loaded with example data.

![FluLens showing a variant's assessment panel](docs/img/variant-panel.png)

---

## Running it

Three ways. They are the same application; pick whichever you like best!

| | Use when |
|---|---|
| [1. In a browser](#1-in-a-browser) | Simplest. Nothing to install, works on any OS |
| [2. As a single file](#2-as-a-single-file) | You want it offline, or on a machine with no internet |
| [3. As a desktop app](#3-as-a-desktop-app) | macOS, and you want it to remember your last run [coming soon]|

### 1. In a browser

Open **<https://flu-crew.github.io/FluLens/>**, click **Open run folder…**, and
choose a Flumina output directory.

Nothing is uploaded. The page reads the folder locally through the browser's own
file picker — see [Your data stays on your machine](#your-data-stays-on-your-machine).

### 2. As a single file

FluLens is one self-contained HTML file with no dependencies and no build step.
Download `flulens.html` from the
[latest release](https://github.com/flu-crew/FluLens/releases/latest) and open it.

```bash
open flulens.html        # macOS
xdg-open flulens.html    # Linux
```

### 3. As a desktop app [coming soon]

Download `FluLens_<version>_universal.dmg` from the
[latest release](https://github.com/flu-crew/FluLens/releases/latest), open it and
drag FluLens to Applications. The build is universal — Apple Silicon and Intel.

The desktop version has real filesystem access, so it reopens your last run
automatically instead of asking every time.

> **If macOS says the app cannot be verified**, the release was not signed. Either
> grab it from a signed release, or allow it once under **System Settings →
> Privacy & Security → Open Anyway**.

Windows and Linux builds are attached to releases too, but the browser version is
the better path on both.

---

## Try it with the example dataset

You do not need a pipeline run to see what FluLens does. `example_run/` is a small
**synthetic** Flumina output: twelve samples, 1,378 calls, all twelve gene
products.

**[▶ Open it live, nothing to download](https://flu-crew.github.io/FluLens/?run=example_run)**

To load it locally instead, which is also the way to check that the desktop app
and the folder picker work before you point them at real data:

| Get it | How |
|---|---|
| [Browse it on GitHub](https://github.com/flu-crew/FluLens/tree/main/example_run) | see the file layout a run is expected to have |
| [`example_run.zip`](https://github.com/flu-crew/FluLens/releases/latest) | attached to every release |
| [Download the whole repository](https://github.com/flu-crew/FluLens/archive/refs/heads/main.zip) | `example_run/` is inside it |

```bash
git clone https://github.com/flu-crew/FluLens.git
# then in FluLens: Open run folder… -> FluLens/example_run
```

It is deliberately built so the controls have something to act on: a library that
fails QC, two segments that never assembled, GATK4 genotype calls with no LoFreq
counterpart, and skewed strand balance on a fraction of the calls.

All twelve samples also ship with **reads**, so the pile-up opens on the example
— click any sample in the `read pile-up` view. The BAMs are capped at ~200×
depth to keep them small (a real `A1` BAM would be 86 MB); allele fractions stay
faithful, only the read count is thinned. `example_run/README.md` has the
reasoning.

---

## What to load into FluLens

A Flumina output directory: the one containing `variant_analysis/`. Everything
except the first file is optional, and FluLens says in the sidebar which of them
it found:

| Path | What it adds |
|---|---|
| `variant_analysis/all_sample_amino_acids.txt` | **required** — the grid itself |
| `reference.fa` | the translated reference row |
| `reference_gtf/*.gtf` | true protein lengths and CDS intervals, so all twelve products appear |
| `variant_analysis/curated_amino_acids.txt` | ▲ ticks marking curated sites |
| `variant_analysis/flumut/markers.tsv` | FluMut marker screening |
| `variant_analysis/flumut_lowfreq/` | markers present *below* consensus |
| `vcf_files/<sample>/lofreq-called-variants.vcf` | strand balance and read counts, read on demand |
| `IRMA_results/<sample>/tables/READ_COUNTS.txt` | the per-segment coverage strip |
| `wfabc*/FIT_results.csv` | selection coefficients and drift tests |

A **metadata CSV** can be loaded separately when the run was configured without
one. The sidebar reports how many samples the join actually matched, which is the
number worth reading: a join that silently matched half your samples looks
identical to one that matched all of them.

If you have no run of your own to hand, the
[example dataset](#try-it-with-the-example-dataset) above populates every one of
these.

---

## What it shows

![The full matrix at genome scale](docs/img/overview.png)

**The matrix.** Every call in the run, positioned by product and codon. Zoom with
the wheel, drag to pan, click a sample name to highlight its row, drag a name to
reorder. Click any header to sort; shift-click to add a second key.

**Variant detail.** Click a cell and FluLens loads that one sample's VCF and shows
the amino-acid change, the raw numbers, allele frequency on a linear or log axis,
strand balance against the *reference* allele's, and a verdict.

**Variant assessment.** Four verdicts — *Looks real*, *Treat with caution*,
*Likely artefact*, *Cannot assess* — each listing its reasons, weighing strand
balance, depth, allele frequency and the number of reads actually supporting the
call. Supporting reads are not depth: a 0.5% call on 15,000× has enormous depth
and may rest on a handful of alt reads.

**Consensus view.** Every sample's own residue at every codon, drawn as
differences from the reference rather than a wash of colour.

**Coverage strip.** Reads recovered per segment, per sample, by decade. The
variant table structurally cannot tell you a segment recovered 5,455 reads rather
than 509,426; this can.

**QC column.** A per-sample verdict from the raw variant table — invariant across
whatever filters you have set, because QC is a fact about the library and not
about the view. Click the mark to override it.

**FluMut markers**, **SNPGenie diversity layers**, **WFABC selection results**, and
**export** to CSV, TSV, TXT, JSON, Markdown or VCF, every file carrying a header
recording the filters that produced it.

---

## Things worth knowing before you trust a number

These are properties of the data, not of this viewer, and they are not visible
from the tables themselves.

**LoFreq and GATK4 do not report the same quantity.** LoFreq's allele frequency is
an allele *fraction*. GATK4's is a *genotype* — a hom-alt call reads 1.0 whatever
the reads say. At one measured site LoFreq said 85.98%, GATK4 said 100%, and the
reads said 85.44%. FluLens reconciles this at load, giving GATK4 rows LoFreq's
fraction wherever both callers found the same change, and flagging the rest as
genotypes. Genotype-only calls are excluded from the consensus, because a genotype
cannot be read as "above 50%".

**Both callers emit a row for the same change**, so the table has two rows per
variant at most sites. Anything counting calls per codon has to dedupe on position
and alternative first.

**FluMut's HA and NA markers are H5/N1-numbered.** `HA1-5` means H5 HA1 numbering
and `NA-1` means N1 NA numbering, so on an H3N2 run those positions are read
against the wrong ruler. The internal genes — PB2, PB1, PA, NP, NS — are
subtype-agnostic and fine. FluLens detects the subtype from the reference segment
names and marks HA/NA findings it cannot confirm. A bare `A_HA` with no subtype
suffix counts as unconfirmable, not as a pass.

**Depth below 100 fakes fixations.** Low template input makes both callers report
false fixed differences, which is why Flumina's own floor is `min_depth 100` and
why any verdict here is downgraded below it.

**The assessment thresholds are absolute and are not the sidebar sliders.** Moving
a slider would otherwise change the verdict that the slider's own filter then
selects on. The panel says where a call sits relative to your current filters
without letting that move the verdict.

**In the SNPGenie layer, 99.1% of cells have one of πN or πS at zero.** At a single
codon in a single sample the value mostly reports *which kind* of difference was
seen, not which exceeds the other. The pattern across codons is the signal.

---

## Your data stays on your machine

FluLens has no server and makes no network requests. The browser version reads
your run folder through the file picker; the desktop version reads it from disk.
Nothing is uploaded, and the hosted page at `flu-crew.github.io` is a static file
that could be saved and run offline with no change in behaviour.

---

## Building from source

The browser version needs no build: `prototypes/flulens.html` *is* the
application. Edit it and reload.

The desktop app is a [Tauri](https://tauri.app) shell around that same file:

```bash
cargo install tauri-cli --version "^2" --locked
rustup target add aarch64-apple-darwin x86_64-apple-darwin
cd desktop && cargo tauri build --target universal-apple-darwin
```

> **The frontend is compiled into the binary.** Editing `prototypes/flulens.html`
> changes nothing the desktop app displays until you rebuild. This is the single
> most common way to lose an hour here.

Further reading, all of it maintainer-facing:

- [`docs/CONTRIBUTING.md`](docs/CONTRIBUTING.md) — how the two run loaders work,
  regenerating the example, and the things that reliably cost people an hour
- [`docs/RELEASING.md`](docs/RELEASING.md) — signing, notarisation and cutting a
  release. Read the DMG section before shipping one: Tauri notarises the app but
  not the disk image, so a green build can still produce a download macOS blocks
- [`desktop/README.md`](desktop/README.md) — the Tauri shell's design

---

## Citing

If FluLens contributed to published work, please cite it — see
[`CITATION.cff`](CITATION.cff). Please cite
[Flumina](https://github.com/flu-crew/Flumina) as well if you used the pipeline
that produced the data.

## License

GPL-3.0-or-later, the same as Flumina. See [`LICENSE`](LICENSE).
