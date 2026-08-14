# Developing FluLens

> Every path and command below is relative to the **repository root**, not to
> `docs/` where this file lives.

## The application is one file

`prototypes/flulens.html` — around 5,000 lines, no dependencies, no build step,
no framework. Open it in a browser and it works. Everything else in this
repository exists to package, document or test that file.

Keep it that way. The single-file property is why FluLens can be emailed to a
collaborator, opened on a locked-down machine, or archived alongside a paper and
still run in ten years.

## Two ways to load a run

Both build the same `files` map of objects with `.text()`, `.bytes(a,b)` and
`.size()`, so every reader below them is shared:

- **The directory picker** (`<input webkitdirectory>`) — the real entry point in
  a browser.
- **`tauriBootstrap()`** — the desktop shell, reading from disk.

There is a third for development only, `devBootstrap()`, which loads a run over
HTTP without a picker. It fires only when the page is **served over http**; under
`file://` it never runs.

### Local development symlinks

`devBootstrap()` reads from `prototypes/run`. Create it yourself — it is
deliberately not committed, because it points into whatever data you happen to
have:

```bash
ln -s /path/to/a/flumina/run   prototypes/run
ln -s /path/to/metadata.csv    prototypes/meta.csv    # optional
```

Then serve the directory and open it. **The pane only accepts `127.0.0.1`** — not
`localhost`, not `file://`:

```bash
python3 -m http.server 8791 --directory prototypes
```

Repointing the symlink is not enough on its own if you have loaded a run already;
the dev fetches would otherwise hit the same URLs and the browser would serve the
previous run's tables from cache. `devBootstrap` fetches with `cache: 'no-store'`
for exactly this reason.

### Static hosting needs a manifest

`devBootstrap` discovers GTFs, per-sample VCFs and read counts by **listing
directories**. GitHub Pages serves no directory index, so a deployed run would
quietly lose the twelve-product grid, the trust panel and the coverage strip
while still looking like it loaded. A run served statically therefore carries
`manifest.json`:

```bash
python3 tools/make_manifest.py example_run
```

Absent a manifest, the autoindex path is used and nothing changes.

## The desktop app embeds the frontend

`beforeBuildCommand` copies `prototypes/flulens.html` into `desktop/dist/` and
Tauri **compiles that directory into the binary**. So:

> Editing `prototypes/flulens.html` changes nothing the desktop app shows, and
> neither does copying it into `dist/` by hand. Every frontend change needs
> `cargo tauri build`.

This has been mistaken for an application bug at least once — the letters scheme
rendering one palette in the app and another in the browser, from a stale
embedded copy.

## Regenerating the example run

`example_run/` is synthetic, generated from a fixed seed against Flumina's own
public reference and ORF code:

```bash
Rscript tools/make_example_run.R --flumina /path/to/Flumina
python3 tools/make_example_bams.py --out example_run     # needs samtools
python3 tools/make_manifest.py example_run
```

It needs R with Biostrings, and a Flumina checkout for `Scripts/fluORFs.R`,
`reference.fa` and `curated_database.csv`. Using the pipeline's own annotation
code is deliberate: a hand-written fixture would drift from Flumina the moment
either changed.

`make_example_bams.py` must run **after** the R script, not before: it reads the
VCFs and IRMA coverage tables that script writes and derives every read from
them, which is what keeps the pile-up agreeing with the table beside it. It only
writes `A4`'s three samples — see `example_run/README.md` for why.

## Things that will bite you

**Declare every persisted preference in the one early block**, above the UI wiring
that restores it from localStorage. Declaring state below its use puts the restore
assignment in a temporal dead zone inside a `try/catch` that swallows the error —
and then either the next unguarded read throws and halts the whole script, or
nothing throws and the setting simply never restores. The second is worse. If
something inexplicable happens at load, check `typeof ready` and `typeof cw`
first; a partially executed script is the signature.

**Changing a default does not change what is already stored.** Read persisted
booleans so that *absent* means "take the current default" and only an explicitly
stored value wins. And a `<select>` falls back to its first option unless
something assigns `.value`, so mark the default option `selected` in the HTML and
keep it in step with the JS default — otherwise the control disagrees with the
rendering, which is worse than either being wrong alone.

**Never do per-row or per-column work inside the draw loop.** The matrix is about
3% dense: iterate the call list, not the grid. Sampling every Nth column — what a
dense-alignment renderer does — draws a fraction of the calls and looks entirely
plausible.

**`requestAnimationFrame` never fires in a hidden tab.** Any "give the browser a
frame to paint" await must race a timeout, or opening a run in a background tab
deadlocks the load with no error and no timeout.

**Nextflow's `work/` sits inside the run directory and competes for every
lookup.** File lookups take the first path that matches, so staged copies are not
inert. Anything under `work/` or `.nextflow/` is excluded when indexing, and the
count is reported in the sidebar so a large exclusion is never mistaken for the
app failing to find things.

**"NA" is a gene name.** Anything round-tripping through a CSV needs
`na.strings = ""`, or every neuraminidase row silently vanishes.

**`cssW` reads 0 in an automated browser pane** until something forces layout,
which collapses `bodyW()` to about 50px and truncates anything computed from the
visible column range. Call `resize()` first when driving the app programmatically.
