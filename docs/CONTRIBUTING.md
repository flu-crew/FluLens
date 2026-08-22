# Developing FluLens

> Every path and command below is relative to the **repository root**, not to
> `docs/` where this file is.

## The application is one file

`flulens.html` — about 5,000 lines, no dependencies, no build step,
and no framework. Open it in a browser and it works. Everything else in this
repository is there to package, document, or test that file.

Keep it this way. Because FluLens is one file, you can email it to a
colleague, open it on a locked-down machine, or store it with a paper. It will
still run in ten years.

## Two ways to load a run

Both make the same `files` map. Each object has `.text()`, `.bytes(a,b)`, and
`.size()`. So all the readers below them are shared:

- **The directory picker** (`<input webkitdirectory>`) — the entry point in
  a browser.
- **`tauriBootstrap()`** — the desktop shell, which reads from disk.

There is a third loader for development only: `devBootstrap()`. It loads a run
over HTTP without a picker. It runs only when the page is **served over http**. Under
`file://`, it never runs.

### Local development symlinks

`devBootstrap()` reads from `prototypes/run`. You must create it. It is
not committed, because it points to your own data:

```bash
ln -s /path/to/a/flumina/run   prototypes/run
ln -s /path/to/metadata.csv    prototypes/meta.csv    # optional
```

Then serve the folder and open it. **The pane accepts only `127.0.0.1`** — not
`localhost`, and not `file://`:

```bash
python3 -m http.server 8791 --directory prototypes
```

If you already loaded a run, a new symlink target is not enough on its own. The dev
fetches use the same URLs, so the browser could serve the tables of the previous run
from the cache. For this reason, `devBootstrap` fetches with `cache: 'no-store'`.

### Static hosting needs a manifest

`devBootstrap` finds GTFs, per-sample VCFs, and read counts when it **lists
directories**. GitHub Pages serves no directory index. So a deployed run would
lose the twelve-product grid, the trust panel, and the coverage strip, but
still look like it loaded. For this reason, a run served statically carries
`manifest.json`:

```bash
python3 tools/make_manifest.py example_run
```

Without a manifest, FluLens uses the autoindex path, and nothing changes.

## The desktop app embeds the frontend

`beforeBuildCommand` copies `flulens.html` into `desktop/dist/`, and
Tauri **compiles that directory into the binary**. So:

> If you edit `flulens.html`, the desktop app shows no change. If you
> copy it into `dist/` by hand, that also shows no change. Every frontend change needs
> `cargo tauri build`.

People have mistaken this for an application bug at least one time. The letters scheme
showed one palette in the app and a different palette in the browser, because the
embedded copy was stale.

## Regenerating the example run

`example_run/` is synthetic. FluLens makes it from a fixed seed against Flumina's own
public reference and ORF code:

```bash
Rscript tools/make_example_run.R --flumina /path/to/Flumina
python3 tools/make_example_bams.py --out example_run     # needs samtools
python3 tools/make_manifest.py example_run
```

It needs R with Biostrings, and a Flumina checkout for `Scripts/fluORFs.R`,
`reference.fa`, and `curated_database.csv`. The scripts use the pipeline's own
annotation code on purpose. A hand-written fixture would drift from Flumina as soon as
one of them changed.

Run `make_example_bams.py` **after** the R script, not before. It reads the
VCFs and IRMA coverage tables that the R script writes, and it makes every read from
them. This keeps the pile-up in agreement with the table beside it. It
writes all twelve samples, capped at ~200× depth (`DEPTH_CAP`), so they ship in
~21 MB — see `example_run/README.md`.

## Problems to avoid

**Declare every persisted preference in the one early block**, above the UI code
that restores it from localStorage. If you declare state below its use, the restore
assignment falls in a temporal dead zone inside a `try/catch` that hides the error.
Then one of two things happens: the next unguarded read throws and stops the whole
script, or nothing throws and the setting never restores. The second one is worse. If
something strange happens at load, check `typeof ready` and `typeof cw`
first. A part-executed script is the sign.

**A change to a default does not change what is already stored.** Read persisted
booleans so that an *absent* value means "use the current default", and only an
explicit stored value wins. A `<select>` falls back to its first option unless
something sets `.value`. So mark the default option `selected` in the HTML and
keep it in step with the JS default. If you do not, the control disagrees with the
view, which is worse than one of them alone being wrong.

**Never do per-row or per-column work inside the draw loop.** The matrix is about
3% dense. Iterate the call list, not the grid. If you sample every Nth column — like
a dense-alignment renderer — you draw a fraction of the calls, and it still looks
correct.

**`requestAnimationFrame` never fires in a hidden tab.** Any await that gives the
browser a frame to paint must also race a timeout. If it does not, a run opened in a
background tab deadlocks the load with no error and no timeout.

**Nextflow's `work/` is inside the run directory and competes for every
lookup.** A file lookup takes the first path that matches, so staged copies are not
inert. FluLens excludes anything under `work/` or `.nextflow/` when it indexes. It
reports the count in the sidebar, so you do not mistake a large exclusion for the
app that fails to find files.

**"NA" is a gene name.** Anything that round-trips through a CSV needs
`na.strings = ""`. If not, every neuraminidase row vanishes.

**`cssW` reads 0 in an automated browser pane** until something forces layout.
That collapses `bodyW()` to about 50px and truncates anything computed from the
visible column range. Call `resize()` first when you drive the app programmatically.
