# FluLens desktop shell

A small Tauri wrapper around `flulens.html`. The app itself does not
change. This shell removes three limits of the browser version.

**Why it exists**

- **The folder picker forgets.** In the browser, you must choose the run directory
  again each session, because the browser cannot persist a `File` handle. The
  shell remembers the last directory and reopens it.
- **`file://` cannot fetch.** So the dev path needs a web server. Under macOS TCC,
  that server could not read `~/Documents` when the repo was there. For this reason,
  development runs against a copy in a scratch directory, not the real
  file. The repo has since moved to `~/Bioinformatics/Github/FluLens`. The
  scratch copy stayed, because a server on the canonical file tests whatever
  saved state the file is in.
- **Ranged BAM reads.** The browser path depends on the dev server to honour
  `Range`. Native file reads have no such dependency.

**What it deliberately does NOT do**

It does not rewrite the app. The shell gives the same `files` map of objects with
`.text()`, `.bytes(a,b)`, and `.size()` that `devBootstrap` fakes and the
directory picker builds for real. So every code path below it is the one that
already works. If this shim and the picker ever disagree, the shim is wrong.

**Build**

    cargo install tauri-cli --version "^2" --locked   # once
    cd desktop && cargo tauri build

The output is `target/release/bundle/macos/FluLens.app`. Drag it to `/Applications`.
`cargo` needs `. "$HOME/.cargo/env"` in `~/.zshrc`. rustup wrote to `~/.profile`,
which zsh does not read, and it could not amend the root-owned `~/.bash_profile`.

**Four things about this build that cost time**

1. **The binary contains the compiled frontend.** If you edit `flulens.html`,
   the app shows no change. If you copy it into `dist/` by hand while the app runs,
   the app also shows no change. Every frontend change needs a rebuild. People
   mistook this for an app bug one time. The letters scheme showed Zappo in the app
   while the browser correctly showed ClustalX, because the embedded copy was stale.
2. **`beforeBuildCommand` runs from the PROJECT ROOT**, the parent of this
   directory — not from here, where `tauri.conf.json` is. So its paths
   read `prototypes/…` and `desktop/dist/…`, and that is correct. `frontendDist`
   on the line above it *is* relative to this file. Two adjacent keys use two
   different bases.
3. **`frontendDist` must hold nothing but the frontend.** At first it was `../prototypes`,
   which carries the `run` symlink to a 28 GB results tree. Tauri walks
   the directory to embed it, and the build was SIGKILLed four times. It looks
   exactly like an OOM kill or a sandbox kill, but it is neither. With `dist/`, it
   compiles in ~17 s.
4. **`bundle.active` must be true to get an icon.** If it is false, the build still
   makes a working binary, so nothing looks wrong. But there is no `.app`, and
   therefore no `Info.plist` and no Dock icon.

**Icon:** `mkicon.py` makes it (stdlib only — no PIL, ImageMagick, or node on
this machine). Then `cargo tauri icon icon.png` expands it into `icons/`.

**Status:** built and running, 2026-08-03 — see the FluLens handoff.
