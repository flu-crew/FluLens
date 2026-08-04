#!/usr/bin/env python3
"""Write manifest.json for a run directory served over static hosting.

FluLens discovers a run by walking it — the GTFs, the per-sample VCFs and the
IRMA read counts are all found by listing directories, not by probing known
paths. A static host (GitHub Pages) serves no directory index, so that walk
returns nothing and a deployed run silently loses the twelve-product grid, the
trust panel and the coverage strip while still looking like it loaded.

manifest.json is a flat array of every path inside the run, relative to its
root. When present, FluLens reads directory listings from it instead.

    python3 tools/make_manifest.py example_run
"""
import json
import os
import sys

SKIP = {".DS_Store", "manifest.json"}


def main() -> int:
    if len(sys.argv) != 2:
        print(__doc__, file=sys.stderr)
        return 2
    root = sys.argv[1]
    if not os.path.isdir(root):
        print(f"not a directory: {root}", file=sys.stderr)
        return 1

    paths = []
    for dirpath, dirnames, filenames in os.walk(root):
        # Nextflow's scratch lives inside the run directory on some layouts and
        # holds a staged copy of every published file. FluLens excludes it when
        # indexing a real run; a manifest that carried it would put the copies
        # back and let them win a lookup.
        dirnames[:] = [d for d in dirnames if d not in (".nextflow", "work")]
        for fn in filenames:
            if fn in SKIP:
                continue
            rel = os.path.relpath(os.path.join(dirpath, fn), root)
            paths.append(rel.replace(os.sep, "/"))

    paths.sort()
    out = os.path.join(root, "manifest.json")
    with open(out, "w") as fh:
        json.dump(paths, fh, indent=0, separators=(",", ":"))
        fh.write("\n")
    print(f"{out}: {len(paths)} paths")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
