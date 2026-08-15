# Example run — synthetic Nanopore data

**Nothing in this directory is a real observation.** The data is synthetic. Do
not cite it or use it to check a result.

## What it is

A five-sample FluPore output directory converted for FluLens. It is small enough
to ship in the repository — 3.6 MB.

| | |
|---|---|
| samples | 5 — `barcode01` through `barcode05` |
| segments | all 8 influenza A segments |
| caller | iVar (the variant caller FluPore uses) |
| reads | single-end ONT, ~60–80× depth per segment |

## How it was made

Synthetic reads were generated from Flumina's reference and run through the
FluPore pipeline. The output was converted into the Flumina-shaped directory
layout that FluLens expects. The BAMs, coverage tables, variant calls, and VCFs
are all present.

## How to use it

Load this folder the same way you load a Flumina example:

- **In the browser:** go to `https://flu-crew.github.io/FluLens/?run=example_run_nanopore`
- **From the file picker:** click *Open run folder…* and select this directory
