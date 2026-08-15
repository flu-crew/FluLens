# Example run — synthetic data

**Nothing in this directory is a real observation.** There is no animal called
`A1`, and no one measured any variant listed here. This directory lets you
open FluLens and see a full screen without a pipeline run first.

Do not cite it, and do not use it to check a result.

## What it is

A twelve-sample Flumina output directory. It is small enough to ship in the repository
— 26 MB. The variant data itself is 476 KB. The rest is the two
parts that must be large to be useful: per-base coverage (4.6 MB) and reads
(21 MB):

| | |
|---|---|
| samples | 12 — four animals (`A1`–`A4`) sampled on days 1, 3 and 5 |
| calls | 1,378 rows — 1,237 LoFreq, 141 GATK4 |
| products | all twelve, including spliced M2 / NEP / PA-X and frameshifted PB1-F2 |
| curated hits | 275 rows against Flumina's own curated site database |
| above consensus | 280 calls over 50% |
| BAMs | all 12 samples, capped at ~200× depth — see below |

## What is real about it

The scaffolding — the part that must stay in step with the pipeline:

- **`reference.fa` is Flumina's own reference**, copied without change.
- **Flumina's `makeGTF.R` makes `reference_gtf/`**, not a person by hand.
- **Flumina's `fluORFs.R` makes the product, CDS, and codon annotation.** It is the same
  code the pipeline uses. So the splice junctions in M2, NEP, and PA-X are the real
  ones, and the column schema cannot drift from real output.

Only the calls are invented — which samples, which positions, which frequencies, and
which depths. `../tools/make_example_run.R` regenerates the whole
directory from a fixed seed, so you can check that claim.

## What it was built to exercise

Each item below is a control. In a demo built from arbitrary data, each control would
have nothing to act on:

- **Both callers at one site.** GATK4 rows carry `af_type=genotype` with an AF of
  1.0, and LoFreq rows carry the fraction. FluLens reconciles 139 of them at load.
  This is why the table has two rows per change at most sites.
- **Two GATK4-only calls**, at depth 4 and 6, with no LoFreq counterpart and so no
  read-level evidence. These get the *Cannot assess* verdict.
- **A failing library.** `A4`'s three samples have median depth 124–147, under
  FluLens's default QC threshold of 200. So they mark ✕. But every LoFreq call
  still clears Flumina's own `min_depth` floor of 100. Zero LoFreq rows are below
  it, the same as real output.
- **Fewer calls in the thin library**, 51–58 against 114–156 elsewhere. If the depth
  alone were scaled, A4 would keep a full set of variants. A reader would notice that
  this is internally inconsistent.
- **Absent segments.** `A4_D3` and `A4_D5` recovered 6 of 8, so you can reach the
  coverage strip's "absent" outline.
- **Strand skew** on ~12% of VCF records. This is what puts calls into *Treat with
  caution* in the trust panel.
- **Recurrent sites**, so the matrix shows columns of calls and not a field of
  single calls.
- **Per-base depth**, in `IRMA_results/<sample>/tables/<SEGMENT>-coverage.txt`. This
  is what the depth profile draws when you open a row in the consensus view. The
  ends taper and the interior varies, because a flat profile would make that panel
  useless. A4's two unassembled segments have no coverage table, the same
  as IRMA would leave them. These tables are 4.6 MB of the 26 MB here. Delete
  them if you need only the variant grid.
- **Reads, for the pile-up** — `BAM_files/<sample>/final_mapped_reads.bam`, from
  `../tools/make_example_bams.py`, for all 12 samples. Every read is
  proper-paired (the set FluLens keeps, and the set that LoFreq's DP4 is counted
  over). Each read carries a quality profile that declines toward its 3′ end, and 12%
  have a soft-clipped tail. The alt allele at each called site is drawn at that
  record's own AF and split across the strands by its own DP4. So the pile-up,
  the frequency in the table, and the strand-balance panel are three views of one
  number, not three guesses. The AF is verified against `samtools
  mpileup` across every called site (median error ~0, mean |error| ~0.01).
  The base qualities are drawn so all three Phred bands have values — about
  4 / 8 / 88 % at ≥30 / 20–29 / <20, against 11 / 3 / 86 on the real run this was
  tuned against. Alt-supporting bases are high and sequencing errors are low. This
  is what makes the pile-up's `fade marks by base quality` control show
  something.

## Depth is capped, so the BAM is a downsample

Every BAM is capped at **~200× per-base depth**, so all twelve ship in ~21 MB.
That is well below the real value: `A1`–`A3` carry LoFreq DP into the thousands
(up to 33,821), and a BAM true to that would be ~86 MB *per sample*. `A4` is
the deliberately thin library (DP 100–764) and is barely capped.

`A4_D3` and `A4_D5` carry 6 of 8 segments, which matches what IRMA assembled for
them. So you can reach the pile-up's absent-segment path from the example.

## What it does not include

No FluMut markers and no WFABC selection results. Both sections handle their own
absence and tell you how to enable them. 

## Note on the reference

Flumina's reference uses bare segment names (`A_HA`, `A_NA`) with no subtype
suffix. So FluLens reports **HA/NA FluMut markers as unconfirmable** on this data.
That is the intended conservative behaviour. FluLens cannot read the subtype from the
names, so it treats an unconfirmable position as a mismatch, not as a pass.
You will see that warning in the sidebar. It is correct.
