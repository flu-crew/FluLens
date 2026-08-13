# Read-level probes

Written 2026-08-13 to settle item 0d — six low-frequency calls IRMA made that
none of LoFreq, iVar or GATK4 reported. They read a run directory in place and
print summaries; none of them writes into the run except `depth_band.py`, which
takes an explicit output path.

All of them need `samtools` on `PATH`. `<run_dir>` is a Flumina output
directory (the one holding `BAM_files/`, `IRMA_results/`, `reference.fa`).

| script | question it answers |
|---|---|
| `pileup_probe.py` | at one position, per allele: count, strand, base quality, MAPQ, distance from the read end, soft-clip and indel proximity |
| `irma_allele_names.py` | which reads does IRMA assign a given allele — feeds the name file below |
| `trace_reads.py` | what did BWA do with those same reads: aligned, soft-clipped, duplicate, elsewhere, unmapped |
| `clip_audit.py` | lay a soft-clipped tail back on the reference and score it — is this sequence BWA declined to align, or the far side of a template switch |
| `junction_test.py` | do the minority reads terminate at the variant more often than the consensus reads do |
| `side_by_side.py` | one read, both alignments, one window of reference — the ground truth behind the statistics |
| `pileup_views.py` | the same position as `samtools depth`, iVar and LoFreq are each shown it |
| `baq_ramp.py` | how much usable depth BAQ costs as a function of distance from a contig edge |
| `depth_band.py` | every position whose raw depth clears `MIN_DEPTH` but whose iVar-visible depth does not |

## Two things that produce convincing wrong answers

(Details of an unpublished run were removed here.)





(Details of an unpublished run were removed here.)





## Reproducing item 0d

```bash
RUN=~/path/to/Analysis_New/WGS

python3 junction_test.py "$RUN"

(Details of an unpublished run were removed here.)



(Details of an unpublished run were removed here.)


python3 depth_band.py "$RUN" "$RUN/depth_profiles/mindepth_blind_band.tsv"
```

(Details of an unpublished run were removed here.)

