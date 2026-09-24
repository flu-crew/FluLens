# Read-level probes

Written 2026-08-13 to settle item 0d — low-frequency calls that IRMA made and
that none of LoFreq, iVar, or GATK4 reported. They read a run directory in place and
print summaries. None of them writes into the run, except `depth_band.py`, which
takes an explicit output path.

All of them need `samtools` on `PATH`. `<run_dir>` is a Flumina output
directory (the one that holds `BAM_files/`, `IRMA_results/`, and `reference.fa`).

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

**Read names.** IRMA appends the fastq comment as `_3:N:0:INDEX` (the `3` is
its marker for a merged pair), and BWA keeps the bare Illumina name. If you compare
them without normalisation, you find none of the reads, and this looks like a result. Every
script here splits on `_` and takes the first field. Illumina names contain no
underscore.

**Negative region starts.** A window built as `pos-400` goes below 1 for a
position near a segment start. Then `samtools view` returns nothing, and the empty
result looks like a clean negative — "BWA does nothing with these reads". For a short
time, it produced exactly that for calls near a segment start. Now both
`clip_audit.py` and `junction_test.py` clamp to 1.

## Reproducing item 0d

```bash
RUN=~/path/to/Analysis_New/WGS
# One call per line, tab-separated: sample, locus, position, ref, alt.
# Keep this file outside the repository.
CALLS=~/path/to/calls.tsv

python3 junction_test.py "$RUN" "$CALLS"

python3 irma_allele_names.py "$RUN/IRMA_results/<sample>/A_PA.bam" A_PA <position> <alt> > /tmp/names.txt
python3 clip_audit.py "$RUN/BAM_files/<sample>/final_mapped_reads.bam" \
    "$RUN/reference.fa" A_PA <position> <alt> /tmp/names.txt

python3 pileup_views.py "$RUN/BAM_files/<sample>/final_mapped_reads.bam" \
    "$RUN/reference.fa" A_PB2:<position>-<position> 30

python3 depth_band.py "$RUN" "$RUN/depth_profiles/mindepth_blind_band.tsv"
```

`depth_band.py` reads every sample, so it takes a long time on a large run.

## Reproducing the LoFreq `-B` measurement

Does `-B` in `lofreq call` recover real variants, or does it make false ones?
`baq_support.py` describes the calls that `-B` adds. `baq_baseline.py` compares
them against the calls LoFreq already makes. This comparison is the only way a support
rate has meaning.

LoFreq is not installed on the laptop, so the calling runs in the image. The
run directory is mounted read-only — nothing here writes into it.

```bash
RUN=~/path/to/Analysis_New/WGS
OUT=/tmp/baq; mkdir -p "$OUT"

# Deterministic spread, not a hand-picked set.
ls "$RUN/BAM_files" | sort | awk 'NR%12==1' > "$OUT/samples.txt"

docker run --rm --platform linux/amd64 \
  -v "$RUN":/data:ro -v "$OUT":/out \
  --entrypoint bash chutter/flumina:latest -c '
while read s; do
  B=/data/BAM_files/$s/final_mapped_reads.bam
  lofreq call    -f /data/reference.fa -o /out/$s.extbaq.vcf $B 2>/dev/null
  lofreq call -B -f /data/reference.fa -o /out/$s.nobaq.vcf  $B 2>/dev/null
done < /out/samples.txt'

python3 baq_support.py  "$RUN" "$OUT" $(tr '\n' ' ' < "$OUT/samples.txt")
python3 baq_baseline.py "$RUN" "$OUT" $(tr '\n' ' ' < "$OUT/samples.txt")
```

The two scripts print the counts that this comparison needs.

`--entrypoint bash` is required: the image's entrypoint is the Flumina launcher.
So a bare `docker run ... lofreq` prints the launcher's help and does not run
LoFreq.

**Indel positions come from iVar's TSV, never from the GATK4 indel VCF.** GATK4
reports almost no indels, because it is a genotype caller and does not see indels
below genotype frequency. If you use it, BAQ looks like it protects against nothing.
