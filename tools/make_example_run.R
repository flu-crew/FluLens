#!/usr/bin/env Rscript
#
# make_example_run.R — build example_run/, the demo dataset shipped with FluLens.
#
# The data is SYNTHETIC. No sample in it is a real animal, and no call in it is a
# real observation. It exists so that someone can open FluLens and see a populated
# screen without first running the pipeline, which is otherwise a hard prerequisite
# for evaluating a viewer.
#
# What is NOT synthetic is the scaffolding: the reference sequence and the curated
# site list are Flumina's own public files, and the product / codon annotation is
# computed by Flumina's own fluORFs.R. So the example exercises the real twelve-ORF
# layout, the real splice junctions, and the real column schema — a hand-written
# fixture would drift from the pipeline the moment either changed.
#
# Usage:
#   Rscript tools/make_example_run.R [--flumina /path/to/Flumina] [--out example_run]
#
# Needs a Flumina checkout for Scripts/fluORFs.R, reference.fa and
# curated_database.csv. You do not need to run this to use the example — the
# output is committed. It is here so the example's provenance is checkable.

suppressWarnings(suppressMessages(library(Biostrings)))

args <- commandArgs(trailingOnly = TRUE)
getarg <- function(flag, default) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) default else args[i + 1L]
}
FLUMINA <- path.expand(getarg("--flumina", "../Flumina"))
OUT     <- path.expand(getarg("--out", "example_run"))

for (f in c("Scripts/fluORFs.R", "reference.fa", "curated_database.csv")) {
  if (!file.exists(file.path(FLUMINA, f)))
    stop("Not a Flumina checkout (missing ", f, "): ", FLUMINA,
         "\nPass --flumina /path/to/Flumina")
}
source(file.path(FLUMINA, "Scripts", "fluORFs.R"))

set.seed(20260803)   # the whole file is reproducible from this line

#############################################
#### Reference and ORFs

ref <- flu_read_fasta(file.path(FLUMINA, "reference.fa"))
segs <- names(ref)
orfs <- lapply(segs, function(s) flu_orfs_for(s, ref[[s]]))
names(orfs) <- segs

cat("reference:", length(segs), "segments,",
    sum(vapply(orfs, length, 1L)), "products\n")

#############################################
#### Samples
#
# Four animals sampled on three days. A time series on purpose: it is what makes
# the sort-by-metadata columns worth clicking, and FluLens keys its WFABC section
# on (animal, locus, position), so a single-timepoint example would leave that
# section explaining its own absence.

animals <- c("A1", "A2", "A3", "A4")
days    <- c(1L, 3L, 5L)
groups  <- c(A1 = "Vaccinated", A2 = "Vaccinated", A3 = "Control", A4 = "Control")

samples <- expand.grid(Animal = animals, Day = days, stringsAsFactors = FALSE)
samples$sample <- sprintf("%s_D%d", samples$Animal, samples$Day)
samples$Group  <- unname(groups[samples$Animal])
samples <- samples[order(samples$Animal, samples$Day), ]
rownames(samples) <- NULL

# A4 is the thin library of the cohort. Its depths are scaled so the MEDIAN lands
# between Flumina's depth-100 floor and FluLens's default QC threshold of 200 —
# every call is one the pipeline would really have emitted, and the sample still
# fails QC. A demo where every sample passes teaches nothing about that column.
depth.scale <- c(A1 = 1.0, A2 = 0.85, A3 = 0.7, A4 = 0.035)

# A thin library also yields FEWER calls, not the same number at lower depth.
# Scaling depth alone would have left A4 with a full complement of variants,
# which is internally inconsistent in a way a reader would notice.
call.scale <- c(A1 = 1.0, A2 = 0.95, A3 = 0.85, A4 = 0.35)

cat("samples:", nrow(samples), "\n")

#############################################
#### Site pool
#
# Calls are drawn from a per-segment pool rather than sampled independently, so
# that sites RECUR across samples. Independent draws give a matrix of singletons,
# which looks plausible in a screenshot and hides the one pattern the view exists
# to show — a column of calls at the same codon.

pool <- do.call(rbind, lapply(segs, function(s) {
  n <- nchar(ref[[s]])
  # roughly one candidate site per 25 nt of segment
  p <- sort(sample(seq(20L, n - 20L), max(8L, round(n / 25))))
  data.frame(seg = s, position = p, stringsAsFactors = FALSE)
}))
# Recurrence weight: most sites appear in a few samples, a handful in most.
pool$weight <- rbeta(nrow(pool), 1.2, 5.5)
pool$weight[sample(nrow(pool), 6)] <- runif(6, 0.75, 0.97)   # the recurrent few

# Seed the curated sites explicitly. Left to chance, a 58-entry database against
# a few hundred random positions produces almost no hits, and the curated ruler
# ticks and the "curated only" site filter are then dead controls in the demo.
cds_to_seg <- function(o)
  unlist(mapply(function(s, e) s:e, o$exons$start, o$exons$end, SIMPLIFY = FALSE))

cur0 <- read.csv(file.path(FLUMINA, "curated_database.csv"),
                 na.strings = "", stringsAsFactors = FALSE)  # "NA" is a gene name
names(cur0)[1] <- "Gene"
seg.of <- vapply(segs, flu_segment_type, "")
seeded <- list()
for (r in seq_len(nrow(cur0))) {
  s <- names(seg.of)[match(cur0$Gene[r], seg.of)]
  if (is.na(s)) next
  o <- Filter(function(x) x$gene == flu_primary_product(cur0$Gene[r]), orfs[[s]])
  if (!length(o)) next
  map <- cds_to_seg(o[[1]])
  i <- (cur0$Amino_Acid[r] - 1L) * 3L + 1L
  if (i > length(map)) next            # position past this reference's protein
  seeded[[length(seeded) + 1L]] <-
    data.frame(seg = s, position = map[i], weight = runif(1, 0.30, 0.65),
               stringsAsFactors = FALSE)
}
pool <- rbind(pool, do.call(rbind, seeded))
pool <- pool[!duplicated(pool[c("seg", "position")]), ]
cat("site pool:", nrow(pool), "sites,", length(seeded), "seeded from the curated database\n")

#############################################
#### Draw calls

mk_alt <- function(b) {
  # transitions ~2x transversions, which is the usual influenza pattern
  ti <- c(A = "G", G = "A", C = "T", T = "C")
  tv <- list(A = c("C", "T"), G = c("C", "T"), C = c("A", "G"), T = c("A", "G"))
  if (runif(1) < 0.67) unname(ti[b]) else sample(tv[[b]], 1)
}

rows <- list()
for (i in seq_len(nrow(samples))) {
  sm <- samples$sample[i]
  ds <- depth.scale[[samples$Animal[i]]]
  cs <- call.scale[[samples$Animal[i]]]
  hit <- pool[runif(nrow(pool)) < pool$weight * cs, , drop = FALSE]
  if (!nrow(hit)) next
  for (j in seq_len(nrow(hit))) {
    s <- hit$seg[j]; p <- hit$position[j]
    rb <- substr(ref[[s]], p, p)
    if (!rb %in% c("A", "C", "G", "T")) next
    # Frequency: a mixture, not one log-normal. The bulk sits under 5% the way
    # real low-frequency calling does, but a pure log-normal puts so little mass
    # past 50% that the consensus view, the 50% threshold and the GATK4
    # corroboration below would all have nothing to act on.
    fq <- if (runif(1) < 0.11) runif(1, 0.55, 0.995)
          else min(0.49, exp(rnorm(1, log(0.018), 1.05)))
    # Floor at 100: Flumina's own min_depth. A real amino-acid table has no call
    # below it, so a demo that carries some would misrepresent the pipeline — and
    # would make FluLens's max-low-depth QC slider look live when on real output
    # it is inert by construction.
    dp <- max(100L, round(rlnorm(1, log(4200), 0.85) * ds))
    rows[[length(rows) + 1L]] <- data.frame(
      sample = sm, seg = s, position = p,
      reference = rb, alternative = mk_alt(rb),
      allele_frequency = round(fq, 6), depth = dp,
      stringsAsFactors = FALSE)
  }
}
calls <- do.call(rbind, rows)
# One row per (sample, site): a site drawn twice is the same change, not two.
calls <- calls[!duplicated(calls[c("sample", "seg", "position")]), ]
cat("draws:", nrow(calls), "\n")

#############################################
#### Annotate with Flumina's own code
#
# One row per (call, product). A position inside the M1/M2 overlap produces two
# rows, which is the pipeline's real behaviour and the reason the twelve-product
# grid has anything in its secondary bands.

ann <- flu_annotate_positions(calls, ref, locus.col = "seg", pos.col = "position")
ann <- ann[!is.na(ann$aa_position), ]
cat("annotated rows:", nrow(ann), "\n")

#############################################
#### Codons and translation

cds.cache <- new.env(parent = emptyenv())
cds_for <- function(seg, product) {
  key <- paste0(seg, "|", product)
  if (!is.null(cds.cache[[key]])) return(cds.cache[[key]])
  o <- Filter(function(x) x$gene == product, orfs[[seg]])[[1]]
  v <- flu_cds_seq(o, ref[[seg]])
  assign(key, v, envir = cds.cache)
  v
}

ann$reference_codon   <- NA_character_
ann$alternative_codon <- NA_character_
for (k in seq_len(nrow(ann))) {
  cds <- cds_for(ann$seg[k], ann$product[k])
  st  <- (ann$aa_position[k] - 1L) * 3L + 1L
  rc  <- substr(cds, st, st + 2L)
  if (nchar(rc) < 3L) next
  ac  <- rc
  substr(ac, ann$codon_position[k], ann$codon_position[k]) <- ann$alternative[k]
  ann$reference_codon[k]   <- rc
  ann$alternative_codon[k] <- ac
}
ann <- ann[!is.na(ann$reference_codon), ]

tr <- function(x) unname(ifelse(is.na(x), NA_character_, GENETIC_CODE[x]))
ann$reference_aa   <- tr(ann$reference_codon)
ann$alternative_aa <- tr(ann$alternative_codon)
ann$aa_changing    <- ifelse(ann$reference_aa == ann$alternative_aa, "NO", "YES")
ann <- ann[!is.na(ann$reference_aa) & !is.na(ann$alternative_aa), ]

#############################################
#### Caller rows
#
# Both callers fire at most sites in a real run, so the table carries two rows per
# change and anything counting calls has to dedupe on position+alternative. The
# example reproduces that, including the part that matters: GATK4's AF is a
# GENOTYPE (1.0 for hom-alt whatever the reads say), not an allele fraction.

ann$method           <- "LoFreq"
ann$locus            <- vapply(ann$seg, flu_segment_type, "")
ann$quality          <- pmin(3000L, round(rlnorm(nrow(ann), log(220), 0.9)))
ann$map_quality      <- "NA"
ann$af_type          <- "fraction"
ann$allele_fraction  <- ann$allele_frequency
ann$variant_id       <- sprintf("%s|%s|%d|%s", ann$sample, ann$seg,
                                ann$position, ann$alternative)

# GATK4 corroboration on the high-frequency calls, carrying a genotype AF.
hi  <- which(ann$allele_frequency > 0.5)
g4  <- ann[hi, , drop = FALSE]
if (nrow(g4)) {
  g4$method           <- "GATK4"
  g4$allele_frequency <- 1.0
  g4$af_type          <- "genotype"
  # allele_fraction stays LoFreq's — that is the column the pipeline added so
  # consumers stop having to infer the reconciliation themselves.
  g4$depth            <- round(g4$depth * runif(nrow(g4), 0.82, 0.98))
}

# A couple of GATK4-only calls in the thinnest library: no LoFreq counterpart, so
# no fraction to borrow and no read-level evidence. These are what FluLens grades
# "Cannot assess" rather than binning as artefacts.
only <- ann[ann$sample == "A4_D5", , drop = FALSE]
only <- head(only[order(-only$depth), , drop = FALSE], 2)
if (nrow(only)) {
  only$method           <- "GATK4"
  only$allele_frequency <- 1.0
  only$af_type          <- "genotype"
  only$allele_fraction  <- NA
  only$depth            <- c(6L, 4L)[seq_len(nrow(only))]
}

tab <- rbind(ann, g4, only)

#############################################
#### Write the amino-acid table

lead <- c("sample", "method", "locus", "product", "product_primary", "position",
          "reference", "alternative", "quality", "depth", "map_quality",
          "allele_frequency", "cds_position", "aa_position", "codon_position",
          "variant_id", "af_type", "allele_fraction", "reference_codon",
          "alternative_codon", "reference_aa", "alternative_aa", "aa_changing")
meta.cols <- c("Animal", "Day", "Group")
tab <- merge(tab, samples[c("sample", meta.cols)], by = "sample", all.x = TRUE)
tab <- tab[order(tab$sample, tab$locus, tab$position, tab$product), ]
tab$product_primary <- ifelse(tab$product_primary, "TRUE", "FALSE")
out <- tab[, c(lead, meta.cols)]

dir.create(file.path(OUT, "variant_analysis"), recursive = TRUE, showWarnings = FALSE)
write.table(out, file.path(OUT, "variant_analysis", "all_sample_amino_acids.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
cat("wrote all_sample_amino_acids.txt:", nrow(out), "rows\n")

#############################################
#### Curated sites
#
# Flumina's own curated_database.csv, joined the way outputSummary.R joins it:
# primary products only, because the database is in primary-product coordinates
# and letting it reach M2 / NEP manufactures coincidental hits.

cur <- cur0   # read once, above, where the sites are seeded into the pool
prim <- out[out$product_primary == "TRUE", ]
cj <- merge(prim, cur, by.x = c("locus", "aa_position"),
            by.y = c("Gene", "Amino_Acid"))
if (nrow(cj)) {
  cj <- cj[order(cj$sample, cj$locus, cj$aa_position), ]
  write.table(cj, file.path(OUT, "variant_analysis", "curated_amino_acids.txt"),
              sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
}
cat("wrote curated_amino_acids.txt:", nrow(cj), "rows\n")

#############################################
#### Per-sample LoFreq VCFs
#
# DP4 is what drives the trust panel, so it has to be internally consistent with
# depth and frequency rather than filler: alt counts must sum to depth*freq, and
# the reference allele needs its own strand split to compare against. A few sites
# are deliberately skewed so "Treat with caution" is reachable in the demo.

lf <- tab[tab$method == "LoFreq", ]
lf <- lf[!duplicated(lf[c("sample", "seg", "position")]), ]
for (sm in samples$sample) {
  v <- lf[lf$sample == sm, , drop = FALSE]
  d <- file.path(OUT, "vcf_files", sm)
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  con <- file(file.path(d, "lofreq-called-variants.vcf"), "w")
  writeLines(c(
    "##fileformat=VCFv4.0",
    "##fileDate=20260803",
    "##source=lofreq call -f reference.fa -o lofreq-called-variants.vcf final_mapped_reads.bam ",
    "##reference=reference.fa",
    "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Raw Depth\">",
    "##INFO=<ID=AF,Number=1,Type=Float,Description=\"Allele Frequency\">",
    "##INFO=<ID=SB,Number=1,Type=Integer,Description=\"Phred-scaled strand bias at this position\">",
    "##INFO=<ID=DP4,Number=4,Type=Integer,Description=\"Counts for ref-forward bases, ref-reverse, alt-forward and alt-reverse bases\">",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"), con)
  if (nrow(v)) {
    v <- v[order(v$seg, v$position), ]
    for (k in seq_len(nrow(v))) {
      dp <- v$depth[k]; af <- v$allele_frequency[k]
      nalt <- max(1L, round(dp * af)); nref <- max(0L, dp - nalt)
      skew <- if (runif(1) < 0.12) runif(1, 0.80, 0.95) else runif(1, 0.42, 0.58)
      af_f <- round(nalt * skew); af_r <- nalt - af_f
      rf_f <- round(nref * runif(1, 0.46, 0.54)); rf_r <- nref - rf_f
      sb <- round(abs(skew - 0.5) * 120)
      writeLines(sprintf("%s\t%d\t.\t%s\t%s\t%d\tPASS\tDP=%d;AF=%.6f;SB=%d;DP4=%d,%d,%d,%d",
                         v$seg[k], v$position[k], v$reference[k], v$alternative[k],
                         v$quality[k], dp, af, sb, rf_f, rf_r, af_f, af_r), con)
    }
  }
  close(con)
}
cat("wrote", nrow(samples), "LoFreq VCFs\n")

#############################################
#### IRMA read counts
#
# Drives the coverage strip. Spread across decades on purpose — a strip where
# every cell is the same shade demonstrates nothing. NS is the weak segment here,
# and A4 is missing two segments outright so the "absent" outline is reachable.

counts.by.sample <- list()   # so the coverage tables below agree with these

for (i in seq_len(nrow(samples))) {
  sm <- samples$sample[i]; ds <- depth.scale[[samples$Animal[i]]]
  d <- file.path(OUT, "IRMA_results", sm, "tables")
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  bias <- c(A_HA = 1.0, A_MP = 1.8, A_NA = 0.35, A_NP = 2.4,
            A_NS = 0.04, A_PA = 0.9, A_PB1 = 0.8, A_PB2 = 0.75)
  cnt <- round(rlnorm(8, log(90000), 0.7) * ds * bias[segs])
  names(cnt) <- segs
  if (sm %in% c("A4_D5", "A4_D3")) cnt[c("A_NS", "A_NA")] <- 0   # absent segments
  counts.by.sample[[sm]] <- cnt
  tot <- sum(cnt)
  ln <- c("Record\tReads\tPatterns\tPairsAndWidows",
          sprintf("1-initial\t%d\tNA\tNA", round(tot * 1.4)),
          sprintf("2-passQC\t%d\t%d\tNA", round(tot * 1.15), round(tot * 0.45)),
          sprintf("3-nomatch\t%d\t%d\tNA", round(tot * 0.02), round(tot * 0.015)),
          sprintf("3-match\t%d\t%d\t%d", tot, round(tot * 0.42), round(tot * 0.52)))
  for (s in segs) if (cnt[[s]] > 0)
    ln <- c(ln, sprintf("4-%s\t%d\t%d\t%d", s, cnt[[s]],
                        round(cnt[[s]] * 0.44), round(cnt[[s]] * 0.51)))
  writeLines(ln, file.path(d, "READ_COUNTS.txt"))
}
cat("wrote", nrow(samples), "READ_COUNTS.txt\n")

#############################################
#### IRMA per-base coverage
#
# tables/<SEGMENT>-coverage.txt, one row per reference base. This is what
# FluLens draws when a sample row is opened, and it is the only per-base depth
# available without a BAM.
#
# Deliberately written LAST. Everything above draws from the same seeded stream,
# so adding a generator earlier in the file would shift every number in the
# variant table and the VCFs; appending here leaves them byte-identical.
#
# Shape matters more than the absolute numbers: flat coverage would make the
# panel pointless. Each segment gets tapered ends, where real amplicon and
# assembly coverage always falls off, plus a couple of smooth interior dips.

for (i in seq_len(nrow(samples))) {
  sm <- samples$sample[i]; ds <- depth.scale[[samples$Animal[i]]]
  cnt <- counts.by.sample[[sm]]
  d <- file.path(OUT, "IRMA_results", sm, "tables")
  for (s in segs) {
    if (cnt[[s]] <= 0) next          # segment never assembled — no table at all
    n <- nchar(ref[[s]])
    # Centre the profile on the depth the variant table reports for this animal,
    # so the two agree rather than telling different stories about one library.
    mu <- 4200 * ds * (0.6 + 0.8 * (cnt[[s]] / max(cnt)))
    x <- seq_len(n) / n
    taper <- pmin(1, pmin(x, 1 - x) / 0.06)          # ends ramp up over ~6%
    wob <- 1 + 0.22 * sin(2 * pi * x * 3 + runif(1, 0, 6)) +
               0.14 * sin(2 * pi * x * 7 + runif(1, 0, 6))
    dep <- pmax(0, round(mu * taper * wob * rlnorm(n, 0, 0.06)))
    base <- strsplit(ref[[s]], "")[[1]]
    writeLines(c(
      paste("Reference_Name", "Position", "Coverage Depth", "Consensus",
            "Deletions", "Ambiguous", "Consensus_Count",
            "Consensus_Average_Quality", sep = "\t"),
      # One decimal on the quality. IRMA writes full float precision; carrying
      # fifteen digits of a number nothing reads costs ~0.5 MB across the example
      # for no information.
      sprintf("%s\t%d\t%d\t%s\t0\t0\t%d\t%.1f", s, seq_len(n), dep, base,
              dep, runif(n, 34, 38))),
      file.path(d, paste0(s, "-coverage.txt")))
  }
}
cat("wrote per-base coverage tables\n")

#############################################
#### Reference and GTF

file.copy(file.path(FLUMINA, "reference.fa"), file.path(OUT, "reference.fa"),
          overwrite = TRUE)

cfg <- file.path(tempdir(), "example_gtf.cfg")
writeLines(c(sprintf("REFERENCE_FILE=%s", normalizePath(file.path(OUT, "reference.fa"))),
             sprintf("OUTPUT_DIRECTORY=%s", normalizePath(OUT))), cfg)
st <- system2("Rscript", c(file.path(FLUMINA, "Scripts", "makeGTF.R"), cfg),
              stdout = TRUE, stderr = TRUE)
cat("makeGTF.R:", if (length(st)) tail(st, 1) else "ok", "\n")

cat("\nexample_run written to:", normalizePath(OUT), "\n")
cat("size:", system2("du", c("-sh", normalizePath(OUT)), stdout = TRUE), "\n")
