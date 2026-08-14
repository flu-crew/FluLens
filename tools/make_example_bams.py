#!/usr/bin/env python3
"""Synthetic per-sample BAMs for example_run/, so the read pile-up has something
to open.

Only A4's three samples get one, and that is a consistency constraint rather
than a shortcut. The pile-up puts the caller's depth and the BAM's depth on
screen at the same time, and this app has been bitten before by two numbers with
the same name disagreeing. A1-A3 carry LoFreq DP of 2,300-4,900 (max 35,140);
a BAM honest to that is ~86 MB for ONE sample. A4's DP is 118-161, which is
1.6 Mbp of alignment per sample and ships.

Everything here is derived from files already in example_run, so the BAM cannot
drift from the tables beside it:

  depth   IRMA's own per-base coverage table, scaled by BWA_FACTOR. The
          pile-up's BAM depth SHOULD sit above the caller's -- the callers
          filter, the BAM does not -- and 1.3 is the ratio the swine run shows.
  alleles every LoFreq call, at its own AF, split across strands by its own
          DP4. The strand-balance panel therefore agrees with the VCF instead
          of being a second opinion about it.
  segments whichever segments IRMA assembled, so A4_D3 and A4_D5 stay at 6 of 8
          and the pile-up's absent-segment path is reachable.

Usage:  python3 tools/make_example_bams.py [--out example_run]
Needs samtools on PATH. Deterministic: same seed, same bytes.
"""

import argparse, os, random, subprocess, sys
from pathlib import Path

SAMPLES = ["A4_D1", "A4_D3", "A4_D5"]
READ_LEN = 150
INSERT_MEAN, INSERT_SD = 300, 45
BWA_FACTOR = 1.3          # BAM depth over IRMA depth; see module docstring
CALLER_HEADROOM = 1.15    # ... and over the caller's own DP, where that is higher
CALL_BUMP_W = 260         # half-width of the bump a call raises, in bases
CLIP_MARGIN = 1.10        # coverage is counted before clipping shortens reads
ERR_RATE = 0.002          # per base, and deliberately low-quality when it fires
CLIP_FRAC = 0.12          # reads carrying a soft-clipped tail
DEL_FRAC = 0.008          # reads carrying a short deletion
SEED = 20260814


def read_fasta(path):
    seqs, name, buf = {}, None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if name:
                    seqs[name] = "".join(buf)
                name, buf = line[1:].split()[0], []
            elif line:
                buf.append(line.upper())
    if name:
        seqs[name] = "".join(buf)
    return seqs


def read_coverage(path):
    """IRMA's per-base depth, as a 1-based list indexed [0] = position 1."""
    depth = []
    with open(path) as fh:
        next(fh)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) > 2:
                depth.append(int(f[2]))
    return depth


def read_vcf(path):
    """{segment: {pos1: (ref, alt, af, p_alt_fwd, p_alt_rev, dp)}}

    The two per-strand probabilities come straight from DP4, so a record with
    17 alt on the forward strand and 2 on the reverse reproduces that skew
    rather than an average of it -- which is the whole point of the strand
    panel."""
    calls = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 8:
                continue
            seg, pos, ref, alt, info = f[0], int(f[1]), f[3], f[4], f[7]
            if len(ref) != 1 or len(alt) != 1:
                continue          # indels are not modelled; none in this set
            kv = dict(p.split("=", 1) for p in info.split(";") if "=" in p)
            af = float(kv.get("AF", 0))
            dp = int(kv.get("DP", 0))
            pf = pr = af
            if "DP4" in kv:
                rf, rr, af_, ar = (int(x) for x in kv["DP4"].split(","))
                if rf + af_ > 0:
                    pf = af_ / (rf + af_)
                if rr + ar > 0:
                    pr = ar / (rr + ar)
            calls.setdefault(seg, {})[pos] = (ref, alt, af, pf, pr, dp)
    return calls


def qual_profile(n, rng):
    """Phred along one read, 5' to 3'. High and flat, then a decline over the
    last third -- the shape every Illumina read has, and the reason the read
    panel says which end is which."""
    out = []
    for i in range(n):
        frac = i / max(1, n - 1)
        base = 38.0 - 7.0 * max(0.0, (frac - 0.6) / 0.4) ** 1.6
        out.append(int(max(20, min(41, rng.gauss(base, 1.6)))))
    return out


def depth_target(irma, calls, seg_len):
    """The per-base depth to aim for.

    IRMA's shape scaled up, floored by what the caller claims. The two sources
    in example_run were generated independently and disagree: 40-50% of A4's
    calls carry a LoFreq DP above IRMA's depth x1.3, by as much as 199x. Left
    alone, the pile-up would print "caller depth 579" over a coverage strip
    reading 140 and look like the depth bug this app has already had twice.
    The BAM has to be the wider number, since the callers filter and it does
    not, so each call raises a triangular bump around itself rather than a
    step -- a rectangular floor would draw a cliff in the coverage strip."""
    target = [d * BWA_FACTOR for d in irma]
    for pos, (_, _, _, _, _, dp) in calls.items():
        need = dp * CALLER_HEADROOM
        p0 = pos - 1
        for p in range(max(0, p0 - CALL_BUMP_W), min(seg_len, p0 + CALL_BUMP_W)):
            v = need * (1 - abs(p - p0) / CALL_BUMP_W)
            if v > target[p]:
                target[p] = v
    return target


def pick_starts(depth, rng, n_frag):
    """Fragment starts drawn in proportion to the target depth, so the profile
    the BAM ends up with has IRMA's shape rather than a flat one."""
    total = sum(depth)
    if total <= 0:
        return []
    cum, acc = [], 0
    for d in depth:
        acc += d
        cum.append(acc / total)
    starts = []
    for _ in range(n_frag):
        u = rng.random()
        lo, hi = 0, len(cum) - 1
        while lo < hi:
            mid = (lo + hi) // 2
            if cum[mid] < u:
                lo = mid + 1
            else:
                hi = mid
        starts.append(lo)
    return starts


def build_mate(seq, ref, start, length, calls, rng, is_rev):
    """One mate. Returns (pos0, cigar, seq, quals) in REFERENCE orientation.

    Everything that makes a read worth looking at individually is decided here:
    which allele it carries at a called site, where the aligner gave up on its
    tail, and what quality it was read at."""
    n = len(ref)
    clip_l = clip_r = 0
    if rng.random() < CLIP_FRAC:
        c = rng.randint(5, 40)
        if rng.random() < 0.5:
            clip_l = c
        else:
            clip_r = c
    aln_len = length - clip_l - clip_r
    if start + aln_len > n:
        aln_len = n - start
    if aln_len < 40:
        return None

    del_at, del_len = -1, 0
    if rng.random() < DEL_FRAC and aln_len > 60:
        del_len = rng.randint(1, 3)
        del_at = rng.randint(20, aln_len - 30)

    bases, quals = [], qual_profile(length, rng)
    if is_rev:
        quals = quals[::-1]        # 5' end is on the right in reference space
    qi = 0
    # soft-clipped bases are real sequence the aligner discarded, so they are
    # mostly wrong against the reference -- that is why they were clipped
    for _ in range(clip_l):
        bases.append(rng.choice("ACGT"))
        qi += 1
    rp = start
    consumed = 0
    while consumed < aln_len:
        if rp == start + del_at and del_len:
            rp += del_len
            consumed += del_len
            continue
        b = ref[rp] if rp < n else "N"
        call = calls.get(rp + 1)
        if call:
            _, alt, _, pf, pr, _ = call
            if rng.random() < (pr if is_rev else pf):
                b = alt
                # Alt-supporting bases are well read -- they have to be, since
                # DP4 counts what cleared LoFreq's BQ floor -- but not uniformly
                # so. A hard randint(31,40) here left the 20-29 Phred band
                # completely EMPTY, which made the pile-up's quality fade look
                # like a two-state switch and gave the quality overlay a band it
                # could never paint.
                quals[qi] = max(quals[qi], int(min(41, max(28, rng.gauss(35, 3)))))
        elif rng.random() < ERR_RATE:
            b = rng.choice([x for x in "ACGT" if x != b])
            # An error is a low-quality base, with a thin tail into the middle
            # band rather than a flat 10-19 block. Tuned against the real swine
            # run, where the three Phred bands hold 11% / 3% / 86% of the
            # mismatches: the middle one is a minority, not a third of them.
            quals[qi] = int(min(29, max(6, rng.gauss(14, 4.5))))
        bases.append(b)
        rp += 1
        consumed += 1
        qi += 1
    for _ in range(clip_r):
        bases.append(rng.choice("ACGT"))
        qi += 1

    m_len = aln_len - del_len
    if del_at >= 0 and del_len:
        cig = []
        if clip_l:
            cig.append(f"{clip_l}S")
        cig.append(f"{del_at}M{del_len}D{m_len - del_at}M")
        if clip_r:
            cig.append(f"{clip_r}S")
        cigar = "".join(cig)
    else:
        cigar = (f"{clip_l}S" if clip_l else "") + f"{m_len}M" + (f"{clip_r}S" if clip_r else "")
    return start, cigar, "".join(bases), quals[: len(bases)]


def sam_records(sample, refs, cov_dir, calls, rng):
    yield "@HD\tVN:1.6\tSO:unsorted"
    for name in refs:
        yield f"@SQ\tSN:{name}\tLN:{len(refs[name])}"
    yield f"@RG\tID:{sample}\tSM:{sample}\tPL:ILLUMINA\tLB:{sample}-synthetic"
    yield "@PG\tID:make_example_bams\tPN:make_example_bams.py\tVN:1"
    yield "@CO\tSYNTHETIC. Generated by tools/make_example_bams.py. Not an observation."

    n_read = 0
    for seg in sorted(refs):
        cov_path = cov_dir / f"{seg}-coverage.txt"
        if not cov_path.exists():
            continue               # IRMA did not assemble it; neither do we
        ref = refs[seg]
        depth = read_coverage(cov_path)[: len(ref)]
        if len(depth) < len(ref):
            depth += [0] * (len(ref) - len(depth))
        seg_calls = calls.get(seg, {})
        target = depth_target(depth, seg_calls, len(ref))
        n_frag = int(sum(target) / (2 * READ_LEN))

        def place(start):
            ins = int(max(READ_LEN + 20, min(600, rng.gauss(INSERT_MEAN, INSERT_SD))))
            a = max(0, min(start, len(ref) - READ_LEN))
            b = max(0, min(a + ins - READ_LEN, len(ref) - READ_LEN))
            return (b, a) if b < a else (a, b)

        frags = [place(s) for s in pick_starts(target, rng, n_frag)]

        # Coordinates first, sequence second, so the depth can be CORRECTED
        # before a single base is built. Sampling starts in proportion to the
        # target does not reproduce the target: every fragment spreads its
        # depth over 150 bases, which convolves the profile and flattens
        # exactly the narrow peaks the calls raise. On the first pass 14-23
        # called sites per sample still came out under the DP their own VCF
        # record claims. Topping up against the REALISED coverage is the only
        # version that can be checked rather than tuned.
        cov = [0] * len(ref)

        def cover(f):
            for s in f:
                for p in range(s, min(len(ref), s + READ_LEN)):
                    cov[p] += 1

        for f in frags:
            cover(f)
        for pos, (_, _, _, _, _, dp) in sorted(seg_calls.items()):
            p0 = pos - 1
            # CLIP_MARGIN: coverage here counts a full READ_LEN, while the
            # aligned block is shorter wherever a tail was clipped.
            need = int(dp * CALLER_HEADROOM * CLIP_MARGIN)
            guard = 0
            while cov[p0] < need and guard < 50000:
                guard += 1
                f = place(rng.randint(max(0, p0 - READ_LEN + 12), p0))
                frags.append(f)
                cover(f)

        for s1, s2 in frags:
            first_rev = rng.random() < 0.5
            m1 = build_mate(seg, ref, s1, READ_LEN, seg_calls, rng, first_rev)
            m2 = build_mate(seg, ref, s2, READ_LEN, seg_calls, rng, not first_rev)
            if not m1 or not m2:
                continue
            n_read += 1
            qname = f"EX-{sample}-{n_read:07d}"
            mapq = 60 if rng.random() < 0.9 else rng.randint(40, 59)
            # 0x1 paired, 0x2 proper pair -- FluLens keeps only proper pairs,
            # because that is the set LoFreq's DP4 is counted over.
            for (pos, cigar, sq, ql), rev, first, mate in (
                (m1, first_rev, True, m2), (m2, not first_rev, False, m1)):
                flag = 0x1 | 0x2 | (0x40 if first else 0x80)
                flag |= 0x10 if rev else 0
                flag |= 0x20 if (not rev) else 0
                # SEQ and QUAL go out in REFERENCE orientation, reverse reads
                # included -- the 0x10 flag is what records the strand, and
                # build_mate has already put the quality decline on the right
                # end for them. Reverse-complementing here (the obvious thing,
                # and what this did first) makes every base of every reverse
                # read disagree with the reference: half the stack turns into
                # noise, and the alt fraction at every called site collapses
                # toward ~0.15 whatever the VCF says.
                tlen = (mate[0] - pos) + READ_LEN
                qual = "".join(chr(33 + q) for q in ql)
                yield (f"{qname}\t{flag}\t{seg}\t{pos + 1}\t{mapq}\t{cigar}\t=\t"
                       f"{mate[0] + 1}\t{tlen if first else -tlen}\t{sq}\t{qual}\tRG:Z:{sample}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="example_run")
    ap.add_argument("--samples", nargs="*", default=SAMPLES)
    args = ap.parse_args()

    root = Path(args.out).resolve()
    refs = read_fasta(root / "reference.fa")
    rng = random.Random(SEED)

    for sample in args.samples:
        cov_dir = root / "IRMA_results" / sample / "tables"
        vcf = root / "vcf_files" / sample / "lofreq-called-variants.vcf"
        if not cov_dir.is_dir() or not vcf.exists():
            print(f"skip {sample}: no coverage tables or no VCF", file=sys.stderr)
            continue
        calls = read_vcf(vcf)
        out_dir = root / "BAM_files" / sample
        out_dir.mkdir(parents=True, exist_ok=True)
        bam = out_dir / "final_mapped_reads.bam"

        view = subprocess.Popen(["samtools", "view", "-b", "-o", str(bam) + ".unsorted", "-"],
                                stdin=subprocess.PIPE, text=True)
        n = 0
        for rec in sam_records(sample, refs, cov_dir, calls, rng):
            view.stdin.write(rec + "\n")
            if not rec.startswith("@"):
                n += 1
        view.stdin.close()
        if view.wait() != 0:
            sys.exit(f"samtools view failed for {sample}")
        subprocess.run(["samtools", "sort", "-o", str(bam), str(bam) + ".unsorted"], check=True)
        os.remove(str(bam) + ".unsorted")
        # FluLens looks for <name>.bai before <name>.bam.bai, and the dev-server
        # path only ever asks for final_mapped_reads.bai -- so name it that.
        subprocess.run(["samtools", "index", str(bam)], check=True)
        os.replace(str(bam) + ".bai", str(out_dir / "final_mapped_reads.bai"))
        print(f"{sample}: {n:,} reads, {bam.stat().st_size / 1e6:.2f} MB")


if __name__ == "__main__":
    main()
