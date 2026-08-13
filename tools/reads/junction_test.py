#!/usr/bin/env python3
"""Do the six IRMA-only calls sit at read/fragment junctions?

For each call, two questions asked of the two alignments:

  IRMA  — how often does a read's aligned block END within 10 bp of the
          variant, for reads carrying the minority allele vs reads carrying
          the consensus allele? A variant that only appears in reads that
          stop right there is a junction artefact, not a variant.

  BWA   — of the reads IRMA counts as the minority allele, how many does BWA
          soft-clip over the position rather than align through it?

Run the two together. MC-522 `A_NA_N2` 104 shows no enrichment on the first
test and is unambiguous on the second.

usage: junction_test.py <run_dir>
"""
import re
import subprocess
import sys
import os

CIG = re.compile(r'(\d+)([MIDNSHP=X])')

CALLS = [
    ('MC-495', 'A_PA',     414,  'A', 'G'),
    ('MC-495', 'A_PA',     1722, 'G', 'A'),
    ('MC-522', 'A_NA_N2',  104,  'T', 'A'),
    ('MC-559', 'A_PB2',    47,   'C', 'A'),
    ('MC-587', 'A_PB1',    842,  'A', 'C'),
    ('MC-695', 'A_MP',     88,   'A', 'G'),
]


def probe(f, t):
    """Return (base, dist_to_nearer_aligned_end, in_softclip)."""
    ops = [(int(n), o) for n, o in CIG.findall(f[5])]
    refp, qp = int(f[3]), 0
    first_ref = last_ref = hit = None
    clip = False
    for i, (n, o) in enumerate(ops):
        if o == 'S':
            fst = i == 0 or (i == 1 and ops[0][1] == 'H')
            start = refp - n if fst else refp
            if start <= t < start + n:
                clip = True
                hit = f[9][qp + (t - start)].upper()
            qp += n
        elif o == 'H':
            pass
        elif o in 'M=X':
            if first_ref is None:
                first_ref = refp
            last_ref = refp + n - 1
            if refp <= t < refp + n:
                hit = f[9][qp + (t - refp)].upper()
                clip = False
            refp += n
            qp += n
        elif o == 'I':
            qp += n
        elif o == 'D':
            last_ref = refp + n - 1
            refp += n
        elif o == 'N':
            refp += n
    if first_ref is None:
        return None, None, clip
    return hit, min(abs(t - first_ref), abs(last_ref - t)), clip


def view(bam, region, flags=None):
    cmd = ['samtools', 'view'] + (['-F', flags] if flags else []) + [bam, region]
    r = subprocess.run(cmd, capture_output=True, text=True)
    return r.stdout.splitlines() if r.returncode == 0 else []


def main():
    run = sys.argv[1]
    print(f'{"sample":<8} {"locus":<9} {"pos":>6} {"chg":>5} '
          f'{"IRMA n_alt":>10} {"alt ends<=10bp":>15} {"ref ends<=10bp":>15} '
          f'{"enrich":>7} {"BWA clips alt":>14} {"BWA aligns alt":>15}')
    for sample, locus, pos, refb, alt in CALLS:
        ibam = os.path.join(run, 'IRMA_results', sample, f'{locus}.bam')
        bbam = os.path.join(run, 'BAM_files', sample, 'final_mapped_reads.bam')
        alt_d, ref_d, alt_names = [], [], set()
        for line in view(ibam, f'{locus}:{pos}-{pos}'):
            f = line.split('\t')
            b, d, _ = probe(f, pos)
            if b is None or d is None:
                continue
            if b == alt:
                alt_d.append(d)
                alt_names.add(f[0].split('_')[0])
            elif b == refb:
                ref_d.append(d)
        if not alt_d:
            print(f'{sample:<8} {locus:<9} {pos:>6}   -- no minority reads found --')
            continue
        a_end = sum(1 for d in alt_d if d <= 10) / len(alt_d)
        r_end = (sum(1 for d in ref_d if d <= 10) / len(ref_d)) if ref_d else 0

        clipped = aligned = 0
        # clamp: a region starting below 1 makes samtools error out, and the
        # empty result reads as "BWA does nothing with these reads"
        lo = max(1, pos - 300)
        for line in view(bbam, f'{locus}:{lo}-{pos + 300}', flags='0x704'):
            f = line.split('\t')
            if f[0] not in alt_names:
                continue
            b, _, clip = probe(f, pos)
            if b != alt:
                continue
            if clip:
                clipped += 1
            else:
                aligned += 1
        enrich = (a_end / r_end) if r_end else float('inf')
        print(f'{sample:<8} {locus:<9} {pos:>6} {refb+">"+alt:>5} '
              f'{len(alt_d):>10} {a_end*100:>14.1f}% {r_end*100:>14.1f}% '
              f'{enrich:>6.1f}x {clipped:>14} {aligned:>15}')


if __name__ == '__main__':
    main()
