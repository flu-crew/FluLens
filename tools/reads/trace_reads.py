#!/usr/bin/env python3
"""Trace the reads IRMA counts as the minority allele into the BWA alignment.

Takes the read names IRMA assigns the minor allele at a position, then finds
every record with those names in the BWA BAM and reports what BWA did with
them: aligned through the position (with which base), soft-clipped over it,
duplicate-flagged, placed elsewhere, or unmapped.

IRMA appends the fastq comment to the read name as `_3:N:0:INDEX` (the 3 is its
marker for a merged pair), and BWA keeps the bare Illumina name, which never
contains an underscore. Comparing them unnormalised finds 0 of 125 and looks
like a real result.

usage: trace_reads.py <irma_bam> <bwa_bam> <contig> <pos1> <minor_allele>
"""
import re
import subprocess
import sys
from collections import Counter

CIG = re.compile(r'(\d+)([MIDNSHP=X])')


def walk(cigar, pos, target):
    """Where does `target` fall in this record?

    Returns (kind, base_index) where kind is one of
    'aligned' | 'deleted' | 'skipped' | 'clipped' | 'outside'.
    base_index indexes SEQ for 'aligned' and 'clipped'.
    """
    ops = [(int(n), o) for n, o in CIG.findall(cigar)]
    refp, qp = pos, 0
    for i, (n, o) in enumerate(ops):
        if o == 'S':
            lo = refp - n if i == 0 else refp
            hi = refp if i == 0 else refp + n
            if lo <= target < hi:
                return 'clipped', qp + (target - lo)
            qp += n
        elif o == 'H':
            pass
        elif o in 'M=X':
            if refp <= target < refp + n:
                return 'aligned', qp + (target - refp)
            refp += n
            qp += n
        elif o == 'I':
            qp += n
        elif o == 'D':
            if refp <= target < refp + n:
                return 'deleted', -1
            refp += n
        elif o == 'N':
            if refp <= target < refp + n:
                return 'skipped', -1
            refp += n
    return 'outside', -1


def irma_minor_names(bam, contig, target, allele):
    out = subprocess.run(['samtools', 'view', bam, f'{contig}:{target}-{target}'],
                         capture_output=True, text=True, check=True).stdout
    names = set()
    for line in out.splitlines():
        f = line.split('\t')
        kind, idx = walk(f[5], int(f[3]), target)
        if kind == 'aligned' and f[9][idx].upper() == allele:
            names.add(f[0].split('_')[0])
    return names


def main():
    irma_bam, bwa_bam, contig, pos, allele = sys.argv[1:6]
    target = int(pos)
    names = irma_minor_names(irma_bam, contig, target, allele.upper())
    print(f'{contig}:{target}  IRMA reads carrying {allele}: {len(names)}')

    # full scan: the mate that does not reach the position, and any record
    # placed on another contig or left unmapped, are part of the answer
    out = subprocess.run(['samtools', 'view', bwa_bam],
                         capture_output=True, text=True, check=True).stdout

    verdict = Counter()
    clipped_base = Counter()
    aligned_base = Counter()
    seen = set()
    for line in out.splitlines():
        f = line.split('\t', 11)
        if f[0] not in names:
            continue
        seen.add(f[0])
        flag = int(f[1])
        tags = []
        if flag & 0x400:
            tags.append('dup')
        if flag & 0x100 or flag & 0x800:
            tags.append('sec/supp')
        if flag & 0x4:
            verdict['unmapped'] += 1
            continue
        if f[2] != contig:
            verdict[f'mapped to {f[2]}'] += 1
            continue
        kind, idx = walk(f[5], int(f[3]), target)
        base = f[9][idx].upper() if idx >= 0 else '-'
        label = kind + (' [' + ','.join(tags) + ']' if tags else '')
        verdict[label] += 1
        if kind == 'clipped':
            clipped_base[base] += 1
        elif kind == 'aligned':
            aligned_base[base] += 1

    print(f'  found in BWA BAM: {len(seen)} of {len(names)} names '
          f'({len(names) - len(seen)} absent entirely)')
    print('  what BWA did with their records:')
    for k, v in verdict.most_common():
        print(f'    {k:<28} {v}')
    if aligned_base:
        print('  base BWA aligned at the position: ' +
              '  '.join(f'{b}={n}' for b, n in aligned_base.most_common()))
    if clipped_base:
        print('  base sitting in the SOFT-CLIPPED tail: ' +
              '  '.join(f'{b}={n}' for b, n in clipped_base.most_common()))


if __name__ == '__main__':
    main()
