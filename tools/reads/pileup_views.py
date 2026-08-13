#!/usr/bin/env python3
"""The same position as three different callers are shown it.

  no BAQ, no overlap removal   — what `samtools depth` and `depth_profiles/`
                                 count, and what item 0d called "BWA depth"
  no BAQ, overlap removal ON   — what iVar is handed, and what its -m floor
                                 and TOTAL_DP are actually measured on
  BAQ ON, overlap removal ON   — roughly what LoFreq works from

The three differ by 25-45% on this library. Verified against iVar's own output:
`TOTAL_DP` equals the middle column exactly.

usage: pileup_views.py <bam> <ref_fa> <contig:start-end> <min_bq>
"""
import subprocess
import sys


def run(bam, ref, region, extra, label, q):
    cmd = (['samtools', 'mpileup', '-aa', '-A', '-d', '0', '-Q', '0',
            '--reference', ref, '-r', region] + extra + [bam])
    o = subprocess.run(cmd, capture_output=True, text=True).stdout.strip()
    if not o:
        print(f'  {label:<36} (no output)')
        return
    f = o.split('\t')
    quals = [ord(c) - 33 for c in f[5]]
    bases = f[4]
    i = k = keep = 0
    alt = {}
    while i < len(bases):
        c = bases[i]
        if c == '^':
            i += 2
            continue
        if c == '$':
            i += 1
            continue
        if c in '+-':
            j = i + 1
            num = ''
            while j < len(bases) and bases[j].isdigit():
                num += bases[j]
                j += 1
            i = j + int(num or 0)
            continue
        if c == '*':
            k += 1
            i += 1
            continue
        if k < len(quals) and quals[k] >= q:
            keep += 1
            u = c.upper()
            if u not in '.,':
                alt[u] = alt.get(u, 0) + 1
        k += 1
        i += 1
    shown = ' '.join(f'{a}={n}' for a, n in sorted(alt.items())) or '-'
    print(f'  {label:<36} mpileup_depth={f[3]:>5}  passing Q{q}={keep:<5} alt={shown}')


def main():
    bam, ref, region, q = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4])
    run(bam, ref, region, ['-B', '-x'], 'no BAQ, no overlap removal', q)
    run(bam, ref, region, ['-B'],       'no BAQ, overlap removal ON (iVar)', q)
    run(bam, ref, region, [],           'BAQ ON, overlap removal ON (LoFreq)', q)


if __name__ == '__main__':
    main()
