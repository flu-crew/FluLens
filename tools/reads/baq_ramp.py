#!/usr/bin/env python3
"""How much usable depth BAQ costs, as a function of distance from a contig edge.

BAQ's HMM has no flanking reference to work with at a contig edge, so it
collapses the qualities of reads hanging off it. On the swine WGS run this
removes about half the Q30-passing bases in the first 50 positions of every
segment, a quarter in the next 50, and nothing past 150 — which is exactly the
shape of LoFreq's call distribution (0 calls in 0-49, 6 in 50-99, 87 in
100-149, 527 in 150-199, pooled over 143 samples).

usage: baq_ramp.py <bam> <ref_fa> <contig> <end_pos>
"""
import subprocess
import sys
import statistics as st


def profile(bam, ref, contig, end, extra):
    cmd = (['samtools', 'mpileup', '-aa', '-A', '-d', '0', '-Q', '0',
            '--reference', ref, '-r', f'{contig}:1-{end}'] + extra + [bam])
    d = {}
    for line in subprocess.run(cmd, capture_output=True, text=True).stdout.splitlines():
        f = line.split('\t')
        if len(f) < 6:
            continue
        q = [ord(c) - 33 for c in f[5]]
        d[int(f[1])] = (len(q), sum(1 for x in q if x >= 30))
    return d


def main():
    bam, ref, contig, end = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4])
    off = profile(bam, ref, contig, end, ['-B'])
    on = profile(bam, ref, contig, end, [])
    print(f'  {"window":<12}{"depth":>8}{"BQ>=30 no BAQ":>16}{"BQ>=30 with BAQ":>18}{"lost":>8}')
    for a, b in ((1, 50), (51, 100), (101, 150), (151, 200), (201, 300), (301, 500)):
        ps = [p for p in off if a <= p <= b]
        if not ps:
            continue
        dep = st.mean(off[p][0] for p in ps)
        o = st.mean(off[p][1] for p in ps)
        n = st.mean(on[p][1] for p in ps)
        print(f'  {f"{a}-{b}":<12}{dep:>8.0f}{o:>16.0f}{n:>18.0f}'
              f'{(100 * (1 - n / o) if o else 0):>7.0f}%')


if __name__ == '__main__':
    main()
