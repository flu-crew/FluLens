#!/usr/bin/env python3
"""Size the band between the depth that gets reported and the depth iVar sees.

`depth_profiles/` and item 0d's "BWA depth" are raw `samtools depth` counts.
iVar is handed a pileup that has already had overlapping mate bases zeroed and
then applies `-q`, so the depth its `-m` floor tests is smaller. Every position
where the raw depth clears 100 but iVar's does not is a position the run calls
adequately covered and no caller ever evaluates.

On the 2026-08-09 swine WGS run: 111,873 positions across 130 of 143 samples.

Writes a per-position table for the band and prints per-sample counts.

usage: depth_band.py <run_dir> <out_tsv> [sample ...]
"""
import os
import subprocess
import sys

MIN_DEPTH = 100
MIN_Q = 30


def visible(bases, quals, minq):
    """Count bases in a pileup column that survive a -q threshold."""
    i = k = keep = 0
    n = len(bases)
    while i < n:
        c = bases[i]
        if c == '^':          # read start, next char is the mapping quality
            i += 2
            continue
        if c == '$':          # read end
            i += 1
            continue
        if c in '+-':         # indel, followed by a length and that many bases
            j = i + 1
            num = ''
            while j < n and bases[j].isdigit():
                num += bases[j]
                j += 1
            i = j + int(num or 0)
            continue
        if c == '*':          # deletion placeholder, consumes a quality
            k += 1
            i += 1
            continue
        if k < len(quals) and (ord(quals[k]) - 33) >= minq:
            keep += 1
        k += 1
        i += 1
    return keep


def sample_scan(bam, ref, emit):
    """Walk the overlap-removed and raw pileups in step; emit band positions."""
    ov = ['samtools', 'mpileup', '-aa', '-A', '-B', '-d', '0', '-Q', '0',
          '--reference', ref, bam]
    raw = ['samtools', 'mpileup', '-aa', '-A', '-B', '-x', '-d', '0', '-Q', '0',
           '--reference', ref, bam]
    with subprocess.Popen(ov, stdout=subprocess.PIPE, text=True,
                          stderr=subprocess.DEVNULL) as p_ov, \
         subprocess.Popen(raw, stdout=subprocess.PIPE, text=True,
                          stderr=subprocess.DEVNULL) as p_raw:
        band = 0
        for l_ov, l_raw in zip(p_ov.stdout, p_raw.stdout):
            f = l_ov.rstrip('\n').split('\t')
            g = l_raw.rstrip('\n').split('\t')
            if len(f) < 6 or len(g) < 6:
                continue
            raw_depth = int(g[3])
            if raw_depth < MIN_DEPTH:
                continue
            vis = visible(f[4], f[5], MIN_Q)
            if vis < MIN_DEPTH:
                band += 1
                emit(f'{f[0]}\t{f[1]}\t{raw_depth}\t{vis}\n')
        return band


def main():
    run, out = sys.argv[1], sys.argv[2]
    ref = os.path.join(run, 'reference.fa')
    bamdir = os.path.join(run, 'BAM_files')
    samples = sys.argv[3:] or sorted(os.listdir(bamdir))

    total_band = 0
    per = []
    with open(out, 'w') as fh:
        fh.write('sample\tcontig\tpos\traw_depth\tivar_visible_depth\n')
        for s in samples:
            bam = os.path.join(bamdir, s, 'final_mapped_reads.bam')
            if not os.path.exists(bam):
                continue
            n = sample_scan(bam, ref, lambda line, s=s: fh.write(f'{s}\t{line}'))
            per.append((s, n))
            total_band += n

    nonzero = sorted(((s, n) for s, n in per if n), key=lambda x: -x[1])
    print(f'samples scanned              : {len(per)}')
    print(f'samples with a band          : {len(nonzero)}')
    print(f'positions raw>=100, iVar<100 : {total_band}')
    if nonzero:
        print('worst samples: ' + '  '.join(f'{s}={n}' for s, n in nonzero[:10]))
    print(f'table written to {out}')


if __name__ == '__main__':
    main()
