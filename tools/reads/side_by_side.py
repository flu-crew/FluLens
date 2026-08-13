#!/usr/bin/env python3
"""Print a few individual reads as BWA lays them down and as IRMA does.

Same read, two alignments, one window of the reference. This is the ground
truth behind every summary statistic in item 0d — and it is what shows, in one
screen, that the read matches the reference through 418 with a single
substitution at 414 and then becomes unrelated sequence at 419.

`.` marks a reference-skip (IRMA merges a read pair into one record with an N
gap between the mates), `-` a deletion.

usage: side_by_side.py <irma_bam> <bwa_bam> <ref_fa> <contig> <pos> <allele> <n_reads> <win>
"""
import re
import subprocess
import sys

CIG = re.compile(r'(\d+)([MIDNSHP=X])')


def load_ref(path, contig):
    seq, on = [], False
    with open(path) as fh:
        for line in fh:
            if line.startswith('>'):
                on = line[1:].split()[0] == contig
            elif on:
                seq.append(line.strip())
    return ''.join(seq).upper()


def render(f, lo, hi):
    """Lay a SAM record over reference window [lo, hi] -> (aligned, clipped)."""
    aligned = [' '] * (hi - lo + 1)
    clipped = [' '] * (hi - lo + 1)
    refp, qp = int(f[3]), 0
    seq = f[9].upper()
    ops = [(int(n), o) for n, o in CIG.findall(f[5])]
    for i, (n, o) in enumerate(ops):
        if o == 'S':
            first = i == 0 or (i == 1 and ops[0][1] == 'H')
            start = refp - n if first else refp
            for k in range(n):
                r = start + k
                if lo <= r <= hi:
                    clipped[r - lo] = seq[qp + k]
            qp += n
        elif o == 'H':
            pass
        elif o in 'M=X':
            for k in range(n):
                r = refp + k
                if lo <= r <= hi:
                    aligned[r - lo] = seq[qp + k]
            refp += n
            qp += n
        elif o == 'I':
            qp += n
        elif o in 'DN':
            fill = '-' if o == 'D' else '.'
            for k in range(n):
                if lo <= refp + k <= hi:
                    aligned[refp + k - lo] = fill
            refp += n
    return ''.join(aligned), ''.join(clipped)


def base_at(f, t):
    refp, qp = int(f[3]), 0
    for n, o in ((int(n), o) for n, o in CIG.findall(f[5])):
        if o == 'S':
            qp += n
        elif o == 'H':
            pass
        elif o in 'M=X':
            if refp <= t < refp + n:
                return f[9][qp + t - refp].upper()
            refp += n
            qp += n
        elif o == 'I':
            qp += n
        elif o in 'DN':
            if refp <= t < refp + n:
                return '*'
            refp += n
    return None


def main():
    irma, bwa, reffa, contig, pos, allele, nreads, win = sys.argv[1:9]
    t, n_show, w = int(pos), int(nreads), int(win)
    lo, hi = max(1, t - w), t + w
    ref = load_ref(reffa, contig)

    iout = subprocess.run(['samtools', 'view', irma, f'{contig}:{t}-{t}'],
                          capture_output=True, text=True, check=True).stdout
    picks = []
    for line in iout.splitlines():
        f = line.split('\t')
        if base_at(f, t) == allele.upper():
            picks.append(f)
        if len(picks) >= n_show:
            break

    bout = subprocess.run(
        ['samtools', 'view', bwa, f'{contig}:{max(1, lo - 200)}-{hi + 200}'],
        capture_output=True, text=True, check=True).stdout
    bwa_by_name = {}
    for line in bout.splitlines():
        f = line.split('\t')
        bwa_by_name.setdefault(f[0], []).append(f)

    print(f'\n### {contig}:{lo}-{hi}   variant at {t} ({ref[t-1]} -> {allele.upper()})')
    print(f'{"reference":>34}  {ref[lo-1:hi]}')
    print(f'{"":>34}  {" " * (t - lo)}^')

    for f in picks:
        name = f[0].split('_')[0]
        a, c = render(f, lo, hi)
        print(f'\n  read {name}')
        print(f'{"IRMA  aligned":>34}  {a}')
        if c.strip():
            print(f'{"IRMA  soft-clipped":>34}  {c}')
        for g in bwa_by_name.get(name, []):
            flag = int(g[1])
            if flag & 0x900:
                continue
            a2, c2 = render(g, lo, hi)
            mate = 'R1' if flag & 0x40 else 'R2'
            dup = ' dup' if flag & 0x400 else ''
            rev = '-' if flag & 0x10 else '+'
            if a2.strip():
                print(f'{f"BWA {mate}{rev}{dup} aligned":>34}  {a2}')
            if c2.strip():
                print(f'{f"BWA {mate}{rev}{dup} soft-clipped":>34}  {c2}')


if __name__ == '__main__':
    main()
