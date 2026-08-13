#!/usr/bin/env python3
"""List the read names IRMA assigns a given allele at a position.

Feeds `clip_audit.py`'s optional name file. Names come out normalised — IRMA's
`_3:N:0:INDEX` comment stripped — so they match the BWA BAM directly.

usage: irma_allele_names.py <irma_bam> <contig> <pos1> <allele> > names.txt
"""
import re
import subprocess
import sys

CIG = re.compile(r'(\d+)([MIDNSHP=X])')


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
    bam, contig, pos, allele = sys.argv[1], sys.argv[2], int(sys.argv[3]), sys.argv[4]
    out = subprocess.run(['samtools', 'view', bam, f'{contig}:{pos}-{pos}'],
                         capture_output=True, text=True, check=True).stdout
    for line in out.splitlines():
        f = line.split('\t')
        if base_at(f, pos) == allele.upper():
            print(f[0].split('_')[0])


if __name__ == '__main__':
    main()
