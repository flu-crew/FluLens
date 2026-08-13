#!/usr/bin/env python3
"""Read-level probe at a single reference position.

Walks the CIGAR of every read overlapping the position and reports, per allele:
count, strand split, base quality, mapping quality, position within the read,
and how close the read's nearest indel / soft-clip sits. Prints a summary only.

Set EXCLUDE_FLAGS=0x704 to see what a caller sees (drops unmapped, secondary,
QC-fail and duplicate records, which is samtools' own default filter). With the
default of 0 you get every record, duplicates included — useful exactly once,
to find out that duplicates are why two depth numbers disagree.

usage: EXCLUDE_FLAGS=0x704 pileup_probe.py <bam> <contig> <pos1> [<pos1> ...]
"""
import os
import subprocess
import sys
from collections import defaultdict
import re
import statistics as st

CIG = re.compile(r'(\d+)([MIDNSHP=X])')
EXCLUDE = os.environ.get('EXCLUDE_FLAGS', '0')


def base_at(read, target):
    """Return (base, bq, qpos, nearest_indel_dist, softclip, readlen, flag)."""
    flag = int(read[1])
    pos = int(read[3])          # 1-based leftmost mapped
    cigar = read[5]
    seq = read[9]
    qual = read[10]
    if cigar == '*':
        return None

    ops = [(int(n), o) for n, o in CIG.findall(cigar)]
    refp = pos
    qp = 0
    hit = None
    indel_ref_positions = []
    softclip = 0
    for n, o in ops:
        if o == 'S':
            softclip += n
            qp += n
        elif o == 'H':
            pass
        elif o in 'M=X':
            if refp <= target < refp + n:
                off = target - refp
                hit = (seq[qp + off], ord(qual[qp + off]) - 33, qp + off)
            refp += n
            qp += n
        elif o == 'I':
            indel_ref_positions.append(refp)
            qp += n
        elif o in 'DN':
            indel_ref_positions.append(refp)
            if refp <= target < refp + n:
                hit = ('*', -1, -1)       # deletion covers the position
            refp += n
    if hit is None:
        return None
    nearest = min((abs(p - target) for p in indel_ref_positions), default=None)
    return hit[0], hit[1], hit[2], nearest, softclip, len(seq), flag


def probe(bam, contig, target):
    cmd = ['samtools', 'view', '-F', EXCLUDE, bam, f'{contig}:{target}-{target}']
    out = subprocess.run(cmd, capture_output=True, text=True, check=True).stdout

    per = defaultdict(list)
    for line in out.splitlines():
        f = line.split('\t')
        r = base_at(f, target)
        if r is None:
            continue
        base, bq, qp, nearest, sc, rlen, flag = r
        mapq = int(f[4])
        rev = bool(flag & 16)
        # distance from the nearer end of the read
        end_dist = min(qp, rlen - 1 - qp) if qp >= 0 else -1
        per[base.upper()].append(
            dict(bq=bq, mapq=mapq, rev=rev, qp=qp, end_dist=end_dist,
                 nearest_indel=nearest, softclip=sc, rlen=rlen))

    total = sum(len(v) for v in per.values())
    if not total:
        print(f'\n=== {contig}:{target}   no reads ===')
        return per
    print(f'\n=== {contig}:{target}   {total} reads overlap ===')
    print(f'{"allele":>6} {"n":>6} {"%":>7} {"fwd/rev":>11} {"meanBQ":>7} '
          f'{"BQ>=30":>7} {"minBQ":>6} {"meanMQ":>7} {"MQ<30":>6} '
          f'{"med_enddist":>12} {"<=10bp end":>11} {"softclip":>9} {"indel<=20":>10}')
    for base in sorted(per, key=lambda b: -len(per[b])):
        v = per[base]
        n = len(v)
        bqs = [x['bq'] for x in v if x['bq'] >= 0]
        mqs = [x['mapq'] for x in v]
        eds = [x['end_dist'] for x in v if x['end_dist'] >= 0]
        fwd = sum(1 for x in v if not x['rev'])
        near_end = sum(1 for x in v if 0 <= x['end_dist'] <= 10)
        sc = sum(1 for x in v if x['softclip'] > 0)
        ind = sum(1 for x in v
                  if x['nearest_indel'] is not None and x['nearest_indel'] <= 20)
        print(f'{base:>6} {n:>6} {100*n/total:>6.2f}% '
              f'{fwd:>5}/{n-fwd:<5} '
              f'{(st.mean(bqs) if bqs else 0):>7.1f} '
              f'{sum(1 for b in bqs if b >= 30):>7} '
              f'{(min(bqs) if bqs else 0):>6} '
              f'{st.mean(mqs):>7.1f} '
              f'{sum(1 for m in mqs if m < 30):>6} '
              f'{(st.median(eds) if eds else 0):>12.0f} '
              f'{near_end:>11} {sc:>9} {ind:>10}')

    print('  surviving gates:')
    for q, mq in ((13, 0), (20, 0), (30, 0), (20, 30), (30, 30)):
        row = []
        kept_total = 0
        for base in sorted(per, key=lambda b: -len(per[b])):
            keep = [x for x in per[base] if x['bq'] >= q and x['mapq'] >= mq]
            kept_total += len(keep)
            if keep:
                row.append(f'{base}={len(keep)}')
        print(f'    BQ>={q:<3} MQ>={mq:<3}  depth={kept_total:<6} ' + '  '.join(row))
    return per


if __name__ == '__main__':
    bam, contig = sys.argv[1], sys.argv[2]
    for p in sys.argv[3:]:
        probe(bam, contig, int(p))
