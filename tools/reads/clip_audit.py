#!/usr/bin/env python3
"""Audit the soft-clipped tails that hide a minority allele from the callers.

For every primary, non-duplicate BWA record whose soft clip covers `target`,
lay the clipped tail back down on the reference at the coordinates it would
occupy and count how well it matches. A tail that matches everywhere except
the variant is sequence BWA declined to align; a tail that is 55-75% mismatched
at high base quality is the far side of a template switch, and the "variant" is
the last base before the junction.

Pass a name file (one read name per line, IRMA's `_3:N:0:` suffix optional) to
restrict the audit to one caller's read set. Without it you get every read at
the position, which is a different and much noisier population.

usage: clip_audit.py <bwa_bam> <ref_fa> <contig> <pos1> <allele> [name_file]
"""
import re
import subprocess
import sys
from collections import Counter
import statistics as st

CIG = re.compile(r'(\d+)([MIDNSHP=X])')


def load_ref(path, contig):
    seq, on = [], False
    with open(path) as fh:
        for line in fh:
            if line.startswith('>'):
                on = line[1:].split()[0] == contig
            elif on:
                seq.append(line.strip())
    return ''.join(seq)


def clip_spans(cigar, pos):
    """Yield (side, qstart, qlen, refstart) for each soft clip.

    refstart is the reference coordinate the first clipped base would occupy.
    """
    ops = [(int(n), o) for n, o in CIG.findall(cigar)]
    refp, qp = pos, 0
    out = []
    for i, (n, o) in enumerate(ops):
        first = i == 0 or (i == 1 and ops[0][1] == 'H')
        last = i == len(ops) - 1 or (i == len(ops) - 2 and ops[-1][1] == 'H')
        if o == 'S':
            if first:
                out.append(('left', qp, n, refp - n))
            elif last:
                out.append(('right', qp, n, refp))
            qp += n
        elif o == 'H':
            pass
        elif o in 'M=X':
            refp += n
            qp += n
        elif o == 'I':
            qp += n
        elif o in 'DN':
            refp += n
    return out


def main():
    bam, reffa, contig, pos, allele = sys.argv[1:6]
    target = int(pos)
    ref = load_ref(reffa, contig)
    keep = None
    if len(sys.argv) > 6:
        with open(sys.argv[6]) as fh:
            keep = {ln.strip().split('_')[0] for ln in fh if ln.strip()}

    # clamp: a region starting below 1 makes samtools return nothing, which
    # reads as "no clipped tails" rather than as the bad query it is
    lo_r = max(1, target - 400)
    out = subprocess.run(
        ['samtools', 'view', '-F', '0x704', bam, f'{contig}:{lo_r}-{target + 400}'],
        capture_output=True, text=True, check=True).stdout

    rows = []
    for line in out.splitlines():
        f = line.split('\t', 11)
        if keep is not None and f[0] not in keep:
            continue
        for side, qs, n, rs in clip_spans(f[5], int(f[3])):
            if not (rs <= target < rs + n):
                continue
            tail = f[9][qs:qs + n].upper()
            quals = [ord(c) - 33 for c in f[10][qs:qs + n]]
            mm, cmp_len, at_target = 0, 0, None
            for k, b in enumerate(tail):
                r = rs + k
                if not (1 <= r <= len(ref)):
                    continue
                cmp_len += 1
                if b != ref[r - 1].upper():
                    mm += 1
                if r == target:
                    at_target = b
            if at_target != allele.upper():
                continue
            rows.append(dict(side=side, cliplen=n, refstart=rs,
                             mm=mm, cmp_len=cmp_len,
                             offset=target - rs,
                             meanq=st.mean(quals) if quals else 0,
                             tlen=abs(int(f[8]))))

    if not rows:
        print('no clipped tails carry that allele')
        return

    print(f'\n=== {contig}:{target}  soft-clipped tails carrying {allele.upper()} '
          f'(primary, non-duplicate): {len(rows)} ===')
    print('  clip side          : ' +
          '  '.join(f'{k}={v}' for k, v in Counter(r["side"] for r in rows).most_common()))
    print(f'  clip length        : median {st.median([r["cliplen"] for r in rows]):.0f} '
          f'(min {min(r["cliplen"] for r in rows)}, max {max(r["cliplen"] for r in rows)})')
    print('  clip begins at ref : ' +
          '  '.join(f'{k}={v}' for k, v in Counter(r["refstart"] for r in rows).most_common(6)))
    print(f'  allele sits {int(st.median([r["offset"] for r in rows]))} bp into the clip (median)')
    print(f'  mean BQ in the clip: {st.mean([r["meanq"] for r in rows]):.1f}')
    print(f'  fragment TLEN      : median {st.median([r["tlen"] for r in rows]):.0f}')
    mmrate = [r['mm'] for r in rows]
    print('  mismatches vs reference across the whole clipped tail:')
    print(f'    median {st.median(mmrate):.0f} of {int(st.median([r["cmp_len"] for r in rows]))} bp; '
          'distribution ' +
          '  '.join(f'{k}mm={v}' for k, v in sorted(Counter(mmrate).items())[:8]))
    only_var = sum(1 for r in rows if r['mm'] == 1)
    print(f'    tails whose ONLY mismatch is the variant itself: {only_var} of {len(rows)}')


if __name__ == '__main__':
    main()
