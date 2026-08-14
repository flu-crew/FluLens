#!/usr/bin/env python3
"""Are LoFreq's -B-only calls real, or is BAQ suppressing them for a reason?

For every call that appears without BAQ and not with it, ask three independent
questions:

  iVar        does the run's own iVar TSV call the same position+ALT? iVar
              already runs with -B, so agreement is partly shared method, but
              its statistics and its PASS test are its own.
  GATK4       does HaplotypeCaller call it? Reported but NOT evidence either
              way at these frequencies: GATK4 is a genotype caller and does not
              call 3% variants at all, so a 0% agreement rate here means
              "wrong instrument", not "unsupported".
  indel       distance to the nearest indel, taken from iVar's TSV and NOT from
              the GATK4 indel VCF. That distinction decides the answer: GATK4
              reports 4 indel records across the whole run, iVar reports 916
              across 129 of 143 samples, because only one of them detects
              indels below genotype frequency. BAQ exists to suppress spurious
              SNPs beside indels, so a -B-only call sitting on one is exactly
              the failure mode, and the GATK4 file would have said there was
              nothing to sit on.

usage: baq_support.py <run_dir> <vcf_dir> <sample> [sample ...]
"""
import os
import sys


def read_vcf(path):
    out = {}
    if not os.path.exists(path):
        return out
    with open(path) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            f = line.rstrip('\n').split('\t')
            if len(f) < 8:
                continue
            info = dict(kv.split('=', 1) for kv in f[7].split(';') if '=' in kv)
            out[(f[0], int(f[1]), f[3], f[4])] = info
    return out


def read_ivar(path):
    """Returns (substitution calls, indel positions by segment)."""
    out, indels = {}, {}
    if not os.path.exists(path):
        return out, indels
    with open(path) as fh:
        hdr = fh.readline().rstrip('\n').split('\t')
        ix = {n: i for i, n in enumerate(hdr)}
        for line in fh:
            f = line.rstrip('\n').split('\t')
            if len(f) < len(hdr):
                continue
            region, pos, alt = f[ix['REGION']], int(f[ix['POS']]), f[ix['ALT']]
            if alt[:1] in '+-':
                indels.setdefault(region, []).append(pos)
            else:
                out[(region, pos, f[ix['REF']], alt)] = (
                    f[ix['ALT_FREQ']], f[ix['TOTAL_DP']], f[ix['PASS']])
    return out, indels


def main():
    run, vcfdir = sys.argv[1], sys.argv[2]
    samples = sys.argv[3:]
    rows = []
    for s in samples:
        base = os.path.join(vcfdir, s)
        with_baq = read_vcf(base + '.extbaq.vcf')
        no_baq = read_vcf(base + '.nobaq.vcf')
        only_b = sorted(set(no_baq) - set(with_baq), key=lambda k: (k[0], k[1]))
        lost = sorted(set(with_baq) - set(no_baq), key=lambda k: (k[0], k[1]))

        rd = os.path.join(run, 'vcf_files', s)
        ivar, indels = read_ivar(os.path.join(rd, 'ivar-called-variants.tsv'))
        gatk = read_vcf(os.path.join(rd, 'gatk4-filtered-snps.vcf'))

        for key in only_b:
            chrom, pos, ref, alt = key
            near = indels.get(chrom, [])
            d = min((abs(pos - p) for p in near), default=None)
            iv = ivar.get(key)
            rows.append((s, chrom, pos, ref, alt,
                         no_baq[key].get('AF', '?'), no_baq[key].get('DP', '?'),
                         'yes' if iv else 'no',
                         'yes' if key in gatk else 'no',
                         'n/a' if d is None else str(d)))
        if lost:
            for key in lost:
                rows.append((s, key[0], key[1], key[2], key[3], 'LOST-BY-B', '', '', '', ''))

    lost = [r for r in rows if r[5] == 'LOST-BY-B']
    rows = [r for r in rows if r[5] != 'LOST-BY-B']
    n = len(rows)
    if not n:
        print('no -B-only calls')
        return

    # Where they sit. The claim -B is meant to fix is that BAQ blinds LoFreq in
    # the first ~100 bases of every segment, so if that is the whole story the
    # recovered calls concentrate there.
    bands = [(0, 49), (50, 99), (100, 149), (150, 199), (200, 10**9)]
    counts = []
    for lo, hi in bands:
        counts.append(sum(1 for r in rows if lo <= r[2] <= hi))
    print('-B-only calls by distance from segment start:')
    print('     0-49   50-99  100-149  150-199    200+')
    print('  ' + ''.join(f'{c:>8}' for c in counts))

    sup_iv = sum(1 for r in rows if r[7] == 'yes')
    d = [int(r[9]) for r in rows if r[9] not in ('n/a', '')]
    near10 = sum(1 for x in d if x <= 10)
    near25 = sum(1 for x in d if x <= 25)
    print(f'\n-B-only calls: {n}   (calls LOST by -B: {len(lost)})')
    print(f'  also called by iVar      : {sup_iv} ({100*sup_iv//n}%)')
    print(f'  within 10 bp of an indel : {near10} ({100*near10//n}%)')
    print(f'  within 25 bp of an indel : {near25} ({100*near25//n}%)')
    if d:
        print(f'  median distance to nearest indel: {sorted(d)[len(d)//2]}')

    risky = sorted((r for r in rows if r[9] not in ('n/a', '') and int(r[9]) <= 10),
                   key=lambda r: int(r[9]))
    if risky:
        print('\n  indel-adjacent (<=10 bp) -- the calls BAQ exists to suppress:')
        print(f'  {"sample":<8} {"segment":<9} {"pos":>6} {"var":<6} {"AF":>8} '
              f'{"DP":>5} {"iVar":>5} {"d":>4}')
        for r in risky:
            print(f'  {r[0]:<8} {r[1]:<9} {r[2]:>6} {r[3]+">"+r[4]:<6} {r[5]:>8} '
                  f'{r[6]:>5} {r[7]:>5} {r[9]:>4}')


if __name__ == '__main__':
    main()
