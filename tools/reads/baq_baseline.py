#!/usr/bin/env python3
"""Compare the calls -B ADDS against the calls LoFreq already made.

A support rate means nothing on its own. If 5% of the -B-only calls are
corroborated by iVar, that is damning only if the calls LoFreq makes WITH BAQ
are corroborated far more often. Same for indel adjacency and allele frequency:
the question is never "is this number low", it is "is it lower than the calls we
already trust".

So every statistic here is computed identically over two sets from the same
samples, the same callers and the same thresholds:

  baseline   calls LoFreq makes with BAQ on -- the run's current output
  -B-only    calls that appear only once BAQ is disabled

usage: baq_baseline.py <run_dir> <vcf_dir> <sample> [sample ...]
"""
import os
import statistics
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from baq_support import read_vcf, read_ivar


def stats(name, calls, ivar_all, indels_all):
    n = len(calls)
    if not n:
        print(f'{name}: none')
        return
    sup = sum(1 for (s, k, info) in calls if k in ivar_all[s])
    dists = []
    for (s, k, info) in calls:
        near = indels_all[s].get(k[0], [])
        if near:
            dists.append(min(abs(k[1] - p) for p in near))
    near10 = sum(1 for d in dists if d <= 10)
    afs = [float(info.get('AF', 0)) for (_s, _k, info) in calls]
    dps = [int(info.get('DP', 0)) for (_s, _k, info) in calls]
    # LoFreq's SB is a phred-scaled Fisher strand-bias p: higher = more skewed.
    sbs = [int(info.get('SB', 0)) for (_s, _k, info) in calls if 'SB' in info]
    print(f'{name:<12} n={n:<6} iVar={100*sup//n:>3}%  '
          f'indel<=10bp={100*near10//n if n else 0:>3}%  '
          f'medAF={statistics.median(afs):.4f}  '
          f'medDP={int(statistics.median(dps)):>5}  '
          f'medSB={int(statistics.median(sbs)) if sbs else "-":>3}')


def main():
    run, vcfdir = sys.argv[1], sys.argv[2]
    samples = sys.argv[3:]
    base, onlyb = [], []
    ivar_all, indels_all = {}, {}
    for s in samples:
        with_baq = read_vcf(os.path.join(vcfdir, s + '.extbaq.vcf'))
        no_baq = read_vcf(os.path.join(vcfdir, s + '.nobaq.vcf'))
        iv, ind = read_ivar(os.path.join(run, 'vcf_files', s,
                                         'ivar-called-variants.tsv'))
        ivar_all[s], indels_all[s] = iv, ind
        for k, info in with_baq.items():
            base.append((s, k, info))
        for k, info in no_baq.items():
            if k not in with_baq:
                onlyb.append((s, k, info))

    print(f'{len(samples)} samples\n')
    stats('baseline', base, ivar_all, indels_all)
    stats('-B-only', onlyb, ivar_all, indels_all)

    # The same comparison restricted to calls LoFreq and iVar could BOTH have
    # seen: iVar never emits a row below its own -m 100 on the filtered depth,
    # so a LoFreq call under thin coverage cannot be corroborated even when it
    # is real, and counting those in the denominator understates both sets.
    print('\nrestricted to DP >= 200 (comfortably above iVar\'s effective floor):')
    stats('baseline', [c for c in base if int(c[2].get('DP', 0)) >= 200],
          ivar_all, indels_all)
    stats('-B-only', [c for c in onlyb if int(c[2].get('DP', 0)) >= 200],
          ivar_all, indels_all)

    # And restricted to where the two callers' frequency ranges overlap.
    print('\nrestricted to DP >= 200 and AF >= 0.02:')
    stats('baseline', [c for c in base if int(c[2].get('DP', 0)) >= 200
                       and float(c[2].get('AF', 0)) >= 0.02], ivar_all, indels_all)
    stats('-B-only', [c for c in onlyb if int(c[2].get('DP', 0)) >= 200
                      and float(c[2].get('AF', 0)) >= 0.02], ivar_all, indels_all)


if __name__ == '__main__':
    main()
