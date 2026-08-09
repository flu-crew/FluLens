#!/usr/bin/env python3
"""Does FluLens' DERIVED consensus agree with IRMA's actual consensus?

FluLens builds its consensus by painting calls above 50% onto the reference --
a re-derivation. IRMA already produced a consensus, and it is the one FluMut
screens. If the two disagree, the app is showing a different consensus from the
one the marker calls were made against.

Derived  : all_sample_amino_acids.txt, af_type=fraction, freq>0.5, max per cell
           (mirrors buildConsensus, including its exclusion of genotype-only AF)
IRMA     : IRMA-consensus-contigs/<sample>.fasta, spliced through the reference
           CDS intervals and translated, compared residue-by-residue to the
           reference protein

Only segments whose contig length equals the reference is compared; a truncated
contig does not share the reference's coordinate frame and is counted separately
rather than guessed at.
"""
import csv, glob, os, re, sys
from collections import defaultdict, Counter

RUN = sys.argv[1]

_B, _AA = 'TCAG', 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
CODONS = {i+j+k: _AA[n] for n, (i, j, k) in
          enumerate((i, j, k) for i in _B for j in _B for k in _B)}
translate = lambda s: ''.join(CODONS.get(s[i:i+3].upper(), 'X') for i in range(0, len(s)-2, 3))
norm_gene = lambda g: re.sub(r'_[HN]?\d+$', '', str(g))


def read_fasta(p):
    d, nm, buf = {}, None, []
    for L in open(p):
        L = L.strip()
        if not L:
            continue
        if L[0] == '>':
            if nm: d[nm] = ''.join(buf)
            nm, buf = L[1:].split()[0], []
        else:
            buf.append(L)
    if nm: d[nm] = ''.join(buf)
    return d


# ---- reference CDS per product (same rule as makeGTF.R / the viewer) ----
cds, seen = {}, set()
gdir = os.path.join(RUN, 'reference_gtf')
for fn in sorted(f for f in os.listdir(gdir) if f.endswith('.gtf')):
    if fn.endswith('combined.gtf'):
        continue
    for line in open(os.path.join(gdir, fn)):
        c = line.rstrip('\n').split('\t')
        if len(c) < 9 or c[2] != 'CDS':
            continue
        m = re.search(r'gene_id "([^"]+)"', c[8])
        if not m:
            continue
        g = norm_gene(m.group(1))
        k = (g, c[0], c[3], c[4])
        if k in seen:
            continue
        seen.add(k)
        cds.setdefault(g, {'seq': c[0], 'ex': []})['ex'].append((int(c[3]), int(c[4])))
for r in cds.values():
    r['ex'].sort()

ref = read_fasta(os.path.join(RUN, 'reference.fa'))
refprot = {g: translate(''.join(ref[r['seq']][s-1:e] for s, e in r['ex']))
           for g, r in cds.items() if r['seq'] in ref}

# ---- derived consensus, mirroring buildConsensus ----
derived, best = {}, {}
aa_path = os.path.join(RUN, 'variant_analysis', 'all_sample_amino_acids.txt')
with open(aa_path) as fh:
    for row in csv.DictReader(fh, delimiter='\t'):
        if row['af_type'] != 'fraction':          # genotype-only AF is not a fraction
            continue
        try:
            f = float(row['allele_frequency'])
        except (TypeError, ValueError):
            continue
        if f <= 0.5:
            continue
        alt = (row.get('alternative_aa') or '').strip()
        if not alt or len(alt) != 1:
            continue
        try:
            pos = int(row['aa_position'])
        except (TypeError, ValueError):
            continue
        k = (row['sample'], row['product'], pos)
        if f > best.get(k, 0):
            best[k] = f
            derived[k] = alt

# ---- IRMA consensus ----
irma = {}
skipped = Counter()
compared_segments = 0
for path in sorted(glob.glob(os.path.join(RUN, 'IRMA-consensus-contigs', '*.fasta'))):
    sample = os.path.basename(path)[:-6]
    contigs = read_fasta(path)
    for g, r in cds.items():
        locus, rseq = r['seq'], ref.get(r['seq'])
        if rseq is None:
            continue
        con = contigs.get(locus)
        if con is None:
            skipped['segment absent from IRMA'] += 1
            continue
        if len(con) != len(rseq):
            skipped['contig truncated (no shared frame)'] += 1
            continue
        compared_segments += 1
        prot = translate(''.join(con[s-1:e] for s, e in r['ex']))
        rp = refprot.get(g, '')
        for i, a in enumerate(prot):
            if i >= len(rp):
                break
            if a != rp[i] and a != 'X':           # X = N in the contig, not a call
                irma[(sample, g, i+1)] = a

# ---- compare ----
dk, ik = set(derived), set(irma)
agree = {k for k in dk & ik if derived[k] == irma[k]}
resid = {k for k in dk & ik if derived[k] != irma[k]}
only_d, only_i = dk - ik, ik - dk

print(f'segments compared: {compared_segments}   skipped: {dict(skipped)}')
print()
print(f'derived consensus residues : {len(dk)}')
print(f'IRMA consensus residues    : {len(ik)}')
print()
print(f'  agree (same cell, same residue) : {len(agree)}')
print(f'  same cell, DIFFERENT residue    : {len(resid)}')
print(f'  in derived only                 : {len(only_d)}')
print(f'  in IRMA only                    : {len(only_i)}')

def show(label, keys, n=8):
    if not keys:
        return
    print(f'\n{label} (first {min(n,len(keys))} of {len(keys)}):')
    for k in sorted(keys)[:n]:
        s, g, p = k
        rp = refprot.get(g, '')
        rr = rp[p-1] if p-1 < len(rp) else '?'
        print(f'  {s:<8} {g:<7} codon {p:<4} ref {rr}  derived {derived.get(k,"-")}  IRMA {irma.get(k,"-")}')

show('SAME CELL, DIFFERENT RESIDUE', resid)
show('IN DERIVED ONLY (app paints a change IRMA does not have)', only_d)
show('IN IRMA ONLY (IRMA has a change the app never paints)', only_i)

by_prod = Counter(g for _, g, _ in only_i)
if by_prod:
    print('\nIRMA-only, by product:', dict(by_prod.most_common()))
