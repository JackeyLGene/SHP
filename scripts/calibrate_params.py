"""SHP parameter calibration: sweep D x percentile for optimal CDS/cDNA separation.

Uses ~300 matched human CDS-cDNA pairs to find optimal hash dimension (D)
and fair-IID percentile (P) that maximize CDS/cDNA structural discrimination.
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

import gzip, statistics, math, itertools
from collections import defaultdict
from shp.core import compute_shp

CDS_PATH = os.environ.get('SHP_CDS', 'g:/GEME/EE/revolution/ensembl_release_115/cds/Homo_sapiens.GRCh38.cds.all.fa.gz')
CDNA_PATH = os.environ.get('SHP_CDNA', 'g:/GEME/EE/revolution/ensembl_release_115/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz')
W = 128; NGRAM = 3; ALPHABET = 'ACGT'

def load_genes(path):
    """Return {gene_symbol: clean_sequence} for protein_coding genes."""
    genes = {}
    with gzip.open(path, 'rt', encoding='utf-8', errors='ignore') as f:
        cur_gene, cur_bio, parts = None, None, None
        for line in f:
            if line.startswith('>'):
                if cur_gene and parts and cur_bio == 'protein_coding':
                    seq = ''.join(parts).upper()
                    clean = ''.join(c for c in seq if c in 'ACGT')
                    if len(clean) >= 300:
                        key = cur_gene.upper()
                        if key and (key not in genes or len(clean) > len(genes.get(key, ''))):
                            genes[key] = clean
                import re
                m_g = re.search(r'gene_symbol:(\S+)', line.strip())
                m_b = re.search(r'gene_biotype:(\S+)', line.strip())
                cur_gene = m_g.group(1) if m_g else None
                cur_bio = m_b.group(1) if m_b else None; parts = []
            else:
                if parts is not None: parts.append(line.strip())
    return genes

print('Loading CDS and cDNA...')
cds = load_genes(CDS_PATH); cdna = load_genes(CDNA_PATH)
common = sorted(set(cds) & set(cdna))
print(f'CDS: {len(cds)} genes, cDNA: {len(cdna)} genes, shared: {len(common)}')

# Sample 300 for speed
import random; random.seed(42)
sample = random.sample(common, min(300, len(common)))
print(f'Calibration sample: {len(sample)} genes')

D_VALUES = [16, 32, 64, 128, 256, 512]
P_VALUES = [90, 95, 97, 98, 99, 99.5]

# Pre-compute fair IID theta0 for each (D, P) pair
print('\nCalibrating theta0 for each (D, P)...')
theta0_map = {}
for D in D_VALUES:
    for P in P_VALUES:
        # Generate fair IID streams and compute percentile threshold
        ds_all = []
        for seed in range(10):
            rng = random.Random(seed)
            seq = ''.join(rng.choice('ACGT') for _ in range(10000))
            result = compute_shp(seq, ngram=NGRAM, dim=D, window=W, theta0=0.0)
            # Extract displacement distribution from compute_shp internals
            # Re-run with zero theta to get raw displacements
            clean = ''.join(c for c in seq.upper() if c in 'ACGT')
            stride = max(1, W // 5)
            from shp.core import _windows, _hash_bucket, _jaccard_distance
            hs = []
            for win in _windows(clean, W, stride):
                kmers = [win[i:i+NGRAM] for i in range(len(win)-NGRAM+1)]
                chroma = {_hash_bucket('C:'+k, D) for k in kmers}
                rhythm = {_hash_bucket('R:'+a+'>'+b, D) for a,b in zip(kmers[:-1], kmers[1:])}
                hs.append(_jaccard_distance(chroma, rhythm))
            ds = sorted([abs(hs[i]-hs[i-1]) for i in range(1, len(hs))])
            idx = max(0, int(P/100 * len(ds)) - 1)
            ds_all.append(ds[idx])
        theta0_map[(D,P)] = statistics.mean(ds_all)
        print(f'  D={D:>4} P={P:>5} -> theta0={theta0_map[(D,P)]:.6f}')

# Sweep
print(f'\n{"="*80}')
print(f'PARAMETER SWEEP: D x Percentile on {len(sample)} CDS-cDNA pairs')
print(f'{"="*80}')
print(f'{"D":>5} {"P":>5} {"theta0":>10} {"CDS_fw":>8} {"cDNA_fw":>8} {"d":>8} {"CDS>cDNA":>7}')
print(f'{"-"*60}')

best_d = -1; best_score = -999; best_D = 64; best_P = 99

for D in D_VALUES:
    for P in P_VALUES:
        t0 = theta0_map[(D,P)]
        cds_fws = []; cdna_fws = []; cdna_gt = 0
        for gene in sample:
            r_cds = compute_shp(cds[gene], ngram=NGRAM, dim=D, window=W, theta0=t0)
            r_cdna = compute_shp(cdna[gene], ngram=NGRAM, dim=D, window=W, theta0=t0)
            cds_fws.append(r_cds.fixed_wit)
            cdna_fws.append(r_cdna.fixed_wit)
            if r_cdna.fixed_wit > r_cds.fixed_wit: cdna_gt += 1

        # Cohen's d: CDS > cDNA is the correct direction
        m_cds = statistics.mean(cds_fws); m_cdna = statistics.mean(cdna_fws)
        pooled = math.sqrt((statistics.stdev(cds_fws)**2 + statistics.stdev(cdna_fws)**2) / 2)
        d = (m_cds - m_cdna) / max(0.0001, pooled)
        n_pos = sum(1 for i in range(len(sample)) if cds_fws[i] > cdna_fws[i])

        print(f'{D:>5} {P:>5} {t0:>10.4f} {m_cds:>8.4f} {m_cdna:>8.4f} {d:>+8.2f} {n_pos:>5}/{len(sample)}')

        if d > best_score:
            best_score = d; best_D = D; best_P = P

print(f'\nOptimal: D={best_D}, P={best_P}, d={best_score:+.2f}')

# Show theta0 at optimal
print(f'theta0(optimal) = {theta0_map[(best_D,best_P)]:.6f}')
print(f'Current default: D=64, P=99, theta0=0.0999')
