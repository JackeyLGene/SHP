"""SHP + LRU on DNA: replace raw cross-harm with cache-based structural tracking.

For each sliding window:
  chroma = set of 3-mer hash bins (which patterns are active)
  rhythm = adjacent 3-mer transition fingerprint (beta-encoded)
  LRU cache: stores (fingerprint -> expected_next, hit_count)
  L1 cross-harm = 1.0 (tension: chroma active but no rhythm match)
                 = 0.0 (harmony: chroma and rhythm aligned via cache hit)
                 = 0.5 (neutral: chroma not active in cache)
  fw_LRU = fraction of transitions where cross-harm is 1.0

Compares to raw SHP fw on the same genes.
"""
import gzip, hashlib, math, random, re, os, sys, csv
from collections import defaultdict, Counter

random.seed(42)

# ── LRU Core (from unified_lru_proof.py) ──
def extended_gcd(a, b):
    if b == 0: return a, 1, 0
    g, x, y = extended_gcd(b, a % b)
    return g, y, x - (a // b) * y

def beta_encode(seq):
    k = len(seq)
    if k == 0: return 0, 1
    n = max(1, math.factorial(k))
    moduli = [n * (i + 1) + 1 for i in range(k)]
    M = 1
    for mod in moduli: M *= mod
    m = 0
    for i, a in enumerate(seq):
        Mi = M // moduli[i]
        _, inv, _ = extended_gcd(Mi, moduli[i])
        m = (m + a * Mi * inv) % M
    return m, n

class DNASHPLRU:
    """DNA SHP + LRU: chroma/rhythm cache with self-referential lookup."""

    def __init__(self, C=8, pattern_len=2, dim=128):
        self.C = C; self.k = pattern_len; self.dim = dim
        self.cache = []
        self.chroma_active = set()
        self.rhythm_history = []
        self.stats = {'total': 0, 'tension': 0, 'harmony': 0, 'neutral': 0,
                      'cache_hits': 0, 'correct_predictions': 0}

    def _hash(self, text):
        return int(hashlib.md5(text.encode()).hexdigest()[:8], 16) % self.dim

    def window_to_pattern(self, win):
        """Extract chroma and rhythm from a DNA window."""
        NG = 3
        kmers = [win[i:i+NG] for i in range(len(win)-NG+1)]
        chroma = {self._hash('C:'+k) for k in kmers}
        rhythm = tuple(self._hash('R:'+kmers[i]+'>'+kmers[i+1])
                       for i in range(len(kmers)-1))
        return chroma, rhythm

    def step(self, win):
        chroma, rhythm = self.window_to_pattern(win)
        self.stats['total'] += 1

        # ── L1: orthogonal projection ──
        chroma_visible = len(chroma & self.chroma_active) > 0 if self.chroma_active else False

        # ── L2-L4: LRU self-lookup ──
        transition_fp = tuple(rhythm[:self.k]) if len(rhythm) >= self.k else None
        cache_hit = False; rhythm_match = False

        if transition_fp and len(transition_fp) == self.k:
            fp_m, fp_n = beta_encode(list(transition_fp))
            fp = (fp_m % 10007, fp_n)

            for i, (cached_fp, _, count) in enumerate(self.cache):
                if cached_fp == fp:
                    cache_hit = True; rhythm_match = True
                    new_count = count + 1
                    entry = self.cache.pop(i)
                    self.cache.insert(0, (fp, entry[1], new_count))
                    break

            if not cache_hit:
                if len(self.cache) >= self.C:
                    self.cache.pop()
                self.cache.insert(0, (fp, None, 1))

        # Update chroma set (LRU eviction)
        for h in chroma:
            self.chroma_active.discard(h)  # refresh
            self.chroma_active.add(h)
        while len(self.chroma_active) > self.C:
            # LRU evict
            if self.rhythm_history:
                old = self.rhythm_history[0]
                if old in self.chroma_active:
                    self.chroma_active.discard(old)

        self.rhythm_history.append(tuple(rhythm[:self.k]))
        if len(self.rhythm_history) > self.C * 10:
            self.rhythm_history = self.rhythm_history[-self.C*5:]

        # L1 cross-harm
        if chroma_visible and not rhythm_match:
            self.stats['tension'] += 1
            crossharm = 1.0
        elif chroma_visible and rhythm_match:
            self.stats['harmony'] += 1
            crossharm = 0.0
        else:
            self.stats['neutral'] += 1
            crossharm = 0.5

        if cache_hit: self.stats['cache_hits'] += 1
        return crossharm

# ── Compute LRU fw for a gene ──
def compute_lru_fw(seq, C=8, W=128, stride=None):
    if stride is None: stride = max(1, W // 5)
    clean = ''.join(c for c in seq.upper() if c in 'ACGT')
    if len(clean) < W: return None

    lru = DNASHPLRU(C=C)
    ch_vals = []
    for start in range(0, len(clean) - W + 1, stride):
        win = clean[start:start+W]
        ch = lru.step(win)
        ch_vals.append(ch)

    fw_lru = lru.stats['tension'] / max(1, lru.stats['total'])
    harmony = lru.stats['harmony'] / max(1, lru.stats['total'])
    hit_rate = lru.stats['cache_hits'] / max(1, lru.stats['total'])
    return fw_lru, harmony, hit_rate, lru.stats['total']

# ── Test on a sample of genes ──
ENSEMBL_DIR = os.environ.get('ENSEMBL_DIR', 'data/ensembl_release_115')
CDS = os.path.join(ENSEMBL_DIR, 'cds', 'Homo_sapiens.GRCh38.cds.all.fa.gz')

def load_genes(path, max_genes=200):
    gene_seqs = {}
    with gzip.open(path, 'rt', encoding='utf-8', errors='ignore') as f:
        cur_gene, cur_bio, parts = None, None, None
        for line in f:
            if line.startswith('>'):
                if cur_gene and parts and cur_bio == 'protein_coding':
                    seq = ''.join(parts).upper()
                    clean = ''.join(c for c in seq if c in 'ACGT')
                    if len(clean) >= 300 and len(clean) % 3 == 0:
                        key = cur_gene.upper()
                        if key and (key not in gene_seqs or len(clean) > len(gene_seqs.get(key, ''))):
                            gene_seqs[key] = clean
                            if len(gene_seqs) >= max_genes: return gene_seqs
                m_g = re.search(r'gene_symbol:(\S+)', line.strip())
                m_b = re.search(r'gene_biotype:(\S+)', line.strip())
                cur_gene = m_g.group(1) if m_g else None
                cur_bio = m_b.group(1) if m_b else None; parts = []
            else:
                if parts is not None: parts.append(line.strip())
    return gene_seqs

print('Loading genes...')
genes = load_genes(CDS, max_genes=100)
print(f'  {len(genes)} genes')

print('\nComputing SHP+LRU fw...')
results = []
for gene, seq in genes.items():
    r = compute_lru_fw(seq)
    if r:
        results.append((gene, *r))

# ── SHP fw for comparison ──
import sys; sys.path.insert(0, 'src')
from shp.core import compute_shp
print('Computing SHP fw (baseline)...')
for i, (gene, fw_l, harm, hit, total) in enumerate(results[:]):
    seq = genes[gene]
    shp_r = compute_shp(seq, dim=128, window=128, theta0=0.078)
    results[i] = (gene, fw_l, harm, hit, total, shp_r.fixed_wit)

print(f'\n{"="*65}')
print(f'SHP+LRU vs Raw SHP on {len(results)} genes')
print(f'{"="*65}')
print(f'  LRU params: C=8, pattern_len=2, D=128, W=128')
print(f'')
fw_lru_vals = [r[1] for r in results]
fw_shp_vals = [r[5] for r in results]
harm_vals = [r[2] for r in results]
hit_vals = [r[3] for r in results]

import statistics
print(f'  {"Metric":<25} {"Mean":>8} {"Median":>8} {"SD":>8}')
print(f'  {"LRU fw (tension)":<25} {statistics.mean(fw_lru_vals):>8.4f} '
      f'{statistics.median(fw_lru_vals):>8.4f} {statistics.stdev(fw_lru_vals):>8.4f}')
print(f'  {"LRU harmony":<25} {statistics.mean(harm_vals):>8.4f} '
      f'{statistics.median(harm_vals):>8.4f} {statistics.stdev(harm_vals):>8.4f}')
print(f'  {"LRU cache hit rate":<25} {statistics.mean(hit_vals):>8.4f} '
      f'{statistics.median(hit_vals):>8.4f} {statistics.stdev(hit_vals):>8.4f}')
print(f'  {"SHP fw (cross-harm)":<25} {statistics.mean(fw_shp_vals):>8.4f} '
      f'{statistics.median(fw_shp_vals):>8.4f} {statistics.stdev(fw_shp_vals):>8.4f}')

# Correlation
def spearman(x, y):
    n = len(x)
    def rk(v):
        idx = sorted(range(n), key=lambda i: v[i])
        r = [0]*n; i = 0
        while i < n:
            j = i;
            while j < n and v[idx[j]] == v[idx[i]]: j += 1
            for k in range(i, j): r[idx[k]] = (i+j-1)/2
            i = j
        return r
    rx, ry = rk(x), rk(y)
    mx, my = statistics.mean(rx), statistics.mean(ry)
    n1 = sum((rx[i]-mx)*(ry[i]-my) for i in range(n))
    d1 = math.sqrt(sum((rx[i]-mx)**2 for i in range(n)) * sum((ry[i]-my)**2 for i in range(n)))
    return n1 / max(0.001, d1)

print(f'\n  rho(LRU_fw, SHP_fw) = {spearman(fw_lru_vals, fw_shp_vals):+.3f}')
print(f'  rho(harmony, hit_rate) = {spearman(harm_vals, hit_vals):+.3f}')

# Top/bottom genes
print(f'\n  Top-5 LRU tension genes:')
idx = sorted(range(len(results)), key=lambda i: -fw_lru_vals[i])[:5]
for i in idx:
    print(f'    {results[i][0]:<15} LRU_fw={fw_lru_vals[i]:.4f} SHP_fw={fw_shp_vals[i]:.4f} '
          f'harmony={harm_vals[i]:.4f} hits={hit_vals[i]:.4f}')

print(f'\n  Bottom-5 LRU tension genes:')
idx = sorted(range(len(results)), key=lambda i: fw_lru_vals[i])[:5]
for i in idx:
    print(f'    {results[i][0]:<15} LRU_fw={fw_lru_vals[i]:.4f} SHP_fw={fw_shp_vals[i]:.4f} '
          f'harmony={harm_vals[i]:.4f} hits={hit_vals[i]:.4f}')

print('\nDone.')
