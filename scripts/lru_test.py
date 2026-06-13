"""LRU quick test: synthetic repetitive vs random DNA streams."""
import math, random, hashlib, statistics

random.seed(42)

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

class LRUStream:
    def __init__(self, C=8, D=128):
        self.C = C; self.D = D
        self.cache = []
        self.chroma_active = set()
        self.stats = {'total': 0, 'tension': 0, 'harmony': 0, 'neutral': 0, 'hits': 0}

    def _hash(self, text):
        return int(hashlib.md5(text.encode()).hexdigest()[:8], 16) % self.D

    def step(self, win):
        kmers = [win[i:i+3] for i in range(len(win)-2)]
        chroma = {self._hash('C:'+k) for k in kmers}
        rhythm_hashes = [self._hash('R:'+kmers[i]+'>'+kmers[i+1]) for i in range(len(kmers)-1)]
        self.stats['total'] += 1

        chroma_visible = len(chroma & self.chroma_active) > 0 if self.chroma_active else False
        fp_m, fp_n = beta_encode(rhythm_hashes[:2]) if len(rhythm_hashes) >= 2 else (0, 1)
        fp = (fp_m % 10007, fp_n)

        hit = False; rhythm_match = False
        for i, (cfp, _, cnt) in enumerate(self.cache):
            if cfp == fp:
                hit = True; rhythm_match = True
                self.cache.pop(i)
                self.cache.insert(0, (fp, None, cnt + 1))
                break
        if not hit:
            if len(self.cache) >= self.C: self.cache.pop()
            self.cache.insert(0, (fp, None, 1))

        for h in chroma: self.chroma_active.add(h)
        while len(self.chroma_active) > self.C * 2:
            self.chroma_active.pop()

        if chroma_visible and not rhythm_match: self.stats['tension'] += 1
        elif chroma_visible and rhythm_match: self.stats['harmony'] += 1
        else: self.stats['neutral'] += 1
        if hit: self.stats['hits'] += 1

# ── Test 1: highly repetitive DNA (AAAAAAAA...)
rep_seq = 'A' * 5000
print('Test 1: Poly-A (highly repetitive)')
lru = LRUStream(C=8)
for start in range(0, len(rep_seq)-128, 25):
    lru.step(rep_seq[start:start+128])
fw = lru.stats['tension'] / lru.stats['total']
print(f'  tension={lru.stats["tension"]} harmony={lru.stats["harmony"]} '
      f'neutral={lru.stats["neutral"]} hits={lru.stats["hits"]} fw={fw:.4f}')

# ── Test 2: random DNA (IID)
rand_seq = ''.join(random.choice('ACGT') for _ in range(5000))
print('\nTest 2: Random IID DNA')
lru2 = LRUStream(C=8)
for start in range(0, len(rand_seq)-128, 25):
    lru2.step(rand_seq[start:start+128])
fw2 = lru2.stats['tension'] / lru2.stats['total']
print(f'  tension={lru2.stats["tension"]} harmony={lru2.stats["harmony"]} '
      f'neutral={lru2.stats["neutral"]} hits={lru2.stats["hits"]} fw={fw2:.4f}')

# ── Test 3: periodic (ATATAT...)
peri_seq = 'AT' * 2500
print('\nTest 3: AT-periodic')
lru3 = LRUStream(C=8)
for start in range(0, len(peri_seq)-128, 25):
    lru3.step(peri_seq[start:start+128])
fw3 = lru3.stats['tension'] / lru3.stats['total']
print(f'  tension={lru3.stats["tension"]} harmony={lru3.stats["harmony"]} '
      f'neutral={lru3.stats["neutral"]} hits={lru3.stats["hits"]} fw={fw3:.4f}')

# ── Test 4: real CDS (first 5000 bp of BRCA2)
print('\nTest 4: Real CDS-like (random but codon-structured)')
cds_like = ''
codons = ['GCT','GCC','GCA','GCG','CGT','CGC','CGA','CGG','AGA','AGG',  # Ala, Arg mix
          'GGT','GGC','GGA','GGG','CCT','CCC','CCA','CCG',              # Gly, Pro
          'TTA','TTG','CTT','CTC','CTA','CTG',                          # Leu
          'GTT','GTC','GTA','GTG','ATT','ATC','ATA']                    # Val, Ile
for _ in range(200):
    cds_like += random.choice(codons)
lru4 = LRUStream(C=8)
for start in range(0, len(cds_like)-128, 25):
    lru4.step(cds_like[start:start+128])
fw4 = lru4.stats['tension'] / lru4.stats['total']
print(f'  tension={lru4.stats["tension"]} harmony={lru4.stats["harmony"]} '
      f'neutral={lru4.stats["neutral"]} hits={lru4.stats["hits"]} fw={fw4:.4f}')

print(f'\nSummary: Repetitive={fw:.4f} Random={fw2:.4f} Periodic={fw3:.4f} CDS={fw4:.4f}')
print('Done.')
