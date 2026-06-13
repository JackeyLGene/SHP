"""Dual-stream LRU: independent chroma cache and rhythm cache.

Chroma LRU: tracks which codon indices are in recent memory.
Rhythm LRU: tracks which transition fingerprints are in recent memory.

Cross-harm = misalignment between the two independent LRU states.
"""
import math, statistics, random, os, sys
from collections import OrderedDict

random.seed(42)
D = 128

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

SYMBOLS = {'A':0,'C':1,'G':2,'T':3}
CODON_TO_IDX = {}
GENETIC_CODE = {}
for codon, aa in {
    'ATA':'I','ATC':'I','ATT':'I','ATG':'M','ACA':'T','ACC':'T','ACG':'T','ACT':'T',
    'AAC':'N','AAT':'N','AAA':'K','AAG':'K','AGC':'S','AGT':'S','AGA':'R','AGG':'R',
    'CTA':'L','CTC':'L','CTG':'L','CTT':'L','CCA':'P','CCC':'P','CCG':'P','CCT':'P',
    'CAC':'H','CAT':'H','CAA':'Q','CAG':'Q','CGA':'R','CGC':'R','CGG':'R','CGT':'R',
    'GTA':'V','GTC':'V','GTG':'V','GTT':'V','GCA':'A','GCC':'A','GCG':'A','GCT':'A',
    'GAC':'D','GAT':'D','GAA':'E','GAG':'E','GGA':'G','GGC':'G','GGG':'G','GGT':'G',
    'TCA':'S','TCC':'S','TCG':'S','TCT':'S','TTC':'F','TTT':'F','TTA':'L','TTG':'L',
    'TAC':'Y','TAT':'Y','TAA':'*','TAG':'*','TGA':'*','TGC':'C','TGT':'C','TGG':'W',
}.items():
    GENETIC_CODE[codon] = aa
    if aa != '*':
        CODON_TO_IDX[codon] = sum(SYMBOLS[x]*(4**(2-j)) for j,x in enumerate(codon))

class DualLRU:
    """Two independent LRU caches. Cross-harm from their misalignment."""

    def __init__(self, C=8):
        self.C = C
        self.chroma_cache = OrderedDict()   # codon_bin -> count
        self.rhythm_cache = OrderedDict()   # fp -> (count, last_codon)
        self.trace = []  # (ch, ch_hit, rh_hit, match) per window

    def step(self, win):
        codons = []
        for i in range(0, len(win)-2, 3):
            c = win[i:i+3]
            if c in CODON_TO_IDX: codons.append(CODON_TO_IDX[c])
        if len(codons) < 3:
            self.trace.append((0.5, False, False, False))
            return 0.5

        # ── Chroma stream: codon bins ──
        chroma_bins = {c % D for c in codons}
        # Update chroma LRU
        for cb in chroma_bins:
            self.chroma_cache.pop(cb, None)
            self.chroma_cache[cb] = self.chroma_cache.get(cb, 0) + 1
            self.chroma_cache.move_to_end(cb, last=False)  # mark as used
        while len(self.chroma_cache) > self.C:
            self.chroma_cache.popitem(last=True)  # LRU evict

        # ── Rhythm stream: transition fingerprints ──
        rhythm_vals = [(codons[i] * 64 + codons[i+1]) % D for i in range(len(codons)-1)]
        fp = None
        if len(rhythm_vals) >= 2:
            fp_m, fp_n = beta_encode(rhythm_vals[:2])
            fp = (fp_m % 10007, fp_n)

        # Update rhythm LRU
        rh_hit = False; rh_last = None
        if fp is not None:
            if fp in self.rhythm_cache:
                rh_hit = True
                cnt, rh_last = self.rhythm_cache[fp]
                self.rhythm_cache.move_to_end(fp, last=False)
                self.rhythm_cache[fp] = (cnt + 1, codons[-1])
            else:
                if len(self.rhythm_cache) >= self.C:
                    self.rhythm_cache.popitem(last=True)
                self.rhythm_cache[fp] = (1, codons[-1])
                self.rhythm_cache.move_to_end(fp, last=False)

        # ── Cross-harm: alignment between streams ──
        ch_hit = len(chroma_bins & set(self.chroma_cache.keys())) > 0
        ch_last = self.rhythm_cache[fp][1] if fp and fp in self.rhythm_cache else None

        # Match: chroma has current codon that matches rhythm's prediction
        match = False
        if ch_hit and rh_hit and ch_last is not None:
            match = (ch_last % D) in chroma_bins

        if ch_hit and rh_hit and match:
            ch = 0.0   # both streams aligned
        elif ch_hit and rh_hit and not match:
            ch = 1.0   # both active but disagree
        elif ch_hit or rh_hit:
            ch = 0.5   # one side active
        else:
            ch = 0.5   # neither active

        self.trace.append((ch, ch_hit, rh_hit, match))
        return ch

    def stats(self):
        if not self.trace: return {}
        ch_vals = [t[0] for t in self.trace]
        m = statistics.mean(ch_vals)
        fw = sum(1 for i in range(1, len(ch_vals)) if abs(ch_vals[i]-ch_vals[i-1]) > 0) / max(1, len(ch_vals)-1)
        ch_hits = sum(1 for t in self.trace if t[1]) / len(self.trace)
        rh_hits = sum(1 for t in self.trace if t[2]) / len(self.trace)
        matches = sum(1 for t in self.trace if t[3]) / len(self.trace)
        return {'mean_h': m, 'fw': fw, 'ch_hit_rate': ch_hits, 'rh_hit_rate': rh_hits, 'match_rate': matches}


# ── Quick test ──
print('Dual-stream LRU validation')
print('='*60)

for label, seq in [
    ('Poly-A repeat', 'A'*5000),
    ('AT periodic', 'AT'*2500),
    ('Random IID', ''.join(random.choice('ACGT') for _ in range(5000))),
    ('Codon-structured', ''.join(random.choice(['GCT','GCC','CGT','CGC','GGT','GGC','CCT','CCC','TTA','TTG']) for _ in range(200))),
]:
    lru = DualLRU(C=8)
    for start in range(0, len(seq)-128, 25):
        lru.step(seq[start:start+128])
    s = lru.stats()
    print(f'  {label:<20} mean_h={s["mean_h"]:.3f} fw={s["fw"]:.3f} '
          f'ch_hit={s["ch_hit_rate"]:.2f} rh_hit={s["rh_hit_rate"]:.2f} match={s["match_rate"]:.2f}')

# ── Real CDS test: Natural vs Synonymous vs Full ──
import gzip, re
from collections import defaultdict

AA_TO_CODONS = defaultdict(list)
for codon, aa in GENETIC_CODE.items():
    if aa != '*': AA_TO_CODONS[aa].append(codon)

def syn_shuffle(seq):
    codons = []
    for i in range(0, len(seq)-2, 3):
        c = seq[i:i+3]
        if c in GENETIC_CODE and GENETIC_CODE[c] != '*':
            codons.append((c, GENETIC_CODE[c]))
    aa_pos = defaultdict(list)
    for pos, (c, aa) in enumerate(codons): aa_pos[aa].append(pos)
    syn_map = {}
    for aa, positions in aa_pos.items():
        choices = [codons[p][0] for p in positions]
        random.shuffle(choices)
        for p, nc in zip(positions, choices): syn_map[p] = nc
    return ''.join(syn_map.get(p, codons[p][0]) for p in range(len(codons)))

def full_shuffle(seq):
    codons = []
    for i in range(0, len(seq)-2, 3):
        c = seq[i:i+3]
        if c in GENETIC_CODE and GENETIC_CODE[c] != '*': codons.append(c)
    random.shuffle(codons); return ''.join(codons)

ENSEMBL_DIR = os.environ.get('ENSEMBL_DIR', 'data/ensembl_release_115')
CDS = os.path.join(ENSEMBL_DIR, 'cds', 'Homo_sapiens.GRCh38.cds.all.fa.gz')

print('\nLoading real CDS...')
gene_seqs = {}
with gzip.open(CDS, 'rt', encoding='utf-8', errors='ignore') as f:
    cur_gene, cur_bio, parts = None, None, None
    for line in f:
        if line.startswith('>'):
            if cur_gene and parts and cur_bio == 'protein_coding':
                seq = ''.join(parts).upper()
                clean = ''.join(c for c in seq if c in 'ACGT')
                if len(clean) >= 300 and len(clean) % 3 == 0:
                    cur_gene = cur_gene.upper() if cur_gene else None
                    if cur_gene and (cur_gene not in gene_seqs or len(clean) > len(gene_seqs[cur_gene])):
                        gene_seqs[cur_gene] = clean
                        if len(gene_seqs) >= 50: break
            m_g = re.search(r'gene_symbol:(\S+)', line.strip())
            m_b = re.search(r'gene_biotype:(\S+)', line.strip())
            cur_gene = m_g.group(1) if m_g else None
            cur_bio = m_b.group(1) if m_b else None; parts = []
        else:
            if parts is not None: parts.append(line.strip())

genes = list(gene_seqs.items())
print(f'  {len(genes)} genes')

print(f'\n{"="*75}')
print(f'REAL CDS: Dual-Stream LRU — Natural vs Synonymous vs Full (C=8)')
print(f'{"="*75}')
print(f'  {"Ablation":<12} {"mean_h":>8} {"fw":>8} {"ch_hit":>7} {"rh_hit":>7} {"match":>7}')
print(f'  {"-"*50}')

for abl_name, abl_fn in [('Natural', lambda s: s), ('Synonymous', syn_shuffle), ('Full', full_shuffle)]:
    all_stats = []
    for gene, seq in genes:
        mod_seq = abl_fn(seq)
        lru = DualLRU(C=8)
        clean = ''.join(c for c in mod_seq if c in 'ACGT')
        for start in range(0, len(clean)-128, max(1, 128//5)):
            lru.step(clean[start:start+128])
        all_stats.append(lru.stats())

    print(f'  {abl_name:<12} '
          f'{statistics.mean([s["mean_h"] for s in all_stats]):>8.4f} '
          f'{statistics.mean([s["fw"] for s in all_stats]):>8.4f} '
          f'{statistics.mean([s["ch_hit_rate"] for s in all_stats]):>7.3f} '
          f'{statistics.mean([s["rh_hit_rate"] for s in all_stats]):>7.3f} '
          f'{statistics.mean([s["match_rate"] for s in all_stats]):>7.3f}')

print('\nDone.')

