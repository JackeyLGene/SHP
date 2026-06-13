"""LRU as Codex: cache the cross-harm TIME SERIES patterns.

SHP produces h_t per window. LRU operates on the h_t sequence to build
a gene-specific structural dynamics codex — which cross-harm patterns recur.

Compare: Natural vs Synonymous CDS — do they build different codexes?
"""
import math, statistics, random, os, sys, gzip, re
from collections import OrderedDict, defaultdict

random.seed(42)
D = 128
W = 128; STRIDE = 25
NG = 3

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
AA_TO_CODONS = defaultdict(list)
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
        AA_TO_CODONS[aa].append(codon)

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

# ── SHP cross-harm ──
def shp_crossharm(win):
    codons = []
    for i in range(0, len(win)-2, 3):
        c = win[i:i+3]
        if c in CODON_TO_IDX: codons.append(CODON_TO_IDX[c])
    if len(codons) < 3: return 0.5
    chroma = {c % D for c in codons}
    rhythm = set()
    for i in range(len(codons)-1):
        rhythm.add((codons[i] * 64 + codons[i+1]) % D)
    union = len(chroma | rhythm)
    return 1.0 - len(chroma & rhythm) / max(1, union)

# ── LRU Codex ──
class CodexLRU:
    """LRU cache on the h_t time series. Cache entries = codex of structural dynamics."""

    def __init__(self, C=16, ctx_len=3):
        self.C = C; self.k = ctx_len
        self.cache = OrderedDict()
        self.h_history = []
        self.stats = {'hits': 0, 'misses': 0, 'total': 0}

    def feed(self, h_val):
        self.h_history.append(h_val)
        self.stats['total'] += 1
        if len(self.h_history) < self.k + 1: return

        # Context: last k h-values, discretized
        ctx = tuple(round(self.h_history[i], 1) for i in range(-self.k-1, -1))
        # Current state transition
        prev = round(self.h_history[-2], 1)
        curr = round(self.h_history[-1], 1)
        transition = (prev, curr)

        # Lookup: does this context->transition exist in cache?
        fp = (beta_encode([hash(str(ctx)) & 0xFFFF])[0] % 10007,
              beta_encode([hash(str(transition)) & 0xFFFF])[0] % 10007)

        if fp in self.cache:
            self.stats['hits'] += 1
            cnt = self.cache[fp]
            self.cache.move_to_end(fp, last=False)
            self.cache[fp] = cnt + 1
        else:
            self.stats['misses'] += 1
            if len(self.cache) >= self.C:
                self.cache.popitem(last=True)
            self.cache[fp] = 1
            self.cache.move_to_end(fp, last=False)

    def hit_rate(self):
        return self.stats['hits'] / max(1, self.stats['total'])

    def cache_size(self):
        return len(self.cache)

    def cache_items(self):
        return list(self.cache.items())[:5]

# ── Load CDS ──
ENSEMBL_DIR = os.environ.get('ENSEMBL_DIR', 'data/ensembl_release_115')
CDS = os.path.join(ENSEMBL_DIR, 'cds', 'Homo_sapiens.GRCh38.cds.all.fa.gz')

print('Loading CDS...')
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
                        if len(gene_seqs) >= 100: break
            m_g = re.search(r'gene_symbol:(\S+)', line.strip())
            m_b = re.search(r'gene_biotype:(\S+)', line.strip())
            cur_gene = m_g.group(1) if m_g else None; cur_bio = m_b.group(1) if m_b else None
            parts = []
        else:
            if parts is not None: parts.append(line.strip())

genes = list(gene_seqs.items())
print(f'  {len(genes)} genes')

# ── Compare codex across ablation ──
print(f'\n{"="*75}')
print(f'LRU CODEX: Structural Dynamics Cache (C=16, ctx=3)')
print(f'{"="*75}')
print(f'  {"Ablation":<12} {"hit_rate":>9} {"cache_size":>10} {"top_fp_examples"}')
print(f'  {"-"*60}')

for abl_name, abl_fn in [('Natural', lambda s: s), ('Synonymous', syn_shuffle), ('Full', full_shuffle)]:
    hit_rates = []; cache_sizes = []
    for gene, seq in genes:
        mod_seq = abl_fn(seq)
        clean = ''.join(c for c in mod_seq if c in 'ACGT')
        codex = CodexLRU(C=16)
        for start in range(0, len(clean)-W+1, STRIDE):
            win = clean[start:start+W]
            h = shp_crossharm(win)
            codex.feed(h)
        hit_rates.append(codex.hit_rate())
        cache_sizes.append(codex.cache_size())

    print(f'  {abl_name:<12} {statistics.mean(hit_rates):>9.4f} {statistics.mean(cache_sizes):>10.1f}')

print('\nDone.')
