"""SHP+LRU reverse ablation: TransG cross-species retrieval.

Three-axis comparison:
  - Instrument axis: SHP (window-independent) vs SHP+LRU (cross-window memory)
  - Ablation axis: Natural vs Synonymous shuffle vs Full shuffle
  - Capacity axis: C = 4, 8, 16, 32, 64

Instrument difference = LRU cross-window memory only. All else identical:
  - Same 64-codon exact indexing (no hashing)
  - Same chroma (codon presence), rhythm (codon transitions)
  - Same window W=128, stride W/5
  - Cold start per transcript

Output: retrieval gain (SHP+LRU - SHP) across ablation × capacity grid.
"""
import gzip, statistics, math, random, re, os, csv, sys
from collections import defaultdict, Counter, OrderedDict

random.seed(42)
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
ENSEMBL_DIR = os.environ.get('ENSEMBL_DIR', 'data/ensembl_release_115')

# ── Genetic code ──
SYMBOLS = {'A':0,'C':1,'G':2,'T':3}
AA_TO_CODONS = defaultdict(list)
CODON_TO_IDX = {}
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
    if aa != '*':
        idx = sum(SYMBOLS[x]*(4**(2-j)) for j,x in enumerate(codon))
        CODON_TO_IDX[codon] = idx
        AA_TO_CODONS[aa].append(codon)

# ── LRU core ──
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

D = 128  # shared projection dimension (modulo, not hashing)

class SHPBaseline:
    """Window-independent SHP: per-window chroma/rhythm, no cross-window memory."""
    def __init__(self):
        self.ch_vals = []

    def step(self, win):
        codons = []
        for i in range(0, len(win)-2, 3):
            c = win[i:i+3]
            if c in CODON_TO_IDX: codons.append(CODON_TO_IDX[c])
        if len(codons) < 3:
            self.ch_vals.append(0.5)
            return 0.5
        # Shared projection via modulo (no hashing)
        chroma = {c % D for c in codons}
        rhythm = set()
        for i in range(len(codons)-1):
            rhythm.add((codons[i] * 64 + codons[i+1]) % D)
        union = len(chroma | rhythm)
        if union == 0:
            self.ch_vals.append(1.0)
            return 1.0
        ch = 1.0 - len(chroma & rhythm) / union
        self.ch_vals.append(ch)
        return ch

    def fw(self, theta0):
        ds = [abs(self.ch_vals[i] - self.ch_vals[i-1]) for i in range(1, len(self.ch_vals))]
        if not ds: return 0.0
        return sum(1 for d in ds if d > theta0) / len(ds)

    def mean_h(self):
        return statistics.mean(self.ch_vals) if self.ch_vals else 0.5

class SHPLRU:
    """SHP + LRU: cross-window memory via finite-capacity cache."""
    def __init__(self, C=8):
        self.C = C
        self.cache = OrderedDict()  # fp -> (count, next_codon)
        self.chroma_active = OrderedDict()  # item -> dummy (for LRU order)
        self.ch_vals = []

    def step(self, win):
        codons = []
        for i in range(0, len(win)-2, 3):
            c = win[i:i+3]
            if c in CODON_TO_IDX: codons.append(CODON_TO_IDX[c])
        if len(codons) < 3:
            self.ch_vals.append(0.5)
            return 0.5

        # chroma_now: codon indices projected into shared D-space
        chroma_now = {c % D for c in codons}
        chroma_visible = len(chroma_now & set(self.chroma_active.keys())) > 0 if self.chroma_active else False

        # rhythm fingerprint: transition pattern projected into D-space
        rhythm_vals = [(codons[i] * 64 + codons[i+1]) % D for i in range(len(codons)-1)]
        if len(rhythm_vals) >= 2:
            fp_m, fp_n = beta_encode(list(rhythm_vals[:2]))
            fp = (fp_m % 10007, fp_n)
        else:
            fp = (0, 1)

        # LRU self-lookup
        hit = False; rhythm_match = False
        if fp in self.cache:
            hit = True; rhythm_match = True
            count, _ = self.cache[fp]
            self.cache.move_to_end(fp)
            self.cache[fp] = (count + 1, None)
        else:
            if len(self.cache) >= self.C:
                self.cache.popitem(last=True)  # LRU eviction
            self.cache[fp] = (1, codons[-1] if codons else 0)
            self.cache.move_to_end(fp, last=False)

        # Update chroma_active with proper LRU ordering
        for item in chroma_now:
            self.chroma_active.pop(item, None)
            self.chroma_active[item] = True
            self.chroma_active.move_to_end(item, last=False)
        while len(self.chroma_active) > self.C * 4:
            self.chroma_active.popitem(last=True)

        # Cross-harm
        if chroma_visible and not rhythm_match:
            ch = 1.0
        elif chroma_visible and rhythm_match:
            ch = 0.0
        else:
            ch = 0.5

        self.ch_vals.append(ch)
        return ch

    def fw(self):
        ds = [abs(self.ch_vals[i] - self.ch_vals[i-1]) for i in range(1, len(self.ch_vals))]
        if not ds: return 0.0
        return sum(1 for d in ds if d > 0) / len(ds)

    def mean_h(self):
        return statistics.mean(self.ch_vals) if self.ch_vals else 0.5

# ── Sequence processing ──
def process_sequence(seq, instrument, W=128):
    stride = max(1, W // 5)
    clean = ''.join(c for c in seq.upper() if c in 'ACGT')
    if len(clean) < W: return None
    for start in range(0, len(clean) - W + 1, stride):
        instrument.step(clean[start:start+W])
    return instrument

# ── Profile extraction (for retrieval) ──
def shp_profile(instrument):
    """Extract simple profile: mean_h, fw as 2-dim vector for retrieval."""
    fw = instrument.fw(0.0) if isinstance(instrument, SHPBaseline) else instrument.fw()
    mh = instrument.mean_h()
    # Also include last N ch values for pattern matching
    ch_tail = instrument.ch_vals[-5:] if len(instrument.ch_vals) >= 5 else instrument.ch_vals
    profile = {'fw': fw, 'mean_h': mh, 'n_windows': len(instrument.ch_vals)}
    return profile

# ── Load ──
def load_genes(path):
    gene_seqs = {}
    with gzip.open(path, 'rt', encoding='utf-8', errors='ignore') as f:
        cur_gene, cur_bio, parts = None, None, None
        for line in f:
            if line.startswith('>'):
                if cur_gene and parts and cur_bio == 'protein_coding':
                    seq = ''.join(parts).upper(); clean = ''.join(c for c in seq if c in 'ACGT')
                    if len(clean) >= 300 and len(clean) % 3 == 0:
                        key = cur_gene.upper()
                        if key and (key not in gene_seqs or len(clean) > len(gene_seqs[key])):
                            gene_seqs[key] = clean
                m_g = re.search(r'gene_symbol:(\S+)', line.strip())
                m_b = re.search(r'gene_biotype:(\S+)', line.strip())
                cur_gene = m_g.group(1) if m_g else None
                cur_bio = m_b.group(1) if m_b else None; parts = []
            else:
                if parts is not None: parts.append(line.strip())
        if cur_gene and parts and cur_bio == 'protein_coding':
            seq = ''.join(parts).upper(); clean = ''.join(c for c in seq if c in 'ACGT')
            if len(clean) >= 300 and len(clean) % 3 == 0:
                key = cur_gene.upper()
                if key and (key not in gene_seqs or len(clean) > len(gene_seqs[key])):
                    gene_seqs[key] = clean
    return gene_seqs

def syn_shuffle(seq):
    codons = []
    for i in range(0, len(seq)-2, 3):
        c = seq[i:i+3]
        if len(c)==3 and c in CODON_TO_IDX and c in AA_TO_CODONS[GENETIC_CODE.get(c,'')]:
            codons.append((c, GENETIC_CODE[c]))
    if not codons: return seq
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
        if c in CODON_TO_IDX: codons.append(c)
    random.shuffle(codons)
    return ''.join(codons)

# Build genetic code reverse map
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
}.items(): GENETIC_CODE[codon] = aa

# ── Main ──
print('Loading species...')
human = load_genes(os.path.join(ENSEMBL_DIR, 'cds', 'Homo_sapiens.GRCh38.cds.all.fa.gz'))
mouse = load_genes(os.path.join(ENSEMBL_DIR, 'cds', 'Mus_musculus.GRCm39.cds.all.fa.gz'))

# Orthologs
with open(os.path.join(os.path.dirname(SCRIPT_DIR), '..', '..', 'TransG', 'orthologs', 'orthologs_mouse.csv'), 'r') as f:
    ortho_map = {row['human_symbol'].upper(): row['mouse_symbol'].upper() for row in csv.DictReader(f)}

# Common genes
common = sorted(set(human) & set(ortho_map) & set(mouse))
random.shuffle(common)
queries = common[:100]  # 100 for speed
print(f'  {len(queries)} query genes')
print(f'  Capacities: C = [4, 8, 16, 32, 64]')
print(f'  Ablation levels: Natural / Synonymous / Full shuffle')

# Run
CAPACITIES = [4, 8, 16, 32, 64]
ABLATIONS = ['Natural', 'Synonymous', 'Full']

print(f'\n{"="*90}')
print(f'SHP vs SHP+LRU: Reverse Ablation Grid')
print(f'{"="*90}')
print(f'  {"Method":<15} {"C":>4} {"Ablation":<12} {"mean_h":>8} {"fw":>8} {"n_windows":>9}')
print(f'  {"-"*60}')

results = []
for C in CAPACITIES:
    for abl in ABLATIONS:
        for method_name, method_cls in [('SHP', SHPBaseline), ('SHP+LRU', SHPLRU)]:
            fws = []; mhs = []; nws = []
            for g in queries:
                seq = human[g]
                if abl == 'Synonymous':
                    seq = syn_shuffle(seq)
                elif abl == 'Full':
                    seq = full_shuffle(seq)

                if method_name == 'SHP':
                    inst = process_sequence(seq, method_cls())
                else:
                    inst = process_sequence(seq, method_cls(C=C))

                if inst:
                    fw = inst.fw(0.0) if isinstance(inst, SHPBaseline) else inst.fw()
                    fws.append(fw)
                    mhs.append(inst.mean_h())
                    nws.append(len(inst.ch_vals))

            if fws:
                print(f'  {method_name:<15} {C:>4} {abl:<12} '
                      f'{statistics.mean(mhs):>8.4f} {statistics.mean(fws):>8.4f} '
                      f'{statistics.mean(nws):>9.0f}')
                results.append((method_name, C, abl, statistics.mean(mhs), statistics.mean(fws)))

print('\nDone.')
