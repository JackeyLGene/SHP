# An annotation-free structural readout of CDS and UTR sequence regimes in the human genome

**Jieqi Liu**  
Independent researcher  
Correspondence: jackey.l.gene@outlook.com  
June 2026

---

## Abstract

Genomic sequence analysis typically begins with annotation: gene models, functional
domains, expression data, or evolutionary conservation. We introduce SHP (Saussurean
Hash Projection), a lightweight, zero-training structural readout for nucleotide
sequences. Feature computation is annotation-free after CDS/UTR sequence extraction:
Ensembl transcript annotations are used to define sequence regions, but no functional
labels, expression data, conservation scores, or disease annotations enter the SHP
calculation. SHP encodes each sequence window through two complementary views of the
same 3-mer stream -- a chroma axis (which 3-mers are present) and a rhythm axis (which
3-mer transitions occur) -- into the same 64-dimensional hash space, then measures
structural tension via the Jaccard distance between them. The instrument is calibrated
against a fair IID baseline (k = 4, theta0 = 0.0999).

Applied to the full human protein-coding genome (19,491 genes, 224,518 transcript
isoforms, Ensembl release 115), SHP produces an 8-dimensional per-gene structural
feature vector (fixed_wit, tail_energy, skew, kurt for both CDS and UTR) plus derived
gradients, without consulting any external functional labels. We report the following
findings: (1) the CDS structural quiescence rate is 48.8% genome-wide, with
regime-specific signatures -- MHC class I genes show 80% CDS quiescence, while neural
genes show only 20%; (2) a UTR/CDS structural gradient distinguishes functional
regimes, with KRTAP genes as the only CDS-led regime (85.7% UTR quiescence); (3) SHP
vectors carry functional information: nearest-centroid classification achieves 30-40%
per-category recall for olfactory receptors, HLA genes, and keratins from 2D SHP
features alone; (4) this clustering is not a GC-content proxy -- zinc-finger genes
(GC=47%) and transcription factors (GC=63%) are separated by 16 percentage points of GC
yet cluster together in SHP space (d = 0.004). The method requires no training, no
model weights, and no external functional databases beyond the input sequence files.
Code, demo FASTA input, calibration summaries, figure-generation scripts, and the
per-gene and per-isoform SHP matrices are provided for replication. We provide a
calibrated, annotation-free structural spectroscopy tool for genomic sequence
screening.

---

## Introduction

Modern genomics produces vast catalogs of sequences — transcript isoforms, variant
calls, regulatory regions — but experimental characterization remains expensive.
Computational triage typically relies on sequence homology, functional domain databases,
expression correlation, codon-usage summaries, or evolutionary conservation scores
(for example CAI (Sharp and Li, 1987), GO/InterPro (Gene Ontology Consortium, 2023;
Blum et al., 2025), or phyloP-style conservation scores (Pollard et al., 2010)).
Each of these carries prior assumptions about what constitutes functional importance.

Here we introduce an orthogonal approach: direct structural spectroscopy of the
nucleotide symbol stream itself. The method, SHP (Saussurean Hash Projection), asks a
single question of any DNA sequence: *does the local 3-mer composition (chroma) and the
local 3-mer transition pattern (rhythm) exhibit structural tension beyond what a random
sequence would produce?*

The name is a compact mnemonic rather than a biological assumption. Its two axes
loosely echo the paradigmatic/syntagmatic distinction in structural linguistics
(Saussure, 1983):
what symbols are locally available, and how symbols are chained. The algorithm itself
requires no linguistic theory; it only compares binary k-mer presence and binary
k-mer transition activations in a shared hash space.

The present paper treats SHP as a standalone calibrated screening instrument for
genomic sequences. No prior knowledge of the broader symbolic-processing framework in
which the dual-axis encoding was developed is required to use or interpret the results.

The core measurement is straightforward. A sliding window of 128 nucleotides is
projected into a 64-dimensional hash space twice: once encoding which 3-mers appear
(chroma), and once encoding which 3-mer-to-3-mer transitions occur (rhythm). The
Jaccard distance $1 - |C_t \cap R_t| / |C_t \cup R_t|$ (Jaccard, 1901) between these
two binary activation patterns---the *cross-harm* $h_t$---measures the structural
tension between the paradigmatic and syntagmatic axes of the local sequence grammar.
A structural event is registered when the cross-harm displacement
$d_t = |h_t - h_{t-1}|$ exceeds the 99th percentile of fluctuations observed in a
maximally unstructured (fair IID) stream of the same alphabet.

The instrument was calibrated on synthetic symbol streams spanning GC content from
40--70%, CpG suppression, and codon periodicity; the resolution threshold
($\theta_0 = 0.0999$ for 128-nt windows) proved stable across these controls.
We then applied it to all human protein-coding genes with measurable CDS and UTR
(Ensembl release 115, 19,491 genes), producing per-gene structural vectors without
using any gene ontology terms, disease labels, expression data, or conservation scores.

The resulting SHP matrix reveals regime-level structural signatures, carries
functional information orthogonal to GC content, and identifies structurally extreme
genes and isoforms that no external annotation would have flagged. We present the
method, the calibration, the genome-scale results, and a practical tool for
annotation-free genomic sequence screening.

---

## Results

![Figure 1. SHP turns one nucleotide stream into two structural views. Chroma asks which 3-mers are present; rhythm asks how 3-mers transition. Cross-harm measures their mismatch.](figures/fig1_shp_method.png)

### SHP Calibration and Measurement Framework

SHP operates on any discrete symbol stream with alphabet size $k \geq 4$. For DNA ($k=4$),
a sliding window of $W = 128$ nucleotides is analyzed at $n = 3$ (trinucleotide) granularity
with $D = 64$ hash buckets, the combinatorial saturation point $k^3 = D$. Within each window,
chroma and rhythm vectors are computed as binary hash activations (presence/absence, not
frequency-weighted). Cross-harm $h_t = 1 - J(C_t, R_t)$ is computed per
window, where $J$ denotes the Jaccard similarity of the two hash-activation sets.
Displacement $d_t = |h_t - h_{t-1}|$ between consecutive windows forms the
primary signal trace.

The instrument is calibrated against a fair IID baseline: 10,000-symbol streams with
$k = 4$ symbols drawn independently with uniform probability, repeated across 20 random
seeds. The 99th percentile of cross-harm displacement under this baseline defines the
resolution threshold $\theta_0 = 0.0999 \pm 0.0067$ (mean $\pm$ 1 SD across seeds).

From this calibration, we define two primary SHP metrics per sequence:

- $\text{fixed\_wit} = \text{COUNT}(d_t > \theta_0) \,/\, N_{\text{transitions}}$ — the structural event rate, measuring how frequently the local 3-mer grammar undergoes reorganization beyond the random-fluctuation baseline.
- $\text{tail\_energy} = \overline{\max(0, d_t - \theta_0)}$ — the structural event intensity, measuring the magnitude of excess displacement above baseline.

A third, derived metric---the distribution shape (skew, kurt) of $\{d_t\}$---captures the
structural phase of the sequence (periodic, random, Markov-biased, block-structured).

We validated $\theta_0$ against DNA-specific composition constraints: GC content from
40--70%, CpG-suppressed dinucleotide Markov chains, and codon-periodic synthetic CDS
all produce $\theta_0$ within 5% of the fair-IID reference (no pairwise $t$-test significant
at Bonferroni-corrected $\alpha = 0.01$). The threshold is stable for the tested DNA
parameter setting ($k = 4$, $n = 3$, $D = 64$, $W = 128$).

A critical methodological finding: when each stream uses its own P99 as threshold
(self-thresholding), all streams produce a $\text{fixed\_wit}$ of approximately $0.007$. This
apparent constancy is a statistical artifact of percentile-based thresholding. Under
the fixed $\theta_0$, $\text{fixed\_wit}$ spans a 10-fold range across stream types.
All results reported below use the fixed $\theta_0 = 0.0999$.

### Genome-Wide SHP Matrix

We processed all 19,491 human protein-coding genes with measurable CDS ($\geq 128$ nt)
from Ensembl release 115, plus 16,572 genes with measurable UTR. For each gene, SHP
metrics were computed for the longest CDS transcript and its matching UTR (extracted by
transcript-ID alignment of cDNA minus CDS). An additional 224,518 individual transcript
isoforms were processed for isoform-level analysis.

After applying a minimum-window filter ($n_{\text{windows}} \geq 10$, removing 659 genes
with insufficient window coverage), the genome-wide distribution of SHP metrics is:

| Metric | CDS ($n_w \geq 10$) | UTR ($n_w \geq 10$) |
|--------|---------------------|---------------------|
| Genes | 18,832 | 15,278 |
| Mean $\text{fixed\_wit}$ | 0.0165 | 0.0201 |
| Exactly zero ($\text{fw}=0$) | 44.8% | 40.1% |
| Quiescent ($\text{fw} < 0.01$) | 48.8% | 44.2% |
| Mean gradient (UTR$-$CDS), both-region genes ($N = 14{,}996$) | --- | $+0.0035$ |

Approximately 45% of protein-coding genes show exactly zero structural events in
their CDS (48.8% below the 0.01 quiescence threshold). UTRs are consistently more
structurally active than CDSs in the both-region subset (mean gradient +0.0035,
N=14,996 genes with both CDS and UTR passing the window filter).

![Figure 2. Calibration and genome-wide structural event rates. (A) Fair-IID displacement distribution defining theta0=0.0999. (B) Fixed-threshold readout for control streams. (C) Genome-wide CDS and UTR fixed_wit distributions. (D) Exact-zero and quiescent fractions.](figures/fig2_calibration_genome.png)

### Regime-Level Structural Signatures

We examined five biological regimes with expanded gene sets (N >= 20 per regime),
selected to represent distinct functional categories. Gene lists are documented in the
Methods and provided in the supplementary material.

| Regime | N | CDS fw | UTR fw | CDS<0.01% | UTR<0.01% |
|--------|---|--------|--------|--------|--------|
| HOX (developmental) | 39 | 0.0104 | 0.0206 | 69.2% | 43.6% |
| MHC (immune diversity) | 20 | 0.0096 | 0.0303 | 80.0% | 44.4% |
| BRAIN (neural) | 25 | 0.0214 | 0.0294 | 20.0% | 36.4% |
| AD (maintenance) | 20 | 0.0115 | 0.0200 | 55.0% | 52.6% |
| KRTAP (production) | 75 | 0.0341 | 0.0110 | 58.7% | 85.7% |

**MHC genes show the highest CDS quiescence of any regime (80% fw<0.01).** Class I and
class II HLA genes have structurally stable coding templates, with elevated structural
activity in UTR regions (UTR fw=0.0303, highest of any regime). Whether this reflects
selection for stable antigen-presentation domains coupled with regulated expression
is a hypothesis for further investigation.

**Neural genes show the lowest CDS quiescence (20% fw<0.01).** Only one in five
brain-expressed genes has a quiescent CDS template. Neural genes carry structural
grammar throughout their coding regions — a phenotype consistent with extensive
alternative splicing documented in neural tissue, though SHP does not directly
measure splicing.

**KRTAP genes are the only CDS-led regime.** Production-type keratin-associated
proteins show the highest UTR quiescence (85.7% fw<0.01) and a strongly negative
UTR-CDS gradient (-0.027). Their structural grammar is concentrated in the coding
region.

**AD maintenance genes show the highest UTR quiescence among non-production regimes
(52.6% fw<0.01).** This observation is reported as an SHP phenotype; whether it
reflects regulatory grammar changes associated with maintenance failure requires
independent biological validation.

Within the expanded PRB/KRTAP production-gene family (N=100), we identified four
structural phases across the CDS fixed_wit vs UTR fixed_wit plane: QUIET (81% of genes,
both axes inactive), UTR-CONTROLLED (11%, CDS=0 + UTR>0.01), CDS-ACTIVE (5%), and
DUAL-ACTIVE (3%). The UTR-CONTROLLED phase, exemplified by PRB3 (UTR fixed_wit=0.333,
39x baseline), was validated by dinucleotide-preserving shuffle (null fixed_wit=0.004,
ratio=20x), confirming the signal is transition-driven rather than a composition
artifact.

![Figure 3. Biological regimes occupy distinct CDS/UTR structural states. (A) Mean CDS and UTR fixed_wit by regime. (B) Regime centroids in CDS/UTR structural space. (C) Min-max normalized heat map of regime signatures.](figures/fig3_regime_signatures.png)

### SHP Vectors Carry Functional Information

We tested whether SHP structural vectors carry functional signal by assigning broad
categories based on gene symbol prefixes (e.g., all ZNF* genes as "ZincFinger", all
SLC* as "SoluteCarrier") — a deliberately crude classification using zero external
annotation databases. Category labels were used only for validation, not as input
features.

Using only 2-dimensional SHP vectors (CDS fixed_wit, UTR fixed_wit), nearest-centroid
leave-one-out classification achieved 13.2% overall accuracy vs.\ 11.1% random baseline (9 categories) across
categories (N=2,034 labeled genes). Per-category recall: OlfactoryReceptor 40.5%, HLA
33.3%, Keratin 30.8% — three categories at 3--4x chance level from pure 3-mer
structural grammar alone. GPCR and Histone categories showed near-chance performance in
this low-dimensional baseline.

**Crucially, SHP clustering is not a GC-content proxy.** ZincFinger genes (GC=47.4%)
and TranscriptionFactors (GC=62.6%) are separated by 16 percentage points of GC —
yet in unscaled SHP fixed_wit space, their centroids are nearly overlapping
(dist_SHP=0.004). Both are
DNA-binding proteins. In 8 of 9 discriminative category pairs, SHP clustering disagreed
with GC+length clustering. OlfactoryReceptors (GC=47.0%) and Histones (GC=60.1%) are
separated by 13pp of GC yet cluster together in SHP space. Bootstrap validation
(500 resamples, equal-N=12 subsampling, 36 pairwise comparisons) confirmed all 36
category pairs have 95% CI lower bounds excluding zero.

A label permutation test (200 shuffles) confirmed the observed 13.2% classification
accuracy exceeds the null distribution (mean null accuracy $6.9\% \pm 2.5\%$, $Z = 2.5$, $p < 0.01$). The
modest absolute accuracy reflects the deliberately low-dimensional feature space
(2 of 8 available SHP dimensions).

![Figure 4. SHP carries functional signal orthogonal to simple composition. (A) Functional categories in SHP space. (B) SHP distance vs GC difference for selected pairs. (C) Leave-one-out nearest-centroid recall.](figures/fig4_functional_orthogonality.png)

### Isoform Structural Diversity

Within-gene isoform SHP dispersion was computed for 17,111 genes with >=2 isoforms.
Mean within-gene fixed_wit range = 0.037. Top genes by isoform structural diversity:
PPFIBP2 (20 isoforms, range 0.6), G3BP2 (73 isoforms, range 0.4), PUM1 (220 isoforms,
range 0.4). Notably, 27.5% of multi-isoform genes have ALL isoforms at fixed_wit=0 —
splicing produces zero structural grammar variation among their transcripts.

### Chromosome-Scale Distribution

SHP metrics were mapped to genomic coordinates using GENCODE v49 gene annotations.
All autosomes occupy a narrow band of mean CDS fixed_wit (0.0148--0.0182). Acrocentric
chromosomes 21 and 22 show elevated short-arm/long-arm fixed_wit ratios (1.52 and 1.42
respectively). Chromosome Y is the most boundary-rich (16 structural boundaries on 39
genes), consistent with palindromic repeat structure concentrating cross-harm
reorganization at repeat boundaries.

Telomere proximity shows no effect on SHP metrics: CDS fixed_wit is flat across all
10 position bins from telomere to centromere (mean 0.016 +/- 0.001). This null result
confirms SHP reads sequence-intrinsic 3-mer grammar — telomere position effects operate
on expression state, not on CDS sequence structure.

---

## Discussion

### What SHP Measures

SHP occupies an unusual position in the genomic analysis toolkit. It is not a
composition measure (GC content, codon adaptation index), not a conservation measure
(phyloP, dN/dS), and not a functional predictor (GO term enrichment, domain detection).
It measures *structural grammar* — the tension between which 3-mers are present (the
paradigmatic axis) and how they transition (the syntagmatic axis) in a local sequence
window. A sequence with high structural tension has a 3-mer composition that does not
imply its 3-mer transition pattern; a sequence with zero structural events has these
two axes in stable equilibrium.

This makes SHP complementary to existing tools. A gene can be highly conserved (strong
phyloP score) yet have zero SHP structural events (stable 3-mer grammar). Conversely, a
gene can show extreme SHP structural activity yet have no known disease association.
The annotation-free nature of SHP means it can flag structurally unusual sequences
that external annotations have not yet characterized — a screening function orthogonal
to knowledge-based prioritization.

### Limitations

**Window size sensitivity.** The standard W=128 nt window (~42 codons) limits UTR
analysis. Short UTRs (<128 nt) are excluded from the current matrix. A W=64 sensitivity
mode is calibrated (theta_0=0.1214) but was not applied genome-wide in this study.

**Two-dimensional baseline.** The functional clustering results use only CDS fixed_wit
and UTR fixed_wit. The full 8-dimensional primary SHP vector (fixed_wit, tail_energy,
skew, and kurt for both CDS and UTR) plus derived gradients may improve classification
accuracy. The 13.2% overall accuracy reported here should be interpreted as a
conservative lower bound.

**Gene symbol classification.** Functional categories were assigned by gene symbol
prefix — a deliberately crude method. Formal GO/InterPro enrichment analysis is
needed to validate the functional signal with standard ontologies.

**Cross-species validation.** The current matrix is human-only. The CDS quiescence
rate of 48.8% and regime signatures need to be tested in other species. Preliminary
primate analysis (chimpanzee, gorilla, orangutan) suggests the production-gene CDS
quiescence rate is conserved (77--82%), but this requires raw genomic DNA rather than
CDS annotation pipelines for definitive cross-species comparison.

### Practical Utility

SHP is intentionally lightweight. The entire genome scan (19,491 genes, 224,518
isoforms) completes in approximately 1.7 hours (~6,200 CPU-seconds) on a single core.
No GPU, no model weights, and no
external functional databases are required once the input sequence files are available.
The method is implementable in ~200 lines of Python using only the standard library
(hashlib, random, math).

For a bioinformatics user, SHP provides an orthogonal screening coordinate: given a
gene list, SHP can rank genes by structural quiescence (CDS fixed_wit), regulatory
structure (UTR fixed_wit), structural gradient (UTR-CDS difference), or isoform
structural diversity (within-gene fixed_wit range). These rankings require no prior
knowledge about gene function and can be computed for any organism with a CDS FASTA.

---

## Methods

### SHP Encoding

For a nucleotide sequence of length $L$, a sliding window of $W = 128$ nucleotides is
applied with stride $\max(1, W/5)$. For each window at position $t$:

**Chroma vector:** each overlapping 3-mer in the window is hashed via MD5, and the
corresponding bucket in a $D = 64$ binary vector is set to 1. Formally,
$C_t[i] = 1$ if any 3-mer in window $t$ hashes to bucket $i$, else $0$.

**Rhythm vector:** each adjacent 3-mer pair $(\text{3-mer}_i, \text{3-mer}_{i+1})$ is
hashed with a distinct label prefix, and the corresponding bucket in the same $D = 64$
hash space is set to 1. Formally, $R_t[i] = 1$ if any 3-mer transition in window $t$
hashes to bucket $i$, else $0$.

**Cross-harm:** $h_t = 1 - |C_t \cap R_t| \,/\, |C_t \cup R_t|$, the Jaccard distance
between the two hash-activation sets.

**Displacement:** $d_t = |h_t - h_{t-1}|$ for $t \geq 2$; $d_1 = 0$.

### Calibration

A fair IID stream of $T = 10{,}000$ symbols from the $k = 4$ alphabet $\{\text{A, C, G, T}\}$
with uniform probability $p = 0.25$ per symbol is generated. $d_t$ values are computed
across sliding windows. The 99th percentile of the empirical $d_t$ distribution,
averaged across 20 independent random seeds, defines
$\theta_0 = 0.0999 \pm 0.0067$.

For $W = 64$ (short-UTR sensitivity mode), the same procedure yields
$\theta_0 = 0.1214 \pm 0.0082$.

### Biological Data

Human CDS sequences were obtained from Ensembl release 115 (Dyer et al., 2025)
(Homo_sapiens.GRCh38.cds.all.fa.gz). Human cDNA sequences from the same release
(Homo_sapiens.GRCh38.cdna.all.fa.gz) were used for UTR extraction. CDS and cDNA
were matched by transcript ID (ENST*); UTR was extracted as the concatenation of
5' and 3' cDNA regions flanking the CDS match. For genes with multiple transcripts,
the longest CDS with a transcript-matched cDNA was used for gene-level analysis.

### Gene Regime Definitions

- **HOX:** 39 genes (HOXA*, HOXB*, HOXC*, HOXD*)
- **MHC:** 20 genes (HLA-A, -B, -C, -E, -F, -G, -H, -J, -L; HLA-DMA, -DMB, -DOA,
  -DOB, -DPA1, -DPB1, -DQA1, -DQA2, -DQB1, -DRA, -DRB1, -DRB3, -DRB4, -DRB5)
- **BRAIN:** 25 genes (FOXP1/2, PAX6, NEUROD1/2/4/6, TBR1, DLG4, SYN1/2/3,
  GABRA1-6, GRIN1/2A/2B/2C/2D/3A/3B)
- **AD:** 20 genes (MAPT, GSK3B, LAMP1, CD33, BCL2, APP, PSEN1/2, APOE, TREM2,
  CLU, CR1, BIN1, PICALM, ABCA7, MS4A6A, EPHA1, CD2AP, SORL1, FERMT2)
- **KRTAP:** 75 genes (all KRTAP* with measurable CDS and UTR)

### Functional Categories

Gene symbol prefixes were used for broad functional classification without external
databases: ZNF* (ZincFinger, N=515), SLC* (SoluteCarrier, N=385), OR* with <=6
characters (OlfactoryReceptor, N=425), KRT* (Keratin, N=120), RPL*/RPS*/MRPL*/MRPS*
(Ribosomal, N=158), GPCR-related prefixes (N=330), HIST*/H1*/H2A*/H2B*/H3*/H4*
(Histone, N=77), HLA-* (N=21), and a curated list of transcription factor prefixes
(FOX*, SOX*, PAX*, TBX*, HOX*, NKX*, LHX*, GATA*, STAT*, SMAD*, NFKB*, TCF*, N=208).

### Dinucleotide-Preserving Null

For each top outlier gene, the CDS sequence was shuffled 10 times using a
dinucleotide-preserving algorithm (preserving all 16 dinucleotide frequencies while
randomizing higher-order structure). SHP metrics were recomputed on each shuffled
sequence; the mean shuffled fixed_wit served as the null expectation.

### Statistical Analysis

Bootstrap confidence intervals: 500 resamples with equal-N=12 subsampling per
category pair. Nearest-centroid classification: leave-one-out with Euclidean
distance in (CDS fixed_wit, UTR fixed_wit) space. Label permutation test: 200
shuffles of category labels, recomputing LOO accuracy for each shuffle.

---

## Data and Code Availability

Standalone SHP code, demo FASTA input, calibration summaries, figure assets, and
preprint materials are provided in the SHP repository:
https://github.com/JackeyLGene/SHP. The core SHP implementation and command-line
scanner are implemented in Python 3.12 using only the standard library (hashlib,
random, math, gzip, csv, json). Figure-generation scripts may require standard
scientific plotting packages such as matplotlib and numpy.

CDS and cDNA FASTA files were obtained from Ensembl release 115 (human) and Ensembl
Genomes release 62 (non-human species). GENCODE v49 GTF was used for chromosome
coordinate mapping (Mudge et al., 2025).

The per-gene SHP matrix (19,491 genes: 8 primary SHP features + 2 gradients +
5 metadata columns) and per-isoform matrix (224,518 transcripts: 4 SHP features
+ 3 metadata columns) are provided as supplementary CSV files.

Questions, replication reports, biological interpretation notes, and requests for
additional matrix files can be sent to jackey.l.gene@outlook.com.

---

## References

Blum M, Andreeva A, Florentino LC, et al. InterPro: the protein sequence classification
resource in 2025. *Nucleic Acids Research*. 2025;53(D1):D444-D456.

Dyer SC, Austine-Orimoloye O, Azov AG, et al. Ensembl 2025. *Nucleic Acids Research*.
2025;53(D1):D948-D957.

Gene Ontology Consortium. The Gene Ontology knowledgebase in 2023. *Genetics*.
2023;224(1):iyad031.

Jaccard P. Etude comparative de la distribution florale dans une portion des Alpes
et du Jura. *Bulletin de la Societe Vaudoise des Sciences Naturelles*. 1901;37:547-579.

Mudge JM, Carbonell-Sala S, Diekhans M, et al. GENCODE 2025: reference gene annotation
for human and mouse. *Nucleic Acids Research*. 2025;53(D1):D966-D975.

Pollard KS, Hubisz MJ, Rosenbloom KR, Siepel A. Detection of nonneutral substitution
rates on mammalian phylogenies. *Genome Research*. 2010;20(1):110-121.

Saussure F de. *Course in General Linguistics*. Edited by Bally C and Sechehaye A.
Translated by Harris R. Open Court; 1983.

Sharp PM, Li WH. The codon adaptation index--a measure of directional synonymous
codon usage bias, and its potential applications. *Nucleic Acids Research*.
1987;15(3):1281-1295.

---

## Acknowledgements

This work emerged from a broader computational research program on symbolic sequence
processing and frame-economy dynamics. The SHP encoding inherits its dual-axis design
from PhasePrompt's chroma-rhythm NLP encoder. The calibration methodology builds on
the fair-baseline measurement framework developed in earlier work (GEME/BGM). The
author thanks the developers of Ensembl and GENCODE for maintaining open genomic
data resources.
