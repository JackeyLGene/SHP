# Figure Captions

## Figure 1. SHP turns one nucleotide stream into two structural views

SHP projects each local nucleotide window into two binary hash activations:
chroma, which records which 3-mers are present, and rhythm, which records which
adjacent 3-mer transitions occur. Cross-harm is the Jaccard distance between
these two views. Consecutive cross-harm displacement defines calibrated
structural events (`fixed_wit`) and excess event intensity (`tail_energy`),
which are assembled into CDS/UTR readout matrices for downstream screening.

## Figure 2. Calibration and genome-wide structural event rates

(A) Fair-IID displacement distribution used to define the fixed threshold
`theta0 = 0.0999`. (B) Fixed-threshold readout separates periodic, fair-IID,
Markov-biased, and block-structured streams. (C) Genome-wide CDS and UTR
`fixed_wit` distributions after the `n_windows >= 10` filter. (D) Exact-zero
and quiescent (`fixed_wit < 0.01`) fractions for CDS and UTR regions.

## Figure 3. Biological regimes occupy distinct CDS/UTR structural states

(A) Mean CDS and UTR `fixed_wit` by regime. (B) Regime centroids in CDS/UTR
structural space; points above the diagonal are UTR-led, points below are
CDS-led. (C) Min-max normalized heat map of regime signatures, with raw values
shown in each cell.

## Figure 4. SHP carries functional signal orthogonal to simple composition

(A) Broad functional gene-symbol categories plotted in CDS/UTR SHP space, with
point size proportional to category size and color indicating CDS GC content.
(B) Selected category pairs show that proximity in SHP space does not reduce to
GC similarity. (C) Leave-one-out nearest-centroid recall using only two SHP
features (`CDS fixed_wit`, `UTR fixed_wit`).
