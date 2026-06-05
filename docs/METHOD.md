# SHP Method Notes

SHP is a calibrated structural readout for symbolic sequences. The current
reference setting is designed for DNA/RNA strings:

```text
alphabet = ACGT
n-gram   = 3
dim      = 64
window   = 128
theta0   = 0.0999
```

## Encoding

For each local window, SHP builds two binary sets.

`chroma` records which k-mers appear in the window. It asks: what symbols are
active here?

`rhythm` records which adjacent k-mer transitions appear in the same window. It
asks: how does local symbolic state move?

Both sets are hashed into the same fixed dimension. The local cross-harm value is
the Jaccard distance between the two sets:

```text
h_t = 1 - |chroma_t intersect rhythm_t| / |chroma_t union rhythm_t|
```

The structural displacement is the absolute change between consecutive windows:

```text
d_t = |h_t - h_(t-1)|
```

## Calibration

The default threshold `theta0 = 0.0999` is the fair-IID baseline used in the
current GeneGrammar experiments for `k=4, n=3, D=64`. A structural event is
counted when:

```text
d_t > theta0
```

Two summary metrics are usually reported:

```text
fixed_wit   = count(d_t > theta0) / count(d_t)
tail_energy = mean(max(0, d_t - theta0))
```

`fixed_wit` measures event frequency. `tail_energy` measures event magnitude.

## Interpretation

SHP is useful when the same symbolic stream can be read along two orthogonal
axes: local presence and local transition. The method is intentionally minimal:
no training, no alignment model, no external labels, and no genetic-code prior.

Recommended claims:

- SHP identifies structural residuals in symbolic sequence streams.
- SHP can be used as a screening coordinate for genes, regions, and transcript
  classes.
- SHP can suggest candidates for biological interpretation or follow-up tests.

Claims to avoid without further evidence:

- disease diagnosis;
- causal mechanism;
- protein function prediction;
- wet-lab validation;
- universal biological law from a single region set.

## Matched Nulls

For publication-level analysis, use matched nulls. At minimum, compare real
segments against shuffled controls that preserve length and nucleotide
composition. For stronger controls, preserve dinucleotide or codon composition
depending on the biological question.

The current CLI computes the readout. It does not yet generate null models.
