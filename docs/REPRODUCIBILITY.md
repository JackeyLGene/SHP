# Reproducibility Checklist

Minimal local check:

```powershell
cd G:\SHP
$env:PYTHONPATH = "src"
python -m shp scan --fasta data\demo.fa --out results\demo_scan.tsv
```

Expected behavior:

- the command writes a TSV file;
- `periodic_acgt` has zero displacement;
- at least one mixed sequence has non-zero `mean_d`;
- records shorter than the selected window return zero windows and zero events.

Publication-level checks:

1. Report the exact sequence source and version, for example GENCODE release.
2. Keep CDS and UTR extraction rules explicit.
3. Use matched nulls for claims about biological regimes.
4. Separate feature computation from biological labels.
5. Treat SHP as a screening coordinate unless independent biological validation
   is available.

Current repository status:

- code path: included;
- toy FASTA: included;
- sample matrix: included;
- paper figures: included;
- full genome-scale matrix: pending external release.
