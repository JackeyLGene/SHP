# Data Notes

This repository intentionally keeps data small.

Included:

- `data/demo.fa`: toy sequences for checking the CLI.
- `data/gene_matrix_sample.csv`: a compact sample of the GeneGrammar matrix.
- `results/biological_calibration.csv`: calibration table copied from the local
  GeneGrammar working tree.
- `results/calibration_summary.json`: calibration summary copied from the local
  GeneGrammar working tree.

Not included:

- the full 19,491-gene matrix;
- raw GENCODE FASTA files;
- large transcript or genome downloads.

Recommended release plan:

1. Keep this repository lightweight and runnable.
2. Attach full matrices to a GitHub release or Zenodo record.
3. Document the exact GENCODE version and extraction scripts when publishing the
   full dataset.

The paper figures in `paper/figures/` are included as assets for review and
discussion. The script in `paper/scripts/make_paper_figures.py` expects the
local data layout used during the GeneGrammar experiments and may require path
adjustment outside the original workspace.
