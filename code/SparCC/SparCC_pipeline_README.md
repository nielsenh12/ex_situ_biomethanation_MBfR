# SparCC Correlation Pipeline (`SparCC.ipynb`)

This notebook computes SparCC pairwise correlations between ASVs from raw
count data, tests them for significance via bootstrap permutation, and
filters down to a final high-confidence correlation network. The core
SparCC algorithm is implemented in `sparcc_core.py` and imported into the
notebook rather than reimplemented there.

All paths below are anonymized placeholders — substitute your own
directories when reproducing this pipeline.

## Method & code reference

- **Method**: Friedman J, Alm EJ (2012). Inferring Correlation Networks
  from Genomic Survey Data. *PLoS Computational Biology*, 8(9): e1002687.
  https://doi.org/10.1371/journal.pcbi.1002687
- **Code**: Core SparCC functions (`sparcc`, `permute_w_replacement`)
  sourced from: https://github.com/bio-developer/sparcc

## Inputs

| File | Description |
|---|---|
| `<data_dir>/raw_counts.tsv` | ASV × sample raw count table (2,571 ASVs, 16 samples including Day 0) |
| `<reference_dir>/iterativeIDs.json` | ASV hash → human-readable iterativeID map |
| `<reference_dir>/table1_codif_output_050526.csv` | Taxonomy table (Kingdom–Species) per ASV |

## Pipeline & imposed thresholds

1. **Rarefaction + prevalence filtering**: drop Day 0, rarefy all
   remaining samples to `RAREFY_DEPTH = 22677` (the minimum per-sample
   total in the raw data), then keep only ASVs present in **≥ 3 samples**
   post-rarefaction.
2. **SparCC correlation estimation**: `sparcc(counts, iters=20, th=0.1,
   xiter=10)` on the rarefied table.
3. **Significance testing**: 100 bootstrap null permutations
   (`permute_w_replacement`), each re-run through `sparcc()` with the same
   parameters, to derive an empirical two-sided p-value per ASV pair.
4. **Multiple-testing correction & filtering**: Benjamini–Hochberg FDR
   correction across all pairs, filtered to **FDR q ≤ 0.01**, then further
   filtered to **|correlation| ≥ 0.5** for the final network.
5. **Annotation**: final network pairs are joined to iterativeIDs and
   taxonomy for both ASVs in each pair.

## Outputs

| File | Description |
|---|---|
| `rarefied_counts.tsv` | Rarefied, prevalence-filtered counts |
| `cor_sparcc_final.tsv` | Real SparCC correlation matrix |
| `permutations_final/permutation_{0..99}.tsv` | Bootstrap null permutations |
| `iters_final/perm_cor_{0..99}.npy` | SparCC correlations per permutation |
| `sparcc_final_results_v2.csv` | All pairs with correlation + empirical p-value |
| `sparcc_final_results_with_fdr.csv` | Same, with BH FDR q-values |
| `sparcc_fdr_significant_pairs.csv` | Pairs passing FDR q ≤ 0.01 |
| `sparcc_fdr_significant_pairs_r_filt.csv` | Final network: FDR q ≤ 0.01 and \|r\| ≥ 0.5 |
| `sparcc_final_iterativeIDs.csv` | Final network annotated with iterativeIDs + taxonomy |
