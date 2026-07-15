# SparCC Correlation Network Construction

This repository rarefies and filters an ASV count table, computes SparCC correlations between ASVs, builds an empirical null distribution via bootstrap permutation, and applies FDR correction to produce a statistically filtered co-occurrence network edge list.

## Method & Software References

**Method:**
Friedman J, Alm EJ (2012). *Inferring Correlation Networks from Genomic Survey Data.* PLoS Computational Biology, 8(9): e1002687. https://doi.org/10.1371/journal.pcbi.1002687

**Code:**
Core SparCC functions (`sparcc`, `permute_w_replacement`) sourced from: https://github.com/bio-developer/sparcc

**Python packages:**
- `pandas`
- `numpy`

## Suggested Repository Structure

```
sparcc-network/
├── README.md
├── code/
│   └── sparcc_core.py          # from bio-developer/sparcc
├── data/
│   └── raw_counts.tsv
└── output/
    ├── rarefied_counts.tsv
    ├── cor_sparcc_final.tsv
    ├── permutations/
    │   └── permutation_0.tsv ... permutation_99.tsv
    ├── iters/
    │   └── perm_cor_0.npy ... perm_cor_99.npy
    ├── sparcc_final_results_v2.csv
    ├── sparcc_final_results_with_fdr.csv
    ├── sparcc_fdr_significant_pairs.csv
    └── sparcc_fdr_significant_pairs_r_filt.csv
```

---

## Pipeline Steps

### Step 1 — Sequencing depth QC
**Purpose:** Sum counts per sample to inspect the sequencing depth distribution and identify a defensible rarefaction depth.
**Input:** `data/raw_counts.tsv`
**Output:** console printout only (no file written)

```python
import pandas as pd

df = pd.read_csv('data/raw_counts.tsv', sep='\t', index_col=0)
totals = df.sum(axis=0)
totals.index = [int(c) for c in totals.index]
print(totals.sort_values())
```

### Step 2 — Day-0 exclusivity check
**Purpose:** Determine how many ASVs are detected only at Day 0 and would therefore be lost once the baseline timepoint is excluded.
**Input:** `data/raw_counts.tsv`
**Output:** console printout only (no file written)

```python
df_no0 = df.drop(columns=['0'])
all_zero = (df_no0.sum(axis=1) == 0).sum()
print(f'ASVs that only ever appeared at Day 0: {all_zero} of {len(df)}')
```

### Step 3 — Rarefaction depth / prevalence threshold check
**Purpose:** Simulate rarefaction at a candidate depth and test how many ASVs would pass various minimum-sample-prevalence cutoffs, to choose a filtering threshold before committing to it.
**Input:** `data/raw_counts.tsv`
**Output:** console printout only (no file written)

```python
import numpy as np

df2 = df.drop(columns=['0'])
RAREFY_DEPTH = 22677
rng = np.random.default_rng(42)

rarefied_check = pd.DataFrame(index=df2.index, columns=df2.columns, dtype=int)
for col in df2.columns:
    counts = df2[col].values.astype(np.int64)
    rarefied_check[col] = rng.multivariate_hypergeometric(counts, RAREFY_DEPTH)

n_samples_present = (rarefied_check > 0).sum(axis=1)
for min_samples in [2, 3, 4, 5, 7, 10]:
    n = (n_samples_present >= min_samples).sum()
    print(f'Present in >= {min_samples} of 13 samples: {n} ASVs')
```

### Step 4 — Rarefaction and prevalence filtering
**Purpose:** Drop Day 0, rarefy remaining samples to a common depth, and retain only ASVs present in ≥3 of 13 samples. Produces the cleaned count table used for correlation analysis.
**Input:** `data/raw_counts.tsv`
**Output:** `output/rarefied_counts.tsv`

```python
df = pd.read_csv('data/raw_counts.tsv', sep='\t', index_col=0)
df.columns = [int(c) for c in df.columns]

df_dropped = df.drop(columns=[0])

RAREFY_DEPTH = 22677
rng = np.random.default_rng(42)

rarefied = pd.DataFrame(index=df_dropped.index, columns=df_dropped.columns, dtype=int)
for col in df_dropped.columns:
    counts = df_dropped[col].values.astype(np.int64)
    sampled = rng.multivariate_hypergeometric(counts, RAREFY_DEPTH)
    rarefied[col] = sampled

n_samples_present = (rarefied > 0).sum(axis=1)
filtered = rarefied[n_samples_present >= 3].copy()

filtered.index.name = 'ASV_ID'
filtered.to_csv('output/rarefied_counts.tsv', sep='\t')
```

### Step 5 — SparCC correlation calculation
**Purpose:** Compute the true SparCC correlation matrix on the rarefied, filtered count table.
**Input:** `output/rarefied_counts.tsv`
**Output:** `output/cor_sparcc_final.tsv`

```python
import sys
sys.path.insert(0, 'code')
from sparcc_core import sparcc

df = pd.read_csv('output/rarefied_counts.tsv', sep='\t', index_col=0)
counts = df.T  # sparcc() expects rows=samples, columns=ASVs

cor_df, cov_df = sparcc(counts, iters=20, th=0.1, xiter=10, verbose=False)
cor_df.to_csv('output/cor_sparcc_final.tsv', sep='\t')
```

### Step 6 — Bootstrap permutation generation
**Purpose:** Generate 100 shuffled-with-replacement versions of the count table to build a null distribution for significance testing (100 permutations is the tool's documented default for empirical p-value estimation).
**Input:** `output/rarefied_counts.tsv`
**Output:** `output/permutations/permutation_{0-99}.tsv` (100 files)

```python
import sys
sys.path.insert(0, 'code')
from sparcc_core import permute_w_replacement

df = pd.read_csv('output/rarefied_counts.tsv', sep='\t', index_col=0)
counts = df.T

N_PERM = 100
for i in range(N_PERM):
    np.random.seed(i)
    perm = permute_w_replacement(counts)
    perm_out = perm.T
    perm_out.index.name = 'ASV_ID'
    perm_out.to_csv(f'output/permutations/permutation_{i}.tsv', sep='\t')
```

### Step 7 — SparCC correlations on permuted data
**Purpose:** Recompute SparCC correlations on each of the 100 permuted count tables, forming the empirical null distribution used for pseudo p-value estimation.
**Input:** `output/permutations/permutation_{0-99}.tsv`
**Output:** `output/iters/perm_cor_{0-99}.npy` (100 files)

```python
import os

perm_dir = 'output/permutations'
out_dir = 'output/iters'
os.makedirs(out_dir, exist_ok=True)

for iter_idx in range(100):
    out_file = f'{out_dir}/perm_cor_{iter_idx}.npy'
    if os.path.exists(out_file):
        continue

    df = pd.read_csv(f'{perm_dir}/permutation_{iter_idx}.tsv', sep='\t', index_col=0)
    counts = df.T

    cor_df, cov_df = sparcc(counts, iters=20, th=0.1, xiter=10, verbose=False)
    np.save(out_file, cor_df.values.astype(np.float32))
```

### Step 8 — Pseudo p-values, FDR correction, and final network filtering
**Purpose:** Compare real SparCC correlations against the permutation-based null distribution to estimate a pseudo p-value per ASV pair, apply Benjamini–Hochberg FDR correction, and filter to the final significant network (FDR q ≤ 0.01, |r| ≥ 0.5).
**Input:** `output/cor_sparcc_final.tsv`, `output/iters/perm_cor_{0-99}.npy`
**Output:** `output/sparcc_fdr_significant_pairs_r_filt.csv`
*(intermediate files `sparcc_final_results_v2.csv`, `sparcc_final_results_with_fdr.csv`, and `sparcc_fdr_significant_pairs.csv` are also written along the way — see summary table below)*

```python
import glob

files_dir = 'output'

cor_df = pd.read_csv(f'{files_dir}/cor_sparcc_final.tsv', sep='\t', index_col=0)
asvs = cor_df.index.tolist()
n = len(asvs)

cor_vals = cor_df.values
iu = np.triu_indices(n, k=1)  # upper triangle = unique ASV pairs
real_cor = cor_vals[iu]

perm_files = sorted(glob.glob(f'{files_dir}/iters/perm_cor_*.npy'),
                     key=lambda x: int(x.split('_')[-1].split('.')[0]))

n_perm = len(perm_files)
perm_array = np.empty((n_perm, len(real_cor)), dtype=np.float32)
for i, f in enumerate(perm_files):
    perm_array[i] = np.load(f)[iu]

# Pseudo p-values: fraction of permutations at least as extreme as the real correlation
abs_real = np.abs(real_cor)
abs_perm = np.abs(perm_array)
ge_count = (abs_perm >= abs_real[np.newaxis, :]).sum(axis=0)
pval = ge_count / n_perm

asv_index = pd.DataFrame({'index': range(n), 'ASV_ID': asvs})
idx_to_id = dict(zip(asv_index['index'], asv_index['ASV_ID']))

edge_df = pd.DataFrame({
    'ASV_1_idx': iu[0], 'ASV_2_idx': iu[1],
    'correlation': real_cor.round(4), 'pvalue': pval
})
edge_df['ASV_1'] = edge_df['ASV_1_idx'].map(idx_to_id)
edge_df['ASV_2'] = edge_df['ASV_2_idx'].map(idx_to_id)
final = edge_df[['ASV_1', 'ASV_2', 'correlation', 'pvalue']].sort_values('correlation', key=abs, ascending=False)
final.to_csv(f'{files_dir}/sparcc_final_results_v2.csv', index=False)

# Benjamini-Hochberg FDR correction
n_tests = len(final)
df_sorted = final.sort_values('pvalue').reset_index(drop=True)
ranks = np.arange(1, n_tests + 1)
bh_q = df_sorted['pvalue'] * n_tests / ranks

bh_q_adj = bh_q.values[::-1]
bh_q_adj = np.minimum.accumulate(bh_q_adj)[::-1]
bh_q_adj = np.clip(bh_q_adj, 0, 1)

df_sorted['fdr_qvalue'] = bh_q_adj
df_sorted.to_csv(f'{files_dir}/sparcc_final_results_with_fdr.csv', index=False)

# Final filtering: FDR q <= 0.01, then |r| >= 0.5
sig = df_sorted[df_sorted['fdr_qvalue'] <= 0.01].copy()
sig.to_csv(f'{files_dir}/sparcc_fdr_significant_pairs.csv', index=False)

sig_r_filt = sig[sig['correlation'].abs() >= 0.5].copy()
sig_r_filt.to_csv(f'{files_dir}/sparcc_fdr_significant_pairs_r_filt.csv', index=False)
```

---

## Overall Input/Output Summary

| Type | File | Description |
|---|---|---|
| Input | `raw_counts.tsv` | Raw ASV × sample count table (tab-delimited), all timepoints including Day 0 |
| Output | `rarefied_counts.tsv` | Rarefied (fixed read depth/sample), prevalence-filtered (≥3 of 13 samples) count table |
| Output | `cor_sparcc_final.tsv` | SparCC correlation matrix computed on real data |
| Output | `permutation_*.tsv` (×100) | Shuffled count tables used to build the null distribution |
| Output | `perm_cor_*.npy` (×100) | SparCC correlation matrices computed on each permuted table |
| Output | `sparcc_final_results_v2.csv` | All pairwise ASV correlations with pseudo p-values |
| Output | `sparcc_final_results_with_fdr.csv` | Same, with Benjamini–Hochberg FDR q-values added |
| Output | `sparcc_fdr_significant_pairs.csv` | Pairs significant at FDR q ≤ 0.01 |
| **Output (final)** | `sparcc_fdr_significant_pairs_r_filt.csv` | Final network: FDR q ≤ 0.01 and \|r\| ≥ 0.5 |

## Notes
- Random seeds are fixed (`np.random.default_rng(42)` for rarefaction; `np.random.seed(i)` per permutation index) for reproducibility.
- Rarefaction depth and the ≥3-of-13-sample prevalence threshold were chosen based on the diagnostic checks in Steps 1–3.
- All file paths above are relative and should be adapted to your local repository layout.
