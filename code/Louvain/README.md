# Positive Co-occurrence Network — Louvain Module Assignment

This README documents the Louvain community (module) detection performed on
the positive-correlation SparCC co-occurrence network, and lists the files
associated with that calculation.

## Input data

| File | Description |
|---|---|
| `sparcc_fdr_significant_pairs_r03.csv` | SparCC pairwise correlations between ASVs that passed FDR-significance filtering. Columns: `ASV_1`, `ASV_2`, `correlation`, `pvalue`, `fdr_qvalue`. 2,218 total significant pairs. |

## Filtering step

The significant-pairs table contains both positive (co-occurrence) and
negative (mutual exclusion) correlations:

- **1,269 positive pairs** (correlation > 0)
- 949 negative pairs (correlation < 0) — **excluded from this calculation**

Only the positive pairs were used to build the network for module detection.
Standard modularity maximization (and the Louvain algorithm built around it)
assumes non-negative edge weights — the null model term at its core compares
observed edge density to expected density under random rewiring, which is not
well-defined when weights can be negative. This is also standard convention
in 16S/microbiome co-occurrence network literature (e.g. SCNIC and similar
tools default to Louvain on a positive-only or magnitude-thresholded
network). Negative correlations were analyzed separately (see
"Related analyses," below) rather than folded into this module assignment.

## Network construction

- **Nodes**: 266 unique ASVs (every ASV appearing in at least one positive
  significant pair)
- **Edges**: 1,269, undirected, weighted by the SparCC correlation
  coefficient (all edge weights strictly > 0 by construction)
- **Connectivity**: 2 connected components — one giant component of 264
  nodes, and one isolated pair of 2 nodes

## Community detection method

- **Algorithm**: Louvain method (Blondel et al., 2008), via
  `networkx.algorithms.community.louvain_communities`
  (networkx 3.6.1's built-in implementation — no external `python-louvain`
  package was available/needed)
- **Parameters**: `weight='weight'` (SparCC correlation as edge weight),
  `resolution=1.0` (default), `seed=42` (fixed for reproducibility — Louvain
  has some randomness in tie-breaking/node visit order, so results can vary
  slightly run-to-run without a fixed seed)

## Results

**11 communities (modules) detected. Modularity Q = 0.6065.**

A modularity above ~0.4 is generally considered evidence of real community
structure, so 0.6065 indicates the positive co-occurrence network partitions
cleanly into modules.

| Module | Size (# ASVs) |
|---|---|
| 0 | 63 |
| 1 | 57 |
| 2 | 52 |
| 3 | 39 |
| 4 | 25 |
| 5 | 8 |
| 6 | 8 |
| 7 | 7 |
| 8 | 3 |
| 9 | 2 |
| 10 | 2 |

Modules 8–10 (sizes 3, 2, 2) are small, largely disconnected fragments rather
than substantial community structure, and may warrant separate treatment or
exclusion in downstream interpretation.

## Output files

| File | Description |
|---|---|
| `asv_community_assignments.csv` | One row per ASV: `ASV` (hash ID), `iterative_ID` (human-readable taxonomic label where available), `community` (module 0–10), `degree` (unweighted degree in the positive network) |

## Related analyses (built on top of this module assignment, not part of it)

These used the same 11-module assignment as an input, but are separate
calculations:

- **z-P (Guimerà & Amaral, 2005) topological role classification** — within-
  module degree z-score and participation coefficient P computed per ASV,
  classifying each node's role (hub, connector, peripheral, etc.) relative to
  these modules.
- **Signed Louvain (exploratory)** — a from-scratch implementation of signed
  modularity (Gómez, Jensen & Arenas, 2009) that includes both positive and
  negative edges in a single joint optimization. This was explored but not
  adopted, per the standard practice in 16S/microbiome literature of keeping
  positive (co-occurrence) and negative (exclusion) analyses separate rather
  than combining them into one clustering.
- **Exclusion-module analysis** — Louvain run separately on the negative-only
  edge subnetwork (weight = |correlation|), plus a module × module matrix
  quantifying mutual exclusion pressure between the 11 positive modules
  found here. This surfaced a strong mutual-exclusion relationship between
  Module 0 and Module 3 specifically.

## Reproducing this calculation

```python
import pandas as pd
import networkx as nx

df = pd.read_csv("sparcc_fdr_significant_pairs_r03.csv")
pos = df[df["correlation"] > 0]

G = nx.Graph()
for _, row in pos.iterrows():
    G.add_edge(row["ASV_1"], row["ASV_2"], weight=row["correlation"])

communities = nx.algorithms.community.louvain_communities(
    G, weight="weight", resolution=1.0, seed=42
)
modularity = nx.algorithms.community.modularity(G, communities, weight="weight")
```
