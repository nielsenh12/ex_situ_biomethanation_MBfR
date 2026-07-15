# Positive Co-occurrence Network — Louvain Module Assignment

This README documents the Louvain community (module) detection performed on
the positive-correlation SparCC co-occurrence network, and lists the files
associated with that calculation.

## Input data

| File | Description |
|---|---|
| `sparcc_fdr_significant_pairs_r_filt.csv` | SparCC pairwise correlations between ASVs that passed FDR-significance filtering. Columns: `ASV_1`, `ASV_2`, `correlation`, `pvalue`, `fdr_qvalue`. 1,577 total significant pairs. |

## Filtering step

The significant-pairs table contains both positive (co-occurrence) and
negative (mutual exclusion) correlations:

- **854 positive pairs** (correlation > 0)
- 723 negative pairs (correlation < 0) — **excluded from this calculation**

Only the positive pairs were used to build the network for module detection.
Standard modularity maximization (and the Louvain algorithm built around it)
assumes non-negative edge weights — the null model term at its core compares
observed edge density to expected density under random rewiring, which is not
well-defined when weights can be negative. This is also standard convention
in 16S/microbiome co-occurrence network literature (e.g. SCNIC and similar
tools default to Louvain on a positive-only or magnitude-thresholded
network). Negative correlations were not folded into this module assignment.

## Network construction

- **Nodes**: 165 unique ASVs (every ASV appearing in at least one positive
  significant pair)
- **Edges**: 854, undirected, weighted by the SparCC correlation
  coefficient (all edge weights strictly > 0 by construction)
- **Connectivity**: 2 connected components — one giant component of 163
  nodes, and one isolated pair of 2 nodes

## Community detection method

- **Algorithm**: Louvain method (Blondel et al., 2008), via
  `networkx.algorithms.community.louvain_communities`
  (networkx 3.6.1's built-in implementation — no external `python-louvain`
  package was available/needed)
- **Parameters**: `weight='weight'` (SparCC correlation as edge weight),
  default resolution (`resolution=1.0`), `seed=42` (fixed for
  reproducibility — Louvain has some randomness in tie-breaking/node visit
  order, so results can vary slightly run-to-run without a fixed seed)

## Results

**9 communities (modules) detected. Modularity Q = 0.5798.**

A modularity above ~0.4 is generally considered evidence of real community
structure, so 0.5798 indicates the positive co-occurrence network partitions
cleanly into modules.

| Module | Size (# ASVs) |
|---|---|
| 0 | 4 |
| 1 | 31 |
| 2 | 20 |
| 3 | 50 |
| 4 | 5 |
| 5 | 9 |
| 6 | 42 |
| 7 | 2 |
| 8 | 2 |

The three smallest modules (0, 7, 8) are not all the same kind of "small,"
and only one is a true disconnected fragment:

| Module | Size | Internal edges | External edges | Notes |
|---|---|---|---|---|
| 7 | 2 | 1 | 0 | Genuinely disconnected — this is the isolated 2-node component (0 edges to the rest of the network). |
| 0 | 4 | 4 | 4 | Embedded in the giant component but weakly cohesive as a module (more external edges than internal). |
| 8 | 2 | 1 | 1 | Embedded in the giant component; not disconnected. |

Modules 4 (size 5, 5 internal / 5 external) and 5 (size 9, 15 internal / 26
external) are also loosely bounded — their external edge counts are
comparable to or exceed their internal edge counts — so they are similarly
weak as standalone modules, despite not being among the three smallest by
size. Modules 1, 2, 3, and 6 are all solidly cohesive (internal edges
substantially outnumber external edges). Only module 7 should be described
as a disconnected fragment; modules 0, 4, 5, and 8 warrant a "weakly-defined
module" caveat rather than a "disconnected" one.

## Output files

| File | Description |
|---|---|
| `asv_community_assignments.csv` | One row per ASV: `ASV` (hash ID), `community` (module 0–8), `degree` (unweighted degree in the positive network). |

## Related analyses (built on top of this module assignment, not part of it)

- **z-P (Guimerà & Amaral, 2005) topological role classification** — within-
  module degree z-score and participation coefficient P computed per ASV,
  classifying each node's role (hub, connector, peripheral, etc.) relative to
  these modules. See `zP_analysis.csv`.

## Reproducing this calculation

```python
import pandas as pd
import networkx as nx

# Paths anonymized — replace with your local input/output directories
input_dir = "<path-to-input-directory>"
output_dir = "<path-to-output-directory>"

# Load the final filtered network (adjust filename/threshold to whatever you settled on)
df = pd.read_csv(f"{input_dir}/sparcc_fdr_significant_pairs_r_filt.csv")

# Positive correlations only
df_positive = df[df["correlation"] > 0].copy()

# Build the graph
G = nx.Graph()
for _, row in df_positive.iterrows():
    G.add_edge(row["ASV_1"], row["ASV_2"], weight=row["correlation"])

print(f"Graph: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")

# Run Louvain community detection
communities = nx.community.louvain_communities(G, weight="weight", seed=42)
print(f"Number of communities found: {len(communities)}")

# Compute modularity score (higher = more well-defined community structure, typically 0.3+ is considered meaningful)
modularity = nx.community.modularity(G, communities, weight="weight")
print(f"Modularity: {modularity:.4f}")

# Build a results table: one row per ASV, with its assigned community and degree
asv_to_community = {}
for comm_id, comm_members in enumerate(communities):
    for asv in comm_members:
        asv_to_community[asv] = comm_id

degree_dict = dict(G.degree())

results = pd.DataFrame({
    "ASV": list(asv_to_community.keys()),
    "community": list(asv_to_community.values()),
})
results["degree"] = results["ASV"].map(degree_dict)
results = results.sort_values(["community", "degree"], ascending=[True, False])

results.to_csv(f"{output_dir}/asv_community_assignments.csv", index=False)
print("Saved asv_community_assignments.csv")
print()
print("Community sizes:")
print(results["community"].value_counts().sort_index())
```

**Console output from this run:**

```
Graph: 165 nodes, 854 edges
Number of communities found: 9
Modularity: 0.5798
Saved asv_community_assignments.csv

Community sizes:
community
0     4
1    31
2    20
3    50
4     5
5     9
6    42
7     2
8     2
Name: count, dtype: int64
```
