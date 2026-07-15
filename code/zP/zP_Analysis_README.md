# z-P Topological Role Analysis — README

This README documents the within-module degree / participation coefficient
(z-P) analysis performed on the positive co-occurrence network and its 11
Louvain modules, following Guimerà & Amaral (2005), "Functional cartography
of complex metabolic networks."

## Purpose

Louvain gives each ASV a module assignment, but says nothing about *how* an
ASV participates in its module — whether it's a densely-connected core
member, a peripheral member, or a bridge connecting multiple modules. The
z-P framework classifies each node's topological role using two metrics
computed relative to the module structure.

## Input

- The positive co-occurrence network (266 ASVs, 1,269 edges) and its 11
  Louvain modules — see the companion Louvain module README/summary report
  for how these were built.

## Metrics

### Within-module degree z-score (z)

For each ASV *i* in module *s*:

```
z_i = (k_i,s − mean(k_s)) / std(k_s)
```

- `k_i,s` = ASV *i*'s **unweighted** degree counting only edges to other
  ASVs in its own module *s* (i.e., how many co-occurrence partners it has
  inside its own module).
- `mean(k_s)` / `std(k_s)` = the mean and population standard deviation of
  that intra-module degree, computed across all ASVs in module *s*.

This measures whether an ASV is unusually well-connected *within its own
module*, relative to other members of that same module. Degree here is
unweighted (edge counts, not correlation-weighted), matching Guimerà &
Amaral's original convention.

### Participation coefficient (P)

For each ASV *i*:

```
P_i = 1 − Σ_s (k_i,s / k_i)²
```

- `k_i` = ASV *i*'s total (unweighted) degree across the whole network.
- `k_i,s` = ASV *i*'s degree restricted to module *s* (summed over all
  modules *s*, including its own).

P is 0 if all of an ASV's connections are inside a single module, and
approaches 1 if its connections are spread evenly across many modules. High
P identifies ASVs that bridge multiple communities rather than sitting
inside just one.

## Role classification thresholds

Each ASV is assigned to one of seven roles based on (z, P):

| Role | Condition |
|---|---|
| Provincial hub | z ≥ 2.5, P ≤ 0.30 |
| Connector hub | z ≥ 2.5, 0.30 < P ≤ 0.75 |
| Kinless hub | z ≥ 2.5, P > 0.75 |
| Ultra-peripheral non-hub | z < 2.5, P ≤ 0.05 |
| Peripheral non-hub | z < 2.5, 0.05 < P ≤ 0.62 |
| Non-hub connector | z < 2.5, 0.62 < P ≤ 0.80 |
| Non-hub kinless | z < 2.5, P > 0.80 |

The z = 2.5 threshold separates "hub" roles (unusually well-connected within
their module) from "non-hub" roles; the P thresholds further separate each
group by how confined vs. spread-out their connections are.

## Results

| Role | Count |
|---|---|
| Ultra-peripheral non-hub | 151 |
| Peripheral non-hub | 111 |
| Non-hub connector | 3 |
| Provincial hub | 1 |
| Connector hub | 0 |
| Kinless hub | 0 |
| Non-hub kinless | 0 |

This distribution follows directly from the network's high modularity
(Q = 0.6065): most ASVs have their connections concentrated within their own
module (low P), so almost everything falls into the two "non-hub, low-P"
categories. Only one ASV qualifies as any kind of hub, and only 3 as
non-hub connectors — these 4 ASVs are the network's only topologically
"unusual" members.

**Notably, all 4 of these role-important ASVs have very low relative
abundance** (all under 0.02%) — a common and often-cited finding in
co-occurrence network analysis: topological importance and abundance are
frequently decoupled, and rare taxa can still play structurally central
roles.

## Output

`zP_analysis.csv` — one row per ASV:

| Column | Description |
|---|---|
| `ASV` | ASV hash ID |
| `iterative_ID` | Human-readable taxonomic label, where available |
| `community` | Louvain module (0–10) |
| `degree` | Total unweighted degree in the positive network |
| `intramodule_degree` | Unweighted degree within its own module (k_i,s) |
| `z_score` | Within-module degree z-score |
| `P_score` | Participation coefficient |
| `role` | One of the 7 categories above |

## Reproducing this calculation

```python
import numpy as np

modules = sorted(set(node_comm.values()))

# intra-module degree
k_i_module = {n: sum(1 for nb in G.neighbors(n) if node_comm[nb] == node_comm[n])
              for n in G.nodes()}

# z-score, computed per module
z_score = {}
for s in modules:
    members = [n for n in G.nodes() if node_comm[n] == s]
    vals = np.array([k_i_module[n] for n in members])
    mean_s, std_s = vals.mean(), vals.std()
    for n in members:
        z_score[n] = (k_i_module[n] - mean_s) / std_s if std_s > 0 else 0.0

# participation coefficient
P_score = {}
for n in G.nodes():
    k_i = G.degree(n)
    if k_i == 0:
        P_score[n] = 0.0
        continue
    s_counts = {}
    for nb in G.neighbors(n):
        s_counts[node_comm[nb]] = s_counts.get(node_comm[nb], 0) + 1
    P_score[n] = 1 - sum((c / k_i) ** 2 for c in s_counts.values())
```

## Limitations

- Degree is unweighted throughout (matching the original Guimerà-Amaral
  convention), so a node with many weak correlations and a node with few
  very strong correlations could receive the same z/P scores. A
  weighted-degree variant was not explored here.
- With only 1 hub and 3 connectors total, sample size for these categories
  is very small — treat any downstream interpretation of "hub behavior" as
  suggestive rather than statistically robust.
