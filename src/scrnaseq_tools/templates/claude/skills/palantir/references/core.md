# Palantir Core Functions (palantir.core)

## run_palantir

**Main Palantir trajectory inference.**

```python
pr_res = palantir.core.run_palantir(
    data,                    # DataFrame or AnnData
    early_cell,              # Start cell barcode
    terminal_states=None,    # Terminal states (pd.Series or None)
    knn=30,                  # K-nearest neighbors
    num_waypoints=1200,      # Number of waypoints
    n_jobs=-1,               # Parallel jobs
    scale_components=True,   # Scale diffusion components
    use_early_cell_as_start=False,  # Use early cell as pseudotime 0
    max_iterations=25,       # Convergence iterations
    eigvec_key='DM_EigenVectors_multiscaled',
    pseudo_time_key='palantir_pseudotime',
    entropy_key='palantir_entropy',
    fate_prob_key='palantir_fate_probabilities',
    save_as_df=None,         # Save as DataFrame vs array
    waypoints_key='palantir_waypoints',
    seed=20
)

# For AnnData, stores results in:
#   - adata.obs['palantir_pseudotime']
#   - adata.obs['palantir_entropy']
#   - adata.obsm['palantir_fate_probabilities']
#   - adata.uns['palantir_waypoints']

# Returns PResults object (or None if AnnData)
```

### Parameters Explained

| Parameter | Description | Default |
|-----------|-------------|---------|
| `early_cell` | Cell barcode to start trajectory from | Required |
| `terminal_states` | pd.Series mapping barcodes to names (None for auto-detect) | None |
| `knn` | Neighbors for Markov chain construction | 30 |
| `num_waypoints` | Waypoints for efficient computation (more = accurate but slower) | 1200 |
| `scale_components` | Scale diffusion components by eigenvalues | True |
| `max_iterations` | Iterations for pseudotime convergence | 25 |
| `use_early_cell_as_start` | Force early cell to have pseudotime 0 | False |

### Terminal States Format

```python
import pandas as pd

# Manual specification
terminal_states = pd.Series(
    ["Erythroid", "Monocyte", "DC"],  # Terminal state names (values)
    index=["cell_ery_001", "cell_mono_001", "cell_dc_001"]  # Cell barcodes (index)
)

# Or let Palantir auto-detect (set to None)
terminal_states = None
```

### Example

```python
import palantir
import pandas as pd

# Define terminal states
terminal_states = pd.Series(
    ['Erythroid', 'Monocyte', 'DC'],
    index=['cell_ery', 'cell_mono', 'cell_dc']
)

# Run Palantir
palantir.core.run_palantir(
    adata,
    early_cell='cell_hsc_001',
    terminal_states=terminal_states,
    num_waypoints=500,
    knn=30
)

# Access results
pseudotime = adata.obs['palantir_pseudotime']
entropy = adata.obs['palantir_entropy']
fate_probs = adata.obsm['palantir_fate_probabilities']
```

---

## identify_terminal_states

**Automatically identify terminal states.**

```python
terminal_cells, boundary_cells = palantir.core.identify_terminal_states(
    ms_data,                 # Multiscale diffusion components
    early_cell,              # Start cell
    knn=30,
    num_waypoints=1200,
    n_jobs=-1,
    max_iterations=25,
    seed=20
)

# Returns:
#   - terminal_cells: Array of terminal state cell indices
#   - boundary_cells: Index of excluded boundary cells
```

---

## PResults Class

**Container for Palantir results (when not using AnnData).**

```python
# Properties
pr_res.pseudotime        # pd.Series of pseudotime values
pr_res.entropy           # pd.Series of entropy values
pr_res.branch_probs      # pd.DataFrame of fate probabilities
pr_res.waypoints         # Array of waypoint indices

# Methods
pr_res.save("results.pkl")  # Save to pickle
pr_res = PResults.load("results.pkl")  # Load from pickle
```

---

## Interpretation Guide

### Pseudotime
- Values from 0 (start) to 1 (terminal)
- Represents ordering along differentiation, not real time
- Cells with similar pseudotime are at similar differentiation stages

### Entropy
- High entropy (~1): Multipotent, can differentiate into multiple fates
- Low entropy (~0): Committed to a specific lineage
- Useful for identifying decision points (branching regions)

### Fate Probabilities
- Matrix of cells × terminal states
- Each row sums to 1
- High probability = likely to differentiate into that lineage
- Balanced probabilities = cell at decision point

---

## Algorithm Overview

1. **Waypoint sampling**: Sample representative cells for efficiency
2. **Pseudotime initialization**: Shortest path from start cell
3. **Markov chain construction**: Build transition probabilities using KNN
4. **Iterative refinement**: Converge pseudotime ordering
5. **Terminal state detection**: Find absorbing states (if not specified)
6. **Fate probability computation**: Absorption probabilities in Markov chain
7. **Entropy calculation**: Shannon entropy of fate probabilities
