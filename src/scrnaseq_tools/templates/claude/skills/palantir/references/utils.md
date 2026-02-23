# Palantir Utility Functions (palantir.utils)

## Diffusion Maps

### run_diffusion_maps

**Compute diffusion maps (core of Palantir).**

```python
dm_res = palantir.utils.run_diffusion_maps(
    data,                    # DataFrame or AnnData (PCA required)
    n_components=10,         # Number of diffusion components
    knn=30,                  # K-nearest neighbors
    alpha=0,                 # Diffusion operator normalization
    seed=0,                  # Random seed
    kernel_backend='scanpy', # Backend: 'scanpy' or 'sklearn'
    pca_key='X_pca',         # Input PCA key
    kernel_key='DM_Kernel',  # Output kernel key
    sim_key='DM_Similarity', # Output similarity key
    eigval_key='DM_EigenValues',    # Output eigenvalue key
    eigvec_key='DM_EigenVectors'    # Output eigenvector key
)

# For AnnData, stores in:
#   - adata.obsm['DM_EigenVectors']
#   - adata.uns['DM_EigenValues']
#   - adata.obsp['DM_Kernel']
```

### determine_multiscale_space

**Determine optimal diffusion component scaling.**

```python
ms_data = palantir.utils.determine_multiscale_space(
    dm_res,                  # Dict or AnnData with diffusion results
    n_eigs=None,             # Number of eigenvalues (auto if None)
    eigval_key='DM_EigenValues',
    eigvec_key='DM_EigenVectors',
    out_key='DM_EigenVectors_multiscaled'
)

# For AnnData, stores in adata.obsm['DM_EigenVectors_multiscaled']
```

### compute_kernel

**Compute diffusion kernel only.**

```python
kernel = palantir.utils.compute_kernel(
    data,                    # DataFrame or AnnData
    knn=30,                  # K-nearest neighbors
    alpha=0,                 # Normalization parameter
    pca_key='X_pca',         # PCA key
    kernel_key='DM_Kernel',  # Output key
    backend=None             # Backend ('scanpy' or 'sklearn')
)
```

---

## Cell Finding Functions

### early_cell

**Find early/progenitor cell of a given type.**

```python
start_cell = palantir.utils.early_cell(
    ad,                      # AnnData
    celltype='HSC',          # Cell type to find
    celltype_column='celltype',  # Column with cell types
    eigvec_key='DM_EigenVectors_multiscaled',
    fallback_seed=None       # Seed for random fallback
)

# Returns cell barcode (str)
# Finds cell of given type at diffusion map extreme
```

### find_terminal_states

**Find terminal states for multiple cell types.**

```python
terminal_states = palantir.utils.find_terminal_states(
    ad,                      # AnnData
    celltypes=['Ery', 'Mono', 'DC'],  # Terminal cell types
    celltype_column='celltype',
    eigvec_key='DM_EigenVectors_multiscaled',
    fallback_seed=None
)

# Returns pd.Series with cell type names as values, barcodes as index
```

### fallback_terminal_cell

**Backup method for finding terminal cells.**

```python
terminal_cell = palantir.utils.fallback_terminal_cell(
    ad,                      # AnnData
    celltype='Erythrocyte',  # Cell type
    celltype_column='anno',  # Cell type column
    eigvec_key='DM_EigenVectors_multiscaled',
    seed=2353                # Random seed
)
```

---

## MAGIC Imputation

### run_magic_imputation

**MAGIC imputation for denoised gene expression.**

```python
imputed = palantir.utils.run_magic_imputation(
    data,                    # Array, DataFrame, or AnnData
    dm_res=None,             # Diffusion results (optional if AnnData)
    n_steps=3,               # Diffusion steps
    sim_key='DM_Similarity', # Similarity matrix key
    expression_key=None,     # Input expression key
    imputation_key='MAGIC_imputed_data',  # Output key
    n_jobs=-1,               # Parallel jobs
    sparse=True,             # Return sparse matrix
    clip_threshold=0.01      # Clip extreme values
)

# For AnnData, stores in adata.layers['MAGIC_imputed_data']
```

**Example:**
```python
# Run MAGIC
palantir.utils.run_magic_imputation(adata)

# Visualize imputed expression
sc.pl.umap(adata, layer='MAGIC_imputed_data', color=['CD34', 'GATA1', 'MPO'])
```

---

## Density and Variability

### run_density

**Estimate cell-state density using Mellon.**

```python
density = palantir.utils.run_density(
    ad,                      # AnnData
    repr_key='DM_EigenVectors',
    density_key='mellon_log_density',
    **kwargs                 # Additional mellon arguments
)

# Stores in adata.obs['mellon_log_density']
```

### run_local_variability

**Compute local gene expression variability.**

```python
local_var = palantir.utils.run_local_variability(
    ad,                      # AnnData
    expression_key='MAGIC_imputed_data',
    distances_key='distances',
    localvar_key='local_variability',
    progress=False,
    eps=1e-16
)

# Stores in adata.varm['local_variability']
```

### run_low_density_variability

**Compute gene variability in low-density regions.**

```python
scores = palantir.utils.run_low_density_variability(
    ad,                      # AnnData
    cell_mask='branch_masks',
    density_key='mellon_log_density',
    localvar_key='local_variability',
    score_key='low_density_gene_variability'
)
```

---

## PCA

### run_pca

**Principal component analysis.**

```python
pca_projections, pca_components = palantir.utils.run_pca(
    data,                    # DataFrame or AnnData
    n_components=300,        # Number of PCs
    use_hvg=True,            # Use highly variable genes only
    pca_key='X_pca'          # Key for storing in AnnData
)

# For AnnData, stores in adata.obsm['X_pca']
```

---

## Preprocessing (palantir.preprocess)

### log_transform

```python
palantir.preprocess.log_transform(
    data,                    # Counts matrix or AnnData
    pseudo_count=0.1         # Pseudocount before log
)

# For AnnData, modifies adata.X in place
# Formula: log2(data + pseudo_count)
```

### filter_counts_data

```python
filtered_data = palantir.preprocess.filter_counts_data(
    data,                    # Counts matrix
    cell_min_molecules=1000, # Minimum molecules per cell
    genes_min_cells=10       # Minimum cells per gene
)
```

### normalize_counts

```python
normalized = palantir.preprocess.normalize_counts(
    data                     # Counts matrix
)
```
