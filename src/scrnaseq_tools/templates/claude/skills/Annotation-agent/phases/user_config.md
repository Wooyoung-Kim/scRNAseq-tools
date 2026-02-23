# Phase 0: User Configuration

Pre-annotation settings. Ask user BEFORE starting annotation.

---

## TF Activity Analysis (Used in Tier 2)

```
TF Activity tool/database/method:

1. Tool:
   [ ] decoupler (Recommended)
   [ ] pySCENIC
   [ ] Other: ___

2. Database:
   [ ] DoRothEA (A, B levels) 
   [ ] DoRothEA (A, B, C levels)
   [ ] CollecTRI - Recommended
   [ ] Custom: ___

3. Method:
   [ ] ULM (Univariate Linear Model) - Recommended
   [ ] MLM (Multivariate Linear Model)
   [ ] VIPER
   [ ] AUCell
```

---

## Pathway Activity Analysis (Used in Tier 3)

```
Pathway Activity tool/database/method:

1. Tool:
   [ ] decoupler (Recommended)
   [ ] GSEApy
   [ ] Other: ___

2. Database:
   [ ] PROGENy (top 500) - Recommended
   [ ] PROGENy (top 100)
   [ ] MSigDB Hallmark
   [ ] KEGG
   [ ] Reactome
   [ ] GO Biological Process
   [ ] Custom: ___

3. Method:
   [ ] mlm - Recommended for decoupler
   [ ] GSEA prerank
   [ ] ORA
   [ ] ssGSEA
   [ ] GSVA
```

---

## Trajectory Analysis (Used in Tier 2)

```
Trajectory tool/settings:

1. Tool:
   [ ] Palantir (Recommended)
   [ ] DPT (Diffusion Pseudotime)
   [ ] Monocle3
   [ ] Slingshot
   [ ] Skip trajectory

2. Root cell:
   [ ] Automatic (early progenitor markers)
   [ ] Manual: ___
```

---

## Subset Re-clustering Preprocessing (Tier 2, 3)

```
When subsetting data for Tier 2/3, how to preprocess before re-clustering:

1. Re-clustering Strategy:
   [ ] Option A: Use existing Harmony reduction (Recommended)
       - Keep batch correction from full data
       - Only re-compute neighbors within subset
       - Fast, consistent with full dataset

   [ ] Option B: New HVG → PCA → No integration
       - Subset-specific HVGs
       - Fresh PCA on subset
       - No batch correction (risk of batch effects)

   [ ] Option C: New HVG → PCA → Harmony → Re-clustering
       - Full re-processing
       - Fresh batch correction on subset
       - Most thorough but may over-correct

2. HVG Settings (if Option B or C):
   Flavor:
   [ ] seurat_v3 (Recommended for UMI data)
   [ ] seurat
   [ ] cell_ranger

   n_top_genes: [3000] (default)

   batch_key: [batch] (for batch-aware HVG selection)

3. Integration Settings (if Option C):
   [ ] Harmony (Recommended)
   [ ] scVI
   [ ] BBKNN

   batch_key: [batch]

4. Neighbors Settings:
   n_neighbors: [15] (default)
   n_pcs: [50] (default, or None to use all PCs)

5. Resolution Strategy:
   [ ] Fixed resolution (user-specified per tier)
   [ ] Multi-resolution scan (Recommended)
       - Run multiple resolutions, evaluate, select best

   If Fixed:
     Tier 2 resolution: [0.8] (default)
     Tier 3 resolution: [1.2] (default)

   If Multi-resolution:
     Resolution range: [0.4, 0.6, 0.8, 1.0, 1.2, 1.5] (default)
     Selection metric:
     [ ] Silhouette score (Recommended)
     [ ] Modularity
     [ ] Manual inspection
```

**Guidance by Subset Size:**
```
| Subset Size | Batch Count | Recommended |
|-------------|-------------|-------------|
| < 5K cells  | Any         | Option A    |
| 5K-20K      | 1-2         | Option A    |
| 5K-20K      | 3+          | Option A/C  |
| > 20K       | Any         | Option B/C  |
```

---

## Default Configuration

If user accepts defaults or skips:

```python
DEFAULT_CONFIG = {
    'tf_activity': {
        'tool': 'decoupler',
        'database': 'collectri',
        'method': 'ulm',
        'split_complexes': False
    },
    'pathway_activity': {
        'tool': 'decoupler',
        'database': 'progeny',
        'method': 'mlm',
        'top_n': 500
    },
    'trajectory': {
        'tool': 'palantir',
        'use': True,
        'root_cell': 'automatic'
    },
    'reclustering': {
        'strategy': 'A',  # 'A', 'B', or 'C'
        'use_rep': 'X_pca_harmony',  # For Option A
        'hvg': {
            'flavor': 'seurat_v3',
            'n_top_genes': 3000,
            'batch_key': 'batch'
        },
        'integration': {
            'tool': 'harmony',
            'batch_key': 'batch'
        },
        'neighbors': {
            'n_neighbors': 15,
            'n_pcs': 50
        },
        'resolution': {
            'strategy': 'multi',  # 'fixed' or 'multi'
            'fixed': {
                'tier2': 0.8,
                'tier3': 1.2
            },
            'multi': {
                'range': [0.4, 0.6, 0.8, 1.0, 1.2, 1.5],
                'metric': 'silhouette'  # 'silhouette', 'modularity', 'manual'
            }
        }
    }
}
```

---

## Configuration Code

```python
def get_analysis_config():
    """
    Store user's analysis configuration.
    Use AskUserQuestion tool to gather preferences.
    """
    config = {
        'tf_activity': {
            'tool': None,
            'database': None,
            'method': None,
            'levels': None
        },
        'pathway_activity': {
            'tool': None,
            'database': None,
            'method': None,
            'top_n': None
        },
        'trajectory': {
            'tool': None,
            'use': True,
            'root_cell': None
        },
        'reclustering': {
            'strategy': None,  # 'A', 'B', or 'C'
            'use_rep': None,   # e.g., 'X_pca_harmony'
            'hvg': {
                'flavor': None,
                'n_top_genes': None,
                'batch_key': None
            },
            'integration': {
                'tool': None,
                'batch_key': None
            }
        }
    }
    return config


def preprocess_subset(subset, config, tier):
    """
    Preprocess subset according to user's re-clustering strategy.
    Call this BEFORE re-clustering in Tier 2/3.
    """
    import scanpy as sc

    strategy = config['reclustering']['strategy']
    neighbors_config = config['reclustering']['neighbors']
    n_neighbors = neighbors_config['n_neighbors']
    n_pcs = neighbors_config['n_pcs']

    if strategy == 'A':
        # Option A: Use existing reduction, re-compute neighbors only
        use_rep = config['reclustering']['use_rep']
        print(f"📊 Re-clustering strategy A: Using existing {use_rep}")
        sc.pp.neighbors(subset, use_rep=use_rep,
                        n_neighbors=n_neighbors, n_pcs=n_pcs)

    elif strategy == 'B':
        # Option B: New HVG → PCA → No integration
        hvg_config = config['reclustering']['hvg']
        print(f"📊 Re-clustering strategy B: New HVG → PCA (no integration)")

        sc.pp.highly_variable_genes(
            subset,
            flavor=hvg_config['flavor'],
            n_top_genes=hvg_config['n_top_genes'],
            batch_key=hvg_config.get('batch_key')
        )
        sc.pp.pca(subset, use_highly_variable=True)
        sc.pp.neighbors(subset, use_rep='X_pca',
                        n_neighbors=n_neighbors, n_pcs=n_pcs)

    elif strategy == 'C':
        # Option C: New HVG → PCA → Integration → Re-clustering
        hvg_config = config['reclustering']['hvg']
        int_config = config['reclustering']['integration']
        print(f"📊 Re-clustering strategy C: New HVG → PCA → {int_config['tool']}")

        sc.pp.highly_variable_genes(
            subset,
            flavor=hvg_config['flavor'],
            n_top_genes=hvg_config['n_top_genes'],
            batch_key=hvg_config.get('batch_key')
        )
        sc.pp.pca(subset, use_highly_variable=True)

        # Integration
        if int_config['tool'] == 'harmony':
            import scanpy.external as sce
            sce.pp.harmony_integrate(subset, key=int_config['batch_key'])
            sc.pp.neighbors(subset, use_rep='X_pca_harmony',
                            n_neighbors=n_neighbors, n_pcs=n_pcs)
        elif int_config['tool'] == 'scvi':
            # scVI integration code
            pass
        elif int_config['tool'] == 'bbknn':
            import bbknn
            bbknn.bbknn(subset, batch_key=int_config['batch_key'],
                        n_pcs=n_pcs, neighbors_within_batch=n_neighbors)

    else:
        raise ValueError(f"Unknown re-clustering strategy: {strategy}")

    return subset


def select_resolution(subset, config, tier, cluster_key_prefix='leiden'):
    """
    Select optimal resolution using fixed or multi-resolution strategy.

    Returns:
        tuple: (best_resolution, cluster_key)
    """
    import scanpy as sc
    from sklearn.metrics import silhouette_score
    import numpy as np

    res_config = config['reclustering']['resolution']
    strategy = res_config['strategy']

    if strategy == 'fixed':
        # Fixed resolution
        resolution = res_config['fixed'][f'tier{tier}']
        cluster_key = f'{cluster_key_prefix}_res{resolution}'
        sc.tl.leiden(subset, resolution=resolution, key_added=cluster_key)
        print(f"📊 Fixed resolution: {resolution} → {subset.obs[cluster_key].nunique()} clusters")
        return resolution, cluster_key

    elif strategy == 'multi':
        # Multi-resolution scan
        resolutions = res_config['multi']['range']
        metric = res_config['multi']['metric']

        print(f"📊 Multi-resolution scan: {resolutions}")
        results = []

        for res in resolutions:
            key = f'{cluster_key_prefix}_res{res}'
            sc.tl.leiden(subset, resolution=res, key_added=key)
            n_clusters = subset.obs[key].nunique()

            if metric == 'silhouette':
                # Skip if only 1 cluster
                if n_clusters <= 1:
                    score = -1
                else:
                    # Use PCA embedding for silhouette
                    use_rep = config['reclustering'].get('use_rep', 'X_pca')
                    if use_rep in subset.obsm:
                        X = subset.obsm[use_rep][:, :50]
                    else:
                        X = subset.obsm['X_pca'][:, :50]
                    labels = subset.obs[key].astype('category').cat.codes
                    score = silhouette_score(X, labels, sample_size=min(5000, X.shape[0]))

            elif metric == 'modularity':
                # Modularity is stored after leiden
                score = subset.uns.get('leiden', {}).get('modularity', 0)

            else:
                score = 0  # Manual selection

            results.append({
                'resolution': res,
                'key': key,
                'n_clusters': n_clusters,
                'score': score
            })
            print(f"  res={res}: {n_clusters} clusters, {metric}={score:.3f}")

        # Select best
        if metric == 'manual':
            print("\n⚠️ Manual selection required. Review clusters and choose resolution.")
            return None, None

        best = max(results, key=lambda x: x['score'])
        print(f"\n✅ Best resolution: {best['resolution']} ({metric}={best['score']:.3f})")

        # Clean up other resolution columns
        for r in results:
            if r['key'] != best['key']:
                del subset.obs[r['key']]

        return best['resolution'], best['key']

    else:
        raise ValueError(f"Unknown resolution strategy: {strategy}")
```

---

## Checklist

- [ ] TF Activity configured
- [ ] Pathway Activity configured
- [ ] Trajectory configured
- [ ] Re-clustering strategy configured (A/B/C)
- [ ] HVG settings confirmed (if B or C): flavor=seurat_v3, n_top_genes=3000
- [ ] Integration settings confirmed (if C)
- [ ] Neighbors settings confirmed: n_neighbors=15, n_pcs=50
- [ ] Resolution strategy configured (fixed/multi)
- [ ] Resolution range confirmed (if multi): [0.4, 0.6, 0.8, 1.0, 1.2, 1.5]
- [ ] Selection metric confirmed (if multi): silhouette/modularity/manual
- [ ] User confirmed all settings
