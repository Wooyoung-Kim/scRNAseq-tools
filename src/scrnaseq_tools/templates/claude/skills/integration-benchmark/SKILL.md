# Integration Benchmark Skill

> `scib-metrics` 기반 batch integration 벤치마킹 가이드

## 개요

Integration 방법들(scVI, Harmony, Scanorama, scANVI 등)의 성능을 **정량적으로** 비교하는 스킬.
두 축으로 평가:
- **Bio conservation** — 생물학적 신호(cell type 분리)가 보존되었는가
- **Batch correction** — 배치 효과가 제거되었는가

---

## 필수 패키지

```bash
pip install scib-metrics scanpy scvi-tools harmonypy scanorama
```

---

## 워크플로우

### Step 1: 데이터 전처리

```python
import scanpy as sc

adata = sc.read_h5ad("input.h5ad")

# raw counts를 위한 layer 확인
assert "counts" in adata.layers, "counts layer 필요"

# HVG 선택 (batch-aware)
adata.raw = adata  # 전체 차원 보존
sc.pp.highly_variable_genes(
    adata,
    flavor="seurat_v3",
    n_top_genes=2000,
    layer="counts",
    batch_key="batch",
    subset=True,
)

# PCA (Unintegrated baseline)
sc.tl.pca(adata, n_comps=30, use_highly_variable=True)
adata.obsm["Unintegrated"] = adata.obsm["X_pca"]
```

> [!IMPORTANT]
> **HVG는 반드시 batch_key를 지정하여 선택**. 배치별 유전자 발현의 편향을 줄입니다.

### Step 2: Integration 방법 실행

각 방법의 latent representation을 `adata.obsm`에 저장합니다.

#### Harmony (PCA 기반, 가장 빠름)

```python
from harmony import harmonize

adata.obsm["Harmony"] = harmonize(
    adata.obsm["X_pca"],
    adata.obs,
    batch_key="batch",
)
```

#### scVI (VAE 기반)

```python
import scvi

scvi.model.SCVI.setup_anndata(adata, layer="counts", batch_key="batch")
vae = scvi.model.SCVI(adata, n_layers=2, n_latent=30, gene_likelihood="nb")
vae.train(max_epochs=300)
adata.obsm["scVI"] = vae.get_latent_representation()
```

#### scANVI (semi-supervised, label 필요)

```python
lvae = scvi.model.SCANVI.from_scvi_model(
    vae,
    adata=adata,
    labels_key="cell_type",  # 알려진 cell type label
    unlabeled_category="Unknown",
)
lvae.train(max_epochs=100, n_samples_per_label=100)
adata.obsm["scANVI"] = lvae.get_latent_representation()
```

#### Scanorama

```python
import numpy as np
import scanorama

batch_cats = adata.obs.batch.cat.categories
adata_list = [adata[adata.obs.batch == b].copy() for b in batch_cats]
scanorama.integrate_scanpy(adata_list)

adata.obsm["Scanorama"] = np.zeros((adata.shape[0], adata_list[0].obsm["X_scanorama"].shape[1]))
for i, b in enumerate(batch_cats):
    adata.obsm["Scanorama"][adata.obs.batch == b] = adata_list[i].obsm["X_scanorama"]
```

### Step 3: Benchmarking

```python
from scib_metrics.benchmark import Benchmarker, BioConservation, BatchCorrection

bm = Benchmarker(
    adata,
    batch_key="batch",
    label_key="cell_type",  # ground truth cell type (필수)
    embedding_obsm_keys=["Unintegrated", "Harmony", "scVI", "scANVI", "Scanorama"],
    bio_conservation_metrics=BioConservation(),
    batch_correction_metrics=BatchCorrection(),
    n_jobs=-1,  # 모든 코어 사용
)
bm.benchmark()
```

### Step 4: 결과 시각화

```python
# 결과 테이블 (scaled)
fig = bm.plot_results_table()
fig.savefig("benchmark_results_scaled.png", dpi=150, bbox_inches="tight")

# 결과 테이블 (raw scores)
fig = bm.plot_results_table(min_max_scale=False)
fig.savefig("benchmark_results_raw.png", dpi=150, bbox_inches="tight")

# DataFrame으로 접근
df = bm.get_results(min_max_scale=False)
df.to_csv("benchmark_results.csv")
print(df.transpose())
```

### Step 5: UMAP 비교 시각화

```python
import matplotlib.pyplot as plt

methods = ["Unintegrated", "Harmony", "scVI", "scANVI", "Scanorama"]
fig, axes = plt.subplots(2, len(methods), figsize=(5 * len(methods), 10))

for i, method in enumerate(methods):
    sc.pp.neighbors(adata, use_rep=method, key_added=f"neighbors_{method}")
    sc.tl.umap(adata, neighbors_key=f"neighbors_{method}")
    adata.obsm[f"X_umap_{method}"] = adata.obsm["X_umap"].copy()

    # Row 1: cell type
    sc.pl.embedding(adata, basis=f"X_umap_{method}", color="cell_type",
                    ax=axes[0, i], show=False, title=f"{method}\n(cell type)",
                    legend_loc="none")
    # Row 2: batch
    sc.pl.embedding(adata, basis=f"X_umap_{method}", color="batch",
                    ax=axes[1, i], show=False, title=f"{method}\n(batch)",
                    legend_loc="none")

plt.tight_layout()
plt.savefig("umap_comparison.png", dpi=150, bbox_inches="tight")
```

---

## 평가 메트릭 해석

### Bio Conservation (높을수록 좋음)

| 메트릭 | 의미 | 해석 |
|--------|------|------|
| **Isolated labels** | 희귀 cell type이 잘 분리되는지 | KMeans F1 on isolated labels |
| **NMI** | 클러스터-실제 label 일치도 | Normalized Mutual Information |
| **ARI** | 클러스터-실제 label 일치도 | Adjusted Rand Index (chance 보정) |
| **Silhouette (label)** | 같은 type은 가까이, 다른 type은 멀리 | -1 ~ 1, 높을수록 좋음 |
| **cLISI** | Cell type별 local mixing | 1에 가까우면 cell type 잘 분리 |

### Batch Correction (높을수록 좋음)

| 메트릭 | 의미 | 해석 |
|--------|------|------|
| **Silhouette (batch)** | 서로 다른 배치가 잘 섞였는지 | 0에 가까우면 잘 섞임 |
| **iLISI** | 배치별 local mixing | 높을수록 배치가 잘 섞임 |
| **KBET** | 배치 비율 보존 | kNN에서 배치 비율이 global과 유사 |
| **Graph connectivity** | 같은 label의 세포들이 연결되는지 | kNN 그래프에서 연결성 |
| **PCR comparison** | PCA 분산에서 배치 효과 제거량 | integration 전후 比較 |

### 종합 결정 기준

```
Total Score = 0.6 × Bio Conservation + 0.4 × Batch Correction
```

> [!TIP]
> - **Bio conservation > Batch correction** — 생물학적 신호 보존이 더 중요
> - scANVI가 일반적으로 **Total Score 최고** (semi-supervised)
> - Harmony가 **속도 대비 성능 최고** (PCA 기반, GPU 불필요)
> - label 정보가 없으면 scVI > Harmony > Scanorama 순서로 시도

---

## 방법 선택 의사결정 트리

```
Label(cell_type) 있는가?
├── YES → scANVI (from_scvi_model 활용, linear_classifier=True 권장)
└── NO
    ├── GPU 있는가?
    │   ├── YES → scVI (n_layers=2, n_latent=30, gene_likelihood="nb")
    │   └── NO
    │       ├── 배치 수 ≤ 10 → Harmony (빠르고 안정적)
    │       └── 배치 수 > 10 → Scanorama (많은 배치에 강함)
    └── 항상 → Unintegrated를 baseline으로 포함하여 비교
```

---

## 체크포인트 저장

```python
# Integration 결과 + 벤치마크 점수 저장
adata.write("integrated_benchmark.h5ad")
df.to_csv("benchmark_scores.csv")
```

---

## References

- [scib-metrics Lung Benchmark](https://scib-metrics.readthedocs.io/en/stable/notebooks/lung_example.html)
- [scVI-tools scanvi_fix Benchmark](https://docs.scvi-tools.org/en/latest/tutorials/notebooks/scrna/scanvi_fix.html)
- [Luecken et al., 2022 (scIB paper)](https://doi.org/10.1038/s41592-021-01336-8)
