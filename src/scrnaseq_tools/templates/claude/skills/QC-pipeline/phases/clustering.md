# Clustering & Resolution Selection

## 개요

Leiden 클러스터링의 resolution 파라미터는 annotation 품질에 직접적인 영향을 미칩니다.
**단일 resolution 대신 multi-resolution 전략**을 사용하여 annotation의 각 tier에 맞는 해상도를 제공합니다.

---

## Multi-resolution 전략 (권장)

```python
import scanpy as sc
from sklearn.metrics import silhouette_score

# 여러 resolution에서 클러스터링
resolutions = [0.3, 0.5, 0.8, 1.0, 1.2, 1.5]
scores = {}

for res in resolutions:
    key = f"leiden_{res}"
    sc.tl.leiden(adata, resolution=res, key_added=key)
    n_clusters = adata.obs[key].nunique()

    # Silhouette score로 품질 평가
    if n_clusters > 1:
        scores[res] = silhouette_score(
            adata.obsm["X_pca_harmony"],  # integration 후 embedding 사용
            adata.obs[key],
            sample_size=min(5000, adata.n_obs),
            random_state=42,
        )
    print(f"  res={res}: {n_clusters} clusters, silhouette={scores.get(res, 'N/A'):.3f}")

best_res = max(scores, key=scores.get)
print(f"\nBest resolution: {best_res} (silhouette={scores[best_res]:.3f})")
```

---

## Resolution 선택 기준

### Silhouette Score

- **값 범위**: -1 ~ 1 (높을수록 좋음)
- **의미**: 같은 클러스터 내 거리 vs 클러스터 간 거리
- **주의**: resolution이 너무 높으면 silhouette가 올라가지만 **over-clustering** 위험

### 데이터 규모별 권장 범위

| 세포 수 | 권장 resolution | 예상 클러스터 수 |
|---------|----------------|----------------|
| < 5,000 | 0.3 – 0.8 | 5 – 15 |
| 5,000 – 30,000 | 0.5 – 1.2 | 10 – 30 |
| > 30,000 | 0.8 – 2.0 | 20 – 50 |

### Over-clustering 방지

```python
# 최소 세포 수 제약
MIN_CELLS_PER_CLUSTER = 50

for res in resolutions:
    key = f"leiden_{res}"
    sizes = adata.obs[key].value_counts()
    small_clusters = (sizes < MIN_CELLS_PER_CLUSTER).sum()
    if small_clusters > 0:
        print(f"⚠️ res={res}: {small_clusters} clusters with <{MIN_CELLS_PER_CLUSTER} cells")
```

---

## Annotation 연결: 3-tier Resolution

Annotation-agent의 3-round 구조에 맞게 3단계 resolution을 준비합니다:

```python
# Tier 1: 대분류 (broad cell types)
sc.tl.leiden(adata, resolution=0.3, key_added="cluster_res_0.3")

# Tier 2: 중분류 (sub-types)
sc.tl.leiden(adata, resolution=0.5, key_added="cluster_res_0.5")

# Tier 3: 소분류 (fine-grained)
sc.tl.leiden(adata, resolution=best_res, key_added=f"cluster_res_{best_res}")

# 기본 cluster 열은 Tier 1으로 설정
adata.obs["cluster"] = adata.obs["cluster_res_0.3"]
```

> **Phase 0 config에서 `resolution.strategy = "multi"`로 설정하면 자동으로 이 과정이 실행됩니다.**

---

## 저장 및 체크포인트

```python
# 체크포인트 저장
adata.write("qc_output/clustered.h5ad")

# 클러스터 정보 기록
for res in resolutions:
    key = f"leiden_{res}"
    if key in adata.obs:
        n = adata.obs[key].nunique()
        print(f"  {key}: {n} clusters")
```
