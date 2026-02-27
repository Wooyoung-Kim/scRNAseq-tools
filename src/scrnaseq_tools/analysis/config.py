"""
Configuration classes for scrna_qc pipeline.
"""

from dataclasses import dataclass, field
from typing import Optional, List, Dict, Any
from pathlib import Path
import yaml


@dataclass
class DoubletConfig:
    """Configuration for doublet detection."""

    enabled: bool = True
    method: str = "scrublet"  # Currently only "scrublet" supported
    expected_doublet_rate: Optional[float] = None  # Auto-estimated if None
    filter_doublets: bool = True  # Remove doublets after detection
    n_prin_comps: int = 30
    min_counts: int = 2
    min_cells: int = 3
    min_gene_variability_pctl: float = 85.0
    n_neighbors: Optional[int] = None  # Auto-selected if None
    threshold: Optional[float] = None  # Auto-selected if None


@dataclass
class AgentConfig:
    """Configuration for AI Agent mode."""

    enabled: bool = True
    auto_detect_batch_key: bool = True
    auto_detect_mt_pattern: bool = True
    auto_select_parameters: bool = True
    generate_structured_report: bool = True
    report_formats: List[str] = field(default_factory=lambda: ["json", "markdown"])
    include_next_steps: bool = True
    verbose_logging: bool = False


@dataclass
class IntegrationConfig:
    """Configuration for integration methods."""
    
    # Which methods to run
    methods: List[str] = field(default_factory=lambda: ["scvi"])
    # Options: ["scvi"], ["harmony"], ["scvi", "harmony"], ["all"]
    
    # scVI parameters
    scvi_n_latent: int = 30
    scvi_n_layers: int = 2
    scvi_max_epochs: int = 200
    scvi_early_stopping: bool = True
    scvi_gene_likelihood: str = "nb"  # "nb" or "zinb"
    
    # Harmony parameters
    harmony_n_pcs: int = 50
    harmony_max_iter: int = 30
    
    # Scanorama parameters (if used)
    scanorama_knn: int = 20
    
    # Benchmarking
    auto_select_best: bool = True
    benchmark_label_key: Optional[str] = None  # For bio conservation metrics
    
    def get_methods_list(self) -> List[str]:
        """Expand 'all' to full method list."""
        if "all" in self.methods:
            return ["scvi", "harmony", "scanorama"]
        return self.methods


@dataclass
class MultiResolutionConfig:
    """Configuration for multi-resolution clustering (for annotation agent)."""

    enabled: bool = False  # Enable multi-resolution clustering
    resolutions: List[float] = field(default_factory=lambda: [0.1, 0.3, 0.5])
    # low: broad cell types, mid: sub-types, high: fine clusters
    auto_select: bool = True  # Auto-select resolutions based on data size
    min_cells_per_cluster: int = 50  # Minimum cells per cluster

    # Output column names
    cluster_key_prefix: str = "cluster_res"  # e.g., cluster_res_0.1, cluster_res_0.3

    def get_resolution_keys(self) -> List[str]:
        """Get cluster column names for each resolution."""
        return [f"{self.cluster_key_prefix}_{r}" for r in self.resolutions]


@dataclass
class PipelineConfig:
    """Configuration for the scRNA-seq QC to Integration pipeline."""

    # ========== QC parameters ==========
    adaptive_qc: bool = True  # Use adaptive MAD-based thresholds
    qc_method: str = "mad"  # "mad" or "quantile"
    qc_nmads: float = 3.0  # Number of MADs for outlier detection
    min_genes: int = 200  # Floor for min_genes
    min_cells: int = 3  # Minimum cells per gene
    max_mt_pct: float = 20.0  # Ceiling for max MT%

    # ========== Normalization ==========
    target_sum: float = 1e4  # Size factor target (counts per 10k)
    n_top_genes: int = 2000  # Number of HVGs
    hvg_flavor: str = "seurat_v3"  # HVG selection method

    # ========== Dimensionality reduction ==========
    n_pcs: int = 50  # Maximum PCs to compute
    n_neighbors: int = 15  # kNN neighbors
    pc_selection_method: str = "elbow"  # "elbow" or "variance_ratio"
    variance_threshold: float = 0.90  # For variance_ratio method

    # ========== Clustering ==========
    resolution: float = 0.5
    clustering_method: str = "leiden"  # "leiden" or "louvain"
    multi_resolution: MultiResolutionConfig = field(default_factory=MultiResolutionConfig)
    
    # ========== Integration ==========
    batch_key: Optional[str] = None  # Column containing batch information
    integration: IntegrationConfig = field(default_factory=IntegrationConfig)

    # ========== Doublet Detection ==========
    doublet: DoubletConfig = field(default_factory=DoubletConfig)

    # ========== Agent Mode ==========
    agent: AgentConfig = field(default_factory=AgentConfig)

    # ========== Reproducibility ==========
    seed: int = 42
    
    # ========== Output ==========
    save_intermediate: bool = True  # Save intermediate h5ad files
    generate_plots: bool = True  # Generate QC and integration plots
    
    @classmethod
    def from_yaml(cls, path: str) -> "PipelineConfig":
        """Load configuration from YAML file."""
        with open(path, 'r') as f:
            config_dict = yaml.safe_load(f) or {}

        # Handle nested IntegrationConfig
        integration_dict = config_dict.pop('integration', {})
        integration_config = IntegrationConfig(**integration_dict) if integration_dict else IntegrationConfig()

        # Handle nested DoubletConfig
        doublet_dict = config_dict.pop('doublet', {})
        doublet_config = DoubletConfig(**doublet_dict) if doublet_dict else DoubletConfig()

        # Handle nested AgentConfig
        agent_dict = config_dict.pop('agent', {})
        agent_config = AgentConfig(**agent_dict) if agent_dict else AgentConfig()

        # Filter to only valid fields
        valid_fields = {f.name for f in cls.__dataclass_fields__.values()}
        filtered = {k: v for k, v in config_dict.items() if k in valid_fields}
        filtered['integration'] = integration_config
        filtered['doublet'] = doublet_config
        filtered['agent'] = agent_config

        return cls(**filtered)
    
    def to_yaml(self, path: str):
        """Save configuration to YAML file."""
        config_dict = {
            'adaptive_qc': self.adaptive_qc,
            'qc_method': self.qc_method,
            'qc_nmads': self.qc_nmads,
            'min_genes': self.min_genes,
            'min_cells': self.min_cells,
            'max_mt_pct': self.max_mt_pct,
            'target_sum': self.target_sum,
            'n_top_genes': self.n_top_genes,
            'hvg_flavor': self.hvg_flavor,
            'n_pcs': self.n_pcs,
            'n_neighbors': self.n_neighbors,
            'pc_selection_method': self.pc_selection_method,
            'variance_threshold': self.variance_threshold,
            'resolution': self.resolution,
            'clustering_method': self.clustering_method,
            'batch_key': self.batch_key,
            'seed': self.seed,
            'save_intermediate': self.save_intermediate,
            'generate_plots': self.generate_plots,
            'integration': {
                'methods': self.integration.methods,
                'scvi_n_latent': self.integration.scvi_n_latent,
                'scvi_n_layers': self.integration.scvi_n_layers,
                'scvi_max_epochs': self.integration.scvi_max_epochs,
                'scvi_early_stopping': self.integration.scvi_early_stopping,
                'scvi_gene_likelihood': self.integration.scvi_gene_likelihood,
                'harmony_n_pcs': self.integration.harmony_n_pcs,
                'harmony_max_iter': self.integration.harmony_max_iter,
                'auto_select_best': self.integration.auto_select_best,
                'benchmark_label_key': self.integration.benchmark_label_key,
            },
            'doublet': {
                'enabled': self.doublet.enabled,
                'method': self.doublet.method,
                'expected_doublet_rate': self.doublet.expected_doublet_rate,
                'filter_doublets': self.doublet.filter_doublets,
                'n_prin_comps': self.doublet.n_prin_comps,
            },
            'agent': {
                'enabled': self.agent.enabled,
                'auto_detect_batch_key': self.agent.auto_detect_batch_key,
                'auto_detect_mt_pattern': self.agent.auto_detect_mt_pattern,
                'auto_select_parameters': self.agent.auto_select_parameters,
                'generate_structured_report': self.agent.generate_structured_report,
                'report_formats': self.agent.report_formats,
            }
        }

        with open(path, 'w') as f:
            yaml.dump(config_dict, f, default_flow_style=False, sort_keys=False)
    
    def with_integration_methods(self, methods: List[str]) -> "PipelineConfig":
        """Return a copy with updated integration methods."""
        from copy import deepcopy
        new_config = deepcopy(self)
        new_config.integration.methods = methods
        return new_config
