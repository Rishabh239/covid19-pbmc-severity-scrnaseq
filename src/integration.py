"""
Integration Module for scRNA-seq Analysis Pipeline
==================================================

Multi-sample integration, clustering, and PBMC cell-type annotation.

When several samples are analyzed together, cells otherwise cluster by sample
(and therefore by condition) rather than by cell type - Lee et al. (2020) saw
exactly this in GSE149689. We correct that with batch integration so downstream
clusters are cell types shared across healthy / mild / severe.

Memory: the heavy steps (ComBat / scale / PCA / integration / clustering) run on
a highly-variable-gene subset (~n_top_genes columns), so peak memory stays small
even for tens of thousands of cells. The full-gene raw counts are preserved on
the returned object (``layers['counts']``) for pseudobulk DE, and log-normalized
values in ``adata.raw`` for marker scoring.

Integration methods (``method=``): 'harmony' (default; calls harmonypy directly),
'combat' (built into scanpy, no compiled deps), or 'bbknn' (needs the bbknn
package). If harmonypy can't be imported, 'harmony' falls back to 'combat'.

scanpy is imported lazily inside functions so the pure-logic helpers remain
importable/testable without a full single-cell stack.

Author: scRNA Analysis Pipeline
License: MIT
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


# Canonical PBMC lineage markers (coarse, transparent auto-annotation).
DEFAULT_PBMC_MARKERS: Dict[str, List[str]] = {
    "CD4 T": ["IL7R", "CD3D", "CD3E", "CD4", "TCF7", "CCR7"],
    "CD8 T": ["CD8A", "CD8B", "CD3D", "GZMK", "CCL5"],
    "NK": ["GNLY", "NKG7", "KLRD1", "KLRF1", "NCAM1"],
    "B": ["MS4A1", "CD79A", "CD79B", "CD19"],
    "Plasma": ["MZB1", "XBP1", "JCHAIN"],
    "CD14 Mono": ["CD14", "LYZ", "S100A8", "S100A9", "VCAN"],
    "CD16 Mono": ["FCGR3A", "MS4A7", "CDKN1C"],
    "DC": ["FCER1A", "CST3", "CLEC10A"],
    "pDC": ["LILRA4", "GZMB", "IL3RA"],
    "Platelet": ["PPBP", "PF4", "ITGA2B"],
}


# ---------------------------------------------------------------------------
# Pure-logic helper (unit-testable without scanpy)
# ---------------------------------------------------------------------------
def assign_labels_from_scores(score_df: pd.DataFrame) -> Dict[str, str]:
    """Map each cluster to its highest-scoring lineage.

    ``score_df``: rows = cluster ids, columns = lineage names, values = mean
    signature score. Returns ``{cluster_id: lineage}``.
    """
    if score_df.shape[1] == 0:
        raise ValueError("score_df has no lineage columns to choose from.")
    return {str(cl): str(score_df.loc[cl].astype(float).idxmax()) for cl in score_df.index}


# ---------------------------------------------------------------------------
# AnnData-facing functions
# ---------------------------------------------------------------------------
def integrate_samples(
    adata,
    batch_key: str = "sample_id",
    method: str = "harmony",
    n_top_genes: int = 2000,
    n_pcs: int = 30,
    n_neighbors: int = 15,
    resolution: float = 1.0,
    min_genes: int = 200,
    max_counts: int = 25000,
    max_mt_pct: float = 10.0,
    results_dir: str = "results",
    run_qc: bool = True,
):
    """QC -> normalize -> batch-aware HVG -> (heavy steps on HVG subset) -> cluster.

    Expects a labeled, raw-count AnnData (e.g. from ``cohort.load_cohort``).
    Returns the full-gene object with ``layers['counts']`` (raw), ``adata.raw``
    (log-norm), ``obs['leiden']`` and ``obsm['X_umap']`` populated.
    """
    import scanpy as sc

    Path(results_dir).mkdir(parents=True, exist_ok=True)

    # Preserve raw counts (sparse) for pseudobulk before any transformation.
    if "counts" not in adata.layers:
        adata.layers["counts"] = adata.X.copy()

    # QC (reuse the pipeline's existing preprocessor for consistency + plots).
    if run_qc:
        try:
            from src.preprocessing import ScRNAPreprocessor
        except Exception:  # pragma: no cover
            from preprocessing import ScRNAPreprocessor
        pre = ScRNAPreprocessor(
            qc_params={"min_genes": min_genes, "max_counts": max_counts,
                       "max_mt_pct": max_mt_pct},
            results_dir=str(Path(results_dir) / "qc"),
        )
        adata = pre.run_qc(adata)
    if "counts" not in adata.layers:
        adata.layers["counts"] = adata.X.copy()

    logger.info("Normalizing (per-cell CP10K + log1p)...")
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.raw = adata  # log-normalized, all genes -> used by DE/marker scoring

    logger.info("Selecting %d batch-aware highly variable genes...", n_top_genes)
    sc.pp.highly_variable_genes(
        adata, n_top_genes=n_top_genes, flavor="seurat", batch_key=batch_key
    )

    method = (method or "harmony").lower()
    if method == "harmony":
        try:
            import harmonypy  # noqa: F401
        except Exception as e:
            logger.warning(
                "harmonypy not importable (%s). Falling back to method='combat' "
                "(no compiled deps). For Harmony: pip install 'harmonypy<2'.", e,
            )
            method = "combat"

    # --- Heavy steps on an HVG-subset copy (keeps memory small) --------------
    # The full-gene raw counts + log-norm stay on `adata`; only cluster labels
    # and embeddings are transferred back afterwards.
    work = adata[:, adata.var["highly_variable"]].copy()
    work.raw = None  # drop the redundant full-gene raw copy from the work object
    logger.info(
        "Integration working set: %d cells x %d HVGs (method=%s)",
        work.n_obs, work.n_vars, method,
    )

    if method == "combat":
        logger.info("ComBat batch correction over '%s'...", batch_key)
        sc.pp.combat(work, key=batch_key)

    sc.pp.scale(work, max_value=10)
    sc.tl.pca(work, n_comps=n_pcs, svd_solver="arpack")

    use_rep = "X_pca"
    if method == "harmony":
        try:
            import harmonypy
            logger.info("Harmony integration over '%s' (direct harmonypy call)...", batch_key)
            ho = harmonypy.run_harmony(work.obsm["X_pca"], work.obs, [batch_key])
            Z = np.asarray(ho.Z_corr)
            if Z.ndim == 2 and Z.shape[0] == work.n_obs:
                emb = Z
            elif Z.ndim == 2 and Z.shape[1] == work.n_obs:
                emb = Z.T
            else:
                raise ValueError(
                    f"unexpected Harmony Z_corr shape {Z.shape} for {work.n_obs} cells"
                )
            work.obsm["X_pca_harmony"] = np.ascontiguousarray(emb)
            use_rep = "X_pca_harmony"
        except Exception as e:
            raise RuntimeError(
                f"Harmony step failed ({e}). Re-run with --integration combat "
                "(built into scanpy, needs no harmonypy)."
            ) from e

    logger.info("Neighbors / UMAP / Leiden (method=%s)...", method)
    if method == "bbknn":
        import scanpy.external as sce
        sce.pp.bbknn(work, batch_key=batch_key, n_pcs=n_pcs)
    else:
        sc.pp.neighbors(work, n_neighbors=n_neighbors, n_pcs=n_pcs, use_rep=use_rep)
    sc.tl.umap(work)
    sc.tl.leiden(work, resolution=resolution, key_added="leiden")
    logger.info("Found %d clusters.", work.obs["leiden"].nunique())

    # --- Transfer results back to the full-gene object -----------------------
    adata.obs["leiden"] = pd.Categorical(work.obs["leiden"].values)
    adata.obsm["X_umap"] = work.obsm["X_umap"]
    adata.obsm["X_pca"] = work.obsm["X_pca"]
    if "X_pca_harmony" in work.obsm:
        adata.obsm["X_pca_harmony"] = work.obsm["X_pca_harmony"]
    del work

    # Diagnostic UMAPs (batch mixing is the thing to eyeball).
    for color in ("leiden", "condition", batch_key):
        if color in adata.obs:
            try:
                sc.pl.umap(adata, color=color, save=f"_{color}.png", show=False)
            except Exception as e:  # pragma: no cover
                logger.warning("Could not plot UMAP(%s): %s", color, e)

    return adata


def annotate_pbmc(
    adata,
    markers: Dict[str, List[str]] = DEFAULT_PBMC_MARKERS,
    groupby: str = "leiden",
    key_added: str = "cell_type",
    use_raw: bool = True,
) -> Tuple[object, Dict[str, str], pd.DataFrame]:
    """Coarsely label clusters by scoring canonical PBMC signatures.

    Returns ``(adata, cluster->label mapping, per-cluster score table)``.
    Treat the labels as a first pass and refine against the marker plots.
    Marker scoring uses ``adata.raw`` (all genes), so markers that aren't highly
    variable are still available.
    """
    import scanpy as sc

    ref_vars = (
        adata.raw.var_names if (use_raw and adata.raw is not None) else adata.var_names
    )
    present = {ct: [g for g in gs if g in ref_vars] for ct, gs in markers.items()}
    present = {ct: gs for ct, gs in present.items() if gs}
    if not present:
        raise ValueError("None of the marker genes are present in the dataset.")

    score_cols = []
    for ct, gs in present.items():
        col = f"_score_{ct}"
        sc.tl.score_genes(adata, gene_list=gs, score_name=col, use_raw=use_raw)
        score_cols.append(col)

    per_cluster = adata.obs.groupby(groupby, observed=True)[score_cols].mean()
    per_cluster.columns = [c.replace("_score_", "") for c in per_cluster.columns]
    mapping = assign_labels_from_scores(per_cluster)

    adata.obs[key_added] = pd.Categorical(
        adata.obs[groupby].astype(str).map(mapping).astype(str)
    )
    adata.obs.drop(columns=score_cols, inplace=True, errors="ignore")

    logger.info("Cluster -> cell type: %s", mapping)
    for ct, n in adata.obs[key_added].value_counts().items():
        logger.info("  %s: %d cells", ct, n)

    try:
        sc.pl.umap(adata, color=key_added, save="_cell_type.png", show=False)
    except Exception as e:  # pragma: no cover
        logger.warning("Could not plot cell-type UMAP: %s", e)

    return adata, mapping, per_cluster


if __name__ == "__main__":
    df = pd.DataFrame(
        {"CD14 Mono": [2.0, 0.1], "CD8 T": [0.2, 1.5], "B": [0.0, 0.3]},
        index=["0", "1"],
    )
    assert assign_labels_from_scores(df) == {"0": "CD14 Mono", "1": "CD8 T"}
    print("integration.py self-test: label assignment OK.")