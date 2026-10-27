"""
Differential State (Pseudobulk DE) Module
=========================================

The core comparative step: severe-vs-mild differential expression *within each
cell type*, done the rigorous way.

Rather than testing per cell (which treats thousands of cells from one patient
as independent replicates and inflates significance), we aggregate raw counts to
one pseudobulk profile per (sample x cell type) and run DESeq2 across samples.
Sample - not cell - is the unit of replication.

``aggregate_counts`` is pure numpy/scipy/pandas and unit-tested here. The DESeq2
step uses pydeseq2 and is exercised on your machine.

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


# ---------------------------------------------------------------------------
# Pure-logic core (unit-testable)
# ---------------------------------------------------------------------------
def aggregate_counts(
    X,
    obs: pd.DataFrame,
    gene_names,
    sample_key: str = "sample_id",
    celltype_key: str = "cell_type",
    condition_key: str = "condition",
    min_cells: int = 10,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Sum raw counts per (sample x cell type) pseudobulk.

    ``X`` is a cells x genes matrix (dense or scipy sparse). Returns
    ``(counts_df [genes x pseudobulk], meta_df [pseudobulk -> sample/celltype/
    condition/n_cells])``. Groups with fewer than ``min_cells`` cells are
    dropped.
    """
    from scipy.sparse import issparse

    obs = obs.reset_index(drop=True)
    genes = list(gene_names)
    data: Dict[str, np.ndarray] = {}
    meta_rows = []

    grouped = obs.groupby([sample_key, celltype_key], observed=True).indices
    for (samp, ct), idx in grouped.items():
        if len(idx) < min_cells:
            continue
        sub = X[idx]
        vec = np.asarray(sub.sum(axis=0)).ravel() if issparse(sub) else np.asarray(sub).sum(axis=0).ravel()
        pbid = f"{samp}|{ct}"
        data[pbid] = vec
        meta_rows.append(
            {
                "pseudobulk_id": pbid,
                sample_key: samp,
                celltype_key: ct,
                condition_key: obs[condition_key].iloc[idx[0]] if condition_key in obs else "NA",
                "n_cells": int(len(idx)),
            }
        )

    if not data:
        raise ValueError(
            f"No (sample, cell type) group reached min_cells={min_cells}. "
            "Lower min_cells or check your labels."
        )

    counts_df = pd.DataFrame(data, index=genes)
    meta_df = pd.DataFrame(meta_rows).set_index("pseudobulk_id")
    counts_df = counts_df[meta_df.index]  # align column order
    return counts_df, meta_df


def top_up_genes(
    res: pd.DataFrame,
    n: int = 200,
    padj_thresh: float = 0.05,
    lfc_thresh: float = 0.0,
) -> List[str]:
    """Significant up-regulated genes from a DESeq2 result (index = gene)."""
    if res is None or res.empty:
        return []
    sig = res[(res["padj"] < padj_thresh) & (res["log2FoldChange"] > lfc_thresh)]
    sig = sig.sort_values("log2FoldChange", ascending=False)
    return sig.head(n).index.astype(str).tolist()


# ---------------------------------------------------------------------------
# AnnData / pydeseq2 functions
# ---------------------------------------------------------------------------
def make_pseudobulk(
    adata,
    sample_key: str = "sample_id",
    celltype_key: str = "cell_type",
    condition_key: str = "condition",
    layer: str = "counts",
    min_cells: int = 10,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Pseudobulk raw counts from an AnnData (wrapper around aggregate_counts)."""
    X = adata.layers[layer] if layer in adata.layers else adata.X
    return aggregate_counts(
        X, adata.obs, adata.var_names,
        sample_key=sample_key, celltype_key=celltype_key,
        condition_key=condition_key, min_cells=min_cells,
    )


def run_deseq2_for_celltype(
    counts_df: pd.DataFrame,
    meta_df: pd.DataFrame,
    celltype: str,
    celltype_key: str = "cell_type",
    condition_key: str = "condition",
    group_test: str = "severe",
    group_ref: str = "mild",
    min_total_count: int = 10,
    min_samples: int = 2,
) -> Optional[pd.DataFrame]:
    """DESeq2 (pydeseq2) for one cell type: ``group_test`` vs ``group_ref``."""
    m = meta_df[
        (meta_df[celltype_key] == celltype)
        & (meta_df[condition_key].isin([group_test, group_ref]))
    ]
    per = m.groupby(condition_key, observed=True).size()
    if per.get(group_test, 0) < min_samples or per.get(group_ref, 0) < min_samples:
        logger.warning(
            "Skipping '%s': need >=%d samples per group (have %s).",
            celltype, min_samples, dict(per),
        )
        return None

    counts = counts_df[m.index].T.round().astype(int)          # samples x genes
    counts = counts.loc[:, counts.sum(axis=0) >= min_total_count]
    metadata = m[[condition_key]].copy()
    metadata[condition_key] = pd.Categorical(
        metadata[condition_key], categories=[group_ref, group_test]
    )

    try:
        from pydeseq2.dds import DeseqDataSet
        from pydeseq2.ds import DeseqStats
    except ImportError:
        logger.error("pydeseq2 not installed. Install with: pip install pydeseq2")
        return None

    try:
        dds = DeseqDataSet(counts=counts, metadata=metadata, design_factors=condition_key)
    except TypeError:  # newer pydeseq2 API
        dds = DeseqDataSet(counts=counts, metadata=metadata, design=f"~{condition_key}")
    dds.deseq2()
    stat = DeseqStats(dds, contrast=[condition_key, group_test, group_ref])
    stat.summary()
    res = stat.results_df.sort_values("padj", na_position="last")
    logger.info(
        "%s: %d genes with padj<0.05", celltype, int((res["padj"] < 0.05).sum())
    )
    return res


def run_all_celltypes(
    adata,
    group_test: str = "severe",
    group_ref: str = "mild",
    celltype_key: str = "cell_type",
    condition_key: str = "condition",
    sample_key: str = "sample_id",
    min_cells: int = 10,
    outdir: str = "results/pseudobulk",
) -> Tuple[Dict[str, pd.DataFrame], pd.DataFrame, pd.DataFrame]:
    """Pseudobulk + DESeq2 across every cell type. Saves one CSV per cell type."""
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)

    counts_df, meta_df = make_pseudobulk(
        adata, sample_key=sample_key, celltype_key=celltype_key,
        condition_key=condition_key, min_cells=min_cells,
    )
    counts_df.to_csv(out / "pseudobulk_counts.csv")
    meta_df.to_csv(out / "pseudobulk_metadata.csv")

    results: Dict[str, pd.DataFrame] = {}
    for ct in sorted(meta_df[celltype_key].astype(str).unique()):
        res = run_deseq2_for_celltype(
            counts_df, meta_df, ct,
            celltype_key=celltype_key, condition_key=condition_key,
            group_test=group_test, group_ref=group_ref,
        )
        if res is not None:
            safe = str(ct).replace(" ", "_").replace("/", "-")
            res.to_csv(out / f"deseq2_{safe}_{group_test}_vs_{group_ref}.csv")
            results[ct] = res
    return results, counts_df, meta_df


def enrich_up_genes(
    gene_list: List[str],
    results_dir: str = "results/enrichment",
    organism: str = "human",
) -> Dict[str, Optional[pd.DataFrame]]:
    """Run GO + KEGG enrichment on an up-gene list via the existing analyzer.

    Note: gseapy's Enrichr backend needs internet access.
    """
    if not gene_list:
        logger.warning("Empty gene list; skipping enrichment.")
        return {"go": None, "kegg": None}
    try:
        from src.enrichment import EnrichmentAnalyzer
    except Exception:  # pragma: no cover
        from enrichment import EnrichmentAnalyzer
    an = EnrichmentAnalyzer(results_dir=results_dir, organism=organism)
    go = an.run_go_enrichment(gene_list)
    kegg = an.run_kegg_enrichment(gene_list)
    if go is not None:
        an.save_results(go, "GO_severe_up.csv")
        an.plot_enrichment_bar(go, title="GO (severe-up)", filename="GO_severe_up.png")
    if kegg is not None:
        an.save_results(kegg, "KEGG_severe_up.csv")
        an.plot_enrichment_bar(kegg, title="KEGG (severe-up)", filename="KEGG_severe_up.png")
    return {"go": go, "kegg": kegg}


if __name__ == "__main__":
    # Self-test of pseudobulk aggregation on a tiny sparse fixture (no scanpy).
    from scipy.sparse import csr_matrix

    # 6 cells x 3 genes; two samples, one cell type
    X = csr_matrix(
        np.array(
            [[5, 0, 1], [3, 0, 0], [4, 1, 0],   # S1 (3 cells)
             [0, 2, 6], [1, 3, 5], [0, 2, 4]],  # S2 (3 cells)
            dtype=float,
        )
    )
    obs = pd.DataFrame(
        {
            "sample_id": ["S1", "S1", "S1", "S2", "S2", "S2"],
            "cell_type": ["Mono"] * 6,
            "condition": ["mild"] * 3 + ["severe"] * 3,
        }
    )
    counts_df, meta_df = aggregate_counts(
        X, obs, ["GENE_A", "GENE_B", "GENE_C"], min_cells=1
    )
    assert counts_df.shape == (3, 2), counts_df.shape
    assert counts_df.loc["GENE_A", "S1|Mono"] == 12  # 5+3+4
    assert counts_df.loc["GENE_C", "S2|Mono"] == 15  # 6+5+4
    assert list(meta_df["condition"]) == ["mild", "severe"]
    assert top_up_genes(
        pd.DataFrame({"log2FoldChange": [2.0, -1.0], "padj": [0.001, 0.2]},
                     index=["UP", "DN"])
    ) == ["UP"]
    print("differential_state.py self-test: pseudobulk aggregation + top_up_genes OK.")
    print(counts_df)
