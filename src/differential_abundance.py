"""
Differential Abundance Module for scRNA-seq Analysis Pipeline
=============================================================

Quantifies how cell-type composition shifts between clinical conditions
(e.g. the monocyte expansion / lymphopenia seen in severe COVID-19).

Approach: compute each sample's cell-type proportions, then compare two
conditions per cell type with a Mann-Whitney U test and BH correction. This is
a transparent, small-cohort-appropriate baseline; for a publication you'd reach
for a compositional model (scCODA) or neighborhood test (Milo), noted in the
docs.

The statistics operate on plain DataFrames (scipy only), so they are fully
unit-testable without scanpy.

Author: scRNA Analysis Pipeline
License: MIT
"""

from __future__ import annotations

import logging
from typing import Optional

import numpy as np
import pandas as pd

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Pure-logic cores (unit-testable)
# ---------------------------------------------------------------------------
def benjamini_hochberg(pvals: pd.Series) -> pd.Series:
    """BH-FDR adjustment that preserves NaNs and the original index."""
    p = pd.Series(pvals, dtype="float64")
    valid = p.dropna()
    if valid.empty:
        return p
    ranked = valid.sort_values()
    n = len(ranked)
    adj = (ranked.values * n) / (np.arange(1, n + 1))
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    adj = np.clip(adj, 0, 1)
    out = p.copy()
    out.loc[ranked.index] = adj
    return out


def proportions_from_obs(
    obs: pd.DataFrame,
    sample_key: str = "sample_id",
    celltype_key: str = "cell_type",
    condition_key: str = "condition",
) -> pd.DataFrame:
    """Long-format per-sample cell-type counts and proportions.

    Columns: [sample_key, celltype_key, n, proportion, condition_key].
    """
    counts = pd.crosstab(obs[sample_key], obs[celltype_key])
    props = counts.div(counts.sum(axis=1), axis=0)

    counts_long = counts.stack().rename("n").reset_index()
    props_long = props.stack().rename("proportion").reset_index()
    merged = counts_long.merge(props_long, on=[sample_key, celltype_key])

    cond = obs.groupby(sample_key, observed=True)[condition_key].first()
    merged[condition_key] = merged[sample_key].map(cond)
    return merged


def test_abundance(
    prop_long: pd.DataFrame,
    group1: str,
    group2: str,
    celltype_key: str = "cell_type",
    condition_key: str = "condition",
    value: str = "proportion",
    min_per_group: int = 2,
) -> pd.DataFrame:
    """Per-cell-type Mann-Whitney U comparison of ``group2`` vs ``group1``.

    ``log2FC`` is log2(median(group2) / median(group1)); a small pseudocount
    avoids division by zero. Returns a table sorted by p-value with BH ``padj``.
    """
    from scipy.stats import mannwhitneyu

    eps = 1e-6
    rows = []
    for ct, sub in prop_long.groupby(celltype_key, observed=True):
        a = sub.loc[sub[condition_key] == group1, value].to_numpy(dtype=float)
        b = sub.loc[sub[condition_key] == group2, value].to_numpy(dtype=float)
        u = p = np.nan
        if len(a) >= min_per_group and len(b) >= min_per_group:
            try:
                u, p = mannwhitneyu(b, a, alternative="two-sided")
            except ValueError:
                u = p = np.nan
        med_a = float(np.median(a)) if len(a) else np.nan
        med_b = float(np.median(b)) if len(b) else np.nan
        rows.append(
            {
                celltype_key: ct,
                f"median_{group1}": med_a,
                f"median_{group2}": med_b,
                "log2FC": float(np.log2((med_b + eps) / (med_a + eps))),
                "U": u,
                "pval": p,
                f"n_{group1}": int(len(a)),
                f"n_{group2}": int(len(b)),
            }
        )
    df = pd.DataFrame(rows)
    if not df.empty:
        df["padj"] = benjamini_hochberg(df["pval"]).values
        df = df.sort_values("pval", na_position="last").reset_index(drop=True)
    return df


# ---------------------------------------------------------------------------
# AnnData-facing wrapper
# ---------------------------------------------------------------------------
def cell_type_proportions(
    adata,
    sample_key: str = "sample_id",
    celltype_key: str = "cell_type",
    condition_key: str = "condition",
) -> pd.DataFrame:
    """Per-sample cell-type proportions from ``adata.obs``."""
    for k in (sample_key, celltype_key, condition_key):
        if k not in adata.obs:
            raise ValueError(f"adata.obs is missing required column '{k}'.")
    return proportions_from_obs(adata.obs, sample_key, celltype_key, condition_key)


if __name__ == "__main__":
    # Self-test of the pure statistics (no scanpy / no AnnData).
    obs = pd.DataFrame(
        {
            "sample_id": [f"S{i}" for i in range(6) for _ in range(100)],
            "condition": (["mild"] * 300) + (["severe"] * 300),
            "cell_type": (
                (["CD14 Mono"] * 20 + ["CD8 T"] * 80) * 3   # mild: 20% mono
                + (["CD14 Mono"] * 60 + ["CD8 T"] * 40) * 3  # severe: 60% mono
            ),
        }
    )
    props = proportions_from_obs(obs)
    res = test_abundance(props, group1="mild", group2="severe")
    mono = res.set_index("cell_type").loc["CD14 Mono"]
    assert mono["log2FC"] > 0, "monocytes should be up in severe in this fixture"
    assert benjamini_hochberg(pd.Series([0.01, np.nan, 0.5])).isna().sum() == 1
    print("differential_abundance.py self-test: proportions + test + BH OK.")
    print(res.to_string(index=False))
