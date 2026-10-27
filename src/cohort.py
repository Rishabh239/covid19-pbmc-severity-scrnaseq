"""
Cohort Module for scRNA-seq Analysis Pipeline
=============================================

Turns the single-sample loader into a case/control cohort workflow.

The existing ``data_loader.load_10x_data`` reads ONE 10x matrix and knows
nothing about experimental groups. For a disease study we need every cell
tagged with its sample of origin and its clinical condition, so we can later
integrate across samples and run condition-vs-condition comparisons.

Written for GSE149689 (Lee et al., Sci. Immunol. 2020 - COVID-19 & influenza
PBMCs), but the mechanics are generic to any combined 10x matrix whose cell
barcodes carry a trailing sample index (``...-1`` / ``...-20``).

Author: scRNA Analysis Pipeline
License: MIT
"""

from __future__ import annotations

import re
import logging
from collections import Counter
from pathlib import Path
from typing import Dict, Optional, Sequence, Union

import numpy as np
import pandas as pd

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# GSE149689 reference metadata
# ---------------------------------------------------------------------------
# Known cohort composition (Lee et al. 2020; ~59,572 PBMCs across 20 samples).
# Used only to *validate* that your suffix->condition map is self-consistent.
EXPECTED_GROUP_COUNTS: Dict[str, int] = {
    "healthy": 4,
    "mild": 4,
    "severe": 6,
    "asymptomatic": 1,
    "flu": 5,
}

# The GEO supplementary matrix concatenates all 20 samples. Each cell barcode
# ends in "-N" (or "_N"), where N is the sample's *position* in the GEO
# "Samples" list. GEO encodes only the sample index in the barcode, NOT the
# condition, so the index->condition map must be read off the sample titles.
#
# HOW TO FILL THIS (do it once, ~5 min):
#   1. Open the GSE149689 GEO page -> "Samples (20)" section.
#   2. GSMs are listed in order: the 1st row is suffix -1, the 2nd is -2, ...
#   3. Read each sample's condition and enter it below, using EXACTLY these
#      labels: healthy / mild / severe / asymptomatic / flu.
#
# FILLED for GSE149689: the GEO titles only give nCoV / Flu / Normal, so the
# COVID severities below come from the paper's Supplementary Table S1 (Lee et
# al. 2020), which lists each patient's clinical "Disease group". Mapping used:
#   severe COVID-19 -> severe   |   mild COVID-19 -> mild
#   mild COVID-19 (asymptomatic) -> asymptomatic
#   severe influenza -> flu     |   Healthy donor -> healthy
# validate_sample_map() confirms these sum to 6 severe / 4 mild / 1 asymptomatic
# / 5 flu / 4 healthy. If you swap in a different matrix, re-verify the order.
SAMPLE_CONDITIONS: Dict[int, Optional[str]] = {
    1:  "severe",        # nCoV 1   (C1)  severe COVID-19
    2:  "mild",          # nCoV 2   (C2)  mild COVID-19
    3:  "flu",           # Flu 1    (F1)  severe influenza
    4:  "flu",           # Flu 2    (F2)  severe influenza
    5:  "healthy",       # Normal 1 (N1)  healthy donor
    6:  "flu",           # Flu 3    (F3)  severe influenza
    7:  "flu",           # Flu 4    (F4)  severe influenza
    8:  "flu",           # Flu 5    (F5)  severe influenza
    9:  "severe",        # nCoV 3   (C3)  severe COVID-19
    10: "severe",        # nCoV 4         severe COVID-19
    11: "mild",          # nCoV 5   (C4)  mild COVID-19
    12: "mild",          # nCoV 6   (C5)  mild COVID-19
    13: "healthy",       # Normal 2 (N2)  healthy donor
    14: "healthy",       # Normal 3 (N3)  healthy donor
    15: "severe",        # nCoV 7   (C6)  severe COVID-19
    16: "severe",        # nCoV 8         severe COVID-19
    17: "severe",        # nCoV 9   (C7)  severe COVID-19
    18: "mild",          # nCoV 10        mild COVID-19
    19: "healthy",       # Normal 4 (N4)  healthy donor
    20: "asymptomatic",  # nCoV 11  (C8)  mild COVID-19 (asymptomatic)
}


_SUFFIX_RE = re.compile(r"[-_](\d+)$")


# ---------------------------------------------------------------------------
# Pure-logic helpers (no scanpy needed - unit-testable on their own)
# ---------------------------------------------------------------------------
def parse_suffix(barcode: str) -> int:
    """Extract the trailing sample index from a barcode.

    'AAACCTGAGAAACGAG-1' -> 1,  'AAAC...._20' -> 20.
    If a barcode carries several trailing indices (e.g. '-1-3'), the LAST one
    is taken, consistent with "the sample index is the final suffix".
    """
    m = _SUFFIX_RE.search(barcode)
    if not m:
        raise ValueError(
            f"Could not parse a sample suffix from barcode {barcode!r}; "
            "expected a trailing '-N' or '_N'."
        )
    return int(m.group(1))


def validate_sample_map(sample_map: Dict[int, Optional[str]]) -> None:
    """Fail loudly unless the suffix->condition map is complete and its group
    sizes match the known GSE149689 composition."""
    unfilled = [k for k, v in sample_map.items() if v in (None, "", "FILL_ME")]
    if unfilled:
        raise ValueError(
            f"SAMPLE_CONDITIONS is still unfilled for suffixes {unfilled}. "
            "Fill them from the GSE149689 GEO 'Samples' list before running."
        )

    counts = Counter(sample_map.values())
    unknown = set(counts) - set(EXPECTED_GROUP_COUNTS)
    if unknown:
        raise ValueError(
            f"Unexpected condition label(s) {sorted(unknown)}. "
            f"Use only {sorted(EXPECTED_GROUP_COUNTS)}."
        )

    mismatches = {
        g: (counts.get(g, 0), n)
        for g, n in EXPECTED_GROUP_COUNTS.items()
        if counts.get(g, 0) != n
    }
    if mismatches:
        raise ValueError(
            "Group counts in your map do not match the known GSE149689 "
            f"composition (got, expected): {mismatches}. "
            "Re-check the sample order on GEO."
        )
    logger.info("Sample map validated: %s", dict(counts))


# ---------------------------------------------------------------------------
# AnnData-facing functions
# ---------------------------------------------------------------------------
def summarize_suffixes(adata) -> pd.DataFrame:
    """Cells per barcode suffix. Cross-check the row count == 20 and the
    per-suffix sizes against the GEO sample list before trusting the map."""
    suffixes = [parse_suffix(bc) for bc in adata.obs_names]
    s = pd.Series(suffixes).value_counts().sort_index()
    return s.rename_axis("suffix").reset_index(name="n_cells")


def label_cells(
    adata,
    sample_map: Dict[int, Optional[str]] = SAMPLE_CONDITIONS,
    validate: bool = True,
):
    """Add 'sample_id' and 'condition' columns to ``adata.obs`` from the suffix map."""
    if validate:
        validate_sample_map(sample_map)

    suffixes = np.array([parse_suffix(bc) for bc in adata.obs_names])
    found = sorted(set(int(x) for x in suffixes))
    expected = sorted(sample_map)
    if found != expected:
        raise ValueError(
            f"Barcode suffixes present {found} do not match the mapped "
            f"suffixes {expected}. Either this isn't the 20-sample combined "
            "matrix, or the barcode format differs from '<bc>-N'."
        )

    adata.obs["sample_id"] = pd.Categorical([f"S{int(n)}" for n in suffixes])
    adata.obs["condition"] = pd.Categorical([sample_map[int(n)] for n in suffixes])
    logger.info("Labeled %d cells across %d samples.", adata.n_obs, len(found))
    logger.info(
        "Cells per condition:\n%s", adata.obs["condition"].value_counts().to_string()
    )
    return adata


def subset_conditions(adata, keep: Sequence[str] = ("healthy", "mild", "severe")):
    """Keep only the requested conditions (drops flu + asymptomatic by default,
    leaving a clean healthy/mild/severe contrast)."""
    if "condition" not in adata.obs:
        raise ValueError("Run label_cells() before subset_conditions().")
    mask = adata.obs["condition"].isin(list(keep)).values
    out = adata[mask].copy()
    out.obs["condition"] = out.obs["condition"].cat.remove_unused_categories()
    logger.info("Subset to %s: %d cells.", list(keep), out.n_obs)
    return out


def per_sample_qc(adata) -> pd.DataFrame:
    """Tidy per-sample table (condition, n_cells, and medians of any standard
    QC columns already present) for a quick cohort QC report."""
    obs = adata.obs
    rows = []
    qc_cols = [
        c
        for c in ("n_genes_by_counts", "n_genes", "total_counts", "pct_counts_mt")
        if c in obs.columns
    ]
    for sid, grp in obs.groupby("sample_id", observed=True):
        row = {
            "sample_id": sid,
            "condition": grp["condition"].iloc[0] if "condition" in grp else "NA",
            "n_cells": len(grp),
        }
        for c in qc_cols:
            row[f"median_{c}"] = float(np.median(grp[c].values))
        rows.append(row)
    return pd.DataFrame(rows).sort_values("sample_id").reset_index(drop=True)


def load_cohort(
    data_dir: Union[str, Path],
    sample_map: Dict[int, Optional[str]] = SAMPLE_CONDITIONS,
    keep: Sequence[str] = ("healthy", "mild", "severe"),
    min_genes: int = 200,
    subset: bool = True,
):
    """High-level entry point: load the combined 10x matrix, label conditions,
    and (optionally) subset to the conditions of interest.

    The pipeline's own loader is imported lazily so this module stays importable
    (and testable) without scanpy installed.
    """
    try:
        from src.data_loader import load_10x_data  # when run as part of the package
    except Exception:  # pragma: no cover - fallback for standalone use
        from data_loader import load_10x_data

    adata = load_10x_data(data_dir, min_genes=min_genes)
    logger.info(
        "Suffix summary before labeling:\n%s",
        summarize_suffixes(adata).to_string(index=False),
    )
    adata = label_cells(adata, sample_map=sample_map)
    if subset:
        adata = subset_conditions(adata, keep=keep)
    return adata


# ---------------------------------------------------------------------------
# Built-in self-test of the pure logic (runs without scanpy or data)
# ---------------------------------------------------------------------------
def _selftest() -> None:
    # parse_suffix
    assert parse_suffix("AAACCTGAGAAACGAG-1") == 1
    assert parse_suffix("TTTTGATCACGTAGGG_20") == 20
    assert parse_suffix("GGG-1-3") == 3
    try:
        parse_suffix("NOSUFFIXHERE")
    except ValueError:
        pass
    else:
        raise AssertionError("parse_suffix should reject a barcode with no suffix")

    # the shipped map is now filled (from Table S1) and must validate
    validate_sample_map(SAMPLE_CONDITIONS)

    # an explicitly empty map must still be rejected
    try:
        validate_sample_map({i: None for i in range(1, 21)})
    except ValueError:
        pass
    else:
        raise AssertionError("empty map should not validate")

    # a *fixture* map with the correct group sizes must pass
    good = {}
    order = (["healthy"] * 4 + ["mild"] * 4 + ["severe"] * 6
             + ["asymptomatic"] * 1 + ["flu"] * 5)
    for i, cond in enumerate(order, start=1):
        good[i] = cond
    validate_sample_map(good)

    # a map with a wrong group size must fail
    bad = dict(good)
    bad[1] = "mild"  # now 5 mild / 3 healthy
    try:
        validate_sample_map(bad)
    except ValueError:
        pass
    else:
        raise AssertionError("wrong group counts should not validate")

    print("cohort.py self-test: all pure-logic checks passed.")


if __name__ == "__main__":
    _selftest()
