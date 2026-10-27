#!/usr/bin/env python
"""
make_figures.py - per-cell-type figures for the cohort DE results
================================================================

Produces up to three figures for one cell type:

  1. Volcano   - from a DESeq2 results CSV (always).
  2. Heatmap   - z-scored pseudobulk expression of the top DE genes, with a
                 condition color bar (needs pseudobulk_counts.csv + metadata).
  3. Pathway   - GO/KEGG enrichment of the up-in-<test> genes via gseapy Enrichr
                 (needs the `gseapy` package and an internet connection).

Each figure is skipped gracefully if its inputs/packages aren't available.
Outputs are named for --celltype, e.g. volcano_CD14_Mono.png.

Note: the pathway step submits only the top `--top-genes` genes by fold change,
not every gene passing the thresholds.

This is the exploratory script. The figures in README.md are built by
make_report.py; the one thing taken from here is the monocyte Enrichr table,
committed at results/enrichment/CD14_Mono/enrichr_severe_up_top15.csv.

Example (all three):
  python make_figures.py \
      --de results/pseudobulk/deseq2_CD14_Mono_severe_vs_mild.csv \
      --pseudobulk results/pseudobulk/pseudobulk_counts.csv \
      --meta results/pseudobulk/pseudobulk_metadata.csv \
      --celltype "CD14 Mono" --test severe --ref mild --outdir figures
"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger("make_figures")

# Genes worth calling out on the monocyte volcano (only labeled if present).
HIGHLIGHT_UP = [
    "S100A8", "S100A9", "S100A12", "SERPINB2", "IL1R2", "IL18R1", "CXCR1",
    "MS4A3", "PRTN3", "MKI67", "HP", "THBS1", "CLU", "IFITM1", "IFITM2",
]
HIGHLIGHT_DOWN = [
    "HLA-DQA1", "HLA-DOA", "GNLY", "GZMB", "IFNG", "KLRB1", "FGFBP2",
]


def _slug(celltype: str) -> str:
    """Filename-safe form of a cell-type label: 'CD14 Mono' -> 'CD14_Mono'."""
    return str(celltype).replace(" ", "_").replace("/", "-")


def _load_de(path: str) -> pd.DataFrame:
    de = pd.read_csv(path, index_col=0)
    lfc = "log2FoldChange" if "log2FoldChange" in de.columns else \
        next(c for c in de.columns if "log2" in c.lower())
    padj = "padj" if "padj" in de.columns else \
        next(c for c in de.columns if "adj" in c.lower() or "fdr" in c.lower())
    de = de.rename(columns={lfc: "log2FoldChange", padj: "padj"})
    return de.dropna(subset=["padj", "log2FoldChange"])


def volcano(de: pd.DataFrame, out: Path, celltype: str, test: str, ref: str,
            lfc_thresh: float = 1.0, padj_thresh: float = 0.05) -> None:
    d = de.copy()
    # floor padj==0 so -log10 is finite; plot those at the top
    tiny = d.loc[d["padj"] > 0, "padj"].min()
    d["padj"] = d["padj"].clip(lower=tiny)
    d["nlp"] = -np.log10(d["padj"])

    up = (d["padj"] < padj_thresh) & (d["log2FoldChange"] >= lfc_thresh)
    dn = (d["padj"] < padj_thresh) & (d["log2FoldChange"] <= -lfc_thresh)
    ns = ~(up | dn)

    fig, ax = plt.subplots(figsize=(7.5, 7))
    ax.scatter(d.loc[ns, "log2FoldChange"], d.loc[ns, "nlp"], s=6,
               c="#c9ccd1", alpha=0.5, linewidths=0, rasterized=True, label="n.s.")
    ax.scatter(d.loc[up, "log2FoldChange"], d.loc[up, "nlp"], s=12,
               c="#c0392b", alpha=0.8, linewidths=0, label=f"up in {test} ({int(up.sum())})")
    ax.scatter(d.loc[dn, "log2FoldChange"], d.loc[dn, "nlp"], s=12,
               c="#2471a3", alpha=0.8, linewidths=0, label=f"up in {ref} ({int(dn.sum())})")

    ax.axhline(-np.log10(padj_thresh), ls="--", lw=0.8, c="0.5")
    ax.axvline(lfc_thresh, ls="--", lw=0.8, c="0.5")
    ax.axvline(-lfc_thresh, ls="--", lw=0.8, c="0.5")

    labels = [g for g in HIGHLIGHT_UP + HIGHLIGHT_DOWN if g in d.index]
    for g in labels:
        x, y = d.loc[g, "log2FoldChange"], d.loc[g, "nlp"]
        ax.annotate(g, (x, y), fontsize=7.5, fontstyle="italic",
                    xytext=(4, 3), textcoords="offset points",
                    color="#7b241c" if x > 0 else "#1a5276")

    ax.set_xlabel(f"log$_2$ fold change  ({test} / {ref})", fontsize=11)
    ax.set_ylabel(r"$-\log_{10}$ adjusted $p$", fontsize=11)
    ax.set_title(f"{celltype}: {test} vs {ref} COVID-19", fontsize=12, fontweight="bold")
    ax.legend(frameon=False, fontsize=9, loc="upper left")
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    fig.tight_layout()
    fig.savefig(out / f"volcano_{_slug(celltype)}.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    logger.info("Wrote %s", out / f"volcano_{_slug(celltype)}.png")


def heatmap(de: pd.DataFrame, pseudobulk_path: str, meta_path: str, out: Path,
            celltype: str, test: str, ref: str, top_n: int = 25) -> None:
    try:
        import seaborn as sns
    except ImportError:
        logger.warning("seaborn not installed (pip install seaborn); skipping heatmap.")
        return

    counts = pd.read_csv(pseudobulk_path, index_col=0)      # genes x pseudobulk_id
    meta = pd.read_csv(meta_path, index_col=0)
    ct_col = "cell_type" if "cell_type" in meta.columns else \
        next(c for c in meta.columns if "cell" in c.lower())
    cond_col = "condition" if "condition" in meta.columns else \
        next(c for c in meta.columns if "cond" in c.lower())

    keep = meta[(meta[ct_col] == celltype) & (meta[cond_col].isin([test, ref]))]
    if keep.empty:
        logger.warning("No pseudobulk samples for %s in {%s,%s}; skipping heatmap.",
                       celltype, test, ref)
        return
    cols = [c for c in keep.index if c in counts.columns]
    mat = counts[cols]

    # top up + top down DE genes present in the matrix
    sig = de[de["padj"] < 0.05]
    up = sig[sig["log2FoldChange"] > 0].sort_values("log2FoldChange", ascending=False)
    dn = sig[sig["log2FoldChange"] < 0].sort_values("log2FoldChange")
    genes = [g for g in list(up.index[:top_n]) + list(dn.index[:top_n]) if g in mat.index]

    # CPM-log then z-score per gene across samples
    cpm = mat.div(mat.sum(axis=0), axis=1) * 1e6
    logmat = np.log1p(cpm.loc[genes])
    z = logmat.sub(logmat.mean(axis=1), axis=0).div(logmat.std(axis=1).replace(0, 1), axis=0)

    cond = keep.loc[cols, cond_col]
    palette = {test: "#c0392b", ref: "#2471a3"}
    col_colors = cond.map(palette)

    g = sns.clustermap(
        z, cmap="RdBu_r", center=0, vmin=-2, vmax=2,
        col_colors=col_colors.values, xticklabels=False,
        yticklabels=True, figsize=(max(6, 0.4 * len(cols)), 0.22 * len(genes) + 2),
        dendrogram_ratio=(0.12, 0.06), cbar_pos=(0.02, 0.83, 0.03, 0.12),
    )
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), fontsize=7, fontstyle="italic")
    handles = [plt.Rectangle((0, 0), 1, 1, color=palette[k]) for k in (test, ref)]
    g.ax_col_dendrogram.legend(handles, [test, ref], frameon=False,
                               loc="center", ncol=2, fontsize=9)
    g.fig.suptitle(f"{celltype}: top DE genes ({test} vs {ref})", y=1.02, fontsize=12)
    g.savefig(out / f"heatmap_{_slug(celltype)}.png", dpi=300, bbox_inches="tight")
    plt.close(g.fig)
    logger.info("Wrote %s", out / f"heatmap_{_slug(celltype)}.png")


def pathway(de: pd.DataFrame, out: Path, celltype: str, test: str,
            top_genes: int = 200,
            padj_thresh: float = 0.05, lfc_thresh: float = 1.0) -> None:
    try:
        import gseapy as gp
    except ImportError:
        logger.warning("gseapy not installed (pip install gseapy); skipping pathway.")
        return

    sig_up = de[(de["padj"] < padj_thresh) & (de["log2FoldChange"] >= lfc_thresh)]
    genes = sig_up.sort_values("log2FoldChange", ascending=False).index.astype(str).tolist()[:top_genes]
    if len(genes) < 5:
        logger.warning("Too few up genes (%d) for enrichment; skipping.", len(genes))
        return

    gene_sets = ["GO_Biological_Process_2021", "MSigDB_Hallmark_2020", "KEGG_2021_Human"]
    try:
        enr = gp.enrichr(gene_list=genes, gene_sets=gene_sets, organism="human",
                         outdir=str(out / "enrichr"))
    except Exception as e:
        logger.warning("Enrichr failed (needs internet): %s; skipping pathway.", e)
        return

    res = enr.results.copy()
    res = res[res["Adjusted P-value"] < 0.1].sort_values("Adjusted P-value").head(15)
    if res.empty:
        logger.warning("No enriched terms at FDR<0.1; skipping plot.")
        return
    res["nlp"] = -np.log10(res["Adjusted P-value"])
    res = res.iloc[::-1]

    fig, ax = plt.subplots(figsize=(8, 0.45 * len(res) + 1.5))
    ax.barh(res["Term"].str.slice(0, 55), res["nlp"], color="#7d3c98")
    ax.set_xlabel(r"$-\log_{10}$ adjusted $p$", fontsize=11)
    ax.set_title(f"Enriched pathways in up-in-{test} monocyte genes", fontsize=12, fontweight="bold")
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    fig.tight_layout()
    fig.savefig(out / f"pathways_{_slug(celltype)}.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    res.to_csv(out / f"pathways_{_slug(celltype)}.csv", index=False)
    logger.info("Wrote %s", out / f"pathways_{_slug(celltype)}.png")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--de", required=True, help="DESeq2 results CSV (gene index).")
    ap.add_argument("--pseudobulk", help="pseudobulk_counts.csv (genes x pseudobulk_id).")
    ap.add_argument("--meta", help="pseudobulk_metadata.csv.")
    ap.add_argument("--celltype", default="CD14 Mono")
    ap.add_argument("--test", default="severe")
    ap.add_argument("--ref", default="mild")
    ap.add_argument("--outdir", default="figures")
    ap.add_argument("--top-genes", type=int, default=200,
                    help="max genes (by fold change) submitted to Enrichr (default: 200)")
    args = ap.parse_args(argv)

    out = Path(args.outdir); out.mkdir(parents=True, exist_ok=True)
    de = _load_de(args.de)

    volcano(de, out, args.celltype, args.test, args.ref)
    if args.pseudobulk and args.meta:
        heatmap(de, args.pseudobulk, args.meta, out, args.celltype, args.test, args.ref)
    else:
        logger.info("No --pseudobulk/--meta given; skipping heatmap.")
    pathway(de, out, args.celltype, args.test, top_genes=args.top_genes)
    logger.info("Done. Figures in %s", out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
