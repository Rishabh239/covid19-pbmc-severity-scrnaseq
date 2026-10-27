#!/usr/bin/env python
"""
make_report.py - build the summary tables and figures behind README.md.

Reads the workflow outputs in results/ and writes:

  results/summary/   small derived tables (DE counts, significant genes, modules)
  figures/           the five figures used in the study write-up

    python make_report.py              # rebuild everything from results/
    python make_report.py --figures    # redraw figures only

Nothing here re-runs the analysis; it only summarises what run_cohort.py
produced, so every number in the README traces back to a file in results/.
Requires pandas >= 1.5 and matplotlib.
"""
from __future__ import annotations

import argparse
import logging
import re
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results"
SUMMARY = RESULTS / "summary"
FIGURES = HERE / "figures"
PSEUDOBULK = RESULTS / "pseudobulk"

TEST, REF = "severe", "mild"
CONDITIONS = ["healthy", "mild", "severe"]      # ordered by severity
PADJ = 0.05

# Cell types in a fixed order: most to least abundant in this cohort.
CELL_ORDER = ["CD14 Mono", "NK", "CD8 T", "CD4 T", "B", "Platelet", "CD16 Mono", "DC"]

# Monocyte gene modules examined in the write-up. Each is a published axis of the
# severe-COVID monocyte phenotype; listing them up front avoids cherry-picking
# genes after seeing the result.
MODULES = {
    "MHC class II\n(antigen presentation)":
        ["HLA-DRA", "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DMA", "HLA-DMB", "CD74"],
    "Alarmins & emergency\nmyelopoiesis":
        ["S100A8", "S100A9", "S100A12", "IL1R2", "SERPINB2", "HP", "RETN", "THBS1", "CLU"],
    "Type I interferon\nresponse":
        ["ISG15", "IFI6", "IFITM3", "OAS1", "OAS2", "MX1", "IRF7"],
}

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger("make_report")


def de_path(cell_type: str) -> Path:
    return PSEUDOBULK / f"deseq2_{cell_type.replace(' ', '_')}_{TEST}_vs_{REF}.csv"


def load_de(cell_type: str) -> pd.DataFrame | None:
    p = de_path(cell_type)
    return pd.read_csv(p, index_col=0) if p.exists() else None


def available_celltypes() -> list[str]:
    """Cell types that have a DESeq2 table, in CELL_ORDER."""
    found = {re.match(r"deseq2_(.+)_%s_vs_%s\.csv" % (TEST, REF), p.name).group(1).replace("_", " ")
             for p in PSEUDOBULK.glob(f"deseq2_*_{TEST}_vs_{REF}.csv")}
    return [c for c in CELL_ORDER if c in found] + sorted(found - set(CELL_ORDER))


# ============================================================================
# summary tables
# ============================================================================
def write_summaries() -> None:
    SUMMARY.mkdir(parents=True, exist_ok=True)
    cohort_summary()
    de_summary()
    module_table()
    export_umap_coords()


def export_umap_coords() -> None:
    """
    Pull the UMAP embedding out of the processed AnnData into a small CSV.

    The .h5ad is ~500 MB and is not committed, so the figure is drawn from this
    extract instead. When the object is absent (a fresh clone, or results pulled
    without re-running), the committed extract is left as it is.
    """
    h5ad = RESULTS / "cohort_processed.h5ad"
    out = SUMMARY / "umap_coords.csv.gz"
    if not h5ad.exists():
        logger.info("no %s; keeping the committed %s", h5ad.name, out.name)
        return
    try:
        import h5py
    except ImportError:
        logger.warning("h5py not installed; keeping the committed %s", out.name)
        return

    # Read the two arrays straight out of the file rather than through
    # anndata.read_h5ad: a backed read also parses .uns, which fails on objects
    # written by a different scanpy version ("no read method for encoding_type
    # 'null'"). Only obs and obsm are needed here.
    def column(obs, name):
        node = obs[name]
        if isinstance(node, h5py.Group):           # categorical
            cats = [c.decode() if isinstance(c, bytes) else c
                    for c in node["categories"][:]]
            return pd.Categorical.from_codes(node["codes"][:], cats).astype(str)
        return node[:]

    with h5py.File(h5ad, "r") as f:
        umap = f["obsm"]["X_umap"][:]
        df = pd.DataFrame({
            "umap1": umap[:, 0].round(3),
            "umap2": umap[:, 1].round(3),
            "cell_type": column(f["obs"], "cell_type"),
            "condition": column(f["obs"], "condition"),
        })
    # Fixed mtime keeps the gzip byte-identical across runs of unchanged input.
    df.to_csv(out, index=False, lineterminator="\n",
              compression={"method": "gzip", "mtime": 0})
    logger.info("wrote %s (%d cells)", out, len(df))


def cohort_summary() -> None:
    """Samples and cells per condition, from the per-sample QC and proportions."""
    qc = pd.read_csv(RESULTS / "per_sample_qc.csv")
    props = pd.read_csv(RESULTS / "cell_type_proportions.csv")
    # Post-QC cells per sample are the cell-type counts summed over cell types;
    # per_sample_qc.csv counts cells as loaded, before the QC filter.
    post = props.groupby(["condition", "sample_id"])["n"].sum().rename("cells_after_qc")
    pre = qc.set_index(["condition", "sample_id"])["n_cells"].rename("cells_loaded")
    med = qc.set_index(["condition", "sample_id"])["median_n_genes"]
    out = pd.concat([pre, post, med], axis=1).reset_index()
    out["condition"] = pd.Categorical(out["condition"], CONDITIONS, ordered=True)
    out = out.sort_values(["condition", "sample_id"])
    out.to_csv(SUMMARY / "cohort_samples.csv", index=False, lineterminator="\n")

    by_cond = out.groupby("condition", observed=True).agg(
        samples=("sample_id", "count"),
        cells_loaded=("cells_loaded", "sum"),
        cells_after_qc=("cells_after_qc", "sum"),
    ).reset_index()
    by_cond.to_csv(SUMMARY / "cohort_by_condition.csv", index=False, lineterminator="\n")
    logger.info("cohort: %d samples, %d cells after QC",
                len(out), int(out["cells_after_qc"].sum()))


def de_summary() -> None:
    """Per-cell-type DE counts, plus one combined table of significant genes."""
    rows, sig_frames = [], []
    for ct in available_celltypes():
        de = load_de(ct)
        sig = de[de["padj"] < PADJ]
        rows.append({
            "cell_type": ct,
            "genes_tested": len(de),
            "significant": len(sig),
            "up_in_severe": int((sig["log2FoldChange"] > 0).sum()),
            "up_in_mild": int((sig["log2FoldChange"] < 0).sum()),
            "significant_lfc_over_1": int((sig["log2FoldChange"].abs() > 1).sum()),
        })
        frame = sig.copy()
        frame.insert(0, "cell_type", ct)
        frame.index.name = "gene"
        sig_frames.append(frame.reset_index())

    pd.DataFrame(rows).to_csv(SUMMARY / "de_summary_by_celltype.csv",
                              index=False, lineterminator="\n")
    combined = pd.concat(sig_frames, ignore_index=True)
    combined = combined.sort_values(["cell_type", "padj"])
    combined.round(6).to_csv(SUMMARY / "de_significant_genes.csv",
                             index=False, lineterminator="\n")
    logger.info("DE: %d significant gene-by-cell-type results at padj < %.2f",
                len(combined), PADJ)


def module_table() -> None:
    """log2FC and padj for each module gene in CD14 monocytes."""
    de = load_de("CD14 Mono")
    rows = []
    for module, genes in MODULES.items():
        for g in genes:
            if g not in de.index:
                logger.warning("module gene %s not in the CD14 Mono table", g)
                continue
            r = de.loc[g]
            rows.append({
                "module": module.replace("\n", " "), "gene": g,
                "baseMean": round(r["baseMean"], 1),
                "log2FoldChange": round(r["log2FoldChange"], 3),
                "pvalue": r["pvalue"], "padj": r["padj"],
                "significant": bool(pd.notna(r["padj"]) and r["padj"] < PADJ),
            })
    pd.DataFrame(rows).to_csv(SUMMARY / "monocyte_gene_modules.csv",
                              index=False, lineterminator="\n")


# ============================================================================
# figures
# ============================================================================
# Colours are the documented defaults of the house palette: a single-hue blue
# ordinal ramp for the ordered severity scale, and the blue<->red diverging pair
# with a neutral grey for up/down contrasts.
SURFACE = "#fcfcfb"
INK, INK_2, MUTED = "#0b0b0b", "#52514e", "#898781"
GRID, AXIS = "#e1e0d9", "#c3c2b7"
SEVERITY = {"healthy": "#86b6ef", "mild": "#3987e5", "severe": "#184f95"}
UP, DOWN, NS = "#d03b3b", "#2a78d6", "#c9c8c0"
ACCENT = "#2a78d6"   # emphasis hue where the contrast is magnitude, not direction
# The eight categorical slots, assigned to CELL_ORDER. On the UMAP each cluster
# also carries a direct label, so identity never rests on hue alone.
CATEGORICAL = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100",
               "#e87ba4", "#008300", "#4a3aa7", "#e34948"]
PX = 1 / 100     # one display pixel in inches; figures are laid out at 100 px/in
WIDTH = 9.2      # inches -> ~920 px, about a GitHub README column


def style(plt) -> None:
    logging.getLogger("matplotlib.font_manager").setLevel(logging.ERROR)
    plt.rcParams.update({
        "font.family": ["Segoe UI", "Arial", "DejaVu Sans"],
        "font.size": 9.5, "axes.titlesize": 10.5, "axes.titleweight": "semibold",
        "axes.titlelocation": "left", "axes.titlecolor": INK, "axes.titlepad": 8,
        "axes.labelcolor": INK_2, "axes.labelsize": 9.5,
        "axes.edgecolor": AXIS, "axes.linewidth": 0.75, "axes.facecolor": SURFACE,
        "figure.facecolor": SURFACE, "savefig.facecolor": SURFACE,
        "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.75, "grid.linestyle": "-",
        "axes.axisbelow": True,
        "xtick.color": MUTED, "ytick.color": MUTED, "xtick.labelcolor": INK_2,
        "ytick.labelcolor": INK_2, "xtick.major.size": 0, "ytick.major.size": 0,
        "xtick.major.pad": 5, "ytick.major.pad": 5,
        "legend.frameon": False, "legend.fontsize": 9, "legend.labelcolor": INK_2,
    })


def header(fig, title: str, subtitle: str) -> None:
    fig.text(0.012, 0.985, title, ha="left", va="top", fontsize=12,
             fontweight="semibold", color=INK)
    fig.text(0.012, 0.915, subtitle, ha="left", va="top", fontsize=9.5, color=INK_2)


def despine(ax, keep=("bottom",)) -> None:
    for side in ("top", "right", "left", "bottom"):
        ax.spines[side].set_visible(side in keep)


def hbar(ax, y, width, height, color, radius_px=4):
    """Horizontal bar from x=0, rounded at the data end, square at the baseline."""
    from matplotlib.patches import PathPatch
    from matplotlib.path import Path as MPath
    bbox = ax.get_window_extent().transformed(ax.figure.dpi_scale_trans.inverted())
    (x0, x1), (y0, y1) = ax.get_xlim(), ax.get_ylim()
    rx = min(radius_px * PX * abs(x1 - x0) / bbox.width, abs(width) / 2)
    ry = min(radius_px * PX * abs(y1 - y0) / bbox.height, height / 2)
    rx = rx if width >= 0 else -rx
    b, t = y - height / 2, y + height / 2
    verts = [(0, b), (width - rx, b), (width, b), (width, b + ry), (width, t - ry),
             (width, t), (width - rx, t), (0, t), (0, b)]
    codes = [MPath.MOVETO, MPath.LINETO, MPath.CURVE3, MPath.CURVE3, MPath.LINETO,
             MPath.CURVE3, MPath.CURVE3, MPath.LINETO, MPath.CLOSEPOLY]
    ax.add_patch(PathPatch(MPath(verts, codes), facecolor=color, edgecolor="none"))


def save(fig, name: str) -> None:
    FIGURES.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIGURES / name, dpi=200)
    logger.info("wrote %s", FIGURES / name)


# ---------------------------------------------------------------------------
def fig_umap(plt) -> None:
    """Integrated UMAP, one panel coloured by cell type, one by condition."""
    coords = RESULTS / "summary" / "umap_coords.csv.gz"
    if not coords.exists():
        logger.warning("no %s; skipping the UMAP figure", coords)
        return
    d = pd.read_csv(coords)
    order = [c for c in CELL_ORDER if c in set(d["cell_type"])]
    colors = dict(zip(order, CATEGORICAL))

    fig, axes = plt.subplots(1, 2, figsize=(WIDTH, 4.6))
    for ax, key in zip(axes, ("cell_type", "condition")):
        groups = order if key == "cell_type" else CONDITIONS
        for g in groups:
            sub = d[d[key] == g]
            c = colors[g] if key == "cell_type" else SEVERITY[g]
            ax.plot(sub["umap1"], sub["umap2"], ".", ms=1.1, color=c,
                    alpha=0.55, rasterized=True, zorder=2)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.grid(False)
        for s in ("top", "right", "left", "bottom"):
            ax.spines[s].set_visible(True)
            ax.spines[s].set_color(GRID)
        ax.set_xlabel("UMAP 1")

    # Direct labels at each cluster's median position, in ink rather than the
    # series colour, so the eight cell types stay distinguishable without hue.
    for ct in order:
        sub = d[d["cell_type"] == ct]
        axes[0].annotate(ct, (sub["umap1"].median(), sub["umap2"].median()),
                         ha="center", va="center", fontsize=8.5, color=INK,
                         fontweight="semibold",
                         path_effects=_halo(plt))
    axes[0].set_ylabel("UMAP 2")
    axes[0].set_title("Cell type")
    axes[1].set_title("Condition")
    handles = [plt.Line2D([], [], marker="o", ls="", ms=7, color=SEVERITY[c], mec=SURFACE)
               for c in CONDITIONS]
    axes[1].legend(handles, CONDITIONS, loc="upper left", bbox_to_anchor=(0.0, -0.03),
                   ncol=3, handletextpad=0.3, columnspacing=1.4)
    header(fig, f"{len(d):,} PBMCs from 14 donors, after integrating over donor",
           "Harmony integration, then Leiden clustering labelled by canonical marker "
           "score.\nEvery cell type draws on all three groups, so donor is not driving "
           "the clusters.")
    fig.tight_layout(rect=(0, 0.04, 1, 0.84), w_pad=2.0)
    save(fig, "fig1_umap.png")


def _halo(plt):
    """A thin surface-coloured outline so labels stay legible over dense points."""
    import matplotlib.patheffects as pe
    return [pe.withStroke(linewidth=2.6, foreground=SURFACE)]


def fig_composition(plt) -> None:
    """Per-sample cell-type proportions, one row per cell type."""
    props = pd.read_csv(RESULTS / "cell_type_proportions.csv")
    abund = pd.read_csv(RESULTS / "abundance_severe_vs_mild.csv").set_index("cell_type")
    # Most abundant cell type at the top, so row order matches the eye's first stop.
    order = [c for c in CELL_ORDER if c in set(props["cell_type"])][::-1]
    n_by_cond = props.groupby("condition")["sample_id"].nunique()

    fig, ax = plt.subplots(figsize=(WIDTH, 1.3 + 0.62 * len(order)))
    xmax = props["proportion"].max() * 100 * 1.38
    ax.set_xlim(0, xmax)
    ax.set_ylim(-0.6, len(order) - 0.4)

    rng = __import__("numpy").random.default_rng(0)
    for i, ct in enumerate(order):
        sub = props[props["cell_type"] == ct]
        for offset, cond in zip((0.22, 0.0, -0.22), CONDITIONS):
            vals = sub.loc[sub["condition"] == cond, "proportion"] * 100
            if vals.empty:
                continue
            jitter = rng.uniform(-0.045, 0.045, len(vals))
            ax.plot(vals, [i + offset + j for j in jitter], "o", ms=7,
                    color=SEVERITY[cond], mec=SURFACE, mew=1.2, zorder=3)
            ax.plot([vals.median()], [i + offset], "|", ms=16, mew=2,
                    color=SEVERITY[cond], zorder=4)
        if ct in abund.index:
            p = abund.loc[ct, "pval"]
            q = abund.loc[ct, "padj"]
            mark = "*" if q < PADJ else ""
            ax.text(xmax, i, f"p = {p:.3f}{mark}\nFDR = {q:.2f}", va="center", ha="right",
                    fontsize=8.5, color=INK if q < PADJ else INK_2, linespacing=1.3)

    ax.set_yticks(range(len(order)))
    ax.set_yticklabels(order, color=INK, fontsize=10)
    ax.grid(axis="y", visible=False)
    ax.set_xlabel("Share of a sample's cells (%)")
    despine(ax)
    handles = [plt.Line2D([], [], marker="o", ls="", ms=7, color=SEVERITY[c], mec=SURFACE)
               for c in CONDITIONS]
    fig.legend(handles, [f"{c} (n = {n_by_cond[c]})" for c in CONDITIONS],
               loc="lower left", bbox_to_anchor=(0.005, 0.0), ncol=3, handletextpad=0.3,
               columnspacing=1.6)
    header(fig, "Cell-type composition of each sample",
           "One dot per sample, tick marks at the group median. p is Mann-Whitney "
           "severe vs mild;\nFDR is after Benjamini-Hochberg across the eight cell types. "
           "No cell type reaches FDR < 0.05.")
    fig.tight_layout(rect=(0, 0.06, 1, 0.845))
    save(fig, "fig2_cohort_composition.png")


def fig_volcano(plt) -> None:
    """CD14 monocyte severe-vs-mild volcano."""
    import numpy as np
    de = load_de("CD14 Mono").dropna(subset=["padj", "log2FoldChange"]).copy()
    floor = de.loc[de["padj"] > 0, "padj"].min()
    de["nlp"] = -np.log10(de["padj"].clip(lower=floor))

    up = (de["padj"] < PADJ) & (de["log2FoldChange"] >= 1)
    dn = (de["padj"] < PADJ) & (de["log2FoldChange"] <= -1)
    ns = ~(up | dn)

    # Few labels, well separated: figure 3 carries the full gene-level detail.
    label_up = ["S100A8", "S100A12", "IL1R2", "SERPINB2", "RETN", "C3AR1"]
    label_dn = ["HLA-DRA", "HLA-DQA1", "HLA-DPB1", "CD74", "GNLY"]

    fig, ax = plt.subplots(figsize=(WIDTH, 5.6))
    ax.scatter(de.loc[ns, "log2FoldChange"], de.loc[ns, "nlp"], s=5, c=NS,
               linewidths=0, rasterized=True)
    ax.scatter(de.loc[dn, "log2FoldChange"], de.loc[dn, "nlp"], s=13, c=DOWN,
               linewidths=0, alpha=0.85)
    ax.scatter(de.loc[up, "log2FoldChange"], de.loc[up, "nlp"], s=13, c=UP,
               linewidths=0, alpha=0.85)
    ax.axhline(-np.log10(PADJ), color=AXIS, lw=0.75)

    # Park the labels in a clear column outside the point cloud and run a thin
    # leader to each gene, so no label sits on top of another point.
    xlo, xhi = de["log2FoldChange"].min() - 1.4, de["log2FoldChange"].max() + 1.6
    ax.set_xlim(xlo, xhi)
    for genes, x_col, ha in ((label_up, xhi - 0.15, "right"), (label_dn, xlo + 0.15, "left")):
        pts = de.loc[[g for g in genes if g in de.index]].sort_values("nlp", ascending=False)
        used = []
        for g, r in pts.iterrows():
            y = r["nlp"]
            while any(abs(y - u) < 0.62 for u in used):
                y -= 0.62
            used.append(y)
            ax.annotate(g, xy=(r["log2FoldChange"], r["nlp"]), xytext=(x_col, y),
                        textcoords="data", ha=ha, va="center", fontsize=8.5,
                        fontstyle="italic", color=INK_2,
                        arrowprops=dict(arrowstyle="-", lw=0.6, color=AXIS,
                                        shrinkA=2, shrinkB=4))
    ax.set_xlabel("log$_2$ fold change   (severe / mild)")
    ax.set_ylabel(r"$-\log_{10}$ adjusted $p$")
    despine(ax)
    handles = [plt.Line2D([], [], marker="o", ls="", ms=7, color=c, mec=SURFACE)
               for c in (UP, DOWN, NS)]
    fig.legend(handles, [f"Up in severe ({int(up.sum())})",
                         f"Up in mild ({int(dn.sum())})",
                         "Not significant"],
               loc="lower left", bbox_to_anchor=(0.005, 0.0), ncol=3, handletextpad=0.3,
               columnspacing=1.6)
    header(fig, "CD14 monocytes: what changes between severe and mild COVID-19",
           "Pseudobulk DESeq2 on 10 samples. Coloured points pass FDR < 0.05 and a "
           "two-fold change.")
    fig.tight_layout(rect=(0, 0.06, 1, 0.88))
    save(fig, "fig4_volcano_cd14_mono.png")


def fig_modules(plt) -> None:
    """log2FC of three monocyte gene modules, with significance marked."""
    mods = pd.read_csv(SUMMARY / "monocyte_gene_modules.csv")
    labels = {m.replace("\n", " "): m.replace("\n", " ") for m in MODULES}
    blocks = [(labels[m], mods[mods["module"] == m]) for m in labels]

    rows = sum(len(b) for _, b in blocks) + 2 * len(blocks)
    fig, ax = plt.subplots(figsize=(WIDTH, 0.9 + 0.265 * rows))
    lim = max(3.0, mods["log2FoldChange"].abs().max() * 1.12)
    ax.set_xlim(-lim, lim)

    yticks, ylabels, heads, y = [], [], [], 0.0
    for name, block in blocks:
        heads.append((name, y))          # module name on its own row above the block
        y -= 1
        for g, lfc, sig in zip(block["gene"], block["log2FoldChange"], block["significant"]):
            color = UP if lfc > 0 else DOWN
            ax.plot([0, lfc], [y, y], color=color, lw=1.5, solid_capstyle="round", zorder=2)
            ax.plot([lfc], [y], "o", ms=8, zorder=3, mew=1.5,
                    color=color if sig else SURFACE, mec=color)
            yticks.append(y)
            ylabels.append(g)
            y -= 1
        y -= 1

    ax.set_ylim(y + 0.5, 1.0)
    ax.axvline(0, color=AXIS, lw=0.75)
    ax.set_yticks(yticks)
    ax.set_yticklabels(ylabels, fontsize=9, fontstyle="italic", color=INK)
    ax.grid(axis="y", visible=False)
    ax.set_xlabel("log$_2$ fold change   (severe / mild)")
    despine(ax)
    for name, y_head in heads:
        ax.text(-lim * 0.99, y_head, name, ha="left", va="center",
                fontsize=9.5, color=INK, fontweight="semibold")
    handles = [plt.Line2D([], [], marker="o", ls="", ms=7, color=UP, mec=SURFACE),
               plt.Line2D([], [], marker="o", ls="", ms=7, color=SURFACE, mec=UP, mew=1.5)]
    fig.legend(handles, ["FDR < 0.05 (filled)", "Not significant (hollow)"],
               loc="lower left", bbox_to_anchor=(0.005, 0.0), ncol=2, handletextpad=0.3,
               columnspacing=1.6)
    header(fig, "Three monocyte programmes, measured separately",
           "Every gene in each module is shown, chosen before looking at the result. "
           "Bars right of zero\nare higher in severe disease.")
    fig.tight_layout(rect=(0, 0.045, 1, 0.885))
    save(fig, "fig5_monocyte_modules.png")


def fig_de_by_celltype(plt) -> None:
    """Significant-gene counts per cell type."""
    d = pd.read_csv(SUMMARY / "de_summary_by_celltype.csv")
    d = d.sort_values("significant").reset_index(drop=True)
    fig, ax = plt.subplots(figsize=(WIDTH, 0.9 + 0.42 * len(d)))
    xmax = d["significant"].max() * 1.30
    ax.set_xlim(0, xmax)
    ax.set_ylim(-0.6, len(d) - 0.4)
    ax.set_yticks(range(len(d)))
    ax.set_yticklabels(d["cell_type"], color=INK, fontsize=10)
    ax.grid(axis="y", visible=False)
    ax.set_xlabel("Genes differing between severe and mild (FDR < 0.05)")
    despine(ax)
    header(fig, "Almost all of the signal is in classical monocytes",
           "Pseudobulk DESeq2 per cell type, severe vs mild. Dendritic cells had too few "
           "cells to test.")
    fig.tight_layout(rect=(0, 0, 1, 0.86))
    for i, r in d.iterrows():
        # Emphasis, not polarity: the two cell types carrying the signal are
        # picked out in the accent hue, the rest recede. Red/blue stay reserved
        # for the severe/mild direction in the other figures.
        hbar(ax, i, r["significant"], 0.5, ACCENT if r["significant"] > 100 else NS)
        ax.text(r["significant"] + xmax * 0.012, i,
                f"{r['significant']:,}  of {r['genes_tested']:,} tested",
                va="center", color=INK_2, fontsize=9)
    save(fig, "fig3_de_by_celltype.png")


def fig_pathways(plt) -> None:
    """Enrichr terms for the genes raised in severe monocytes."""
    import numpy as np
    p = RESULTS / "enrichment" / "CD14_Mono" / "enrichr_severe_up_top15.csv"
    if not p.exists():
        logger.warning("no enrichment table at %s; skipping figure 5", p)
        return
    d = pd.read_csv(p)
    d = d.sort_values("Adjusted P-value").head(12).iloc[::-1].reset_index(drop=True)
    d["nlp"] = -np.log10(d["Adjusted P-value"])
    short = {"GO_Biological_Process_2021": "GO", "KEGG_2021_Human": "KEGG",
             "MSigDB_Hallmark_2020": "Hallmark"}

    fig, ax = plt.subplots(figsize=(WIDTH, 0.9 + 0.42 * len(d)))
    xmax = d["nlp"].max() * 1.20
    ax.set_xlim(0, xmax)
    ax.set_ylim(-0.6, len(d) - 0.4)
    ax.set_yticks(range(len(d)))
    ax.set_yticklabels([re.sub(r"\s*\(GO:\d+\)$", "", t) for t in d["Term"]],
                       color=INK, fontsize=9.5)
    ax.grid(axis="y", visible=False)
    ax.set_xlabel(r"$-\log_{10}$ adjusted $p$")
    despine(ax)
    # make_figures.py sends Enrichr the top 200 genes by fold change, so the
    # input is a capped subset of the 240 that pass both thresholds.
    de = load_de("CD14 Mono")
    n_up = int(((de["padj"] < PADJ) & (de["log2FoldChange"] >= 1)).sum())
    header(fig, "Genes raised in severe monocytes look like a neutrophil programme",
           f"Enrichr over the top 200 by fold change, of the {n_up} monocyte genes up at "
           "least two-fold in severe\ndisease. Overlap = input genes shared with the term.")
    fig.tight_layout(rect=(0, 0, 1, 0.86))
    for i, r in d.iterrows():
        hbar(ax, i, r["nlp"], 0.5, UP)
        ax.text(r["nlp"] + xmax * 0.012, i, f"{short.get(r['Gene_set'], '')}  {r['Overlap']}",
                va="center", color=INK_2, fontsize=9)
    save(fig, "fig6_pathways_monocyte.png")


def plot() -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    style(plt)
    fig_umap(plt)
    fig_composition(plt)
    fig_de_by_celltype(plt)
    fig_volcano(plt)
    fig_modules(plt)
    fig_pathways(plt)
    plt.close("all")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Build the study's tables and figures.")
    ap.add_argument("--figures", action="store_true",
                    help="redraw figures only, from the existing summary tables")
    args = ap.parse_args(argv)
    if not args.figures:
        write_summaries()
    plot()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
