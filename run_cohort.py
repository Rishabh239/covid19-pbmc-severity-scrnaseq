#!/usr/bin/env python
"""
run_cohort.py - End-to-end COVID-19 severe-vs-mild PBMC cohort workflow
=======================================================================

Chains the cohort extension into one run:

    load + label (GSE149689)  ->  Harmony integration + clustering
      ->  PBMC cell-type annotation  ->  differential abundance
      ->  pseudobulk DESeq2 (severe vs mild, per cell type)  ->  enrichment

Before running, fill ``SAMPLE_CONDITIONS`` in ``src/cohort.py`` from the
GSE149689 GEO "Samples" list (the loader refuses to proceed otherwise).

Example:
    python run_cohort.py --data data/GSE149689 --outdir results/covid_cohort

Place the three GEO supplementary files in --data, renamed to:
    barcodes.tsv.gz  features.tsv.gz  matrix.mtx.gz
"""

import argparse
import logging
import sys
from pathlib import Path

import pandas as pd

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger("run_cohort")


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="COVID-19 severe-vs-mild PBMC cohort workflow (GSE149689).",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    p.add_argument("--data", "-i", required=True,
                   help="Directory with combined 10x files (barcodes/features/matrix).")
    p.add_argument("--outdir", "-o", default="results/covid_cohort",
                   help="Output directory.")
    p.add_argument("--keep", nargs="+", default=["healthy", "mild", "severe"],
                   help="Conditions to retain (default: healthy mild severe).")
    p.add_argument("--test", default="severe", help="Test condition for DE.")
    p.add_argument("--ref", default="mild", help="Reference condition for DE.")
    p.add_argument("--resolution", type=float, default=1.0, help="Leiden resolution.")
    p.add_argument("--n-top-genes", type=int, default=2000, help="HVGs for integration.")
    p.add_argument("--integration", default="harmony",
                   choices=["harmony", "combat", "bbknn"],
                   help="Batch integration method (default: harmony; "
                        "combat needs no compiled deps).")
    p.add_argument("--min-cells", type=int, default=10,
                   help="Min cells per (sample x cell type) pseudobulk.")
    p.add_argument("--no-enrichment", action="store_true",
                   help="Skip enrichment (which needs internet for Enrichr).")
    return p


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)

    # Imported here so --help works even before the single-cell stack is installed.
    from src.cohort import load_cohort, per_sample_qc
    from src.integration import integrate_samples, annotate_pbmc
    from src.differential_abundance import cell_type_proportions, test_abundance
    from src.differential_state import run_all_celltypes, top_up_genes, enrich_up_genes
    from src.data_loader import save_h5ad

    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)

    logger.info("[1/6] Load + label cohort")
    adata = load_cohort(args.data, keep=tuple(args.keep))
    per_sample_qc(adata).to_csv(out / "per_sample_qc.csv", index=False)

    logger.info("[2/6] Harmony integration + clustering")
    adata = integrate_samples(
        adata, method=args.integration, resolution=args.resolution,
        n_top_genes=args.n_top_genes, results_dir=str(out / "integration"),
    )

    logger.info("[3/6] PBMC cell-type annotation")
    adata, mapping, scores = annotate_pbmc(adata)
    scores.to_csv(out / "celltype_scores.csv")
    pd.Series(mapping, name="cell_type").to_csv(out / "cluster_to_celltype.csv")
    save_h5ad(adata, out / "cohort_processed.h5ad")

    logger.info("[4/6] Differential abundance (%s vs %s)", args.test, args.ref)
    props = cell_type_proportions(adata)
    props.to_csv(out / "cell_type_proportions.csv", index=False)
    abund = test_abundance(props, group1=args.ref, group2=args.test)
    abund.to_csv(out / f"abundance_{args.test}_vs_{args.ref}.csv", index=False)

    logger.info("[5/6] Pseudobulk DESeq2 (%s vs %s)", args.test, args.ref)
    results, counts_df, meta_df = run_all_celltypes(
        adata, group_test=args.test, group_ref=args.ref,
        min_cells=args.min_cells, outdir=str(out / "pseudobulk"),
    )

    logger.info("[6/6] Enrichment on severe-up genes")
    if args.no_enrichment:
        logger.info("  Skipped (--no-enrichment).")
    else:
        for ct, res in results.items():
            genes = top_up_genes(res)
            if genes:
                safe = str(ct).replace(" ", "_").replace("/", "-")
                enrich_up_genes(genes, results_dir=str(out / "enrichment" / safe))

    logger.info("Done. Outputs written to %s", out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
