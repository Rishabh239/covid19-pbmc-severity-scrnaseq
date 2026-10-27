# Running the workflow

How to reproduce the analysis in [`../README.md`](../README.md) from the raw GEO
matrix. For what the results mean, read the README; this file is setup and
operation only.

## 1. Get the data

**GSE149689** — Lee et al., *Science Immunology* 2020, "Immunophenotyping of
COVID-19 and influenza highlights the role of type I interferons in development
of severe COVID-19". Open access on GEO (CC BY 4.0), no data-access application
needed.

The series ships one combined 10x matrix holding all 20 samples; each barcode
ends in a sample index (`-1` … `-20`).

1. Open the
   [GSE149689 GEO page](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE149689)
   and download the supplementary single-cell matrix (barcodes, features, matrix).
2. Put the three files in `data/GSE149689/` and rename them to exactly:

   ```
   barcodes.tsv.gz   features.tsv.gz   matrix.mtx.gz
   ```

   Strip the `GSE149689_` prefix. If the features file is named `genes.tsv.gz`,
   rename it to `features.tsv.gz`.

## 2. The sample-to-condition map

This is already done for GSE149689 and needs no action.

GEO's sample titles only say nCoV / Flu / Normal, so the COVID **severities**
were read from the paper's Supplementary Table S1 (each patient's clinical
"Disease group") and written into `SAMPLE_CONDITIONS` in
[`../src/cohort.py`](../src/cohort.py). It validates to 6 severe, 4 mild,
1 asymptomatic, 5 influenza, 4 healthy.

If you point the loader at a different matrix, re-derive the map from that
dataset's GEO order. `validate_sample_map()` refuses to run until the group
counts are consistent, and `summarize_suffixes(adata)` prints cells per suffix
so you can confirm the sample count.

## 3. Install

```bash
pip install -r requirements.txt
```

**Windows note.** `requirements.txt` pins `harmonypy<2` deliberately. harmonypy
2.x uses a C++/Armadillo backend with pre-built wheels only for Linux and macOS,
so on Windows pip tries to compile it and fails looking for BLAS. The pinned
pure-Python build installs cleanly and behaves the same.

If Harmony still will not install, use an integrator that needs no compilation:

```bash
python run_cohort.py --data data/GSE149689 --integration combat
```

`--integration` accepts `harmony` (default), `combat` (built into scanpy) or
`bbknn`. If `harmony` is requested but harmonypy cannot be imported, the
pipeline falls back to `combat` automatically.

## 4. Run

```bash
python run_cohort.py --data data/GSE149689 --outdir results
```

The published results used the defaults: Harmony integration, Leiden resolution
1.0, 2,000 highly variable genes, minimum 10 cells per pseudobulk, severe vs
mild contrast, healthy/mild/severe retained. Expect roughly 30–60 minutes for
37k cells on a laptop.

Useful flags:

| Flag | Default | Effect |
|---|---|---|
| `--keep` | `healthy mild severe` | conditions to retain |
| `--test` / `--ref` | `severe` / `mild` | the DE contrast |
| `--resolution` | `1.0` | Leiden clustering granularity |
| `--n-top-genes` | `2000` | highly variable genes used for integration |
| `--min-cells` | `10` | minimum cells per donor × cell type pseudobulk |
| `--integration` | `harmony` | `harmony`, `combat` or `bbknn` |
| `--no-enrichment` | off | skip Enrichr, which needs internet |

## 5. Rebuild the report

```bash
python make_report.py              # summary tables + the five README figures
python make_report.py --figures    # figures only, from the committed tables
```

`make_report.py` reads `results/` and writes `results/summary/` and `figures/`.
It re-runs nothing, so it is safe to use on the committed results without the
raw data present.

The monocyte pathway figure (figure 6) is the one part not built from
`run_cohort.py` output alone. It reads an Enrichr table produced by
[`../make_figures.py`](../make_figures.py), committed at
`results/enrichment/CD14_Mono/enrichr_severe_up_top15.csv`. To regenerate it:

```bash
python make_figures.py \
    --de results/pseudobulk/deseq2_CD14_Mono_severe_vs_mild.csv \
    --pseudobulk results/pseudobulk/pseudobulk_counts.csv \
    --meta results/pseudobulk/pseudobulk_metadata.csv \
    --celltype "CD14 Mono" --test severe --ref mild --outdir figures
```

That writes `figures/pathways_CD14_Mono.csv` (plus a volcano and heatmap of its
own, separate from the README figures). Copy the CSV over the committed table to
refresh it:

```bash
cp figures/pathways_CD14_Mono.csv results/enrichment/CD14_Mono/enrichr_severe_up_top15.csv
```

Enrichr receives only the top genes by fold change — 200 by default, set with
`--top-genes`. That cap is why figure 6 says "top 200 of 240".

## 6. What each step writes

| Step | Module | Output under `--outdir` |
|---|---|---|
| Load and label | `src/cohort.py` | `per_sample_qc.csv` |
| Integrate and cluster | `src/integration.py` | `integration/qc/*.png`, UMAP embeddings |
| Annotate | `src/integration.py` | `cluster_to_celltype.csv`, `celltype_scores.csv`, `cohort_processed.h5ad` |
| Differential abundance | `src/differential_abundance.py` | `cell_type_proportions.csv`, `abundance_severe_vs_mild.csv` |
| Pseudobulk DE | `src/differential_state.py` | `pseudobulk/deseq2_<cell type>_severe_vs_mild.csv`, `pseudobulk_counts.csv`, `pseudobulk_metadata.csv` |
| Enrichment | `src/enrichment.py` | `enrichment/<cell type>/GO_severe_up.csv`, `KEGG_severe_up.csv` |

`cohort_processed.h5ad` is about 500 MB and is not committed — it exceeds
GitHub's file size limit. Everything else in `results/` is in the repository.

## 7. Notes and caveats

- **Annotation is coarse.** Cell types come from a marker-score pass over Leiden
  clusters. Review the labels against the marker plots in `results/qc/` before
  trusting anything downstream.
- **Differential abundance** uses a nonparametric test on a small number of
  donors. Compositional models (scCODA) or neighbourhood tests (Milo) are more
  rigorous.
- **Enrichment** goes through gseapy's Enrichr backend, which needs internet.
- **pydeseq2** needs at least 2 donors per group per cell type and enough cells
  per pseudobulk (`--min-cells`). Dendritic cells fail this in the severe group,
  so they are not tested.
