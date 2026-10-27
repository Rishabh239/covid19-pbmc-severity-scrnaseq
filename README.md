# What separates severe from mild COVID-19 in blood immune cells

A single-cell RNA-seq reanalysis of **GSE149689** — 37,088 PBMCs from 14 donors
(4 healthy, 4 mild COVID-19, 6 severe COVID-19). The question: **when one virus
leaves most people with a mild illness but sends a minority into critical care,
what does the difference look like in circulating immune cells?**

This is an *in silico* reanalysis of public data. No wet-lab work was done here,
and the dataset is not mine — it was published by
[Lee et al., *Science Immunology* 2020](https://doi.org/10.1126/sciimmunol.abd1554).
What is mine is the workflow, the analysis choices, and the write-up below.

Every number here comes from a file in [`results/`](results/), rebuilt by
[`make_report.py`](make_report.py).

## Summary

**Classical monocytes carry almost all of the signal.** Of the eight cell types
tested, CD14 monocytes gave 588 genes differing between severe and mild disease
at FDR < 0.05; the next-largest was platelets at 149, and four cell types gave
fewer than five.

Within those monocytes, two programmes move in opposite directions:

| Programme | Direction in severe | Significant? |
|---|---|---|
| MHC class II / antigen presentation (*HLA-DRA*, *CD74*, …) | **down**, up to 12-fold | yes, all 7 genes at FDR < 0.05 |
| Alarmins & emergency myelopoiesis (*S100A8/9/12*, *IL1R2*, *RETN*, …) | **up**, up to 37-fold | yes, all 9 genes at FDR < 0.05 |
| Type I interferon response (*ISG15*, *IFI6*, *MX1*, …) | up | **no** — all 7 genes FDR > 0.29 |

That combination — an antigen-presentation-low, alarmin-high monocyte — is the
phenotype described in severe COVID-19 and in bacterial sepsis, where it is read
as a marker of a dysfunctional innate response rather than a stronger one.

**Three results worth stating plainly:**

1. **No shift in cell-type composition survives multiple-testing correction.**
   Dendritic cells and CD8 T cells are both lower in severe disease
   (p = 0.014 and 0.019, Mann-Whitney), but with 4 mild against 6 severe donors,
   neither clears FDR < 0.05 across eight cell types. The direction matches the
   literature; this cohort is too small to call it.
2. **The type I interferon signature did not reproduce.** All seven interferon-
   stimulated genes trend upward in severe disease, and none reaches
   significance. This was the headline of the original study, so its absence
   here is a genuine negative result, discussed in [section 3](#3-discussion).
3. **The genes raised in severe monocytes are a neutrophil programme.**
   Neutrophil degranulation and neutrophil-mediated immunity dominate the
   enrichment at FDR ≈ 2 × 10⁻²², alongside complement, coagulation and
   neutrophil extracellular trap formation.

![Cell-type composition of each sample](figures/fig2_cohort_composition.png)

## 1. Data and methods

### Dataset

| | |
|---|---|
| Accession | [GSE149689](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE149689) (Lee et al., *Sci Immunol* 2020) |
| Assay | 10x Genomics Chromium 3′ v3, PBMCs |
| Donors used | 14 — healthy 4, mild COVID-19 4, severe COVID-19 6 |
| Donors excluded | 1 asymptomatic and 5 severe influenza, to keep a clean severity contrast |
| Cells | 59,963 loaded → **37,088 after QC**, across 22,104 genes |
| Severity labels | from the paper's Supplementary Table S1; GEO titles only distinguish nCoV / Flu / Normal |

Cells per group after QC: healthy 14,495, mild 14,085, severe 8,508. The severe
group has more donors but fewer cells, which matters for the pseudobulk step
below.

### Workflow

```
load + label (GSE149689)  ->  Harmony integration over donor  ->  Leiden clustering
   ->  marker-based PBMC annotation  ->  differential abundance (severe vs mild)
   ->  pseudobulk DESeq2 per cell type  ->  Enrichr pathway enrichment
```

| Step | Module | Method |
|---|---|---|
| Load and label | `src/cohort.py` | per-barcode donor and condition from the 10x suffix |
| Integrate | `src/integration.py` | Harmony over `sample_id`, 2,000 HVGs |
| Cluster and annotate | `src/integration.py` | Leiden (resolution 1.0) → 26 clusters → 8 PBMC types by canonical marker score |
| Composition | `src/differential_abundance.py` | per-donor proportions, Mann-Whitney, Benjamini-Hochberg across cell types |
| Expression | `src/differential_state.py` | pseudobulk per donor × cell type, DESeq2 (pydeseq2), severe vs mild |
| Pathways | `src/enrichment.py`, `make_figures.py` | Enrichr over genes up ≥ 2-fold in severe |

Run as:

```bash
python run_cohort.py --data data/GSE149689 --outdir results
```

### Why pseudobulk rather than per-cell testing

Treating thousands of cells from one donor as independent replicates inflates
significance, because cells from the same person are not independent
observations — the pseudoreplication problem. Summing counts to one profile per
donor × cell type makes the **donor** the unit of replication, which is current
best practice for case/control single-cell designs. The cost is power: the
monocyte comparison rests on 4 mild against 6 severe donors, not thousands of
cells.

## 2. Results

### 2.1 Cell types

![Integrated UMAP by cell type and condition](figures/fig1_umap.png)

Leiden clustering of the Harmony-integrated data gave 26 clusters, labelled by
canonical marker score into 8 PBMC populations:

| Cell type | Cells | Cell type | Cells |
|---|---:|---|---:|
| CD14 monocytes | 10,221 | B | 4,547 |
| NK | 6,956 | Platelet | 3,190 |
| CD8 T | 5,296 | CD16 monocytes | 1,559 |
| CD4 T | 4,843 | Dendritic cells | 476 |

Annotation is a marker-score first pass, not a curated reference mapping — see
[Limitations](#4-limitations).

### 2.2 Composition shifts do not survive correction

Figure 1 shows every donor's cell-type composition. The largest effects are a
loss of dendritic cells (log₂FC −3.42) and of CD8 T cells (−1.60) in severe
disease, both nominally significant:

| Cell type | log₂FC severe/mild | p | FDR |
|---|---:|---:|---:|
| Dendritic cells | −3.42 | 0.014 | 0.076 |
| CD8 T | −1.60 | 0.019 | 0.076 |
| CD16 monocytes | −2.33 | 0.114 | 0.229 |
| Platelet | +1.71 | 0.067 | 0.178 |
| CD14 monocytes | +0.57 | 0.352 | 0.403 |

Nothing reaches FDR < 0.05. Dendritic-cell loss and T-cell lymphopenia are both
well described in severe COVID-19, and the directions here match — but a
Mann-Whitney test on 4 against 6 donors cannot establish them, and reporting
them as established would overstate what this cohort shows.

### 2.3 Expression changes concentrate in classical monocytes

![Genes differing between severe and mild, per cell type](figures/fig3_de_by_celltype.png)

CD14 monocytes gave 588 genes at FDR < 0.05 (353 up in severe, 235 up in mild).
Platelets gave 149, overwhelmingly down in severe — but platelet transcriptomes
in PBMC preparations are heavily shaped by ambient RNA and sample handling, so I
have not interpreted them. Dendritic cells could not be tested: only one severe
donor had enough cells to form a pseudobulk profile.

![CD14 monocyte volcano](figures/fig4_volcano_cd14_mono.png)

### 2.4 Two programmes move, a third does not

Three gene modules were chosen *before* inspecting the results, each a published
axis of the severe-COVID monocyte phenotype, and every gene in each is shown:

![Three monocyte gene modules](figures/fig5_monocyte_modules.png)

- **MHC class II falls.** All seven genes drop significantly: *HLA-DQA1* −3.57,
  *HLA-DPB1* −2.60, *HLA-DPA1* −2.29, *HLA-DRA* −1.91, *CD74* −1.75,
  *HLA-DMA* −1.50, *HLA-DMB* −1.48 (all FDR < 0.05). Reduced monocyte HLA-DR is
  one of the most reproducible markers of severe COVID-19 and of sepsis.
- **Alarmins and emergency-myelopoiesis genes rise.** All nine go up
  significantly: *IL1R2* +5.21, *SERPINB2* +4.77, *THBS1* +4.23, *HP* +4.12,
  *CLU* +3.97, *RETN* +3.18, *S100A8* +2.58, *S100A12* +1.82, *S100A9* +1.63.
  *RETN*, *HP*, *IL1R2* and *S100A8/9/12* together describe the immature,
  immunosuppressive monocyte state reported in sepsis.
- **Type I interferon does not separate the groups.** All seven interferon-
  stimulated genes trend up — *ISG15* +1.80, *IFITM3* +1.30, *IFI6* +1.29 — and
  none is significant (smallest FDR 0.29, for *OAS1*).

*TNF* and *IL1B* both trend **downward** in severe monocytes (−1.57 and −1.60,
FDR 0.23 and 0.53), the opposite of the cytokine-storm picture, though neither
is significant either.

### 2.5 The severe monocyte signature is a neutrophil programme

![Enriched pathways](figures/fig6_pathways_monocyte.png)

Neutrophil degranulation, neutrophil-mediated immunity and neutrophil activation
all return FDR ≈ 2 × 10⁻²², with 40 of the input genes in each term. Complement,
coagulation and neutrophil extracellular trap formation follow. The monocytes of
severely ill patients are expressing a granulocyte-like programme — consistent
with the release of immature myeloid cells from the marrow under emergency
myelopoiesis, which is the standard reading of this signature.

## 3. Discussion

The monocyte result is strong and internally consistent: antigen presentation
down, alarmins up, a neutrophil-like programme on top, all from a design that
treats the donor as the replicate. It reproduces a phenotype reported
independently in severe COVID-19 and in bacterial sepsis, which is reassuring
because nothing in this workflow was tuned to find it.

**The interferon result deserves more attention than the positive findings.**
The original study's headline was a type I interferon response distinguishing
severe COVID-19. This reanalysis does not reproduce it as significant. Several
explanations are compatible with the data, and I cannot separate them here:

- **Power.** Six severe donors against four mild, with severe pseudobulks built
  from a median of 299 cells against 895 for mild. Interferon-stimulated genes
  are also highly variable between donors, which is exactly the case where small
  n fails.
- **Timing.** The interferon response peaks early and decays; the severity
  contrast depends on when each donor was sampled relative to symptom onset.
  Days-since-onset is in the paper's Table S1 but is not modelled here.
- **Contrast.** This compares severe against *mild COVID-19*. The original
  analysis makes more of the comparison against healthy donors and against
  influenza, where interferon differences are larger.
- **Annotation.** A coarse marker-score label for "CD14 Mono" merges states that
  a finer annotation would separate, diluting a signal confined to a subset.

The honest conclusion is that this cohort supports the monocyte dysfunction
axis and is **underpowered for the interferon axis** — not that the interferon
finding is wrong.

## 4. Limitations

- **Small cohort.** 4 mild against 6 severe donors for every comparison. The
  composition tests in particular are underpowered, and none survives FDR
  correction.
- **Unbalanced depth.** Severe samples yielded far fewer cells (8,508 across 6
  donors) than mild (14,085 across 4). Severe monocyte pseudobulks rest on as
  few as 141 cells.
- **Coarse annotation.** Cell types come from a marker-score pass over Leiden
  clusters, not a curated reference. Monocyte subsets, T-cell states and
  plasmablasts are not resolved, and clusters were not manually reviewed
  against marker plots.
- **No covariates.** Age, sex, days since symptom onset and treatment are all in
  the source metadata and none is modelled, despite each plausibly affecting
  these genes.
- **Composition testing is unsophisticated.** A Mann-Whitney test on proportions
  ignores their compositional nature; scCODA or Milo would be more appropriate.
- **Platelet results are not interpreted.** Ambient RNA makes them unreliable in
  PBMC data.
- **Enrichment is one-sided and capped.** Only genes up in severe were tested,
  and only the top 200 by fold change of the 240 qualifying genes were submitted
  to Enrichr. Genes down in severe were not tested for enrichment.
- **The healthy group is loaded but barely used.** All reported contrasts are
  severe vs mild; healthy donors appear in Figure 1 and in the integration only.

## 5. Next steps

1. Model days-since-onset and age as covariates in the DESeq2 design — the most
   likely explanation for the missing interferon signal.
2. Add the severe-vs-healthy and mild-vs-healthy contrasts, and the influenza
   donors as a different-virus comparison.
3. Replace marker-score annotation with reference mapping against a curated PBMC
   atlas, then re-test within monocyte subsets.
4. Re-test composition with scCODA, which handles the compositional constraint
   properly.
5. Score interferon activity per cell as a module score rather than gene by
   gene, which is better powered than seven separate tests.

## Reproduce

```bash
pip install -r requirements.txt

# 1. Download the GSE149689 supplementary matrix into data/GSE149689/, renamed to
#    barcodes.tsv.gz, features.tsv.gz, matrix.mtx.gz  (see docs/workflow.md)
# 2. Run the analysis (~30-60 min on a laptop for 37k cells)
python run_cohort.py --data data/GSE149689 --outdir results

# 3. Rebuild the summary tables and figures in this README
python make_report.py
```

`make_report.py --figures` redraws the figures from the committed tables without
re-running anything. Full setup notes, the sample-to-condition mapping and
Windows-specific caveats are in [`docs/workflow.md`](docs/workflow.md).

## Files

| Path | Contents |
|---|---|
| `results/per_sample_qc.csv` | cells and median genes per donor, as loaded |
| `results/cell_type_proportions.csv` | per-donor cell-type counts and proportions |
| `results/abundance_severe_vs_mild.csv` | composition test, with BH-corrected FDR |
| `results/cluster_to_celltype.csv`, `celltype_scores.csv` | cluster annotation and the marker scores behind it |
| `results/pseudobulk/deseq2_<cell type>_severe_vs_mild.csv` | full DESeq2 results, 7 cell types |
| `results/pseudobulk/pseudobulk_counts.csv`, `_metadata.csv` | the donor × cell-type count matrix fed to DESeq2 |
| `results/enrichment/<cell type>/` | Enrichr output per cell type |
| `results/summary/` | derived tables behind this README (DE counts, significant genes, gene modules, cohort totals) |
| `results/summary/umap_coords.csv.gz` | UMAP coordinates and labels per cell, so the UMAP redraws without the 527 MB object |
| `results/qc/` | UMAPs by cell type, condition, donor and cluster, plus QC violins |
| `figures/` | the six figures above |
| `src/`, `run_cohort.py` | the analysis code |

The processed AnnData object (`cohort_processed.h5ad`, 527 MB) is **not** in this
repository — it exceeds GitHub's file limit. `run_cohort.py` regenerates it from
the GEO matrix.

## Licence and disclaimer

Code is MIT-licensed (see [LICENSE](LICENSE)). The GSE149689 data belongs to its
original authors and is used here under its GEO terms; please cite
[Lee et al. 2020](https://doi.org/10.1126/sciimmunol.abd1554) for the data.

Research and portfolio use only. These results are a reanalysis by one person
and have not been peer reviewed. Nothing here is clinical guidance.
