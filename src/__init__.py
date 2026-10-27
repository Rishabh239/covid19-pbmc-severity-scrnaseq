"""
scRNA-seq Analysis Pipeline

A comprehensive pipeline for single-cell RNA sequencing analysis including
preprocessing, clustering, differential expression, pathway enrichment,
and machine learning classification.

Modules:
    - data_loader: Load 10x Genomics data
    - preprocessing: QC, normalization, clustering
    - differential_expression: Marker gene analysis
    - enrichment: GO and KEGG pathway analysis
    - classification: ML cell type classification
    - cli: Command-line interface

Author: Rishabh Tripathi
License: MIT
"""

__version__ = "1.0.0"

from .data_loader import load_10x_data, load_h5ad, save_h5ad
from .preprocessing import ScRNAPreprocessor, preprocess_adata
from .differential_expression import DifferentialExpression, find_markers
from .enrichment import EnrichmentAnalyzer, run_enrichment
from .classification import CellTypeClassifier

# Cohort (case/control) extension
from .cohort import load_cohort, label_cells, subset_conditions, SAMPLE_CONDITIONS
from .integration import integrate_samples, annotate_pbmc
from .differential_abundance import cell_type_proportions, test_abundance
from .differential_state import make_pseudobulk, run_all_celltypes, top_up_genes

__all__ = [
    'load_10x_data',
    'load_h5ad',
    'save_h5ad',
    'ScRNAPreprocessor',
    'preprocess_adata',
    'DifferentialExpression',
    'find_markers',
    'EnrichmentAnalyzer',
    'run_enrichment',
    'CellTypeClassifier',
    # cohort extension
    'load_cohort',
    'label_cells',
    'subset_conditions',
    'SAMPLE_CONDITIONS',
    'integrate_samples',
    'annotate_pbmc',
    'cell_type_proportions',
    'test_abundance',
    'make_pseudobulk',
    'run_all_celltypes',
    'top_up_genes',
]
