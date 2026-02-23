## `snATAC-seq/`
Contains scripts for single-cell chromatin accessibility preprocessing, peak calling, regulatory analysis, and integration with RNA-defined tumor meta-programs.

### `01-Pre-processing.R`
- Builds a consensus peak set across samples from per-sample BED files.
- Quantifies chromatin accessibility using fragment files.
- Merges all samples into a single Signac/Seurat object with batch, patient, and week metadata.
- Performs QC filtering (TSS enrichment, nucleosome signal, FRiP, blacklist ratio).
- Runs the standard scATAC workflow (TF-IDF, LSI, clustering, and t-SNE visualization).

### `02-TME_annotation.R`
- Imports snRNA-seq–derived cell-type annotations and transfers them to the scATAC object by barcode matching.
- Computes gene activity scores and performs label transfer from the snRNA-seq reference to predict cell identities in ATAC.
- Refines tumor microenvironment (TME) cell-type annotations using automated references (SingleR/celldex) and ATAC clustering.
- Generates coverage and marker accessibility plots for key lineage and state-defining loci.

### `03-Metaprogram_assignment.R`
- Transfers RNA-derived malignant meta-program labels to scATAC malignant cells via label transfer (CCA anchors) using gene activity profiles.
- Performs malignant-only snATAC dimensional reduction and clustering (TF-IDF/LSI) for state-level visualization.
- Stores per-cell meta-program assignments in the combined ATAC object for downstream regulatory analyses.

### `04-Peak_calling.R`
- Performs MACS2 peak calling on the annotated snATAC dataset using fine-grained groups (cell types and malignant meta-programs).
- Builds a new peak-by-cell matrix from the called peaks and re-runs the standard snATAC workflow (TF-IDF/LSI, clustering, t-SNE).
- Generates dedicated TME and malignant subsets for downstream accessibility and state-specific analyses.

### `05-Motif_analysis.R`
- Adds TF motif annotations to malignant snATAC peaks (HOCOMOCO v13/JASPAR-format motifs) and computes motif enrichment in differentially accessible peaks across meta-programs.
- Builds TF-by-metaprogram enrichment matrices and performs TF footprinting for selected regulators.
- Integrates ENCODE cCRE annotations (PLS/pELS/dELS) to link motif-bearing peaks to target genes and runs pathway enrichment on linked genes.
