## `snATAC-seq/`
Contains scripts for single-cell chromatin accessibility preprocessing, peak calling, regulatory analysis, and integration with RNA-defined tumor meta-programs.

### `01-Peak_calling.R`
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
