## `snRNA-seq/`
Contains scripts for single-nucleus RNA-seq preprocessing, malignant cell identification, tumor microenvironment annotation, extraction of malignant transcriptional programs, meta-program definition and assignment, functional characterization, regulatory network inference, and validation across external single-cell and bulk melanoma cohorts.

### `01-Individual_preprocessing.R`
- Example of per-sample snRNA-seq preprocessing, including QC filtering, normalization, dimensional reduction, and clustering.
- Identifies and annotates malignant cells using SCEVAN-based CNV inference.

### `02-Merge_objects.R`
- Merges per-sample Seurat objects into a single dataset for downstream malignant/TME analyses.

### `03-TME_annotation.R`
- Annotates non-malignant tumor microenvironment compartments within the merged dataset.
- Produces curated cell-type labels used for downstream stratified analyses and visualization.

### `04-Run_NMF.R`
- Runs per-sample Non-negative Matrix Factorization (NMF) on malignant cells to extract transcriptional programs.
- Exports program loadings and gene weights for meta-program consolidation.

### `05-Metaprograms.R`
- Consolidates per-sample NMF programs into consensus malignant meta-programs shared across samples.
- Generates meta-program gene signatures used throughout downstream analyses.

### `06-Metaprogram_assignment.R`
- Scores individual malignant cells for each meta-program signature and assigns a dominant meta-program per cell.
- Produces the meta-program annotations used for state-resolved differential and regulatory analyses.

### `07-DEGs_memento.ipynb`
- Performs differential expression analyses to quantify state-specific expression changes.
- Outputs ranked gene lists and summary statistics for enrichment testing and interpretation.

### `08-GSEA_DEGs.R`
- Runs Gene Set Enrichment Analysis (GSEA) on ranked differential expression results.
- Summarizes enriched pathways and functional themes across malignant states/meta-programs.

### `09-CellChat.R`
- Infers ligand–receptor-mediated cell–cell communication networks using CellChat.
- Compares signaling interactions across conditions/states and identifies key sender/receiver cell types.

### `10-Melanoma_bulk.R`
- Validates meta-program signatures in external bulk melanoma transcriptomic datasets.
- Assesses clinical relevance of malignant states identified.

### `11-Pozniak_et_al.R`
- Validates and projects meta-program signatures in an independent melanoma single-cell cohort (Pozniak et al., 2024).
- Quantifies cross-cohort reproducibility of malignant states and associated transcriptional programs.

### `12-GRN_Pozniak.R`
- Prepares expression matrices, TF lists, and required inputs for gene regulatory network (GRN) inference in the external cohort.
- Exports formatted files for pySCENIC and downstream regulon analysis.

### `13-Script_pySCENIC.sbatch`
- SLURM script to run the pySCENIC pipeline (GRN inference, motif enrichment, regulon scoring).
  
### `14-NFATC2.R`
- Performs targeted analysis of NFATC2 regulatory activity across malignant states.
- Links NFATC2-associated programs/regulons to resistant tumor phenotypes and state-specific signatures.
