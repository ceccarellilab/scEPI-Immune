## `snRNA_snATAC/`
Contains integrative multiome analyses combining gene expression and chromatin accessibility to characterize tumor cell states, regulatory programs, and transcription factor activity across malignant cell states.

### `01-Co-embedding_preprocessing.R`
- Builds per-sample multiome Seurat objects (RNA + ATAC) from 10x outputs, merges them, and applies joint RNA/ATAC QC filtering.
- Integrates curated metadata (cell-type and metaprogram annotations) from snRNA-seq and snATAC-seq analyses into the multiome object.
- Computes RNA-only, ATAC-only, and weighted multi-modal (WNN) embeddings, and saves malignant and TME multiome subsets for downstream co-embedding and regulatory analyses.

### `02-SCENICplus_pyciTopic.py`
- Creates a pycisTopic object from a peak-by-cell count matrix and associated annotations/metadata.
- Performs topic modeling across multiple topic numbers, selects the optimal model, and stores it in the cisTopic object.
- Binarizes topics and identifies candidate enhancer regions and differentially accessible regions (DARs) by meta-program.
- Exports topic/DAR region sets as BED files for downstream SCENIC+ regulatory network inference.

### `03-Create_cisTarget_db.sh`
- Generates custom cisTarget motif-ranking databases from consensus peak regions by creating padded FASTA sequences and scanning them against curated motif collections.
- Builds SCENIC+-compatible motif databases for downstream cis-regulatory enrichment and regulon inference.

### `04-snakemake_config.yaml`
- Defines input/output paths and analysis parameters for the SCENIC+ Snakemake workflow.

### `05-Save_SCENICplus_output.ipynb`
- Organizes and exports regulons, motif enrichments, and regulatory network objects for downstream interpretation and visualization.

### `06-SCENICplus_downstream.R`
- Imports SCENIC+ eRegulon metadata, AUCell matrices (gene- and region-based), and RSS specificity scores.
- Adds eRegulon AUCell scores as new Seurat assays and computes scaled meta-program–level averages for genes, regions, and TF expression.
- Ranks and prioritizes transcription factors per meta-program based on differential TF expression, gene-based AUCell, and region-based AUCell activity.
- Generates integrative heatmap visualizations of regulon activity (AUC), specificity (RSS), regulon size, and overlap with meta-program gene signatures.

### `07-Overlap_SCENICplusEncode.R`
- Intersects SCENIC+ eRegulon regions from prioritized meta-program–specific TFs with ENCODE hg38 cCRE annotations (pELS and dELS) using `GRanges`.
- Quantifies the percentage of regulon regions overlapping protein-coding–linked cCREs per TF and meta-program, and exports a consolidated overlap table.

### `08-Perturbation_regressors.py`
- Trains SCENIC+ gene expression regressors from scRNA counts using TF–target links derived from eRegulon metadata.
- Saves fitted models for downstream in silico perturbation simulations.

### `09-perturbation_iteration.py`
- Performs in silico TF knockdown simulations using pre-trained regressors and exports per-TF perturbation results.

### `10-perturbation_sbatch.sbatch`
- SLURM submission script for parallel transcription factor perturbation analyses.

### `11-MP4_perturbation.R`
- Aggregates TF perturbation outputs across iterations and restricts analysis to cells assigned to the MP_4 (neural crest-like) meta-program.
- Computes per-TF mean log2 fold-change of MP_4 signature genes relative to baseline (iteration 0) to quantify perturbation impact.
