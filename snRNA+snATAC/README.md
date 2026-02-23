## `snRNA+snATAC/`
Contains integrative multiome analyses combining gene expression and chromatin accessibility to characterize tumor cell states, regulatory programs, and transcription factor activity across malignant cell states.

### `01-Co-embedding_preprocessing.R`
- Builds per-sample multiome Seurat objects (RNA + ATAC) from 10x outputs, merges them, and applies joint RNA/ATAC QC filtering.
- Integrates curated metadata (cell-type and metaprogram annotations) from snRNA-seq and snATAC-seq analyses into the multiome object.
- Computes RNA-only, ATAC-only, and weighted multi-modal (WNN) embeddings, and saves malignant and TME multiome subsets for downstream co-embedding and regulatory analyses.

### `02-SCENICplus_pyciTopic.py`
- Creates a pycisTopic object from a peak-by-cell count matrix and associated annotations/metadata.
- Runs MALLET LDA topic modeling across multiple topic numbers, selects the optimal model, and stores it in the cisTopic object.
- Binarizes topics and identifies candidate enhancer regions and differentially accessible regions (DARs) by meta-program.
- Exports topic/DAR region sets as BED files for downstream SCENIC+ regulatory network inference.

