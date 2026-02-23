# scEPI-Immune

### Abstract
Melanoma plasticity drives immune evasion and therapy resistance through dynamic cell-state transitions beyond genetic alterations. Although epigenetic remodeling is central to this process, its impact under therapeutic pressure remains unclear. We profiled longitudinal biopsies from melanoma patients treated in the phase Ib NIBIT-M4 epi-immunotherapy trial (NCT02608437, DNMT1 inhibitor plus anti-CTLA4) using single-cell multiome and spatial transcriptomics. Seven malignant meta-programs were identified, including a rare Wnt/β-catenin melanocytic state and a de-differentiated neural crest–like state enriched in non-responders. Spatial analyses showed that homotypic clustering stabilizes resistant programs, with neural crest–like cells forming compact niches. Responders displayed enrichment of antigen presentation/interferon program and coordinated T and B cell expansion, whereas non-responders retained stable neural crest–like clusters. Epigenetic therapy reactivated transposable elements, priming innate immunity and enhancing immunogenicity. NFATC2 emerged as a master regulator of neural crest–like states and resistance; its perturbation promoted differentiation and immunogenicity. These findings define mechanisms of resistance and nominate β-catenin and NFATC2 as therapeutic vulnerabilities.

----

### 10X Multiome GEX+ATAC
To characterize the regulatory basis of melanoma cell-state plasticity under therapy, we performed single-nucleus multiome sequencing to jointly profile gene expression (GEX) and chromatin accessibility (ATAC) from longitudinal melanoma biopsies (13 samples from 5 patients; R n=2, NR n=3), enabling integrated analysis of transcriptional programs and their underlying regulatory landscapes.

### Visium HD
To analyze how melanoma cells interact with their microenvironment and how cellular neighborhoods shape tumor heterogeneity and progression, we investigated the heterogeneity and spatial organization of metastatic melanoma lesions under therapy in nine samples from four patients (R n=2, and NR n=2) using 10x Genomics VisiumHD profiling.

----

###  `snRNA-seq/`
Contains all scripts for single-nucleus RNA-seq preprocessing, malignant state identification, meta-program extraction, functional characterization, regulatory inference, and cross-cohort validation.

#### `01-Individual_preprocessing.R`
Example of per-sample single-nucleus RNA-seq preprocessing, quality control, normalization, and initial clustering. Malignant cell annotation with SCEVAN.

#### `02-Merge_objects.R`
Merging of individual Seurat objects into a unified dataset for downstream analysis.

#### `03-TME_annotation.R`
Annotation of non-malignant compartments within the tumor microenvironment.

#### `04-Run_NMF.R`
Per-sample Non-negative Matrix Factorization (NMF) to extract malignant transcriptional programs.

#### `05-Metaprograms.R`
Identification of malignant transcriptional meta-programs.

#### `06-Metaprogram_assignment.R`
Per-cell meta-program scoring and assignment across samples.

#### `07-DEGs_memento.ipynb`
Differential expression analysis using Memento for state-specific contrasts.

#### `08-GSEA_DEGs.R`
Gene Set Enrichment Analysis of ranked differential expression results.

#### `09-CellChat.R`
Inference and comparison of ligand–receptor-mediated cell–cell communication networks.

#### `10-Melanoma_bulk.R`
Projection and validation of single-cell meta-programs in bulk melanoma transcriptomic datasets.

#### `11-Pozniak_et_al.R`
Cross-cohort validation and projection of meta-programs in independent melanoma datasets.

#### `12-GRN_Pozniak.R`
Preparation of expression matrices and transcription factor sets for gene regulatory network inference.

#### `13-Script_pySCENIC.sbatch`
SLURM submission script for pySCENIC-based gene regulatory network reconstruction.

#### `14-NFATC2.R`
Targeted analysis of NFATC2 regulatory activity and its association with resistant tumor states.

--

### `snATAC-seq/`
Contains scripts for single-cell chromatin accessibility preprocessing, peak calling, regulatory analysis, and integration with RNA-defined tumor meta-programs.

#### `01-Pre-processing.R`
Single-cell ATAC-seq preprocessing, quality control, dimensional reduction, and chromatin accessibility profiling.

#### `02-TME_annotation.R`
Annotation of microenvironmental compartments based on chromatin accessibility signatures.

#### `03-Metaprogram_assignment.R`
Projection and assignment of RNA-defined malignant meta-programs onto ATAC profiles.

#### `04-Peak_calling.R`
Genome-wide peak calling and construction of consensus accessible chromatin regions after cell type annotation.

#### `05-Motif_analysis.R`
Transcription factor motif enrichment and regulatory activity analysis across tumor states.

--

### `snRNA+snATAC/`
Contains integrative multiome analyses combining gene expression and chromatin accessibility to characterize tumor cell states, regulatory programs, and transcription factor activity across malignant cell states.

#### `01-Co-embedding_preprocessing.R`
Preprocessing and construction of joint RNA–ATAC co-embeddings for SCENIC+ input.

#### `02-SCENICplus_pyciTopic.py`
Execution of SCENIC+ workflow using pycisTopic for topic modeling and regulatory network inference.

#### `03-Create_cisTarget_db.sh`
Generation of custom cisTarget databases required for motif enrichment and regulon detection.

#### `04-snakemake_config.yaml`
Configuration file specifying parameters and resources for the SCENIC+ Snakemake pipeline.

#### `05-Save_SCENICplus_output.ipynb`
Processing and export of SCENIC+ output objects for downstream analyses.

#### `06-SCENICplus_downstream.R`
Downstream analysis of SCENIC+ regulons and transcription factor activity across tumor states.

#### `07-Overlap_SCENICplusEncode.R`
Overlap analysis between SCENIC+ regulatory regions and ENCODE reference datasets.

#### `08-Perturbation_regressors.py`
Construction of regression models for transcription factor perturbation simulations.

#### `09-perturbation_iteration.py`
Iterative in silico perturbation of regulatory networks to assess state transitions.

#### `10-perturbation_sbatch.sbatch`
SLURM submission script for perturbation simulations.

#### `11-MP4_perturbation.R`
Analysis and visualization of perturbation effects on the MP4 (neural crest–like) malignant meta-program.

--

### Visium HD analysis
This folder contains the complete workflow used to analyze high-resolution 10x Genomics VisiumHD melanoma samples and to integrate tumor metaprograms using Optimal Transport (OT).

## Notebooks

### `Annotations/`
For each sample, describe the annotation of cells in tumor vs non-tumor

### `1-Segmentation.ipynb`
Performs nuclei segmentation on high-resolution H&E images to obtain single-cell spatial coordinates. Outputs segmented cells and spatially indexed gene expression matrices.

### `2-OT_mapping.ipynb`
Applies Optimal Transport (OT) to transfer tumor metaprogram annotations and cell type labels from single-cell RNA-seq data to VisiumHD spatial cells based on shared gene-expression profiles.
Implements:
- Construction of gene-expression cost matrices  
- Entropic regularized OT  
- Transport-based metaprogram assignment 

### `3-Infiltration.ipynb`
Quantifies immune cell infiltration within tumor regions. Computes infiltration metrics and compares responders vs non-responders.

### `4-TSPS_correlations.ipynb`
Computes correlations between Tumor State Proximity Score (TSPS) between dispearsed and clustered cells within a metaprogram.

### `5-TSPS.ipynb`
Calculates the Tumor State Proximity Score (TSPS) for each metaprogram to quantify spatial clustering and domain compactness.

### `6-Spatial_Interactions.ipynb`
Analyzes heterotypic spatial interactions between malignant and non-malignant cells.  

### `Xenium_Pozniak/`
Validation analyses using independent Xenium spatial transcriptomics data.  
Includes OT-based transport of metaprograms and spatial validation of localization patterns.





