# scEPI-Immune

### Abstract
Melanoma plasticity drives immune evasion and therapy resistance through dynamic cell-state transitions beyond genetic alterations. Although epigenetic remodeling is central to this process, its impact under therapeutic pressure remains unclear. We profiled longitudinal biopsies from melanoma patients treated in the phase Ib NIBIT-M4 epi-immunotherapy trial (NCT02608437, DNMT1 inhibitor plus anti-CTLA4) using single-cell multiome and spatial transcriptomics. Seven malignant meta-programs were identified, including a rare Wnt/β-catenin melanocytic state and a de-differentiated neural crest–like state enriched in non-responders. Spatial analyses showed that homotypic clustering stabilizes resistant programs, with neural crest–like cells forming compact niches. Responders displayed enrichment of antigen presentation/interferon program and coordinated T and B cell expansion, whereas non-responders retained stable neural crest–like clusters. Epigenetic therapy reactivated transposable elements, priming innate immunity and enhancing immunogenicity. NFATC2 emerged as a master regulator of neural crest–like states and resistance; its perturbation promoted differentiation and immunogenicity. These findings define mechanisms of resistance and nominate β-catenin and NFATC2 as therapeutic vulnerabilities.

## 10X Multiome GEX+ATAC
Single-nucleus multiome sequencing was performed to jointly profile gene expression (GEX) and chromatin accessibility (ATAC) from longitudinal melanoma biopsies.

### snRNA

#### 01-Individual_preprocessing.R
Example of per-sample single-nucleus RNA-seq preprocessing, quality control, normalization, and initial clustering. Malignant cell annotation with SCEVAN.

#### 02-Merge_objects.R
Merging of individual Seurat objects into a unified dataset for downstream analysis.

#### 03-TME_annotation.R
Annotation of non-malignant compartments within the tumor microenvironment.

#### 04-Run_NMF.R
Per-sample Non-negative Matrix Factorization (NMF) to extract malignant transcriptional programs.

#### 05-Metaprograms.R
Identification of malignant transcriptional meta-programs.

#### 06-Metaprogram_assignment.R
Per-cell meta-program scoring and assignment across samples.

#### 07-DEGs_memento.ipynb
Differential expression analysis using Memento for state-specific contrasts.

#### 08-GSEA_DEGs.R
Gene Set Enrichment Analysis of ranked differential expression results.

#### 09-CellChat.R
Inference and comparison of ligand–receptor-mediated cell–cell communication networks.

#### 10-Melanoma_bulk.R
Projection and validation of single-cell meta-programs in bulk melanoma transcriptomic datasets.

#### 11-Pozniak_et_al.R
Cross-cohort validation and projection of meta-programs in independent melanoma datasets.

#### 12-GRN_Pozniak.R
Preparation of expression matrices and transcription factor sets for gene regulatory network inference.

#### 13-Script_pySCENIC.sbatch
SLURM submission script for pySCENIC-based gene regulatory network reconstruction.

#### 14-NFATC2.R
Targeted analysis of NFATC2 regulatory activity and its association with resistant tumor states.


### Visium HD
To analyze how melanoma cells interact with their microenvironment and how cellular neighborhoods shape tumor heterogeneity and progression, we investigated the heterogeneity and spatial organization of metastatic melanoma lesions under therapy in nine samples from four patients (R n=2, and NR n=2) using 10x Genomics VisiumHD profiling.

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





