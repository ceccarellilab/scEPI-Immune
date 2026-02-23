# scEPI-Immune

### Abstract
Melanoma plasticity drives immune evasion and therapy resistance through dynamic cell-state transitions beyond genetic alterations. Although epigenetic remodeling is central to this process, its impact under therapeutic pressure remains unclear. We profiled longitudinal biopsies from melanoma patients treated in the phase Ib NIBIT-M4 epi-immunotherapy trial (NCT02608437, DNMT1 inhibitor plus anti-CTLA4) using single-cell multiome and spatial transcriptomics. Seven malignant meta-programs were identified, including a rare Wnt/β-catenin melanocytic state and a de-differentiated neural crest–like state enriched in non-responders. Spatial analyses showed that homotypic clustering stabilizes resistant programs, with neural crest–like cells forming compact niches. Responders displayed enrichment of antigen presentation/interferon program and coordinated T and B cell expansion, whereas non-responders retained stable neural crest–like clusters. Epigenetic therapy reactivated transposable elements, priming innate immunity and enhancing immunogenicity. NFATC2 emerged as a master regulator of neural crest–like states and resistance; its perturbation promoted differentiation and immunogenicity. These findings define mechanisms of resistance and nominate β-catenin and NFATC2 as therapeutic vulnerabilities.

***

### 10X Multiome GEX+ATAC
To characterize the regulatory basis of melanoma cell-state plasticity under therapy, we performed single-nucleus multiome sequencing to jointly profile gene expression (GEX) and chromatin accessibility (ATAC) from longitudinal melanoma biopsies (13 samples from 5 patients; R n=2, NR n=3), enabling integrated analysis of transcriptional programs and their underlying regulatory landscapes.

### Visium HD
To analyze how melanoma cells interact with their microenvironment and how cellular neighborhoods shape tumor heterogeneity and progression, we investigated the heterogeneity and spatial organization of metastatic melanoma lesions under therapy in nine samples from four patients (R n=2, and NR n=2) using 10x Genomics VisiumHD profiling.

***

### `snRNA-seq/`
Contains the complete single-nucleus RNA-seq analytical workflow, including preprocessing and quality control, malignant cell identification, tumor microenvironment annotation, extraction and assignment of malignant transcriptional meta-programs, differential expression and functional enrichment analyses, inference of gene regulatory and cell–cell communication networks, and cross-cohort validation in independent single-cell and bulk melanoma datasets.

***

### `snATAC-seq/`
Contains the single-nucleus chromatin accessibility analysis workflow, including preprocessing and quality control, cell-type and tumor state annotation, projection of RNA-defined malignant meta-programs onto ATAC profiles, genome-wide peak calling and consensus region construction, and transcription factor motif enrichment analyses to characterize state-specific regulatory programs.

***

### `snRNA+snATAC/`
Contains the integrative multiome analysis framework combining single-nucleus RNA-seq and ATAC-seq to define malignant cell states, reconstruct gene regulatory networks using SCENIC+, characterize transcription factor activity and cis-regulatory landscapes, benchmark regulatory regions against ENCODE references, and perform in silico transcription factor perturbation simulations to evaluate meta-program stability and state transitions.

***

#### `Figure_script.R`
Contains the original code used to generate main and supplementary figures.

***

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





