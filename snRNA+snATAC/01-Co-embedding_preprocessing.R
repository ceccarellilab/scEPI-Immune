###################
#### Libraries ####
###################
library(Seurat)
library(Signac)
library(dplyr)
library(ggplot2)
library(SingleCellExperiment)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomicRanges)
library(parallel)

setwd("/home3/ciervo/scMULTIOME/Analisi/")
# dir.create("scRNA_scATAC/RData/")

source('colori_finali.R')

# Load objects
# Batch 1
path_batch1 = "/home3/adefalco/scRNA+ATACT_NIBIT/Batch_1/"
batch1 = list.dirs(path_batch1, recursive = F, full.names = F)[1:3]

# Batch 2
path_batch2 = "/home3/adefalco/scRNA+ATACT_NIBIT/Batch_2/"
batch2 = list.dirs(path_batch2, recursive = F, full.names = F)

# Batch 3
path_batch3 = "/home3/adefalco/scRNA+ATACT_NIBIT/Batch_3/"
batch3 = list.dirs(path_batch3, recursive = F, full.names = F)

# Batch 4
path_batch4 = "/home3/ciervo/scMULTIOME/Batch_4/"
batch4 = list.dirs(path_batch4, recursive = F, full.names = F)

# Batch 5
path_batch5 = "/home3/ciervo/scMULTIOME/Batch_5/"
batch5 = list.dirs(path_batch5, recursive = F, full.names = F)[1]

### 1. Read in peak sets
path.atac = c(paste0(path_batch1, batch1, "/outs/"),
              paste0(path_batch2, batch2, "/outs/"),
              paste0(path_batch3, batch3, "/outs/"),
              paste0(path_batch4, batch4, "/outs/"),
              paste0(path_batch5, batch5, "/outs/")
)
names(path.atac) <- c(batch1, batch2, batch3, batch4, batch5)
# rm(path_batch1, path_batch2, path_batch3, path_batch4, path_batch5, batch1, batch2, batch3, batch4, batch5)

path.atac
GR.atac <- mclapply(path.atac, function(x){
  peaks <- read.table(
    file = paste0(x, "/atac_peaks.bed"),
    col.names = c("chr", "start", "end"))
  makeGRangesFromDataFrame(peaks)
},
mc.cores = length(path.atac)
)
names(GR.atac) = names(path.atac)

seu.multiome.peaks <- reduce(unlist(as(GR.atac, "GRangesList")))

# Filter out bad peaks based on length
peakwidths <- width(seu.multiome.peaks)
seu.multiome.peaks <- seu.multiome.peaks[peakwidths  < 10000 & peakwidths > 20]
seu.multiome.peaks

# Create barcodes
allBarcodes <- mclapply(path.atac, function(x) {
  read.table(paste0(x, "/filtered_feature_bc_matrix/barcodes.tsv.gz"), header = FALSE)$V1
}, mc.cores = length(path.atac))

# Create fragment objects
allFragments <- mclapply(names(path.atac), function(x){
  CreateFragmentObject(
    path = paste0(path.atac[x], "atac_fragments.tsv.gz"),
    cells = allBarcodes[[x]]
  )
}, mc.cores = length(path.atac))
names(allFragments) <- names(path.atac)

# Quantify peaks in each dataset
allPeaksCounts <- mclapply(names(path.atac), function(x){
  FeatureMatrix(
    fragments = allFragments[[x]],
    features = seu.multiome.peaks,
    cells = allBarcodes[[x]])
}, mc.cores = length(path.atac))
names(allPeaksCounts) <- names(path.atac)

allMetadata <- mclapply(path.atac, function(x) {
  y = as.data.frame(readr::read_csv(paste0(x, "per_barcode_metrics.csv"))[, c("barcode", "atac_peak_region_fragments", "atac_fragments")])
  rownames(y) = y[, 1]
  return(y)
}, mc.cores = length(path.atac))
names(allMetadata) = names(path.atac)

# Create the objects
# RNA 
allCountsGE <- mclapply(path.atac, function (x){
  counts <- Read10X(paste0(x, "filtered_feature_bc_matrix"))
  # create a Seurat seu.multiome containing the RNA data
  countsGE <- CreateSeuratObject(
    counts = counts$`Gene Expression`,
    assay = "RNA",
    min.features = 500, min.cells = 3
  )
  countsGE
}, mc.cores = length(path.atac))
names(allCountsGE) <- names(path.atac)

annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevels(annotation) <- paste0('chr', seqlevels(annotation))

seu.finale <- mclapply(1:length(allCountsGE), function(x){
  allCountsGE[[x]][["ATAC"]] <- CreateChromatinAssay(
    counts = allPeaksCounts[[x]], sep = c(":", "-"),
    fragments = allFragments[[x]], 
    annotation = annotation,
    meta.data = allBarcodes[[x]]
  )
  # x$Sample.ID <- x
  allCountsGE[[x]]
}, mc.cores = length(path.atac))
names(seu.finale) = names(path.atac)

seu.finale = mclapply(1:length(seu.finale), function(x) {
  seu.finale[[x]] = AddMetaData(seu.finale[[x]], allMetadata[[x]])
}, mc.cores = length(seu.finale))
names(seu.finale) = names(path.atac)

rm(allBarcodes, allFragments, allPeaksCounts, allCountsGE, allMetadata, annotation, seu.multiome.peaks, GR.atac, peakwidths)
save(seu.finale, file = "scRNA_scATAC/RData/multiome.RData")

#### Merge dei SeuObj
load("scRNA_scATAC/RData/multiome.RData")

head(colnames(seu.finale$Pat02_W0))
names(seu.finale) = gsub("Pat3", "Pat03", names(seu.finale))
names(seu.finale) = gsub("Pat6", "Pat06", names(seu.finale))
names(seu.finale) = gsub("W0_1", "W0", names(seu.finale))

seu.multiome = merge(
  x = seu.finale[[1]],
  y = seu.finale[2:length(seu.finale)],
  add.cell.ids = names(seu.finale)
)

seu.multiome@assays$RNA
seu.multiome@assays$ATAC

DefaultAssay(seu.multiome)
dim(seu.multiome) # 36,601 genes x 88,834 cells
DefaultAssay(seu.multiome) = 'ATAC'
dim(seu.multiome) # 263,139 peaks x 88,834 cells
DefaultAssay(seu.multiome) = 'RNA'
rm(seu.finale)

save(seu.multiome, file = "scRNA_scATAC/RData/seu_multiome.RData")

# Filtering merged object based on filters applied separatelly ----
seu.multiome$percent.mt = PercentageFeatureSet(seu.multiome, pattern = "^MT-")

DefaultAssay(seu.multiome) = "ATAC"

seu.multiome <- NucleosomeSignal(seu.multiome)
seu.multiome <- TSSEnrichment(seu.multiome)

seu.multiome$pct_reads_in_peaks <- seu.multiome$atac_peak_region_fragments / seu.multiome$atac_fragments * 100
seu.multiome$blacklist_ratio <- FractionCountsInRegion(
  object = seu.multiome, 
  assay = "ATAC",
  regions = blacklist_hg38_unified
)

seu.multiome$orig.ident = "N/A"
seu.multiome$orig.ident[grep("Pat02", rownames(seu.multiome@meta.data))] = "Pat02"
seu.multiome$orig.ident[grep("Pat03", rownames(seu.multiome@meta.data))] = "Pat03"
seu.multiome$orig.ident[grep("Pat06", rownames(seu.multiome@meta.data))] = "Pat06"
seu.multiome$orig.ident[grep("Pat13", rownames(seu.multiome@meta.data))] = "Pat13"
seu.multiome$orig.ident[grep("Pat14", rownames(seu.multiome@meta.data))] = "Pat14"
table(seu.multiome$orig.ident)

seu.multiome$Week = "N/A"
seu.multiome$Week[grep("W0", rownames(seu.multiome@meta.data))] = "W0"
seu.multiome$Week[grep("W4", rownames(seu.multiome@meta.data))] = "W4"
seu.multiome$Week[grep("W12", rownames(seu.multiome@meta.data))] = "W12"
table(seu.multiome$Week)

seu.multiome$Patient = seu.multiome$orig.ident
seu.multiome$Patient = gsub("Pat", "Patient ", seu.multiome$Patient)
table(seu.multiome$Patient)

seu.multiome$orig.ident = paste0(seu.multiome$orig.ident, "_", seu.multiome$Week)
table(seu.multiome$orig.ident)

seu.multiome$Batch = "Batch 3"
seu.multiome$Batch[seu.multiome$orig.ident %in% c("Pat02_W0", "Pat02_W4", "Pat02_W12")] = "Batch 1"
seu.multiome$Batch[seu.multiome$orig.ident %in% c("Pat13_W4")] = "Batch 2"
seu.multiome$Batch[seu.multiome$orig.ident %in% c("Pat14_W12", "Pat03_W4", "Pat06_W0", "Pat06_W4")] = "Batch 4"
seu.multiome$Batch[seu.multiome$orig.ident %in% c("Pat14_W4")] = "Batch 5"
table(seu.multiome$Batch)
save(seu.multiome, file = "scRNA_scATAC/RData/seu_multiome.RData")

# -- #
load("scRNA_scATAC/RData/seu_multiome.RData")
dim(seu.multiome) # 88,834
table(seu.multiome$orig.ident)

seu.multiome$Passed.QC = 'Not passed QC'
seu.multiome$Passed.QC[seu.multiome$nCount_RNA < 50000 & seu.multiome$percent.mt < 5 & seu.multiome$nFeature_RNA < 7500] = 'Passed RNA QC'
seu.multiome$Passed.QC[seu.multiome$Passed.QC %in% 'Passed RNA QC' & seu.multiome$nCount_ATAC > 2500 & seu.multiome$nCount_ATAC < 20000 & seu.multiome$pct_reads_in_peaks > 20 & seu.multiome$blacklist_ratio < 0.05 & seu.multiome$nucleosome_signal < 1.5 & seu.multiome$TSS.enrichment > 2] = 'Passed RNA and ATAC QC'
seu.multiome$Passed.QC[seu.multiome$Passed.QC %in% 'Not passed QC' & seu.multiome$nCount_ATAC > 2500 & seu.multiome$nCount_ATAC < 20000 & seu.multiome$pct_reads_in_peaks > 20 & seu.multiome$blacklist_ratio < 0.05 & seu.multiome$nucleosome_signal < 1.5 & seu.multiome$TSS.enrichment > 2] = 'Passed ATAC QC'
table(seu.multiome$Passed.QC)

VlnPlot(seu.multiome, features = c('nCount_RNA', 'nFeature_RNA', 'percent.mt', 'nCount_ATAC', 'nFeature_ATAC', 'pct_reads_in_peaks', 'blacklist_ratio', 'nucleosome_signal', 'TSS.enrichment'), group.by = 'Passed.QC', pt.size = 0)
p1 = VlnPlot(seu.multiome, 'nCount_RNA', pt.size = 0, group.by = 'Passed.QC') + NoLegend() + geom_hline(yintercept = 50000, linetype = 'dashed', color = 'grey80')
p2 = VlnPlot(seu.multiome, 'nFeature_RNA', pt.size = 0, group.by = 'Passed.QC') + NoLegend() + geom_hline(yintercept = 7500, linetype = 'dashed', color = 'grey80')
p3 = VlnPlot(seu.multiome, 'percent.mt', pt.size = 0, group.by = 'Passed.QC') + NoLegend() + geom_hline(yintercept = 5, linetype = 'dashed', color = 'grey80')
cowplot::plot_grid(p1, p2, p3, ncol = 1)

p1 = VlnPlot(seu.multiome, 'nCount_ATAC', pt.size = 0, group.by = 'Passed.QC') + NoLegend() + geom_hline(yintercept = c(2500, 20000), linetype = 'dashed', color = 'grey80')
p2 = VlnPlot(seu.multiome, 'pct_reads_in_peaks', pt.size = 0, group.by = 'Passed.QC') + NoLegend() + geom_hline(yintercept = 20, linetype = 'dashed', color = 'grey80')
p3 = VlnPlot(seu.multiome, 'blacklist_ratio', pt.size = 0, group.by = 'Passed.QC') + NoLegend() + geom_hline(yintercept = 0.05, linetype = 'dashed', color = 'grey80')
p4 = VlnPlot(seu.multiome, 'nucleosome_signal', pt.size = 0, group.by = 'Passed.QC') + NoLegend() + geom_hline(yintercept = 1.5, linetype = 'dashed', color = 'grey80')
p5 = VlnPlot(seu.multiome, 'TSS.enrichment', pt.size = 0, group.by = 'Passed.QC') + NoLegend() + geom_hline(yintercept = 2, linetype = 'dashed', color = 'grey80')

save(seu.multiome, file = 'scRNA_scATAC/RData/seu_multiome.RData')

# Add metadata ----
load("scRNA_scATAC/RData/seu_multiome.RData")
load("scRNA/RData/merged_object.RData", verbose = TRUE)
load("scATAC/RData/scATAC_macs2.RData", verbose = TRUE)
colnames(seu@meta.data)
colnames(combined.new@meta.data)
df.metadata = rbind(seu@meta.data[, c('Barcode', 'Cell_annotation', 'Metaprogram_assignment')], combined.new@meta.data[, c('Barcode', 'Cell_annotation', 'Metaprogram_assignment')])
df.metadata = df.metadata[!duplicated(df.metadata$Barcode), ]
save(df.metadata, file = 'scRNA_scATAC/RData/df_metadata_all.RData')

rownames(df.metadata) = gsub('Pat03_W0_1', 'Pat03_W0', rownames(df.metadata))
sum(rownames(df.metadata) %in% rownames(seu.multiome@meta.data)) # 76491

seu.multiome = AddMetaData(seu.multiome, metadata = df.metadata)
View(seu.multiome@meta.data)
rm(meta.data, p1, p2, p3, p4, p5)

save(seu.multiome, file = 'scRNA_scATAC/RData/seu_multiome.RData')

# Subset multiome dataset ----
load("scRNA_scATAC/RData/seu_multiome.RData")
table(seu.multiome$orig.ident, seu.multiome$Metaprogram_assignment)
# seu.multiome = subset(seu.multiome, Passed.QC %in% 'Not passed QC', invert = TRUE)
# dim(seu.multiome) # 76,491 cells

table(seu.multiome$Passed.QC)
seu.multiome = subset(seu.multiome, Passed.QC %in% 'Passed RNA and ATAC QC')
dim(seu.multiome) # 38,029 cells

# RNA processing ----
DefaultAssay(seu.multiome) = "RNA"
# seu.multiome = JoinLayers(seu.multiome)
seu.multiome = seu.multiome %>% 
  NormalizeData() %>% 
  FindVariableFeatures(nfeatures = 3000) %>% 
  ScaleData(features = rownames(seu.multiome)) %>% 
  RunPCA()

nPCA(seu.multiome)

seu.multiome <- seu.multiome %>% 
  FindNeighbors(reduction = 'pca', dims = 1:14) %>% 
  FindClusters(resolution = 0.1) %>% 
  RunTSNE(dims = 1:14, reduction.name = 'tsne_RNA', reduction.key = 'tsneRNA_')

# p1 = DimPlot(seu.multiome, reduction = 'tsne_RNA')
DimPlot(seu.multiome, reduction = 'tsne_RNA', group.by = 'RNA_snn_res.0.1') + 
  DimPlot(seu.multiome, group.by = 'orig.ident', reduction = 'tsne_RNA') + 
  DimPlot(seu.multiome, group.by = 'Passed.QC', reduction = 'tsne_RNA') +
  DimPlot(seu.multiome, group.by = 'Cell_annotation', reduction = 'tsne_RNA', cols = ccolors)
# DimPlot(seu.multiome, group.by = 'Cell_annotation', reduction = 'tsne_ATAC', cols = ccolors)

FeaturePlot(seu.multiome, features = c('nCount_RNA', 'nFeature_RNA', 'percent.mt'), ncol = 3)

DimPlot(seu.multiome, group.by = 'Cell_annotation', reduction = 'tsne_RNA', cols = ccolors) + 
  DimPlot(seu.multiome, group.by = 'RNA_snn_res.0.1', reduction = 'tsne_RNA')

table(seu.multiome$Cell_annotation, seu.multiome$RNA_snn_res.0.1)

# ATAC processing ----
DefaultAssay(seu.multiome) = "ATAC"

seu.multiome <- seu.multiome %>% 
  RunTFIDF() %>% 
  FindTopFeatures(min.cutoff = "q75") %>% 
  RunSVD() %>% 
  FindNeighbors(reduction = 'lsi', dims = 2:50) %>% 
  FindClusters(resolution = 0.1, algorithm = 3) %>% 
  RunTSNE(dims = 2:50, reduction.name = 'tsne_ATAC', reduction = 'lsi', reduction.key = 'tsneATAC_')

DimPlot(seu.multiome, reduction = 'tsne_ATAC') + DimPlot(seu.multiome, reduction = 'tsne_ATAC', group.by = 'orig.ident')
DimPlot(seu.multiome, reduction = 'tsne_ATAC', group.by = 'Cell_annotation', cols = ccolors)

DimPlot(seu.multiome, reduction = 'tsne_ATAC', group.by = 'ATAC_snn_res.0.1') +
  DimPlot(seu.multiome, group.by = 'orig.ident', reduction = 'tsne_ATAC') + 
  # DimPlot(seu.multiome, group.by = 'Passed.QC', reduction = 'tsne_ATAC') +
  DimPlot(seu.multiome, group.by = 'Cell_annotation', reduction = 'tsne_ATAC', cols = ccolors)

# build a joint neighbor graph using both assays
seu.multiome <- FindMultiModalNeighbors(
  object = seu.multiome,
  reduction.list = list("pca", "lsi"), 
  dims.list = list(1:14, 2:50),
  modality.weight.name = c("RNA.weight", "ATAC.weight"),
  verbose = TRUE
)

# build a joint tSNE visualization
# DefaultAssay(seu.multiome) = 'RNA'
seu.multiome <- RunTSNE(
  object = seu.multiome,
  nn.name = "weighted.nn",
  assay = "RNA",
  verbose = TRUE, 
  reduction.name = 'tsne_multimodal',
)

# Clustering multimodal 
seu.multiome <- FindClusters(object = seu.multiome, resolution = 0.1, graph.name="wsnn", algorithm = 3)

# seu.multiome@reductions$tsne_multimodal

DimPlot(seu.multiome, reduction = 'tsne_multimodal', group.by = "Patient") +
  DimPlot(seu.multiome, reduction = 'tsne_multimodal', group.by = "Week", cols = week_cols) +
  DimPlot(seu.multiome, reduction = 'tsne_multimodal', group.by = "Cell_annotation", cols = ccolors)

seu.multiome$Week[seu.multiome$Week == 'W0'] = 'Week 0'
seu.multiome$Week[seu.multiome$Week == 'W4'] = 'Week 4'
seu.multiome$Week[seu.multiome$Week == 'W12'] = 'Week 12'

seu.multiome$Week = factor(seu.multiome$Week, levels = c('Week 0', 'Week 4', 'Week 12'))
p1 = DimPlot(seu.multiome, reduction = 'tsne_multimodal', group.by = 'Patient') + ggtitle('Patient')
p2 = DimPlot(seu.multiome, reduction = 'tsne_multimodal', group.by = 'Week', cols = week_cols) + ggtitle('Week')
p3 = DimPlot(seu.multiome, reduction = 'tsne_multimodal', group.by = 'Cell_annotation', cols = ccolors) + ggtitle('Cell type')

cowplot::plot_grid(p1, p2, p3, ncol = 3, align = 'hv')

table(seu.multiome$Metaprogram_assignment, seu.multiome$Cell_annotation, exclude = NULL)

save(seu.multiome, file = 'scRNA_scATAC/RData/seu_multiome_filtered.RData')
rm(p1, p2, p3, df.metadata)

# Malignant subset ----
table(seu.multiome$Cell_annotation, seu.multiome$Metaprogram_assignment, exclude = NULL)
malignant.multiome = subset(seu.multiome, Cell_annotation %in% 'Malignant')
dim(malignant.multiome) # 64532

# ATAC
DefaultAssay(malignant.multiome) = 'ATAC'
any(rowSums(malignant.multiome) == 0) # Check if some peaks has 0 values: FALSE

malignant.multiome = malignant.multiome %>% 
  RunTFIDF() %>% 
  FindTopFeatures(min.cutoff = "q75") %>% 
  RunSVD() %>% 
  FindNeighbors(reduction = 'lsi', dims = 2:50) %>% 
  FindClusters(resolution = 0.1, algorithm = 3) %>% 
  RunTSNE(dims = 2:50, reduction.name = 'tsne_ATAC', reduction = 'lsi', reduction.key = 'tsneATAC_')

# RNA
DefaultAssay(malignant.multiome) = 'RNA'
malignant.multiome = malignant.multiome %>% 
  NormalizeData() %>% 
  FindVariableFeatures(nfeatures = 3000) %>% 
  ScaleData(features = rownames(malignant.multiome)) %>% 
  RunPCA()

ndims = nPCA(malignant.multiome) # 13
malignant.multiome = malignant.multiome %>% 
  FindNeighbors(reduction = 'pca', dims = 1:ndims) %>% 
  FindClusters(resolution = 0.1) %>% 
  RunTSNE(dims = 1:ndims, reduction.name = 'tsne_RNA', reduction.key = 'tsneRNA_')

# Multiome
# build a joint neighbor graph using both assays
malignant.multiome <- FindMultiModalNeighbors(
  object = malignant.multiome,
  reduction.list = list("pca", "lsi"), 
  dims.list = list(1:ndims, 2:50),
  modality.weight.name = c("RNA.weight", "ATAC.weight"),
  verbose = TRUE
)

# build a joint tSNE visualization
# DefaultAssay(malignant.multiome) = 'RNA'
malignant.multiome <- RunTSNE(
  object = seu.multiome,
  nn.name = "weighted.nn",
  assay = "RNA",
  verbose = TRUE, 
  reduction.name = 'tsne_multimodal',
)

malignant.multiome$Week = factor(malignant.multiome$Week, levels = names(week_cols))
p1 = DimPlot(malignant.multiome, reduction = 'tsne_multimodal', group.by = 'Patient') + ggtitle('Patient')
p2 = DimPlot(malignant.multiome, reduction = 'tsne_multimodal', group.by = 'Week', cols = week_cols) + ggtitle('Week')
p3 = DimPlot(malignant.multiome, reduction = 'tsne_multimodal', group.by = 'Metaprogram_assignment', cols = colori_mp) + ggtitle('Metaprograms')

cowplot::plot_grid(p1, p2, p3, ncol = 3, align = 'hv')
save(malignant.multiome, file = 'scRNA_scATAC/RData/malignant_subset.RData')
rm(malignant.multiome, p1, p2, p3)

# TME subset ----
table(seu.multiome$Cell_annotation, seu.multiome$Metaprogram_assignment, exclude = NULL)
tme.multiome = subset(seu.multiome, Cell_annotation %in% 'Malignant', invert = TRUE)
tme.multiome = subset(tme.multiome, Passed.QC %in% 'Passed RNA and ATAC QC')
dim(tme.multiome) # 3207
DefaultAssay(tme.multiome) = 'ATAC'

tme.multiome = tme.multiome %>% 
  RunTFIDF() %>% 
  FindTopFeatures(min.cutoff = "q90") %>% 
  RunSVD()

tme.multiome = tme.multiome %>% 
  FindNeighbors(reduction = 'lsi', dims = 2:15) %>% 
  FindClusters(resolution = 0.1, algorithm = 3) %>% 
  RunTSNE(dims = 2:15, reduction.name = 'tsne_ATAC', reduction = 'lsi', reduction.key = 'tsneATAC_')

DefaultAssay(tme.multiome) = 'RNA'
tme.multiome = tme.multiome %>% 
  NormalizeData() %>% 
  FindVariableFeatures(nfeatures = 2000) %>% 
  ScaleData(features = rownames(tme.multiome), vars.to.regress = 'percent.mt') %>% 
  RunPCA()

ndims = nPCA(tme.multiome) # 17
tme.multiome = tme.multiome %>% 
  FindNeighbors(reduction = 'pca', dims = 1:ndims) %>% 
  FindClusters(resolution = 0.1) %>% 
  RunTSNE(dims = 1:ndims, reduction.name = 'tsne_RNA', reduction.key = 'tsneRNA_', check_duplicates = FALSE)

# build a joint neighbor graph using both assays
tme.multiome <- FindMultiModalNeighbors(
  object = tme.multiome,
  reduction.list = list("pca", "lsi"), 
  dims.list = list(1:17, 2:15),
  modality.weight.name = c("RNA.weight", "ATAC.weight"),
  verbose = TRUE
)

# build a joint tSNE visualization
# DefaultAssay(malignant.multiome) = 'RNA'
tme.multiome <- RunTSNE(
  object = tme.multiome,
  nn.name = "weighted.nn",
  assay = "ATAC",
  verbose = TRUE, 
  reduction.name = 'tsne_multimodal',
  check_duplicates = FALSE
)

tme.multiome$Week = factor(tme.multiome$Week, levels = names(week_cols))
tme.multiome$Cell_annotation = factor(tme.multiome$Cell_annotation, levels = c(names(colors_tme), 'Cycling cells'))
p1 = DimPlot(tme.multiome, reduction = 'tsne_multimodal', group.by = 'Patient') + ggtitle('Patient')
p2 = DimPlot(tme.multiome, reduction = 'tsne_multimodal', group.by = 'Week', cols = week_cols) + ggtitle('Week')
p3 = DimPlot(tme.multiome, reduction = 'tsne_multimodal', group.by = 'Cell_annotation', cols = colors_tme, na.value = 'grey95') + ggtitle('Cell type')

cowplot::plot_grid(p1, p2, p3, ncol = 3, align = 'hv')
save(tme.multiome, file = 'scRNA_scATAC/RData/tme_subset.RData')

# tme.filtered
DefaultAssay(tme.multiome.filtered) = 'ATAC'
tme.multiome.filtered = tme.multiome.filtered %>% 
  RunTFIDF() %>% 
  FindTopFeatures(min.cutoff = "q90") %>% 
  RunSVD()
length(TopFeatures(tme.multiome.filtered))

tme.multiome.filtered = tme.multiome.filtered %>% 
  FindNeighbors(reduction = 'lsi', dims = 2:10) %>% 
  FindClusters(resolution = 0.1, algorithm = 3) %>% 
  RunTSNE(dims = 2:10, reduction.name = 'tsne_ATAC', reduction = 'lsi', reduction.key = 'tsneATAC_')

DefaultAssay(tme.multiome.filtered) = 'RNA'
tme.multiome.filtered = tme.multiome.filtered %>% 
  NormalizeData() %>% 
  FindVariableFeatures(nfeatures = 2000) %>% 
  ScaleData(features = rownames(tme.multiome.filtered)) %>% 
  RunPCA()

ndims = nPCA(tme.multiome.filtered) # 17
tme.multiome.filtered = tme.multiome.filtered %>% 
  FindNeighbors(reduction = 'pca', dims = 1:ndims) %>% 
  FindClusters(resolution = 0.1) %>% 
  RunTSNE(dims = 1:ndims, reduction.name = 'tsne_RNA', reduction.key = 'tsneRNA_', check_duplicates = FALSE)

# build a joint neighbor graph using both assays
tme.multiome.filtered <- FindMultiModalNeighbors(
  object = tme.multiome.filtered,
  reduction.list = list("pca", "lsi"), 
  dims.list = list(1:13, 2:10),
  modality.weight.name = c("RNA.weight", "ATAC.weight"),
  verbose = TRUE
)

# build a joint tSNE visualization
# DefaultAssay(tme.multiome.filtered) = 'RNA'
tme.multiome.filtered <- RunTSNE(
  object = tme.multiome.filtered,
  nn.name = "weighted.nn",
  assay = "RNA",
  verbose = TRUE, 
  reduction.name = 'tsne_multimodal',
  check_duplicates = FALSE
)

DimPlot(tme.multiome.filtered, reduction = 'tsne_multimodal', group.by = 'Patient') + ggtitle('Patient') +
  DimPlot(tme.multiome.filtered, reduction = 'tsne_multimodal', group.by = 'Week', cols = week_cols) + ggtitle('Week') +
  DimPlot(tme.multiome.filtered, reduction = 'tsne_multimodal', group.by = 'Cell_annotation', cols = colors_tme, na.value = 'grey95') + ggtitle('Cell type')

tme.multiome$Week = factor(tme.multiome$Week, levels = names(week_cols))
tme.multiome$Cell_annotation = factor(tme.multiome$Cell_annotation, levels = c(names(colors_tme), 'Cycling cells'))
p1 = DimPlot(tme.multiome, reduction = 'tsne_multimodal', group.by = 'Patient') + ggtitle('Patient')
p2 = DimPlot(tme.multiome, reduction = 'tsne_multimodal', group.by = 'Week', cols = week_cols) + ggtitle('Week')
p3 = DimPlot(tme.multiome, reduction = 'tsne_multimodal', group.by = 'Cell_annotation', cols = colors_tme, na.value = 'grey95') + ggtitle('Cell type')

cowplot::plot_grid(p1, p2, p3, ncol = 3, align = 'hv')
save(tme.multiome.filtered, file = 'scRNA_scATAC/RData/tme_subset.RData')
