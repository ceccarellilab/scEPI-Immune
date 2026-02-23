library(Seurat)
library(Signac)
library(dplyr)
library(ggplot2)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomicRanges)
library(parallel)
setwd("/home3/ciervo/scMULTIOME/Analisi/")
# dir.create("scATAC")
# dir.create("scATAC/RData")

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
rm(path_batch1, path_batch2, path_batch3, path_batch4, path_batch5, batch1, batch2, batch3, batch4, batch5)

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

combined.peaks <- reduce(unlist(as(GR.atac, "GRangesList")))

# Filter out bad peaks based on length
peakwidths <- width(combined.peaks)
combined.peaks <- combined.peaks[peakwidths  < 10000 & peakwidths > 20]
combined.peaks

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
    features = combined.peaks,
    cells = allBarcodes[[x]])
}, mc.cores = length(path.atac))
names(allPeaksCounts) <- names(path.atac)

annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevels(annotation) <- paste0('chr', seqlevels(annotation))

seu.final <- mclapply(names(allPeaksCounts), function(x){
  tmp <- CreateChromatinAssay(
    counts = allPeaksCounts[[x]], sep = c(":", "-"),
    fragments = allFragments[[x]], annotation = annotation, 
    meta.data = allBarcodes[[x]]
  )
  # allPeaksCounts[[x]]$Sample.ID <- x
  tmp
}, mc.cores = length(path.atac))
names(seu.final) = names(path.atac)

seu.final = mclapply(names(seu.final), function(x){
  CreateSeuratObject(
    counts = seu.final[[x]],
    assay = "ATAC"
  )
}, mc.cores = length(path.atac))
names(seu.final) = names(path.atac)

allMetadata <- mclapply(path.atac, function(x) {
  y = as.data.frame(readr::read_csv(paste0(x, "per_barcode_metrics.csv"))[, c("barcode", "atac_peak_region_fragments", "atac_fragments")])
  rownames(y) = y[, 1]
  return(y)
}, mc.cores = length(path.atac))
names(allMetadata) = names(path.atac)

seu.final = mclapply(1:length(seu.final), function(x) {
  seu.final[[x]] = AddMetaData(seu.final[[x]], allMetadata[[x]])
}, mc.cores = length(seu.final))
names(seu.final) = names(path.atac)

rm(allBarcodes, allFragments, allPeaksCounts, allMetadata, annotation, combined.peaks, GR.atac, peakwidths)

# save(seu.final, file = "scATAC/RData/scATAC.RData")

# merge all datasets, adding a cell ID to make sure cell names are unique
# load("scATAC/RData/scATAC.RData")
combined <- merge(
  x = seu.final[[1]],
  y = seu.final[2:length(seu.final)],
  add.cell.ids = names(seu.final)
)

combined$Barcode = rownames(combined@meta.data)
combined$OriginalBarcode = gsub(pattern = paste0(names(path.atac), collapse = "_|"), replacement = "", combined$Barcode)
combined$orig.ident = unlist(parallel::mclapply(1:length(combined$Barcode), FUN = 
                                                  function(x) paste0(strsplit(x = combined$Barcode, "_")[[x]][1:2], collapse = "_"), 
                                                mc.cores = 70)
)
View(combined@meta.data)
# save(combined, file = "Analisi_scRNA_scATAC/scATAC_combined.RData")

combined$Patient = combined$orig.ident
combined$Patient = gsub("Pat", "Patient ", combined$Patient)
combined$Patient = gsub("_W0|_W4|_W12", "", combined$Patient)
combined$Patient = gsub("6", "06", combined$Patient)
combined$Patient = gsub("3", "03", combined$Patient)
table(combined$Patient)

combined$Week = combined$orig.ident
combined$Week = gsub("Pat02_|Pat13_|Pat14_|Pat3_|Pat6_", "", combined$Week)
combined$Week = gsub("W", "Week ", combined$Week)
table(combined$Week)

combined$Batch = "Batch 3"
combined$Batch[combined$orig.ident %in% batch1] = "Batch 1"
combined$Batch[combined$orig.ident %in% batch2] = "Batch 2"
combined$Batch[combined$orig.ident %in% batch4] = "Batch 4"
combined$Batch[combined$orig.ident %in% batch5] = "Batch 5"
table(combined$Batch, exclude = NULL)

colnames(combined@meta.data)

save(combined, file = "scATAC/RData/scATAC_combined.RData")

# Filtering
combined <- NucleosomeSignal(combined, n = NULL)
# save(combined, file = "scATAC/RData/scATAC_combined.RData")

combined <- TSSEnrichment(combined)
combined$pct_reads_in_peaks <- combined$atac_peak_region_fragments / combined$atac_fragments * 100
combined$blacklist_ratio <- FractionCountsInRegion(
  object = combined, 
  assay = "ATAC",
  regions = blacklist_hg38_unified
)

# QC metrics 
DensityScatter(combined, x = 'nCount_ATAC', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)

VlnPlot(object = combined,
        features = c("nCount_ATAC", "nFeature_ATAC", "TSS.enrichment", "nucleosome_signal", 
                     "blacklist_ratio","pct_reads_in_peaks"),
        pt.size = 0.00005
)

VlnPlot(object = combined,
        features = c("nCount_ATAC", "nFeature_ATAC", "TSS.enrichment", "nucleosome_signal", 
                     "blacklist_ratio","pct_reads_in_peaks"),
        pt.size = 0, group.by = 'orig.ident'
)

combined$nucleosome_group = ifelse(combined$nucleosome_signal > 1.5, 'high', 'low')

FragmentHistogram(object = combined, group.by = 'nucleosome_signal')

dim(combined)
# 263,139  88,834

summary(combined$nCount_ATAC)
summary(combined$nFeature_ATAC)
dim(combined)
# 263,139  88,834

combined <- subset(
  x = combined,
  subset = nCount_ATAC > 2500 &
    nCount_ATAC < 20000 &
    pct_reads_in_peaks > 20 &
    blacklist_ratio < 0.05 &
    nucleosome_signal < 1.5 &
    TSS.enrichment > 2
)

dim(combined)
# 263,139  55,850

VlnPlot(object = combined,
        features = c("nCount_ATAC", "nFeature_ATAC", "TSS.enrichment", "nucleosome_signal"),
        ncol = 4,
        pt.size = 0,
        group.by = "Cell_annotation"
)
dim(combined)
any(rowSums(GetAssayData(object = combined, layer = "counts")) == 0) # FALSE
save(combined, file = "scATAC/RData/scATAC_combined.RData")

# scATAC workflow ----
load(file = "scATAC/RData/scATAC_combined.RData")

DefaultAssay(combined) <- 'ATAC'
combined <- RunTFIDF(combined)
combined <- FindTopFeatures(combined, min.cutoff = 'q5')
combined <- RunSVD(object = combined)

DepthCor(combined, n = 50)

combined <- FindNeighbors(object = combined, reduction = 'lsi', dims = 2:50)
combined <- FindClusters(object = combined, algorithm = 3, resolution = 0.8)
combined <- RunTSNE(object = combined, reduction = 'lsi', dims = 2:50)

DimPlot(object = combined, label = TRUE, reduction = "tsne")

DimPlot(object = combined, label = FALSE, group.by = "Patient") /
  DimPlot(object = combined, label = FALSE, reduction = "tsne", group.by = "Week", cols = week_cols)
save(combined, file = "scATAC/RData/scATAC_combined.RData")
