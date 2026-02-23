library(Seurat)
library(Signac)
library(dplyr)
library(ggplot2)
library(SingleCellExperiment)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomicRanges)
library(parallel)
library(future)
setwd("/home3/ciervo/scMULTIOME/Analisi/")

# Peak calling after cell type annotation ----
load('/home3/ciervo/scMULTIOME/Analisi/scATAC/RData/scATAC_combined.RData')

combined$Fine_annotation = combined$Metaprogram_assignment
combined$Fine_annotation[is.na(combined$Fine_annotation)] = combined$Cell_annotation[is.na(combined$Fine_annotation)]
table(combined$Fine_annotation, exclude = NULL)

peaks = CallPeaks(combined, group.by = "Fine_annotation", macs2.path = "/storage/qnap_home/ciervo/miniconda3/envs/macs2/bin/macs2")
head(peaks)
peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_hg38_unified, invert = TRUE)
head(peaks)
save(peaks, file = '/home3/ciervo/scMULTIOME/Analisi/scATAC/RData/peak_calling_celltype.RData')

# load('/home3/ciervo/scMULTIOME/Analisi/scATAC/RData/peak_calling_celltype.RData')
# load('/home/noviello/TE_utilities/Genome_class_hg38_extended.RDa')

plan('multisession', workers = 50)
options(future.globals.maxSize = 50 * 1024 ^ 3) # for 50 Gb RAM
combined.new = FeatureMatrix(fragments = Fragments(combined), features = peaks)
plan('sequential')

combined.new <- CreateChromatinAssay(combined.new, fragments = Fragments(combined), annotation = Annotation(combined))
combined.new <- CreateSeuratObject(combined.new, assay = "ATAC", meta.data = combined@meta.data)

DefaultAssay(combined.new) <- 'ATAC'
combined.new <- RunTFIDF(combined.new)
combined.new <- FindTopFeatures(combined.new, min.cutoff = 'q75')
combined.new <- RunSVD(object = combined.new)

DepthCor(combined.new, n = 50)

combined.new <- FindNeighbors(object = combined.new, reduction = 'lsi', dims = 2:50)
combined.new <- FindClusters(object = combined.new, algorithm = 3, resolution = 0.5)
combined.new <- RunTSNE(object = combined.new, reduction = 'lsi', dims = 2:50)

DimPlot(object = combined.new, label = TRUE)

DimPlot(object = combined.new, label = FALSE, group.by = "Patient") /
  DimPlot(object = combined.new, label = FALSE, reduction = "tsne", group.by = "Week", cols = week_cols)

DimPlot(object = combined.new, label = FALSE, reduction = "tsne", group.by = 'Cell_annotation', cols = ccolors)

save(combined.new, file = "/home3/ciervo/scMULTIOME/Analisi/scATAC/RData/scATAC_macs2.RData")

# TME subset ----
load(file = "scATAC/RData/scATAC_macs2.RData")
table(combined.new$Cell_annotation, exclude = NULL)
tme = subset(combined.new, Cell_annotation %in% c("Cycling cells", "Malignant"), invert = TRUE)

DefaultAssay(tme) <- 'ATAC'
tme <- RunTFIDF(tme)
tme <- FindTopFeatures(tme, min.cutoff = 'q75')
tme <- RunSVD(object = tme)

# DepthCor(tme, n = 50)

tme <- FindNeighbors(object = tme, reduction = 'lsi', dims = 2:25)
tme <- FindClusters(object = tme, algorithm = 3, resolution = 0.5)
tme <- RunTSNE(object = tme, reduction = 'lsi', dims = 2:25)

DimPlot(tme)

p1 = DimPlot(object = tme, label = FALSE, reduction = "tsne", group.by = 'Cell_annotation', cols = colors_tme) + NoLegend()
p2 = DimPlot(object = tme, label = FALSE, group.by = "Patient") / DimPlot(object = tme, label = FALSE, reduction = "tsne", group.by = "Week", cols = week_cols)

cowplot::plot_grid(p1, p2, rel_widths = c(.7, .6))

# Accessibility dotplot
DefaultAssay(tme)
future::plan('multisession', workers = 20)
options(future.globals.maxSize = 50 * 1024 ^ 3) # for 50 Gb RAM
gene.activities <- GeneActivity(tme)
future::plan('sequential')

tme[['Activity']] <- CreateAssayObject(counts = gene.activities)
DefaultAssay(tme) = 'Activity'
tme <- NormalizeData(object = tme, 
                     assay = 'Activity', 
                     normalization.method = 'LogNormalize',
                     scale.factor = median(tme$nCount_Activity)
)
rm(gene.activities)

DefaultAssay(tme) = 'Activity'
save(tme, file = 'scATAC/RData/scATAC_tme.RData')

genes_to_plot = c("MS4A1", 'BANK1', "MZB1", 'JCHAIN', # Plasma cells; # B cells
                  'CD4', 'FOXP1', 'IL7R', # CD4 & CD4 Tcm
                  'LAG3', 'HAVCR2', 'PDCD1', # Exhausted 
                  "CD8A", 'CD8B', "GZMK", 'IL9R', 'RORC', # CD8 & CD8 Tcm
                  "FOXP3", 'IL2RA', # Tregs
                  'IL10', 'CD86', # M1
                  'MARCO', 'CTSD',  # M2
                  "PDGFRA", 'PDPN', "ACTA2", "RGS5", # iCAF; myoCAF
                  "PECAM1", 'VWF', # Endothelial cells
                  "KRT14", 'KRT1' # Keratinocytes
)

tme$Cell_annotation = factor(tme$Cell_annotation, levels = names(colors_tme))
DotPlot(tme, features = genes_to_plot, group.by = 'Cell_annotation') + 
  RotatedAxis() + coord_flip() +
  viridis::scale_color_viridis()

# Metaprogram subset ----
load("/home3/ciervo/scMULTIOME/Analisi/scATAC/RData/scATAC_macs2.RData")
malignant = subset(combined.new, Cell_annotation %in% "Malignant")

DefaultAssay(malignant) <- 'ATAC'
malignant <- RunTFIDF(malignant)
malignant <- FindTopFeatures(malignant, min.cutoff = 'q75')
malignant <- RunSVD(object = malignant)

# DepthCor(malignant, n = 50)

malignant <- FindNeighbors(object = malignant, reduction = 'lsi', dims = 2:50)
malignant <- FindClusters(object = malignant, algorithm = 3, resolution = 0.5)
malignant <- RunTSNE(object = malignant, reduction = 'lsi', dims = 2:50)

DimPlot(malignant)

p1 = DimPlot(object = malignant, label = FALSE, reduction = "tsne", group.by = 'Metaprogram_assignment', cols = colori_mp) + NoLegend()
p2 = DimPlot(object = malignant, label = FALSE, group.by = "Patient") / DimPlot(object = malignant, label = FALSE, reduction = "tsne", group.by = "Week", cols = week_cols)

cowplot::plot_grid(p1, p2, rel_widths = c(.7, .6))
rm(p1, p2)

save(malignant, file = 'scATAC/RData/malignant_subset.RData')
