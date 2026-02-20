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

# Annotate MPs by transfering label from RNA to ATAC
load(file = "scATAC/RData/scATAC_combined.RData")
load(file = "scRNA/RData/malignant_subset.RData")

dim(seu) # 53,702 cells
table(seu$assegnazione_filtrata)

DefaultAssay(seu)
table(combined$Cell_annotation, exclude = NULL)
malignant = subset(combined, Cell_annotation %in% "Malignant")

# future::plan('multisession', workers = 20)
# options(future.globals.maxSize = 50 * 1024 ^ 3) # for 50 Gb RAM
# gene.activities <- GeneActivity(combined)
# future::plan('sequential')
# 
# combined[['Activity']] <- CreateAssayObject(counts = gene.activities)

DefaultAssay(malignant) = 'Activity'
FeaturePlot(malignant, features = c('nCount_ATAC', 'nCount_Activity', 'nFeature_ATAC', 'nFeature_Activity'))
FeaturePlot(malignant, features = c('PMEL', 'MLANA', 'PTPRC', 'CD8A', 'GZMA', 'CD4', 'FOXP3', 'MS4A1', 'MRC1', 'PDGFRA'), ncol = 3)

# DotPlot(malignant, features = c('PMEL', 'MLANA', 'PTPRC', 'CD8A', 'GZMA', 'CD4', 'FOXP3', 'MS4A1', 'MRC1', 'PDGFRA')) + RotatedAxis()

DefaultAssay(malignant) <- 'ATAC'
malignant <- RunTFIDF(malignant)
malignant <- FindTopFeatures(malignant, min.cutoff = 'q90')
malignant <- RunSVD(object = malignant)

DepthCor(malignant, n = 50)

malignant <- FindNeighbors(object = malignant, reduction = 'lsi', dims = 2:50)
malignant <- FindClusters(object = malignant, algorithm = 3, resolution = 0.8)
malignant <- RunTSNE(object = malignant, reduction = 'lsi', dims = 2:50)

DefaultAssay(malignant) = 'Activity'
malignant <- NormalizeData(object = malignant, 
                           assay = 'Activity', 
                           normalization.method = 'LogNormalize',
                           # normalization.method = 'RC',
                           # scale.factor = 1e6
                           scale.factor = median(combined$nCount_Activity)
)

transfer.anchors <- FindTransferAnchors(
  reference = seu,
  query = malignant,
  reduction = 'cca'
)

predicted.labels <- TransferData(
  anchorset = transfer.anchors,
  refdata = seu$assegnazione_filtrata,
  weight.reduction = malignant[['lsi']],
  dims = 2:50
)

malignant = AddMetaData(malignant, metadata = predicted.labels$predicted.id, col.name = "predicted.id")
DefaultAssay(malignant) = 'ATAC'
DimPlot(malignant, group.by = "Patient") + 
  DimPlot(malignant, group.by = "Week", cols = week_cols) +
  DimPlot(malignant, group.by = "Cell_annotation", cols = ccolors) +
  DimPlot(malignant, group.by = "predicted.id", cols = c(colori_mp, "nc" = "grey90")) 

table(malignant$predicted.id, exclude = NULL)

DefaultAssay(malignant) = "ATAC"
CoveragePlot(malignant, region = StringToGRanges("chr14-32937952-32940741"), extend.upstream = 1500, extend.downstream = 1500, group.by = 'predicted.id')
CoveragePlot(malignant, region = "MMP16", extend.upstream = 1500, extend.downstream = 1500, group.by = 'predicted.id')


tmp = malignant$predicted.id
table(malignant$predicted.id)
combined$Metaprogram_assignment = NA
combined$Metaprogram_assignment[names(tmp)] = tmp
annotation_rna = seu$Metaprogram_assignment
annotation_rna = annotation_rna[names(annotation_rna) %in% rownames(malignant@meta.data)]
combined$Metaprogram_assignment[names(annotation_rna)] = annotation_rna
table(combined$Cell_annotation, combined$Metaprogram_assignment, exclude = NULL)

combined$Metaprogram_assignment[names(annotation_rna)] = annotation_rna

save(combined, file = "scATAC/RData/scATAC_combined.RData")