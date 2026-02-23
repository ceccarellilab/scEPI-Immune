library(Seurat)
library(Signac)
library(dplyr)
library(ggplot2)
library(SingleCellExperiment)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomicRanges)
library(parallel)
library(SingleR)
library(celldex)
setwd("/home3/ciervo/scMULTIOME/Analisi/")

## Add metadata from snRNA-seq (TME annotation) ----
View(combined@meta.data)
combined$Barcode = gsub("Pat3", "Pat03", combined$Barcode)
combined$Barcode = gsub("Pat6", "Pat06", combined$Barcode)
combined = RenameCells(combined, new.names = combined$Barcode)

load("scRNA/RData/metadata.RData")
head(rownames(meta.data))
head(colnames(combined))

View(meta.data[meta.data$Patient %in% "Patient 03", ])
View(combined@meta.data[combined$Patient %in% "Patient 03", ])
View(combined@meta.data[combined$Patient %in% "Patient 06", ])

combined = AddMetaData(combined, metadata = meta.data[, "Cell_annotation", drop=FALSE])
View(combined@meta.data)

DimPlot(combined, group.by = c("seurat_clusters", "Patient", "Batch"))
DimPlot(combined, group.by = "Cell_annotation", cols = ccolors)
save(combined, file = "scATAC/RData/scATAC_combined.RData")

### Promoters ----
load("scATAC/RData/scATAC_combined.RData")
DefaultAssay(combined) = 'ATAC'

promoter_pmel <- StringToGRanges("chr12-55965958-55966297")
promoter_mitf = StringToGRanges('chr3-69739190-69739539')
promoter_mlana = StringToGRanges('chr9-5890707-5891049')
promoter_ptprc = StringToGRanges('chr1-198638886-198639230')
promoter_ms4a1 = StringToGRanges('chr11-60455505-60455853')
promoter_itgam = StringToGRanges('chr16-31259761-31260109')
promoter_vwf = StringToGRanges('chr12-6124697-6125029')
promoter_pdgfra = StringToGRanges('chr4-54228962-54229306')

CoveragePlot(combined, region = c('PMEL'), extend.upstream = 1500, extend.downstream = 1500, region.highlight = promoter_pmel)
CoveragePlot(combined, region = c('MITF'), extend.upstream = 1500, extend.downstream = 1500, region.highlight = promoter_mitf)
CoveragePlot(combined, region = c('MLANA'), extend.upstream = 1500, extend.downstream = 1500, region.highlight = promoter_mlana)
CoveragePlot(combined, region = c('PTPRC'), extend.upstream = 1500, extend.downstream = 1500, region.highlight = promoter_ptprc)
CoveragePlot(combined, region = c('MS4A1'), extend.upstream = 1500, extend.downstream = 1500, region.highlight = promoter_ms4a1)
CoveragePlot(combined, region = c('ITGAM'), extend.upstream = 1500, extend.downstream = 1500, region.highlight = promoter_itgam)
CoveragePlot(combined, region = c('VWF'), extend.upstream = 1500, extend.downstream = 1500, region.highlight = promoter_vwf)
CoveragePlot(combined, region = c('PDGFRA'), extend.upstream = 1500, extend.downstream = 1500, region.highlight = promoter_pdgfra)
CoveragePlot(combined, region = c('KRT14'), extend.upstream = 1500, extend.downstream = 1500)
CoveragePlot(combined, region = c('RGS5'), extend.upstream = 1500, extend.downstream = 1500)

save(combined, file = 'scATAC/RData/scATAC_combined.RData')

# Annotation of cell types ----
table(combined$Cell_annotation, combined$seurat_clusters)
combined$Cell_annotation[is.na(combined$Cell_annotation)] = "To annotate"
combined$Cell_annotation[combined$seurat_clusters %in% c(0, 1, 2, 3, 4, 5, 6, 9, 11, 12, 14, 18, 19) & 
                           combined$Cell_annotation %in% 'To annotate'] = "Malignant"
save(combined, file = "scATAC/RData/scATAC_combined.RData")

rm(promoter_inos, promoter_itgam, promoter_mitf, promoter_mlana, promoter_ms4a1, promoter_pdgfra, promoter_pmel, promoter_ptprc, promoter_vwf)

# predict labels using scRNA seq ----
load("scRNA/RData/merged_object.RData", verbose = TRUE)
dim(seu) # 58,670 cells
table(seu$Cell_annotation)

DefaultAssay(seu)

future::plan('multisession', workers = 20)
options(future.globals.maxSize = 50 * 1024 ^ 3) # for 50 Gb RAM
gene.activities <- GeneActivity(combined)
future::plan('sequential')

combined[['Activity']] <- CreateAssayObject(counts = gene.activities)
DefaultAssay(combined) = 'Activity'
combined <- NormalizeData(object = combined, 
                          assay = 'Activity', 
                          normalization.method = 'LogNormalize',
                          scale.factor = median(combined$nCount_Activity)
)
rm(gene.activities)

DefaultAssay(combined) = 'Activity'
FeaturePlot(combined, features = c('nCount_ATAC', 'nCount_Activity', 'nFeature_ATAC', 'nFeature_Activity'))
FeaturePlot(combined, features = c('PMEL', 'MLANA', 'PTPRC', 'CD8A', 'GZMA', 'CD4', 'FOXP3', 'MS4A1', 'MRC1', 'PDGFRA'), ncol = 3)

DotPlot(combined, features = c('PMEL', 'MLANA', 'PTPRC', 'CD8A', 'GZMA', 'CD4', 'FOXP3', 'MS4A1', 'MRC1', 'PDGFRA')) + RotatedAxis()

transfer.anchors <- FindTransferAnchors(
  reference = seu,
  query = combined,
  reduction = 'cca'
)

predicted.labels <- TransferData(
  anchorset = transfer.anchors,
  refdata = seu$Cell_annotation,
  weight.reduction = combined[['lsi']],
  dims = 2:50
)

combined = AddMetaData(combined, metadata = predicted.labels$predicted.id, col.name = "predicted.id")

DimPlot(combined, label = TRUE) + 
  DimPlot(combined, group.by = "Patient") + 
  DimPlot(combined, group.by = "predicted.id", cols = ccolors) + 
  DimPlot(combined, group.by = "Cell_annotation", cols = ccolors)
DimPlot(combined, group.by = "Week", cols = week_cols)

save(combined, file = "scATAC/RData/scATAC_combined.RData")

# SingleR ----
load('scATAC/RData/scATAC_combined.RData')
colnames(combined@meta.data)
DefaultAssay(combined) = 'Activity'

seu_sce <- as.matrix(combined@assays$Activity@data)

bped <- BlueprintEncodeData()
pred_bped_main <- SingleR(test = seu_sce, ref = bped, labels = bped$label.main)
pruneScores(pred_bped_main)
combined[['celltype_bped_main']] <- pred_bped_main$pruned.labels
pred_bped_fine <- SingleR(test = seu_sce, ref = bped, labels = bped$label.fine)
pruneScores(pred_bped_fine)
combined[['celltype_bped_fine']] <- pred_bped_fine$pruned.labels

table(combined$celltype_bped_fine, combined$Cell_annotation)

hpcad = celldex::HumanPrimaryCellAtlasData()
pred_hpcad_main <- SingleR(test = seu_sce, ref = hpcad, labels = hpcad$label.main)
combined[['celltype_hpcad_main']] <- pred_hpcad_main$pruned.labels
pred_hpcad_fine <- SingleR(test = seu_sce, ref = hpcad, labels = hpcad$label.fine)
combined[['celltype_hpcad_fine']] <- pred_hpcad_fine$pruned.labels

save(combined, file = 'scATAC/RData/scATAC_combined.RData')

table(combined$celltype_hpcad_fine, combined$Cell_annotation)

monaco = celldex::MonacoImmuneData()
pred_monaco_main <- SingleR(test = seu_sce, ref = monaco, labels = monaco$label.main, num.threads = 10)
combined[['celltype_monaco_main']] <- pred_monaco_main$pruned.labels
pred_monaco_fine <- SingleR(test = seu_sce, ref = monaco, labels = monaco$label.fine, num.threads = 25)
combined[['celltype_monaco_fine']] <- pred_monaco_fine$pruned.labels

table(combined$celltype_monaco_fine, combined$Cell_annotation)
DimPlot(combined, group.by = c('Cell_annotation', 'celltype_monaco_fine'))
DimPlot(combined, group.by = 'Cell_annotation', cols = ccolors)

save(combined, file = 'scATAC/RData/scATAC_combined.RData')

# TME ----
load("scATAC/RData/scATAC_combined.RData")
DefaultAssay(combined) = 'ATAC'
dim(combined) # 55,850

table(combined$Cell_annotation, exclude = NULL)

tme = subset(combined, subset = Cell_annotation %in% c('Malignant', 'Cycling cells'), invert = TRUE)
dim(tme) # 10,778

tme <- RunTFIDF(tme)
tme <- FindTopFeatures(tme, min.cutoff = 'q90')
tme <- RunSVD(object = tme)

DepthCor(tme, n = 50)

tme <- FindNeighbors(object = tme, reduction = 'lsi', dims = 2:30)
tme <- FindClusters(object = tme, algorithm = 3, resolution = 1.2)
tme <- RunTSNE(object = tme, reduction = 'lsi', dims = 2:30)

DimPlot(tme, label = TRUE) + DimPlot(tme, group.by = 'orig.ident')
DimPlot(tme, group.by = 'predicted.id', cols = ccolors) + DimPlot(tme, group.by = 'Cell_annotation', cols = ccolors)
DimPlot(tme, group.by = 'Batch') + DimPlot(tme, group.by = 'Patient')

DimPlot(tme, label = T) + DimPlot(tme, group.by = 'orig.ident') +
  DimPlot(tme, group.by = 'Cell_annotation', cols = ccolors) + DimPlot(tme, group.by = 'celltype_monaco_fine') + NoLegend()

table(tme$celltype_bped_fine, tme$ATAC_snn_res.1.2)
tme$celltype_bped_main[tme$ATAC_snn_res.1.2 %in% 0] = 'T cells - CD8'
tme$celltype_bped_main[tme$ATAC_snn_res.1.2 %in% 1] = 'Macrophages'
tme$celltype_bped_main[tme$ATAC_snn_res.1.2 %in% 2] = 'T cells - CD8'
tme$celltype_bped_main[tme$ATAC_snn_res.1.2 %in% 3] = 'Macrophages'
tme$celltype_bped_main[tme$ATAC_snn_res.1.2 %in% 4] = 'Melanocytes'
tme$celltype_bped_main[tme$ATAC_snn_res.1.2 %in% 5] = 'T cells - CD4'
tme$celltype_bped_main[tme$ATAC_snn_res.1.2 %in% 6] = 'Endothelial cells'
tme$celltype_bped_main[tme$ATAC_snn_res.1.2 %in% 7] = 'Macrophages'
tme$celltype_bped_main[tme$ATAC_snn_res.1.2 %in% 8] = 'T cells - CD4'
tme$celltype_bped_main[tme$ATAC_snn_res.1.2 %in% 9] = 'Macrophages'
tme$celltype_bped_main[tme$ATAC_snn_res.1.2 %in% 10] = 'CAF'
tme$celltype_bped_main[tme$ATAC_snn_res.1.2 %in% 11] = 'T cells - CD8'
tme$celltype_bped_main[tme$ATAC_snn_res.1.2 %in% 12] = 'Macrophages'
tme$celltype_bped_main[tme$ATAC_snn_res.1.2 %in% 13] = 'T cells - Tregs'
tme$celltype_bped_main[tme$ATAC_snn_res.1.2 %in% 14] = 'T cells - CD8'
tme$celltype_bped_main[tme$ATAC_snn_res.1.2 %in% 15] = 'Macrophages'
tme$celltype_bped_main[tme$ATAC_snn_res.1.2 %in% 16] = 'Plasma cells'
tme$celltype_bped_main[tme$ATAC_snn_res.1.2 %in% 17] = 'T cells - CD8'
tme$celltype_bped_main[tme$ATAC_snn_res.1.2 %in% 18] = 'CAF'
tme$celltype_bped_main[tme$ATAC_snn_res.1.2 %in% 19] = 'B cells'
tme$celltype_bped_main[tme$ATAC_snn_res.1.2 %in% 20] = 'CAF'
tme$celltype_bped_main[tme$ATAC_snn_res.1.2 %in% 21] = 'Keratinocytes'
tme$celltype_bped_main[tme$ATAC_snn_res.1.2 %in% 22] = 'Macrophages'

DimPlot(tme, group.by = 'celltype_bped_main')

DimPlot(tme, label = T)
table(tme$celltype_monaco_fine, tme$ATAC_snn_res.1.2)
tme$celltype_monaco_main[tme$ATAC_snn_res.1.2 %in% 0] = 'T cells - CD8'
tme$celltype_monaco_main[tme$ATAC_snn_res.1.2 %in% 1] = 'Macrophages'
tme$celltype_monaco_main[tme$ATAC_snn_res.1.2 %in% 2] = 'T cells - CD8'
tme$celltype_monaco_main[tme$ATAC_snn_res.1.2 %in% 3] = 'T cells - CD8'
tme$celltype_monaco_main[tme$ATAC_snn_res.1.2 %in% 4] = 'Melanocytes'
tme$celltype_monaco_main[tme$ATAC_snn_res.1.2 %in% 5] = 'T cells - CD4'
tme$celltype_monaco_main[tme$ATAC_snn_res.1.2 %in% 6] = 'Endothelial cells'
tme$celltype_monaco_main[tme$ATAC_snn_res.1.2 %in% 7] = 'Macrophages'
tme$celltype_monaco_main[tme$ATAC_snn_res.1.2 %in% 8] = 'T cells - CD4'
tme$celltype_monaco_main[tme$ATAC_snn_res.1.2 %in% 9] = 'Macrophages'
tme$celltype_monaco_main[tme$ATAC_snn_res.1.2 %in% 10] = 'CAF'
tme$celltype_monaco_main[tme$ATAC_snn_res.1.2 %in% 11] = 'T cells - CD8'
tme$celltype_monaco_main[tme$ATAC_snn_res.1.2 %in% 12] = 'Macrophages'
tme$celltype_monaco_main[tme$ATAC_snn_res.1.2 %in% 13] = 'T cells - Tregs'
tme$celltype_monaco_main[tme$ATAC_snn_res.1.2 %in% 14] = 'T cells - CD8'
tme$celltype_monaco_main[tme$ATAC_snn_res.1.2 %in% 15] = 'Macrophages'
tme$celltype_monaco_main[tme$ATAC_snn_res.1.2 %in% 16] = 'Plasma cells'
tme$celltype_monaco_main[tme$ATAC_snn_res.1.2 %in% 17] = 'T cells - CD8'
tme$celltype_monaco_main[tme$ATAC_snn_res.1.2 %in% 18] = 'CAF'
tme$celltype_monaco_main[tme$ATAC_snn_res.1.2 %in% 19] = 'B cells'
tme$celltype_monaco_main[tme$ATAC_snn_res.1.2 %in% 20] = 'CAF'
tme$celltype_monaco_main[tme$ATAC_snn_res.1.2 %in% 21] = 'Keratinocytes'
tme$celltype_monaco_main[tme$ATAC_snn_res.1.2 %in% 22] = 'Macrophages'

DimPlot(tme, label = T, group.by = 'celltype_monaco_main') + DimPlot(tme, label = T, group.by = 'celltype_bped_main')

table(tme$celltype_bped_main, tme$ATAC_snn_res.1.2)

tme$celltype_bped_main[tme$ATAC_snn_res.1.2 %in% c(0, 2, 11, 14, 17)] = 'CD8'
tme$celltype_bped_main[tme$ATAC_snn_res.1.2 %in% c(5, 8, 13)] = 'CD4'
tme$celltype_bped_main[tme$ATAC_snn_res.1.2 %in% c(19)] = 'B cells'
tme$celltype_bped_main[tme$ATAC_snn_res.1.2 %in% c(6)] = 'Endothelial cells'
tme$celltype_bped_main[tme$ATAC_snn_res.1.2 %in% c(21)] = 'Keratinocytes'
tme$celltype_bped_main[tme$ATAC_snn_res.1.2 %in% c(10, 20, 18)] = 'Fibroblasts'
tme$celltype_bped_main[tme$ATAC_snn_res.1.2 %in% c(16)] = 'Plasma cells'
tme$celltype_bped_main[tme$ATAC_snn_res.1.2 %in% c(1, 3, 7, 9, 12, 15, 22)] = 'Macrophages'
tme$celltype_bped_main[tme$ATAC_snn_res.1.2 %in% c(4)] = 'Malignant'

DimPlot(tme, group.by = 'celltype_bped_main', cells = names(tmp))

table(tme$celltype_bped_fine, tme$ATAC_snn_res.1.2)

# tme$Cell_annotation[tme$Cell_annotation %in% 'T cells' & tme$ATAC_snn_res.1.2 %in% c(2, 11)] = 'T cells - CD8'
# tme$Cell_annotation[tme$Cell_annotation %in% 'T cells' & tme$ATAC_snn_res.1.2 %in% c(0, 17)] = 'T cells - CD8 Tcm'
# tme$Cell_annotation[tme$Cell_annotation %in% 'T cells' & tme$ATAC_snn_res.1.2 %in% c(5)] = 'T cells - CD4'
# tme$Cell_annotation[tme$Cell_annotation %in% 'T cells' & tme$ATAC_snn_res.1.2 %in% c(8)] = 'T cells - CD4 Tcm'
# tme$Cell_annotation[tme$Cell_annotation %in% 'T cells' & tme$ATAC_snn_res.1.2 %in% c(5, 8)] = 'T cells - CD4'

x <- subset(tme, celltype_bped_main %in% 'Malignant', invert = TRUE)
x <- RunTFIDF(x)
x <- FindTopFeatures(x, min.cutoff = 'q90')
x <- RunSVD(object = x)

DepthCor(x, n = 50)

x <- FindNeighbors(object = x, reduction = 'lsi', dims = 2:30)
x <- FindClusters(object = x, algorithm = 3, resolution = 1.2)
x <- RunTSNE(object = x, reduction = 'lsi', dims = 2:30)
DimPlot(x) + DimPlot(x, group.by = 'celltype_bped_main')
VlnPlot(tme, features = c('nCount_ATAC', 'nFeature_ATAC'), ncol = 1)

DimPlot(x, cells = names(tmp), label = TRUE) + DimPlot(x, group.by = 'Cell_annotation', cells = names(tmp), cols = ccolors) +
  DimPlot(x, group.by = 'orig.ident', cells = names(tmp))

table(x$Cell_annotation, x$ATAC_snn_res.1.2)
x$Cell_annotation[x$Cell_annotation %in% c("Macrophages", "Undefined", "CAF", "T cells") & x$ATAC_snn_res.1.2 %in% c(1, 3, 5, 11, 14, 22)] = "M1 macrophages"
x$Cell_annotation[x$Cell_annotation %in% c("Macrophages", "Undefined", "CAF", "T cells") & x$ATAC_snn_res.1.2 %in% c(8)] = "M2 macrophages"
x$Cell_annotation[x$Cell_annotation %in% c("Macrophages", "Undefined", "CAF", "T cells") & x$ATAC_snn_res.1.2 %in% c(9, 20)] = "iCAF"
x$Cell_annotation[x$Cell_annotation %in% c("Macrophages", "Undefined", "CAF", "T cells") & x$ATAC_snn_res.1.2 %in% c(17)] = "myoCAF"
x$Cell_annotation[x$Cell_annotation %in% c("Macrophages", "Undefined", "CAF", "T cells") & x$ATAC_snn_res.1.2 %in% c(0, 16, 18)] = "T cells - CD8 Tcm"
x$Cell_annotation[x$Cell_annotation %in% c("Macrophages", "Undefined", "CAF", "T cells") & x$ATAC_snn_res.1.2 %in% c(2, 10)] = "T cells - CD8"
x$Cell_annotation[x$Cell_annotation %in% c("Macrophages", "Undefined", "CAF", "T cells") & x$ATAC_snn_res.1.2 %in% c(4)] = "T cells - CD4"
x$Cell_annotation[x$Cell_annotation %in% c("Macrophages", "Undefined", "CAF", "T cells") & x$ATAC_snn_res.1.2 %in% c(7)] = "T cells - CD4 Tcm"
x$Cell_annotation[x$Cell_annotation %in% c("Macrophages", "Undefined", "CAF", "T cells") & x$ATAC_snn_res.1.2 %in% c(13)] = "T cells - Exhausted T cell"
x$Cell_annotation[x$Cell_annotation %in% c("Macrophages", "Undefined", "CAF", "T cells") & x$ATAC_snn_res.1.2 %in% c(12)] = "T cells - Tregs"
x$Cell_annotation[x$Cell_annotation %in% c("Macrophages", "Undefined", "CAF", "T cells") & x$ATAC_snn_res.1.2 %in% c(19)] = "B cells"
x$Cell_annotation[x$Cell_annotation %in% c("Macrophages", "Undefined", "CAF", "T cells") & x$ATAC_snn_res.1.2 %in% c(6)] = "Endothelial cells"
x$Cell_annotation[x$Cell_annotation %in% c("Macrophages", "Undefined", "CAF", "T cells") & x$ATAC_snn_res.1.2 %in% c(21)] = "Keratinocytes"
x$Cell_annotation[x$Cell_annotation %in% c("Macrophages", "Undefined", "CAF", "T cells") & x$ATAC_snn_res.1.2 %in% c(15)] = "Plasma cells"

x$Cell_annotation[x$Cell_annotation %in% c(1, 3, 5, 14)] = 'M2 macrophages'
x$Cell_annotation[x$Cell_annotation %in% c(8, 11)] = 'M1 macrophages'

x = AddMetaData(x, metadata = tmp, col.name = 'Cell_annotation')

DimPlot(x, group.by = "Cell_annotation", cols = ccolors)

table(x$Cell_annotation, x$orig.ident) 
table(x$Cell_annotation)
tmp = x$Cell_annotation[x$Cell_annotation %in% c('M1 macrophages', 'M2 macrophages')]

# DefaultAssay(x) = 'Activity'
# markers.genescore = FindAllMarkers(x, only.pos = TRUE)
# markers.genescore = markers.genescore[markers.genescore$p_val_adj < 0.05, ]
# markers.genescore %>% group_by(cluster) %>%
#   dplyr::filter(avg_log2FC > 1) %>%
#   slice_head(n = 10) %>%
#   ungroup() -> top10

malignant = tme$Barcode[tme$celltype_bped_main %in% "Malignant"]
cell_annot = x$Cell_annotation

# combined$Cell_annotation[names(malignant)] = "Malignant"
# combined$Cell_annotation[names(cell_annot)] = cell_annot

combined = AddMetaData(combined, metadata = tmp, col.name = 'Cell_annotation')

DimPlot(tme, group.by = 'Cell_annotation', cols = ccolors)

DimPlot(combined, group.by = "Cell_annotation", cols = ccolors)
save(combined, file = "scATAC/RData/scATAC_combined.RData")

DefaultAssay(combined) = 'ATAC'

DimPlot(combined, group.by = 'Patient') +
  DimPlot(combined, group.by = "Week", cols = week_cols)

DimPlot(combined, group.by = "Cell_annotation", cols = ccolors)

promoter_pmel <- StringToGRanges("chr12-55965958-55966297")
promoter_mitf = StringToGRanges('chr3-69739190-69739539')
# promoter_axl = StringToGRanges('chr19-41219054-41219327')
promoter_mlana = StringToGRanges('chr9-5890707-5891049')
promoter_inos = StringToGRanges('chr17-27793684-27794270')
promoter_ptprc = StringToGRanges('chr1-198638886-198639230')
promoter_ms4a1 = StringToGRanges('chr11-60455505-60455853')
# M1
promoter_cd86 = StringToGRanges('chr1-206772302-206772648')
# M2
promoter_itgam = StringToGRanges('chr16-31259761-31260109')

promoter_vwf = StringToGRanges('chr12-6124697-6125029')
promoter_pdgfra = StringToGRanges('chr4-54228962-54229306')
promoter_krt14 = StringToGRanges("chr17-41586835-41587180")
promoter_mzb1 = StringToGRanges("chr5-139389796-139390129")
promoter_rgs5 = StringToGRanges("chr1-163203032-163203373")

CoveragePlot(combined, region = c('PMEL'), region.highlight = promoter_pmel, extend.upstream = 1500, extend.downstream = 1500, group.by = 'Cell_annotation', peaks = TRUE)
CoveragePlot(combined, region = c('MITF'), extend.upstream = 1500, extend.downstream = 1500, region.highlight = promoter_mitf, group.by = 'Cell_annotation')
# CoveragePlot(combined, region = c('AXL'), extend.upstream = 1500, extend.downstream = 1500, region.highlight = promoter_axl, group.by = "Cell_annotation")
CoveragePlot(combined, region = c('MLANA'), extend.upstream = 1500, extend.downstream = 1500, region.highlight = promoter_mlana, group.by = "Cell_annotation")
CoveragePlot(combined, region = c('PTPRC'), extend.upstream = 1500, extend.downstream = 1500, region.highlight = promoter_ptprc, group.by = "Cell_annotation")
CoveragePlot(combined, region = c('MS4A1'), extend.upstream = 1500, extend.downstream = 1500, region.highlight = promoter_ms4a1, group.by = "Cell_annotation")
CoveragePlot(combined, region = c('MZB1'), extend.upstream = 1500, extend.downstream = 1500, region.highlight = promoter_mzb1, group.by = "Cell_annotation")

CoveragePlot(combined, region = c('IL10'), extend.upstream = 1500, extend.downstream = 1500, region.highlight = promoter_cd86, group.by = "Cell_annotation")
CoveragePlot(combined, region = c('ITGAM'), extend.upstream = 1500, extend.downstream = 1500, region.highlight = promoter_itgam, group.by = "Cell_annotation")

tmp = subset(combined, orig.ident %in% 'Pat15_W12')
CoveragePlot(tmp, region = c('MITF'), extend.upstream = 1500, extend.downstream = 1500, group.by = "Cell_annotation")
CoveragePlot(tmp, region = c('KRT14'), extend.upstream = 1500, extend.downstream = 1500, region.highlight = promoter_krt14, group.by = "Cell_annotation")

CoveragePlot(tmp, region = c('GZMA'), extend.upstream = 1500, extend.downstream = 1500, group.by = "Cell_annotation")

CoveragePlot(combined, region = c('KRT14'), extend.upstream = 1500, extend.downstream = 1500, region.highlight = promoter_krt14, group.by = "Cell_annotation")
CoveragePlot(combined, region = c('VWF'), extend.upstream = 1500, extend.downstream = 1500, region.highlight = promoter_vwf, group.by = "Cell_annotation")
CoveragePlot(combined, region = c('PDGFRA'), extend.upstream = 1500, extend.downstream = 1500, region.highlight = promoter_pdgfra, group.by = "Cell_annotation")
CoveragePlot(combined, region = c('RGS5'), extend.upstream = 1500, extend.downstream = 1500, region.highlight = promoter_rgs5, group.by = "Cell_annotation")
# CoveragePlot(combined, region = c('NOS2'), extend.upstream = 1500, extend.downstream = 1500, region.highlight = promoter_inos)
CoveragePlot(combined, region = c('CD4'), extend.upstream = 1500, extend.downstream = 1500, group.by = "Cell_annotation", region.highlight = StringToGRanges("chr12-6789331-6789678"))
CoveragePlot(combined, region = c('CD8A'), extend.upstream = 1500, extend.downstream = 1500, group.by = "Cell_annotation", region.highlight = StringToGRanges("chr2-86790816-86791161"))
CoveragePlot(combined, region = c('FOXP3'), extend.upstream = 1500, extend.downstream = 1500, group.by = "Cell_annotation", region.highlight = StringToGRanges("chrX-49264658-49265000"))
CoveragePlot(combined, region = c('LAG3'), extend.upstream = 1500, extend.downstream = 1500, group.by = "Cell_annotation", region.highlight = StringToGRanges("chr12-6772261-6772587"))

DefaultAssay(combined) = "Activity"
DotPlot(combined, features = c("PMEL", "MITF", "MLANA",
                               "PTPRC", "MS4A1", "MZB1",
                               "ITGAM", "VWF", "PDGFRA", "RGS5",
                               "CD4", "CD8A", "FOXP3", "LAG3"), group.by = "Cell_annotation") + RotatedAxis()
DefaultAssay(combined) = "ATAC"
save(combined, file = 'scATAC/RData/scATAC_combined.RData')