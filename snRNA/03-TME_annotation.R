###################
#### Libraries ####
###################
library(Seurat)
library(SeuratData)
library(dplyr)
library(parallel)
# setwd("/home3/ciervo/scMULTIOME/Analisi/")
# dir.create("/home3/ciervo/scMULTIOME/Analisi/scRNA/RData")
source("/home/caruso/scProject/NYnontumor/package_code/functionScRNAseq.R")
source('colori_finali.R')
source("/home/ciervo/Functions/mww_celltype.R")

############################
#### Load merged object ####
############################
load("scRNA/RData/merged_object.RData")
tme = subset(seu, Malignant %in% "Non malignant")
rm(seu)

tme = tme %>% 
  NormalizeData() %>% 
  FindVariableFeatures(nfeatures = 3000) %>% 
  ScaleData(features = rownames(tme)) %>% 
  RunPCA()

ElbowPlot(tme, ndims = 50) # 15

tme = tme %>% 
  FindNeighbors(dims = 1:15) %>% 
  FindClusters(resolution = 0.1) %>% 
  RunTSNE(dims = 1:15)

DimPlot(tme, reduction = "tsne", label = T) + 
  DimPlot(tme, reduction = "tsne", group.by = "Patient") + 
  DimPlot(tme, reduction = "tsne", group.by = "Batch")
FeaturePlot(tme, features = c("COL1A1", "FOXP3", "CD8A", "CD4", "KRT14", "PECAM1", "MS4A1", "MZB1"), label = TRUE)

##############################
#### DEGs Seurat clusters ####
##############################
future::plan("multisession", workers = 20)
markers = FindAllMarkers(tme, only.pos = TRUE)
future::plan("sequential")

markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10
DoHeatmap(tme, features = top10$gene) + NoLegend()

tme$Cell_annotation = "N/A"
tme$Cell_annotation[tme$RNA_snn_res.0.1 %in% c(0)] = "T cells"
tme$Cell_annotation[tme$RNA_snn_res.0.1 %in% c(1, 2)] = "Myeloid cells"
tme$Cell_annotation[tme$RNA_snn_res.0.1 %in% c(3)] = "Endothelial cells"
tme$Cell_annotation[tme$RNA_snn_res.0.1 %in% c(4)] = "Keratinocytes"
tme$Cell_annotation[tme$RNA_snn_res.0.1 %in% c(5, 7)] = "CAF"
tme$Cell_annotation[tme$RNA_snn_res.0.1 %in% c(9)] = "Plasma cells"
tme$Cell_annotation[tme$RNA_snn_res.0.1 %in% c(6, 8)] = "Malignat"
tme$Cell_annotation[tme$RNA_snn_res.0.1 %in% c(10)] = "B cells"
table(tme$Cell_annotation)

FeatureScatter(tme, feature1 = "S.Score", "G2M.Score") + 
  geom_vline(xintercept = 0.2, linetype = "dashed", color = "grey70") + 
  geom_hline(yintercept = 0.2, linetype = "dashed", color = "grey70")
tme$Cell_annotation[tme$S.Score > 0.2 | tme$G2M.Score > 0.2] = "Cycling cells"
cell_annotation = tme$Cell_annotation
save(cell_annotation, file = "scRNA/RData/tme_cell_annotation.RData")

#####################################################
#### Adding main annotation to the merged object ####
#####################################################
load("scRNA/RData/merged_object.RData")
load("scRNA/RData/tme_cell_annotation.RData")

cell_annotation[cell_annotation %in% "Malignat"] = "Malignant"
seu = AddMetaData(seu, metadata = cell_annotation, col.name = "Cell_annotation")
seu$Cell_annotation[is.na(seu$Cell_annotation)] = seu$Malignant[is.na(seu$Cell_annotation)]

DimPlot(seu, group.by = "Cell_annotation")
save(seu, file = "scRNA/RData/merged_object.RData")

#################################
#### Refining TME annotation ####
#################################
load("scRNA/RData/merged_object.RData")
tme = subset(seu, Cell_annotation %in% c("Malignant", "Cycling cells"), invert = TRUE)
dim(tme) # 31,093 x 4,783
rm(seu)

tme = tme %>% 
  NormalizeData() %>% 
  FindVariableFeatures(nfeatures = 3000) %>% 
  ScaleData(features = rownames(tme)) %>% 
  RunPCA()
ElbowPlot(tme, ndims = 50) # 12

tme = tme %>% 
  FindNeighbors(dims = 1:12) %>% 
  RunTSNE(dims = 1:12)

# Identify myeliod cell types
tme = FindClusters(tme, resolution = 0.5)
DimPlot(tme, reduction = "tsne", label = T) | DimPlot(tme, group.by = "Cell_annotation")

clusters = c(1, 3, 4, 6, 12)
marker.list <- mclapply(clusters, FUN = function(x) {
  markers <- FindMarkers(tme, ident.1 = x, assay = "RNA")
  markers <- markers[order(markers$avg_log2FC, decreasing = TRUE), , drop = FALSE]
  genes <- markers$avg_log2FC
  names(genes) <- rownames(markers)
  return(genes)
},
mc.cores = length(clusters)
)
names(marker.list) <- paste0("markers.cluster", clusters)
lapply(marker.list, function(x) names(head(x, 20)))

# MWW on ranked list
load("/home/caruso/Analisi2021/signatures/scTHI_c8_signatures_968.RData", verbose = T)
signature_Colors = signature_Colors[signature_Colors$ALLPhenotypeFinal %in% c("MacrophagesM1", "MacrophagesM2"), ]
signature = signature[signature_Colors$ALLPhenotypeOri]
source("/home/ciervo/EPICA/NIBIT/signature_macrophages.R") # Macrophages signatures

signature = c(signature, macrophages)

tmp = MultirankedMww(rankedLists = marker.list, geneSet = signature, minLenGeneSet = 1, ncore = 32)
tmp = lapply(tmp, FUN = function(x) x[x$qValue < 0.05 & x$logit2_NES > 0.58, ])
tmp = lapply(tmp, FUN = function(x) x[order(x$logit2_NES, decreasing = T), ])
lapply(tmp, function(x) rownames(x[1:4, ]))

tme$Cell_annotation[tme$RNA_snn_res.0.5 %in% c(4, 3, 12)] = "M1 macrophages"
tme$Cell_annotation[tme$RNA_snn_res.0.5 %in% c(1, 6)] = "M2 macrophages"


DimPlot(tme, group.by = "Cell_annotation", cols = ccolors) + ggtitle('Cell annotation')
rm(macrophages, marker.list, signature_Colors, signature, tmp)

# tme = AddModuleScore(tme, features = list(M1 = c("TNF", "NOS2", "IL1B", "IL6", "CXCL9", "CXCL10", "CD80", "CD86", "HLA-DRA", "PTGS2"),
#                                           M2 = c("ARG1", "CD206", "CD163", "IL10", "TGFB1", "CCL17", "CCL18", "CCL22", "PPARG", "RETNLA"))
# )

# CAFs are forming different clusters
DimPlot(tme, reduction = "tsne", label = T) + 
  DimPlot(tme, group.by = "Cell_annotation", reduction = "tsne")
tme = FindClusters(tme, resolution = 0.7) 

clusters = c(8, 9, 11)
marker.list <- mclapply(clusters, FUN = function(x) {
  markers <- FindMarkers(tme, ident.1 = x, assay = "RNA")
  markers <- markers[order(markers$avg_log2FC, decreasing = TRUE), , drop = FALSE]
  genes <- markers$avg_log2FC
  names(genes) <- rownames(markers)
  return(genes)
},
mc.cores = length(clusters)
)
names(marker.list) <- paste0("markers.cluster", clusters)
lapply(marker.list, function(x) names(head(x, 20)))
FeaturePlot(tme, features = c("PDGFRA", "ACTA2", "RGS5"))

tme$Cell_annotation[tme$RNA_snn_res.0.7 %in% 9] = "myoCAF"
tme$Cell_annotation[tme$RNA_snn_res.0.7 %in% c(8, 11)] = "iCAF"

DimPlot(tme, group.by = "Cell_annotation")

# T cells
DimPlot(tme, reduction = "tsne", label = T) + 
  DimPlot(tme, group.by = "Cell_annotation", reduction = "tsne")

tme = FindClusters(tme, resolution = 5) # increasing resolution to distinguish T cell subtypes
DimPlot(tme, reduction = "tsne", label = T)

clusters <- c(0, 1, 2, 4, 5, 6, 7)
marker.list <- mclapply(clusters, FUN = function(x) {
  markers <- FindMarkers(tme, ident.1 = x, assay = "RNA")
  markers <- markers[order(markers$avg_log2FC, decreasing = TRUE), , drop = FALSE]
  genes <- markers$avg_log2FC
  names(genes) <- rownames(markers)
  return(genes)
  },
  mc.cores = length(clusters)
)
names(marker.list) <- paste0("markers.cluster", clusters)

source("/home/caruso/scProject/NYnontumor/package_code/functionScRNAseq.R")
load("/home/caruso/Analisi2021/signatures/scTHI_c8_signatures_968.RData", verbose = T)
signature_Colors = signature_Colors[signature_Colors$ALLPhenotypeFinal %in% c("AnergicTcell", "CD4", "CD8", "ExhaustedTcell",
                                                                              "NaturalKiller", "TgammaDelta", "Th1", "Tcell",
                                                                              "Tcm", "TDoubleNeg", "Tem", "TFH", "TgammaDelta",
                                                                              "Th17", "Th2", "Th22", "Th9", "Tregs",
                                                                              "TNaturalKiller"), ]
signature_Colors = signature_Colors[1:90, ]
signature = signature[signature_Colors$ALLPhenotypeOri]
signature = signature[-grep("xCellOld", names(signature))]

tmp = MultirankedMww(rankedLists = marker.list, geneSet = signature, minLenGeneSet = 5, ncore = 32)
tmp = lapply(tmp, FUN = function(x) x[x$qValue < 0.1, ])
tmp = lapply(tmp, FUN = function(x) x[order(x$logit2_NES, decreasing = T), ])
lapply(tmp, function(x) rownames(x[1:4, ]))

tme$Cell_annotation[tme$RNA_snn_res.5 %in% 0] = "T cells - CD8"
tme$Cell_annotation[tme$RNA_snn_res.5 %in% 1] = "T cells - CD4 Tcm"
tme$Cell_annotation[tme$RNA_snn_res.5 %in% 4] = "T cells - Exhausted T cell"
tme$Cell_annotation[tme$RNA_snn_res.5 %in% 5] = "T cells - CD4"
tme$Cell_annotation[tme$RNA_snn_res.5 %in% 6] = "T cells - Tregs"
tme$Cell_annotation[tme$RNA_snn_res.5 %in% 7] = "T cells - CD8 Tcm"

DimPlot(tme, group.by = "Cell_annotation", reduction = "tsne")
# save(tme, file = "scRNA/RData/tme_subset.RData")

cell_annotation = tme$Cell_annotation
save(cell_annotation, file = "scRNA/RData/fine_cell_annotation.RData")

####################################################
#### Adding refined cell types to seurat object ####
####################################################
load("scRNA/RData/merged_object.RData")

tme = subset(seu, Cell_annotation %in% c("Malignant", "Cycling cells"), invert = TRUE)
rm(seu, cell_annotation)

tme = tme %>% 
  NormalizeData() %>% 
  FindVariableFeatures(nfeatures = 3000) %>% 
  ScaleData(features = rownames(tme)) %>%
  RunPCA() # 12

tme = tme %>% 
  FindNeighbors(dims = 1:12) %>% 
  RunTSNE(dims = 1:12)

llevels = names(colors_tme)
tme$Cell_annotation = factor(tme$Cell_annotation, levels = llevels)
DimPlot(tme, group.by = "Cell_annotation", cols = colors) + ggtitle("Cell type")

FeaturePlot(tme, features = c("PDGFRA", "ACTA2", "RGS5", # iCAF; myoCAF
                              "PTPRC", "CD8A", "CD4", "GZMA", "FOXP3", # T cells  
                              "ITGAM", "MRC1",  # Myeloid lineage 
                              "KRT14", "EPCAM", # Keratinocytes
                              "PECAM1", # Endothelial cells
                              "MS4A1", "MZB1" # Plasma cells; # B cells
), ncol = 3)

names(ccolors)[!names(ccolors) %in% c('Malignant', 'Cycling cells')]
# tme$Cell_annotation = factor(tme$Cell_annotation, names(ccolors)[!names(ccolors) %in% c('Malignant', 'Cycling cells')])
DotPlot(tme, features = c('PDGFRA', 'ACTA2', 'RGS5', 
                          'PTPRC', 'CD8A', 'CD4', 'GZMA', 'FOXP3', 'IL2RA', 
                          'NCR1', 'GNLY', 
                          'LAG3', 'CTLA4', 'HAVCR2', 'PDCD1', 
                          'ITGAM', 'MRC1', 
                          'KRT14', 
                          'PECAM1', 'VWF', 
                          'MS4A1', 'MZB1'), group.by = 'Cell_annotation') + RotatedAxis()
seu$Cell_annotation[is.na(seu$Cell_annotation)] = seu$Malignant[is.na(seu$Cell_annotation)]
save(seu, file = "scRNA/RData/merged_object.RData")

