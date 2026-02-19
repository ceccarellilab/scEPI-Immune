###################
#### Libraries ####
###################
library(Seurat)
library(SeuratData)
library(dplyr)

# setwd("/home3/ciervo/scMULTIOME/Analisi/")
# dir.create("/home3/ciervo/scMULTIOME/Analisi/scRNA/RData")
source('colori_finali.R')

######################
#### Load objects ####
######################
Path_samples <- read.table("/home3/ciervo/scMULTIOME/Path_samples", quote = "'", comment.char = "")$V1 %>% as.list()
seu_list = lapply(Path_samples, function(x){
  get(load(x))
})
names(seu_list) = lapply(Path_samples, function(x){
  gsub(".RData", "", strsplit(x, split = "/")[[1]][9])
  })
rm(path, patients, weeks)

# Merge
seu = merge(seu_list[[1]], seu_list[2:13], merge.data = TRUE)
head(seu)
tail(seu)
rm(seu_list)
dim(seu)
seu$Patient = gsub(" - W0| - W0_1| - W4| - W12", "", seu$Patient)

#########################
#### Seurat workflow ####
#########################
seu = JoinLayers(seu)
seu = seu %>% 
  NormalizeData() %>% 
  FindVariableFeatures(nfeatures = 3000) %>% 
  ScaleData(features = rownames(seu)) %>% 
  RunPCA()

ElbowPlot(seu, ndims = 50) # 22

seu = seu %>% 
  FindNeighbors(dims = 1:22) %>% 
  FindClusters(resolution = 0.2) %>% 
  RunTSNE(dims = 1:22)

DimPlot(seu, reduction = "tsne", label = T) + 
  DimPlot(seu, reduction = "tsne", group.by = "Patient") + 
  DimPlot(seu, reduction = "tsne", group.by = "class") + 
  DimPlot(seu, reduction = "tsne", group.by = "Batch")
table(seu$Batch, seu$orig.ident)

seu$Patient = gsub(" - W0| - W4| - W12| - W0_1", "", seu$Patient)

DimPlot(seu, reduction = "tsne", group.by = "Patient") + 
  DimPlot(seu, reduction = 'tsne', group.by = 'Week', cols = week_cols) + 
  DimPlot(seu, reduction = "tsne", group.by = "Batch") +
  DimPlot(seu, reduction = "tsne", group.by = "Malignant", cols = colors_malignant)

seu = CellCycleScoring(seu, s.features = cc.genes.updated.2019$s.genes, g2m.features = cc.genes.updated.2019$g2m.genes)
FeaturePlot(seu, features = c("PMEL", "PTPRC"), blend = TRUE)

# Adding response metadata
seu$Responder = "NR"
seu$Responder[seu$Patient %in% c("Patient 03", "Patient 14")] = "R"

save(seu, file = "scRNA/RData/merged_object.RData")