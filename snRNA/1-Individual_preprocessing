###################
#### Libraries ####
###################
library(Seurat)
library(dplyr)
library(SingleCellExperiment)
library(SCEVAN)

# setwd("/home/ciervo/EPICA/NIBIT/")
# dir.create("Individual_analysis")
# dir.create("Individual_analysis/Patient_02")
# dir.create("Individual_analysis/Patient_02/Week0")

setwd("Individual_analysis/Patient_02/Week0/")

######################
#### Load objects ####
######################
seu = Read10X(data.dir =  "/home3/adefalco/scRNA+ATACT_NIBIT/Batch_1/Pat02_W0/outs/filtered_feature_bc_matrix")$`Gene Expression`
seu = CreateSeuratObject(seu, min.features = 500, min.cells = 3, project = "Pat02_W0")
seu = RenameCells(seu, add.cell.id = "Pat02_W0")
dim(seu)

seu$percent.mt = PercentageFeatureSet(seu, pattern = "^MT-")
seu$Patient = "Patient 02 - W0"
seu$Week = "Week 0"
seu$Batch = "Batch 1"

###################
#### Filtering ####
###################
seu = subset(seu, subset = nCount_RNA < 50000 & 
               percent.mt < 5 &
               nFeature_RNA < 7500)
dim(seu) 
writeLines(rownames(seu@meta.data), con = "allCells.txt")

#################
#### SCEVAN  ####
#################
cells = read.table("allCells.txt")$V1
cells = gsub("Pat02_W0_", "", cells)

dir.create("SCEVAN_results")
setwd("SCEVAN_results")
counts <- Read10X(data.dir = "/home3/adefalco/scRNA+ATACT_NIBIT/Batch_1/Pat02_W0/outs/filtered_feature_bc_matrix")$`Gene Expression`
counts <- counts[, cells]
res_SCEVAN <-  pipelineCNA(counts, SUBCLONES = TRUE, par_cores = 60)
save(res_SCEVAN, file = "Pat02_W0.RData")
 
# Add SCEVAN results
setwd("/home/ciervo/EPICA/NIBIT/Individual_analysis/Patient_02/Week0/")
load("SCEVAN_results/Pat02_W0.RData", verbose = TRUE)
rownames(res_SCEVAN) = paste0("Pat02_W0_", rownames(res_SCEVAN))

seu = AddMetaData(seu, res_SCEVAN)
rm(res_SCEVAN)

##########################
#### Seurat workflow  ####
##########################
seu = seu %>% 
  NormalizeData() %>% 
  FindVariableFeatures() %>% 
  ScaleData(features = rownames(seu)) %>% 
  RunPCA()
ElbowPlot(seu, ndims = 50)

seu = seu %>% 
  FindNeighbors(dims = 1:13) %>% 
  RunTSNE(dims = 1:13) %>% 
  # RunUMAP(dims = 1:13) %>% 
  FindClusters(resolution = 0.1)

##############
#### Plot ####
##############
png("tsne.png", width = 5, height = 5, units = "in", res = 600)
DimPlot(seu, reduction = "tsne")
dev.off()

png("markers.png", width = 10, height = 5, units = "in", res = 600)
FeaturePlot(seu, features = c("PMEL", "MITF", "PTPRC", "CD8A"), reduction = "tsne")
dev.off()

DimPlot(seu)
FeaturePlot(seu, features = "percent.mt")

png("scevan.png", width = 10, height = 5, units = "in", res = 600)
DimPlot(seu, reduction = "tsne") + DimPlot(seu, group.by = "class", reduction = "tsne", cols = c("grey80", "forestgreen", "red"))
dev.off() 

save(seu, file = "Pat02_W0.RData")
