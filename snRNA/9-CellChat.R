###################
#### Libraries ####
###################
library(Seurat)
library(dplyr)
library(ggplot2)
library(CellChat)
library(ComplexHeatmap)
library(circlize)

setwd("/home3/ciervo/scMULTIOME/Analisi/")
source("colori_finali.R")

load('scRNA/RData/merged_object.RData')
table(seu$Metaprogram_assignment, exclude = NULL)
table(seu$Metaprogram_assignment, seu$Cell_annotation, exclude = NULL)
seu$CellChat = seu$Metaprogram_assignment
table(seu$CellChat, exclude = NULL)

seu$CellChat[is.na(seu$CellChat)] = seu$Cell_annotation[is.na(seu$CellChat)]
table(seu$CellChat, seu$Cell_annotation)

sort(table(seu$CellChat, exclude = NULL))
seu = subset(seu, CellChat %in% c('Cycling cells', "non-classified"), invert = TRUE)
dim(seu)

table(seu$Responder, exclude = NULL, seu$orig.ident)

# seu = subset(seu, Responder %in% 'R')
seu = subset(seu, Responder %in% 'NR')
table(seu$CellChat)

seu = NormalizeData(seu)

data.input <- seu[["RNA"]]$data # normalized data matrix
labels <- seu$CellChat
table(labels, exclude = NULL)

meta <- data.frame(samples = seu$orig.ident, labels = seu$CellChat, row.names = rownames(seu@meta.data)) # create a dataframe of the cell labels
meta$samples = factor(meta$samples)
meta$labels = factor(meta$labels)

data.input = data.input[, rownames(meta)]

# Creating CellChat object ----
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels")
cellchat <- setIdent(cellchat, ident.use = "labels")
levels(cellchat@idents) 
groupSize <- as.numeric(table(cellchat@idents)) 
groupSize

CellChatDB <- CellChatDB.human 
showDatabaseCategory(CellChatDB)
cellchat@DB <- CellChatDB

cellchat <- subsetData(cellchat)

# dir.create('scRNA/CellChat')

#### Pre-processing ####
future::plan("multisession", workers = 50)
options(future.globals.maxSize= 5000*1024^2)

print('identifyOverExpressedGenes')
cellchat <- identifyOverExpressedGenes(cellchat, do.fast = FALSE)
print('identifyOverExpressedInteractions')
cellchat <- identifyOverExpressedInteractions(cellchat)
print('projectData')

#### Infer cell-cell communication ####
print('computeCommunProb')
cellchat <- computeCommunProb(cellchat, type = "triMean")
print('filterCommunication')
cellchat <- filterCommunication(cellchat, min.cells = 10)

#### Infer cell-cell communication at pathway level ####
print('computeCommunProbPathway')
cellchat <- computeCommunProbPathway(cellchat)
print('aggregateNet')
cellchat <- aggregateNet(cellchat)

# save(cellchat, file = "scRNA/CellChat/cellChat_object_Responder.RData")
save(cellchat, file = "scRNA/CellChat/cellChat_object_NonResponder.RData")

####################
#### RESPONDERS ####
####################
load("scRNA/CellChat/cellChat_object_Responder.RData")
df.net <- subsetCommunication(cellchat)
# write.table(df.net, file = "scRNA/CellChat/output_cellChat_Responder.txt", row.names = FALSE, quote = FALSE)

#### CellChat plots
groupSize <- as.numeric(table(cellchat@idents))
colors_cellchat = c(colors_tme, colori_mp)
colors_cellchat = colors_cellchat[levels(cellchat@idents)]
par(mfrow = c(1,2), xpd = TRUE)
options(repr.plot.width = 40, repr.plot.height = 40)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge = F, title.name = "Number of interactions", vertex.label.cex = .8, color.use = colors_cellchat)
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength", vertex.label.cex = .8, color.use = colors_cellchat)

mat <- cellchat@net$weight
mp = names(colori_mp)
par(mfrow = c(2,4), xpd=TRUE, cex = 1, mar = c(2, 2, 2, 2))
for (MP in mp) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[MP, ] <- mat[MP, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = MP, edge.label.cex = 3,
                   vertex.label.cex = .7,
                   color.use = colors_cellchat, arrow.size = 0.05, arrow.width = 0.05, remove.isolate = TRUE)
}

tme = setdiff(names(colors_cellchat), mp)
par(mfrow = c(3,5), xpd=TRUE)
for (cell in tme) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[cell, ] <- mat[cell, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = cell, vertex.label.cex = .8, color.use = colors_cellchat, arrow.size = 0.05, arrow.width = 0.05)
}

#### Plot MHC-I and MHC-II in R samples
sort(cellchat@netP$pathways)
levels(cellchat@idents)

netVisual_aggregate(cellchat, signaling = 'MHC-I', arrow.size = 0.5,
                    layout = "circle", vertex.label.cex = 0.7, 
                    color.use = colors_cellchat, remove.isolate = TRUE
                    )
netVisual_aggregate(cellchat, signaling = 'MHC-II', arrow.size = 0.5,
                    layout = "circle", vertex.label.cex = 0.7,
                    color.use = colors_cellchat, remove.isolate = TRUE
                    )

########################
#### NON RESPONDERS ####
########################
load("scRNA/CellChat/cellChat_object_NonResponder.RData")
# df.net <- subsetCommunication(cellchat)
# write.table(df.net, file = "scRNA/CellChat/output_cellChat_NonResponder.txt", row.names = FALSE, quote = FALSE)

#### CellChat plots
groupSize <- as.numeric(table(cellchat@idents))
colors_cellchat = c(colors_tme, colori_mp)
colors_cellchat = colors_cellchat[levels(cellchat@idents)]
par(mfrow = c(1,2), xpd = TRUE)
options(repr.plot.width = 40, repr.plot.height = 40)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge = F, title.name = "Number of interactions", vertex.label.cex = .8, color.use = colors_cellchat)
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength", vertex.label.cex = .8, color.use = colors_cellchat)

par(mfrow = c(1, 1), xpd = TRUE)
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle", color.use = colors_cellchat)

# Contribution of each LR pair
netAnalysis_contribution(cellchat, signaling = 'PTN')

# Extract ligand-receptor interaction
pairLR <- extractEnrichedLR(cellchat, signaling = 'PTN', geneLR.return = FALSE)
LR.show <- pairLR$interaction_name # show one ligand-receptor pair

par(mfrow = c(2,2), xpd = TRUE, cex = 1)
options(repr.plot.width = 20, repr.plot.height = 20)
netVisual_individual(cellchat, signaling = pathways.show,  pairLR.use = LR.show[c(1,3)], 
                     color.use = colors_cellchat, remove.isolate = TRUE, 
                     arrow.size = 0.5, nCol = 1)