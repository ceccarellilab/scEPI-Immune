###################################
#### NMF on individual samples ####
###################################
library(Seurat)
library(NMF)
library(doMC)
library(parallel)
setwd("/home3/ciervo/scMULTIOME/Analisi/")
dir.create("scRNA/NMF")

load("scRNA/RData/merged_object.RData")
seu = subset(seu, Cell_annotation %in% 'Malignant')
dim(seu) # 31,093 x 53,702

sort(table(seu$orig.ident))
table(seu$orig.ident)
# min 277 max 8,523
# Pat03_W12 discarded for the low number of tumor cells

# Save gene expression matrices
seu.list = SplitObject(seu, split.by = 'orig.ident')
seu.list$Pat03_W12 = NULL

dir.create('scRNA/NMF/Matrices')
mclapply(seu.list, FUN = function(x){
  sample = unique(x$orig.ident)
  
  counts <- as.matrix(x@assays$RNA$counts)
  counts <- rowMeans(counts)
  counts <- sort(counts, decreasing = T)
  genes.filter <- names(counts)[1:7000]
  seu = subset(x, features = genes.filter)
  seu = NormalizeData(seu, verbose = F)
  mat = as.matrix(seu@assays$RNA$data)
  mat = mat - rowMeans(mat)
  
  mat[mat < 0] = 0
  tmp = rowSums(mat)
  sum(tmp == 0)
  mat = mat[tmp > 0, ]

  print(dim(mat))
  save(mat, file = paste0("scRNA/NMF/Matrices/", sample, ".RData"))
}, mc.cores = length(seu.list)
)

# NMF ----
dir.create("scRNA/NMF/NMF_list")
registerDoMC(4)

samples = list.files("scRNA/NMF/Matrices/", full.names = F)
samples = gsub(".RData", "", samples)
f
oreach(s = 1:length(samples))  %dopar% {
  message(samples[s])
  load(paste0("scRNA/NMF/Matrices/", samples[s], ".RData"))
  res.list = mclapply(4:9, function(r) nmf(mat, rank = r, nrun = 10, seed = 123, .opt = "vp5") , 
                      mc.cores = 6)
  names(res.list) = paste0("Rank_", 4:9)
  save(res.list, file = paste0("scRNA/NMF/NMF_list/Rank_list_10runs_", samples[s], ".RData"))
}

# Saving H and W matrices ----
dir.create("scRNA/NMF/Factors")
dir.create("scRNA/NMF/Factors_H")

# H matries
mat = NULL
lf = list.files("scRNA/NMF/NMF_list/", pattern = "10runs")
for(file in lf){
  load(paste0("scRNA/NMF/NMF_list/", file))
  sample = paste0(strsplit(gsub(".RData", "", file), split = "_")[[1]][4], '_', strsplit(gsub(".RData", "", file), split = "_")[[1]][5])
  for(i in 1:length(res.list)){
    rank = names(res.list)[i]
    H = res.list[[i]]@fit@H
    rownames(H) = paste0(sample, "_", rank, ".", 1:nrow(H))
    mat = rbind(mat, H)
  }
  save(mat, file = paste0("scRNA/NMF/Factors_H/", sample, "_factors.RData"))
  mat = NULL
}

# W matrices
mat = NULL
lf = list.files("scRNA/NMF/NMF_list/", pattern = "10runs")
for(file in lf){
  load(paste0("scRNA/NMF/NMF_list/", file))
  sample = paste0(strsplit(gsub(".RData", "", file), split = "_")[[1]][4], '_', strsplit(gsub(".RData", "", file), split = "_")[[1]][5])
  for(i in 1:length(res.list)){
    rank = names(res.list)[i]
    W = res.list[[i]]@fit@W
    colnames(W) = paste0(sample, "_", rank, ".", 1:ncol(W))
    mat = cbind(mat, W)
  }
  save(mat, file = paste0("scRNA/NMF/Factors/", sample, "_factors.RData"))
  mat = NULL
}


