############################################################
########## Metaprogram clustering as described by ##########
### Gavish et al, 2023 - DOI: 10.1038/s41586-023-06130-4 ###
############################################################
library(reshape2)
library(NMF)
library(ggplot2)
library(RColorBrewer)
library(viridis)
library(ComplexHeatmap)
library(circlize)
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
library(gprofiler2)
custom_magma <- c(colorRampPalette(c("white", rev(magma(323, begin = 0.15))[1]))(10), rev(magma(323, begin = 0.18)))

setwd("/home3/ciervo/scMULTIOME/Analisi/")
source("/home/ciervo/HeadAndNeck/Analisi_03aprile24/NMF/Robust_nmf.R")   
source('colori_finali.R')

## Input:
# Genes_nmf_w_basis is a list in which each entry contains NMF gene-scores of a single sample. In our study we ran NMF using ranks 4-9 on the top 7000 genes in each sample. Hence each entry in Genes_nmf_w_basis is a matrix with 7000 rows (genes) X 39 columns (NMF programs)  
# For the code below to run smoothly, please use the following nomenclature:
# 1) End entry names in Genes_nmf_w_basis (i.e. each sample name) with '_rank4_9_nruns10.RDS' 
# 2) End each matrix column name with an extension that represents the NMF rank and program index. for example '_rank4_9_nrun10.RDS.4.1' to represent the first NMF program obtained using rank=4, or '_rank4_9_nrun10.RDS.6.5' to represent the fifth NMF program obtained using rank=6    
# See Genes_nmf_w_basis_example.RDS for an example 
# We define MPs in 2 steps:
# 1) The function robust_nmf_programs.R performs filtering, so that programs selected for defining MPs are:
#    i) Robust                - recur in more that one rank within the sample 
#    ii) Non-redundant        - once a NMF program is selected, other programs within the sample that are similar to it are removed
#    iii) Not sample specific - has similarity to NMF programs in other samples 
# ** Please see https://github.com/gabrielakinker/CCLE_heterogeneity for more details on how to define robust NMF programs 
# 2) Selected NMFs are then clustered iteratively, as described in Figure S1. At the end of the process, each cluster generates a list of the 50 genes (i.e. the MP) that represent the NMF programs that contributed to the cluster. Notably, not all initially selected NMFs end up participating in a cluster  

# ---------------------------------------------------------------------------------------------------- #
# Select NMF programs ----
# ----------------------------------------------------------------------------------------------------
## Parameters 
ngenes = 50
intra_min_parameter <- floor(ngenes * 0.7)
intra_max_parameter <- floor(ngenes * 0.2)
inter_min_parameter <- floor(ngenes * 0.2)

lf = list.files("scRNA/NMF/Factors/", all.files = FALSE, recursive = FALSE, 
                full.names = FALSE, pattern = "_factors")
llist = list()
for(file in lf){
  load(paste0("scRNA/NMF/Factors/", file))
  sample = gsub("_factors.RData", "", file)
  llist[[sample]] = mat
  rm(mat)
}

# Get top 50 genes for each NMF program 
nmf_programs <- lapply(llist, function(x) apply(x, 2, function(y) names(sort(y, decreasing = TRUE))[1:ngenes]))
nmf_programs <- lapply(nmf_programs, toupper) ## convert all genes to uppercase 

tmp <- lapply(nmf_programs, function(x) apply(x, 2, function(y) sum(grepl(pattern = "MT-|RPL|RPS", x = y)))) # discard programs with high RPL-S/MT- genes
tmp = lapply(tmp, function(x) names(x)[x <= ngenes * 0.2])
tmp = unlist(lapply(tmp, function(x) Reduce(c, x)))
nmf_programs = lapply(nmf_programs, function(x) x[, colnames(x) %in% tmp])

# For each sample, select robust NMF programs (i.e. observed using different ranks in the same sample), remove redundancy due to multiple ranks, and apply a filter based on the similarity to programs from other samples. 
nmf_filter_ccle <- robust_nmf_programs(nmf_programs, intra_min = intra_min_parameter, intra_max = intra_max_parameter, 
                                       inter_filter = FALSE, inter_min = inter_min_parameter)
nmf_programs <- lapply(nmf_programs, function(x) x[, is.element(colnames(x), nmf_filter_ccle), drop = FALSE])
nmf_programs <- do.call(cbind, nmf_programs)
dim(nmf_programs)
# 75 programs

# calculate similarity between programs
nmf_intersect <- apply(nmf_programs, 2, function(x) apply(nmf_programs , 2, function(y) length(intersect(x, y)))) 

# hierarchical clustering of the similarity matrix 
nmf_intersect_hc <- hclust(as.dist(ngenes-nmf_intersect), method = "average")
nmf_intersect_hc <- reorder(as.dendrogram(nmf_intersect_hc), colMeans(nmf_intersect))
nmf_intersect <- nmf_intersect[order.dendrogram(nmf_intersect_hc), order.dendrogram(nmf_intersect_hc)]

# ----------------------------------------------------------------------------------------------------
# Cluster selected NMF programs to generate MPs
# ----------------------------------------------------------------------------------------------------
### Parameters for clustering
Min_intersect_initial <- 10    # the minimal intersection cutoff for defining the first NMF program in a cluster
Min_intersect_cluster <- 10    # the minimal intersection cutoff for adding a new NMF to the forming cluster
Min_group_size        <- 2   # the minimal group size to consider for defining the first NMF program in a cluster

Sorted_intersection       <-  sort(apply(nmf_intersect , 2, function(x) (length(which(x >= Min_intersect_initial))-1)), decreasing = TRUE)

Cluster_list              <- list()   ### Every entry contains the NMFs of a chosen cluster
MP_list                   <- list()
k                         <- 1
Curr_cluster              <- c()

nmf_intersect_original    <- nmf_intersect

while (Sorted_intersection[1] >= Min_group_size) {  
  
  Curr_cluster <- c(Curr_cluster , names(Sorted_intersection[1]))
  
  ### intersection between all remaining NMFs and Genes in MP 
  Genes_MP <- nmf_programs[,names(Sorted_intersection[1])] # Genes in the forming MP are first chosen to be those in the first NMF. Genes_MP always has only 50 genes and evolves during the formation of the cluster
  nmf_programs <- nmf_programs[,-match(names(Sorted_intersection[1]) , colnames(nmf_programs))]  # remove selected NMF
  Intersection_with_Genes_MP <- sort(apply(nmf_programs, 2, function(x) length(intersect(Genes_MP,x))) , decreasing = TRUE) # intersection between all other NMFs and Genes_MP  
  NMF_history <- Genes_MP  # has genes in all NMFs in the current cluster, for redefining Genes_MP after adding a new NMF 
  
  ### Create gene list is composed of intersecting genes (in descending order by frequency). When the number of genes with a given frequency span bewond the 50th genes, they are sorted according to their NMF score.    
  while ( Intersection_with_Genes_MP[1] >= Min_intersect_cluster) {  
    
    Curr_cluster  <- c(Curr_cluster , names(Intersection_with_Genes_MP)[1])
    
    Genes_MP_temp   <- sort(table(c(NMF_history , nmf_programs[,names(Intersection_with_Genes_MP)[1]])), decreasing = TRUE)   ## Genes_MP is newly defined each time according to all NMFs in the current cluster 
    Genes_at_border <- Genes_MP_temp[which(Genes_MP_temp == Genes_MP_temp[ngenes])]   ### genes with overlap equal to the 50th gene
    
    if (length(Genes_at_border)>1){
      ### Sort last genes in Genes_at_border according to maximal NMF gene scores
      ### Run across all NMF programs in Curr_cluster and extract NMF scores for each gene
      Genes_curr_NMF_score <- c()
      for (i in Curr_cluster) {
        curr_study <- paste0(strsplit(i , "_")[[1]][[1]], "_", strsplit(i , "_")[[1]][[2]])
        Q <- llist[[curr_study]][match(names(Genes_at_border),toupper(rownames(llist[[curr_study]])))[!is.na(match(names(Genes_at_border),toupper(rownames(llist[[curr_study]]))))], i] 
        names(Q) <- names(Genes_at_border[!is.na(match(names(Genes_at_border),toupper(rownames(llist[[curr_study]]))))])  ### sometimes when adding genes the names do not appear 
        Genes_curr_NMF_score <- c(Genes_curr_NMF_score,  Q )
      }
      Genes_curr_NMF_score_sort <- sort(Genes_curr_NMF_score , decreasing = TRUE)
      Genes_curr_NMF_score_sort <- Genes_curr_NMF_score_sort[unique(names(Genes_curr_NMF_score_sort))]   
      
      Genes_MP_temp <- c(names(Genes_MP_temp[which(Genes_MP_temp > Genes_MP_temp[ngenes])]) , names(Genes_curr_NMF_score_sort))
      
    } else {
      Genes_MP_temp <- names(Genes_MP_temp)[1:ngenes]
      
    }
    
    NMF_history <- c(NMF_history , nmf_programs[,names(Intersection_with_Genes_MP)[1]]) 
    Genes_MP <- Genes_MP_temp[1:ngenes]
    
    nmf_programs <- nmf_programs[,-match(names(Intersection_with_Genes_MP)[1] , colnames(nmf_programs))]  # remove selected NMF
    
    Intersection_with_Genes_MP <- sort(apply(nmf_programs, 2, function(x) length(intersect(Genes_MP,x))) , decreasing = TRUE) # intersection between all other NMFs and Genes_MP  
    
  }
  
  Cluster_list[[paste0("Cluster_",k)]] <- Curr_cluster
  MP_list[[paste0("MP_",k)]] <- Genes_MP
  k <- k+1
  
  nmf_intersect <- nmf_intersect[-match(Curr_cluster,rownames(nmf_intersect) ) , -match(Curr_cluster,colnames(nmf_intersect) ) ]  # Remove current chosen cluster
  
  Sorted_intersection <-  sort(apply(nmf_intersect , 2, function(x) (length(which(x>=Min_intersect_initial))-1)  ) , decreasing = TRUE)   # Sort intersection of remaining NMFs not included in any of the previous clusters
  
  Curr_cluster <- c()
  print(dim(nmf_intersect)[2])
}

####  Sort Jaccard similarity plot according to new clusters:
inds_sorted <- c()

for (j in 1:length(Cluster_list)){
  inds_sorted <- c(inds_sorted , match(Cluster_list[[j]] , colnames(nmf_intersect_original)))
}

inds_new <- c(inds_sorted, which(is.na( match(1:dim(nmf_intersect_original)[2],inds_sorted)))) ### clustered NMFs will appear first, and the latter are the NMFs that were not clustered
df = nmf_intersect_original[inds_new, inds_new]
df = df[nrow(df):1, ]
ccol = colorRampPalette(custom_magma)

### Heatmap ####
Heatmap(df, name = "Jaccard\nsimilarity",
        cluster_columns = F, 
        show_row_names = F,
        cluster_column_slices = F,
        cluster_row_slices = F,
        cluster_rows = F, 
        show_column_names = T,
        col = ccol(75),
        width = unit(20, "cm"),
        height = unit(20, "cm")
)

save(Cluster_list, file = 'scRNA/NMF/Program_cluster_list.RData')

#### Adding row and column track
annot_df = data.frame("Row_programs" = rownames(df),
                      "Col_programs" = colnames(df))

annot_df$MP = ""
annot_df$MP[annot_df$Row_programs %in% Cluster_list$Cluster_1] = "MP 1" 
annot_df$MP[annot_df$Row_programs %in% Cluster_list$Cluster_2] = "MP 2" 
annot_df$MP[annot_df$Row_programs %in% Cluster_list$Cluster_3] = "MP 3" 
annot_df$MP[annot_df$Row_programs %in% Cluster_list$Cluster_4] = "MP 4" 
annot_df$MP[annot_df$Row_programs %in% Cluster_list$Cluster_5] = "MP 5" 
annot_df$MP[annot_df$Row_programs %in% Cluster_list$Cluster_6] = "MP 6"
annot_df$MP[annot_df$Row_programs %in% Cluster_list$Cluster_7] = "MP 7"

annot_df$MP_cols = "white"
annot_df$MP_cols[annot_df$MP %in% "MP 1"] = colori_mp[1]
annot_df$MP_cols[annot_df$MP %in% "MP 2"] = colori_mp[2]
annot_df$MP_cols[annot_df$MP %in% "MP 3"] = colori_mp[3]
annot_df$MP_cols[annot_df$MP %in% "MP 4"] = colori_mp[4]
annot_df$MP_cols[annot_df$MP %in% "MP 5"] = colori_mp[5]
annot_df$MP_cols[annot_df$MP %in% "MP 6"] = colori_mp[6]
annot_df$MP_cols[annot_df$MP %in% "MP 7"] = colori_mp[7]

annot_df = annot_df[annot_df$MP %in% paste0('MP ', 1:7), ]
df = df[annot_df$Row_programs, annot_df$Row_programs]

mp <- annot_df$MP
names(mp) <- annot_df$Row_programs
mp_col <- unique(annot_df$MP_cols)
names(mp_col) <- unique(annot_df$MP)

row_ha = rowAnnotation(Metaprograms = mp,
                       col = list(Metaprograms = mp_col
                       ), na_col = "white", show_annotation_name = FALSE,
                       annotation_name_gp = gpar(fontsize = 10))
# column_ha 
# View(annot_df)

annot_df$Patient = unlist(lapply(1:length(annot_df$Col_programs), FUN = function(x) strsplit(annot_df$Col_programs[x], '_')[[1]][[1]]))
table(annot_df$Patient)
annot_df$Patient_cols = "white"
annot_df$Patient_cols[annot_df$Patient %in% "Pat02"] = "#F8766D"
annot_df$Patient_cols[annot_df$Patient %in% "Pat03"] = "#A3A500"
annot_df$Patient_cols[annot_df$Patient %in% "Pat06"] = "#00BF7D"
annot_df$Patient_cols[annot_df$Patient %in% "Pat14"] = "#00B0F6"
annot_df$Patient_cols[annot_df$Patient %in% "Pat15"] = "#E76BF3"

pat <- annot_df$Patient
names(pat) <- annot_df$Col_programs
pat_col <- unique(annot_df$Patient_cols)
names(pat_col) <- unique(annot_df$Patient)

column_ha = columnAnnotation(Patient = pat, 
                             col = list(Patient = pat_col), 
                             na_col = "white", show_annotation_name = FALSE,
                             annotation_name_gp = gpar(fontsize = 10))

df = df[, rev(colnames(df))]
ccol = colorRampPalette(custom_magma)
h1 = Heatmap(df, name = "Jaccard\nsimilarity",
             # top_annotation = column_ha,
             left_annotation = row_ha, 
             cluster_columns = F, 
             show_row_names = F,
             cluster_column_slices = F,
             cluster_row_slices = F,
             cluster_rows = F, 
             show_column_names = F,
             col = ccol(75),
             width = unit(15, "cm"),
             height = unit(15, "cm")
)
h1 

pdf(file = 'scRNA/NMF_heatmap.pdf', width = 10, height = 10)
h1
dev.off()

save(df, MP_list, file = 'scRNA/NMF/metaprograms.RData')
write.table(x = MP_list, file = "scRNA/NMF/Metaprograms.txt", sep = "\t", row.names = F, col.names = T)

# Enrichment MP programs 50 genes ----
# Load MP 
load("scRNA/NMF/metaprograms.RData")
MP = as.data.frame(MP_list)

llist = list('GO' = read.gmt('MSigDb_11Feb2025/c5.go.v2024.1.Hs.symbols.gmt'),
             'REAC' = read.gmt('MSigDb_11Feb2025/c2.cp.reactome.v2024.1.Hs.symbols.gmt'),
             'WIKIPATHWAYS' = read.gmt('MSigDb_11Feb2025/c2.cp.wikipathways.v2024.1.Hs.symbols.gmt'),
             'KEGG' = read.gmt('MSigDb_11Feb2025/c2.cp.kegg_legacy.v2024.1.Hs.symbols.gmt'),
             'HALLMARKS' = read.gmt('MSigDb_11Feb2025/h.all.v2024.1.Hs.symbols.gmt')
             
)

signatures <- Reduce(rbind, llist)
signatures <- signatures[-grep("GOCC|GOMF", signatures$term), ]
rm(llist)

MP_top50_enrichment = lapply(1:7, function(x){
  res = enricher(gene = MP[, x],
                 TERM2GENE = signatures,
                 minGSSize = 15, 
                 maxGSSize = 250,
                 pvalueCutoff = 1, 
                 pAdjustMethod = "fdr"
  )
})
names(MP_top50_enrichment) = paste0('MP', 1:7)
save(MP_top50_enrichment, file = 'scRNA/NMF/Enrichment_top50.RData')