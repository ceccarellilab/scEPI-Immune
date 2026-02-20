###################
#### Libraries ####
###################
library(Seurat)
library(SeuratData)
library(dplyr)
library(ggplot2)
library(clusterProfiler)
library(enrichplot)

setwd('/home3/ciervo/scMULTIOME/Analisi/')
source('colori_finali.R')

llist = list('GO' = read.gmt('MSigDb_11Feb2025/c5.go.v2024.1.Hs.symbols.gmt'),
             'REAC' = read.gmt('MSigDb_11Feb2025/c2.cp.reactome.v2024.1.Hs.symbols.gmt'),
             'WIKIPATHWAYS' = read.gmt('MSigDb_11Feb2025/c2.cp.wikipathways.v2024.1.Hs.symbols.gmt'),
             'KEGG' = read.gmt('MSigDb_11Feb2025/c2.cp.kegg_legacy.v2024.1.Hs.symbols.gmt'),
             'HALLMARKS' = read.gmt('MSigDb_11Feb2025/h.all.v2024.1.Hs.symbols.gmt')
             
)

signatures <- Reduce(rbind, llist)
signatures <- signatures[-grep("GOCC|GOMF", signatures$term), ]
rm(llist)

GSEA_list = vector('list', 7)
names(GSEA_list) = paste0('MP_', 1:7)

set.seed(1234)
for(i in 1:7){
  ans <- as.data.frame(readr::read_csv(paste0("scRNA/memento/MP_", i, ".csv")))
  ranked <- ans$de_coef
  names(ranked) <- ans$gene
  ranked <- sort(ranked, decreasing = TRUE)
  head(ranked)
  res <- clusterProfiler::GSEA(geneList = ranked, 
                               minGSSize = 15, 
                               maxGSSize = 250, eps = 0,
                               pvalueCutoff = 0.05, 
                               TERM2GENE = signatures)
  pdf(file = paste0('scRNA/memento/GSEA_MP_', i, '.pdf'), width = 10, height = 10)
  plot(enrichplot::dotplot(res, x = 'NES', split = '.sign'))
  dev.off()
  GSEA_list[[i]] = res
}
save(GSEA_list, file = 'scRNA/memento/result_GSEA.RData')
