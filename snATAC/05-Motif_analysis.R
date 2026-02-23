setwd('/home3/ciervo/scMULTIOME/Analisi/')
library(Seurat)
library(Signac)
library(dplyr)
library(ggplot2)
library(SingleCellExperiment)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomicRanges)
library(parallel)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(motifmatchr)
library(JASPAR2024)
library(TFBSTools)
library(chromVAR)
library(ggseqlogo)
library(motifbreakR)
library(universalmotif)
source('colori_finali.R')

# Add motif to the object ----
# Get a list of motif position frequency matrices from the JASPAR database
# jaspar <- JASPAR2024()
# sq24 <- RSQLite::dbConnect(RSQLite::SQLite(), db(jaspar))
# motifs24 <- TFBSTools::getMatrixSet(sq24, list(species = "Homo sapiens", collection = "CORE"))

# Hocomoco v13 motifs
load('scATAC/RData/malignant_subset.RData')
DefaultAssay(malignant)

malignant = subset(malignant, Metaprogram_assignment %in% 'non-classified', invert = TRUE)
DefaultAssay(malignant)

motifs <- universalmotif::read_meme(file = 'scATAC/H13CORE_meme_format.meme')
head(motifs, 5)
names(motifs) = lapply(1:length(motifs), function(x) motifs[[x]]@name)
motifs = lapply(motifs, function(x) universalmotif::convert_motifs(x, class = 'TFBSTools-PFMatrix'))
motifs.list <- do.call(PFMatrixList, motifs)
rm(motifs)

# Scan the DNA sequence of each peak for the presence of each motif
malignant <- AddMotifs(malignant, genome = BSgenome.Hsapiens.UCSC.hg38, pfm = motifs.list)
# save(malignant, file = 'scATAC/RData/malignant_subset_motifs.RData')

# Transcription factors ----
load('scATAC/RData/malignant_subset_motifs.RData')
Idents(malignant) = 'Metaprogram_assignment'
DefaultAssay(malignant) = 'ATAC'

DA.peaks = mclapply(1:7, FUN = function(x){
  metaprogram = paste0('MP_', x)
  future::plan('sequential')
  top.da.peak <- FindMarkers(malignant, ident.1 = metaprogram)
  top.da.peak = rownames(top.da.peak[top.da.peak$p_val_adj < 0.05 & top.da.peak$avg_log2FC >= 0.58, ])
  return(top.da.peak)
}, mc.cores = 7)
names(DA.peaks) = names(colori_mp)

set.seed(12345)
TF.of.interest = mclapply(1:7, FUN = function(x){
  metaprogram = names(colori_mp)[x]
  top.da.peak = DA.peaks[[metaprogram]]
  open.peaks <- AccessiblePeaks(malignant, idents = metaprogram)
  
  # match the overall GC content in the peak set
  meta.feature <- GetAssayData(malignant, assay = "ATAC", layer = "meta.features")
  peaks.matched <- MatchRegionStats(
    meta.feature = meta.feature[open.peaks, ],
    query.feature = meta.feature[top.da.peak, ],
    n = 50000
  )
  
  # test enrichment (hypergeometric test)
  enriched.motifs <- FindMotifs(
    object = malignant,
    features = top.da.peak, 
    background = peaks.matched
  )
  
  # Annotation dei geni in cui cade il motif
  motif.use <- enriched.motifs[enriched.motifs$fold.enrichment >= 1.3, ] # Motivo di interesse
  motif.use$TF.name <- sapply(motif.use$motif.name, function(x) strsplit(x, '\\.')[[1]][[1]])
  motif.use = motif.use[order(motif.use$fold.enrichment, decreasing = TRUE), ]
  motif.use = motif.use[!duplicated(motif.use$TF.name), ]
  motif.use$MP = metaprogram
  motif.use = motif.use[c(1:30, which(motif.use$TF.name %in% 'MITF')), ]
}, mc.cores = 7)
names(TF.of.interest) = names(colori_mp)
TF.of.interest[['MP_5']] = na.omit(TF.of.interest[['MP_5']])
save(TF.of.interest, file = 'scATAC/TFOI_MPs.RData')
load('scATAC/TFOI_MPs.RData')


tmp = lapply(TF.of.interest, function(x) x$TF.name)
tmp = Reduce(c, tmp) %>% unique()

# Selected TF
tf_vector <- c(
  "MITF", "KMT2B", "EGR1", "EGR2", "HINFP", "E2F4", "TFDP1", 
  "MYCN", "GCM2", "E2F3", "E2F1", "KLF3", "KLF16", "NRF1", 
  "CUX2", "CUX1", "MBD1", "JUND", "NFE2", "BACH2", "FOSL2", 
  "JUN", "JUNB", "ATF3", "FOSL1", "FOSB", "FOS", "BATF", 
  "BATF3", "BACH1", "NF2L2", "NF2L1", "XBP1", "MAFG", "ATF4", 
  "ZN441", "NFAC1", "NFAT5", "NFAC3", "NFAC4", "SP1", "NFAC2", 
  "IRF7", "IRF5", "IRF2", "IRF8", "IRF4", "IRF9", "IRF6", "PRDM1", 
  "IRF3", "IRF1", "STAT1", "STAT2", "NFKB1", "NFKB2", "RELB", 
  "TF65", "MAX", "MLX", "ARNT2", "MLXPL", "NPAS2", "MNT", "LEF1", 
  "TF7L2", "TF7L1", "TCF7", "GLIS3"
)
setdiff(tmp, tf_vector)
# tmp = tmp[tmp %in% f_vector]

enrichment_mat = matrix(data = NA, nrow = length(tf_vector), ncol = 7)
colnames(enrichment_mat) = names(TF.of.interest)
rownames(enrichment_mat) = tf_vector

for(tf in rownames(enrichment_mat)){
  for(mp in colnames(enrichment_mat)){
    if(tf %in% TF.of.interest[[mp]][, 'TF.name']){
      enrichment_mat[tf, mp] = TF.of.interest[[mp]][ TF.of.interest[[mp]]$TF.name %in% tf, 'fold.enrichment']
    }else{
      enrichment_mat[tf, mp] = 0
    }
  }
}

# enrichment_mat = enrichment_mat[order(enrichment_mat[, 1], decreasing = TRUE), ]
head(enrichment_mat)
range(enrichment_mat) # 0 4
# quantile(enrichment_mat, c(.95))
enrichment_mat[enrichment_mat > 3] = 3

tf_family = readxl::read_xlsx('scATAC/TF_family.xlsx')
table(tf_family$`TF family`) == 1
tf_family$`TF family`[tf_family$`TF family` %in% c('GCM', 'MBD', 'Nrf1', 'zf-CXXC')] = 'Others'
families_colors = pals::brewer.paired(length(unique(tf_family$`TF family`)))
names(families_colors) = sort(unique(tf_family$`TF family`))
tf_family$TF.family.col = factor(tf_family$`TF family`, levels = names(families_colors), labels = families_colors)

save(enrichment_mat, tf_family, file = 'scATAC/Heatmap_motif.RData')

# Footprinting ----
load('scATAC/RData/malignant_subset_motifs.RData')
Idents(malignant) = 'Metaprogram_assignment'
table(Idents(malignant), exclude = NULL)
DefaultAssay(malignant) = 'ATAC'

load('scATAC/TFOI_MPs.RData')

# tf_vector <- c(
#   "MITF", "KMT2B", "EGR1", "EGR2", "HINFP", "E2F4", "TFDP1", 
#   "MYCN", "GCM2", "E2F3", "E2F1", "KLF3", "KLF16", "NRF1", 
#   "CUX2", "CUX1", "MBD1", "JUND", "NFE2", "BACH2", "FOSL2", 
#   "JUN", "JUNB", "ATF3", "FOSL1", "FOSB", "FOS", "BATF", 
#   "BATF3", "BACH1", "NF2L2", "NF2L1", "XBP1", "MAFG", "ATF4", 
#   "ZN441", "NFAC1", "NFAT5", "NFAC3", "NFAC4", "SP1", "NFAC2", 
#   "IRF7", "IRF5", "IRF2", "IRF8", "IRF4", "IRF9", "IRF6", "PRDM1", 
#   "IRF3", "IRF1", "STAT1", "STAT2", "NFKB1", "NFKB2", "RELB", 
#   "TF65", "MAX", "MLX", "ARNT2", "MLXPL", "NPAS2", "MNT", "LEF1", 
#   "TF7L2", "TF7L1", "TCF7", "GLIS3"
# )

motif.name = Reduce(rbind, TF.of.interest)
motif.name = motif.name[motif.name$TF.name %in% tf_vector, ]
motif.name = motif.name$motif %>% unique

tmp = Motifs(malignant)
tmp@motif.names[grep('^ATF', tmp@motif.names)]

motif.name = c(motif.name, c('MYC.H13CORE.0.P.B', 'AP2A.H13CORE.0.PSM.A', 'PPARG.H13CORE.0.P.B', 'NFAC2.H13CORE.0.P.B'))

malignant <- Footprint(
  object = malignant,
  motif.name = c('EGR1.H13CORE.0.PS.A', 'ATF3.H13CORE.0.P.B'),
  genome = BSgenome.Hsapiens.UCSC.hg38, 
  in.peaks = TRUE
)
save(malignant, file = 'scATAC/RData/malignant_subset_motifs.RData')

source('/home/ciervo/Functions/scATAC_function.R')
PlotFootprint_mod(malignant, features = c('IRF1.H13CORE.0.P.B'), 
                  group.by = 'Metaprogram_assignment', 
                  group.colors = colori_mp)

# Annotation with ENCODE regions pELS, dELS, PLS ----
# ENCODE annotation  
# H13CORE_motifs_annotation <- readr::read_delim("scATAC/H13CORE_motifs_annotation.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
load('scATAC/RData/malignant_subset_motifs.RData')
# load('/home/noviello/TE_utilities/Genome_class_hg38_extended.RDa')

# cCRE encode
load('/home/noviello/TE_utilities/cCRE_ENCODE.RDa', verbose = TRUE)
CCRE_links = vroom::vroom("/home/noviello/TE_utilities/ENCODE/V4-hg38.Gene-Links.eQTLs.txt",col_names = F)
head(cCRE_ENCODE)
CCRE_links[1:5, ]
length(intersect(CCRE_links$X1, cCRE_ENCODE$X4)) # 0
length(intersect(CCRE_links$X1, cCRE_ENCODE$X5)) # 765816
colnames(mcols(cCRE_ENCODE))[2] = 'cCRE_link'
colnames(CCRE_links)[1] = 'cCRE_link'

signatures <- list('GO' = read.gmt(gmtfile = "MSigDb_11Feb2025/c5.go.v2024.1.Hs.symbols.gmt"),
                   'KEGG' = read.gmt(gmtfile = "MSigDb_11Feb2025/c2.cp.kegg_legacy.v2024.1.Hs.symbols.gmt"),
                   'REAC' = read.gmt(gmtfile = "MSigDb_11Feb2025/c2.cp.reactome.v2024.1.Hs.symbols.gmt"),
                   'HALL' = read.gmt(gmtfile = "MSigDb_11Feb2025/h.all.v2024.1.Hs.symbols.gmt")
)
signatures$GO = signatures$GO[grep('^GOBP', signatures$GO$term), ]
signatures = Reduce(rbind, signatures)

Idents(malignant) = 'Metaprogram_assignment'
DefaultAssay(malignant) = 'ATAC'

DA.peaks = mclapply(1:7, FUN = function(x){
  metaprogram = names(colori_mp)[x]
  future::plan('sequential')
  top.da.peak <- FindMarkers(malignant, ident.1 = metaprogram)
  top.da.peak = rownames(top.da.peak[top.da.peak$p_val_adj < 0.05 & top.da.peak$avg_log2FC >= 0.58, ])
  return(top.da.peak)
}, mc.cores = 5)
names(DA.peaks) = names(colori_mp)
mean(sapply(DA.peaks, length))

TF.of.interest = mclapply(1:7, FUN = function(x){
  metaprogram = names(colori_mp)[x]
  top.da.peak = DA.peaks[[metaprogram]]
  open.peaks <- AccessiblePeaks(malignant, idents = metaprogram)
  
  # match the overall GC content in the peak set
  meta.feature <- GetAssayData(malignant, assay = "ATAC", layer = "meta.features")
  peaks.matched <- MatchRegionStats(
    meta.feature = meta.feature[open.peaks, ],
    query.feature = meta.feature[top.da.peak, ],
    n = 50000
  )
  
  # test enrichment (hypergeometric test)
  enriched.motifs <- FindMotifs(
    object = malignant,
    features = top.da.peak, 
    background = peaks.matched
  )
  
  # Annotation dei geni in cui cade il motif
  motif.use <- enriched.motifs$motif[enriched.motifs$fold.enrichment >= 1.3] # Motivo di interesse
  return(motif.use)
}, mc.cores = 5)
names(TF.of.interest) = names(colori_mp)
tmp = lapply(TF.of.interest, function(x) sapply(x, function(y) strsplit(y, '\\.')[[1]][[1]]))
tmp = Reduce(c, tmp) %>% unique()

CCRE_links$cCRE_link

Genes_PLS = mclapply(1:7, FUN = function(x){
  metaprogram = names(colori_mp)[x]
  motif.use = TF.of.interest[[metaprogram]]
  
  # Annotation dei geni in cui cade il motif
  peaks.with.motif <- as.matrix(GetMotifData(object = malignant)[DA.peaks[[metaprogram]], motif.use]) # Matrix T/F dei picchi che hanno il motivo
  interest.peaks <- rownames(peaks.with.motif)[apply(peaks.with.motif, 1, function(x) any(x == TRUE))]
  
  regions = GRanges(seqnames = sapply(strsplit(interest.peaks,"-"), function(x) x[1]),
                    ranges = IRanges(start = as.numeric(sapply(strsplit(interest.peaks,"-"), function(x) x[2])),
                                     end = as.numeric(sapply(strsplit(interest.peaks,"-"), function(x) x[3]))))
  
  hg38_PromEnh = cCRE_ENCODE[cCRE_ENCODE$X6 %in% c('PLS')]
  subset.hg38 = subsetByOverlaps(x = hg38_PromEnh, ranges = regions)
  
  subset.hg38 = merge(subset.hg38, CCRE_links, by = 'cCRE_link')
  subset.hg38 = subset.hg38[!is.na(subset.hg38$X3), ]
  subset.hg38$region = paste0(subset.hg38$seqnames, '_', subset.hg38$start, "_", subset.hg38$end)
  subset.hg38 = distinct(subset.hg38, region, X3, .keep_all = TRUE)
  subset.hg38 = subset.hg38[subset.hg38$X4.y %in% 'protein_coding', ]
  # subset.hg38 = subset.hg38[!is.na(subset.hg38$gene_name), ]
  
  # rm(meta.feature, peaks.matched, enriched.motifs, peaks.with.motif, interest.peaks, regions, hg38_PromEnh)
  
  # Enrichment
  DEGs <- as.data.frame(readr::read_csv(paste0("scRNA/memento/", metaprogram, ".csv")))
  DEGs$padjust = p.adjust(DEGs$de_pval)
  DEGs_genes = DEGs$gene[DEGs$de_coef >= log(2) & DEGs$padjust < 0.05]
  genes = intersect(unique(subset.hg38$X3), DEGs_genes)
  return(genes)
}, mc.cores = 7)
names(Genes_PLS) = names(colori_mp)

Genes_pELS = mclapply(1:7, FUN = function(x){
  metaprogram = names(colori_mp)[x]
  motif.use = TF.of.interest[[metaprogram]]
  
  # Annotation dei geni in cui cade il motif
  peaks.with.motif <- as.matrix(GetMotifData(object = malignant)[DA.peaks[[metaprogram]], motif.use]) # Matrix T/F dei picchi che hanno il motivo
  interest.peaks <- rownames(peaks.with.motif)[apply(peaks.with.motif, 1, function(x) any(x == TRUE))]
  
  regions = GRanges(seqnames = sapply(strsplit(interest.peaks,"-"), function(x) x[1]),
                    ranges = IRanges(start = as.numeric(sapply(strsplit(interest.peaks,"-"), function(x) x[2])),
                                     end = as.numeric(sapply(strsplit(interest.peaks,"-"), function(x) x[3]))))
  
  hg38_PromEnh = cCRE_ENCODE[cCRE_ENCODE$X6 %in% c('pELS')]
  subset.hg38 = subsetByOverlaps(x = hg38_PromEnh, ranges = regions)
  
  subset.hg38 = merge(subset.hg38, CCRE_links, by = 'cCRE_link')
  subset.hg38 = subset.hg38[!is.na(subset.hg38$X3), ]
  subset.hg38$region = paste0(subset.hg38$seqnames, '_', subset.hg38$start, "_", subset.hg38$end)
  subset.hg38 = distinct(subset.hg38, region, X3, .keep_all = TRUE)
  subset.hg38 = subset.hg38[subset.hg38$X4.y %in% 'protein_coding', ]
  # subset.hg38 = subset.hg38[!is.na(subset.hg38$gene_name), ]
  
  # rm(meta.feature, peaks.matched, enriched.motifs, peaks.with.motif, interest.peaks, regions, hg38_PromEnh)
  
  # Enrichment
  DEGs <- as.data.frame(readr::read_csv(paste0("scRNA/memento/", metaprogram, ".csv")))
  DEGs$padjust = p.adjust(DEGs$de_pval)
  DEGs_genes = DEGs$gene[DEGs$de_coef >= log(2) & DEGs$padjust < 0.05]
  genes = intersect(unique(subset.hg38$X3), DEGs_genes)
  return(genes)
}, mc.cores = 7)
names(Genes_pELS) = names(colori_mp)

Genes_dELS = mclapply(1:7, FUN = function(x){
  metaprogram = names(colori_mp)[x]
  motif.use = TF.of.interest[[metaprogram]]
  
  # Annotation dei geni in cui cade il motif
  peaks.with.motif <- as.matrix(GetMotifData(object = malignant)[DA.peaks[[metaprogram]], motif.use]) # Matrix T/F dei picchi che hanno il motivo
  interest.peaks <- rownames(peaks.with.motif)[apply(peaks.with.motif, 1, function(x) any(x == TRUE))]
  
  regions = GRanges(seqnames = sapply(strsplit(interest.peaks,"-"), function(x) x[1]),
                    ranges = IRanges(start = as.numeric(sapply(strsplit(interest.peaks,"-"), function(x) x[2])),
                                     end = as.numeric(sapply(strsplit(interest.peaks,"-"), function(x) x[3]))))
  
  hg38_PromEnh = cCRE_ENCODE[cCRE_ENCODE$X6 %in% c('dELS')]
  subset.hg38 = subsetByOverlaps(x = hg38_PromEnh, ranges = regions)
  
  subset.hg38 = merge(subset.hg38, CCRE_links, by = 'cCRE_link')
  subset.hg38 = subset.hg38[!is.na(subset.hg38$X3), ]
  subset.hg38$region = paste0(subset.hg38$seqnames, '_', subset.hg38$start, "_", subset.hg38$end)
  subset.hg38 = distinct(subset.hg38, region, X3, .keep_all = TRUE)
  subset.hg38 = subset.hg38[subset.hg38$X4.y %in% 'protein_coding', ]
  # subset.hg38 = subset.hg38[!is.na(subset.hg38$gene_name), ]
  
  # rm(meta.feature, peaks.matched, enriched.motifs, peaks.with.motif, interest.peaks, regions, hg38_PromEnh)
  
  # Enrichment
  DEGs <- as.data.frame(readr::read_csv(paste0("scRNA/memento/", metaprogram, ".csv")))
  DEGs$padjust = p.adjust(DEGs$de_pval)
  DEGs_genes = DEGs$gene[DEGs$de_coef >= log(2) & DEGs$padjust < 0.05]
  genes = intersect(unique(subset.hg38$X3), DEGs_genes)
  return(genes)
}, mc.cores = 7)
names(Genes_dELS) = names(colori_mp)

DF.list.PLS = mclapply(1:7, FUN = function(x){
  metaprogram = names(colori_mp)[x]
  genes = Genes_PLS[[metaprogram]]
  res <- clusterProfiler::enricher(genes,
                                   TERM2GENE = signatures,
                                   minGSSize = 10,
                                   maxGSSize = 250,
                                   pAdjustMethod = "fdr",
                                   pvalueCutoff = 1,
  )
  df = res@result
  df = df[df$pvalue < 0.1, ]
  df = df[order(df$p.adjust, decreasing = FALSE), ]
  df$MSigDB = 'GO:BP'
  df$MSigDB[grep("KEGG_", df$Description)] = 'KEGG'
  df$MSigDB[grep("REACTOME_", df$Description)] = 'REACTOME'
  df$MSigDB[grep("HALLMARK_", df$Description)] = 'HALLMARK'
  df$Label = gsub("GOBP_|KEGG_|REACTOME_|HALLMARK_", "", df$Description)
  df$Label = gsub("_", " ", df$Label)
  # df = df[c(which(df$MSigDB %in% 'GO:BP')[1:50], which(df$MSigDB %in% 'REACTOME')[1:50], 
  #           which(df$MSigDB %in% 'KEGG')[1:50], which(df$MSigDB %in% 'HALLMARK')[1:50]), ]
  # df = na.omit(df)
  df = df[order(df$p.adjust, decreasing = FALSE), ]
  df$MSigDB = factor(df$MSigDB, levels = c('GO:BP', 'REACTOME', 'KEGG', 'HALLMARK'))
  df$size <- cut(
    df$p.adjust,
    breaks = c(-Inf, 0.01, 0.05, Inf),
    labels = c("p.adjust < 0.01", "0.01 < p.adjust < 0.05", "p.adjust > 0.05"),
    right = FALSE
  )
  df$size <- factor(df$size, 
                    levels = c("p.adjust < 0.01", "0.01 < p.adjust < 0.05", "p.adjust > 0.05"), 
                    labels = c(1,2,3)) 
  df$Metaprogram = metaprogram
  return(df)
}, mc.cores = 7)
names(DF.list.PLS) = names(colori_mp)

DF.list.pELS = mclapply(1:7, FUN = function(x){
  metaprogram = names(colori_mp)[x]
  genes = Genes_pELS[[metaprogram]]
  res <- clusterProfiler::enricher(genes,
                                   TERM2GENE = signatures,
                                   minGSSize = 10,
                                   maxGSSize = 250,
                                   pAdjustMethod = "fdr",
                                   pvalueCutoff = 1,
  )
  df = res@result
  df = df[df$pvalue < 0.1, ]
  df = df[order(df$p.adjust, decreasing = FALSE), ]
  df$MSigDB = 'GO:BP'
  df$MSigDB[grep("KEGG_", df$Description)] = 'KEGG'
  df$MSigDB[grep("REACTOME_", df$Description)] = 'REACTOME'
  df$MSigDB[grep("HALLMARK_", df$Description)] = 'HALLMARK'
  df$Label = gsub("GOBP_|KEGG_|REACTOME_|HALLMARK_", "", df$Description)
  df$Label = gsub("_", " ", df$Label)
  # df = df[c(which(df$MSigDB %in% 'GO:BP')[1:50], which(df$MSigDB %in% 'REACTOME')[1:50], 
  #           which(df$MSigDB %in% 'KEGG')[1:50], which(df$MSigDB %in% 'HALLMARK')[1:50]), ]
  # df = na.omit(df)
  df = df[order(df$p.adjust, decreasing = FALSE), ]
  df$MSigDB = factor(df$MSigDB, levels = c('GO:BP', 'REACTOME', 'KEGG', 'HALLMARK'))
  df$size <- cut(
    df$p.adjust,
    breaks = c(-Inf, 0.01, 0.05, Inf),
    labels = c("p.adjust < 0.01", "0.01 < p.adjust < 0.05", "p.adjust > 0.05"),
    right = FALSE
  )
  df$size <- factor(df$size, 
                    levels = c("p.adjust < 0.01", "0.01 < p.adjust < 0.05", "p.adjust > 0.05"), 
                    labels = c(1,2,3)) 
  df$Metaprogram = metaprogram
  return(df)
}, mc.cores = 7)
names(DF.list.pELS) = names(colori_mp)

DF.list.dELS = mclapply(1:7, FUN = function(x){
  metaprogram = names(colori_mp)[x]
  genes = Genes_dELS[[metaprogram]]
  res <- clusterProfiler::enricher(genes,
                                   TERM2GENE = signatures,
                                   minGSSize = 10,
                                   maxGSSize = 250,
                                   pAdjustMethod = "fdr",
                                   pvalueCutoff = 1,
  )
  df = res@result
  df = df[df$pvalue < 0.1, ]
  df = df[order(df$p.adjust, decreasing = FALSE), ]
  df$MSigDB = 'GO:BP'
  df$MSigDB[grep("KEGG_", df$Description)] = 'KEGG'
  df$MSigDB[grep("REACTOME_", df$Description)] = 'REACTOME'
  df$MSigDB[grep("HALLMARK_", df$Description)] = 'HALLMARK'
  df$Label = gsub("GOBP_|KEGG_|REACTOME_|HALLMARK_", "", df$Description)
  df$Label = gsub("_", " ", df$Label)
  # df = df[c(which(df$MSigDB %in% 'GO:BP')[1:50], which(df$MSigDB %in% 'REACTOME')[1:50], 
  #           which(df$MSigDB %in% 'KEGG')[1:50], which(df$MSigDB %in% 'HALLMARK')[1:50]), ]
  # df = na.omit(df)
  df = df[order(df$p.adjust, decreasing = FALSE), ]
  df$MSigDB = factor(df$MSigDB, levels = c('GO:BP', 'REACTOME', 'KEGG', 'HALLMARK'))
  df$size <- cut(
    df$p.adjust,
    breaks = c(-Inf, 0.01, 0.05, Inf),
    labels = c("p.adjust < 0.01", "0.01 < p.adjust < 0.05", "p.adjust > 0.05"),
    right = FALSE
  )
  df$size <- factor(df$size, 
                    levels = c("p.adjust < 0.01", "0.01 < p.adjust < 0.05", "p.adjust > 0.05"), 
                    labels = c(1,2,3)) 
  df$Metaprogram = metaprogram
  return(df)
}, mc.cores = 7)
names(DF.list.dELS) = names(colori_mp)

save(DA.peaks, TF.of.interest, Genes_PLS, Genes_dELS, Genes_pELS, DF.list.PLS, DF.list.dELS, DF.list.pELS, file = 'scATAC/RData/Results/DF.list_TFenrichment_ENCODE.RData')
openxlsx::write.xlsx(DF.list.PLS, file = 'scATAC/RData/Results/TF_promoter_enrichment_ENCODE.xlsx')
openxlsx::write.xlsx(DF.list.dELS, file = 'scATAC/RData/Results/TF_distal_enhancer_enrichment_ENCODE.xlsx')
openxlsx::write.xlsx(DF.list.pELS, file = 'scATAC/RData/Results/TF_proximal_enhancer_enrichment_ENCODE.xlsx')
