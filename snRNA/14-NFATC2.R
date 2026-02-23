###################
#### Libraries ####
###################
setwd("/home3/ciervo/scMULTIOME/Analisi")
# BiocManager::install('saezlab/decoupleR')
library(decoupleR)
library(tidyr)
library(dplyr)
library(parallel)
library(tibble)
library(caret)
library(pROC)
library(ggplot2)

# 1. NFATC2 activity in bulk melanoma ----
load('Melanoma_Bulk/dataframe_bulk_all.RData')

summary(mmat$purity)
mmat = mmat[mmat$purity >= 0.3, ]
mmat$Treatment[mmat$Dataset %in% c('vanallen', 'amato')] = 'PRE'
table(mmat$Dataset, mmat$Therapy, exclude = NULL)
table(mmat$Treatment, exclude = NULL)
# mmat = mmat[!is.na(mmat$Treatment), ]
# mmat = mmat[!(mmat$Dataset %in% 'amato'), ]
samples = rownames(mmat)

llist = list.files('scRNA/RData/Melanoma_bulk/Matrices',full.names = T)[2:9]
load('scRNA/Pozniak_SCENIC/Regulon_SCENIC.RData', verbose = TRUE)

create_df <- function(tf_list) {
  do.call(rbind, lapply(names(tf_list), function(tf_name) {
    tf_data <- tf_list[[tf_name]]
    pos_df <- if (!is.null(tf_data)) {
      data.frame(source = tf_name, target = tf_data, mor = 1, stringsAsFactors = FALSE)
    } else NULL
  }))
}
net = create_df(regulons_list)
net <- net %>% distinct(source, target, mor, .keep_all = TRUE)

# decoupleR - Univariate Linear Model
prova_decoupleR = mclapply(llist, function(x) {
  expr = get(load(x))
  expr = expr[, colnames(expr) %in% samples]
  
  regulons_in <- subset(net, target %in% rownames(expr))
  
  ul <- run_ulm(
    mat = expr,
    network = regulons_in,
    .source = "source", .target = "target"
  )
  
  prova_wide <- ul %>%
    dplyr::select(source, condition, score) %>%
    pivot_wider(names_from = condition, values_from = score) %>% 
    as.data.frame()
  rownames(prova_wide) = prova_wide$source
  prova_wide$source = NULL
  return(prova_wide)
  
}, mc.cores = length(llist))

# Activity - expression
# NFATC2 expression
nfatc2_expr = lapply(llist, function(x) {
  expr = get(load(x))
  expr = expr['NFATC2', colnames(expr) %in% samples]
  return(expr)
})
nfatc2_expr = Reduce(c, nfatc2_expr)
nfatc2_expr = scale(nfatc2_expr)

# Activity
nfatc2 = lapply(prova_decoupleR, function(x) x['NFATC2', ])
nfatc2 = Reduce(c, nfatc2) %>% unlist()
nfatc2[1:5]
sum(is.na(nfatc2))
nfatc2 = nfatc2[!is.na(nfatc2)]
nfatc2 = scale(nfatc2)

# rownames(nfatc2) = gsub('X', '', rownames(nfatc2))
# names(irf1) = gsub('X', '', names(irf1))
# names(mitf) = gsub('X', '', names(mitf))
# names(lef1) = gsub('X', '', names(lef1))

sum(rownames(mmat) %in% rownames(nfatc2))
keep = intersect(rownames(nfatc2), rownames(mmat))
mmat = mmat[keep, ]
nfatc2 = nfatc2[keep, ]
identical(rownames(mmat), names(nfatc2))
mmat$NFATC2_activity = nfatc2

# rownames(nfatc2_expr) = gsub('X', '', rownames(nfatc2_expr))
identical(rownames(mmat), rownames(nfatc2_expr))
sum(rownames(mmat) %in% rownames(nfatc2_expr))
nfatc2_expr = nfatc2_expr[rownames(mmat), , drop = FALSE]
mmat$NFATC2_expression = nfatc2_expr
mmat$Net_ = mmat$NFATC2_activity - mmat$NFATC2_expression

mmat$Treatment[mmat$Treatment %in% c('EDT', 'ON')] = 'POST'
colnames(mmat)
table(mmat$Treatment, mmat$Dataset, exclude = NULL)

save(mmat, file = 'Melanoma_Bulk/BulkICI_NFACT2_activity_PozReg.RData')
# 2. NFATC2 KD ----
# https://doi.org/10.1038/s41388-019-0729-2; https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE101323
library(dplyr)
library(readr)
library(limma)
library(ggplot2)
library(data.table)
library(edgeR)
library(sva)
library(yaGST)
library(ComplexHeatmap)
library(circlize)
source("/home/caruso/scProject/NYnontumor/package_code/functionScRNAseq.R")
source("colori_finali.R")

x <- read.ilmn("Melanoma_Bulk/GSE101323_nfatc2_ko/GSE101323_Non-normalized_data.txt.gz",probeid = "ID_REF", annotation = 'SYMBOL')

pData = data.frame('Sample' = colnames(x$E), 'Type' = c(rep('Ctrl_shRNA 1', 3), rep('Ctrl_shRNA 2', 3), rep('NFATC2_shRNA 1', 3), rep('NFATC2_shRNA 2', 3)), row.names = colnames(x$E))
pData$Type.col = 'lightblue'
pData$Type.col[pData$Type %in% 'NFATC2_shRNA 1'] = 'steelblue'
pData$Type.col[pData$Type %in% 'Ctrl_shRNA 2'] = 'peachpuff'
pData$Type.col[pData$Type %in% 'NFATC2_shRNA 2'] = 'darkorange'
pData$Condition = 'Ctrl'
pData$Condition[grep('NFATC2', pData$Type)] = 'KD'
pData$sh = '1'
pData$sh[grep('shRNA 2', pData$Type)] = '2'
pData$Batch = paste0(pData$Condition, '_', pData$sh)
batch <- as.numeric(as.factor(pData$sh))

# Pre normalization
expr = x$E
boxplot(log2(expr+1), col = pData$Type.col)
pca <- prcomp(t(log2(expr+1)))
pca_data <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2])
pca_data <- merge(pca_data, pData, by = "row.names", all = T)
rownames(pca_data) <- pca_data$Row.names
pca_data$Row.names <- NULL

ggplot(data = pca_data, aes(x = PC1, y = PC2, col = Type.col)) +
  scale_colour_identity("Type", breaks = pca_data$Type.col, 
                        labels = pca_data$Type, guide = "legend") + 
  geom_point(size = 4) +
  theme_bw() + 
  labs() +
  geom_hline(yintercept = 0, color = "grey80", linetype = 'dashed') + 
  geom_vline(xintercept = 0, color = "grey80", linetype = 'dashed')

# Normalization
y <- neqc(x)
expr_norm = y$E
# expr_batch <- removeBatchEffect(expr_norm, batch)
# expr_batch <- sva::ComBat(expr_norm, batch) # uso la prima correzione
boxplot(expr_norm, col = pData$Type.col)
pca <- prcomp(t(expr_norm))
# pca <- prcomp(t(log2(expr_batch+1)))
pca_data <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2])
pca_data <- merge(pca_data, pData, by = "row.names", all = T)
rownames(pca_data) <- pca_data$Row.names
pca_data$Row.names <- NULL

ggplot(data = pca_data, aes(x = PC1, y = PC2, col = Type.col)) +
  scale_colour_identity("Type", breaks = pca_data$Type.col, 
                        labels = pca_data$Type, guide = "legend") + 
  geom_point(size = 4) +
  theme_bw() + 
  labs() +
  geom_hline(yintercept = 0, color = "grey80", linetype = 'dashed') + 
  geom_vline(xintercept = 0, color = "grey80", linetype = 'dashed')

annotation_df = read_delim("Melanoma_Bulk/GSE101323_nfatc2_ko/GPL10558_HumanHT-12_V4_0_R1_15002873_B.txt", 
                           delim = "\t", escape_double = FALSE, 
                           trim_ws = TRUE)
annotation_df <- as.data.frame(annotation_df)

annotation_df <- annotation_df[annotation_df$Probe_Id %in% rownames(expr), ]
annotation_df <- annotation_df[!is.na(annotation_df$Symbol) & annotation_df$Symbol != "", ]
expr[annotation_df$Symbol %in% 'NFATC2', ]
expr <- expr[annotation_df$Probe_Id, ]
expr_df <- cbind(annotation_df["Symbol"], expr)
expr <- expr_df %>%
  group_by(Symbol, .drop = TRUE) %>%
  summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)), .groups = "drop")
expr = as.data.frame(expr)
rownames(expr) <- expr$Symbol
expr$Symbol <- NULL
expr[1:5, ]
expr['NFATC2', ]
expr['HLA-B', ]
expr_raw = expr

# expr_norm = expr_batch
expr_norm <- expr_norm[annotation_df$Probe_Id, ]
expr_df <- cbind(annotation_df["Symbol"], expr_norm)
expr_norm <- expr_df %>%
  group_by(Symbol, .drop = TRUE) %>%
  summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>% 
  as.data.frame()
rownames(expr_norm) <- expr_norm$Symbol
expr_norm$Symbol <- NULL
expr_norm[1:5, ]
expr_norm['MITF', ]
expr_norm['NFATC2', ]
expr_norm['HLA-B', ]

save(expr_raw, expr_norm, pData, file = 'Melanoma_Bulk/GSE101323_nfatc2_ko/expr_data.RData')

# DEGs ----
load('Melanoma_Bulk/GSE101323_nfatc2_ko/expr_data.RData')

# sh1 
pData_sh1 <- pData[grep('1', pData$Type), ]
ddata <- expr_norm[, rownames(pData_sh1)]

groups <- rep("KD", ncol(ddata))
names(groups) <- colnames(ddata)
groups[names(groups) %in% rownames(pData_sh1)[grep('Ctrl', pData_sh1$Type)]] <- "Ctrl"
groups <- as.factor(groups)

design <- model.matrix(~ 0 + groups)
colnames(design) <- levels(groups)

fit <- lmFit(ddata, design)
contrast <- makeContrasts(KD - Ctrl, levels=design)
fit2 <- contrasts.fit(fit, contrast)
fit2 <- eBayes(fit2)
ans <- topTable(fit2, number = nrow(ddata))
ans$DEGs <- ""
ans$DEGs[ans$logFC < -1.5 & ans$adj.P.Val < 0.05] <- "Down regulated"
ans$DEGs[ans$logFC > 1.5 & ans$adj.P.Val < 0.05] <- "Up regulated"
table(ans$DEGs)

ans$Label = ""
ans$Label[ans$DEGs == "Up regulated"] = rownames(ans)[ans$DEGs == "Up regulated"]
ans$Label[ans$DEGs == "Down regulated"] = rownames(ans)[ans$DEGs == "Down regulated"]
save(ans, file = 'Melanoma_Bulk/GSE101323_nfatc2_ko/DEGs_sh1.RData')

# ssMwwGST using the top 50 NMF genes and Tsoi signatures ----
DEGs = parallel::mclapply(1:7, FUN = function(x) {
  metaprogram = names(colori_mp)[x]
  tmp = as.data.frame(readr::read_csv(paste0("scRNA/memento/", metaprogram, ".csv")))
  tmp$padjust = p.adjust(tmp$de_pval)
  tmp = tmp[order(tmp$de_coef, decreasing = TRUE), ]
  tmp_genes = tmp$gene[tmp$de_coef >= log(2) & tmp$padjust < 0.05][1:50]
}, mc.cores = 7)
names(DEGs) = names(colori_mp)
sapply(DEGs, length)

# ssGSEA MP and Tsoi ----
load('Melanoma_Bulk/GSE101323_nfatc2_ko/expr_data.RData')

expr_norm = expr_norm[, rownames(pData)[pData$sh %in% '1']]
expr_norm_log = log2(expr_norm+1)

ssMwwGst(geData = expr_norm, geneSet = DEGs, ncore = 12, minLenGeneSet = 15, filename = 'Melanoma_Bulk/GSE101323_nfatc2_ko/ssMWWGst_top50up')

load('Melanoma_Bulk/GSE101323_nfatc2_ko/ssMWWGst_top50up_MWW.RData')
load('Melanoma_Bulk/GSE101323_nfatc2_ko/expr_data.RData', verbose = TRUE)

range(NES)
summary(NES)
quantile(NES, c(0.1, 0.9))
NES[NES < quantile(NES, 0.1)] = quantile(NES, 0.1)
NES[NES > quantile(NES, 0.9)] = quantile(NES, 0.9)
range(NES)

samples = pData$Sample[pData$sh %in% '1']
samples_col = pData$Type.col[pData$sh %in% '1']
names(samples_col) = pData$Sample[pData$sh %in% '1']

column_ha = HeatmapAnnotation(Condition = samples,
                              col = list(Condition = samples_col), na_col = "white", 
                              annotation_name_side = "left", show_legend = FALSE, 
                              show_annotation_name = FALSE)
col_fun = colorRamp2(c(min(NES), 0,  max(NES)), c("blue", "white", "red")) ### range di colori dell'heatmap

h1 = Heatmap(NES, 
             # width = unit(5, "in"), height = unit(5, "in"),
             heatmap_legend_param = list(legend_direction = "horizontal"),
             row_labels = c('Cell cylce', 'Melanocytic I', 'EMT/Hypoxia', 'Neural crest-like', 'Antigen presentation/Interferon', 'Melanocytic II', 'Wnt/B-catenin'), 
             name = "NES", col = col_fun, # column_split = ssplit,
             top_annotation = column_ha, column_names_rot = 30,
             cluster_rows = F, show_column_names = FALSE, 
             row_title = "", cluster_columns = F)
h1
lgd = Legend(labels = c("Ctrl", "shNFATC2"), 
             title = "Condition", 
             legend_gp = gpar(fill = c('lightblue', 'steelblue'))
)
draw(h1, annotation_legend_list = list(lgd), heatmap_legend_side = "bottom")

# Tsoi
load("/home/caruso/Analisi2021/Maio_EPICA/dati_clinici/Anichini_signatures/Melanoma_Anichini_Signatures.RData", verbose = TRUE)
rm(AntiPD1_response, EMT_and_MELANOMA_diff_MARKERS, Pozniak_immune, Tirosh_signatures)

ssMwwGst(geData = expr_norm, geneSet = TSOI_differentiation_Stage, ncore = 12, minLenGeneSet = 15, filename = 'Melanoma_Bulk/GSE101323_nfatc2_ko/Tsoi')
load('Melanoma_Bulk/GSE101323_nfatc2_ko/Tsoi_MWW.RData')

range(NES)
NES[NES < -1] = -1
NES[NES > 1] = 1

samples = pData$Sample[c(1:3, 7:9)]
samples_col = pData$Type.col[c(1:3, 7:9)]
names(samples_col) = pData$Sample[c(1:3, 7:9)]
ssplit = pData$sh[c(1:3, 7:9)]

column_ha = HeatmapAnnotation(Condition = samples,
                              col = list(Condition = samples_col), na_col = "white",
                              annotation_name_side = "left", show_legend = F)
col_fun = colorRamp2(c(min(NES), mean(NES),  max(NES)), c("blue", "white", "red")) ### range di colori dell'heatmap

h1 = Heatmap(NES, 
             # width = unit(5, "in"), height = unit(5, "in"),
             name = "NES", col = col_fun, 
             column_split = ssplit,
             top_annotation = column_ha, column_names_rot = 30,
             cluster_rows = F,
             row_title = "", cluster_columns = F)
lgd = Legend(labels = c("Ctrl 1", "KO 1" #, 
                        # 'Ctrl 2', 'KO 2'
), 
title = "Condition", 
legend_gp = gpar(fill = c('lightblue', 'steelblue'
                          # , 'peachpuff', 'darkorange'
))
)
draw(h1, annotation_legend_list = lgd)

NES_long = reshape2::melt(NES)
NES_long$Condition = "Ctrl"
NES_long$Condition[grep("G|H|I", NES_long$Var2)] = "KD"

NES_mean <- NES_long %>%
  group_by(Var1, Condition) %>%
  summarise(
    mean_value = mean(value, na.rm = TRUE),
    .groups = "drop"
  )

NES_mean = NES_mean[NES_mean$Condition %in% "Ctrl", ]
order = NES_mean$Var1[order(NES_mean$mean_value, decreasing = TRUE)]
unique(order)

NES_long$Var1 = factor(NES_long$Var1, levels = order)
NES_long$value = as.numeric(NES_long$value)

levels(NES_long$Var1) = c("Neural crest-like", "Neural crest-like\nTransitory", "Transitory", 
                          "Undifferentiated\nNeural crest-like", "Undifferentiated", 
                          "Transitory\nMelanocytic", "Melanocytic")

ggplot(NES_long, aes(x = Condition, y = value, fill = Var1)) +
  geom_boxplot(
    aes(linetype = Condition),
    width = 0.3,
    alpha = 1,
    position = position_dodge(width = 1),
    outliers = FALSE
  ) +
  scale_linetype_manual(
    values = c("Ctrl" = "solid", "KD" = "twodash")
  ) +
  scale_fill_manual(values = nes_cols) +
  theme_minimal(base_size = 11) +
  theme(
    strip.background = element_rect(fill = "#f0f0f0", color = NA),
    strip.text = element_text(size = 10, face = "bold", colour = "black"),
    axis.text.x = element_text(size = 9, face = "bold"),
    axis.text.y = element_text(size = 9),
    axis.title.y = element_text(size = 10, face = "bold"),
    axis.line = element_line(colour = "black"),
    panel.grid = element_blank(),
    legend.position = "none",
    panel.spacing.x = unit(0.2, "lines"),
    panel.spacing.y = unit(0.5, "lines")
  ) +
  facet_grid(~Var1, scales = "free_x", axes = "all_x") +
  # stat_compare_means(
  #   aes(group = Condition),
  #   method = "t.test",
  #   # comparisons = list(c("Ctrl", "KD")),
  #   label = "p.signif",
  #   label.x.npc = "center",
  #   # label.y.npc = 0.9,
  #   hide.ns = FALSE
  # ) +
  labs(
    x = "",
    y = "NES",
    fill = ""
  )
