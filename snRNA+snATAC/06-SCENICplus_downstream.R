# Add AUCell scores SCENICplus ----
library(Seurat)
library(SeuratData)
library(Signac)
library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)

setwd('/home3/ciervo/scMULTIOME/Analisi/')
source('colori_finali.R')

# Load obj ----
load('scRNA_scATAC/RData/seu_multiome_filtered.RData')
seu.multiome = subset(seu.multiome, Metaprogram_assignment %in% paste0('MP_', 1:7))

# Output from Downstream_analysis_SCENICplus.ipynb
# Load eRegulon_metadata ----
eRegulon_metadata <- as.data.frame(read_csv("scRNA_scATAC/SCENICplus/CSV_output/eRegulon_metadata.csv"))
# length(unique(eRegulon_metadata$eRegulon_name)) # 188
# eRegulon_metadata[1:5, ]
# summary(sapply(eRegulon_metadata$Gene_signature_name, function(x) as.numeric(gsub('g)', '', strsplit(x, '\\(')[[1]][2])) ) )
# summary(sapply(eRegulon_metadata$Region_signature_name, function(x) as.numeric(gsub('r)', '', strsplit(x, '\\(')[[1]][2])) ) )

eRegulon_metadata = eRegulon_metadata[eRegulon_metadata$regulation == 1, ]

# Load eRegulon_AUC ----
eRegulon_filtered_AUC <- as.data.frame(read_csv("scRNA_scATAC/SCENICplus/CSV_output/eRegulon_filtered_AUC.csv"))
rownames(eRegulon_filtered_AUC) = eRegulon_filtered_AUC$Cell
head(rownames(eRegulon_filtered_AUC))
rownames(eRegulon_filtered_AUC) = gsub('___cisTopic', '', rownames(eRegulon_filtered_AUC))
eRegulon_filtered_AUC$Cell = NULL
eRegulon_filtered_AUC = as.matrix(eRegulon_filtered_AUC)
eRegulon_filtered_AUC = t(eRegulon_filtered_AUC)
eRegulon_filtered_AUC[1:5, 1:5]

setdiff(colnames(eRegulon_filtered_AUC), colnames(seu.multiome))
colnames(eRegulon_filtered_AUC) = gsub('Pat03_W0_1', 'Pat03_W0', colnames(eRegulon_filtered_AUC))

eRegulon_filtered_AUC_regions = eRegulon_filtered_AUC[grep(pattern = 'r\\)', rownames(eRegulon_filtered_AUC)), ]
eRegulon_filtered_AUC_regions = eRegulon_filtered_AUC_regions[unique(eRegulon_metadata$Region_signature_name), ]
eRegulon_filtered_AUC_genes = eRegulon_filtered_AUC[grep(pattern = 'g\\)', rownames(eRegulon_filtered_AUC)), ]
eRegulon_filtered_AUC_genes = eRegulon_filtered_AUC_genes[unique(eRegulon_metadata$Gene_signature_name), ]

setdiff(colnames(eRegulon_filtered_AUC_genes), colnames(seu.multiome))
rownames(eRegulon_filtered_AUC_genes) = gsub('_', '\\.', rownames(eRegulon_filtered_AUC_genes))
rownames(eRegulon_filtered_AUC_regions) = gsub('_', '\\.', rownames(eRegulon_filtered_AUC_regions))

# Add AUC genes and regions to seu obj ----
seu.multiome[['eRegulon_gene']] = CreateAssayObject(data = eRegulon_filtered_AUC_genes)
seu.multiome[['eRegulon_region']] = CreateAssayObject(data = eRegulon_filtered_AUC_regions)

DefaultAssay(seu.multiome) = 'eRegulon_gene'
seu.multiome = ScaleData(seu.multiome, features = rownames(seu.multiome))
seu.multiome@assays$eRegulon_gene$scale.data[1:5, 1:5]

mean_scaled_auc_genes = AverageExpression(seu.multiome, return.seurat = FALSE, group.by = 'Metaprogram_assignment', assays = 'eRegulon_gene')
mean_scaled_auc_genes = as.matrix(mean_scaled_auc_genes$eRegulon_gene)
mean_scaled_auc_genes = t(scale(t(mean_scaled_auc_genes)))

DefaultAssay(seu.multiome) = 'eRegulon_region'
seu.multiome = ScaleData(seu.multiome, features = rownames(seu.multiome))
seu.multiome@assays$eRegulon_region$scale.data[1:5, 1:5]

mean_scaled_auc_regions = AverageExpression(seu.multiome, return.seurat = FALSE, group.by = 'Metaprogram_assignment', assays = 'eRegulon_region')
mean_scaled_auc_regions = as.matrix(mean_scaled_auc_regions$eRegulon_region)
mean_scaled_auc_regions = t(scale(t(mean_scaled_auc_regions)))

DefaultAssay(seu.multiome) = 'RNA'
seu.multiome = seu.multiome %>% 
  NormalizeData() #%>% 
# ScaleData(features = rownames(seu.multiome))

length(unique(eRegulon_metadata$TF)) # 122
length(intersect(unique(eRegulon_metadata$TF), rownames(seu.multiome))) # 122
mean_scaled_tf_expr = AverageExpression(seu.multiome, return.seurat = FALSE, group.by = 'Metaprogram_assignment', assays = 'RNA', features = unique(eRegulon_metadata$TF))
mean_scaled_tf_expr = as.matrix(mean_scaled_tf_expr$RNA)
mean_scaled_tf_expr = t(scale(t(mean_scaled_tf_expr)))

# Prioritizing TFs ----
# TFs were subsequently ranked based on:
# (1) TF expression in a given MP relative to all other MPs,
# (2) eRegulon AUCell scores on gene expression in a given MP relative to all other MPs,
# (3) eRegulon AUCell scores on chromatin accessibility in a given MP relative to all other MPs.
# Only TFs ranked within the top 30 across all three criteria were retained for downstream analyses.

for(i in nrow(mean_scaled_auc_genes)) {
  for(j in ncol(mean_scaled_auc_genes)){
    mean_scaled_auc_genes[i, j] = mean_scaled_auc_genes[i, j] - mean(mean_scaled_auc_genes[i, -j])
  }
}
Top_30_Regulon_Genes = apply(mean_scaled_auc_genes, 2, function(x) head(names(sort(x, decreasing = TRUE)), 30) )

for(i in nrow(mean_scaled_auc_regions)) {
  for(j in ncol(mean_scaled_auc_regions)){
    mean_scaled_auc_regions[i, j] = mean_scaled_auc_regions[i, j] - mean(mean_scaled_auc_regions[i, -j])
  }
}
Top_30_Regulon_Regions = apply(mean_scaled_auc_regions, 2, function(x) head(names(sort(x, decreasing = TRUE)), 30) )

for(i in nrow(mean_scaled_tf_expr)) {
  for(j in ncol(mean_scaled_tf_expr)){
    mean_scaled_tf_expr[i, j] = mean_scaled_tf_expr[i, j] - mean(mean_scaled_tf_expr[i, -j])
  }
}
Top_30_TF_expr = apply(mean_scaled_tf_expr, 2, function(x) head(names(sort(x, decreasing = TRUE)), 30) )

intersect_MPs = vector('list', 7)
names(intersect_MPs) = paste0('MP-', 1:7)

for(mp in names(intersect_MPs)){
  tmp = data.frame('Expr_gene' = Top_30_TF_expr[, mp], 'AUC_gene' = Top_30_Regulon_Genes[, mp], 'AUC_region' = Top_30_Regulon_Regions[, mp])
  tmp = cbind(tmp, tf = unlist(lapply(tmp[,1], function(x) ifelse(any(grepl(x, tmp[, 2])) & any(grepl(x, tmp[, 3])), 'yes', 'no') ) ) )
  tmp = tmp[tmp$tf %in% 'yes', 'Expr_gene']
  intersect_MPs[[mp]] = tmp
}
intersect_MPs

save(intersect_MPs, eRegulon_filtered_AUC_genes, eRegulon_filtered_AUC_regions, mean_scaled_auc_genes, mean_scaled_auc_regions, file = 'scRNA_scATAC/SCENICplus/downstream_scenicplus.RData')

# writeLines(unlist(MP_TF[[1]]), con = 'scRNA_scATAC/SCENICplus/tf_list_MP1.txt')
# writeLines(unlist(MP_TF[[2]]), con = 'scRNA_scATAC/SCENICplus/tf_list_MP2.txt')
# writeLines(unlist(MP_TF[[3]]), con = 'scRNA_scATAC/SCENICplus/tf_list_MP3.txt')
# writeLines(unlist(MP_TF[[4]]), con = 'scRNA_scATAC/SCENICplus/tf_list_MP4.txt')
# writeLines(unlist(MP_TF[[5]]), con = 'scRNA_scATAC/SCENICplus/tf_list_MP5.txt')
# writeLines(unlist(MP_TF[[6]]), con = 'scRNA_scATAC/SCENICplus/tf_list_MP6.txt')
# writeLines(unlist(MP_TF[[7]]), con = 'scRNA_scATAC/SCENICplus/tf_list_MP7.txt')

# Genes AUC plot ----
load('scRNA_scATAC/SCENICplus/downstream_scenicplus.RData')
RSS_scores <- as.matrix(read_csv("scRNA_scATAC/SCENICplus/CSV_output/RSS_scores.csv"))
rownames(RSS_scores) = c('MP_2', 'MP_7', 'MP_5', 'MP_1', 'MP_4', 'MP_6', "MP_3")
RSS_scores = RSS_scores[paste0('MP_', 1:7), ]
colnames(RSS_scores)

tf_to_plot = Reduce(c, intersect_MPs) %>% unique()

tmp = unique(eRegulon_metadata$Gene_signature_name[eRegulon_metadata$TF %in% tf_to_plot])
RSS_scores = RSS_scores[, tmp]

Rss_melt = reshape2::melt(t(RSS_scores))
colnames(Rss_melt) = c('Gene_signature_name', 'MP', 'RSS')
expr = reshape2::melt(mean_scaled_auc_genes)

rel_data = expr
colnames(rel_data) = c('Gene_signature_name', "MP", 'value')
eRegulon_metadata = eRegulon_metadata[eRegulon_metadata$TF %in% tf_to_plot, ]
eRegulon_metadata = eRegulon_metadata[!duplicated(eRegulon_metadata$TF), ]
setdiff(eRegulon_metadata$Gene_signature_name, tmp) # 0

rel_data$MP = gsub('-', '_', rel_data$MP)
rel_data$Gene_signature_name = gsub('\\.', '_', rel_data$Gene_signature_name)
rel_data = rel_data[rel_data$Gene_signature_name %in% tmp, ]

rel_data = merge(rel_data, eRegulon_metadata[, c('TF', 'eRegulon_name', 'Gene_signature_name')], by = 'Gene_signature_name', all.x = TRUE, all.y = FALSE)
rel_data = merge(rel_data, Rss_melt, by = c('Gene_signature_name', 'MP'), all = TRUE)
rel_data$eRegulon_name = as.character(rel_data$eRegulon_name)
rel_data$eRegulon_name = sapply(rel_data$eRegulon_name, function(x) paste0(strsplit(x, '_')[[1]][1], '(+)') )

range(rel_data$value)
range(rel_data$RSS)

rel_data$value[rel_data$value < -1.5] = -1.5
rel_data$value[rel_data$value > 1.5] = 1.5

rel_data[1:15, ]
rel_data$Target_genes = sapply(rel_data$Gene_signature_name, function(x) strsplit(x, '\\(')[[1]][2])
rel_data$Target_genes = gsub('g)', '', rel_data$Target_genes)
rel_data$Target_genes = log2(as.numeric(rel_data$Target_genes))

rel_data$MP = factor(rel_data$MP, names(colori_mp))

rel_data <- rel_data %>%
  mutate(MP_num = as.numeric(gsub("MP_", "", MP)))

regulon_order <- rel_data %>%
  group_by(eRegulon_name) %>%
  summarize(max_MP = MP_num[which.max(value)], .groups = "drop") %>%
  arrange(max_MP) %>%
  mutate(order = row_number())

rel_data <- rel_data %>%
  left_join(regulon_order, by = "eRegulon_name") %>%
  mutate(eRegulon_name = reorder(eRegulon_name, order))

rel_data$MP <- factor(rel_data$MP, levels = paste0("MP_", 1:7))

range(rel_data$RSS)

p1 = ggplot(data = rel_data, mapping = aes_string(x = 'MP', y = 'eRegulon_name')) +
  geom_tile(mapping = aes_string(fill = 'value')) +
  geom_point(mapping = aes_string(size = 'RSS'), colour="black",pch=21, fill = 'black') +
  scale_radius(range = c(0.3, 5)) +
  scale_fill_distiller(palette = "RdYlBu", limits=c(-1.5,1.5)) +
  theme_pubr() +
  labs(fill = 'Regulon AUC', size = 'Gene AUC RSS') +
  theme( # axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(angle = 60, hjust = 1),
    plot.margin = margin(0, 5, 0, 0),
    axis.title.y = element_blank(),
    axis.ticks.x = element_blank(), 
    axis.title.x = element_blank(), 
    axis.text.x = element_blank(),
    legend.position = 'right') 
p1

p2 = ggplot(rel_data, aes(y = eRegulon_name, x = Target_genes)) +
  geom_col(fill = "grey30") +
  theme_pubr() +
  labs(x = 'log2[target genes]') +
  theme(axis.line.x = element_blank(),
        axis.ticks.x = element_blank(), 
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        plot.margin = margin(0, 5, 0, 0))

(p2/p1 + plot_layout(heights = c(1, 3))) & coord_flip()

# Regions AUC plot ----
eRegulon_metadata <- as.data.frame(read_csv("scRNA_scATAC/SCENICplus/CSV_output/eRegulon_metadata.csv"))
eRegulon_metadata = eRegulon_metadata[eRegulon_metadata$regulation == 1, ]

RSS_scores <- as.matrix(read_csv("scRNA_scATAC/SCENICplus/CSV_output/RSS_scores.csv"))
rownames(RSS_scores) = c('MP_2', 'MP_7', 'MP_5', 'MP_1', 'MP_4', 'MP_6', "MP_3")
RSS_scores = RSS_scores[paste0('MP_', 1:7), ]
colnames(RSS_scores)

tmp = unique(eRegulon_metadata$Region_signature_name[eRegulon_metadata$TF %in% tf_to_plot])
RSS_scores = RSS_scores[, tmp]

Rss_melt = reshape2::melt(t(RSS_scores))
colnames(Rss_melt) = c('Region_signature_name', 'MP', 'RSS')
expr = reshape2::melt(mean_scaled_auc_regions)

rel_data = expr
colnames(rel_data) = c('Region_signature_name', "MP", 'value')
eRegulon_metadata = eRegulon_metadata[eRegulon_metadata$TF %in% tf_to_plot, ]
eRegulon_metadata = eRegulon_metadata[!duplicated(eRegulon_metadata$TF), ]
setdiff(eRegulon_metadata$Region_signature_name, tmp) # 0

rel_data$MP = gsub('-', '_', rel_data$MP)
rel_data$Region_signature_name = gsub('\\.', '_', rel_data$Region_signature_name)
rel_data = rel_data[rel_data$Region_signature_name %in% tmp, ]

rel_data = merge(rel_data, eRegulon_metadata[, c('TF', 'eRegulon_name', 'Region_signature_name')], by = 'Region_signature_name', all.x = TRUE, all.y = FALSE)
rel_data = merge(rel_data, Rss_melt, by = c('Region_signature_name', 'MP'), all = TRUE)
rel_data$eRegulon_name = as.character(rel_data$eRegulon_name)
rel_data$eRegulon_name = sapply(rel_data$eRegulon_name, function(x) paste0(strsplit(x, '_')[[1]][1], '(+)') )

range(rel_data$value)
range(rel_data$RSS)

rel_data$value[rel_data$value < -1.5] = -1.5
rel_data$value[rel_data$value > 1.5] = 1.5

rel_data[1:15, ]
rel_data$Target_genes = sapply(rel_data$Region_signature_name, function(x) strsplit(x, '\\(')[[1]][2])
rel_data$Target_genes = gsub('r)', '', rel_data$Target_genes)
rel_data$Target_genes = log2(as.numeric(rel_data$Target_genes))

rel_data$MP = factor(rel_data$MP, names(colori_mp))

rel_data <- rel_data %>%
  mutate(MP_num = as.numeric(gsub("MP_", "", MP)))

regulon_order <- rel_data %>%
  group_by(eRegulon_name) %>%
  summarize(max_MP = MP_num[which.max(value)], .groups = "drop") %>%
  arrange(max_MP) %>%
  mutate(order = row_number())

rel_data <- rel_data %>%
  left_join(regulon_order, by = "eRegulon_name") %>%
  mutate(eRegulon_name = reorder(eRegulon_name, order))

rel_data$MP <- factor(rel_data$MP, levels = paste0("MP_", 1:7))

range(rel_data$RSS)

p1 = ggplot(data = rel_data, mapping = aes_string(x = 'MP', y = 'eRegulon_name')) +
  geom_tile(mapping = aes_string(fill = 'value')) +
  geom_point(mapping = aes_string(size = 'RSS'), colour="black",pch=21, fill = 'black') +
  scale_radius(range = c(0.3, 5)) +
  scale_fill_gradient2(low = 'blue', mid = 'grey90', high = 'yellow', limits=c(-1.5,1.5)) +
  # scale_fill_distiller(palette = "YlGnBu", limits=c(-1.5,1.5)) +
  theme_pubr() +
  labs(fill = 'Regulon AUC', size = 'Region AUC RSS') +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(angle = 60, hjust = 1),
        plot.margin = margin(0, 5, 0, 0),
        legend.position = 'right') 

p2 = ggplot(rel_data, aes(y = eRegulon_name, x = Target_genes)) +
  geom_col(fill = "grey30") +
  theme_pubr() +
  labs(x = 'log2[target regions]') +
  theme(axis.line.x = element_blank(),
        axis.ticks.x = element_blank(), 
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        plot.margin = margin(0, 5, 0, 0))

(p2/p1 + plot_layout(heights = c(1, 3))) & coord_flip()

# Overlap between MP signatures (50 genes) and Regulon per each MP ----
load('scRNA_scATAC/SCENICplus/downstream_scenicplus.RData', verbose = TRUE)
rm(eRegulon_filtered_AUC_genes, eRegulon_filtered_AUC_regions, mean_scaled_auc_genes, mean_scaled_auc_regions)

load('scRNA/NMF/metaprograms.RData')
rm(df)

eRegulon_metadata <- as.data.frame(read_csv("scRNA_scATAC/SCENICplus/CSV_output/eRegulon_metadata.csv"))
eRegulon_metadata = eRegulon_metadata[eRegulon_metadata$regulation == 1, ]
names(intersect_MPs) = gsub('-', '_', names(intersect_MPs))

overlap_df = matrix(NA, nrow = 7, ncol = length(unique(Reduce(c, intersect_MPs))))
rownames(overlap_df) = names(intersect_MPs)
colnames(overlap_df) = unique(Reduce(c, intersect_MPs))

for(tf in colnames(overlap_df)){
  for(mp in rownames(overlap_df)){
    regulon_genes = eRegulon_metadata[eRegulon_metadata$TF %in% tf, 'Gene']
    mp_genes = MP_list[[mp]]
    overlap_df[mp, tf] = length(intersect(regulon_genes, mp_genes))
  }
}

overlap_df = t(overlap_df)
overlap_df_melted = melt(overlap_df)
overlap_df_melted$Label = overlap_df_melted$value
overlap_df_melted$Label[overlap_df_melted$Label == 0] = ''

order_tf = rel_data$eRegulon_name
order_tf = gsub('\\(\\+\\)', '', order_tf)

overlap_df_melted$Var1 = factor(overlap_df_melted$Var1, levels = order_tf)

p3 = ggplot(overlap_df_melted, aes(x=Var2, y=Var1, fill=value)) +
  geom_tile(color="white") +
  geom_text(aes(label=Label), color="black") +
  scale_fill_gradient(low="white", high="#006766", limits = c(0,max(overlap_df_melted$value)) ) +
  labs(x="", y="", fill="Overlap") +
  theme_pubr() +
  theme(plot.margin = margin(0, 5, 0, 0),
        # panel.grid = element_blank(),
        # text = element_text(size = 15, face = 'bold'),
        # legend.text = element_text(size = 10), 
        # legend.title = element_text(size = 10),
        # legend.direction = 'horizontal', 
        legend.position = 'right', 
        # legend.title.position = 'top',
        axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(angle = 60, hjust = 1)
  )
(p2/p1/p3 + plot_layout(heights = c(1, 3, 3))) & coord_flip()