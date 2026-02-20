########################
##### Figure Script ####
########################
# setwd('/home3/ciervo/scMULTIOME/Analisi/')
# dir.create('Figures')

###################
#### Libraries ####
###################
library(Seurat)
library(SeuratObject)
library(Signac)
library(parallel)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(viridis)
library(ComplexHeatmap)
library(circlize)
library(clusterProfiler)
library(readr)
library(reshape2)
library(patchwork)
library(survival)
library(survminer)
library(caret)
library(pROC)
library(tibble)
library(scales)
library(stringr)

source('colori_finali.R')

##################
#### Figure 1 ####
##################
# dir.create('Figures/Fig1')
# 1B tSNE malignant/non malignant ----
load('scRNA_scATAC/RData/seu_multiome_filtered.RData')
dim(seu.multiome)
# table(seu.multiome$Metaprogram_assignment)
# Reductions(seu.multiome)

seu.multiome$Malignant = factor(seu.multiome$Malignant, levels = names(colors_malignant))
p = DimPlot(seu.multiome, cols = colors_malignant, group.by = 'Malignant', reduction = 'tsne_multimodal', pt.size = .5) + 
  ggtitle('Malignant')
p
ggsave(filename = 'Figures/Fig1/1b_tSNE_malignant.pdf', plot = p, device = 'pdf', width = 10, height = 8,
       dpi = 600, units = 'in', bg = 'white')

# 1C: tSNE TME ----
load('scRNA/RData/merged_object.RData')
tme = subset(seu, Cell_annotation %in% c('Malignant'), invert = TRUE)
dim(tme)
table(tme$Cell_annotation, exclude = NULL)

tme = tme %>% 
  NormalizeData() %>% 
  FindVariableFeatures(nfeatures = 3000) %>% 
  ScaleData(features = rownames(tme)) %>% 
  RunPCA() %>% 
  FindNeighbors(dims = 1:12) %>% 
  RunTSNE(dims = 1:12)

tme$Cell_annotation = factor(tme$Cell_annotation, levels = names(colors_tme))
p = DimPlot(tme, cols = colors_tme, group.by = 'Cell_annotation', reduction = 'tsne', pt.size = .5) + 
  ggtitle('Cell type')
p
ggsave(filename = 'Figures/Fig1/1c_tSNE_TME_RNA.pdf', plot = p, device = 'pdf', width = 10, height = 8,
       dpi = 600, units = 'in', bg = 'white')

# 1D: Markers TME ----
# Heatmap expression
load('scRNA/RData/merged_object.RData')
seu = subset(seu, Cell_annotation %in% c('Cycling cells', 'Malignant'), invert = TRUE)
table(seu$Cell_annotation, exclude = NULL)

genes_to_plot = c("MS4A1", 'PAX5', # 'BANK1', 
                  "MZB1", 'JCHAIN', # Plasma cells; # B cells
                  'CD4', 'IL7R', # CD4 & CD4 Tcm #'FOXP1', 
                  'LAG3', 'HAVCR2', 'PDCD1', 'TIGIT', # Exhausted 
                  "CD8A", 'CD8B', "GZMK", 'IL9R', 'RORC', # CD8 & CD8 Tcm
                  "FOXP3", 'IL2RA', # Tregs
                  'IL1B', 'CD86', # M1
                  'MARCO', 'CTSD',  # M2
                  "PDGFRA", 'PDPN', "ACTA2", "RGS5", # iCAF; myoCAF
                  "PECAM1", 'VWF', # Endothelial cells
                  "KRT14", 'KRT1' # Keratinocytes
)

mat = matrix(NA, length(genes_to_plot), 14)
rownames(mat) = genes_to_plot
colnames(mat) = unique(seu$Cell_annotation)
table(seu$Cell_annotation)

for (type in colnames(mat)) {
  x = rownames(seu@meta.data)[seu$Cell_annotation %in% type]
  tmp = seu@assays$RNA$data[genes_to_plot, x]
  # tmp = seu@assays$Imputed_Magic$data[genes_to_plot, x]
  mat[, type] = apply(tmp, 1, function(x) mean(x))
}

mat_scaled = t(scale(t(mat)))
range(mat_scaled)
# quantile(mat_scaled, c(.02, .98))
# mat_scaled[mat_scaled < -2] = -2
# mat_scaled[mat_scaled > 2] = 2

metaprogram <- factor(colnames(mat_scaled), levels = names(colors_tme))
# names(metaprogram) <- rownames(seu@meta.data)
metaprogram_colors <- colors_tme

column_ha = HeatmapAnnotation(Metaprograms = metaprogram,
                              col = list(Metaprograms = metaprogram_colors),
                              annotation_name_side = "left", 
                              show_annotation_name = FALSE, annotation_label = 'Cell type', show_legend = F)
# split_rows = genes_to_plot
split_rows = c(rep('B cells', 2), 
               rep('Plasme cells', 2),
               rep('CD4', 2), 
               rep('Exhausted', 4),
               rep('CD8', 5),
               rep('Tregs', 2),
               rep('M1', 2), 
               rep('M2', 2),
               rep('CAF', 4), 
               rep('EC', 2), 
               rep('Keratinocytes', 2)) 
split_rows = factor(split_rows, levels = unique(split_rows))

ccol = circlize::colorRamp2(c(-2, 0, 2), c('blue', 'white', 'red'))
ggap = c(0.5, 1.5, 0.5, 1.5, 1.5, 0.5, 1.5, 1.5, 0.5, 1.5, 0.5, 1.5, 1.5, 1.5)

lgd = Legend(col_fun = ccol, title = "Average\nexpression", at = c(-2, 0, 2), 
             direction = "horizontal")
ht <- Heatmap(mat_scaled, name = 'Average\nexpression',
              heatmap_legend_param = list(legend_direction = "horizontal"),
              col = ccol,
              cluster_rows = FALSE,
              cluster_columns = FALSE, 
              row_split = split_rows, 
              row_title = NULL,
              show_column_names = FALSE, 
              column_title = NULL,
              column_title_rot = 90,
              column_gap = unit(ggap, "mm"),
              gap = unit(1.5, 'mm'),
              top_annotation = column_ha, 
              column_split = metaprogram,
              row_names_side = "left"
)
ht

pdf(file = 'Figures/Fig1/1d_Heatmap_TME.pdf', width = 6, height = 9)
draw(ht, heatmap_legend_side = "bottom")
dev.off()

# Dotplot ATAC
load('scATAC/RData/scATAC_tme.RData')
tme$Cell_annotation = factor(tme$Cell_annotation, levels = names(colors_tme))

p = DotPlot(tme, features = rev(genes_to_plot), group.by = 'Cell_annotation') + 
  labs(x = "", y = "") +
  guides(size = guide_legend(title = "% accessibility", title.position = "top", position = 'bottom', direction = 'horizontal'), 
         color = guide_colorbar(title = "Scaled activity", title.position = "top", position = 'bottom', direction = 'horizontal')) + 
  RotatedAxis() + coord_flip() +
  viridis::scale_color_viridis() +
  theme(legend.title.align = 0.5,
        legend.spacing.x = unit(5, "lines"),
        axis.text.y = element_text(hjust = 1, vjust = 0.5, angle = 0))
p
ggsave(filename = 'Figures/Fig1/1d_Dotplot_ATAC_tme.pdf', plot = p, device = 'pdf', width = 5, height = 7,
       dpi = 600, units = 'in', bg = 'white')

# 1E: Distribution TME----
# Barplot
load('scRNA_scATAC/RData/seu_multiome_filtered.RData')
df.tme = seu.multiome@meta.data[!seu.multiome$Cell_annotation %in% c('Malignant'), ]
df.tme = df.tme[!is.na(df.tme$Cell_annotation), ]

sample_props <- df.tme %>%
  group_by(orig.ident, Week, Responder, Cell_annotation) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(orig.ident) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

avg_props <- sample_props %>%
  group_by(Week, Responder, Cell_annotation) %>%
  summarise(mean_prop = mean(prop), .groups = "drop")

norm_props <- avg_props %>%
  group_by(Week, Responder) %>%
  mutate(percent = 100 * mean_prop / sum(mean_prop)) %>%
  ungroup()

norm_props$Cell_annotation = factor(norm_props$Cell_annotation, names(colors_tme))

p = ggplot(norm_props, aes(x = Week, y = percent, fill = Cell_annotation)) +
  geom_bar(stat = "identity", position = "stack", color = "white", width = 0.9) +
  facet_grid(~ Responder, scales = "free_x", space = "fixed") +
  scale_fill_manual(values = colors_tme) +
  labs(
    x = "Week",
    y = "Cell type composition (%)",
    fill = NULL
  ) +
  scale_y_continuous(
    expand = c(0, 0),
    limits = c(0, 100)
  ) +
  cowplot::theme_cowplot() +
  theme(
    axis.title = element_text(size = 14, face = "plain", color = "black"),
    axis.text = element_text(size = 12, color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.line = element_line(color = "black"),
    
    strip.background = element_blank(),
    strip.text = element_text(size = 14, face = "bold", color = "black"),
    
    legend.position = "right",
    legend.text = element_text(size = 12),
    legend.key.size = unit(0.6, "cm"),
    
    panel.spacing = unit(0.5, "lines"),
    plot.margin = margin(10, 10, 10, 10)
  )
p

ggsave(filename = 'Figures/Fig1/1e_Barplot_TME.pdf', plot = p, device = 'pdf', width = 8, height = 8,
       dpi = 600, units = 'in', bg = 'white')

################################
#### Supplementary Figure 1 ####
################################
# S1B: tSNE Week ----
seu.multiome$Week = factor(seu.multiome$Week, levels = names(week_cols))
p = DimPlot(seu.multiome, cols = week_cols, group.by = 'Week', reduction = 'tsne_multimodal', pt.size = .5) + 
  ggtitle('Week')
p
ggsave(filename = 'Figures/Fig1/Supp_1b_tSNE_week.pdf', plot = p, device = 'pdf', width = 10, height = 8,
       dpi = 600, units = 'in', bg = 'white')

# S1C: PD1 signalling CD4-CD8 ----
load('scRNA/RData/merged_object.RData')
seu = subset(seu, Cell_annotation %in% c('T cells - CD4', 'T cells - CD4 Tcm', 'T cells - CD8', 'T cells - CD8 Tcm'))
table(seu$Cell_annotation, exclude = NULL)

load('scATAC/RData/scATAC_tme.RData')
tme = subset(tme, Cell_annotation %in% c('T cells - CD4', 'T cells - CD4 Tcm', 'T cells - CD8', 'T cells - CD8 Tcm'))
table(tme$Cell_annotation, exclude = NULL)

# Enrichment
signatures <- read.gmt(gmtfile = "MSigDb_11Feb2025/c2.cp.reactome.v2024.1.Hs.symbols.gmt")
pd1_signalling = signatures$gene[signatures$term %in% 'REACTOME_PD_1_SIGNALING']
ctla4_signalling = signatures$gene[signatures$term %in% 'REACTOME_CTLA4_INHIBITORY_SIGNALING']

seu$Test = seu$Cell_annotation
seu$Test[grep('CD4', seu$Test)] = 'T cells - CD4'
seu$Test[grep('CD8', seu$Test)] = 'T cells - CD8'
seu = AddModuleScore(seu, features = list(PD1 = pd1_signalling), name = 'PD', seed = 123, search = T)
seu = AddModuleScore(seu, features = list(CTLA4 = ctla4_signalling), name = 'CTLA4_', seed = 123, search = T)

Idents(seu) = 'Test'
table(Idents(seu))
seu$Week = factor(seu$Week, names(week_cols))

df_long <- seu@meta.data %>%
  select(PD1, CTLA4_1, Test, Week)
df_long$Week = factor(df_long$Week, levels = names(week_cols))

p1 = ggplot(df_long[df_long$Test %in% 'T cells - CD4', ], aes(x = Week, y = PD1, fill = Week)) +
  geom_violin(trim = FALSE, color = "black", alpha = 0.8) +
  geom_boxplot(width = 0.1, outlier.shape = NA, color = "black", fill = "white", lwd = 0.3) +
  # geom_jitter(width = 0.15, size = 0.6, alpha = 0.5, shape = 21, stroke = 0.2, color = "black", fill = "grey90") +
  scale_fill_manual(values = week_cols) +
  # geom_text(data = pvalues, aes(x = 1.5, y = 4.5, label = paste0('p=', round(p, 3))), inherit.aes = FALSE, size = 3.5) +
  # geom_text(data = pvalues, aes(x = 1.5, y = 4.5, label = p_label), inherit.aes = FALSE, size = 5) +
  # ggh4x::facet_nested_wrap(~Gene + Week, scales = "free_y", nrow = 2) +
  stat_compare_means(label = 'p.signif', comparisons = list('C1' = c('Week 0', 'Week 4'), 'C2' = c('Week 0', 'Week 12')), label.y = c(0.55, 0.65)) +
  theme_minimal(base_size = 11) +
  labs(
    y = "Module score",
    x = NULL,
    title = 'Gene expression'
  ) +
  theme(
    strip.background = element_rect(fill = "#f0f0f0", color = NA),
    strip.text = element_text(size = 10, face = "bold", colour = 'black'),
    axis.text.x = element_text(size = 9, face = "bold"),
    axis.text.y = element_text(size = 9),
    axis.title.y = element_text(size = 10, face = "bold"),
    axis.line = element_line(colour = 'black'),
    panel.grid = element_blank(),
    legend.position = "none", 
    panel.spacing.x = unit(0.2, "lines"),
    panel.spacing.y = unit(0.5, "lines")
  )
p1

p2 = ggplot(df_long[df_long$Test %in% 'T cells - CD8', ], aes(x = Week, y = PD1, fill = Week)) +
  geom_violin(trim = FALSE, color = "black", alpha = 0.8) +
  geom_boxplot(width = 0.1, outlier.shape = NA, color = "black", fill = "white", lwd = 0.3) +
  # geom_jitter(width = 0.15, size = 0.6, alpha = 0.5, shape = 21, stroke = 0.2, color = "black", fill = "grey90") +
  scale_fill_manual(values = week_cols) +
  # geom_text(data = pvalues, aes(x = 1.5, y = 4.5, label = paste0('p=', round(p, 3))), inherit.aes = FALSE, size = 3.5) +
  # geom_text(data = pvalues, aes(x = 1.5, y = 4.5, label = p_label), inherit.aes = FALSE, size = 5) +
  # ggh4x::facet_nested_wrap(~Gene + Week, scales = "free_y", nrow = 2) +
  stat_compare_means(label = 'p.signif', comparisons = list('C1' = c('Week 0', 'Week 4'), 'C2' = c('Week 0', 'Week 12')), label.y = c(1, 1.1)) +
  theme_minimal(base_size = 11) +
  labs(
    y = "Module score",
    x = NULL,
    title = 'Gene expression'
  ) +
  theme(
    strip.background = element_rect(fill = "#f0f0f0", color = NA),
    strip.text = element_text(size = 10, face = "bold", colour = 'black'),
    axis.text.x = element_text(size = 9, face = "bold"),
    axis.text.y = element_text(size = 9),
    axis.title.y = element_text(size = 10, face = "bold"),
    axis.line = element_line(colour = 'black'),
    panel.grid = element_blank(),
    legend.position = "none", 
    panel.spacing.x = unit(0.2, "lines"),
    panel.spacing.y = unit(0.5, "lines")
  )
p2

DefaultAssay(tme) = 'Activity'
tme$Test = tme$Cell_annotation
tme$Test[grep('CD4', tme$Test)] = 'T cells - CD4'
tme$Test[grep('CD8', tme$Test)] = 'T cells - CD8'
tme = AddModuleScore(tme, features = list(PD1 = pd1_signalling), name = 'PD', seed = 123, search = T)

Idents(tme) = 'Test'
table(Idents(tme))
tme$Week = factor(tme$Week, names(week_cols))

df_long <- tme@meta.data %>%
  select(PD1, Test, Week)
df_long$Week = factor(df_long$Week, levels = names(week_cols))

p3 = ggplot(df_long[df_long$Test %in% 'T cells - CD4', ], aes(x = Week, y = PD1, fill = Week)) +
  geom_violin(trim = FALSE, color = "black", alpha = 0.8) +
  geom_boxplot(width = 0.1, outlier.shape = NA, color = "black", fill = "white", lwd = 0.3) +
  # geom_jitter(width = 0.15, size = 0.6, alpha = 0.5, shape = 21, stroke = 0.2, color = "black", fill = "grey90") +
  scale_fill_manual(values = week_cols) +
  # geom_text(data = pvalues, aes(x = 1.5, y = 4.5, label = paste0('p=', round(p, 3))), inherit.aes = FALSE, size = 3.5) +
  # geom_text(data = pvalues, aes(x = 1.5, y = 4.5, label = p_label), inherit.aes = FALSE, size = 5) +
  # ggh4x::facet_nested_wrap(~Gene + Week, scales = "free_y", nrow = 2) +
  stat_compare_means(label = 'p.signif', comparisons = list('C1' = c('Week 0', 'Week 4'), 'C2' = c('Week 0', 'Week 12')), label.y = c(0.45, 0.55)) +
  theme_minimal(base_size = 11) +
  labs(
    y = "Module score",
    x = NULL,
    title = 'Accessibility'
  ) +
  theme(
    strip.background = element_rect(fill = "#f0f0f0", color = NA),
    strip.text = element_text(size = 10, face = "bold", colour = 'black'),
    axis.text.x = element_text(size = 9, face = "bold"),
    axis.text.y = element_text(size = 9),
    axis.title.y = element_text(size = 10, face = "bold"),
    axis.line = element_line(colour = 'black'),
    panel.grid = element_blank(),
    legend.position = "none", 
    panel.spacing.x = unit(0.2, "lines"),
    panel.spacing.y = unit(0.5, "lines")
  )
p3

p4 = ggplot(df_long[df_long$Test %in% 'T cells - CD8', ], aes(x = Week, y = PD1, fill = Week)) +
  geom_violin(trim = FALSE, color = "black", alpha = 0.8) +
  geom_boxplot(width = 0.1, outlier.shape = NA, color = "black", fill = "white", lwd = 0.3) +
  # geom_jitter(width = 0.15, size = 0.6, alpha = 0.5, shape = 21, stroke = 0.2, color = "black", fill = "grey90") +
  scale_fill_manual(values = week_cols) +
  # geom_text(data = pvalues, aes(x = 1.5, y = 4.5, label = paste0('p=', round(p, 3))), inherit.aes = FALSE, size = 3.5) +
  # geom_text(data = pvalues, aes(x = 1.5, y = 4.5, label = p_label), inherit.aes = FALSE, size = 5) +
  # ggh4x::facet_nested_wrap(~Gene + Week, scales = "free_y", nrow = 2) +
  stat_compare_means(label = 'p.signif', comparisons = list('C1' = c('Week 0', 'Week 4'), 'C2' = c('Week 0', 'Week 12')), label.y = c(0.6, 0.7)) +
  theme_minimal(base_size = 11) +
  labs(
    y = "Module score",
    x = NULL,
    title = 'Accessibility'
  ) +
  theme(
    strip.background = element_rect(fill = "#f0f0f0", color = NA),
    strip.text = element_text(size = 10, face = "bold", colour = 'black'),
    axis.text.x = element_text(size = 9, face = "bold"),
    axis.text.y = element_text(size = 9),
    axis.title.y = element_text(size = 10, face = "bold"),
    axis.line = element_line(colour = 'black'),
    panel.grid = element_blank(),
    legend.position = "none", 
    panel.spacing.x = unit(0.2, "lines"),
    panel.spacing.y = unit(0.5, "lines")
  )
p4

p = cowplot::plot_grid(p1, p3, p2, p4, align = 'hv', ncol = 2)

ggsave(filename = 'Figures/Fig1/Supp_1c_Vln_PD1_TCells.pdf', plot = p, device = 'pdf', width = 9, height = 6,
       dpi = 600, units = 'in', bg = 'white')

# S1D: ICI markers Weeks in T cells ----
df = seu@meta.data[seu$Cell_annotation %in% c('T cells - CD4', 'T cells - CD4 Tcm', 'T cells - Exhausted T cell', 'T cells - CD8', 'T cells - CD8 Tcm', 'T cells - Tregs'), c("Week","Responder","Cell_annotation")]
exh_genes = c('PDCD1', 'CTLA4', 'LAG3', 'HAVCR2',  'TIGIT', 'BTLA')
tmp = t(seu@assays$RNA$data[exh_genes, rownames(df)])
identical(rownames(df), rownames(tmp))
df = cbind(df, tmp)
rm(tmp)
head(df)
#                             Week          Responder   Cell_annotation       PDCD1 CTLA4 LAG3  HAVCR2    TIGIT    BTLA
# Pat02_W0_AAAGCTTGTTAACACG-1 Week 0        NR          T cells - CD4 Tcm     0     0    0      0         1.559528 0.00000
# Pat02_W0_ACAACACTCCGGCTAA-1 Week 0        NR          T cells - CD4         0     0    0      0         0.000000 0.00000
# Pat02_W0_ACATTAGTCTAGCGAT-1 Week 0        NR          T cells - CD8 Tcm     0     0    0      0         2.517899 0.00000
# Pat02_W0_AGAACCGCACCTCGCT-1 Week 0        NR          T cells - CD4 Tcm     0     0    0      0         0.000000 1.37434
# Pat02_W0_AGCTTCCTCAGCACGC-1 Week 0        NR          T cells - CD4 Tcm     0     0    0      0         0.000000 0.00000
# Pat02_W0_AGGATCCGTTAATGCG-1 Week 0        NR          T cells - CD4         0     0    0      0         0.000000 0.00000

seu$Responder = factor(seu$Responder, c('R', 'NR'))
seu$Week = factor(seu$Week, levels = names(week_cols))
VlnPlot(object = seu, features = exh_genes, cols = week_cols, group.by = 'Week', split.by = 'Responder', split.plot = TRUE, ncol = 5)

df_long <- df %>%
  pivot_longer(cols = c(PDCD1, CTLA4, LAG3, HAVCR2, TIGIT, BTLA),
               names_to = "Gene",
               values_to = "Expression")
df_long$Responder = factor(df_long$Responder, levels = c('R', 'NR'))
df_long$Week = factor(df_long$Week, levels = names(week_cols))
head(df_long)
#   Week   Responder Cell_annotation   Gene         Expression
# 1 Week 0 NR        T cells - CD4 Tcm PDCD1        0   
# 2 Week 0 NR        T cells - CD4 Tcm CTLA4        0   
# 3 Week 0 NR        T cells - CD4 Tcm LAG3         0   
# 4 Week 0 NR        T cells - CD4 Tcm HAVCR2       0   
# 5 Week 0 NR        T cells - CD4 Tcm TIGIT        1.56
# 6 Week 0 NR        T cells - CD4 Tcm BTLA         0   

pvalues <- df_long %>%
  group_by(Gene, Week) %>%
  summarise(
    p = tryCatch(t.test(Expression ~ Responder)$p.value, error = function(e) NA),
    .groups = "drop"
  ) %>%
  mutate(
    p_label = case_when(
      p < 0.001 ~ "***",
      p < 0.01  ~ "**",
      p < 0.05  ~ "*",
      TRUE ~ "ns"
    )
  )

# Violin plot
p = ggplot(df_long, aes(x = Responder, y = Expression, fill = Responder)) +
  geom_violin(trim = FALSE, color = "black", alpha = 0.8) +
  geom_boxplot(width = 0.1, outlier.shape = NA, color = "black", fill = "white", lwd = 0.3) +
  geom_jitter(width = 0.15, size = 0.6, alpha = 0.5, shape = 21, stroke = 0.2, color = "black", fill = "grey90") +
  scale_fill_manual(values = c('NR' = 'cadetblue', 'R' = 'salmon')) +
  # geom_text(data = pvalues, aes(x = 1.5, y = 4.5, label = paste0('p=', round(p, 3))), inherit.aes = FALSE, size = 3.5) +
  geom_text(data = pvalues, aes(x = 1.5, y = 4.5, label = p_label), inherit.aes = FALSE, size = 5) +
  ggh4x::facet_nested_wrap(~Gene + Week, scales = "free_y", nrow = 2) +
  # stat_compare_means(label = 'p.signif', ref.group = 'NR', hide.ns = TRUE) +
  theme_minimal(base_size = 11) +
  labs(
    y = "Expression Level",
    x = NULL,
    title = NULL
  ) +
  theme(
    strip.background = element_rect(fill = "#f0f0f0", color = NA),
    strip.text = element_text(size = 10, face = "bold", colour = 'black'),
    axis.text.x = element_text(size = 9, face = "bold"),
    axis.text.y = element_text(size = 9),
    axis.title.y = element_text(size = 10, face = "bold"),
    axis.line = element_line(colour = 'black'),
    panel.grid = element_blank(),
    legend.position = "none",
    panel.spacing.x = unit(0.2, "lines"),
    panel.spacing.y = unit(0.5, "lines")
  )
p

ggsave(filename = 'Figures/Fig1/Supp_1d_IC_WeekR.pdf', plot = p, device = 'pdf', width = 7.5, height = 6,
       dpi = 600, units = 'in', bg = 'white')

##################
#### Figure 2 ####
##################
# dir.create('Figures/Fig2')
# 2A: Jaccard Similarity NMF ----
custom_magma <- c(colorRampPalette(c("white", rev(magma(323, begin = 0.15))[1]))(10), rev(magma(323, begin = 0.18)))

load("scRNA/NMF/metaprograms.RData")
load("scRNA/NMF/Program_cluster_list.RData")
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

pdf(file = 'Figures/Fig2/2a_NMF_heatmap.pdf', width = 10, height = 10)
h1
dev.off()

# 2B: Heatmap NMF Scores ----
lf = list.files("scRNA/NMF/Factors/", all.files = F, recursive = F, full.names = F,
                pattern = "_factors")
llist = list()
for(file in lf){
  load(paste0("scRNA/NMF/Factors/", file))
  sample = gsub("_factors.RData", "", file)
  llist[[sample]] = mat
  rm(mat)
}

load('scRNA/NMF/Program_cluster_list.RData')
program_list = Reduce(c, Cluster_list)
load('scRNA/NMF/metaprograms.RData')
ggenes = unique(Reduce(c, MP_list))

llist_2 = lapply(llist, function(x) x[, which(colnames(x) %in% program_list) ])
tmp = llist_2[[1]]
dim(tmp)
llist_2[[1]] = NULL

for(i in 1:length(llist_2)){
  tmp = merge(x = tmp, y = llist_2[[i]], by = 'row.names', all = TRUE)
  rownames(tmp) = tmp$Row.names
  tmp$Row.names = NULL
}

range(tmp, na.rm = TRUE)
tmp[is.na(tmp)] = 0
length(intersect(rownames(tmp), ggenes))
setdiff(ggenes, rownames(tmp))
# tmp = tmp / rowMeans(tmp)

ccols = colorRamp2(c(0, 20, 60), c('grey90', '#add8e6', '#006766'))

program_order = lapply(1:7, function(x) rep(paste0('MP_', x), length(Cluster_list[[x]]) ))
program_order = Reduce(c, program_order)
names(program_order) = unlist(Cluster_list)

gene_order = unlist(lapply(1:7, function(x) rep(paste0('MP_', x), 50)))
names(gene_order) = unlist(MP_list)
tmp = tmp[names(gene_order), names(program_order)]
tmp = as.matrix(tmp)

column_ha = HeatmapAnnotation(Metaprograms = program_order,
                              col = list(Metaprograms = colori_mp),
                              annotation_name_side = "left", 
                              show_annotation_name = FALSE, annotation_label = 'Cell type', show_legend = F)

row_ha = rowAnnotation(Metaprograms = gene_order,
                       col = list(Metaprograms = colori_mp),
                       show_annotation_name = FALSE, annotation_label = 'Cell type', show_legend = F)

h1 = Heatmap(tmp, col = ccols,
             name = "NMF score", 
             heatmap_legend_param = list(legend_direction = "horizontal"),
             #gap = unit(1.6, 'mm'),
             
             show_row_names = FALSE,
             cluster_rows = FALSE,
             row_split = gene_order, 
             cluster_row_slices = FALSE,
             left_annotation = row_ha, 
             row_title = NULL,
             
             show_column_names = FALSE,
             cluster_columns = FALSE,
             column_split = program_order,
             cluster_column_slices = FALSE,
             top_annotation = column_ha,
             column_title = NULL
)

pdf(file = 'Figures/Fig2/2b_Heatmap_NMF_score.pdf', width = 4, height = 5)
draw(h1, heatmap_legend_side = "bottom")
dev.off()

# 2C: barplot ----
load('scRNA_scATAC/RData/seu_multiome_filtered.RData')
df = seu.multiome@meta.data[seu.multiome$Metaprogram_assignment %in% paste0('MP_', 1:7), ]
table(df$Responder, df$Metaprogram_assignment)

sample_props <- df %>%
  group_by(orig.ident, Week, Responder, Metaprogram_assignment) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(orig.ident) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

table(sample_props$orig.ident, sample_props$Metaprogram_assignment)

# MP4 in NR vs R
x = sample_props$prop[sample_props$Responder %in% 'NR' & sample_props$Metaprogram_assignment %in% 'MP_4']
x = c(x, 0)
y = sample_props$prop[sample_props$Responder %in% 'R' & sample_props$Metaprogram_assignment %in% 'MP_4']
y = c(y, 0, 0, 0)
wilcox.test(x, y)

# MP5 in NR w12 vs R w12
x = sample_props$prop[sample_props$Responder %in% 'NR' & sample_props$Metaprogram_assignment %in% 'MP_5']
x = c(x, 0)
y = sample_props$prop[sample_props$Responder %in% 'R' & sample_props$Metaprogram_assignment %in% 'MP_5']
# y = c(y, 0, 0)
wilcox.test(x, y)

avg_props <- sample_props %>%
  group_by(Week, Responder, Metaprogram_assignment) %>%
  summarise(mean_prop = mean(prop), .groups = "drop")

norm_props <- avg_props %>%
  group_by(Week, Responder) %>%
  mutate(percent = 100 * mean_prop / sum(mean_prop)) %>%
  ungroup()

norm_props$Metaprogram_assignment = factor(norm_props$Metaprogram_assignment, names(colori_mp))

p1 = ggplot(norm_props, aes(x = Week, y = percent, fill = Metaprogram_assignment)) +
  geom_bar(stat = "identity", position = "stack", color = "white", width = 0.9) +
  facet_grid(~ Responder, scales = "free_x", space = "fixed") +
  scale_fill_manual(values = colori_mp, labels = c('Cell cycle', 'Melanocytic', 'EMT / Response to stress', 'Neural crest-like',
                                                   'Antigen presentation / Interferon', 'Pigmentation', 'Undefined'),
                    guide = guide_legend(nrow = 2)) +
  labs(
    x = "Week",
    y = "Cell type composition (%)",
    fill = NULL
  ) +
  scale_y_continuous(
    expand = c(0, 0),
    limits = c(0, 100.1)
  ) +
  cowplot::theme_cowplot() +
  theme(
    axis.title = element_text(size = 14, face = "plain", color = "black"),
    axis.text = element_text(size = 12, color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.line = element_line(color = "black"),
    
    strip.background = element_blank(),
    strip.text = element_text(size = 14, face = "bold", color = "black"),
    
    legend.position = "bottom",
    legend.text = element_text(size = 12),
    legend.key.size = unit(0.6, "cm"),
    
    panel.spacing = unit(0.5, "lines"),  # Less space between R and NR
    plot.margin = margin(10, 10, 10, 10)
  )
p1

ggsave(filename = 'Figures/Fig2/2c_Barplot_MP.pdf', plot = p1, device = 'pdf', width = 8, height = 8,
       dpi = 600, units = 'in', bg = 'white')

# 2D: Motif enrichment and functionality analysis ----
# Heatmap TF motif enrichment
load('scATAC/Heatmap_motif.RData', verbose = TRUE)
# Changing the TF names to gene names
identical(rownames(enrichment_mat), tf_family$`Nome TF`) # TRUE
rownames(enrichment_mat) = tf_family$`Gene name`

ccol = circlize::colorRamp2(c(0, 2, 2.5, 3), c('white', '#8bafd0', '#827ab3', '#1f1d35'))
mat = enrichment_mat

col_ha = columnAnnotation(Metaprogram = colnames(mat),
                          col = list(Metaprogram = colori_mp),
                          show_annotation_name = FALSE,
                          annotation_label = 'TF Family',
                          show_legend = FALSE)

TF_family_col = tf_family$TF.family.col
names(TF_family_col) = tf_family$`Gene name`

row_ha = rowAnnotation(TF_family = tf_family$`Gene name`, 
                       col = list(TF_family = TF_family_col), 
                       show_annotation_name = FALSE,
                       annotation_label = 'TF Family',
                       show_legend = FALSE)

h1 = Heatmap(mat, name = 'enrichment',
             # column_gap = unit(ggap, "mm"),
             col = ccol,
             cluster_rows = FALSE,
             cluster_columns = FALSE, 
             # row_split = ssplit, 
             row_title = NULL,
             left_annotation = row_ha,
             show_column_names = FALSE, 
             top_annotation = col_ha, 
             # column_title = '77 top-ranked TF in malignant cells',
             column_title_gp=grid::gpar(fontsize=20, fontface='bold')
)
h1

lgd_list = list(Legend(labels = names(families_colors), 
                       legend_gp = gpar(fill = families_colors), 
                       direction = "horizontal",
                       title = "TF Family"),
                Legend(labels = names(colori_mp), 
                       legend_gp = gpar(fill = colori_mp),
                       direction = "horizontal",
                       title = 'Metaprograms')
)


pdf(file = 'Figures/Fig2/2e_Heatmap_TF_Enrichemnt.pdf', width = 7, height = 17)
draw(h1, annotation_legend_list = lgd_list, annotation_legend_side = "bottom", heatmap_legend_side = "bottom", merge_legends = TRUE)
dev.off()

# 2E: Bubbleplot ENCODE regions ----
# Functionality analysis 
load("scATAC/RData/Results/DF.list_TFenrichment_ENCODE.RData", verbose = TRUE)

df.pELS = Reduce(rbind, DF.list.pELS)
df.pELS$cre_type = 'pELS'
df.dELS = Reduce(rbind, DF.list.dELS)
df.dELS$cre_type = 'dELS'
df.PLS = Reduce(rbind, DF.list.PLS)
df.PLS$cre_type = 'PLS'

df = rbind(df.PLS, df.pELS, df.dELS)
category = openxlsx::read.xlsx('scATAC/Category.xlsx')
category = category[!duplicated(category$ID), ]
category$Category[1:6] = 'Genome Integrity and\nCell Cycle Control'
category$Category[c(10, 41, 44, 46:48)] = 'WNT/B-Catenin\nPathway Regulation'
# category = rbind(category, c("HALLMARK_WNT_BETA_CATENIN_SIGNALING", 'Signalling and ion transport'))
# category = rbind(category, c("KEGG_WNT_SIGNALING_PATHWAY", 'Signalling and ion transport'))
# openxlsx::write.xlsx(category, file = 'scATAC/Category.xlsx')

df = df[df$ID %in% category$ID, ]
df = merge(df, category, by = 'ID', all.x = TRUE)
unique(df$Category)
df$Category[df$Category %in% 'Neuronal development'] = 'Neuronal\ndevelopment'
# df$Category[df$Category %in% 'Metabolic processes'] = 'Metabolic\nprocesses' 
df$Category[df$Category %in% 'EMT and Response to stress'] = 'EMT and\nResponse to stress' 
df$Category[df$Category %in% 'Signalling and ion transport'] = 'WNT/β-Catenin\nPathway Regulation' 
# df$Category[df$Category %in% "Cell cycle and DNA repair"] = 'Genome Integrity and\nCell Cycle Control'
# df$Category[df$Category %in% 'Epithelial development'] = 'Epithelial\ndevelopment' 
df$Category[df$Category %in% 'Immune and inflammation'] = 'Inflammatory and\nAntigen-Specific Immunity' 
df$Category[df$Category %in% 'Morphogenesis and development'] = 'Morphogenesis and\ndevelopment' 
df$Category[df$Category %in% 'Transcription and translation'] = 'Transcription and\ntranslation' 
unique(df$Category)

df$Label[df$Label %in% 'ANTIGEN PROCESSING AND PRESENTATION OF EXOGENOUS ANTIGEN'] = 'ANTIGEN PROCESSING/PRESENTATION OF EXOG. ANTIGEN'
df$Label[df$Label %in% 'ANTIGEN PROCESSING AND PRESENTATION OF ENDOGENOUS ANTIGEN'] = 'ANTIGEN PROCESSING/PRESENTATION OF ENDOG. ANTIGEN'
df$Label[df$Label %in% 'DEGRADATION OF BETA CATENIN BY THE DESTRUCTION COMPLEX'] = 'DEG. OF BETA CATENIN BY THE DESTRUCTION COMPLEX'

df = df[!df$Label %in% c('RIBOSOME', 'NEGATIVE REGULATION OF CELLULAR SENESCENCE', 'EUKARYOTIC TRANSLATION INITIATION',
                         'ESTABLISHMENT OF CELL POLARITY', 'BASEMENT MEMBRANE ORGANIZATION',
                         'REPRESSION OF WNT TARGET GENES', 'MET RECEPTOR RECYCLING'), ]

unique(df$Category)
df$Category = factor(df$Category, levels = c(unique(df$Category)[3], 
                                             unique(df$Category)[4],
                                             unique(df$Category)[5],
                                             unique(df$Category)[6],
                                             unique(df$Category)[2],
                                             unique(df$Category)[1])
)

df$size <- cut(
  df$pvalue,
  breaks = c(-Inf, 0.01, 0.05, 0.1),
  labels = c("pvalue < 0.01", "0.01 < pvalue < 0.05", "0.05 < pvalue < 0.1"),
  right = FALSE
)
table(df$size, exclude = NULL)
df$size <- factor(df$size, 
                  levels = c("pvalue < 0.01", "0.01 < pvalue < 0.05", "0.05 < pvalue < 0.1"), 
                  labels = c(1,2,3)) 

df$Metaprogram = factor(df$Metaprogram, levels = names(colori_mp), labels = c('Cell cycle', 'Melanocytic I', 'Hypoxia/\nEMT', 
                                                                              'Neural\ncrest-like', 
                                                                              'Antigen\npresentation/\nInterferon',
                                                                              'Melanocytic II', 'WNT/\nB-Catenin')
)
names(colori_mp) = labels = c('Cell cycle', 'Melanocytic I', 'Hypoxia/\nEMT', 
                              'Neural\ncrest-like', 'Antigen\npresentation/\nInterferon',
                              'Melanocytic II', 'WNT/\nB-Catenin')

df$Label[df$Label %in% "REGULATION OF MITF M DEPENDENT GENES INVOLVED IN CELL CYCLE AND PROLIFERATION"] = "REGULATION OF MITF M DEPENDENT GENES INVOLVED IN\nCELL CYCLE AND PROLIFERATION"
df = df[!df$Label %in% 'ANTIGEN PROCESSING AND PRESENTATION',]

# 12x9
ggplot(df, aes(x = cre_type, y = Label, fill = Metaprogram, size = size)) +
  geom_point(shape = 21, color = "black", stroke = 0.2, alpha = 0.7) +  # bordo nero sottile
  scale_size_manual(
    values = c(4,2,1),
    labels = c(expression("p-value" < 0.01), 
               expression("0.01 < p-value < 0.05"), 
               expression("0.05 < p-value < 0.1"))
  ) +
  scale_fill_manual(values = colori_mp) +
  facet_grid(~Category~Metaprogram, scales = 'free', space = 'free', switch = "y") +
  labs(x = NULL, y = NULL) +
  guides(
    size = guide_legend(title = expression(italic("p")*"-value"), reverse = TRUE, override.aes = list(alpha = 1)),
    fill = guide_none()
  ) +
  scale_x_discrete(expand = c(1, 1)) +  # Riduce lo spazio
  scale_y_discrete(position = "right") +
  theme_minimal(base_size = 17) +
  theme_minimal(base_size = 17) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(color = "black", size = 9),
    axis.text.x = element_text(angle = 45, hjust = 1),  
    strip.text = element_text(face = "bold", size = 9),
    strip.background = element_rect(fill = "grey90", color = NA),
    legend.text = element_text(size = 9),
    legend.title = element_text(size = 12, face = "bold"),
    legend.position = "bottom",
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0, face = "bold", size = 9),
    strip.text.x = element_text(size = 7),
    plot.margin = margin(10, 10, 10, 10),
    panel.spacing.y = unit(1, "lines"),
    panel.spacing.x = unit(0.2, "lines")
  )
ggsave('Figures/Fig2/2e_TF_Genes_EnhancerPromoter_Enrichment.pdf', width = 15, height = 10)

################################
#### Supplementary Figure 2 ####
################################
# S2A: Correlation MPs and overlap ----
# Correlation
load(file = 'scRNA/RData/malignant_subset.RData')
colnames(seu@meta.data)
seu@meta.data = seu@meta.data[, c(1:20, 40:47)]
colnames(seu@meta.data) = gsub('Cluster', 'MP_', colnames(seu@meta.data))

cor_result <- cor(seu@meta.data[, paste0('MP_', 1:7)], method = "pearson")
cor_result_table = cor_result
write.table(cor_result_table, file = 'MP_correlation.txt')
# print(paste("Pearson correlation:", round(cor_result, 3)))
colnames(cor_result)

# heatmap
cor_result = melt(cor_result)

p = ggplot(cor_result, aes(x=Var2, y=Var1, fill=value)) +
  geom_tile(color="white") +
  # geom_text(aes(label=round(value, 3)), color="black") +
  scale_fill_gradient2(low = "#6AACD0", mid = "white", high="#E58267", limits = c(-1,1), breaks = c(-1, 0, 1)) +
  labs(x="", y="", fill="Correlation") +
  theme_minimal() +
  theme(axis.text = element_text(size = 10), 
        axis.text.x = element_text(hjust = 1, angle = 30),
        panel.grid = element_blank(),
        text = element_text(size = 15, face = 'bold'),
        legend.text = element_text(size = 10), 
        legend.title = element_text(size = 10, face = 'bold') #,
        # legend.direction = 'horizontal', legend.position = 'bottom', legend.title.position = 'top'
  )
p

ggsave(filename = 'Figures/Fig2/Supp_2a_MP_correlation.pdf', plot = p, device = 'pdf', width = 5.5, height = 4,
       dpi = 600, units = 'in', bg = 'white')

# overlap genes between mps
MP = as.data.frame(MP_list)
overlap = apply(MP , 2, function(x) apply(MP , 2, function(y) length(intersect(x,y)))) 

overlap <- melt(overlap)
overlap$Label = overlap$value
overlap$Label[overlap$Label == 0] = ''
# overlap$Label[overlap$Label == 50] = ''

p = ggplot(overlap, aes(x=Var2, y=Var1, fill=value)) +
  geom_tile(color="white") +
  geom_text(aes(label=Label), color="black") +
  scale_fill_gradient(low="white", high="#006766", limits = c(0,50), breaks = c(0, 25, 50)) +
  labs(x="", y="", fill="Gene overlap") +
  theme_minimal() +
  theme(axis.text = element_text(size = 10), 
        axis.text.x = element_text(hjust = 1, angle = 30),
        panel.grid = element_blank(),
        text = element_text(size = 15, face = 'bold'),
        legend.text = element_text(size = 10), 
        legend.title = element_text(size = 10, face = 'bold') #,
        # legend.direction = 'horizontal', legend.position = 'bottom', legend.title.position = 'top'
  )
ggsave(filename = 'Figures/Fig2/Supp_2a_Gene_overlap.pdf', plot = p, device = 'pdf', width = 5.5, height = 4,
       dpi = 600, units = 'in', bg = 'white')

# S2B: Heatmap AddModuleScores ----
load('scRNA/RData/malignant_subset.RData')
colnames(seu@meta.data)
df = seu@meta.data[, c('Metaprogram_assignment', paste0('Cluster', 1:7))]
rm(seu)

colnames(df)[2:ncol(df)] = paste0('MP_', 1:7)
order_cells = lapply(names(colori_mp), function(x) {
  tmp = df[df$Metaprogram_assignment %in% x, ]
  tmp = tmp[order(tmp[, x], decreasing = TRUE), ]
  rownames(tmp)
})

mat = df[Reduce(c, order_cells), 2:8]
range(mat)

ccols = colorRamp2(c(-1, 0, 1), c('blue', 'grey90', 'red'))

column_ha = HeatmapAnnotation(Metaprograms = names(colori_mp),
                              col = list(Metaprograms = colori_mp),
                              annotation_name_side = "left", 
                              show_annotation_name = FALSE, annotation_label = 'Cell type', show_legend = F)

cell_order = lapply(1:7, function(x) rep(paste0('MP_', x), length(order_cells[[x]])) )
cell_order = Reduce(c, cell_order)
names(cell_order) = Reduce(order_cells)

row_ha = rowAnnotation(Metaprograms = cell_order,
                       col = list(Metaprograms = colori_mp),
                       show_annotation_name = FALSE, annotation_label = 'Cell type', show_legend = F)

h1 = Heatmap(as.matrix(mat), col = ccols,
             name = "MP score", 
             heatmap_legend_param = list(legend_direction = "horizontal", at = seq(-1, 1)),
             
             show_row_names = FALSE,
             cluster_rows = FALSE,
             row_split = cell_order,
             cluster_row_slices = FALSE,
             left_annotation = row_ha,
             row_title = NULL,
             row_gap = unit(0.4, 'mm'),
             
             # column_gap = unit(1, 'mm'),
             # column_split = names(colori_mp),
             show_column_names = FALSE,
             cluster_columns = FALSE,
             cluster_column_slices = FALSE,
             top_annotation = column_ha,
             column_title = NULL
)

pdf(file = 'Figures/Fig2/Supp_2b_Heatmap_AddModuleScore.pdf', width = 3, height = 5)
draw(h1, heatmap_legend_side = "bottom")
dev.off()

# S2C: Epigenetic plasticity ----
# Barplot correct prediction
ATAC_correct_per_class <- predicted.labels %>%
  dplyr::select(True_label, predicted.id) %>%
  dplyr::mutate(Correct = (True_label == predicted.id))

ATAC_correct_per_class

df = ATAC_correct_per_class %>%
  group_by(True_label) %>%
  summarise(perc_correct = mean(Correct) * 100) %>% 
  as.data.frame()
df$True_label = factor(df$True_label, names(celltype))
p2 = ggplot(df, aes(x = True_label, y = perc_correct, fill = True_label)) +
  geom_bar(stat = "identity", position = "stack", color = "white", width = 0.9) +
  scale_fill_manual(values = celltype) +
  labs(
    x = "Cell types",
    y = "Correct prediction (%)",
    fill = NULL
  ) +
  cowplot::theme_cowplot() +
  theme(
    axis.title = element_text(size = 14, face = "plain", color = "black"),
    axis.text = element_text(size = 12, color = "black"),
    axis.text.x = element_text(size = 12, color = "black", angle = 30, hjust = 1),
    
    axis.ticks = element_line(color = "black"),
    axis.line = element_line(color = "black"),
    
    strip.background = element_blank(),
    strip.text = element_text(size = 14, face = "bold", color = "black"),
    
    legend.position = "none",
    # legend.text = element_text(size = 12),
    # legend.key.size = unit(0.6, "cm"),
    
    panel.spacing = unit(0.5, "lines"),  # Less space between R and NR
    plot.margin = margin(10, 10, 10, 10)
  )

ggsave(filename = 'Figures/Supp_2C_barplot_epigeneticplasticity.pdf', plot = p2, device = 'pdf', width = 8, height = 5,
       dpi = 600, units = 'in', bg = 'white')

# Heatmap percentage of prediction score
score_cols <- grep("^prediction.score", colnames(predicted.labels), value = TRUE)
score_cols = score_cols[1:13]

df_mean <- predicted.labels %>%
  group_by(True_label) %>%
  summarise(across(all_of(score_cols), mean, .names = "{.col}"))

colnames(df_mean) <- c("True_label", str_replace(colnames(df_mean)[-1], "prediction.score.", ""))
mat <- as.matrix(df_mean[,-1])
colnames(mat) = gsub('mean_', '', colnames(mat))
colnames(mat) = gsub('\\.', ' ', colnames(mat))
colnames(mat) = gsub('B Plasma cells', 'B/Plasma cells', colnames(mat))

rownames(mat) <- df_mean$True_label

setdiff(names(celltype),colnames(mat))
mat = mat[rev(names(celltype)), names(celltype)]

col_fun <- colorRamp2(c(0, 0.5, 1), c("#f7f7f7", "#35b779", "#31688e"))

p3 = Heatmap(
  mat,
  name = "Prob.",
  col = col_fun,
  cluster_rows = FALSE, cluster_columns = FALSE,
  row_title = "Predicted cell type probability (ATAC)",
  column_title = "Cell type annotation (RNA)",
  heatmap_legend_param = list(title = "Probability"), column_names_rot = 30
)

pdf(file = 'Figures/Supp_2C_PrbabilityHeatmap_epigeneticplasticity.pdf', width = 8, height = 6)
draw(p3)
dev.off()

# S2D: Correlation with signatures Tsoi, Tirosh, Pozniak, Baron, Hoek, Rambow, Varfaillie, Wouters, Soldatov ----
load('Supp_tables/Signatures_MP_correlation.RData')
load("/home/caruso/Analisi2021/Maio_EPICA/dati_clinici/Anichini_signatures/Melanoma_Anichini_Signatures.RData", verbose = TRUE)
rm(AntiPD1_response, EMT_and_MELANOMA_diff_MARKERS, Pozniak_immune)
Tirosh_signatures$TIR_T = NULL
Tirosh_signatures$TIR_B = NULL
Tirosh_signatures$TIR_MACRO = NULL
Tirosh_signatures$TIR_ENDO = NULL
Tirosh_signatures$TIR_CAF = NULL
load("/home/ciervo/EPICA/NIBIT/scRNA/Signatures_melanoma_Pozniak_100.RData", verbose = T)

colnames(cor_result)
mat = cor_result[1:7, 8:ncol(cor_result)]
mp_order = c('MP_1', 'MP_5', 'MP_2', 'MP_6', 'MP_7', 'MP_4', 'MP_3')

mat = mat[mp_order, ]
# ssplit = c(rep('Tsoi et al.', length(TSOI_differentiation_Stage)), 
#            rep('Tirosh et al.', length(Tirosh_signatures)),
#            rep('Pozniak et al.', length(signature)))
# names(ssplit) = c(names(TSOI_differentiation_Stage), names(Tirosh_signatures), names(signature))

# heatmap
color_heatmap = colorRamp2(c(-1, -.5, -.01, 0, .01, .5, 1), c('#053061', '#6AACD0', '#FFFFFF', '#FFFFFF', '#FFFFFF', '#E58267', '#67001F') ) 
row_ha = rowAnnotation(Metaprogram = mp_order, 
                       col = list(Metaprogram = colori_mp), 
                       na_col = "white", show_annotation_name = FALSE,
                       annotation_name_gp = gpar(fontsize = 10)
)
range(mat)
colnames(mat)

# Correlation between Baron and MP
mat_tmp = mat[, grep("Baron", colnames(mat))]
mat_tmp = mat_tmp[, c("Baron_Mature_melanocytes", "Baron_Stress-like", "Baron_Neural_crest")]
h1 <- Heatmap(mat_tmp, 
              col = color_heatmap,
              name = 'Correlation',
              cluster_rows = FALSE, 
              show_row_names = FALSE,
              cluster_columns = FALSE, 
              show_column_names = TRUE, 
              left_annotation = row_ha,
              column_names_rot = 45, 
              column_title = "Baron et al.",
              column_title_gp = gpar(fontface = 'bold', fontsize = 15), 
              column_names_gp = gpar(fontsize = 12), 
              column_labels = c("Mature melanocytic", "Stress-like", "Neural crest-like"),
              show_column_dend = FALSE
)
h1

# Correlation between Hoek and MP
mat_tmp = mat[, grep("Hoek", colnames(mat))]
mat_tmp
h2 <- Heatmap(mat_tmp, 
              col = color_heatmap,
              name = 'Correlation',
              cluster_rows = FALSE, 
              show_row_names = FALSE,
              cluster_columns = FALSE, 
              show_column_names = TRUE,
              left_annotation = row_ha,
              column_names_rot = 45,
              column_title = "Hoek et al.",
              column_title_gp = gpar(fontface = 'bold', fontsize = 15), 
              column_names_gp = gpar(fontsize = 12), 
              column_labels = c("Proliferative", "Invasive"),
              show_column_dend = FALSE
)
h2

colnames(mat)

# Correlation between Rambow and MP
mat_tmp = mat[, grep("Rambow", colnames(mat))]
mat_tmp
h3 <- Heatmap(mat_tmp, 
              col = color_heatmap,
              name = 'Correlation',
              cluster_rows = FALSE, 
              show_row_names = FALSE,
              cluster_columns = FALSE, 
              show_column_names = TRUE,
              left_annotation = row_ha,
              column_names_rot = 45, 
              column_title = "Rambow et al.",
              column_title_gp = gpar(fontface = 'bold', fontsize = 15),
              column_names_gp = gpar(fontsize = 12), 
              column_labels = c('Pigmentation', 'SMC', 'NCSC', 'Invasive')
              
)
h3

colnames(mat)

# Correlation between Verfaillie and MP
mat_tmp = mat[, grep("Verfaillie", colnames(mat))]
mat_tmp = mat_tmp[, c("Verfaillie_Proliferative", "Verfaillie_Invasive")]
h4 <- Heatmap(mat_tmp, 
              col = color_heatmap,
              name = 'Correlation',
              cluster_rows = FALSE, 
              show_row_names = FALSE,
              cluster_columns = FALSE, 
              show_column_names = TRUE,
              left_annotation = row_ha,
              column_names_rot = 45, 
              column_title = "Verfaillie et al.",
              column_title_gp = gpar(fontface = 'bold', fontsize = 15),
              column_names_gp = gpar(fontsize = 12), 
              column_labels = c('Proliferative', 'Invasive')
              
)
h4
colnames(mat)

# Correlation between Wouters and MP
mat_tmp = mat[, grep("Wouters", colnames(mat))]
mat_tmp = mat_tmp[, c("Wouters_Melanocytic", "Wouters_Intermediate", "Wouters_Mesenchymal-like")]
h5 <- Heatmap(mat_tmp, 
              col = color_heatmap,
              name = 'Correlation',
              cluster_rows = FALSE, 
              show_row_names = FALSE,
              cluster_columns = FALSE, 
              show_column_names = TRUE,
              left_annotation = row_ha,
              column_names_rot = 45, 
              column_title = "Wouters et al.",
              column_title_gp = gpar(fontface = 'bold', fontsize = 15),
              column_names_gp = gpar(fontsize = 12), 
              column_labels = c("Melanocytic", "Intermediate", "Mesenchymal-like")
              
)
h5

# heatmap
load("Supp_tables/Signatures_MP_correlation.RData")
load("/home/caruso/Analisi2021/Maio_EPICA/dati_clinici/Anichini_signatures/Melanoma_Anichini_Signatures.RData", verbose = TRUE)
rm(AntiPD1_response, EMT_and_MELANOMA_diff_MARKERS, Pozniak_immune)
Tirosh_signatures$TIR_T = NULL
Tirosh_signatures$TIR_B = NULL
Tirosh_signatures$TIR_MACRO = NULL
Tirosh_signatures$TIR_ENDO = NULL
Tirosh_signatures$TIR_CAF = NULL
load("/home/ciervo/EPICA/NIBIT/scRNA/Signatures_melanoma_Pozniak_100.RData", verbose = T)

colnames(cor_result)
mat = cor_result[1:7, 8:ncol(cor_result)]
mp_order = c('MP_1', 'MP_5', 'MP_2', 'MP_6', 'MP_7', 'MP_4', 'MP_3')

mat = mat[mp_order, ]
# color_heatmap = colorRamp2(c(-1, -.5, -.1, 0, .1, .5, 1), c('#053061', '#6AACD0', '#FFFFFF', '#FFFFFF', '#FFFFFF', '#E58267', '#67001F') )
row_ha = rowAnnotation(Metaprogram = mp_order, 
                       col = list(Metaprogram = colori_mp), 
                       na_col = "white", show_annotation_name = FALSE,
                       annotation_name_gp = gpar(fontsize = 10)
)
range(mat)
colnames(mat)
mat_tmp = mat[, names(signature)]
h6 <- Heatmap(mat_tmp, 
              col = color_heatmap,
              name = 'Correlation',
              cluster_rows = FALSE, 
              show_row_names = FALSE,
              cluster_columns = TRUE, 
              show_column_names = TRUE, 
              left_annotation = row_ha,
              column_names_rot = 45, 
              column_title = "Pozniak et al.",
              column_title_gp = gpar(fontface = 'bold', fontsize = 15), 
              column_names_gp = gpar(fontsize = 12), 
              column_labels = c('Neural crest-like', 'INF-a/b response', 'Melanocytic',
                                'Mesenchymal', 'Mitochondrial', 'Antigen presentation',
                                'Mitotic', 'Stress (p53)', 'Stress (hypoxia)'),
              show_column_dend = FALSE
)

mat_tmp = mat[, rev(names(TSOI_differentiation_Stage))]

# Correlation between Tsoi and MP
h7 <- Heatmap(mat_tmp, 
              col = color_heatmap,
              name = 'Correlation',
              cluster_rows = FALSE, 
              show_row_names = FALSE,
              cluster_columns = FALSE, 
              show_column_names = TRUE,
              left_annotation = row_ha,
              column_names_rot = 45,
              column_title = "Tsoi et al.",
              column_title_gp = gpar(fontface = 'bold', fontsize = 15), 
              column_names_gp = gpar(fontsize = 12), 
              column_labels = c('Melanocytic', 'Melanocytic\nTransitory',
                                'Transitory', 'Transitory\nNeural crest-like',
                                'Neural crest-like', 'Neural crest-like\nUndifferentiated',
                                'Undifferentiated'),
              show_column_dend = FALSE
)

# Correlation between Tirosh and MP
mat_tmp = mat[, c('TIR_MELA_CCYCLE', 'TIR_MITF', 'TIR_MELA', 'TIR_AXL')]

h8 <- Heatmap(mat_tmp, 
              col = color_heatmap,
              name = 'Correlation',
              cluster_rows = FALSE, 
              show_row_names = FALSE,
              cluster_columns = FALSE, 
              show_column_names = TRUE,
              left_annotation = row_ha,
              column_names_rot = 45, 
              column_title = "Tirosh et al.",
              column_title_gp = gpar(fontface = 'bold', fontsize = 15),
              column_names_gp = gpar(fontsize = 12), 
              column_labels = c('Melanoma\ncell cycle', 'MITF-program', 'Melanoma', 'AXL-program')
              
)

# Soldotov
load('Revisioni/Correlation_signatures_Soldatov.RData')
colnames(cor_result)
mat = cor_result[1:7, 8:ncol(cor_result)]
mp_order = c('MP_1', 'MP_5', 'MP_2', 'MP_6', 'MP_7', 'MP_4', 'MP_3')

mat = mat[mp_order, ]
# ssplit = c(rep('Tsoi et al.', length(TSOI_differentiation_Stage)), 
#            rep('Tirosh et al.', length(Tirosh_signatures)),
#            rep('Pozniak et al.', length(signature)))
# names(ssplit) = c(names(TSOI_differentiation_Stage), names(Tirosh_signatures), names(signature))

# heatmap
color_heatmap = colorRamp2(c(-1, -.5, -.01, 0, .01, .5, 1), c('#053061', '#6AACD0', '#FFFFFF', '#FFFFFF', '#FFFFFF', '#E58267', '#67001F') ) 
row_ha = rowAnnotation(Metaprogram = mp_order, 
                       col = list(Metaprogram = colori_mp), 
                       na_col = "white", show_annotation_name = FALSE,
                       annotation_name_gp = gpar(fontsize = 10)
)
range(mat)
colnames(mat)

# Correlation between Baron and MP
h9 <- Heatmap(mat, 
              col = color_heatmap,
              name = 'Correlation',
              cluster_rows = FALSE, 
              show_row_names = FALSE,
              cluster_columns = FALSE, 
              show_column_names = TRUE, 
              left_annotation = row_ha,
              column_names_rot = 45, 
              column_title = "Soldotov et al.",
              column_title_gp = gpar(fontface = 'bold', fontsize = 15), 
              column_names_gp = gpar(fontsize = 12), 
              column_labels = c('Cranial Neural Crest', 'Trunk Neural Crest', 'Pre delamination',
                                'Delamination', 'Mesenchymal', 'Melanoblast'),
              show_column_dend = FALSE
)


pdf(file = 'Figures/Fig2/Supp_2c_Correlation_heatmap_signatures.pdf', width = 22, height = 4)
h7+h8+h6+h1+h2+h3+h4+h5+h9
dev.off()

# S2E: Stemness ----
load("scRNA/RData/malignant_subset.RData")
df_cell <- seu@meta.data[, c("Metaprogram_assignment", "StemSig_top100", "orig.ident")] %>%
  as.data.frame() %>%
  dplyr::rename(sample_id = orig.ident)
df_cell[1:5, ]
df_cell$Metaprogram_assignment <- factor(
  df_cell$Metaprogram_assignment,
  levels = names(colori_mp),
  labels = c("Cell Cycle", "Melanocytic I", "Hypoxia/EMT",
             "Neural crest-like", "Antigen presentation/\nInterferon",
             "Melanocytic II", "Wnt/B-catenin")
)
df_cell[1:5, ]
colori_mp_nn <- colori_mp
names(colori_mp_nn) <- levels(df_cell$Metaprogram_assignment)

my_comp <- list(
  "MP1MP5" = c(names(colori_mp_nn)[1], names(colori_mp_nn)[5]),
  "MP2MP5" = c(names(colori_mp_nn)[2], names(colori_mp_nn)[5]),
  "MP3MP5" = c(names(colori_mp_nn)[3], names(colori_mp_nn)[5]),
  "MP4MP5" = c(names(colori_mp_nn)[4], names(colori_mp_nn)[5]),
  "MP6MP5" = c(names(colori_mp_nn)[6], names(colori_mp_nn)[5]),
  "MP7MP5" = c(names(colori_mp_nn)[7], names(colori_mp_nn)[5]),
  "MP1MP4" = c(names(colori_mp_nn)[1], names(colori_mp_nn)[4]),
  "MP2MP4" = c(names(colori_mp_nn)[2], names(colori_mp_nn)[4]),
  "MP3MP4" = c(names(colori_mp_nn)[3], names(colori_mp_nn)[4]),
  "MP6MP4" = c(names(colori_mp_nn)[6], names(colori_mp_nn)[4]),
  "MP7MP4" = c(names(colori_mp_nn)[7], names(colori_mp_nn)[4])
)

pb = df_cell %>%
  group_by(sample_id, Metaprogram_assignment) %>%
  summarise(StemSig_top100 = median(StemSig_top100, na.rm = TRUE), .groups = "drop") %>% 
  as.data.frame()

res_list = vector("list", length(my_comp))
nm_vec = names(my_comp)

for (i in 1:length(my_comp)) {
  
  comp = my_comp[[i]]
  nm   = nm_vec[i]
  
  tmp = pb %>%
    dplyr::filter(Metaprogram_assignment %in% comp) %>%
    dplyr::mutate(Metaprogram_assignment = droplevels(Metaprogram_assignment))
  
  n1 = sum(tmp$Metaprogram_assignment == comp[1])
  n2 = sum(tmp$Metaprogram_assignment == comp[2])
  
  pval = wilcox.test(StemSig_top100 ~ Metaprogram_assignment, data = tmp)$p.value
  
  res_list[[i]] = data.frame(
    comparison = nm,
    group1 = comp[1],
    group2 = comp[2],
    p = pval,
    stringsAsFactors = FALSE
  )
}

res_tests = dplyr::bind_rows(res_list)

y_max = max(pb$StemSig_top100, na.rm = TRUE)
step_abs = 0.06 * diff(range(pb$StemSig_top100, na.rm = TRUE))

p_anno = res_tests %>%
  mutate(
    y.position = y_max + row_number() * step_abs,
    label = dplyr::case_when(
      is.na(p)   ~ NA_character_,
      p < 0.001  ~ "***",
      p < 0.01   ~ "**",
      p < 0.05   ~ "*",
      TRUE       ~ "ns"
    )
  )

p <- ggplot(pb, aes(x = Metaprogram_assignment, y = StemSig_top100, fill = Metaprogram_assignment)) +
  geom_boxplot(width = 0.3, alpha = 0.6, outliers = TRUE) +
  scale_fill_manual(values = colori_mp_nn) +
  ggpubr::stat_pvalue_manual(
    p_anno,
    xmin = "group1", xmax = "group2",
    y.position = "y.position",
    label = "label",
    tip.length = 0.01,
    size = 4,
    inherit.aes = FALSE
  ) +
  coord_cartesian(clip = "off") +   # <<< CHIAVE
  labs(x = "", y = "Stemness Module Score") +
  theme_minimal(base_size = 11) +
  theme(
    plot.margin = margin(t = 20, r = 5, b = 5, l = 5),  # spazio sopra
    axis.text.x = element_text(size = 9, face = "bold", angle = 45, hjust = 1),
    axis.text.y = element_text(size = 9),
    axis.title.y = element_text(size = 10, face = "bold"),
    axis.line = element_line(colour = "black"),
    panel.grid = element_blank(),
    legend.position = "none"
  )
p

ggsave(filename = 'Revisioni/boxplot_stemness.pdf', plot = p, device = 'pdf', height = 5, width = 6,
       dpi = 600, units = 'in', bg = 'white')

# S2F-G: CellChat ----
# R
load("scRNA/CellChat/cellChat_object_Responder.RData")
groupSize <- as.numeric(table(cellchat@idents))
colors_cellchat = c(colors_tme, colori_mp)
colors_cellchat = colors_cellchat[levels(cellchat@idents)]
levels(cellchat@idents)

# MHC-I e MHC-II
sapply(c('MHC-I', 'MHC-II'), function(i) {
  par(mfrow = c(1, 1), xpd = TRUE, cex = 1, mar = c(0, 1, 0, 1))
  options(repr.plot.width = 4, repr.plot.height = 4)
  pdf(file = paste0('Figures/Fig2/Supp_2e_CellChat_', i, '_Responder.pdf'), width = 4, height = 6)
  netVisual_aggregate(cellchat, signaling = i, arrow.size = 0.5,
                      layout = "circle", vertex.label.cex = 0.7,
                      color.use = colors_cellchat, remove.isolate = TRUE)
  dev.off()
}
)

# NR
load("scRNA/CellChat/cellChat_object_NonResponder.RData")
groupSize <- as.numeric(table(cellchat@idents))
colors_cellchat = c(colors_tme, colori_mp)
colors_cellchat = colors_cellchat[levels(cellchat@idents)]
levels(cellchat@idents)

pairLR <- extractEnrichedLR(cellchat, signaling = 'PTN', geneLR.return = FALSE)
LR.show <- pairLR$interaction_name # show one ligand-receptor pair
LR.show = LR.show[c(1,3)]

sapply(LR.show, function(i) {
  pdf(file = paste0('Figures/Fig2/Supp_2d_CellChat_', i, '_NonResponder.pdf'), width = 4, height = 6)
  netVisual_individual(cellchat, signaling = 'PTN',  pairLR.use = i, 
                       color.use = colors_cellchat, remove.isolate = TRUE, arrow.size = 0.5)
  dev.off()
}
)

################################
#### Supplementary Figure 3 ####
################################
# S3A: Footprinting ----
source('/home/ciervo/Functions/scATAC_function.R')
load('scATAC/RData/malignant_subset_motifs.RData')
load('scATAC/TFOI_MPs.RData')

tf_vector = c('E2F3', 'KLF3', 'MITF', 'KMT2B', 'FOS', 'JUN', 'NFAC2', 'SP1', 'IRF1', 'STAT1', 'ARNT2', 'MAX', 'LEF1', 'TF7L2')

motif.name = Reduce(rbind, TF.of.interest)
motif.name = motif.name[motif.name$TF.name %in% tf_vector, ]
motif.name = motif.name$motif %>% unique
motif.name = c(motif.name, c('MYC.H13CORE.0.P.B', 'AP2A.H13CORE.0.PSM.A', 'PPARG.H13CORE.0.P.B', 'NFAC2.H13CORE.0.P.B', 'E2F1.H13CORE.0.P.B'))

p1 = PlotFootprint_mod(malignant, features = c('E2F1.H13CORE.0.P.B', 'EGR1.H13CORE.0.PS.A' #MP1
                                               # 'MITF.H13CORE.0.P.B', 'AP2A.H13CORE.0.PSM.A' #MP2
                                               # 'FOS.H13CORE.0.P.B', 'ATF3.H13CORE.0.P.B' #MP3
                                               # 'NFAC2.H13CORE.0.P.B', 'ZEB1.H13CORE.0.P.B' #MP4
                                               # 'IRF1.H13CORE.0.P.B', 'STAT1.H13CORE.1.P.B' #MP5
                                               # 'MAX.H13CORE.2.S.C', 'MYC.H13CORE.0.P.B' #MP6
                                               # 'LEF1.H13CORE.0.PSM.A', 'TF7L2.H13CORE.0.P.B' #MP7
),
group.by = 'Metaprogram_assignment', group.colors = colori_mp, label.top = 0,
plot.title = c('E2F1', 'EGR1'
               # 'MITF', 'TFAP2A'
               # 'FOS', 'ATF3'
               # 'NFATC2', 'ZEB1'
               # 'IRF1', 'STAT1'
               # 'MAX', 'MYC'
               # 'LEF1', 'TCF7L2'
),
priority.group = 'MP_1'
)
p2 = MotifPlot(
  object = malignant,
  motifs = c( 'E2F1.H13CORE.0.P.B', 'EGR1.H13CORE.0.PS.A' ), 
  use.names = FALSE
)

pdf(file = 'Figures/Fig2/Supp_3b_MP1_Footprinting.pdf', width = 8, height = 5)
cowplot::plot_grid(p1, p2, nrow = 2, rel_heights = c(0.8, 0.5))
dev.off()

# S3B: SCENIC results ----
# Genes AUC plot
load('scRNA_scATAC/SCENICplus/downstream_scenicplus.RData')
RSS_scores <- as.matrix(read_csv("scRNA_scATAC/SCENICplus/CSV_output/RSS_scores.csv"))
rownames(RSS_scores) = c('MP_2', 'MP_7', 'MP_5', 'MP_1', 'MP_4', 'MP_6', "MP_3")
RSS_scores = RSS_scores[paste0('MP_', 1:7), ]
colnames(RSS_scores)
eRegulon_metadata <- as.data.frame(read_csv("scRNA_scATAC/SCENICplus/CSV_output/eRegulon_metadata.csv"))
eRegulon_metadata = eRegulon_metadata[eRegulon_metadata$regulation == 1, ]

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
regulon_order$eRegulon_name_TF = gsub(pattern = '\\(\\+\\)', replacement = '', regulon_order$eRegulon_name)

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

# overlap 
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

order_tf = regulon_order$eRegulon_name_TF
# order_tf = gsub('\\(\\+\\)', '', order_tf)

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

p = (p2/p1/p3 + plot_layout(heights = c(1, 3, 3))) & coord_flip()

ggsave(filename = 'Figures/Fig2/Supp_3a_SCENIC_Heatmap.pdf', plot = p, device = 'pdf', width = 20, height = 10, 
       dpi = 600, units = 'in', bg = 'white')

##################
#### Figure 6 ####
##################
# dir.create('Figures/Fig6')
# 6A: Bulk Melanoma (ICI-response) ----
load('Melanoma_Bulk/dataframe_bulk_all.RData')
summary(mmat$purity)
mmat = mmat[mmat$purity >= 0.3, ]
# mmat$Treatment[mmat$Dataset %in% c('vanallen', 'amato')] = 'PRE'
table(mmat$Dataset, mmat$Therapy, exclude = NULL)
mmat$Treatment[mmat$Treatment %in% 'EDT'] = 'POST'
mmat$Treatment[mmat$Treatment %in% 'ON'] = 'POST'
mmat$Treatment[is.na(mmat$Treatment)] = 'Unknown'
mmat$Treatment[mmat$Treatment %in% ''] = 'Unknown'
table(mmat$Treatment, mmat$Dataset)
colnames(mmat)
mmat_treatment = reshape2::melt(mmat[, c(1:8, 10)])
mmat_treatment$response_NR = factor(mmat_treatment$response_NR, levels = c("R", "N"))
table(mmat$Treatment, mmat$Dataset, exclude = NULL)
table(mmat$Therapy, mmat$Dataset, exclude = NULL)

mmat_treatment = reshape2::melt(mmat[, c(1:8, 10)])
table(mmat_treatment$Treatment, exclude = NULL)
mmat_treatment = mmat_treatment[!mmat_treatment$Treatment %in% 'Unknown', ]

mmat_treatment$response_NR = factor(mmat_treatment$response_NR, levels = c("R", "N"), labels = c('R', 'NR'))
mmat_treatment$Treatment = factor(mmat_treatment$Treatment, levels = c("PRE", "POST"))
table(mmat_treatment$response_NR, mmat_treatment$Treatment)
#     PRE   POST
# R   854   210
# NR  1134  616

p = ggplot(data = mmat_treatment, aes(x = response_NR, y = value, fill = response_NR)) + 
  geom_boxplot(width = 0.3, alpha = 0.6, position = position_dodge(width = 1), outliers = FALSE) +
  theme_minimal(base_size = 11) +
  theme(
    strip.background = element_rect(fill = "#f0f0f0", color = NA),
    strip.text = element_text(size = 10, face = "bold", colour = 'black'),
    axis.text.x = element_text(size = 9, face = "bold"),
    axis.text.y = element_text(size = 9),
    axis.title.y = element_text(size = 10, face = "bold"),
    axis.line = element_line(colour = 'black'),
    panel.grid = element_blank(),
    legend.position = "none",
    panel.spacing.x = unit(0.2, "lines"),
    panel.spacing.y = unit(0.5, "lines")
    ) + 
  ylim(c(0, 0.5)) +
  facet_grid(~Treatment~variable, scales = 'free_x', axes = "all_x") +
  stat_compare_means(
    method = "wilcox.test", 
    label.x.npc = 'center', 
    label.y = 0.45, 
    hide.ns = TRUE
    ) +
  scale_fill_manual(
    values = c('NR' = '#5F9EA0', 'R' = '#FA8072')
    ) +
  labs(
    x = "",
    y = "Cell fraction",
    fill = "Response"
    )

ggsave(filename = 'Figures/Fig6/6a_Bulk_ICI.pdf', plot = p, device = 'pdf', width = 9, height = 5,
       dpi = 600, units = 'in', bg = 'white')

# 6B: Forestplot MPs ----
# PRE
meta = mmat[mmat$Treatment %in% c('PRE'), ]
meta = meta[!is.na(meta$`vital status`), ]
meta = meta[!meta$`vital status` %in% '', ]
sum(is.na(meta$`overall survival (days)`))
meta = meta[!is.na(meta$`overall survival (days)`), ]
# View(table(meta$patient_name, meta$Dataset))

meta$`OS months` = meta$`overall survival (days)` / 30.44

table(meta$`vital status`, exclude = NULL)

meta$Sstate[meta$`vital status` %in% c('0:LIVING', 'Alive')] = 0
meta$Sstate[meta$`vital status` %in% c('1:DECEASED', 'Dead')] = 1
table(meta$`vital status`, meta$Sstate, exclude = NULL)

meta$`OS months`
meta$MP1_hl = ifelse(meta$MP_1 >= median(meta$MP_1), 'High', 'Low')
meta$MP2_hl = ifelse(meta$MP_2 >= median(meta$MP_2), 'High', 'Low')
meta$MP3_hl = ifelse(meta$MP_3 >= median(meta$MP_3), 'High', 'Low')
meta$MP4_hl = ifelse(meta$MP_4 >= median(meta$MP_4), 'High', 'Low')
meta$MP5_hl = ifelse(meta$MP_5 >= median(meta$MP_5), 'High', 'Low')
meta$MP6_hl = ifelse(meta$MP_6 >= median(meta$MP_6), 'High', 'Low')
meta$MP7_hl = ifelse(meta$MP_7 >= median(meta$MP_7), 'High', 'Low')

meta$MP1_hl = factor(meta$MP1_hl, levels = c("High", "Low"))
meta$MP2_hl = factor(meta$MP2_hl, levels = c("High", "Low"))
meta$MP3_hl = factor(meta$MP3_hl, levels = c("High", "Low"))
meta$MP4_hl = factor(meta$MP4_hl, levels = c("High", "Low"))
meta$MP5_hl = factor(meta$MP5_hl, levels = c("High", "Low"))
meta$MP6_hl = factor(meta$MP6_hl, levels = c("High", "Low"))
meta$MP7_hl = factor(meta$MP7_hl, levels = c("High", "Low"))

cox_model <- coxph(Surv(`OS months`, Sstate) ~ MP1_hl+MP2_hl+MP3_hl+MP4_hl+MP5_hl+MP6_hl+MP7_hl, data = meta)
summary(cox_model)
p1 = ggforest(cox_model, data = meta, main = 'Hazard ratio - Pre treatment')
p1

# POST
meta = mmat[mmat$Treatment %in% c('POST'), ]
meta = meta[!is.na(meta$`vital status`), ]
meta = meta[!meta$`vital status` %in% '', ]
sum(is.na(meta$`overall survival (days)`))
meta = meta[!is.na(meta$`overall survival (days)`), ]
# View(table(meta$patient_name, meta$Dataset))

meta$`OS months` = meta$`overall survival (days)` / 30.44

table(meta$`vital status`, exclude = NULL)

meta$Sstate[meta$`vital status` %in% c('0:LIVING', 'Alive')] = 0
meta$Sstate[meta$`vital status` %in% c('1:DECEASED', 'Dead')] = 1
table(meta$`vital status`, meta$Sstate, exclude = NULL)

meta$`OS months`
meta$MP1_hl = ifelse(meta$MP_1 >= median(meta$MP_1), 'High', 'Low')
meta$MP2_hl = ifelse(meta$MP_2 >= median(meta$MP_2), 'High', 'Low')
meta$MP3_hl = ifelse(meta$MP_3 >= median(meta$MP_3), 'High', 'Low')
meta$MP4_hl = ifelse(meta$MP_4 >= median(meta$MP_4), 'High', 'Low')
meta$MP5_hl = ifelse(meta$MP_5 >= median(meta$MP_5), 'High', 'Low')
meta$MP6_hl = ifelse(meta$MP_6 >= median(meta$MP_6), 'High', 'Low')
meta$MP7_hl = ifelse(meta$MP_7 >= median(meta$MP_7), 'High', 'Low')

meta$MP1_hl = factor(meta$MP1_hl, levels = c("High", "Low"))
meta$MP2_hl = factor(meta$MP2_hl, levels = c("High", "Low"))
meta$MP3_hl = factor(meta$MP3_hl, levels = c("High", "Low"))
meta$MP4_hl = factor(meta$MP4_hl, levels = c("High", "Low"))
meta$MP5_hl = factor(meta$MP5_hl, levels = c("High", "Low"))
meta$MP6_hl = factor(meta$MP6_hl, levels = c("High", "Low"))
meta$MP7_hl = factor(meta$MP7_hl, levels = c("High", "Low"))

cox_model <- coxph(Surv(`OS months`, Sstate) ~ MP1_hl+MP2_hl+MP3_hl+MP4_hl+MP5_hl+MP6_hl+MP7_hl, data = meta)
summary(cox_model)
p2 = ggforest(cox_model, data = meta, main = 'Hazard ratio - Post treatment')
p2

pplot = cowplot::plot_grid(p1, p2, align = 'hv', nrow = 2)
ggsave(filename = 'Figures/Fig6/6b_Forestplot_MPs.pdf', plot = pplot, device = 'pdf', height = 8, width = 7,
       dpi = 600, units = 'in', bg = 'white')

# 6C: MP4 TF Perturbation ----
load('scRNA_scATAC/SCENICplus/feather_output/MP_4/DF_ScoreTF_MP4_AfterPerturbation.RData')
# Label
label_df <- plot_df %>%
  group_by(TF) %>%
  filter(Iteration == max(Iteration)) %>%
  ungroup()

# iter_max <- max(plot_df$Iteration)
p = ggplot(plot_df, aes(x = Iteration, y = mean_log2FC, group = TF)) +
  geom_line(aes(color = TF), size = 0.5) +
  geom_point(size = 1.6, shape = 21, stroke = 0.3, color = "black", fill = "white") +
  geom_text_repel(
    data = label_df,
    aes(label = TF, color = TF),
    nudge_x = 0.4,             
    direction = "y",
    hjust = 0,
    size = 3,
    fontface = "bold",
    segment.color = "black",
    segment.linetype = "dashed",
    segment.size = 0.3,
    segment.ncp = 5,
    box.padding = 0.4,
    point.padding = 0.3,
    max.overlaps = Inf,
    show.legend = FALSE
  ) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "black", linewidth = 0.3) +
  scale_x_continuous(
    breaks = seq(0, max(plot_df$Iteration)),
    limits = c(0, max(plot_df$Iteration) + 1),
    expand = c(0.01, 0.01)
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  theme_classic(base_size = 11) +
  labs(
    x = "Perturbation iteration",
    y = expression("Average score")
  ) +
  theme(
    axis.text = element_text(size = 9),
    axis.text.x = element_text(size = 9, face = "bold"),
    axis.title.y = element_text(size = 10, face = "bold"),
    axis.title.x = element_text(size = 10),
    axis.line = element_line(color = "black"),
    panel.grid = element_blank(),
    legend.position = "none",
    plot.title = element_blank(),
    plot.margin = margin(5, 15, 5, 5)
  )

ggsave(filename = 'Figures/PerturbationIter_MP4_score.pdf', plot = p, device = 'pdf', width = 7.5, height = 6,
       dpi = 600, units = 'in', bg = 'white')

# 6D: NFATC2 activity ----
# ROC curve
load('Melanoma_Bulk/BulkICI_NFACT2_activity_PozReg.RData')
df = mmat[!is.na(mmat$Treatment), ]
df$Treatment = factor(df$Treatment, levels = c('PRE', 'POST'))
df$response_NR = factor(df$response_NR, levels = c('R', 'N'), labels = c('R', 'NR'))

fit_and_roc <- function(df, tf = "Net_", label = "PRE") {
  df <- df %>%
    mutate(Response = factor(response_NR, levels = c("R","N"),
                             labels = c("R","NR")),
           NFATC2_activity = as.numeric(NFATC2_activity),
           Net_ = as.numeric(Net_))
  
  form <- as.formula(paste("Response ~", tf))
  
  ctrl <- trainControl(
    method = "LOOCV",
    classProbs = TRUE,
    summaryFunction = twoClassSummary
  )
  
  model <- train(
    form,
    data = df,
    method = "glm",
    family = binomial(link = "logit"),
    trControl = ctrl,
    metric = "ROC"
  )
  
  prob_NR <- predict(model, type = "prob")[, "NR"]
  
  roc_obj <- pROC::roc(response = df$Response, predictor = prob_NR, quiet = TRUE)
  auc_val <- pROC::auc(roc_obj)
  list(roc_obj = roc_obj, auc = as.numeric(auc_val))
}

ddf_pre  <- mmat %>% filter(Treatment %in% "PRE")
res_pre  <- fit_and_roc(ddf_pre,  tf = "Net_", label = "PRE")
auc_PRE  <- round(res_pre$auc, 2)
res_pre = res_pre$roc_obj
ci.auc(res_pre)
# 95% CI: 0.5286-0.6611 (DeLong)

ddf_post <- mmat %>% filter(Treatment %in% "POST")
res_post <- fit_and_roc(ddf_post, tf = "Net_", label = "POST")
auc_POST <- round(res_post$auc, 2)
res_post = res_post$roc_obj
ci.auc(res_post)
# 95% CI: 0.5846-0.8146 (DeLong)

p = ggroc(list(PRE = res_pre, POST = res_post), legacy.axes = TRUE, size = 1.4) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey90") +
  coord_equal() +
  theme_minimal() +
  labs(title = "NFATC2 activity - PRE vs POST",
       x = "1 - Specificity", y = "Sensitivity",
       color = "") +
  scale_color_manual(values = c("PRE" = "#009688", "POST" = "#EF5350")) +
  scale_linetype_manual(breaks = c("PRE","POST"), values = c("PRE" = "dashed", "POST" = "solid")) +
  annotate("text", x = 0.2, y = 0.15,
           label = paste0("AUC PRE = ", auc_PRE, ' (CI: 0.53-0.66)'),
           color = "#009688", hjust = 0) +
  annotate("text", x = 0.2, y = 0.05,
           label = paste0("AUC POST = ", auc_POST, ' (CI: 0.58-0.81)'),
           color = "#EF5350", hjust = 0) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "bottom")

ggsave(filename = 'Figures/Fig6/6c_NFATC2_prediction.pdf', plot = p, device = 'pdf', height = 3.5, width = 5,
       dpi = 600, units = 'in', bg = 'white')

# Boxplot
p = ggplot(df, aes(x = response_NR, y = Net_, fill = response_NR)) +
  geom_boxplot(width = 0.3,  alpha = 0.6, position = position_dodge(width = 1), outliers = FALSE) +
  facet_wrap(~Treatment) +
  stat_compare_means(label.x = 1.5, label = 'p.signif') +
  theme_minimal(base_size = 11) +
  labs(y = 'Post transcriptional activity', x = '') +
  scale_fill_manual(values = c('NR' = '#5F9EA0', 'R' = '#FA8072')) +
  theme(
    strip.background = element_rect(fill = "#f0f0f0", color = NA),
    strip.text = element_text(size = 10, face = "bold", colour = 'black'),
    axis.text.x = element_text(size = 9, face = "bold"),
    axis.text.y = element_text(size = 9),
    axis.title.y = element_text(size = 10, face = "bold"),
    axis.line = element_line(colour = 'black'),
    panel.grid = element_blank(),
    legend.position = "none",
    panel.spacing.x = unit(0.2, "lines"),
    panel.spacing.y = unit(0.5, "lines")
  )

ggsave(filename = 'Figures/Fig6/6c_NFATC2_activity.pdf', plot = p, device = 'pdf', height = 3.5, width = 5,
       dpi = 600, units = 'in', bg = 'white')

# 6E: Metaprogram genes after NFATC2-KD ----
load('Melanoma_Bulk/GSE101323_nfatc2_ko/expr_data.RData')

selected_genes = list('Cell cycle' = c("CHEK1", "POLE", "BRCA2", "CENPF", "KIF23"),
                      'Melanocytic I' = c("TYRP1", "SALL4", "KCNJ13", "TSPAN10", "SERPINF1"),
                      'EMT/Hypoxia' = c("DDIT4", "VEGFA", "HK2", "PFKFB4", "CDKN2A"),
                      'Neural crest-like' = c("PRICKLE2", "SEMA3D", "COL1A2", "SDC2", "TFAP2B"),
                      'Antigen presentation/Interferon' = c("HLA-B", "TAP1", "B2M", "STAT1", "IFIH1"),
                      'Melanocytic II' = c("RUNX1", "SLCO3A1", "SYT14", "ABCB5", "POSTN"),
                      'WNT/B-catenin' = c("APCDD1", "ARHGAP29", "PRRX1", "SLC24A4", "ONECUT1")
)
unlist_selected_genes = Reduce(c, selected_genes)
length(intersect(rownames(expr_norm), unlist_selected_genes))

mat = expr_norm[unlist_selected_genes, pData$Sample[pData$sh %in% '1']]
mat_scaled = t(scale(t(mat)))
range(mat_scaled)
mat_scaled[mat_scaled < -2] = -2
mat_scaled[mat_scaled > 2] = 2

samples = pData$Sample[pData$sh %in% '1']
samples_col = pData$Type.col[pData$sh %in% '1']
names(samples_col) = pData$Sample[pData$sh %in% '1']

col_fun = colorRamp2(c(-2, 0, 2), c("blue", "grey90", "red")) 

split_rows = c(rep('Cell cycle', 5), 
               rep('Melanocytic I', 5),
               rep('EMT/Hypoxia', 5), 
               rep('Neural crest-like', 5),
               rep('Antigen presentation/Interferon', 5),
               rep('Melanocytic II', 5),
               rep('WNT/B-catenin', 5))

split_rows = factor(split_rows, levels = unique(split_rows))

mp = as.character(split_rows)
names(mp) = unlist_selected_genes
tmp_col = colori_mp
names(tmp_col) = names(selected_genes)

column_ha = HeatmapAnnotation(MP = mp,
                              col = list(MP = tmp_col), na_col = "white",
                              annotation_name_side = "left", show_legend = F, show_annotation_name = F)
row_ha = rowAnnotation(Condition = anno_simple(samples, col = samples_col), show_annotation_name = F)
mat_scaled = t(mat_scaled)

h1 = Heatmap(
  mat_scaled, name = 'Scaled expression',
  heatmap_legend_param = list(legend_direction = "horizontal"),
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE, 
  # row_split = split_rows, 
  column_split = split_rows,
  row_title = NULL, 
  left_annotation = row_ha,
  show_row_names = FALSE,
  show_column_names = TRUE, 
  column_title = NULL, column_names_rot = 45,
  # column_title_rot = 45, 
  gap = unit(1.5, 'mm'),
  top_annotation = column_ha, 
  row_names_side = "left"
)

h1

lgd = packLegend(Legend(labels = c("Ctrl", "shNFATC2"),
                        title = "Condition",
                        legend_gp = gpar(fill = c('lightblue', 'steelblue')), 
                        direction = 'horizontal'),
                 Legend(labels = names(tmp_col),
                        title = "MP",
                        legend_gp = gpar(fill = tmp_col),
                        direction = 'horizontal'),
                 direction = "horizontal",
                 gap = unit(6, "mm") 
)

pdf(file = 'Figures/Fig6/6e_Heatmap_MP_marker.pdf', width = 9, height = 4)
draw(h1, annotation_legend_list = lgd, heatmap_legend_side = "bottom", merge_legends = TRUE)
dev.off()

################################
#### Supplementary Figure 8 ####
################################
# S6A: ROC % metaprograms: Pozniak and NIBIT-M4 ----
load('pozniak_tmp_roc.RData', verbose = TRUE)
load('nibit_tmp_roc.RData', verbose = TRUE)
rm(df, sample_props)

# MP4 - PRE and POST
NIBIT_M4$PRE
Pozniak_MP_4$PRE
auc_nibit  <- round(as.numeric(auc(NIBIT_M4$PRE)), 2)
auc_pozniak <- round(as.numeric(auc(Pozniak_MP_4$PRE)), 2)
pROC::ci.auc(NIBIT_M4$PRE)
pROC::ci.auc(Pozniak_MP_4$PRE)

p1 <- ggroc(list(`NIBIT-M4` = NIBIT_M4$PRE, `Pozniak et al.` = Pozniak_MP_4$PRE), legacy.axes = TRUE, size = 1.4) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey90") +
  coord_equal() +
  theme_minimal() +
  labs(
    title = "MP4 (PRE ICI)",
    x = "1 - Specificity", 
    y = "Sensitivity",
    color = ""
    ) +
  scale_color_manual(
    values = c("NIBIT-M4" = "#6A3D9A", "Pozniak et al." = "#DAA520")
    ) +
  annotate(
    "text", x = 0.25, y = 0.15,
    label = paste0("AUC NIBIT-M4 = ", auc_nibit, ' (CI: 0.26-1)'),
    color = "#6A3D9A", hjust = 0
    ) +
  annotate(
    "text", x = 0.25, y = 0.05,
    label = paste0("AUC Pozniak et al. = ", auc_pozniak, ' (CI: 0.43-0.95)'),
    color = "#DAA520", hjust = 0
           ) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
    )
p1

NIBIT_M4$POST
Pozniak_MP_4$POST
auc_nibit  <- round(as.numeric(auc(NIBIT_M4$POST)), 2)
auc_pozniak <- round(as.numeric(auc(Pozniak_MP_4$POST)), 2)
pROC::ci.auc(NIBIT_M4$POST)
pROC::ci.auc(Pozniak_MP_4$POST)

p2 <- ggroc(list(`NIBIT-M4` = NIBIT_M4$POST, `Pozniak et al.` = Pozniak_MP_4$POST), legacy.axes = TRUE, size = 1.4) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey90") +
  coord_equal() +
  theme_minimal() +
  labs(
    title = "MP4 (PRE ICI)",
    x = "1 - Specificity", 
    y = "Sensitivity",
    color = ""
  ) +
  scale_color_manual(
    values = c("NIBIT-M4" = "#6A3D9A", "Pozniak et al." = "#DAA520")
    ) +
  annotate(
    "text", x = 0.25, y = 0.15,
    label = paste0("AUC NIBIT-M4 = ", auc_nibit, ' (CI: 1-1)'),
    color = "#6A3D9A", hjust = 0
    ) +
  annotate(
    "text", x = 0.25, y = 0.05,
    label = paste0("AUC Pozniak et al. = ", auc_pozniak, ' (CI: 0.56-1)'),
    color = "#DAA520", hjust = 0
    ) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
    )
p2

pplot = cowplot::plot_grid(p1, p2, align = 'hv')
ggsave(filename = 'Figures/Fig6/Supp_6a_ROC_percentage_MP4.pdf', plot = pplot, device = 'pdf', width = 7, height = 3.5,
       dpi = 600, units = 'in', bg = 'white')

# MP5 - PRE and POST
NIBIT_M5$PRE
Pozniak_MP_5$PRE
auc_nibit  <- round(as.numeric(auc(NIBIT_M5$PRE)), 2)
auc_pozniak <- round(as.numeric(auc(Pozniak_MP_5$PRE)), 2)
pROC::ci.auc(NIBIT_M5$PRE)
pROC::ci.auc(Pozniak_MP_5$PRE)

p1 <- ggroc(list(`NIBIT-M4` = NIBIT_M5$PRE, `Pozniak et al.` = Pozniak_MP_5$PRE), legacy.axes = TRUE, size = 1.4) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey90") +
  coord_equal() +
  theme_minimal() +
  labs(
    title = "MP5 (PRE ICI)",
    x = "1 - Specificity", y = "Sensitivity",
    color = ""
    ) +
  scale_color_manual(
    values = c("NIBIT-M4" = "#6A3D9A", "Pozniak et al." = "#DAA520")
    ) +
  annotate(
    "text", x = 0.25, y = 0.15,
    label = paste0("AUC NIBIT-M4 = ", auc_nibit, ' (CI: 1-1)'),
    color = "#6A3D9A", hjust = 0
    ) +
  annotate(
    "text", x = 0.25, y = 0.05,
    label = paste0("AUC Pozniak et al. = ", auc_pozniak, ' (CI: 0.48-0.97)'),
    color = "#DAA520", hjust = 0
    ) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
    )
p1

NIBIT_M5$POST
Pozniak_MP_5$POST
auc_nibit  <- round(as.numeric(auc(NIBIT_M5$POST)), 2)
auc_pozniak <- round(as.numeric(auc(Pozniak_MP_5$POST)), 2)
pROC::ci.auc(NIBIT_M5$POST)
pROC::ci.auc(Pozniak_MP_5$POST)

p2 <- ggroc(list(`NIBIT-M4` = NIBIT_M5$POST, `Pozniak et al.` = Pozniak_MP_5$POST), legacy.axes = TRUE, size = 1.4) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey90") +
  coord_equal() +
  theme_minimal() +
  labs(
    title = "MP5 (POST ICI)",
    x = "1 - Specificity", y = "Sensitivity",
    color = ""
    ) +
  scale_color_manual(
    values = c("NIBIT-M4" = "#6A3D9A", "Pozniak et al." = "#DAA520")
    ) +
  annotate(
    "text", x = 0.25, y = 0.15,
    label = paste0("AUC NIBIT-M4 = ", auc_nibit, ' (CI: 0.47-1)'),
    color = "#6A3D9A", hjust = 0
    ) +
  annotate(
    "text", x = 0.25, y = 0.05,
    label = paste0("AUC Pozniak et al. = ", auc_pozniak, ' (CI: 0.50-1)'),
    color = "#DAA520", hjust = 0
    ) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
    )
p2

pplot = cowplot::plot_grid(p1, p2, align = 'hv')
ggsave(filename = 'Figures/Fig6/Supp_6a_ROC_percentage_MP5.pdf', plot = pplot, device = 'pdf', width = 7, height = 3.5,
       dpi = 600, units = 'in', bg = 'white')

# S8B: Survival (MP4) ----
load('Melanoma_Bulk/dataframe_bulk_all.RData')
summary(mmat$purity)
mmat = mmat[mmat$purity >= 0.3, ]
# mmat$Treatment[mmat$Dataset %in% c('vanallen', 'amato')] = 'PRE'
table(mmat$Dataset, mmat$Therapy, exclude = NULL)
mmat$Treatment[mmat$Treatment %in% 'EDT'] = 'POST'
mmat$Treatment[mmat$Treatment %in% 'ON'] = 'POST'
mmat$Treatment[is.na(mmat$Treatment)] = 'Unknown'
mmat$Treatment[mmat$Treatment %in% ''] = 'Unknown'
table(mmat$Treatment, mmat$Dataset)
colnames(mmat)

# PRE
meta = mmat[mmat$Treatment %in% c('PRE'), ]
meta = meta[!is.na(meta$`vital status`), ]
meta = meta[!meta$`vital status` %in% '', ]
sum(is.na(meta$`overall survival (days)`))
meta = meta[!is.na(meta$`overall survival (days)`), ]
# View(table(meta$patient_name, meta$Dataset))

meta$`OS months` = meta$`overall survival (days)` / 30.44

table(meta$`vital status`, exclude = NULL)
meta$Sstate[meta$`vital status` %in% c('0:LIVING', 'Alive')] = 0
meta$Sstate[meta$`vital status` %in% c('1:DECEASED', 'Dead')] = 1
table(meta$`vital status`, meta$Sstate, exclude = NULL)

meta$`OS months`
meta$MP4_hl = ifelse(meta$MP_4 >= median(meta$MP_4), 'High', 'Low')
meta$MP4_hl = factor(meta$MP4_hl, levels = c("High", "Low"))

fit <- survfit(Surv(`OS months`, Sstate) ~ MP4_hl, data = meta)
fit
survdiff(Surv(`OS months`, Sstate) ~ MP4_hl, data = meta)

## change palette
p <- ggsurvplot(fit,
                pval = TRUE, conf.int = FALSE,
                pval.method = FALSE,
                pval.coord = c(1.2, 0.25), # p-value location
                pval.size = 4, # adjust p-value size
                title = "Before Treatment: MP4",
                risk.table = F, 
                risk.table.col = "strata",
                linetype = 1, 
                xlab = "OS (months)",
                # ggtheme = theme_bw(), # Change ggplot2 theme
                palette = c("#558B2F", "#AED581"),
                legend = "none",
                # legend = c(.9, .12),
                # legend.labs = levels(meta$MP4_hl),
                surv.scale = "percent",
)

p1 = ggpar(p, 
           font.main = c(12, "bold"),
           font.x = c(12, "bold"),
           font.y = c(12, "bold"),
           font.caption = c(12, "bold"), 
           font.legend = c(12),
           font.tickslab = c(8))

p1$plot + theme(legend.title = element_text(face = "bold", size = 12),
                legend.text = element_text(size = 10))

# Post treatment 
meta = mmat[mmat$Treatment %in% c('POST'), ]
meta = meta[!is.na(meta$`vital status`), ]
meta = meta[!meta$`vital status` %in% '', ]
sum(is.na(meta$`overall survival (days)`))
meta = meta[!is.na(meta$`overall survival (days)`), ]
# View(table(meta$patient_name, meta$Dataset))

meta$`OS months` = meta$`overall survival (days)` / 30.44

table(meta$`vital status`, exclude = NULL)
table(meta$M4_grade)

meta$Sstate[meta$`vital status` %in% c('0:LIVING', 'Alive')] = 0
meta$Sstate[meta$`vital status` %in% c('1:DECEASED', 'Dead')] = 1
table(meta$`vital status`, meta$Sstate, exclude = NULL)

meta$`OS months`
meta$MP4_hl = ifelse(meta$MP_4 >= median(meta$MP_4), 'High', 'Low')
meta$MP4_hl = factor(meta$MP4_hl, levels = c("High", "Low"))

fit <- survfit(Surv(`OS months`, Sstate) ~ MP4_hl, data = meta)
fit
survdiff(Surv(`OS months`, Sstate) ~ MP4_hl, data = meta)

p <- ggsurvplot(fit,
                pval = TRUE, conf.int = FALSE,
                pval.method = FALSE,
                pval.coord = c(1.2, 0.25), # p-value location
                pval.size = 4, # adjust p-value size
                title = "Post Treatment: MP4",
                risk.table = F, 
                risk.table.col = "strata",
                linetype = 1, 
                xlab = "OS (months)",
                # ggtheme = theme_bw(), # Change ggplot2 theme
                palette = c("#558B2F", "#AED581"),
                # legend = "none",
                legend = c(.8, .7),
                legend.labs = levels(meta$MP4_hl),
                surv.scale = "percent"
)

p2 = ggpar(p, 
           font.main = c(12, "bold"),
           font.x = c(12, "bold"),
           font.y = c(12, "bold"),
           font.caption = c(12, "bold"), 
           font.legend = c(12),
           font.tickslab = c(8))

p2$plot + theme(legend.title = element_text(face = "bold", size = 12),
                legend.text = element_text(size = 10))

comb = p1$plot / p2$plot
comb

ggsave(filename = 'Figures/Fig6/Supp6c_Survival_MP4.pdf', plot = comb, device = 'pdf', height = 7, width = 5,
       dpi = 600, units = 'in', bg = 'white')

# S8C: Multivariate Cox ----
# Cox multivariata con clinical (Stage e Gender)
# Load pData
load('Melanoma_Bulk/GSE78220/pData_Hugo.RData')
load('Melanoma_Bulk/GSE91061/pData_Riaz.RData')
load("Melanoma_Bulk/Nathanson_2017/pData_Nathanson.RData")
load("Melanoma_Bulk/PRJEB23709/pData_Gide.RData")
load('Melanoma_Bulk/phs000452/pData_Liu.RData')
load("Melanoma_Bulk/GSE145996/pData_Amato.RData")
load("Melanoma_Bulk/pData_Snyder.RData", verbose = TRUE)
updated = read_tsv(file = "Revisioni/Updates_pdata_bulk/skcm_mskcc_2014_clinical_data.tsv")
updated = updated[, c("Patient ID", "Age at Diagnosis", "Sex")]
colnames(pData_Snyder_melanoma)[2] = "Patient ID"
pData_Snyder_melanoma = merge(pData_Snyder_melanoma, updated, "Patient ID")
rm(updated)
load('Melanoma_Bulk/pData_VanAllen.RData', verbose = TRUE)
updated = readxl::read_xlsx("Revisioni/Updates_pdata_bulk/Van_Allen_pData.xlsx") %>% as.data.frame()
colnames(updated)[1] = "Patient ID"
colnames(pData_VanAllen_melanoma)[2] = "Patient ID"
pData_VanAllen_melanoma = merge(pData_VanAllen_melanoma, updated, "Patient ID")
rm(updated)
load("Melanoma_Bulk/GSE168204/pData_Du.RData", verbose = TRUE)
pData_Du_melanoma$Gender = NA

colnames(pData_Snyder_melanoma)[46:47] = c("age_start", "Gender")
colnames(pData_Snyder_melanoma)[3] = "sample_id"
colnames(pData_VanAllen_melanoma)[c(3, 47:48)] = c("sample_id", "age_start", "Gender")
colnames(pData_Du_melanoma)[c(1, 8)] = c("sample_id", "age_start")

pData_clinical = rbind(pData_Hugo_melanoma[, c("sample_id", "Gender", "age_start")],
                       pData_Riaz_melanoma[, c("sample_id", "Gender", "age_start")],
                       pData_Nathanson_melanoma[, c("sample_id", "Gender", "age_start")],
                       pData_Gide_melanoma[, c("sample_id", "Gender", "age_start")],
                       pData_Liu_melanoma[, c("sample_id", "Gender", "age_start")],
                       pData_Amato_melanoma[, c("sample_id", "Gender", "age_start")],
                       pData_Snyder_melanoma[, c("sample_id", "Gender", "age_start")],
                       pData_VanAllen_melanoma[, c("sample_id", "Gender", "age_start")],
                       pData_Du_melanoma[, c("sample_id", "Gender", "age_start")]
)
rm(pData_Amato_melanoma, pData_Du_melanoma, pData_Gide_melanoma, pData_Hugo_melanoma, pData_Liu_melanoma, 
   pData_Nathanson_melanoma, pData_Riaz_melanoma, pData_Snyder_melanoma, pData_VanAllen_melanoma)

load('Melanoma_Bulk/dataframe_bulk_all.RData')
summary(mmat$purity)
mmat = mmat[mmat$purity >= 0.3, ]
mmat$Treatment[mmat$Dataset %in% c('vanallen', 'amato')] = 'PRE'
table(mmat$Dataset, mmat$Therapy, exclude = NULL)
mmat$Treatment[mmat$Treatment %in% 'EDT'] = 'POST'
mmat$Treatment[mmat$Treatment %in% 'ON'] = 'POST'
mmat$Treatment[is.na(mmat$Treatment)] = 'Unknown'
mmat$Treatment[mmat$Treatment %in% ''] = 'Unknown'
table(mmat$Treatment, mmat$Dataset)
colnames(mmat)

pData_clinical = pData_clinical[pData_clinical$sample_id %in% intersect(pData_clinical$sample_id, rownames(mmat)), ]
pData_clinical = pData_clinical %>% as.data.frame()
rownames(pData_clinical) = pData_clinical$sample_id

mmat$sample_id = rownames(mmat)
mmat = merge(mmat, pData_clinical, by = "sample_id")
mmat[1:5, ]

table(mmat$Gender)
mmat$Gender[mmat$Gender %in% "female"] = "Female"
mmat$Gender[mmat$Gender %in% "male"] = "Male"
table(mmat$Gender, mmat$Treatment, exclude = NULL)

# PRE
meta = mmat[mmat$Treatment %in% c('PRE'), ]
meta = meta[!is.na(meta$age_start), ]
meta = meta[!is.na(meta$Gender), ]

meta = meta[!is.na(meta$`vital status`), ]
meta = meta[!meta$`vital status` %in% '', ]
sum(is.na(meta$`overall survival (days)`))
meta = meta[!is.na(meta$`overall survival (days)`), ]
nrow(meta) # 143 samples

meta$`OS months` = meta$`overall survival (days)` / 30.44

table(meta$`vital status`, exclude = NULL)
meta$Sstate[meta$`vital status` %in% c('0:LIVING', 'Alive')] = 0
meta$Sstate[meta$`vital status` %in% c('1:DECEASED', 'Dead')] = 1
table(meta$`vital status`, meta$Sstate, exclude = NULL)

meta$`OS months`
meta$`Cell Cycle` = ifelse(meta$MP_1 >= median(meta$MP_1), 'High', 'Low')
meta$`Melanocytic I` = ifelse(meta$MP_2 >= median(meta$MP_2), 'High', 'Low')
meta$`Hypoxia/EMT` = ifelse(meta$MP_3 >= median(meta$MP_3), 'High', 'Low')
meta$`Neural crest-like` = ifelse(meta$MP_4 >= median(meta$MP_4), 'High', 'Low')
meta$`Antigen presentation/\nInterferon` = ifelse(meta$MP_5 >= median(meta$MP_5), 'High', 'Low')
meta$`Melanocytic II` = ifelse(meta$MP_6 >= median(meta$MP_6), 'High', 'Low')
meta$`WNT/β-Catenin` = ifelse(meta$MP_7 >= median(meta$MP_7), 'High', 'Low')

meta$`Cell Cycle` = factor(meta$`Cell Cycle`, levels = c("High", "Low"))
meta$`Melanocytic I` = factor(meta$`Melanocytic I`, levels = c("High", "Low"))
meta$`Hypoxia/EMT` = factor(meta$`Hypoxia/EMT`, levels = c("High", "Low"))
meta$`Neural crest-like` = factor(meta$`Neural crest-like`, levels = c("High", "Low"))
meta$`Antigen presentation/\nInterferon` = factor(meta$`Antigen presentation/\nInterferon`, levels = c("High", "Low"))
meta$`Melanocytic II` = factor(meta$`Melanocytic II`, levels = c("High", "Low"))
meta$`WNT/β-Catenin` = factor(meta$`WNT/β-Catenin`, levels = c("High", "Low"))

colnames(meta)[26] = "Age"

cox <- coxph(
  Surv(`overall survival (days)`, Sstate) ~
    `Cell Cycle` + `Melanocytic I` + `Hypoxia/EMT` + `Neural crest-like` +
    `Antigen presentation/\nInterferon` + `Melanocytic II` + `WNT/β-Catenin` +
    Gender + 
    Age
  ,
  data = meta
)
summary(cox)
ggforest(cox, data = meta)
cox$coefficients
cox$coefficients[5]
names(cox$coefficients)[5]
names(cox$coefficients)[5] = "`Antigen presentation/\nInterferon`Low"
names(cox$coefficients)[5]
summary(cox)

p1 = ggforest(cox, data = meta, main = 'Hazard ratio - Pre treatment')

# POST
meta = mmat[mmat$Treatment %in% c('POST'), ]
meta = meta[!is.na(meta$age_start), ]
meta = meta[!is.na(meta$Gender), ]

meta = meta[!is.na(meta$`vital status`), ]
meta = meta[!meta$`vital status` %in% '', ]
sum(is.na(meta$`overall survival (days)`))
meta = meta[!is.na(meta$`overall survival (days)`), ]
nrow(meta) # 29 samples

meta$`OS months` = meta$`overall survival (days)` / 30.44

table(meta$`vital status`, exclude = NULL)
meta$Sstate[meta$`vital status` %in% c('0:LIVING', 'Alive')] = 0
meta$Sstate[meta$`vital status` %in% c('1:DECEASED', 'Dead')] = 1
table(meta$`vital status`, meta$Sstate, exclude = NULL)

meta$`OS months`
meta$`Cell Cycle` = ifelse(meta$MP_1 >= median(meta$MP_1), 'High', 'Low')
meta$`Melanocytic I` = ifelse(meta$MP_2 >= median(meta$MP_2), 'High', 'Low')
meta$`Hypoxia/EMT` = ifelse(meta$MP_3 >= median(meta$MP_3), 'High', 'Low')
meta$`Neural crest-like` = ifelse(meta$MP_4 >= median(meta$MP_4), 'High', 'Low')
meta$`Antigen presentation/\nInterferon` = ifelse(meta$MP_5 >= median(meta$MP_5), 'High', 'Low')
meta$`Melanocytic II` = ifelse(meta$MP_6 >= median(meta$MP_6), 'High', 'Low')
meta$`WNT/β-Catenin` = ifelse(meta$MP_7 >= median(meta$MP_7), 'High', 'Low')

meta$`Cell Cycle` = factor(meta$`Cell Cycle`, levels = c("High", "Low"))
meta$`Melanocytic I` = factor(meta$`Melanocytic I`, levels = c("High", "Low"))
meta$`Hypoxia/EMT` = factor(meta$`Hypoxia/EMT`, levels = c("High", "Low"))
meta$`Neural crest-like` = factor(meta$`Neural crest-like`, levels = c("High", "Low"))
meta$`Antigen presentation/\nInterferon` = factor(meta$`Antigen presentation/\nInterferon`, levels = c("High", "Low"))
meta$`Melanocytic II` = factor(meta$`Melanocytic II`, levels = c("High", "Low"))
meta$`WNT/β-Catenin` = factor(meta$`WNT/β-Catenin`, levels = c("High", "Low"))

colnames(meta)[26] = "Age"

cox <- coxph(
  Surv(`overall survival (days)`, Sstate) ~
    `Cell Cycle` + `Melanocytic I` + `Hypoxia/EMT` + `Neural crest-like` +
    `Antigen presentation/\nInterferon` + `Melanocytic II` + `WNT/β-Catenin` +
    Gender + 
    Age
  ,
  data = meta
)
summary(cox)
ggforest(cox, data = meta)

cox$coefficients
cox$coefficients[5]
names(cox$coefficients)[5]
names(cox$coefficients)[5] = "`Antigen presentation/\nInterferon`Low"
names(cox$coefficients)[5]
summary(cox)

p2 = ggforest(cox, data = meta, main = 'Hazard ratio - Post treatment')
p2

pplot = cowplot::plot_grid(p1, p2, align = 'hv', nrow = 2)
pplot

ggsave(filename = "Revisioni/forest_plot_Rev.pdf", 
       plot = pplot, device = 'pdf', height = 8, width = 7,
       dpi = 600, units = 'in', bg = 'white')

################################
#### Supplementary Figure 9 ####
################################
# S9B: 
load('Melanoma_Bulk/GSE101323_nfatc2_ko/Tsoi_MWW.RData')
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
  facet_grid(~ Var1, scales = "free_x", axes = "all_x") +
  labs(
    x = "",
    y = "NES",
    fill = ""
  )

ggsave(filename = 'Revisioni/tsoi_deconv_2.pdf', device = 'pdf', width = 12, height = 4,
       dpi = 600, units = 'in', bg = 'white')

# S9C: NES heatmap ----
load('Melanoma_Bulk/GSE101323_nfatc2_ko/ssMWWGst_top50up_MWW.RData')
load('Melanoma_Bulk/GSE101323_nfatc2_ko/expr_data.RData')

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
col_fun = colorRamp2(c(min(NES), 0,  max(NES)), c("blue", "white", "red"))

h1 = Heatmap(NES, 
             # width = unit(5, "in"), height = unit(5, "in"),
             heatmap_legend_param = list(legend_direction = "horizontal"),
             row_labels = c('Cell cylce', 'Melanocytic I', 'EMT/Hypoxia', 'Neural crest-like', 
                            'Antigen presentation/Interferon', 'Melanocytic II', 'Wnt/B-catenin'), 
             name = "NES", col = col_fun, # column_split = ssplit,
             top_annotation = column_ha, column_names_rot = 30,
             cluster_rows = F, show_column_names = FALSE, 
             row_title = "", cluster_columns = F)
h1

lgd = Legend(labels = c("Ctrl", "shNFATC2"),
             title = "Condition",
             legend_gp = gpar(fill = c('lightblue', 'steelblue'))
)

pdf(file = 'Figures/Fig6/Supp_6b_Heatmap_MP_NES.pdf', width = 5, height = 5)
draw(h1, heatmap_legend_side = "bottom")
dev.off()

