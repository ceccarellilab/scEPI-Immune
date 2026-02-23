####################################################
# Transposable Elements (TE) Multiome Analysis 
####################################################

# --- 0. Setup & Dependencies ---
library(vroom)
library(dplyr)
library(tibble)
library(Seurat)
library(scran)
library(ggplot2)
library(ggrepel)
library(ComplexHeatmap)
library(circlize)
library(reshape2)
library(tidyverse)
library(scales)
library(grid)
library(RColorBrewer)

set.seed(123)

# Define relative paths for reproducibility
base_dir    <- "."
data_dir    <- file.path(base_dir, "data")
results_dir <- file.path(base_dir, "results")
figures_dir <- file.path(results_dir, "figures")
objects_dir <- file.path(results_dir, "objects")

# Create directory structure
dirs_to_create <- c(results_dir, figures_dir, objects_dir, file.path(results_dir, "DEG"))
lapply(dirs_to_create, dir.create, showWarnings = FALSE, recursive = TRUE)


# --- 1. Utility Functions ---

# Function to read and merge MATES matrices
# TE locus-level matrices were produced using the MATES pipeline prior to downstream analysis.
# See: https://github.com/mcgilldinglab/MATES

read_mates_matrix <- function(sample, modality, mates_dir) {
  
  path_u <- file.path(mates_dir, paste0(sample, "_", modality), "Unique_TE_MTX.csv")
  path_m <- file.path(mates_dir, paste0(sample, "_", modality), "Multi_TE_MTX.csv")
  
  # Load unique and multi-mapping matrices
  mat_u <- vroom(path_u) %>% column_to_rownames(colnames(.)[1]) %>% as.matrix()
  mat_m <- vroom(path_m) %>% column_to_rownames(colnames(.)[1]) %>% as.matrix()
  
  # Handle potential differences in cell/gene sets
  genes <- union(rownames(mat_u), rownames(mat_m))
  cells <- union(colnames(mat_u), colnames(mat_m))
  
  m_u <- matrix(0, nrow = length(genes), ncol = length(cells), dimnames = list(genes, cells))
  m_m <- matrix(0, nrow = length(genes), ncol = length(cells), dimnames = list(genes, cells))
  
  m_u[rownames(mat_u), colnames(mat_u)] <- mat_u
  m_m[rownames(mat_m), colnames(mat_m)] <- mat_m
  
  # Combine and prefix cell names with sample ID
  out <- m_u + m_m
  colnames(out) <- paste0(sample, "_", colnames(out))
  
  return(out)
}

# Generic function to merge matrices by union of IDs
merge_union <- function(x, y) {
  full_join(
    rownames_to_column(as.data.frame(x), "ID"),
    rownames_to_column(as.data.frame(y), "ID"),
    by = "ID"
  ) %>%
    replace(is.na(.), 0) %>%
    column_to_rownames("ID")
}

# Seurat Object creation with Scran normalization
create_seurat_TE <- function(count_matrix, metadata) {
  
  # Filter common cells
  common_cells <- intersect(colnames(count_matrix), rownames(metadata))
  count_matrix <- count_matrix[, common_cells]
  
  sce <- SingleCellExperiment(assays = list(counts = count_matrix))
  sce <- sce[rowSums(counts(sce)) > 0, ]
  sce <- sce[, colSums(counts(sce)) > 3]
  
  # Normalization
  sce <- computeSumFactors(sce)
  sce <- logNormCounts(sce, log = FALSE)
  
  obj <- CreateSeuratObject(
    counts = assay(sce, "counts"),
    data   = log(normcounts(sce) + 1),
    meta.data = metadata[colnames(sce), ]
  )
  
  return(obj)
}


# --- 2. Data Loading & Processing ---

mates_dir <- file.path(data_dir, "raw/MATES/result_MTX")
samples   <- unique(gsub("_scRNA|_scATAC", "", list.files(mates_dir)))

# Process RNA and ATAC TE matrices
mat_RNA_list  <- lapply(samples, read_mates_matrix, modality = "scRNA", mates_dir = mates_dir)
mat_ATAC_list <- lapply(samples, read_mates_matrix, modality = "scATAC", mates_dir = mates_dir)

merged_TE_RNA  <- Reduce(merge_union, mat_RNA_list)
merged_TE_ATAC <- Reduce(merge_union, mat_ATAC_list)

# Load reference single-cell object for metadata
load(file.path(data_dir, "merged_object.RData")) 
df_metadata <- seu@meta.data

# Create Seurat objects for TEs
s_obj_RNA  <- create_seurat_TE(merged_TE_RNA, df_metadata)
s_obj_ATAC <- create_seurat_TE(merged_TE_ATAC, df_metadata)

save(s_obj_RNA, s_obj_ATAC, file = file.path(objects_dir, "s_obj_TE_combined.RData"))


# --- 3. Figure 5B: Treatment Induced TEs (Week 4 vs Week 0) ---

run_scoreMarkers <- function(obj, celltype_filter) {
  obj_sub <- subset(obj, Malignant %in% celltype_filter & Week %in% c("Week 0", "Week 4"))
  Idents(obj_sub) <- "Week"
  res <- scoreMarkers(obj_sub@assays$RNA$data, groups = Idents(obj_sub))
  return(as.data.frame(res$`Week 4`))
}

# Differential analysis
RNA_mal  <- run_scoreMarkers(s_obj_RNA,  "Malignant")
ATAC_mal <- run_scoreMarkers(s_obj_ATAC, "Malignant")
RNA_imm  <- run_scoreMarkers(s_obj_RNA,  "Non malignant")
ATAC_imm <- run_scoreMarkers(s_obj_ATAC, "Non malignant")

# Process for plotting
merged_df <- merge(ATAC_imm, RNA_imm, by = "row.names", suffixes = c("_ATAC","_RNA"))
colnames(merged_df)[1] <- "Gene"

# Load TE annotation 
# for human_TEs.csv see: https://github.com/mcgilldinglab/MATES
te_anno_path <- file.path(data_dir, "human_TEs.csv")
TE_anno <- vroom::vroom(te_anno_path) %>%
  distinct(TE_Name, .keep_all = TRUE) %>%
  mutate(TE_Name = gsub("_", "-", TE_Name)) %>%
  column_to_rownames("TE_Name")

# Mapping annotation
merged_df$TE_Fam <- TE_anno[merged_df$Gene, "TE_Fam"]
merged_df$TE_Fam_lab <- sapply(strsplit(merged_df$TE_Fam, "/"), `[`, 1)
merged_df$TE_Fam_lab <- ifelse(grepl("SINE|LINE|LTR", merged_df$TE_Fam_lab),
                               gsub("\\?", "", merged_df$TE_Fam_lab), "Others")
merged_df$TE_Fam_lab <- factor(merged_df$TE_Fam_lab, levels = c("SINE", "LINE", "LTR", "Others"))

# Calculate combined rank
merged_df$Mean_combined_Rank <- rowMeans(merged_df[, c("rank.logFC.cohen_ATAC", "rank.logFC.cohen_RNA")])
merged_df$Label <- ifelse(grepl("SINE|LINE|LTR", merged_df$TE_Fam), merged_df$Gene, "")

# Visualization: Figure 5B
custom_colors <- c("SINE" = "#E41A1C", "LINE" = "#4DAF4A", "LTR" = "#377EB8", "Others" = "gray70")

fig5b <- ggplot(merged_df, aes(x = rank.logFC.cohen_ATAC, y = rank.logFC.cohen_RNA,
                               size = Mean_combined_Rank, color = TE_Fam_lab)) +
  geom_point(alpha = 0.8) +
  scale_color_manual(values = custom_colors) +
  scale_y_reverse(limits = c(1000, 0)) +
  scale_x_reverse(limits = c(1000, 0)) +
  scale_size_continuous(range = c(3, 1)) +
  geom_text_repel(aes(label = Label), size = 3, max.overlaps = 50, fontface = "bold.italic") +
  theme_minimal() +
  labs(x = "Rank logFC (Accessibility)", y = "Rank logFC (Expression)", color = "TE Family") +
  guides(size = "none")

ggsave(file.path(figures_dir, "Figure5B_Treatment_induced_TEs.pdf"), fig5b, width = 6, height = 6)


# --- 4. Metaprogram (MP) Analysis - Figure 5D & S7 ---

run_MP_score <- function(obj) {
  obj_sub <- subset(obj, Malignant == "Malignant" & MP %in% paste0("MP_", 1:7))
  Idents(obj_sub) <- "MP"
  scoreMarkers(obj_sub@assays$RNA$data, groups = Idents(obj_sub))
}

MP_RNA  <- run_MP_score(s_obj_RNA)
MP_ATAC <- run_MP_score(s_obj_ATAC)

# Helper for heatmap matrix
prep_rank_mat <- function(datalist) {
  map_dfr(names(datalist), function(nm) {
    d <- as.data.frame(datalist[[nm]])
    tibble(Gene = rownames(d), MP = nm, Rank = d$rank.AUC)
  }) %>%
    pivot_wider(names_from = MP, values_from = Rank) %>%
    column_to_rownames("Gene") %>%
    as.matrix()
}

m_rna  <- prep_rank_mat(MP_RNA)
m_atac <- prep_rank_mat(MP_ATAC)
common  <- intersect(rownames(m_rna), rownames(m_atac))

# Visualization: Supplementary Figure S7
col_fun <- colorRamp2(c(1, 50, 100), c("#9E0142", "#F7F7F7", "#313695"))
ht_opt(legend_title_side = "topcenter")

ht_list <- Heatmap(m_atac[common,], name = "ATAC Rank", col = col_fun) +
  Heatmap(m_rna[common,], name = "RNA Rank", col = col_fun)

pdf(file.path(figures_dir, "TE_MP_Heatmap_S7.pdf"), width = 8, height = 10)
draw(ht_list, merge_legends = TRUE)
dev.off()


# Visualization: Figure 5D
prepare_MP_df <- function(datalist){
  
  df <- merge(mp_list[["MP_4"]],
              mp_list[["MP_2"]],
              by="row.names",
              suffixes=c("_MP4","_MP2"))
  
  out <- data.frame(
    Gene = df$Row.names,
    x = df$rank.logFC.cohen_MP2,
    y = df$rank.logFC.cohen_MP4
  )
  
  out$TE_Class <- TE_anno[out$Gene,"class"]
  return(out)
}

df_rna  <- prepare_MP_df(MP_RNA)
df_atac <- prepare_MP_df(MP_ATAC)

p_rna <- ggplot(df_rna, aes(x,y,color=TE_Class)) +
  geom_point(alpha=0.5) +
  theme_minimal() +
  labs(title="RNA MP2 vs MP4")

p_atac <- ggplot(df_atac, aes(x,y,color=TE_Class)) +
  geom_point(alpha=0.5) +
  theme_minimal() +
  labs(title="ATAC MP2 vs MP4")

ggsave(file.path(figures_dir,"MP2_vs_MP4_RNA_Figure5D.pdf"),
       p_rna, width=5, height=5)

ggsave(file.path(figures_dir,"MP2_vs_MP4_ATAC_Figure5D.pdf"),
       p_atac, width=5, height=5)



# --- 5. TE-TAA Analysis (ANTARES Integration) ---
# all_candidates.RData has been generated using ANTARES pipeline
# See: https://github.com/ceccarellilab/ANTARES
# HLA_patients.RData contains HLA info for each patient 

# Load ANTARES results
load(file.path(data_dir, "ANTARES/all_candidates.RData"))
load(file.path(data_dir, "ANTARES/HLA_patients.RData"))

all_df <- merge(all, HLA_patients, by = "HLA") %>%
  mutate(features_locus = paste0(
    gsub('\\:', '-', gsub('\\*', '', genomic_localization_RE)), '-', RE_name
  ))

# Extract expression data
get_te_taa_data <- function(obj, features) {
  mat  <- FetchData(obj, vars = features)
  meta <- obj@meta.data[rownames(mat), c('Metaprogram_assignment', 'Patient', 'Week')]
  cbind(meta, mat) %>%
    filter(Metaprogram_assignment %in% paste0("MP_", 1:7)) %>%
    rownames_to_column("Cell")
}

# [ s_obj_RNA_locus and s_obj_ATAC_locus came from MATES loci quantification]

exp_data  <- get_te_taa_data(s_obj_RNA_locus, unique(all_df$features_locus))
atac_data <- get_te_taa_data(s_obj_ATAC_locus, unique(all_df$features_locus))

meta_exp_df  <- melt(exp_data)
colnames(meta_exp_df)[5] <- "features_locus"

meta_atac_df <- melt(atac_data)
colnames(meta_atac_df)[5] <- "features_locus"

# Antigenic annotation
meta_exp_df$is_antigenic  <- 0
meta_atac_df$is_antigenic <- 0

tmp_antigenic <- unique(paste(all_df$Patient, all_df$features_locus))

meta_exp_df$is_antigenic[
  paste(meta_exp_df$Patient, meta_exp_df$features_locus) %in% tmp_antigenic
] <- 1

meta_atac_df$is_antigenic[
  paste(meta_atac_df$Patient, meta_atac_df$features_locus) %in% tmp_antigenic
] <- 1

# Keep common cells
common_cells <- intersect(unique(meta_atac_df$Cell), unique(meta_exp_df$Cell))

meta_atac_df <- meta_atac_df[meta_atac_df$Cell %in% common_cells, ]
colnames(meta_atac_df)[6] <- "ATAC"

meta_exp_df <- meta_exp_df[meta_exp_df$Cell %in% common_cells, ]
colnames(meta_exp_df)[6] <- "RNA"

meta_all <- meta_atac_df %>%
  dplyr::select(ATAC, Cell) %>%
  inner_join(meta_exp_df, by = "Cell")

# Compute proportions
summary_df <- meta_all %>%
  group_by(Metaprogram_assignment, Week) %>%
  summarise(
    total_expressed_antigenic = n_distinct(Cell[RNA != 0 & ATAC != 0 & is_antigenic == 1]),
    total_expressed           = n_distinct(Cell[RNA != 0 & ATAC != 0]),
    .groups = 'drop'
  )

total_df <- meta_exp_df %>%
  group_by(Metaprogram_assignment, Week) %>%
  summarise(total_cells = n_distinct(Cell), .groups = 'drop')

summary_df <- left_join(summary_df, total_df, by = c("Metaprogram_assignment", "Week"))

summary_df$Proportion <- (summary_df$total_expressed / summary_df$total_cells) * 100

# Plot Figure 5E - Immunoediting
colori_mp <- brewer.pal(7, "Dark2")

p_immuno <- ggplot(summary_df, aes(x = Week, y = Proportion, color = Metaprogram_assignment)) +
  geom_line(aes(group = Metaprogram_assignment), color = "gray80", size = 1) +
  geom_point(size = 4) +
  scale_color_manual(values = colori_mp) +
  labs(y = "Proportion (%) immunoedited cells", x = NULL, color = "Metaprogram") +
  theme_minimal() +
  theme(panel.grid.major.x = element_blank())

ggsave(file.path(figures_dir, "TE_TAA_immunoediting.pdf"), p_immuno, width = 6, height = 4)