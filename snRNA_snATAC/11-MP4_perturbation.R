library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(arrow)
library(scran)
library(readr)
library(parallel)

setwd('/home3/ciervo/scMULTIOME/Analisi/')
source('colori_finali.R')

# TF perturbation  ----
# Load metadata
load('scRNA_scATAC/RData/metadata.RData')
df_metadata = metadata[metadata$Metaprogram_assignment %in% paste0('MP_', 1:7), ]
rm(metadata)
df_metadata$Cell = rownames(df_metadata)

# Load metaprogram signature genes
metaprogram <- "MP_4"
load('scRNA/NMF/metaprograms.RData')
signature_genes <- MP_list[[metaprogram]]

# Directory with all TF folders
base_path <- 'scRNA_scATAC/SCENICplus/feather_output/MP_4'
TF_folders <- list.dirs(path = base_path, full.names = TRUE, recursive = FALSE)
tf_folder = TF_folders[1]
all_results <- mclapply(TF_folders, function(tf_folder) {
  tf_name <- basename(tf_folder)
  message("Processing TF: ", tf_name)
  
  # Read feather files per perturbation
  feather_files <- list.files(path = tf_folder, pattern = "*.feather", full.names = TRUE)
  llist <- lapply(feather_files, read_feather)
  llist <- lapply(llist, as.data.frame)
  names(llist) <- paste0("Perturbation_", seq_along(llist) - 1)
  
  llist <- lapply(llist, function(x) {
    x$index <- gsub("___cisTopic", "", x$index)
    rownames(x) <- x$index
    x$index <- NULL
    rownames(x) <- gsub("Pat03_W0_1", "Pat03_W0", rownames(x))
    return(x)
  })
  
  # Select cells belonging to MP_4
  selected_cells <- df_metadata %>%
    filter(Metaprogram_assignment == metaprogram) %>%
    pull(Cell)
  
  # Intersect genes of MP4
  signature_genes_present <- intersect(signature_genes, colnames(llist[[1]]))
  if (length(signature_genes_present) < 3) {
    warning("TF ", tf_name, " has too few genes in signature. Skipping.")
    return(NULL)
  }
  
  # Average expression
  averaged_list <- list()
  for (i in seq_along(llist)) {
    df_iter <- llist[[i]]
    cells_in_both <- intersect(rownames(df_iter), selected_cells)
    
    averaged <- df_iter[cells_in_both, signature_genes_present, drop = FALSE] %>%
      summarise(across(everything(), mean, na.rm = TRUE))
    
    averaged$Iteration <- i - 1
    averaged_list[[i]] <- averaged
  }
  
  # Baseline
  baseline_df <- averaged_list[[1]][, signature_genes_present, drop = FALSE]
  
  # Log2 fold change
  log2fc_df <- bind_rows(lapply(seq_along(averaged_list), function(i) {
    df <- averaged_list[[i]]
    log2fc <- df[signature_genes_present] - baseline_df
    log2fc$Iteration <- i - 1
    return(log2fc)
  }))
  
  # Mean log2FC per iteration
  log2fc_summary <- log2fc_df %>%
    rowwise() %>%
    mutate(mean_log2FC = mean(c_across(all_of(signature_genes_present)), na.rm = TRUE)) %>%
    ungroup() %>%
    dplyr::select(Iteration, mean_log2FC) %>%
    mutate(TF = tf_name)
  
  return(log2fc_summary)
}, mc.cores = 10)
names(all_results) <- basename(TF_folders)[seq_along(all_results)]

# Combine all TFs
plot_df <- bind_rows(all_results)
save(plot_df, file = 'scRNA_scATAC/SCENICplus/feather_output/MP_4/DF_ScoreTF_MP4_AfterPerturbation.RData')
