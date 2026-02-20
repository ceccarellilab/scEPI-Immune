# SCENIC using Pozniak et al., 2024 dataset ----
# devtools::install_github("aertslab/SCENIC") 
# remotes::install_github("aertslab/SCopeLoomR")
# packageVersion("SCENIC")
# nBiocManager::install("RcisTarget")
library(Seurat)
library(SCENIC)
library(SingleCellExperiment)
library(ggplot2)
library(SCopeLoomR)
library(RcisTarget)
library(dplyr)
library(data.table)
library(tidyr)
library(stringr)
library(readr)
library(purrr)
source('/home/ciervo/Functions/Choose_nPCA.R')

setwd('/home3/ciervo/scMULTIOME/Analisi/')
# dir.create('scRNA/Pozniak_SCENIC')

# human:
# system("wget https://resources.aertslab.org/cistarget/databases/old/homo_sapiens/hg38/refseq_r80/mc9nr/gene_based/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather")
# system("wget https://resources.aertslab.org/cistarget/databases/old/homo_sapiens/hg38/refseq_r80/mc9nr/gene_based/hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather")
# system("wget https://resources.aertslab.org/cistarget/motif2tf/motifs-v9-nr.hgnc-m0.001-o0.0.tbl")
# system("wget https://resources.aertslab.org/cistarget/tf_lists/allTFs_hg38.txt")
# mc9nr: Motif collection version 9: 24k motifs

# Load data ----
load("scRNA/RData/Pozniak_scRNA.RData", verbose = TRUE)
dim(pozniak)
table(pozniak$Metaprogram_assignment, pozniak$Response)
pozniak = subset(pozniak, Metaprogram_assignment %in% "non-classified", invert = TRUE)
dim(pozniak)

annotation = pozniak@meta.data[ , c("Malignant_clusters", "Metaprogram_assignment")]
annotation$Barcode = rownames(annotation)
table(annotation$Metaprogram_assignment, annotation$Malignant_clusters)
write.table(annotation, file = 'scRNA/Pozniak_SCENIC/Metaprogram_annotation.txt', row.names = FALSE)
DimPlot(pozniak, group.by = 'Metaprogram_assignment')

# Prepare data to compute pyscenic
exprMat <- as.matrix(pozniak@assays$RNA$counts)
dim(exprMat) # 22897 genes x 12321 cells
cell.info <- pozniak@meta.data

c <- floor(ncol(exprMat) * 0.01) # 127
nonzero <- exprMat > 0L
keep_genes <- rowSums(as.matrix(nonzero)) >= c
exprMat <- exprMat[keep_genes, ]
dim(exprMat) # 15359 x 12321

# Further exclude mitochondrial genes
idx <- grep("^MT-", rownames(exprMat))
rownames(exprMat)[idx]
exprMat <- exprMat[-idx, ]
dim(exprMat) # 15346 x 12321

# Further exclude ribosomal genes
idx <- grep("^RP[LS]", rownames(exprMat))
rownames(exprMat)[idx]
exprMat <- exprMat[-idx, ]
dim(exprMat) # 15250 x 12321

# Extract default embedding (e.g. UMAP or PCA coordinates)
default.umap <- Embeddings(pozniak, reduction = "umap")
default.umap.name <- "UMAP"

SCopeLoomR::build_loom(file.name = "scRNA/Pozniak_SCENIC/pySCENIC_input_malignant.loom",
                       dgem = exprMat, 
                       default.embedding = default.umap,
                       default.embedding.name = "UMAP",
                       title = "Filtered RNA assay for pySCENIC",
                       genome = "Homo sapiens")
SCopeLoomR::add_annotated_clustering(loom = loom, group = cell.info$Metaprogram_assignment)
SCopeLoomR::close_loom(loom = loom)

write.table(exprMat, file = 'scRNA/Pozniak_SCENIC/pySCENIC_input_malignant.tsv', quote=FALSE, sep='\t', col.names = NA)

# Run pySCENIC using sbatch file "Script_pySCENIC.sbatch"

# Refine pySCENIC regulons by weight cutoff (per-TF 90th percentile) ----
# To obtain a high-confidence network, the initial regulons were refined by applying a weight-based cutoff:
# for each TF, the 90th percentile of its target connection weights was used as a threshold, and only
# interactions with weight >= this value were retained.

path_regs <- "scRNA/Pozniak_SCENIC/regulons.csv"               # regulon -> targets
path_adj  <- "scRNA/Pozniak_SCENIC/expr_mat.adjacencies.tsv"   # GRNBoost2: TF, target, weight/importance
out_csv   <- "scRNA/Pozniak_SCENIC/regulons_90pct.csv"         # output

min_targets     <- 10 
weight_quantile <- 0.90 

regs_raw <- fread(path_regs) %>% 
  as.data.frame()
colnames(regs_raw)[1] <- "regulon"
regs_raw[1:5, ]

extract_gene_names <- function(x) {
  if (is.na(x) || x == "") return(character(0))
  stringr::str_extract_all(x, "(?<=[\"'])[A-Za-z0-9._-]+(?=[\"'])")[[1]]
}

regs_edges <- regs_raw %>% 
  mutate(
    tf = str_replace(regulon, "_(motif|track).*", ""),
    targets_list = map(TargetGenes, extract_gene_names)
  ) %>% 
  tidyr::unnest_longer(targets_list, values_to = "target") %>% 
  transmute(tf = str_trim(tf), target = str_trim(target)) %>% 
  filter(!is.na(target), target != "") %>% 
  distinct(tf, target)
head(regs_edges)

# GRNBoost2
adj <- fread(path_adj) %>% 
  as.data.frame()
names(adj) = tolower(names(adj))
adj[1:5, ]
#         TF  target importance
# 1 HNRNPA1  IMPDH2   17.74031
# 2     UBB    SKP1   17.38152
# 3 HNRNPA1   BANF1   15.47248
# 4   UQCRB    TBCA   14.41454
# 5 HNRNPA1 SLC25A3   14.35478

adj <- adj %>% 
  mutate(across(c(TF,target), as.character)) %>% 
  distinct(TF, target, .keep_all = TRUE)
regs_edges[1:5, ]
#   tf     target 
# 1 ARID3A IL17B  
# 2 ARID3A ZNF449 
# 3 ARID3A SLC30A1
# 4 ARID3A DACT3  
# 5 ARID3A ATXN7 
colnames(adj)[1] = "tf"

edges <- regs_edges %>% 
  inner_join(adj, by = c("tf","target"))

tf_thr <- edges %>% 
  group_by(tf) %>% 
  summarise(w_thr = quantile(importance, probs = weight_quantile, na.rm = TRUE), .groups = "drop")

edges_f <- edges %>% 
  inner_join(tf_thr, by = "tf") %>% 
  filter(importance >= w_thr)

if(min_targets > 0){
  keep_tf <- edges_f %>% 
    count(tf, name = "n") %>% 
    filter(n >= min_targets) %>% pull(tf)
  edges_f <- edges_f %>% 
    filter(tf %in% keep_tf)
}

regulons_final <- edges_f %>% 
  arrange(tf, desc(importance)) %>% 
  group_by(tf) %>% 
  summarise(n_targets = n_distinct(target),
            targets = paste(unique(target), collapse = ","),
            .groups = "drop")

regulons_list <- setNames(
  lapply(strsplit(regulons_final$targets, ","), trimws),
  regulons_final$tf
)
regulons_list[1:2]
# $ATF3
# [1] "EGR1"       "STX1A"      "VIL1"       "SEMA4A"     "DUSP1"      "JUN"        "MMP23B"     "AVPR1A"     "SERTAD1"    "DGKH"       "EGR3"       "MCTP2"      "COL13A1"   
# [14] "DES"        "IER5"       "RGS2"       "SESN2"      "HSPA6"      "PALM2"      "PER1"       "BTBD19"     "NFKBIZ"     "MIR181A1HG" "LIPG"       "MCL1"       "ZCCHC14"   
# [27] "PTCH2"      "MAFF"       "EGR2"       "GLIS3"      "JUND"       "PMAIP1"     "IER3"       "CSRNP1"     "NR4A1"      "CXCL1"      "LHX2"       "SPEG"       "MKNK2"     
# [40] "PPP1R15A"  
# 
# $ATF4
# [1] "ST13"       "SARS"       "SOX10"      "PPARA"      "DDIT3"      "CEBPG"      "PCK2"       "NFE2L2"     "CIR1"       "GARS"       "CHAC1"      "ATF5"       "MYLIP"     
# [14] "ADM2"       "RBFOX2"     "TOX4"       "TXN2"       "SESN2"      "OSGIN1"     "HSD17B14"   "C6orf48"    "PATZ1"      "AARS"       "PBX2"       "DEDD2"      "YARS"      
# [27] "BEX2"       "C2orf49"    "SOX4"       "REXO4"      "STK19"      "TRIB3"      "PHF1"       "TNRC6B"     "NFIL3"      "ZC3H8"      "ZNF12"      "TRIB2"      "MIEF1"     
# [40] "PGPEP1"     "EP300"      "EIF4ENIF1"  "ZNF426"     "PLEKHO1"    "ZCCHC3"     "OCEL1"      "PPIG"       "SYF2"       "TSGA10"     "SYS1"       "LIPT1"      "GPT2"      
# [53] "CCNT2"      "TRIM16"     "ZFP90"      "LMBR1L"     "GGA1"       "KDM3A"      "DDX27"      "AMDHD2"     "PISD"       "UCHL1"      "TOB2"       "TCF25"      "TARS"      
# [66] "ZRSR2"      "NATD1"      "CDC42EP4"   "C8orf33"    "SUPT5H"     "EIF1"       "ZNF274"     "PEX26"      "GRK3"       "ARFGAP1"    "ZFAND2A"    "ZNF623"     "NARS"      
# [79] "COPS2"      "DNAJB2"     "AFF3"       "FXYD1"      "CMTM3"      "SRF"        "XPOT"       "KLC4"       "ERN1"       "ABCG1"      "KLHL24"     "DLGAP1-AS2" "ZNF397"    
# [92] "CWC22"      "EEF1D"

regulons_list["NFATC2"]
# $NFATC2
# [1] "EXPH5"     "MAP2"      "SEMA3D"    "FZD7"      "TFAP2C"    "DSTYK"     "MGLL"      "CAV1"      "SSBP2"     "ABCA5"     "RTN4RL1"   "LPAR1"     "EXTL1"     "PCSK6"    
# [15] "SCUBE2"    "CSPG4"     "CYFIP2"    "ST8SIA1"   "CHN2"      "EVI2A"     "TFAP2B"    "MATN2"     "TMEM168"   "SORBS1"    "APCDD1"    "SLC35F1"   "PRKCA"     "S100A4"  
# [29] "SPPL3"     "SPRY4"     "ABHD2"     "SLC1A4"    "SLCO4A1"   "CELF2"     "PSMG3-AS1" "ADAMTS5"   "ABCA10"    "MTUS1"     "WDR63"     "FAM149A"   "TIAM2"     "KANK1"    
# [43] "DMTN"      "ATP8A2"    "SGCD"      "COPZ2"     "PRELP"     "MAP3K3"    "APBA2"     "HOXB-AS1"  "SH3D19"    "FGFR1"     "KCNK2"     "SPARCL1"   "FGF1"      "VAMP5"   

save(regulons_list, file = 'scRNA/Pozniak_SCENIC/Regulon_SCENIC.RData')
