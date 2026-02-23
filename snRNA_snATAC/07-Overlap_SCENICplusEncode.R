###################
#### Libraries ####
###################
library(Seurat)
library(SeuratData)
library(Signac)
library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)
library(GenomicRanges)
library(S4Vectors)
library(BSgenome.Hsapiens.UCSC.hg38)
setwd('/home3/ciervo/scMULTIOME/Analisi/')
source('colori_finali.R')

# Loading cCRE encode ----
load('/home/noviello/TE_utilities/cCRE_ENCODE.RDa', verbose = TRUE)
CCRE_links = vroom::vroom("/home/noviello/TE_utilities/ENCODE/V4-hg38.Gene-Links.eQTLs.txt",col_names = F)
head(cCRE_ENCODE)
CCRE_links[1:5, ]
length(intersect(CCRE_links$X1, cCRE_ENCODE$X4)) # 0
length(intersect(CCRE_links$X1, cCRE_ENCODE$X5)) # 765816
colnames(mcols(cCRE_ENCODE))[2] = 'cCRE_link'
colnames(CCRE_links)[1] = 'cCRE_link'

CCRE_links$cCRE_link

# Loading prioritized TF per MP ----
TF = readxl::read_xlsx("Revisioni/SupplementaryTable7.xlsx", sheet = 2) %>% as.list()
TF = lapply(TF, na.omit)

# Loading eRegulon_metadata ----
eRegulon_metadata <- as.data.frame(read_csv("scRNA_scATAC/SCENICplus/CSV_output/eRegulon_metadata.csv"))
eRegulon_metadata = eRegulon_metadata[eRegulon_metadata$regulation == 1, ]

MP_list <- list(
  MP_1 = TF$`Cell cycle`,
  MP_2 = TF$`Melanocytic I`,
  MP_3 = TF$`Hypoxia/EMT`,
  MP_4 = TF$`Neural crest-like`,
  MP_5 = TF$`Antigen presentation/Interferon`,
  MP_6 = TF$`Melanocytic II`,
  MP_7 = TF$`WNT/B-catenin`
)

hg38_PromEnh_pELS <- cCRE_ENCODE[cCRE_ENCODE$X6 %in% "pELS"]
hg38_PromEnh_dELS <- cCRE_ENCODE[cCRE_ENCODE$X6 %in% "dELS"]

MP_results <- list()
for (mp_name in names(MP_list)) {
  mp_tfs <- unique(na.omit(MP_list[[mp_name]]))
  
  MP_df <- data.frame(
    row.names = mp_tfs,
    pELS = rep(NA_real_, length(mp_tfs)),
    dELS = rep(NA_real_, length(mp_tfs))
  )
  
  for (tf in mp_tfs) {
    message(paste0("[", mp_name, "] SCENIC+ regions per ", tf))
    regulon_regions <- eRegulon_metadata[eRegulon_metadata$TF %in% tf, "Region"]
    regulon_regions <- gsub(":", "-", regulon_regions)

    regions <- GRanges(
      seqnames = sapply(strsplit(regulon_regions, "-"), function(x) x[1]),
      ranges = IRanges(
        start = as.numeric(sapply(strsplit(regulon_regions, "-"), function(x) x[2])),
        end   = as.numeric(sapply(strsplit(regulon_regions, "-"), function(x) x[3]))
      )
    )
    message(paste0("[", mp_name, "] Overlap con pELS / dELS"))
    
    subset.hg38_pELS <- subsetByOverlaps(x = hg38_PromEnh_pELS, ranges = regions)
    subset.hg38_dELS <- subsetByOverlaps(x = hg38_PromEnh_dELS, ranges = regions)
    
    subset.hg38_pELS <- merge(subset.hg38_pELS, CCRE_links, by = "cCRE_link")
    subset.hg38_dELS <- merge(subset.hg38_dELS, CCRE_links, by = "cCRE_link")
    
    subset.hg38_pELS <- subset.hg38_pELS[!is.na(subset.hg38_pELS$X3), ]
    subset.hg38_dELS <- subset.hg38_dELS[!is.na(subset.hg38_dELS$X3), ]
    
    # pELS
    if (nrow(subset.hg38_pELS) > 0) {
      subset.hg38_pELS$region <- paste0(subset.hg38_pELS$seqnames, "_", subset.hg38_pELS$start, "_", subset.hg38_pELS$end)
      subset.hg38_pELS <- dplyr::distinct(subset.hg38_pELS, region, X3, .keep_all = TRUE)
      subset.hg38_pELS <- subset.hg38_pELS[subset.hg38_pELS$X4.y %in% "protein_coding", ]
      perc_pELS <- round(length(unique(subset.hg38_pELS$region)) / length(unique(regulon_regions)) * 100, 2)
    } else {
      perc_pELS <- 0
    }
    
    # dELS
    if (nrow(subset.hg38_dELS) > 0) {
      subset.hg38_dELS$region <- paste0(subset.hg38_dELS$seqnames, "_", subset.hg38_dELS$start, "_", subset.hg38_dELS$end)
      subset.hg38_dELS <- dplyr::distinct(subset.hg38_dELS, region, X3, .keep_all = TRUE)
      subset.hg38_dELS <- subset.hg38_dELS[subset.hg38_dELS$X4.y %in% "protein_coding", ]
      perc_dELS <- round(length(unique(subset.hg38_dELS$region)) / length(unique(regulon_regions)) * 100, 2)
    } else {
      perc_dELS <- 0
    }
    
    MP_df[tf, ] <- c(perc_pELS, perc_dELS)
    message(paste0("[", mp_name, "] ", tf, ": Done"))
  }
  
  MP_df$TF <- rownames(MP_df)
  MP_df$MP <- mp_name
  
  MP_results[[mp_name]] <- MP_df
}

MP_df_all <- do.call(rbind, MP_results)

DF_long <- MP_df_all %>%
  pivot_longer(cols = c("pELS", "dELS"), names_to = "variable", values_to = "value")

save(DF_long, file = "Revisioni/df_overlap_ENDODE_SCENIC.RData")
