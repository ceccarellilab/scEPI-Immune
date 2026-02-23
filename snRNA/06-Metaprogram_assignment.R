###################
#### Libraries ####
###################
library(Seurat)
library(dplyr)
library(ggplot2)
source('colori_finali.R')

simplicity_score = function(df, percentage = 0.1, ncores = 20) {
  cclass = parallel::mclapply(1:nrow(df), FUN = function(x) {
    tmp = unlist(df[x,])
    ifelse( (sort(tmp, decreasing = TRUE)[1] -  sort(tmp, decreasing = TRUE)[2])/sort(tmp, decreasing = TRUE)[1] >= percentage,
            names(sort(tmp, decreasing = TRUE))[1], 
            'non-classified')
  }, 
  mc.cores = ncores
  )
  cclass = unlist(cclass)
  print(sort(table(cclass), decreasing = TRUE))
  names(cclass) = rownames(df)
  return(cclass)
}

# Load seu object 
load('scRNA/RData/merged_object.RData')
load('scRNA/NMF/metaprograms.RData', verbose = TRUE)
rm(df)
seu = subset(seu, Cell_annotation %in% 'Malignant')
names(MP_list)

# AddModuleScore with the signatures ----
seu = AddModuleScore(seu, features = MP_list, ctrl = 50, seed = 1234)
View(seu@meta.data)
filter_1 = apply(seu@meta.data[, paste0('Cluster', 1:7)], 1, function(x) ifelse(max(x) > 0.3, names(which.max(x)), 'non-classified') )
tmp = seu@meta.data[names(filter_1)[!(filter_1 %in% 'non-classified')], paste0('Cluster', 1:7)]
filter_2 = simplicity_score(tmp, ncores = 30)
head(filter_2, 12)

seu = AddMetaData(seu, filter_2, col.name = 'Metaprogram_assignment')
seu$Metaprogram_assignment[is.na(seu$Metaprogram_assignment)] = 'non-classified'
table(seu$Metaprogram_assignment)

seu$Metaprogram_assignment[seu$Metaprogram_assignment %in% 'Cluster1'] = 'MP_1'
seu$Metaprogram_assignmnetn[seu$Metaprogram_assignment %in% 'Cluster2'] = 'MP_2'
seu$Metaprogram_assignmnetn[seu$Metaprogram_assignment %in% 'Cluster3'] = 'MP_3'
seu$Metaprogram_assignmnetn[seu$Metaprogram_assignment %in% 'Cluster4'] = 'MP_4'
seu$Metaprogram_assignment[seu$Metaprogram_assignment %in% 'Cluster5'] = 'MP_5'
seu$Metaprogram_assignment[seu$Metaprogram_assignment %in% 'Cluster6'] = 'MP_6'
seu$Metaprogram_assignment[seu$Metaprogram_assignment %in% 'Cluster7'] = 'MP_7'

save(seu, file = 'scRNA/RData/malignant_subset.RData')

assignment = seu@meta.data[, c(paste0('Cluster', 1:7), 'Metaprogram_assignment')]
save(assignment, file = 'scRNA/RData/assegnazione_metaprograms_SimplicityScore.RData')

##############################################
#### Correlation with Melanoma Signatures ####
##############################################
load("scRNA/RData/malignant_subset.RData")

# AddModuleScore signatures di Ponziak et al., 2024 DOI: 0.1016/j.cell.2023.11.037
load("/home/ciervo/EPICA/NIBIT/scRNA/Signatures_melanoma_Pozniak_100.RData", verbose = T)
seu = AddModuleScore(seu, features = signature, seed = 1234)
colnames(seu@meta.data)[30:ncol(seu@meta.data)] = names(signature)

# AddModuleScore signatures di Tsoi et al., 2018 DOI: 10.1016/j.ccell.2018.03.017
load("/home/caruso/Analisi2021/Maio_EPICA/dati_clinici/Anichini_signatures/Melanoma_Anichini_Signatures.RData", verbose = TRUE)
rm(AntiPD1_response, EMT_and_MELANOMA_diff_MARKERS, Pozniak_immune)

seu = AddModuleScore(seu, features = TSOI_differentiation_Stage, seed = 1234)
colnames(seu@meta.data)[39:45] = names(TSOI_differentiation_Stage)

# AddModuleScore signatures di Tirosh et al., 2016  DOI: 10.1126/science.aad0501
names(Tirosh_signatures)
Tirosh_signatures = Tirosh_signatures[1:4]

seu = AddModuleScore(seu, features = Tirosh_signatures, seed = 1234)
colnames(seu@meta.data)[46:49] = names(Tirosh_signatures)

scores_df = seu@meta.data[, c(21:27, 30:ncol(seu@meta.data))]
head(scores_df)
colnames(scores_df)

cor_result <- cor(scores_df, method = "pearson")

# dir.create('Supp_tables')
save(cor_result, file = 'Supp_tables/Signatures_MP_correlation.RData')

# AddModuleScore signatures di Baron et al., 2020 DOI: 10.1016/j.cels.2020.08.018
load("Revisioni/Signatures_Melanoma/Baron_signatures.RData", verbose = TRUE)
seu = AddModuleScore(seu, features = Baron_signatures, seed = 1234)
colnames(seu@meta.data)[50:52] = paste0("Baron_", names(Baron_signatures))
rm(Baron_signatures)

# AddModuleScore signatures di Hoek et al., 2008 DOI: 10.1158/0008-5472.CAN-07-2491
load("Revisioni/Signatures_Melanoma/Hoek_signatures.RData", verbose = TRUE)
seu = AddModuleScore(seu, features = Hoek_signatures, ctrl = max(unlist(lapply(Hoek_signatures, length))), seed = 1234)
colnames(seu@meta.data)
colnames(seu@meta.data)[53:54] = paste0("Hoek_", names(Hoek_signatures))
rm(Hoek_signatures)

# AddModuleScore signatures di Rambow et al., 2019 DOI: 10.1101/gad.329771.119
load("Revisioni/Signatures_Melanoma/Rambow_Signature.RData", verbose = TRUE)
seu = AddModuleScore(seu, features = Rambow_signatures, ctrl = max(unlist(lapply(Rambow_signatures, length))), seed = 1234)
colnames(seu@meta.data)
colnames(seu@meta.data)[55:58] = paste0("Rambow_", names(Rambow_signatures))
rm(Rambow_signatures, Rambow_sigs)

# AddModuleScore signatures di Verfaillie et al., 2015 DOI: 10.1038/ncomms7683
load("Revisioni/Signatures_Melanoma/Verfaillie_signatures.RData", verbose = TRUE)
seu = AddModuleScore(seu, features = Verfaillie_signatures, ctrl = max(unlist(lapply(Verfaillie_signatures, length))), seed = 1234)
colnames(seu@meta.data)
colnames(seu@meta.data)[59:60] = paste0("Verfaillie_", names(Verfaillie_signatures))
rm(Verfaillie_signatures)

# AddModuleScore signatures di Wouters et al., 2020 DOI: 10.1038/s41556-020-0547-3
load("Revisioni/Signatures_Melanoma/Wouters_signatures.RData", verbose = TRUE)
seu = AddModuleScore(seu, features = Wouters_signatures, seed = 1234)
colnames(seu@meta.data)
colnames(seu@meta.data)[61:63] = paste0("Wouters_", names(Wouters_signatures))
rm(Wouters_signatures)

scores_df = seu@meta.data[, c(21:27, 50:63)]
cor_result <- cor(scores_df, method = "pearson")
save(cor_result, file = 'Revisioni/Correlation_signatures_mel.RData')

# AddModuleScore signatures di Soldatov et al., 2019 DOI: 10.1126/science.aas9536
load("Revisioni/Signatures_Melanoma/Soldatov_signatures.RData", verbose = TRUE)
seu = AddModuleScore(seu, features = Soldatov_signatures, ctrl = 50, seed = 1234, search =TRUE)
colnames(seu@meta.data)
colnames(seu@meta.data)[64:69] = paste0("Soldatov_", names(Soldatov_signatures))
rm(Soldatov_signatures)
scores_df = seu@meta.data[, c(21:27, 64:69)]
cor_result <- cor(scores_df, method = "pearson")
save(cor_result, file = 'Revisioni/Correlation_signatures_Soldatov.RData')

##################
#### Stemness ####
##################
load("scRNA/RData/malignant_subset.RData")

stemness_signature = as.data.frame(TCGAbiolinks::SC_PCBC_stemSig)
stemness_signature = stemness_signature[stemness_signature$`TCGAbiolinks::SC_PCBC_stemSig` >= 0, , drop = FALSE]
stemness_signature_top100 = stemness_signature[order(stemness_signature$`TCGAbiolinks::SC_PCBC_stemSig`, decreasing = T), , drop = FALSE]
stemness_signature_top100 = stemness_signature_top100[1:100, ,drop = FALSE]

Signatures_Stemness = list("StemSig" = rownames(stemness_signature),
                           "StemSig_top100" = rownames(stemness_signature_top100),
                           "EctoStemSig" = rownames(ecto_stemness_signature),
                           "EctoStemSig_top100" = rownames(ecto_stemness_signature_top100))
# View(seu@meta.data)
seu = AddModuleScore(seu, features = Signatures_Stemness, name = "SigStem", search = TRUE, seed = 1234)
colnames(seu@meta.data)[30:33] = names(Signatures_Stemness)
save(seu, file = "scRNA/RData/malignant_subset.RData")
