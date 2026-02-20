library(Seurat)
library(dplyr)
library(ggplot2)
library(caret)
library(pROC)
source('colori_finali.R')

setwd('/home3/ciervo/scMULTIOME/Analisi/')

# Load malignant subset by Pozniak et al., 2024
pozniak = readRDS('/home/mbesharat/Pozniak/Pozniak_cellannot/Malignant_MP_simp_score_annot_MP_auto_annot.Rds')
dim(pozniak)
View(pozniak@meta.data)

DefaultAssay(pozniak) = 'RNA'
sort(table(pozniak$sample_ID))
# scrCMA068  scrCMA044  scrCMA131  scrCMA121  scrCMA093  scrCMA130  scrCMA120  scrCMA119  scrCMA055  scrCMA064  scrCMA049  scrCMA041 sc5rCMA188  scrCMA094 sc5rCMA074 sc5rCMA070 
# 10         11         22         34         35         42         48         50         62         65        103        111        127        147        162        168 
# scrCMA046  scrCMA076  scrCMA072  scrCMA050  scrCMA048  scrCMA090  scrCMA089  scrCMA063 sc5rCMA136  scrCMA077  scrCMA054  scrCMA087 sc5rCMA144  scrCMA038  scrCMA040 sc5rCMA141 
# 186        189        203        220        367        384        398        402        514        541        573        613        659        811        939       1120 
# scrCMA091 sc5rCMA149  scrCMA036 
# 1273       1552       2067 

length(table(pozniak$sample_ID)) # 35 samples
length(which(table(pozniak$sample_ID) > quantile(as.numeric(table(pozniak$sample_ID)), .1))) # 31 samples

# Filter out samples with malignant cell counts <= the 10th percentile
pozniak = subset(pozniak, sample_ID %in% names(which(table(pozniak$sample_ID) > quantile(as.numeric(table(pozniak$sample_ID)), .1) )) )
dim(pozniak)

pozniak = pozniak %>%
  NormalizeData() %>%
  FindVariableFeatures(nfeatures = 3000) %>%
  ScaleData()

# Load metaprogram genes 
load('scRNA/NMF/metaprograms.RData')
rm(df)
pozniak = AddModuleScore(pozniak, features = MP_list, name = 'MP_', ctrl = 50, seed = 123)

# boxplot simplicity score ----
df = pozniak@meta.data
df$MP_annotation = apply(df[, paste0('MP_', 1:7)], 1, function(x) ifelse(max(x) > 0.05, names(which.max(x)), 'non-classified'))
pozniak = AddMetaData(pozniak, metadata = df$MP_annotation, 'Metaprogram_assignment')
save(pozniak, file = 'scRNA/RData/Pozniak_scRNA.RData')

#### 
load("scRNA/RData/Pozniak_scRNA.RData", verbose = TRUE)
df = pozniak@meta.data
sample_props <- df %>%
  count(patient_ID, Timepoint, Response, Metaprogram_assignment, name = "n") %>%
  tidyr::complete(
    tidyr::nesting(patient_ID, Timepoint, Response),
    Metaprogram_assignment,
    fill = list(n = 0)
  ) %>%
  group_by(patient_ID, Timepoint, Response) %>%
  mutate(
    total = sum(n),
    prop  = ifelse(total > 0, n / total, 0)
  ) %>%
  select(-total) %>%
  ungroup() %>%
  as.data.frame()


# LOOCV function
fit_roc <- function(df) {
  ctrl <- trainControl(method = "LOOCV", classProbs = TRUE,
                       summaryFunction = twoClassSummary)
  fit <- train(Response ~ prop,
               data = df,
               method = "glm",
               family = binomial(),
               trControl = ctrl,
               metric = "ROC")
  p_nr <- predict(fit, type = "prob")[,"NR"]                  # prob. classe positiva
  roc(df$Response, p_nr, levels = c("R","NR"), quiet = TRUE)
}

# MP4
mmat = sample_props[sample_props$Metaprogram_assignment %in% 'MP_4', ]
mmat$Response = ifelse(mmat$Response %in% "NonResponders", "NR", "R")
df_pre  <- filter(mmat, Timepoint == "BT")
df_post <- filter(mmat, Timepoint == "OT")

roc_PRE  <- fit_roc(df_pre)
roc_POST <- fit_roc(df_post)

auc_PRE  <- as.numeric(round(auc(roc_PRE), 2))
auc_POST <- as.numeric(round(auc(roc_POST), 2))
Pozniak_MP_4 = list(PRE = roc_PRE, POST = roc_POST)

# MP5
mmat = sample_props[sample_props$Metaprogram_assignment %in% 'MP_5', ]
mmat$Response = ifelse(mmat$Response %in% "NonResponders", "NR", "R")
df_pre  <- filter(mmat, Timepoint == "BT")
df_post <- filter(mmat, Timepoint == "OT")

roc_PRE  <- fit_roc(df_pre)
roc_POST <- fit_roc(df_post)

auc_PRE  <- as.numeric(round(auc(roc_PRE), 2))
auc_POST <- as.numeric(round(auc(roc_POST), 2))
Pozniak_MP_5 = list(PRE = roc_PRE, POST = roc_POST)

save(df, sample_props, Pozniak_MP_4, Pozniak_MP_5, file = 'pozniak_tmp_roc.RData')

# NIBIT-M4 cohort
load('scRNA_scATAC/RData/seu_multiome_filtered.RData')
df = seu.multiome@meta.data[seu.multiome$Metaprogram_assignment %in% paste0('MP_', 1:7), ]
table(df$Responder, df$Metaprogram_assignment)
rm(seu.multiome)

sample_props <- df %>%
  count(orig.ident, Week, Responder, Metaprogram_assignment, name = "n") %>%
  complete(
    nesting(orig.ident, Week, Responder),
    Metaprogram_assignment,
    fill = list(n = 0)
  ) %>%
  group_by(orig.ident, Week, Responder) %>%
  mutate(
    total = sum(n),
    prop  = ifelse(total > 0, n / total, 0)
  ) %>%
  select(-total) %>%
  ungroup() %>%
  as.data.frame()

table(sample_props$orig.ident, sample_props$Metaprogram_assignment)

# MP_4
mmat <- sample_props %>%
  filter(Metaprogram_assignment == "MP_4") %>%
  mutate(
    Treatment = ifelse(Week == "Week 0", "PRE", "POST"),
    prop      = as.numeric(as.character(prop)),
    Responder = factor(Responder, levels = c("R","NR"))  
  ) %>%
  tidyr::drop_na(prop)
mmat$prop = as.numeric(mmat$prop)

df_pre <- filter(mmat, Treatment == "PRE")
df_post <- filter(mmat, Treatment == "POST")

roc_PRE  <- fit_roc(df_pre)
roc_POST <- fit_roc(df_post)

auc_PRE  <- as.numeric(auc(roc_PRE))
auc_POST <- as.numeric(auc(roc_POST))

NIBIT_M4 = list(PRE = roc_PRE, POST = roc_POST)

# MP_5
mmat <- sample_props %>%
  filter(Metaprogram_assignment == "MP_5") %>%
  mutate(
    Treatment = ifelse(Week == "Week 0", "PRE", "POST"),
    prop      = as.numeric(as.character(prop)),
    Responder = factor(Responder, levels = c("R","NR"))  
  ) %>%
  tidyr::drop_na(prop)
mmat$prop = as.numeric(mmat$prop)

df_pre <- filter(mmat, Treatment == "PRE")
df_post <- filter(mmat, Treatment == "POST")

roc_PRE  <- fit_roc(df_pre)
roc_POST <- fit_roc(df_post)

auc_PRE  <- as.numeric(auc(roc_PRE))
auc_POST <- as.numeric(auc(roc_POST))

NIBIT_M5 = list(PRE = roc_PRE, POST = roc_POST)

save(df, sample_props, NIBIT_M4, NIBIT_M5, file = 'nibit_tmp_roc.RData')
