# Malanoma bulk - Response (downloaded from TIGER http://tiger.canceromics.org/#/immuneResponse) ----
setwd('/home3/ciervo/scMULTIOME/Analisi/')
library(estimate)
library(readr)
require(org.Hs.eg.db)

source('colori_finali.R')
source('/home/ciervo/Functions/quantileNormalization.R')
tpm_from_rpkm <- function(x) {
  rpkm.sum <- colSums(x)
  return(t(t(x) / (1e-06 * rpkm.sum)))
}

# GSE78220 (Hugo et al., 2016 Cell) ----
Hugo_melanoma = as.data.frame(readRDS('Melanoma_Bulk/GSE78220/Melanoma-GSE78220.Response.Rds'))
pData_Hugo_melanoma = readRDS('Melanoma_Bulk/GSE78220/Melanoma-GSE78220.Response_Clinical.Rds')

Hugo_melanoma[1:5, 1:5]
rownames(Hugo_melanoma) = Hugo_melanoma$GENE_SYMBOL
Hugo_melanoma$GENE_SYMBOL = NULL

Hugo_melanoma = as.matrix(Hugo_melanoma)
range(colSums(Hugo_melanoma))
boxplot(log2(Hugo_melanoma+1)[, 1:5])

Hugo_melanoma = t(quantileNormalization(t(Hugo_melanoma)))
boxplot(log2(Hugo_melanoma+1)[, 1:5])



df_bulk <- data.frame(GeneSymbol = rownames(Hugo_melanoma), Hugo_melanoma)
df_bulk <- rbind(c("GeneSymbol", colnames(Hugo_melanoma)), df_bulk)
df_bulk[1:10, 1:5]

dir.create('scRNA/RData/Melanoma_bulk')
write.table(df_bulk, quote = FALSE, row.names = F, col.names = F, sep = "\t", 
            file = "scRNA/RData/Melanoma_bulk/RNAbulk_Hugo_CIBERSORT.txt")

Hugo_melanoma = log2(Hugo_melanoma + 1)
range(Hugo_melanoma)

# dir.create('scRNA/RData/Melanoma_bulk/Matrices')
save(Hugo_melanoma, file = 'scRNA/RData/Melanoma_bulk/Matrices/Hugo.RData')

# dir.create("scRNA/RData/Melanoma_bulk/ESTIMATE")
Hugo_melanoma <- Hugo_melanoma[, pData_Hugo_melanoma$sample_id]

write.table(Hugo_melanoma, sep = "\t", file = "scRNA/RData/Melanoma_bulk/ESTIMATE/inputEstimate_Hugo.txt", quote = F)
filterCommonGenes(input.f = "scRNA/RData/Melanoma_bulk/ESTIMATE/inputEstimate_Hugo.txt", output.f = "scRNA/RData/Melanoma_bulk/ESTIMATE/estimate_Hugo.gct", id = "GeneSymbol")
estimateScore("scRNA/RData/Melanoma_bulk/ESTIMATE/estimate_Hugo.gct", output.ds = "scRNA/RData/Melanoma_bulk/ESTIMATE/estimateScore_Hugo.gct", platform  = "illumina")

## Loading estimate results
tableScores <-read.table(file = "scRNA/RData/Melanoma_bulk/ESTIMATE/estimateScore_Hugo.gct", quote="", header = TRUE, sep = "\t", stringsAsFactors = FALSE, skip = 2)
rownames(tableScores) <- tableScores$NAME
tableScores$NAME <- NULL
tableScores$Description <- NULL
tableScores <- t(as.matrix(tableScores))
identical(rownames(tableScores), pData_Hugo_melanoma$sample_id)

pData_Hugo_melanoma = cbind(pData_Hugo_melanoma, tableScores)
pData_Hugo_melanoma$purity = cos(0.6049872018 + 0.0001467884 * pData_Hugo_melanoma$ESTIMATEScore)
save(pData_Hugo_melanoma, file = 'Melanoma_Bulk/GSE78220/pData_Hugo.RData')

# Plot 
load('Melanoma_Bulk/GSE78220/pData_Hugo.RData')
clinical = pData_Hugo_melanoma
rm(pData_Hugo_melanoma)
CIBERSORTx_Hugo_Results <- read_delim("Melanoma_Bulk/CIBERSORTx_Hugo_Results.txt", 
                                      delim = "\t", escape_double = FALSE, 
                                      trim_ws = TRUE)

mat = as.data.frame(CIBERSORTx_Hugo_Results[, 1:8])
rownames(mat) = mat$Mixture
rownames(clinical) = clinical$sample_id

mat = mat[rownames(clinical), ]
identical(rownames(clinical), rownames(mat))

colnames(clinical)
mat = cbind(mat, clinical[, c("sample_id", "response_NR", 'purity')])

for(barcode in mat$sample_id){
  message(barcode)
  ppurity = mat[mat$sample_id %in% barcode, 'purity']
  tmp2 <- as.matrix(mat[mat$sample_id %in% barcode, 2:8])
  tmp2 = tmp2*ppurity
  mat[mat$sample_id %in% barcode, 2:8] = tmp2
}
rowSums(mat[, 2:8])
# mat$TME = 1-mat$purity

# llevels = mat$Mixture
mat$Mixture = NULL
mat$purity = NULL
mat$sample_id = NULL

mmat = reshape2::melt(mat, c("response_NR"))
table(mmat$response_NR)
mmat$response = factor(mmat$response_NR, levels = c("R", "N"), labels = c('R', 'NR'))

p1 = ggplot(data = mmat, aes(x = response, y = value, fill = response)) + 
  # geom_violin(alpha=0.7) +
  geom_boxplot(width=0.3,  alpha = 0.6, position = position_dodge(width = 1)) +
  theme_classic(base_size = 14) + 
  theme(
    axis.text.y = element_text(size = 10),         # Show y-axis text
    axis.text.x = element_text(size = 10),         # Show x-axis text
    axis.title.x = element_blank(),                # Remove x-axis title
    # axis.title.y = element_blank(),                # Remove y-axis title
    strip.background = element_blank(),            # Remove strip background
    # strip.text.y = element_blank(),                # Remove strip text
    legend.position = "none",                       # Adjust legend position
    # panel.border = element_rect(colour = "black", size = 0.5), # Add a border to the plot
    panel.spacing = unit(0.05, "lines"),
    plot.title = element_text(hjust = 0.5, size = 14), # Center title
    plot.subtitle = element_text(hjust = 0.5, size = 10) # Center subtitle
  ) + 
  facet_grid(~variable, scales = 'free_x') +
  stat_compare_means(method = "wilcox.test", label = "p.format", label.x.npc = 'center', label.y = max(as.numeric(mmat$value), na.rm = TRUE) * 1.05) +
  scale_fill_manual(
    values = c('NR' = '#5F9EA0', 'R' = '#FA8072'),
    labels = c('NR', 'R')) +
  labs(
    title = "Hugo et al. 2016",
    subtitle = "anti-PD1",
    x = "",
    y = "Cell fraction",
    fill = "Response"
  )
p1

# GSE91061 (Riaz et al., 2017 Cell) ----
Riaz_melanoma = as.data.frame(readRDS('Melanoma_Bulk/GSE91061/Melanoma-GSE91061.Response.Rds'))
pData_Riaz_melanoma = readRDS('Melanoma_Bulk/GSE91061/Melanoma-GSE91061.Response_Clinical.Rds')

Riaz_melanoma[1:5, 1:5]
rownames(Riaz_melanoma) = Riaz_melanoma$GENE_SYMBOL
Riaz_melanoma$GENE_SYMBOL = NULL

Riaz_melanoma = as.matrix(Riaz_melanoma)
range(Riaz_melanoma)
boxplot(log2(Riaz_melanoma + 1)[, 1:5])
Riaz_melanoma = t(quantileNormalization(t(Riaz_melanoma)))
boxplot(log2(Riaz_melanoma + 1)[, 1:5])

df_bulk <- data.frame(GeneSymbol = rownames(Riaz_melanoma), Riaz_melanoma)
df_bulk <- rbind(c("GeneSymbol", colnames(Riaz_melanoma)), df_bulk)
df_bulk[1:10, 1:5]

# dir.create('scRNA/RData/Melanoma_bulk')
write.table(df_bulk, quote = FALSE, row.names = F, col.names = F, sep = "\t", 
            file = "scRNA/RData/Melanoma_bulk/RNAbulk_Riaz_CIBERSORT.txt")

Riaz_melanoma = log2(Riaz_melanoma + 1)
range(Riaz_melanoma)

save(Riaz_melanoma, file = 'scRNA/RData/Melanoma_bulk/Matrices/Riaz.RData')

dim(Riaz_melanoma)
# dir.create("scRNA/RData/Melanoma_bulk/ESTIMATE")
Riaz_melanoma <- Riaz_melanoma[, pData_Riaz_melanoma$sample_id]

write.table(Riaz_melanoma, sep = "\t", file = "scRNA/RData/Melanoma_bulk/ESTIMATE/inputEstimate_Riaz.txt", quote = F)
filterCommonGenes(input.f = "scRNA/RData/Melanoma_bulk/ESTIMATE/inputEstimate_Riaz.txt", output.f = "scRNA/RData/Melanoma_bulk/ESTIMATE/estimate_Riaz.gct", id = "GeneSymbol")
estimateScore("scRNA/RData/Melanoma_bulk/ESTIMATE/estimate_Riaz.gct", output.ds = "scRNA/RData/Melanoma_bulk/ESTIMATE/estimateScore_Riaz.gct", platform  = "illumina")

## Loading estimate results
tableScores <-read.table(file = "scRNA/RData/Melanoma_bulk/ESTIMATE/estimateScore_Riaz.gct", quote="", header = TRUE, sep = "\t", stringsAsFactors = FALSE, skip = 2)
rownames(tableScores) <- tableScores$NAME
tableScores$NAME <- NULL
tableScores$Description <- NULL
tableScores <- t(as.matrix(tableScores))
identical(rownames(tableScores), pData_Riaz_melanoma$sample_id)

pData_Riaz_melanoma = cbind(pData_Riaz_melanoma, tableScores)
pData_Riaz_melanoma$purity = cos(0.6049872018 + 0.0001467884 * pData_Riaz_melanoma$ESTIMATEScore)
range(pData_Riaz_melanoma$purity)
save(pData_Riaz_melanoma, file = 'Melanoma_Bulk/GSE91061/pData_Riaz.RData')

# Plot ----
library(readr)
load('Melanoma_Bulk/GSE91061/pData_Riaz.RData')
clinical = pData_Riaz_melanoma
rm(pData_Riaz_melanoma)
CIBERSORTx_Riaz_Results <- read_delim("Melanoma_Bulk/CIBERSORTx_Riaz_Results.txt", 
                                      delim = "\t", escape_double = FALSE, 
                                      trim_ws = TRUE)

mat = as.data.frame(CIBERSORTx_Riaz_Results[, 1:8])
rownames(mat) = mat$Mixture
rownames(clinical) = clinical$sample_id

mat = mat[rownames(clinical), ]
identical(rownames(clinical), rownames(mat))

colnames(clinical)
mat = cbind(mat, clinical[, c("sample_id", "response_NR", 'purity')])

for(barcode in mat$sample_id){
  message(barcode)
  ppurity = mat[mat$sample_id %in% barcode, 'purity']
  tmp2 <- as.matrix(mat[mat$sample_id %in% barcode, 2:8])
  tmp2 = tmp2*ppurity
  mat[mat$sample_id %in% barcode, 2:8] = tmp2
}
rowSums(mat[, 2:8])
# mat$TME = 1-mat$purity

# llevels = mat$Mixture
mat$Mixture = NULL
mat$purity = NULL
mat$sample_id = NULL

mmat = reshape2::melt(mat, c("response_NR"))
table(mmat$response_NR)
mmat = mmat[!mmat$response_NR %in% 'UNK', ]
mmat$response = factor(mmat$response_NR, levels = c("R", "N"), labels = c('R', 'NR'))

p2 = ggplot(data = mmat, aes(x = response, y = value, fill = response)) + 
  # geom_violin(alpha=0.7) +
  geom_boxplot(width=0.3,  alpha = 0.6, position = position_dodge(width = 1)) +
  theme_classic(base_size = 14) + 
  theme(
    axis.text.y = element_text(size = 10),         # Show y-axis text
    axis.text.x = element_text(size = 10),         # Show x-axis text
    axis.title.x = element_blank(),                # Remove x-axis title
    # axis.title.y = element_blank(),                # Remove y-axis title
    strip.background = element_blank(),            # Remove strip background
    # strip.text.y = element_blank(),                # Remove strip text
    legend.position = "none",                       # Adjust legend position
    # panel.border = element_rect(colour = "black", size = 0.5), # Add a border to the plot
    panel.spacing = unit(0.05, "lines"),
    plot.title = element_text(hjust = 0.5, size = 14), # Center title
    plot.subtitle = element_text(hjust = 0.5, size = 10) # Center subtitle
  ) + 
  facet_grid(~variable, scales = 'free_x') +
  stat_compare_means(method = "wilcox.test", label = "p.format", label.x.npc = 'center', label.y = max(as.numeric(mmat$value), na.rm = TRUE) * 1.05) +
  scale_fill_manual(
    values = c('NR' = '#5F9EA0', 'R' = '#FA8072'),
    labels = c('NR', 'R')) +
  labs(
    title = "Riaz et al. 2017",
    subtitle = "anti-PD1",
    x = "",
    y = "Cell fraction",
    fill = "Response"
  )
p2

# Nathanson (Nathanson et al., 2017 Cancer Immunol Res) ----
Nathanson_melanoma = as.data.frame(readRDS('Melanoma_Bulk/Nathanson_2017/Melanoma-Nathanson_2017.Response.Rds'))
pData_Nathanson_melanoma = readRDS('Melanoma_Bulk/Nathanson_2017/Melanoma-Nathanson_2017.Response_clinical.Rds')

Nathanson_melanoma[1:5, 1:5]
rownames(Nathanson_melanoma) = Nathanson_melanoma$GENE_SYMBOL
Nathanson_melanoma$GENE_SYMBOL = NULL

Nathanson_melanoma = as.matrix(Nathanson_melanoma)
Nathanson_melanoma = na.omit(Nathanson_melanoma)
range(Nathanson_melanoma)

Nathanson_melanoma = tpm_from_rpkm(Nathanson_melanoma)
range(Nathanson_melanoma)

Nathanson_melanoma = as.matrix(Nathanson_melanoma)
range(Nathanson_melanoma)
boxplot(log2(Nathanson_melanoma+1)[, 1:5])
Nathanson_melanoma = t(quantileNormalization(t(Nathanson_melanoma)))
boxplot(log2(Nathanson_melanoma+1)[, 1:5])


df_bulk <- data.frame(GeneSymbol = rownames(Nathanson_melanoma), Nathanson_melanoma)
df_bulk <- rbind(c("GeneSymbol", colnames(Nathanson_melanoma)), df_bulk)
df_bulk[1:10, 1:5]

# dir.create('scRNA/RData/Melanoma_bulk')
write.table(df_bulk, quote = FALSE, row.names = F, col.names = F, sep = "\t", 
            file = "scRNA/RData/Melanoma_bulk/RNAbulk_Nathanson_CIBERSORT.txt")

Nathanson_melanoma = log2(Nathanson_melanoma + 1)
range(Nathanson_melanoma)

save(Nathanson_melanoma, file = 'scRNA/RData/Melanoma_bulk/Matrices/Nathanson.RData')

dim(Nathanson_melanoma)
# dir.create("scRNA/RData/Melanoma_bulk/ESTIMATE")
pData_Nathanson_melanoma = pData_Nathanson_melanoma[pData_Nathanson_melanoma$sample_id %in% intersect(colnames(Nathanson_melanoma), pData_Nathanson_melanoma$sample_id), ]
table(pData_Nathanson_melanoma$response_NR)

Nathanson_melanoma <- Nathanson_melanoma[, pData_Nathanson_melanoma$sample_id]

write.table(Nathanson_melanoma, sep = "\t", file = "scRNA/RData/Melanoma_bulk/ESTIMATE/inputEstimate_Nathanson.txt", quote = F)
filterCommonGenes(input.f = "scRNA/RData/Melanoma_bulk/ESTIMATE/inputEstimate_Nathanson.txt", output.f = "scRNA/RData/Melanoma_bulk/ESTIMATE/estimate_Nathanson.gct", id = "GeneSymbol")
estimateScore("scRNA/RData/Melanoma_bulk/ESTIMATE/estimate_Nathanson.gct", output.ds = "scRNA/RData/Melanoma_bulk/ESTIMATE/estimateScore_Nathanson.gct", platform  = "illumina")

## Loading estimate results
tableScores <-read.table(file = "scRNA/RData/Melanoma_bulk/ESTIMATE/estimateScore_Nathanson.gct", quote="", header = TRUE, sep = "\t", stringsAsFactors = FALSE, skip = 2)
rownames(tableScores) <- tableScores$NAME
tableScores$NAME <- NULL
tableScores$Description <- NULL
tableScores <- t(as.matrix(tableScores))
identical(rownames(tableScores), paste0('X', pData_Nathanson_melanoma$sample_id))

pData_Nathanson_melanoma = cbind(pData_Nathanson_melanoma, tableScores)
pData_Nathanson_melanoma$purity = cos(0.6049872018 + 0.0001467884 * pData_Nathanson_melanoma$ESTIMATEScore)
range(pData_Nathanson_melanoma$purity)
save(pData_Nathanson_melanoma, file = 'Melanoma_Bulk/Nathanson_2017/pData_Nathanson.RData')

# Plot ----
library(readr)
load('Melanoma_Bulk/Nathanson_2017/pData_Nathanson.RData')
clinical = pData_Nathanson_melanoma
rm(pData_Nathanson_melanoma)
CIBERSORTx_Nathanson_Results <- read_delim("Melanoma_Bulk/CIBERSORTx_Nathanson_Results.txt", 
                                           delim = "\t", escape_double = FALSE, 
                                           trim_ws = TRUE)

mat = as.data.frame(CIBERSORTx_Nathanson_Results[, 1:8])
rownames(mat) = mat$Mixture
rownames(clinical) = clinical$sample_id

mat = mat[rownames(clinical), ]
identical(rownames(clinical), rownames(mat))

colnames(clinical)
mat = cbind(mat, clinical[, c("sample_id", "response_NR", 'purity')])

for(barcode in mat$sample_id){
  message(barcode)
  ppurity = mat[mat$sample_id %in% barcode, 'purity']
  tmp2 <- as.matrix(mat[mat$sample_id %in% barcode, 2:8])
  tmp2 = tmp2*ppurity
  mat[mat$sample_id %in% barcode, 2:8] = tmp2
}
rowSums(mat[, 2:8])
# mat$TME = 1-mat$purity

# llevels = mat$Mixture
mat$Mixture = NULL
mat$purity = NULL
mat$sample_id = NULL

mmat = reshape2::melt(mat, c("response_NR"))
table(mmat$response_NR)
# mmat = mmat[!mmat$response_NR %in% 'UNK', ]
mmat$response = factor(mmat$response_NR, levels = c("R", "N"), labels = c('R', 'NR'))

p3 = ggplot(data = mmat, aes(x = response, y = value, fill = response)) + 
  # geom_violin(alpha=0.7) +
  geom_boxplot(width=0.3,  alpha = 0.6, position = position_dodge(width = 1)) +
  theme_classic(base_size = 14) + 
  theme(
    axis.text.y = element_text(size = 10),         # Show y-axis text
    axis.text.x = element_text(size = 10),         # Show x-axis text
    axis.title.x = element_blank(),                # Remove x-axis title
    # axis.title.y = element_blank(),                # Remove y-axis title
    strip.background = element_blank(),            # Remove strip background
    # strip.text.y = element_blank(),                # Remove strip text
    legend.position = "none",                       # Adjust legend position
    # panel.border = element_rect(colour = "black", size = 0.5), # Add a border to the plot
    panel.spacing = unit(0.05, "lines"),
    plot.title = element_text(hjust = 0.5, size = 14), # Center title
    plot.subtitle = element_text(hjust = 0.5, size = 10) # Center subtitle
  ) + 
  facet_grid(~variable, scales = 'free_x') +
  stat_compare_means(method = "wilcox.test", label = "p.format", label.x.npc = 'center', label.y = max(as.numeric(mmat$value), na.rm = TRUE) * 1.05) +
  scale_fill_manual(
    values = c('NR' = '#5F9EA0', 'R' = '#FA8072'),
    labels = c('NR', 'R')) +
  labs(
    title = "Nathanson et al. 2017",
    subtitle = "anti-CTLA4",
    x = "",
    y = "Cell fraction",
    fill = "Response"
  )
p3
table(clinical$Therapy)

# Gide (Gide et al., 2019 Cancer Cell) ----
Gide_melanoma = as.data.frame(readRDS('Melanoma_Bulk/PRJEB23709/Melanoma-PRJEB23709.Response.Rds'))
pData_Gide_melanoma = readRDS('Melanoma_Bulk/PRJEB23709/Melanoma-PRJEB23709.Response_clinical.Rds')

Gide_melanoma[1:5, 1:5]
rownames(Gide_melanoma) = Gide_melanoma$GENE_SYMBOL
Gide_melanoma$GENE_SYMBOL = NULL

Gide_melanoma = as.matrix(Gide_melanoma)
Gide_melanoma = na.omit(Gide_melanoma)
range(Gide_melanoma)

table(pData_Gide_melanoma$response_NR)
# N  R 
# 42 49

Gide_melanoma = as.matrix(Gide_melanoma)
range(Gide_melanoma)
boxplot(log2(Gide_melanoma+1)[, 1:5])
Gide_melanoma = t(quantileNormalization(t(Gide_melanoma)))
boxplot(log2(Gide_melanoma+1)[, 1:5])


df_bulk <- data.frame(GeneSymbol = rownames(Gide_melanoma), Gide_melanoma)
df_bulk <- rbind(c("GeneSymbol", colnames(Gide_melanoma)), df_bulk)
df_bulk[1:10, 1:5]

# dir.create('scRNA/RData/Melanoma_bulk')
write.table(df_bulk, quote = FALSE, row.names = F, col.names = F, sep = "\t", 
            file = "scRNA/RData/Melanoma_bulk/RNAbulk_Gide_CIBERSORT.txt")

Gide_melanoma = log2(Gide_melanoma + 1)
range(Gide_melanoma)

save(Gide_melanoma, file = 'scRNA/RData/Melanoma_bulk/Matrices/Gide.RData')

# dir.create("scRNA/RData/Melanoma_bulk/ESTIMATE")
dim(Gide_melanoma)
Gide_melanoma <- Gide_melanoma[, pData_Gide_melanoma$sample_id]

write.table(Gide_melanoma, sep = "\t", file = "scRNA/RData/Melanoma_bulk/ESTIMATE/inputEstimate_Gide.txt", quote = F)
filterCommonGenes(input.f = "scRNA/RData/Melanoma_bulk/ESTIMATE/inputEstimate_Gide.txt", output.f = "scRNA/RData/Melanoma_bulk/ESTIMATE/estimate_Gide.gct", id = "GeneSymbol")
estimateScore("scRNA/RData/Melanoma_bulk/ESTIMATE/estimate_Gide.gct", output.ds = "scRNA/RData/Melanoma_bulk/ESTIMATE/estimateScore_Gide.gct", platform  = "illumina")

## Loading estimate results
tableScores <-read.table(file = "scRNA/RData/Melanoma_bulk/ESTIMATE/estimateScore_Gide.gct", quote="", header = TRUE, sep = "\t", stringsAsFactors = FALSE, skip = 2)
rownames(tableScores) <- tableScores$NAME
tableScores$NAME <- NULL
tableScores$Description <- NULL
tableScores <- t(as.matrix(tableScores))
identical(rownames(tableScores), pData_Gide_melanoma$sample_id)

pData_Gide_melanoma = cbind(pData_Gide_melanoma, tableScores)
pData_Gide_melanoma$purity = cos(0.6049872018 + 0.0001467884 * pData_Gide_melanoma$ESTIMATEScore)
save(pData_Gide_melanoma, file = 'Melanoma_Bulk/PRJEB23709/pData_Gide.RData')

# Plot ----
load('Melanoma_Bulk/PRJEB23709/pData_Gide.RData')
clinical = pData_Gide_melanoma
rm(pData_Gide_melanoma)
CIBERSORTx_Gide_Results <- read_delim("Melanoma_Bulk/CIBERSORTx_Gide_Results.txt", 
                                      delim = "\t", escape_double = FALSE, 
                                      trim_ws = TRUE)

mat = as.data.frame(CIBERSORTx_Gide_Results[, 1:8])
rownames(mat) = mat$Mixture
rownames(clinical) = clinical$sample_id

mat = mat[rownames(clinical), ]
identical(rownames(clinical), rownames(mat))

colnames(clinical)
mat = cbind(mat, clinical[, c("sample_id", "response_NR", 'purity')])

for(barcode in mat$sample_id){
  message(barcode)
  ppurity = mat[mat$sample_id %in% barcode, 'purity']
  tmp2 <- as.matrix(mat[mat$sample_id %in% barcode, 2:8])
  tmp2 = tmp2*ppurity
  mat[mat$sample_id %in% barcode, 2:8] = tmp2
}
rowSums(mat[, 2:8])
# mat$TME = 1-mat$purity

# llevels = mat$Mixture
mat$Mixture = NULL
mat$purity = NULL
mat$sample_id = NULL

mmat = reshape2::melt(mat, c("response_NR"))
table(mmat$response_NR)
# mmat = mmat[!mmat$response_NR %in% 'UNK', ]
mmat$response = factor(mmat$response_NR, levels = c("R", "N"), labels = c('R', 'NR'))

p4 = ggplot(data = mmat, aes(x = response, y = value, fill = response)) + 
  # geom_violin(alpha=0.7) +
  geom_boxplot(width=0.3,  alpha = 0.6, position = position_dodge(width = 1)) +
  theme_classic(base_size = 14) + 
  theme(
    axis.text.y = element_text(size = 10),         # Show y-axis text
    axis.text.x = element_text(size = 10),         # Show x-axis text
    axis.title.x = element_blank(),                # Remove x-axis title
    # axis.title.y = element_blank(),                # Remove y-axis title
    strip.background = element_blank(),            # Remove strip background
    # strip.text.y = element_blank(),                # Remove strip text
    legend.position = "none",                       # Adjust legend position
    # panel.border = element_rect(colour = "black", size = 0.5), # Add a border to the plot
    panel.spacing = unit(0.05, "lines"),
    plot.title = element_text(hjust = 0.5, size = 14), # Center title
    plot.subtitle = element_text(hjust = 0.5, size = 10) # Center subtitle
  ) + 
  facet_grid(~variable, scales = 'free_x') +
  stat_compare_means(method = "wilcox.test", label = "p.format", label.x.npc = 'center', label.y = max(as.numeric(mmat$value), na.rm = TRUE) * 1.05) +
  scale_fill_manual(
    values = c('NR' = '#5F9EA0', 'R' = '#FA8072'),
    labels = c('NR', 'R')) +
  labs(
    title = "Gide et al. 2019",
    subtitle = "anti-CTLA4+anti-PD1 / anti-PD1",
    x = "",
    y = "Cell fraction",
    fill = "Response"
  )
p4
table(clinical$Therapy)

# phs000452 (Liu et al. 2019, Nat Med) ----
Liu_melanoma = as.data.frame(readRDS('Melanoma_Bulk/phs000452/Melanoma-phs000452.Response.Rds'))
pData_Liu_melanoma = readRDS('Melanoma_Bulk/phs000452/Melanoma-phs000452.Response_Clinical.Rds')

Liu_melanoma[1:5, 1:5]
rownames(Liu_melanoma) = Liu_melanoma$GENE_SYMBOL
Liu_melanoma$GENE_SYMBOL = NULL

Liu_melanoma = as.matrix(Liu_melanoma)
range(Liu_melanoma)
boxplot(log2(Liu_melanoma+1)[, 1:5])
Liu_melanoma = t(quantileNormalization(t(Liu_melanoma)))
boxplot(log2(Liu_melanoma+1)[, 1:5])


df_bulk <- data.frame(GeneSymbol = rownames(Liu_melanoma), Liu_melanoma)
df_bulk <- rbind(c("GeneSymbol", colnames(Liu_melanoma)), df_bulk)
df_bulk[1:10, 1:5]

# dir.create('scRNA/RData/Melanoma_bulk')
write.table(df_bulk, quote = FALSE, row.names = F, col.names = F, sep = "\t", 
            file = "scRNA/RData/Melanoma_bulk/RNAbulk_Liu_CIBERSORT.txt")

Liu_melanoma = log2(Liu_melanoma + 1)
range(Liu_melanoma)

save(Liu_melanoma, file = 'scRNA/RData/Melanoma_bulk/Matrices/Liu.RData')

# dir.create("scRNA/RData/Melanoma_bulk/ESTIMATE")
Liu_melanoma <- Liu_melanoma[, pData_Liu_melanoma$sample_id]

write.table(Liu_melanoma, sep = "\t", file = "scRNA/RData/Melanoma_bulk/ESTIMATE/inputEstimate_Liu.txt", quote = F)
filterCommonGenes(input.f = "scRNA/RData/Melanoma_bulk/ESTIMATE/inputEstimate_Liu.txt", output.f = "scRNA/RData/Melanoma_bulk/ESTIMATE/estimate_Liu.gct", id = "GeneSymbol")
estimateScore("scRNA/RData/Melanoma_bulk/ESTIMATE/estimate_Liu.gct", output.ds = "scRNA/RData/Melanoma_bulk/ESTIMATE/estimateScore_Liu.gct", platform  = "illumina")

## Loading estimate results
tableScores <-read.table(file = "scRNA/RData/Melanoma_bulk/ESTIMATE/estimateScore_Liu.gct", quote="", header = TRUE, sep = "\t", stringsAsFactors = FALSE, skip = 2)
rownames(tableScores) <- tableScores$NAME
tableScores$NAME <- NULL
tableScores$Description <- NULL
tableScores <- t(as.matrix(tableScores))
identical(rownames(tableScores), pData_Liu_melanoma$sample_id)

pData_Liu_melanoma = cbind(pData_Liu_melanoma, tableScores)
pData_Liu_melanoma$purity = cos(0.6049872018 + 0.0001467884 * pData_Liu_melanoma$ESTIMATEScore)

pData_Liu_melanoma = pData_Liu_melanoma[grep('^P', pData_Liu_melanoma$patient_name), ]
rownames(pData_Liu_melanoma) = pData_Liu_melanoma$sample_id

save(pData_Liu_melanoma, file = 'Melanoma_Bulk/phs000452/pData_Liu.RData')

# Plot ----
library(readr)
load('Melanoma_Bulk/phs000452/pData_Liu.RData')
clinical = pData_Liu_melanoma
rm(pData_Liu_melanoma)
CIBERSORTx_Liu_Results <- read_delim("Melanoma_Bulk/CIBERSORTx_Liu_Results.txt", 
                                     delim = "\t", escape_double = FALSE, 
                                     trim_ws = TRUE)

mat = as.data.frame(CIBERSORTx_Liu_Results[, 1:8])
rownames(mat) = mat$Mixture
rownames(clinical) = clinical$sample_id

mat = mat[rownames(clinical), ]
identical(rownames(clinical), rownames(mat))

colnames(clinical)
mat = cbind(mat, clinical[, c("sample_id", "response_NR", 'purity')])

for(barcode in mat$sample_id){
  message(barcode)
  ppurity = mat[mat$sample_id %in% barcode, 'purity']
  tmp2 <- as.matrix(mat[mat$sample_id %in% barcode, 2:8])
  tmp2 = tmp2*ppurity
  mat[mat$sample_id %in% barcode, 2:8] = tmp2
}
rowSums(mat[, 2:8])
# mat$TME = 1-mat$purity

# llevels = mat$Mixture
mat$Mixture = NULL
mat$purity = NULL
mat$sample_id = NULL

mmat = reshape2::melt(mat, c("response_NR"))
table(mmat$response_NR)
mmat$response = factor(mmat$response_NR, levels = c("R", "N"), labels = c('R', 'NR'))

p5 = ggplot(data = mmat, aes(x = response, y = value, fill = response)) + 
  # geom_violin(alpha=0.7) +
  geom_boxplot(width=0.3,  alpha = 0.6, position = position_dodge(width = 1)) +
  theme_classic(base_size = 14) + 
  theme(
    axis.text.y = element_text(size = 10),         # Show y-axis text
    axis.text.x = element_text(size = 10),         # Show x-axis text
    axis.title.x = element_blank(),                # Remove x-axis title
    # axis.title.y = element_blank(),                # Remove y-axis title
    strip.background = element_blank(),            # Remove strip background
    # strip.text.y = element_blank(),                # Remove strip text
    legend.position = "none",                       # Adjust legend position
    # panel.border = element_rect(colour = "black", size = 0.5), # Add a border to the plot
    panel.spacing = unit(0.05, "lines"),
    plot.title = element_text(hjust = 0.5, size = 14), # Center title
    plot.subtitle = element_text(hjust = 0.5, size = 10) # Center subtitle
  ) + 
  facet_grid(~variable, scales = 'free_x') +
  stat_compare_means(method = "wilcox.test", label = "p.format", label.x.npc = 'center', label.y = max(as.numeric(mmat$value), na.rm = TRUE) * 1.05) +
  scale_fill_manual(
    values = c('NR' = '#5F9EA0', 'R' = '#FA8072'),
    labels = c('NR', 'R')) +
  labs(
    title = "Liu et al. 2019",
    subtitle = "anti-PD1",
    x = "",
    y = "Cell fraction",
    fill = "Response"
  )
p5
table(clinical$Therapy)

# GSE145996 (Amato et al. 2020, Cancers (Basel)) ----
Amato_melanoma = as.data.frame(readRDS('Melanoma_Bulk/GSE145996/Melanoma_GSE145996.Response.Rds'))
pData_Amato_melanoma = readRDS('Melanoma_Bulk/GSE145996/Melanoma_GSE145996.Response_Clinical.Rds')

Amato_melanoma[1:5, 1:5]
rownames(Amato_melanoma) = Amato_melanoma$GENE_SYMBOL
Amato_melanoma$GENE_SYMBOL = NULL

Amato_melanoma = as.matrix(Amato_melanoma)
Amato_melanoma = na.omit(Amato_melanoma)
range(Amato_melanoma)
boxplot(log2(Amato_melanoma+1))
Amato_melanoma = t(quantileNormalization(t(Amato_melanoma)))
boxplot(log2(Amato_melanoma+1))


df_bulk <- data.frame(GeneSymbol = rownames(Amato_melanoma), Amato_melanoma)
df_bulk <- rbind(c("GeneSymbol", colnames(Amato_melanoma)), df_bulk)
df_bulk[1:10, 1:5]

# dir.create('scRNA/RData/Melanoma_bulk')
write.table(df_bulk, quote = FALSE, row.names = F, col.names = F, sep = "\t", 
            file = "scRNA/RData/Melanoma_bulk/RNAbulk_Amato_CIBERSORT.txt")

Amato_melanoma = log2(Amato_melanoma + 1)
range(Amato_melanoma)

save(Amato_melanoma, file = 'scRNA/RData/Melanoma_bulk/Matrices/Amato.RData')

# dir.create("scRNA/RData/Melanoma_bulk/ESTIMATE")
Amato_melanoma <- Amato_melanoma[, pData_Amato_melanoma$sample_id]

write.table(Amato_melanoma, sep = "\t", file = "scRNA/RData/Melanoma_bulk/ESTIMATE/inputEstimate_Amato.txt", quote = F)
filterCommonGenes(input.f = "scRNA/RData/Melanoma_bulk/ESTIMATE/inputEstimate_Amato.txt", output.f = "scRNA/RData/Melanoma_bulk/ESTIMATE/estimate_Amato.gct", id = "GeneSymbol")
estimateScore("scRNA/RData/Melanoma_bulk/ESTIMATE/estimate_Amato.gct", output.ds = "scRNA/RData/Melanoma_bulk/ESTIMATE/estimateScore_Amato.gct", platform  = "illumina")

## Loading estimate results
tableScores <-read.table(file = "scRNA/RData/Melanoma_bulk/ESTIMATE/estimateScore_Amato.gct", quote="", header = TRUE, sep = "\t", stringsAsFactors = FALSE, skip = 2)
rownames(tableScores) <- tableScores$NAME
tableScores$NAME <- NULL
tableScores$Description <- NULL
tableScores <- t(as.matrix(tableScores))
identical(rownames(tableScores), pData_Amato_melanoma$sample_id)

pData_Amato_melanoma = cbind(pData_Amato_melanoma, tableScores)
pData_Amato_melanoma$purity = cos(0.6049872018 + 0.0001467884 * pData_Amato_melanoma$ESTIMATEScore)
save(pData_Amato_melanoma, file = 'Melanoma_Bulk/GSE145996/pData_Amato.RData')

# Plot ----
library(readr)
load('Melanoma_Bulk/GSE145996/pData_Amato.RData')
clinical = pData_Amato_melanoma
rm(pData_Amato_melanoma)
CIBERSORTx_Amato_Results <- read_delim("Melanoma_Bulk/CIBERSORTx_Amato_Results.txt", 
                                       delim = "\t", escape_double = FALSE, 
                                       trim_ws = TRUE)

mat = as.data.frame(CIBERSORTx_Amato_Results[, 1:8])
rownames(mat) = mat$Mixture
rownames(clinical) = clinical$sample_id

mat = mat[rownames(clinical), ]
identical(rownames(clinical), rownames(mat))

colnames(clinical)
mat = cbind(mat, clinical[, c("sample_id", "response_NR", 'purity')])

for(barcode in mat$sample_id){
  message(barcode)
  ppurity = mat[mat$sample_id %in% barcode, 'purity']
  tmp2 <- as.matrix(mat[mat$sample_id %in% barcode, 2:8])
  tmp2 = tmp2*ppurity
  mat[mat$sample_id %in% barcode, 2:8] = tmp2
}
rowSums(mat[, 2:8])
# mat$TME = 1-mat$purity

# llevels = mat$Mixture
mat$Mixture = NULL
mat$purity = NULL
mat$sample_id = NULL

mmat = reshape2::melt(mat, c("response_NR"))
table(mmat$response_NR)
mmat$response = factor(mmat$response_NR, levels = c("R", "N"), labels = c('R', 'NR'))

p6 = ggplot(data = mmat, aes(x = response, y = value, fill = response)) + 
  # geom_violin(alpha=0.7) +
  geom_boxplot(width=0.3,  alpha = 0.6, position = position_dodge(width = 1)) +
  theme_classic(base_size = 14) + 
  theme(
    axis.text.y = element_text(size = 10),         # Show y-axis text
    axis.text.x = element_text(size = 10),         # Show x-axis text
    axis.title.x = element_blank(),                # Remove x-axis title
    # axis.title.y = element_blank(),                # Remove y-axis title
    strip.background = element_blank(),            # Remove strip background
    # strip.text.y = element_blank(),                # Remove strip text
    legend.position = "none",                       # Adjust legend position
    # panel.border = element_rect(colour = "black", size = 0.5), # Add a border to the plot
    panel.spacing = unit(0.05, "lines"),
    plot.title = element_text(hjust = 0.5, size = 14), # Center title
    plot.subtitle = element_text(hjust = 0.5, size = 10) # Center subtitle
  ) + 
  facet_grid(~variable, scales = 'free_x') +
  stat_compare_means(method = "wilcox.test", label = "p.format", label.x.npc = 'center', label.y = max(as.numeric(mmat$value), na.rm = TRUE) * 1.05) +
  scale_fill_manual(
    values = c('NR' = '#5F9EA0', 'R' = '#FA8072'),
    labels = c('NR', 'R')) +
  labs(
    title = "Amato et al. 2020",
    # subtitle = "Pozniak et al.",
    x = "",
    y = "Cell fraction",
    fill = "Response"
  )
p6
table(clinical$Therapy)

# cbioportal (Snyder et al. 2014, N Engl J Med) ----
load("/home/noviello/EPICA_Analysis/NIBITM4_scRNAseq/figure_paper/dataset_bulk_ICI.RDa", verbose = TRUE)
rm(RNA_seq_NORM, RNA_seq_noNORM)

MSKCC_RNA_Seq_expression_median <- read.delim("/home/noviello/EPICA_Analysis/NIBITM4_paper/Figures/dataset/MSKCC_RNA_Seq_expression_median.txt")

MSKCC_RNA_Seq_expression_median <- MSKCC_RNA_Seq_expression_median[,c("Entrez_Gene_Id", intersect(colnames(MSKCC_RNA_Seq_expression_median),clinical$Patient.ID))]
range(MSKCC_RNA_Seq_expression_median[,2:ncol(MSKCC_RNA_Seq_expression_median)])
MSKCC_RNA_Seq_expression_median[1:5, 1:5]
MSKCC_RNA_Seq_expression_median = MSKCC_RNA_Seq_expression_median[!duplicated(MSKCC_RNA_Seq_expression_median$Entrez_Gene_Id), ]

symbols = mapIds(org.Hs.eg.db, keys = as.character(MSKCC_RNA_Seq_expression_median$Entrez_Gene_Id), keytype = "ENTREZID", column = "SYMBOL")
symbols = symbols[!is.na(symbols)]
symbols = symbols[!duplicated(symbols)]

rownames(MSKCC_RNA_Seq_expression_median) = MSKCC_RNA_Seq_expression_median$Entrez_Gene_Id
MSKCC_RNA_Seq_expression_median = MSKCC_RNA_Seq_expression_median[names(symbols), ]
rownames(MSKCC_RNA_Seq_expression_median) = symbols
MSKCC_RNA_Seq_expression_median$Entrez_Gene_Id = NULL
range(MSKCC_RNA_Seq_expression_median)

MSKCC_RNA_Seq_expression_median = MSKCC_RNA_Seq_expression_median[-which(rowSums(is.na(MSKCC_RNA_Seq_expression_median)) != 0),]
range(MSKCC_RNA_Seq_expression_median)

MSKCC_RNA_Seq_expression_median = as.matrix(MSKCC_RNA_Seq_expression_median)
MSKCC_RNA_Seq_expression_median = na.omit(MSKCC_RNA_Seq_expression_median)
range(MSKCC_RNA_Seq_expression_median)
boxplot(log2(MSKCC_RNA_Seq_expression_median+1))
MSKCC_RNA_Seq_expression_median = t(quantileNormalization(t(MSKCC_RNA_Seq_expression_median)))
boxplot(log2(MSKCC_RNA_Seq_expression_median+1))


df_bulk <- data.frame(GeneSymbol = rownames(MSKCC_RNA_Seq_expression_median), MSKCC_RNA_Seq_expression_median)
df_bulk <- rbind(c("GeneSymbol", colnames(MSKCC_RNA_Seq_expression_median)), df_bulk)
df_bulk[1:10, 1:5]

# dir.create('scRNA/RData/Melanoma_bulk')
write.table(df_bulk, quote = FALSE, row.names = F, col.names = F, sep = "\t", 
            file = "scRNA/RData/Melanoma_bulk/RNAbulk_Snyder_CIBERSORT.txt")

MSKCC_RNA_Seq_expression_median = log2(MSKCC_RNA_Seq_expression_median + 1)
range(MSKCC_RNA_Seq_expression_median)

save(MSKCC_RNA_Seq_expression_median, file = 'scRNA/RData/Melanoma_bulk/Matrices/Snyder.RData')

## Loading estimate results
colnames(MSKCC_RNA_Seq_expression_median)
pData_Snyder_melanoma = clinical[clinical$Sample.ID %in% colnames(MSKCC_RNA_Seq_expression_median), ]
pData_Snyder_melanoma$purity = cos(0.6049872018 + 0.0001467884 * pData_Snyder_melanoma$ESTIMATE.Score)
save(pData_Snyder_melanoma, file = 'Melanoma_Bulk/pData_Snyder.RData')
# Plot ----
library(readr)
load('Melanoma_Bulk/pData_Snyder.RData')
clinical = pData_Snyder_melanoma
rm(pData_Snyder_melanoma)
CIBERSORTx_Snyder_Results <- read_delim("Melanoma_Bulk/CIBERSORTx_Snyder_Results.txt", 
                                        delim = "\t", escape_double = FALSE, 
                                        trim_ws = TRUE)

mat = as.data.frame(CIBERSORTx_Snyder_Results[, 1:8])
rownames(mat) = mat$Mixture
rownames(clinical) = clinical$Sample.ID

mat = mat[rownames(clinical), ]
identical(rownames(clinical), rownames(mat))

colnames(clinical)
mat = cbind(mat, clinical[, c("Sample.ID", "Responder", 'purity')])

for(barcode in mat$Sample.ID){
  message(barcode)
  ppurity = mat[mat$Sample.ID %in% barcode, 'purity']
  tmp2 <- as.matrix(mat[mat$Sample.ID %in% barcode, 2:8])
  tmp2 = tmp2*ppurity
  mat[mat$Sample.ID %in% barcode, 2:8] = tmp2
}
rowSums(mat[, 2:8])
# mat$TME = 1-mat$purity

# llevels = mat$Mixture
mat$Mixture = NULL
mat$purity = NULL
mat$Sample.ID = NULL

mmat = reshape2::melt(mat, c("Responder"))
table(mmat$Responder)
mmat$response = factor(mmat$Responder, levels = c("R", "NR"), labels = c('R', 'NR'))

p7 = ggplot(data = mmat, aes(x = response, y = value, fill = response)) + 
  # geom_violin(alpha=0.7) +
  geom_boxplot(width=0.3,  alpha = 0.6, position = position_dodge(width = 1)) +
  theme_classic(base_size = 14) + 
  theme(
    axis.text.y = element_text(size = 10),         # Show y-axis text
    axis.text.x = element_text(size = 10),         # Show x-axis text
    axis.title.x = element_blank(),                # Remove x-axis title
    # axis.title.y = element_blank(),                # Remove y-axis title
    strip.background = element_blank(),            # Remove strip background
    # strip.text.y = element_blank(),                # Remove strip text
    legend.position = "none",                       # Adjust legend position
    # panel.border = element_rect(colour = "black", size = 0.5), # Add a border to the plot
    panel.spacing = unit(0.05, "lines"),
    plot.title = element_text(hjust = 0.5, size = 14), # Center title
    plot.subtitle = element_text(hjust = 0.5, size = 10) # Center subtitle
  ) + 
  facet_grid(~variable, scales = 'free_x') +
  stat_compare_means(method = "wilcox.test", label = "p.format", label.x.npc = 'center', label.y = max(as.numeric(mmat$value), na.rm = TRUE) * 1.05) +
  scale_fill_manual(
    values = c('NR' = '#5F9EA0', 'R' = '#FA8072'),
    labels = c('NR', 'R')) +
  labs(
    title = "Snyder et al. 2014",
    # subtitle = "Pozniak et al.",
    x = "",
    y = "Cell fraction",
    fill = "Response"
  )
p7


# cbioportal (Van Allen et al. 2015, Science) ----
load("/home/noviello/EPICA_Analysis/NIBITM4_scRNAseq/figure_paper/dataset_bulk_ICI.RDa", verbose = TRUE)
rm(RNA_seq_NORM, RNA_seq_noNORM)

DFCI_RNA_Seq_expression_median <- read.delim("/home/noviello/EPICA_Analysis/NIBITM4_paper/Figures/dataset/DFCI_RNA_Seq_expression_median.txt")

DFCI_RNA_Seq_expression_median <- DFCI_RNA_Seq_expression_median[,c("Entrez_Gene_Id", intersect(colnames(DFCI_RNA_Seq_expression_median),clinical$Patient.ID))]
range(DFCI_RNA_Seq_expression_median[,2:ncol(DFCI_RNA_Seq_expression_median)])
DFCI_RNA_Seq_expression_median[1:5, 1:5]
DFCI_RNA_Seq_expression_median = DFCI_RNA_Seq_expression_median[!duplicated(DFCI_RNA_Seq_expression_median$Entrez_Gene_Id), ]

require(org.Hs.eg.db)
symbols = mapIds(org.Hs.eg.db, keys = as.character(DFCI_RNA_Seq_expression_median$Entrez_Gene_Id), keytype = "ENTREZID", column = "SYMBOL")
symbols = symbols[!is.na(symbols)]
symbols = symbols[!duplicated(symbols)]

rownames(DFCI_RNA_Seq_expression_median) = DFCI_RNA_Seq_expression_median$Entrez_Gene_Id
DFCI_RNA_Seq_expression_median = DFCI_RNA_Seq_expression_median[names(symbols), ]
rownames(DFCI_RNA_Seq_expression_median) = symbols
DFCI_RNA_Seq_expression_median$Entrez_Gene_Id = NULL
range(DFCI_RNA_Seq_expression_median)

DFCI_RNA_Seq_expression_median = DFCI_RNA_Seq_expression_median[-which(rowSums(is.na(DFCI_RNA_Seq_expression_median)) != 0),]
range(DFCI_RNA_Seq_expression_median)

DFCI_RNA_Seq_expression_median = as.matrix(DFCI_RNA_Seq_expression_median)
DFCI_RNA_Seq_expression_median = na.omit(DFCI_RNA_Seq_expression_median)
range(DFCI_RNA_Seq_expression_median)
boxplot(log2(DFCI_RNA_Seq_expression_median+1))
DFCI_RNA_Seq_expression_median = t(quantileNormalization(t(DFCI_RNA_Seq_expression_median)))
boxplot(log2(DFCI_RNA_Seq_expression_median+1))


df_bulk <- data.frame(GeneSymbol = rownames(DFCI_RNA_Seq_expression_median), DFCI_RNA_Seq_expression_median)
df_bulk <- rbind(c("GeneSymbol", colnames(DFCI_RNA_Seq_expression_median)), df_bulk)
df_bulk[1:10, 1:5]

# dir.create('scRNA/RData/Melanoma_bulk')
write.table(df_bulk, quote = FALSE, row.names = F, col.names = F, sep = "\t", 
            file = "scRNA/RData/Melanoma_bulk/RNAbulk_VanAllen_CIBERSORT.txt")

DFCI_RNA_Seq_expression_median = log2(DFCI_RNA_Seq_expression_median + 1)
range(DFCI_RNA_Seq_expression_median)

save(DFCI_RNA_Seq_expression_median, file = 'scRNA/RData/Melanoma_bulk/Matrices/VanAllen.RData')

## Loading estimate results
colnames(DFCI_RNA_Seq_expression_median)
pData_VanAllen_melanoma = clinical[clinical$Sample.ID %in% colnames(DFCI_RNA_Seq_expression_median), ]
pData_VanAllen_melanoma$purity = cos(0.6049872018 + 0.0001467884 * pData_VanAllen_melanoma$ESTIMATE.Score)
save(pData_VanAllen_melanoma, file = 'Melanoma_Bulk/pData_VanAllen.RData')

# Plot ----
library(readr)
load('Melanoma_Bulk/pData_VanAllen.RData')
clinical = pData_VanAllen_melanoma
rm(pData_VanAllen_melanoma)
CIBERSORTx_VanAllen_Results <- read_delim("Melanoma_Bulk/CIBERSORTx_VanAllen_Results.txt", 
                                          delim = "\t", escape_double = FALSE, 
                                          trim_ws = TRUE)

mat = as.data.frame(CIBERSORTx_VanAllen_Results[, 1:8])
rownames(mat) = mat$Mixture
rownames(clinical) = clinical$Sample.ID

mat = mat[rownames(clinical), ]
identical(rownames(clinical), rownames(mat))

colnames(clinical)
mat = cbind(mat, clinical[, c("Sample.ID", "Responder", 'purity')])

for(barcode in mat$Sample.ID){
  message(barcode)
  ppurity = mat[mat$Sample.ID %in% barcode, 'purity']
  tmp2 <- as.matrix(mat[mat$Sample.ID %in% barcode, 2:8])
  tmp2 = tmp2*ppurity
  mat[mat$Sample.ID %in% barcode, 2:8] = tmp2
}
rowSums(mat[, 2:8])
# mat$TME = 1-mat$purity

# llevels = mat$Mixture
mat$Mixture = NULL
mat$purity = NULL
mat$Sample.ID = NULL

mmat = reshape2::melt(mat, c("Responder"))
table(mmat$Responder)
mmat$response = factor(mmat$Responder, levels = c("R", "NR"), labels = c('R', 'NR'))

p7 = ggplot(data = mmat, aes(x = response, y = value, fill = response)) + 
  # geom_violin(alpha=0.7) +
  geom_boxplot(width=0.3,  alpha = 0.6, position = position_dodge(width = 1)) +
  theme_classic(base_size = 14) + 
  theme(
    axis.text.y = element_text(size = 10),         # Show y-axis text
    axis.text.x = element_text(size = 10),         # Show x-axis text
    axis.title.x = element_blank(),                # Remove x-axis title
    # axis.title.y = element_blank(),                # Remove y-axis title
    strip.background = element_blank(),            # Remove strip background
    # strip.text.y = element_blank(),                # Remove strip text
    legend.position = "none",                       # Adjust legend position
    # panel.border = element_rect(colour = "black", size = 0.5), # Add a border to the plot
    panel.spacing = unit(0.05, "lines"),
    plot.title = element_text(hjust = 0.5, size = 14), # Center title
    plot.subtitle = element_text(hjust = 0.5, size = 10) # Center subtitle
  ) + 
  facet_grid(~variable, scales = 'free_x') +
  stat_compare_means(method = "wilcox.test", label = "p.format", label.x.npc = 'center', label.y = max(as.numeric(mmat$value), na.rm = TRUE) * 1.05) +
  scale_fill_manual(
    values = c('NR' = '#5F9EA0', 'R' = '#FA8072'),
    labels = c('NR', 'R')) +
  labs(
    title = "Van Allen et al. 2017",
    # subtitle = "Pozniak et al.",
    x = "",
    y = "Cell fraction",
    fill = "Response"
  )
p7

# GSE168204 (Du et al. 2021, Nat. Commun) ----
Du_melanoma = as.data.frame(read_delim("Melanoma_Bulk/GSE168204/GSE168204_norm_counts_TPM_GRCh38.p13_NCBI.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE))
annot = as.data.frame(read_delim("Melanoma_Bulk/GSE168204/Human.GRCh38.p13.annot.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE))
Du_melanoma = merge(Du_melanoma, annot[, c('GeneID', 'Symbol')], by = 'GeneID')
Du_melanoma = Du_melanoma[!duplicated(Du_melanoma$Symbol), ]
rownames(Du_melanoma) = Du_melanoma$Symbol
Du_melanoma$GeneID = NULL
Du_melanoma$Symbol = NULL
# Du_melanoma = as.matrix(log2(Du_melanoma + 1))

pData_Du_melanoma = as.data.frame(readxl::read_excel("Melanoma_Bulk/GSE168204/pData.xlsx"))
rownames(pData_Du_melanoma) = pData_Du_melanoma$Accession
pData_Du_melanoma = pData_Du_melanoma[colnames(Du_melanoma), ]

Du_melanoma[1:5, 1:5]

Du_melanoma = as.matrix(Du_melanoma)
Du_melanoma = na.omit(Du_melanoma)
range(Du_melanoma)
boxplot(log2(Du_melanoma+1))
Du_melanoma = t(quantileNormalization(t(Du_melanoma)))
boxplot(log2(Du_melanoma+1))


df_bulk <- data.frame(GeneSymbol = rownames(Du_melanoma), Du_melanoma)
df_bulk <- rbind(c("GeneSymbol", colnames(Du_melanoma)), df_bulk)
df_bulk[1:10, 1:5]

# dir.create('scRNA/RData/Melanoma_bulk')
write.table(df_bulk, quote = FALSE, row.names = F, col.names = F, sep = "\t", 
            file = "scRNA/RData/Melanoma_bulk/RNAbulk_Du_CIBERSORT.txt")

Du_melanoma = log2(Du_melanoma + 1)
range(Du_melanoma)

save(Du_melanoma, file = 'scRNA/RData/Melanoma_bulk/Matrices/Du.RData')

# dir.create("scRNA/RData/Melanoma_bulk/ESTIMATE")
Du_melanoma <- Du_melanoma[, pData_Du_melanoma$Accession]

write.table(Du_melanoma, sep = "\t", file = "scRNA/RData/Melanoma_bulk/ESTIMATE/inputEstimate_Du.txt", quote = F)
filterCommonGenes(input.f = "scRNA/RData/Melanoma_bulk/ESTIMATE/inputEstimate_Du.txt", output.f = "scRNA/RData/Melanoma_bulk/ESTIMATE/estimate_Du.gct", id = "GeneSymbol")
estimateScore("scRNA/RData/Melanoma_bulk/ESTIMATE/estimate_Du.gct", output.ds = "scRNA/RData/Melanoma_bulk/ESTIMATE/estimateScore_Du.gct", platform  = "illumina")

## Loading estimate results
tableScores <-read.table(file = "scRNA/RData/Melanoma_bulk/ESTIMATE/estimateScore_Du.gct", quote="", header = TRUE, sep = "\t", stringsAsFactors = FALSE, skip = 2)
rownames(tableScores) <- tableScores$NAME
tableScores$NAME <- NULL
tableScores$Description <- NULL
tableScores <- t(as.matrix(tableScores))
identical(rownames(tableScores), pData_Du_melanoma$Accession)

pData_Du_melanoma = cbind(pData_Du_melanoma, tableScores)
pData_Du_melanoma$purity = cos(0.6049872018 + 0.0001467884 * pData_Du_melanoma$ESTIMATEScore)
colnames(pData_Du_melanoma)[colnames(pData_Du_melanoma) == 'Antibody'] = 'Therapy'
colnames(pData_Du_melanoma)[colnames(pData_Du_melanoma) == 'Response immune checkpoint blockade therapy'] = 'Responder_NR'

save(pData_Du_melanoma, file = 'Melanoma_Bulk/GSE168204/pData_Du.RData')

# Plot ----
library(readr)
load('Melanoma_Bulk/pData_Snyder.RData')
clinical = pData_Snyder_melanoma
rm(pData_Snyder_melanoma)
CIBERSORTx_Snyder_Results <- read_delim("Melanoma_Bulk/CIBERSORTx_Snyder_Results.txt", 
                                        delim = "\t", escape_double = FALSE, 
                                        trim_ws = TRUE)

mat = as.data.frame(CIBERSORTx_Snyder_Results[, 1:8])
rownames(mat) = mat$Mixture
rownames(clinical) = clinical$Sample.ID

mat = mat[rownames(clinical), ]
identical(rownames(clinical), rownames(mat))

colnames(clinical)
mat = cbind(mat, clinical[, c("Sample.ID", "Responder", 'purity')])

for(barcode in mat$Sample.ID){
  message(barcode)
  ppurity = mat[mat$Sample.ID %in% barcode, 'purity']
  tmp2 <- as.matrix(mat[mat$Sample.ID %in% barcode, 2:8])
  tmp2 = tmp2*ppurity
  mat[mat$Sample.ID %in% barcode, 2:8] = tmp2
}
rowSums(mat[, 2:8])
# mat$TME = 1-mat$purity

# llevels = mat$Mixture
mat$Mixture = NULL
mat$purity = NULL
mat$Sample.ID = NULL

mmat = reshape2::melt(mat, c("Responder"))
table(mmat$Responder)
mmat$response = factor(mmat$Responder, levels = c("R", "NR"), labels = c('R', 'NR'))

p7 = ggplot(data = mmat, aes(x = response, y = value, fill = response)) + 
  # geom_violin(alpha=0.7) +
  geom_boxplot(width=0.3,  alpha = 0.6, position = position_dodge(width = 1)) +
  theme_classic(base_size = 14) + 
  theme(
    axis.text.y = element_text(size = 10),         # Show y-axis text
    axis.text.x = element_text(size = 10),         # Show x-axis text
    axis.title.x = element_blank(),                # Remove x-axis title
    # axis.title.y = element_blank(),                # Remove y-axis title
    strip.background = element_blank(),            # Remove strip background
    # strip.text.y = element_blank(),                # Remove strip text
    legend.position = "none",                       # Adjust legend position
    # panel.border = element_rect(colour = "black", size = 0.5), # Add a border to the plot
    panel.spacing = unit(0.05, "lines"),
    plot.title = element_text(hjust = 0.5, size = 14), # Center title
    plot.subtitle = element_text(hjust = 0.5, size = 10) # Center subtitle
  ) + 
  facet_grid(~variable, scales = 'free_x') +
  stat_compare_means(method = "wilcox.test", label = "p.format", label.x.npc = 'center', label.y = max(as.numeric(mmat$value), na.rm = TRUE) * 1.05) +
  scale_fill_manual(
    values = c('NR' = '#5F9EA0', 'R' = '#FA8072'),
    labels = c('NR', 'R')) +
  labs(
    title = "Snyder et al. 2014",
    # subtitle = "Pozniak et al.",
    x = "",
    y = "Cell fraction",
    fill = "Response"
  )
p7
# All dataset together ----
# Hugo
# Checked GEO, tiger and original paper: OK!
load('Melanoma_Bulk/GSE78220/pData_Hugo.RData')
clinical = pData_Hugo_melanoma
rm(pData_Hugo_melanoma)
CIBERSORTx_Hugo_Results <- read_delim("Melanoma_Bulk/CIBERSORTx_Hugo_Results.txt", 
                                      delim = "\t", escape_double = FALSE, 
                                      trim_ws = TRUE)

mat = as.data.frame(CIBERSORTx_Hugo_Results[, 1:8])
rownames(mat) = mat$Mixture
rownames(clinical) = clinical$sample_id

mat = mat[rownames(clinical), ]
identical(rownames(clinical), rownames(mat))

colnames(clinical)
mat = cbind(mat, clinical[, c("sample_id", "response_NR", 'purity', 'Treatment', "overall survival (days)", "vital status", 'Therapy', 'patient_name', 'Total Mutation')])

for(barcode in mat$sample_id){
  message(barcode)
  ppurity = mat[mat$sample_id %in% barcode, 'purity']
  tmp2 <- as.matrix(mat[mat$sample_id %in% barcode, 2:8])
  tmp2 = tmp2*ppurity
  mat[mat$sample_id %in% barcode, 2:8] = tmp2
}
rowSums(mat[, 2:8])
# mat$TME = 1-mat$purity

# llevels = mat$Mixture
mat$Mixture = NULL
# mat$purity = NULL
# rownames(mat) = mat$sample_id
mat$sample_id = NULL
hugo = mat

# Riaz
# Checked GEO, tiger and original paper: not ok; combining the info
load('Melanoma_Bulk/GSE91061/pData_Riaz.RData')
clinical = pData_Riaz_melanoma
# adding the response for some of the samples uknown
clinical$response[clinical$patient_name %in% paste0('Pt', c(81,80,109))] = 'PD'
clinical$response_NR[clinical$patient_name %in% paste0('Pt', c(81,80,109))] = 'N'
clinical$response[clinical$patient_name %in% paste0('Pt', c(105,69,35))] = 'PRCR'
clinical$response_NR[clinical$patient_name %in% paste0('Pt', c(105,69,35))] = 'R'

rm(pData_Riaz_melanoma)
CIBERSORTx_Riaz_Results <- read_delim("Melanoma_Bulk/CIBERSORTx_Riaz_Results.txt", 
                                      delim = "\t", escape_double = FALSE, 
                                      trim_ws = TRUE)

mat = as.data.frame(CIBERSORTx_Riaz_Results[, 1:8])
rownames(mat) = mat$Mixture
rownames(clinical) = clinical$sample_id

mat = mat[rownames(clinical), ]
identical(rownames(clinical), rownames(mat))

colnames(clinical)
mat = cbind(mat, clinical[, c("sample_id", "response_NR", 'purity', 'Treatment', "overall survival (days)", "vital status", 'Therapy', 'patient_name', 'Total Mutation')])

for(barcode in mat$sample_id){
  message(barcode)
  ppurity = mat[mat$sample_id %in% barcode, 'purity']
  tmp2 <- as.matrix(mat[mat$sample_id %in% barcode, 2:8])
  tmp2 = tmp2*ppurity
  mat[mat$sample_id %in% barcode, 2:8] = tmp2
}
rowSums(mat[, 2:8])
# mat$TME = 1-mat$purity

# llevels = mat$Mixture
mat$Mixture = NULL
# mat$purity = NULL
# rownames(mat) = mat$sample_id
mat$sample_id = NULL
riaz = mat

# Nathanson
# Checked tiger and original paper: OK!
load('Melanoma_Bulk/Nathanson_2017/pData_Nathanson.RData')
clinical = pData_Nathanson_melanoma
rm(pData_Nathanson_melanoma)
CIBERSORTx_Nathanson_Results <- read_delim("Melanoma_Bulk/CIBERSORTx_Nathanson_Results.txt", 
                                           delim = "\t", escape_double = FALSE, 
                                           trim_ws = TRUE)

mat = as.data.frame(CIBERSORTx_Nathanson_Results[, 1:8])
rownames(mat) = mat$Mixture
rownames(clinical) = clinical$sample_id

mat = mat[rownames(clinical), ]
identical(rownames(clinical), rownames(mat))

colnames(clinical)
mat = cbind(mat, clinical[, c("sample_id", "response_NR", 'purity', 'Treatment', "overall survival (days)", "vital status", 'Therapy', 'patient_name', 'Total Mutation')])

for(barcode in mat$sample_id){
  message(barcode)
  ppurity = mat[mat$sample_id %in% barcode, 'purity']
  tmp2 <- as.matrix(mat[mat$sample_id %in% barcode, 2:8])
  tmp2 = tmp2*ppurity
  mat[mat$sample_id %in% barcode, 2:8] = tmp2
}
rowSums(mat[, 2:8])
# mat$TME = 1-mat$purity

# llevels = mat$Mixture
mat$Mixture = NULL
# mat$purity = NULL
# rownames(mat) = mat$sample_id
mat$sample_id = NULL

nathanson = mat

# Gide
# Checked tiger and original paper: OK!
load('Melanoma_Bulk/PRJEB23709/pData_Gide.RData')
clinical = pData_Gide_melanoma
rm(pData_Gide_melanoma)
CIBERSORTx_Gide_Results <- read_delim("Melanoma_Bulk/CIBERSORTx_Gide_Results.txt", 
                                      delim = "\t", escape_double = FALSE, 
                                      trim_ws = TRUE)

mat = as.data.frame(CIBERSORTx_Gide_Results[, 1:8])
rownames(mat) = mat$Mixture
rownames(clinical) = clinical$sample_id

mat = mat[rownames(clinical), ]
identical(rownames(clinical), rownames(mat))

colnames(clinical)
mat = cbind(mat, clinical[, c("sample_id", "response_NR", 'purity', 'Treatment', "overall survival (days)", "vital status", 'Therapy', 'patient_name', 'Total Mutation')])

for(barcode in mat$sample_id){
  message(barcode)
  ppurity = mat[mat$sample_id %in% barcode, 'purity']
  tmp2 <- as.matrix(mat[mat$sample_id %in% barcode, 2:8])
  tmp2 = tmp2*ppurity
  mat[mat$sample_id %in% barcode, 2:8] = tmp2
}
rowSums(mat[, 2:8])
# mat$TME = 1-mat$purity

# llevels = mat$Mixture
mat$Mixture = NULL
# mat$purity = NULL
# rownames(mat) = mat$sample_id
mat$sample_id = NULL

gide = mat

# Liu
# Checked phs000452 and original paper: OK!
load('Melanoma_Bulk/phs000452/pData_Liu.RData')
clinical = pData_Liu_melanoma
rm(pData_Liu_melanoma)
CIBERSORTx_Liu_Results <- read_delim("Melanoma_Bulk/CIBERSORTx_Liu_Results.txt", 
                                     delim = "\t", escape_double = FALSE, 
                                     trim_ws = TRUE)

mat = as.data.frame(CIBERSORTx_Liu_Results[, 1:8])
rownames(mat) = mat$Mixture

clinical = clinical[grep('^P', clinical$patient_name), ]
rownames(clinical) = clinical$sample_id

clinical_paper = openxlsx::read.xlsx(xlsxFile = 'Melanoma_Bulk/phs000452/41591_2019_654_MOESM4_ESM.xlsx', sheet = 1)
rownames(clinical_paper) = c('N/A', clinical_paper$`Supplemental.Table.1:.Genomic.and.Clinical.Patient.and.Tumor.Characteristics`[2:nrow(clinical_paper)])
clinical_paper[1,1] = 'Patients'
colnames(clinical_paper) = clinical_paper[1, ]
clinical_paper = clinical_paper[2:nrow(clinical_paper), ]
clinical_paper = clinical_paper[grep('^P', clinical_paper$Patients), ]

clinical$patient_name = gsub('_T_M|_T_P|_T', '', clinical$patient_name)
clinical_paper = clinical_paper[clinical_paper$Patients %in% clinical$patient_name, ]
clinical_paper = clinical_paper[clinical$patient_name, ]
clinical_paper = clinical_paper[, c(1, 42:44)]
clinical_paper$Treatment = NA
clinical_paper$Treatment[clinical_paper$`biopsyContext (1=Pre-Ipi; 2=On-Ipi; 3=Pre-PD1; 4=On-PD1)` %in% 1] = 'PRE'
clinical_paper$Treatment[clinical_paper$`biopsyContext (1=Pre-Ipi; 2=On-Ipi; 3=Pre-PD1; 4=On-PD1)` %in% 2] = 'ON'
clinical_paper$Treatment[clinical_paper$`biopsyContext (1=Pre-Ipi; 2=On-Ipi; 3=Pre-PD1; 4=On-PD1)` %in% 3 & clinical_paper$daysBiopsyAfterIpiStart %in% c('na', 'unk')] = 'PRE'
clinical_paper$Treatment[clinical_paper$`biopsyContext (1=Pre-Ipi; 2=On-Ipi; 3=Pre-PD1; 4=On-PD1)` %in% 3 & as.numeric(clinical_paper$daysBiopsyAfterIpiStart) < 0] = 'PRE'
clinical_paper$Treatment[clinical_paper$`biopsyContext (1=Pre-Ipi; 2=On-Ipi; 3=Pre-PD1; 4=On-PD1)` %in% 3 & as.numeric(clinical_paper$daysBiopsyAfterIpiStart) > 0] = 'ON'
clinical_paper$Treatment[clinical_paper$`biopsyContext (1=Pre-Ipi; 2=On-Ipi; 3=Pre-PD1; 4=On-PD1)` %in% 4] = 'ON'

clinical_paper$Therapy = NA
clinical_paper$Therapy[clinical_paper$`biopsyContext (1=Pre-Ipi; 2=On-Ipi; 3=Pre-PD1; 4=On-PD1)` %in% 1] = 'anti-CTLA4'
clinical_paper$Therapy[clinical_paper$`biopsyContext (1=Pre-Ipi; 2=On-Ipi; 3=Pre-PD1; 4=On-PD1)` %in% 2] = 'anti-CTLA4'
clinical_paper$Therapy[clinical_paper$`biopsyContext (1=Pre-Ipi; 2=On-Ipi; 3=Pre-PD1; 4=On-PD1)` %in% 3 & clinical_paper$daysBiopsyAfterIpiStart %in% c('na', 'unk')] = 'anti-CTLA4+anti-PD1'
clinical_paper$Therapy[clinical_paper$`biopsyContext (1=Pre-Ipi; 2=On-Ipi; 3=Pre-PD1; 4=On-PD1)` %in% 3 & as.numeric(clinical_paper$daysBiopsyAfterIpiStart) < 0] = 'anti-CTLA4+anti-PD1'
clinical_paper$Therapy[clinical_paper$`biopsyContext (1=Pre-Ipi; 2=On-Ipi; 3=Pre-PD1; 4=On-PD1)` %in% 3 & as.numeric(clinical_paper$daysBiopsyAfterIpiStart) > 0] = 'anti-CTLA4'
clinical_paper$Therapy[clinical_paper$`biopsyContext (1=Pre-Ipi; 2=On-Ipi; 3=Pre-PD1; 4=On-PD1)` %in% 4] = 'anti-CTLA4+anti-PD1'

table(clinical_paper$Treatment, clinical_paper$Therapy, exclude = NULL)

colnames(clinical)
colnames(clinical_paper)[1] = 'patient_name'

clinical = merge(clinical[, c(1:4, 6:12, 14:21)], clinical_paper[, c(1, 5:6)], by = c('patient_name'))
rownames(clinical) = clinical$sample_id

mat = mat[rownames(clinical), ]
identical(rownames(clinical), rownames(mat))

colnames(clinical)
mat = cbind(mat, clinical[, c("sample_id", "response_NR", 'purity', 'Treatment', "overall survival (days)", "vital status", 'Therapy', 'patient_name', 'Total Mutation')])

for(barcode in mat$sample_id){
  message(barcode)
  ppurity = mat[mat$sample_id %in% barcode, 'purity']
  tmp2 <- as.matrix(mat[mat$sample_id %in% barcode, 2:8])
  tmp2 = tmp2*ppurity
  mat[mat$sample_id %in% barcode, 2:8] = tmp2
}
rowSums(mat[, 2:8])
# mat$TME = 1-mat$purity

# llevels = mat$Mixture
mat$Mixture = NULL
# mat$purity = NULL
# rownames(mat) = mat$sample_id
mat$sample_id = NULL

liu = mat

# Amato
load('Melanoma_Bulk/GSE145996/pData_Amato.RData')
clinical = pData_Amato_melanoma
rm(pData_Amato_melanoma)
CIBERSORTx_Amato_Results <- read_delim("Melanoma_Bulk/CIBERSORTx_Amato_Results.txt", 
                                       delim = "\t", escape_double = FALSE, 
                                       trim_ws = TRUE)

mat = as.data.frame(CIBERSORTx_Amato_Results[, 1:8])
rownames(mat) = mat$Mixture
rownames(clinical) = clinical$sample_id

mat = mat[rownames(clinical), ]
identical(rownames(clinical), rownames(mat))

colnames(clinical)
mat = cbind(mat, clinical[, c("sample_id", "response_NR", 'purity', 'Treatment', "overall survival (days)", "vital status", 'Therapy', 'patient_name', 'Total Mutation')])

for(barcode in mat$sample_id){
  message(barcode)
  ppurity = mat[mat$sample_id %in% barcode, 'purity']
  tmp2 <- as.matrix(mat[mat$sample_id %in% barcode, 2:8])
  tmp2 = tmp2*ppurity
  mat[mat$sample_id %in% barcode, 2:8] = tmp2
}
rowSums(mat[, 2:8])
# mat$TME = 1-mat$purity

# llevels = mat$Mixture
mat$Mixture = NULL
# mat$purity = NULL
# rownames(mat) = mat$sample_id
mat$sample_id = NULL

amato = mat

# Snyder
load('Melanoma_Bulk/pData_Snyder.RData')
clinical = pData_Snyder_melanoma
rm(pData_Snyder_melanoma)
CIBERSORTx_Snyder_Results <- read_delim("Melanoma_Bulk/CIBERSORTx_Snyder_Results.txt", 
                                        delim = "\t", escape_double = FALSE, 
                                        trim_ws = TRUE)

mat = as.data.frame(CIBERSORTx_Snyder_Results[, 1:8])
rownames(mat) = mat$Mixture
rownames(clinical) = clinical$Sample.ID

mat = mat[rownames(clinical), ]
identical(rownames(clinical), rownames(mat))

colnames(clinical)
clinical$Therapy = 'anti-CTLA4'
mat = cbind(mat, clinical[, c(3, 44, 9:10, 45, 46, 2, 6)])

for(barcode in mat$Sample.ID){
  message(barcode)
  ppurity = mat[mat$Sample.ID %in% barcode, 'purity']
  tmp2 <- as.matrix(mat[mat$Sample.ID %in% barcode, 2:8])
  tmp2 = tmp2*ppurity
  mat[mat$Sample.ID %in% barcode, 2:8] = tmp2
}
rowSums(mat[, 2:8])
# mat$TME = 1-mat$purity

# llevels = mat$Mixture
mat$Mixture = NULL
# mat$purity = NULL
# rownames(mat) = mat$Sample.ID
mat$Sample.ID = NULL

snyder = mat
snyder$Treatment = 'POST'
snyder$Treatment[snyder$Patient.ID %in% c('CR6126', 'LSD0167', 'NR3549', 'NR5784', 'SD1494', 'SD2056')] = 'PRE'

# Van Allen
load('Melanoma_Bulk/pData_VanAllen.RData')
clinical = pData_VanAllen_melanoma
rm(pData_VanAllen_melanoma)
CIBERSORTx_VanAllen_Results <- read_delim("Melanoma_Bulk/CIBERSORTx_VanAllen_Results.txt", 
                                          delim = "\t", escape_double = FALSE, 
                                          trim_ws = TRUE)

mat = as.data.frame(CIBERSORTx_VanAllen_Results[, 1:8])
rownames(mat) = mat$Mixture
rownames(clinical) = clinical$Sample.ID

mat = mat[rownames(clinical), ]
identical(rownames(clinical), rownames(mat))

colnames(clinical)
clinical$Therapy = 'anti-CTLA4'
mat = cbind(mat, clinical[, c(3, 44, 9:10, 45, 46, 2, 6)])

for(barcode in mat$Sample.ID){
  message(barcode)
  ppurity = mat[mat$Sample.ID %in% barcode, 'purity']
  tmp2 <- as.matrix(mat[mat$Sample.ID %in% barcode, 2:8])
  tmp2 = tmp2*ppurity
  mat[mat$Sample.ID %in% barcode, 2:8] = tmp2
}
rowSums(mat[, 2:8])
# mat$TME = 1-mat$purity

# llevels = mat$Mixture
mat$Mixture = NULL
# mat$purity = NULL
# rownames(mat) = mat$Sample.ID
mat$Sample.ID = NULL

vanallen = mat

# Du
load('Melanoma_Bulk/GSE168204/pData_Du.RData')
clinical = pData_Du_melanoma
rm(pData_Du_melanoma)
CIBERSORTx_Du_Results <- read_delim("Melanoma_Bulk/CIBERSORTx_Du_Results.txt", 
                                    delim = "\t", escape_double = FALSE, 
                                    trim_ws = TRUE)

mat = as.data.frame(CIBERSORTx_Du_Results[, 1:8])
rownames(mat) = mat$Mixture
rownames(clinical) = clinical$Accession

mat = mat[rownames(clinical), ]
identical(rownames(clinical), rownames(mat))

colnames(clinical)
mat = cbind(mat, clinical[, c("Accession", "Responder_NR", 'Treatment state', 'purity', 'Therapy')])
colnames(mat)[9:10] = c('sample_id', 'response_NR')

for(barcode in mat$sample_id){
  message(barcode)
  ppurity = mat[mat$sample_id %in% barcode, 'purity']
  tmp2 <- as.matrix(mat[mat$sample_id %in% barcode, 2:8])
  tmp2 = tmp2*ppurity
  mat[mat$sample_id %in% barcode, 2:8] = tmp2
}

rowSums(mat[, 2:8])
# mat$TME = 1-mat$purity

# llevels = mat$Mixture
mat$Mixture = NULL
# mat$purity = NULL
# rownames(mat) = mat$sample_id
mat$sample_id = NULL
du = mat

# All dataset together
colnames(hugo)
colnames(riaz)
colnames(gide)
colnames(nathanson)
colnames(liu)
colnames(amato)
colnames(snyder)
colnames(snyder)[8] = 'response_NR'
snyder$Treatment = rep(NA, nrow(snyder))
snyder$Overall.Survival..Months. = floor(snyder$Overall.Survival..Months. * 30.4375)
colnames(snyder)[9] = 'overall survival (days)'
colnames(snyder)[10] = 'vital status'
colnames(snyder)[13] = 'patient_name'
colnames(snyder)[14] = 'Total Mutation'
snyder = snyder[, colnames(amato)]
colnames(vanallen)
colnames(vanallen)[8] = 'response_NR'
vanallen$Treatment = rep('PRE', nrow(vanallen))
vanallen$Overall.Survival..Months. = floor(vanallen$Overall.Survival..Months. * 30.4375)
colnames(vanallen)[9] = 'overall survival (days)'
colnames(vanallen)[10] = 'vital status'
colnames(vanallen)[13] = 'patient_name'
colnames(vanallen)[14] = 'Total Mutation'
vanallen = vanallen[, colnames(snyder)]
colnames(du)
colnames(du)[9] = 'Treatment'
du$`overall survival (days)` = rep(NA, nrow(du))
du$`vital status` = rep(NA, nrow(du))
du$patient_name = rep(NA, nrow(du))
du$`Total Mutation` = rep(NA, nrow(du))
du = du[, colnames(snyder)]

hugo$Dataset = rep('Hugo', nrow(hugo))
riaz$Dataset = rep('riaz', nrow(riaz))
nathanson$Dataset = rep('nathanson', nrow(nathanson))
gide$Dataset = rep('gide', nrow(gide))
liu$Dataset = rep('liu', nrow(liu))
amato$Dataset = rep('amato', nrow(amato))
snyder$Dataset = rep('snyder', nrow(snyder))
vanallen$Dataset = rep('vanallen', nrow(vanallen))
du$Dataset = rep('du', nrow(du))

mmat = rbind(hugo, riaz, nathanson, gide, liu, amato, snyder, vanallen, du)
table(mmat$response_NR, mmat$Therapy)
mmat$Therapy[mmat$Therapy %in% 'anti-CTLA-4'] = 'anti-CTLA4'
mmat$Therapy[mmat$Therapy %in% 'anti-PD-1'] = 'anti-PD1'
mmat$Therapy[mmat$Therapy %in% 'anti-CTLA-4+anti-PD-1'] = 'anti-CTLA4+anti-PD1'
mmat$Therapy[mmat$Therapy %in% 'anti-PD1+anti-CTLA4'] = 'anti-CTLA4+anti-PD1'
mmat$response_NR[mmat$response_NR %in% 'NR'] = 'N'
mmat = mmat[!mmat$response_NR %in% 'UNK', ]
mmat = mmat[!mmat$Therapy %in% c('anti-PDL1'), ]
table(mmat$response_NR, mmat$Therapy, exclude = NULL)

save(mmat, file = 'Melanoma_Bulk/dataframe_bulk_all.RData')