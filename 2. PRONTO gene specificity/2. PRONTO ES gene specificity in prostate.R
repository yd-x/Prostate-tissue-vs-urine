# This script is for looking for the specificity in prostate and prostate cancer
# in HPA database. The normal tissue data is the "RNA consensus tissue gene data"
# generated from HPA and GTEx through HPA group's own normalization pipeline. 
# The cancer data is from TCGA data downloaded from HPA website.

# Package preparation -----------------------------------------------------

setwd("G:/Yuandong Xing/OneDrive - Queen's University/project codes and results/2. PRONTO gene specificity")
setwd("E:/study/Queens/OneDrive - Queen's University/project codes and results/2. PRONTO gene specificity")

library(rio)
library(biomaRt)
library(openxlsx)
library(dplyr)
library(Vennerable)
library(reshape2)
library(stringr)
library(ggpubr)

load("PRONTO ES gene specificity in prostate.RData")


# Import and extract the data ---------------------------------------------

## PRONTO Epi and Str genes

# Import genes and protein expression check results
PRONTO.epi <- read.xlsx(xlsxFile = "PRONTO Gene Classification.xlsx",
                        sheet = "Epi in IHC images")
PRONTO.str <- read.xlsx(xlsxFile = "PRONTO Gene Classification.xlsx",
                        sheet = "Str in IHC images")
PRONTO.epi <- PRONTO.epi[
  ,c('Gene','ID','Normal.Prostate.Specific.protein',
     'Prostate.cancer.specific.protein.expression')]
PRONTO.str <- PRONTO.str[
  ,c('Gene','ID','Normal.Prostate.Specific.protein',
     'Prostate.cancer.specific.protein.expression')]

# Convert gene symbols to Ensembl IDs
# ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
# PRONTO.epi.ensemble <- 
#   getBM(mart = ensembl,filters = "entrezgene_id",values = PRONTO.epi$ID,
#         attributes = c("ensembl_gene_id","entrezgene_id","hgnc_symbol"))
# PRONTO.str.ensemble <- 
#   getBM(mart = ensembl,filters = "entrezgene_id",values = PRONTO.str$ID,
#         attributes = c("ensembl_gene_id","entrezgene_id","hgnc_symbol"))
# duplicated(PRONTO.epi.ensemble$entrezgene_id)
# duplicated(PRONTO.str.ensemble$entrezgene_id)

# There are some manual modifications in Ensembl ID list. For the convenience of
# repeating the codes, save them into an xlsx file
# write.xlsx(list(PRONTO.epi.ensemble,PRONTO.str.ensemble),
#            sheetName = c("Epi","Str"),file = "PRONTO ES genes Ensemble ID.xlsx",
#            overwrite = T)
PRONTO.epi.ensemble <- read.xlsx("PRONTO ES genes Ensemble ID.xlsx",
                                 sheet = "Epi")
PRONTO.str.ensemble <- read.xlsx("PRONTO ES genes Ensemble ID.xlsx",
                                 sheet = "Str")
PRONTO.ensemble <- rbind(PRONTO.epi.ensemble,PRONTO.str.ensemble)

## HPA and GTEx consensus data

HPA.Consensus <- import(file = "rna_tissue_consensus.tsv.zip")

# Extract only PRONTO genes
HPA.Consensus <- HPA.Consensus[
  HPA.Consensus$Gene %in% PRONTO.ensemble$ensembl_gene_id,]

# Rearrange the HPA Consensus data.
# The original expression data frame downloaded from HPA website has 5 columns:
# ensemble gene ID, gene symbol, sample ID, the organ/tissue type where tissues 
# are from, nTPM value.
# After rearrangement, an expression data frame, whose rows are samples and columns
# are 1. gene expression FPKM values, 2. the organ where tissues are from, is 
# generated.
colnames(HPA.Consensus) <- c("Gene","Gene_name","Tissue","nTPM")
rearranged.HPA.Consensus <- dcast(HPA.Consensus,Gene_name+Gene ~ Tissue)
rownames(rearranged.HPA.Consensus) <- rearranged.HPA.Consensus$Gene_name
rearranged.HPA.Consensus$Gene_name <- NULL

# Extract Only urinary tract organs from the original data
rearranged.HPA.Consensus <- 
  rearranged.HPA.Consensus[,c("prostate","kidney","urinary bladder")]

## TCGA data

HPA.TCGA <- import(file = "rna_cancer_sample.tsv.zip")

# Extract Only urinary tract organs from the original data
# unique(HPA.TCGA$Cancer)
# The following cancer types are of interest:
# PRAD: prostate adenocarcinoma
# BLCA: bladder urothelial carcinoma
# KICH: Kidney chromophobe cell carcinoma
# KIRC: kidney renal clear cell carcinoma
# KIRP: Kidney renal papillary cell carcinoma
HPA.TCGA <- 
  HPA.TCGA[HPA.TCGA$Cancer %in% c("PRAD","BLCA","KICH","KIRC","KIRP"),]

# Extract only PRONTO genes
HPA.TCGA <- HPA.TCGA[HPA.TCGA$Gene %in% PRONTO.ensemble$ensembl_gene_id,]

# Rearrange HPA TCGA data.
# The original expression data frame downloaded from HPA website has 4 columns:
# ensemble gene ID, sample ID, the organ/tissue type where tissues are from, FPKM 
# value.
# After rearrangement, an expression data frame, whose rows are samples and columns
# are 1. gene expression FPKM values, 2. the organ where tissues are from, is 
# generated.
rearranged.HPA.TCGA <- dcast(HPA.TCGA,Sample ~ Gene)
rearranged.HPA.TCGA$cancer.type <- rearranged.HPA.TCGA$Cancer
rownames(rearranged.HPA.TCGA) <- rearranged.HPA.TCGA$Sample
rearranged.HPA.TCGA$Sample <- NULL

# Add gene symbols and cancer types.
colnames(rearranged.HPA.TCGA) <- 
  PRONTO.ensemble$hgnc_symbol[match(colnames(rearranged.HPA.TCGA),
                                    PRONTO.ensemble$ensembl_gene_id)]
rearranged.HPA.TCGA$cancer.type <- 
  HPA.TCGA$Cancer[match(rownames(rearranged.HPA.TCGA),HPA.TCGA$Sample)]


# Gene expression and fold changes in different organs --------------------

# Summarize gene expression in prostate compared to that in other urinary tract 
# organs in TCGA and HPA Consensus data

## HPA Consensus data

# Log-transform the data
rearranged.HPA.Consensus <- log2(rearranged.HPA.Consensus+1)

# Separate stromal and epithelial genes and calculate minimum fold change
rearranged.HPA.Consensus.epi <- rearranged.HPA.Consensus[
  rownames(rearranged.HPA.Consensus) %in% PRONTO.epi.ensemble$hgnc_symbol,]
rearranged.HPA.Consensus.str <- rearranged.HPA.Consensus[
  rownames(rearranged.HPA.Consensus) %in% PRONTO.str.ensemble$hgnc_symbol,]
# 2 Epi genes are left.
# rownames(rearranged.HPA.Consensus)[
#   !rownames(rearranged.HPA.Consensus) %in% 
#     c(rownames(rearranged.HPA.Consensus.epi),rownames(rearranged.HPA.Consensus.str))]
# Add C16orf70 and FAM122C to Epi ensembl ID list.
# C16orf70 is an alias of PHAF1. FAM122C is an alias of PABIR3.
# PRONTO.epi.ensemble <- edit(PRONTO.epi.ensemble)
# rearranged.HPA.Consensus.epi <- rearranged.HPA.Consensus[
#   rownames(rearranged.HPA.Consensus) %in% PRONTO.epi.ensemble$hgnc_symbol,]
HPA.Consensus.epi.minFC <- apply(rearranged.HPA.Consensus.epi,MARGIN = 1,FUN = function(x){
  x["prostate"]-max(x[c("kidney","urinary bladder")])}) %>% sort(decreasing = T)
HPA.Consensus.epi.minFC <- 2^HPA.Consensus.epi.minFC
HPA.Consensus.str.minFC <- apply(rearranged.HPA.Consensus.str,MARGIN = 1,FUN = function(x){
  x["prostate"]-max(x[c("kidney","urinary bladder")])}) %>% sort(decreasing = T)
HPA.Consensus.str.minFC <- 2^HPA.Consensus.str.minFC

## TCGA data. 

# Log-transform the data and calculate mean difference.
rearranged.HPA.TCGA[,colnames(rearranged.HPA.TCGA)!="cancer.type"] <- 
  log2(rearranged.HPA.TCGA[,colnames(rearranged.HPA.TCGA)!="cancer.type"]+1)
HPA.TCGA.mean <- 
  data.frame(summarise_all(group_by(rearranged.HPA.TCGA,cancer.type),mean))
rownames(HPA.TCGA.mean) <- HPA.TCGA.mean$cancer.type
HPA.TCGA.mean$cancer.type <- NULL
HPA.TCGA.mean <- data.frame(t(HPA.TCGA.mean))

# Separate stromal and epithelial genes and calculate minimum fold change
HPA.TCGA.epi.mean <- HPA.TCGA.mean[
  rownames(HPA.TCGA.mean) %in% PRONTO.epi.ensemble$hgnc_symbol,]
HPA.TCGA.str.mean <- HPA.TCGA.mean[
  rownames(HPA.TCGA.mean) %in% PRONTO.str.ensemble$hgnc_symbol,]
# 1 Epi gene was left
# rownames(HPA.TCGA.mean)[!rownames(HPA.TCGA.mean) %in% 
#                           c(rownames(HPA.TCGA.epi.mean),rownames(HPA.TCGA.str.mean))]
# Add NKX3.1 to Epi ensembl ID list
# PRONTO.epi.ensemble <- edit(PRONTO.epi.ensemble)
# HPA.TCGA.epi.mean <- HPA.TCGA.mean[
#   rownames(HPA.TCGA.mean) %in% PRONTO.epi.ensemble$hgnc_symbol,]
HPA.TCGA.epi.minFC <- apply(HPA.TCGA.epi.mean,MARGIN = 1,FUN = function(x){
  x["PRAD"]-max(x[c("BLCA","KICH","KIRC","KIRP")])}) %>% sort(decreasing = T)
HPA.TCGA.epi.minFC <- 2^HPA.TCGA.epi.minFC
HPA.TCGA.str.minFC <- apply(HPA.TCGA.str.mean,MARGIN = 1,FUN = function(x){
  x["PRAD"]-max(x[c("BLCA","KICH","KIRC","KIRP")])}) %>% sort(decreasing = T)
HPA.TCGA.str.minFC <- 2^HPA.TCGA.str.minFC

# Visualize an example
rearranged.HPA.Consensus.long <- 
  melt(cbind(rearranged.HPA.Consensus,gene=rownames(rearranged.HPA.Consensus)),
       id.vars = "gene",variable.name = "organ",
       value.name = "median.expression")
ggbarplot(
  rearranged.HPA.Consensus.long[rearranged.HPA.Consensus.long$gene=="AZGP1",],
  x = "organ",y = "median.expression",xlab = F,ylab = F,
  title = "median expression",fill = "tomato3")


# Integrate the gene expression and protein expression results ------------

# Str genes
selected.str.RNA <- PRONTO.str.ensemble$entrezgene_id[
  match(names(HPA.Consensus.str.minFC)[HPA.Consensus.str.minFC>1],
        PRONTO.str.ensemble$hgnc_symbol)]
selected.str.protein <- 
  PRONTO.str$ID[PRONTO.str$Normal.Prostate.Specific.protein=="Y" | 
                  PRONTO.str$Prostate.cancer.specific.protein.expression=="Y"]
selected.str <- union(selected.str.protein,selected.str.RNA)
names(selected.str) <- PRONTO.str$Gene[match(selected.str,PRONTO.str$ID)]

# Epi genes
selected.epi.RNA <- PRONTO.epi.ensemble$entrezgene_id[
  match(names(HPA.Consensus.epi.minFC)[HPA.Consensus.epi.minFC>1],
        PRONTO.epi.ensemble$hgnc_symbol)]
selected.epi.protein <- 
  PRONTO.epi$ID[PRONTO.epi$Normal.Prostate.Specific.protein!="N" & 
                  PRONTO.epi$Prostate.cancer.specific.protein.expression!="N"]
selected.epi <- union(selected.epi.protein,selected.epi.RNA)
names(selected.epi) <- PRONTO.epi$Gene[match(selected.epi,PRONTO.epi$ID)]


# Write the results for publishing. ---------------------------------------

# Check if the symbols match first
PRONTO.epi$Gene %in% names(HPA.Consensus.epi.minFC)
PRONTO.epi$Gene %in% names(HPA.TCGA.epi.minFC)
PRONTO.str$Gene %in% names(HPA.Consensus.str.minFC)
PRONTO.str$Gene %in% names(HPA.TCGA.str.minFC)
# Some epithelial genes have different aliases:
# "C16orf70" "FAM122C"  "NKX3-1"   "SQRDL"

# Change their symbols
PRONTO.epi$Gene[!PRONTO.epi$Gene %in% rownames(HPA.TCGA.epi.mean)] <-
  PRONTO.epi.ensemble$hgnc_symbol[
    match(PRONTO.epi$ID[!PRONTO.epi$Gene %in% rownames(HPA.TCGA.epi.mean)],
          PRONTO.epi.ensemble$entrezgene_id)]

# Change "NKX3.1" symbol in HPA.TCGA.epi.minFC into "NKX3-1"
names(HPA.TCGA.epi.minFC)[names(HPA.TCGA.epi.minFC)=="NKX3.1"] <- "NKX3-1"

# Change symbol "C16orf70" "FAM122C" in HPA.Consensus.epi.minFC into "PHAF1"
# "PABIR3"
names(HPA.Consensus.epi.minFC)[
  names(HPA.Consensus.epi.minFC) %in% c("C16orf70","FAM122C")] <- 
  c("PHAF1","PABIR3")

# Change row names of PRONTO.epi and PRONTO.str into gene symbols
rownames(PRONTO.epi) <- PRONTO.epi$Gene
rownames(PRONTO.str) <- PRONTO.str$Gene

# Change detailed description of protein staining (i.e., "a little higher...")
# into "Y"
PRONTO.epi$Normal.Prostate.Specific.protein[
  str_detect(string = PRONTO.epi$Normal.Prostate.Specific.protein,
             pattern = "A litle Higher than\nkidney and bladder")] <- "Y"
PRONTO.epi$Prostate.cancer.specific.protein.expression[
  str_detect(string = PRONTO.epi$Prostate.cancer.specific.protein.expression,
             pattern = "A litle Higher than\nrenal and bladder cancer")] <- "Y"

# Order the gene symbols
HPA.Consensus.epi.minFC <- data.frame(
  HPA.Consensus.epi.minFC,Genes=names(HPA.Consensus.epi.minFC)) %>%
  arrange(Genes,.by_group = T)
HPA.TCGA.epi.minFC <- data.frame(
  HPA.TCGA.epi.minFC,Genes=names(HPA.TCGA.epi.minFC)) %>%
  arrange(Genes,.by_group = T)
HPA.Consensus.str.minFC <- data.frame(
  HPA.Consensus.str.minFC,Genes=names(HPA.Consensus.str.minFC)) %>%
  arrange(Genes,.by_group = T)
HPA.TCGA.str.minFC <- data.frame(
  HPA.TCGA.str.minFC,Genes=names(HPA.TCGA.str.minFC)) %>%
  arrange(Genes,.by_group = T)
PRONTO.epi <- arrange(PRONTO.epi,Gene,.by_group = T)
PRONTO.str <- arrange(PRONTO.str,Gene,.by_group = T)

PRONTO.epi.specificity <- cbind(
  HPA.Consensus.epi.minFC,HPA.TCGA.epi.minFC,PRONTO.epi)
PRONTO.str.specificity <- cbind(
  HPA.Consensus.str.minFC,HPA.TCGA.str.minFC,PRONTO.str)
PRONTO.epi.specificity$Genes <- NULL
PRONTO.epi.specificity$Genes <- NULL
PRONTO.epi.specificity$Gene <- NULL
PRONTO.epi.specificity$ID <- NULL
PRONTO.str.specificity$Genes <- NULL
PRONTO.str.specificity$Genes <- NULL
PRONTO.str.specificity$Gene <- NULL
PRONTO.str.specificity$ID <- NULL

write.xlsx(list(PRONTO.epi.specificity,PRONTO.str.specificity),
           "results/selected genes prostate specificity.xlsx",overwrite = T,
           colNames=T,rowNames=T,
           sheetName = c("epi.specificity","str.specificity"))
