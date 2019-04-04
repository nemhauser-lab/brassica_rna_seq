# comparison to Law et. al. 2018
# NOTE: this script references lists from analysis_v3, run that file first.

library(rprojroot)
library(openxlsx)

rootDir <- find_root(is_rstudio_project)
lawS1Path <- file.path(rootDir, "data", "Law_et_al_2018",
                       "PP2018-RA-00062DR1_Supplemental_Table_1.xlsx")

# getting the data -------------------------------------------------------------

df <- read.xlsx(xlsxFile=lawS1Path, sheet=1, skipEmptyRows=FALSE, startRow=2)
col <- colnames(df)
col <- gsub("\\.-\\.", "_", col)
col <- gsub("\\.", "_", col)
colnames(df) <- col
lawDF <- df

# calculate fold changes
lawDF$IDL_D3_vs_D0 <- log(lawDF$IDL_D3/lawDF$T0, 2)
lawDF$IDL_D6_vs_D0 <- log(lawDF$IDL_D6/lawDF$T0, 2)
lawDF$DP_D3_vs_D0 <- log(lawDF$DP_D3/lawDF$T0, 2)
lawDF$DP_D6_vs_D0 <- log(lawDF$DP_D6/lawDF$T0, 2)

# WT S vs P
SvP_DF <- allGenesWT_SvP
SvP_DF$Gene_ID <- row.names(SvP_DF)
SvP_DF <- dplyr::left_join(SvP_DF, edinburghGeneInfo, by="Gene_ID")

commonGenes <- lawDF[lawDF$Locus_ID %in% SvP_DF$A._thaliana_best_hit, "Locus_ID"]

# combined dataframe
CommonDF <- dplyr::left_join(lawDF, SvP_DF[,c(1:8)], by=c("Locus_ID" = "A._thaliana_best_hit"))
# sum(is.na(CommonDF$IDL_D3_vs_D0))
CommonDF <- CommonDF[!is.na(CommonDF$IDL_D3_vs_D0), ]
CommonDF <- CommonDF[!is.na(CommonDF$logFC), ]
CommonDF$logFC_WT <- CommonDF$logFC
CommonDF$logFC_phyB <- allGenesphyB_SvP[CommonDF$Gene_ID, "logFC"]

# standard deviations of idl and dp columns
sdD3Idl <- sd(CommonDF$IDL_D3_vs_D0)
sdD6Idl <- sd(CommonDF$IDL_D6_vs_D0)
sdD3Dp <- sd(CommonDF$DP_D3_vs_D0)
sdD6Dp <- sd(CommonDF$DP_D6_vs_D0)

# Analysis/comparisons ---------------------------------------------------------

ggplot(CommonDF, aes(x=IDL_D3_vs_D0, y=logFC_WT)) + geom_point(stroke=0, alpha=0.5)
ggplot(CommonDF, aes(x=IDL_D6_vs_D0, y=logFC_WT)) + geom_point(stroke=0, alpha=0.5)
ggplot(CommonDF, aes(x=IDL_D3_vs_D0, y=IDL_D6_vs_D0)) + geom_point(stroke=0, alpha=0.5)

# correlation between log fold changes in different sets.

cor(CommonDF[, c("IDL_D3_vs_D0",
              "IDL_D6_vs_D0",
              "DP_D3_vs_D0",
              "DP_D6_vs_D0",
              "logFC_WT")], method="pearson")

cor(CommonDF[, c("IDL_D3_vs_D0",
              "IDL_D6_vs_D0",
              "DP_D3_vs_D0",
              "DP_D6_vs_D0",
              "logFC_WT")], method="spearman")


ggplot(CommonDF, aes(x=DP_D3_vs_D0, y=logFC_WT)) + geom_point(stroke=0, alpha=0.5)
ggplot(CommonDF, aes(x=DP_D6_vs_D0, y=logFC_WT)) + geom_point(stroke=0, alpha=0.5)
ggplot(CommonDF, aes(x=DP_D3_vs_D0, y=DP_D6_vs_D0)) + geom_point(stroke=0, alpha=0.5)

ggplot(CommonDF, aes(x=IDL_D6_vs_D0, y=DP_D6_vs_D0)) + geom_point(stroke=0, alpha=0.5)

ggplot(lawDF, aes(x=IDL_D3_vs_D0, y=IDL_D6_vs_D0)) + geom_point()
plot(CommonDF$IDL_D3_vs_D0, CommonDF$logFC)

nrow(subset(SvP_DF, FDR<0.01 & abs(logFC)>1 & (A._thaliana_best_hit %in% commonGenes)
            , select=A._thaliana_best_hit))


sum(lawDF$Locus_ID[1:1000] %in% subset(SvP_DF, FDR<0.01&abs(logFC)>1)$A._thaliana_best_hit)

nrow(subset(CommonDF, FDR<0.01 & abs(logFC)>1))
nrow(subset(CommonDF, abs(IDL_D3_vs_D0)>2))
nrow(subset(CommonDF, FDR<0.01 & abs(logFC)>1 & abs(IDL_D3_vs_D0)>2))

sdD3Idl <- sd(CommonDF$IDL_D3_vs_D0)
sdD6Idl <- sd(CommonDF$IDL_D6_vs_D0)
sdD3Dp <- sd(CommonDF$DP_D3_vs_D0)
sdD6Dp <- sd(CommonDF$DP_D6_vs_D0)


length(unique(CommonDF$Locus_ID))
length(unique(subset(CommonDF, FDR<2 & abs(logFC)>1)$Locus_ID))


# D3 vs D0 IDL only
nrow(subset(CommonDF, abs(IDL_D3_vs_D0)>2 & abs(DP_D3_vs_D0)<2))
nrow(subset(CommonDF, abs(IDL_D3_vs_D0)>2 & abs(DP_D3_vs_D0)<2 & FDR<0.01 & abs(logFC)>1))

length(unique(subset(CommonDF, abs(IDL_D3_vs_D0)>1 & abs(DP_D3_vs_D0)<1)$Locus_ID))
length(unique(subset(CommonDF, abs(IDL_D3_vs_D0)>1 & abs(DP_D3_vs_D0)<1 & FDR<2 & abs(logFC)>1)$Locus_ID))

# D3 vs D0 DP only
nrow(subset(CommonDF, abs(DP_D3_vs_D0)>2 & abs(IDL_D3_vs_D0)<2))
nrow(subset(CommonDF, abs(DP_D3_vs_D0)>2 & abs(IDL_D3_vs_D0)<2 & FDR<0.01 & abs(logFC)>1))

length(unique(subset(CommonDF, abs(DP_D3_vs_D0)>1 & abs(IDL_D3_vs_D0)<1)$Locus_ID))
length(unique(subset(CommonDF, abs(DP_D3_vs_D0)>1 & abs(IDL_D3_vs_D0)<1 & FDR<2 & abs(logFC)>1)$Locus_ID))

# D6 vs D0 IDL only
nrow(subset(CommonDF, abs(IDL_D6_vs_D0)>2 & abs(DP_D6_vs_D0)<2))
nrow(subset(CommonDF, abs(IDL_D6_vs_D0)>2 & abs(DP_D6_vs_D0)<2 & FDR<0.01 & abs(logFC)>2))

length(unique(subset(CommonDF, abs(IDL_D6_vs_D0)>1 & abs(DP_D6_vs_D0)<1)$Locus_ID))
length(unique(subset(CommonDF, abs(IDL_D6_vs_D0)>1 & abs(DP_D6_vs_D0)<1 & FDR<2 & abs(logFC)>1)$Locus_ID))

# D6 vs D0 DP only
nrow(subset(CommonDF, abs(DP_D6_vs_D0)>2 & abs(IDL_D6_vs_D0)<2))
nrow(subset(CommonDF, abs(DP_D6_vs_D0)>2 & abs(IDL_D6_vs_D0)<2 & FDR<0.01 & abs(logFC)>1))

length(unique(subset(CommonDF, abs(DP_D6_vs_D0)>1 & abs(IDL_D6_vs_D0)<1)$Locus_ID))
length(unique(subset(CommonDF, abs(DP_D6_vs_D0)>1 & abs(IDL_D6_vs_D0)<1 & FDR<2 & abs(logFC)>1)$Locus_ID))


# ------------------------------------------------------------------------------
# Law 2018 supplemental table 2

lawS2Path <- file.path(rootDir, "data", "Law_et_al_2018",
                       "PP2018-RA-00062DR1_Supplemental_Table_2.xlsx")

df <- read.xlsx(xlsxFile=lawS2Path, sheet=1, skipEmptyRows=FALSE, startRow=2,
                cols=1:8)

# Table S2 uses cell color to indicate if genes were up or down regulated in
# the liturature, we are hard coding those values below.
litValSign <- c(1, -1, 1, -1, 1, 1, -1, 1, 1, 1, 1, -1, -1, -1, -1, -1, -1, -1,
                -1, 1, -1, 1, -1, -1, -1, 1, -1, -1, 1, 1, 1, -1, 1, -1, -1, -1)

df$litSign <- litValSign

rownames(df) <- df$AGI

temp6 <- subset(CommonDF, Locus_ID %in% df$AGI)

temp6$litSign <- df[temp6$Locus_ID, "litSign"]
temp6$logFC_phyB <- allGenesphyB_SvP[temp6$Gene_ID, "logFC"]

sum(df$`IDL_D6-D3` * df$litSign > 0)/nrow(df)
# [1] 0.75  ???

# fraction of genes similar between WT S-P and lit values
sum(temp6$litSign * temp6$logFC > 0)/nrow(temp6)
# 0.7755

# fraction of genes similar between phyB S-P and lit values
sum(temp6$litSign * temp6$logFC_phyB > 0)/nrow(temp6)
# 0.7347

# binomial significance test
binom.test(32,49,0.5,alternative="greater")

sum(temp6$DP_D6_vs_D0 * temp6$litSign > 0)
sum(temp6$IDL_D6_vs_D0 * temp6$litSign > 0)



