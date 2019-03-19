# comparison to Law et. al. 2018
# NOTE: this script references lists from analysis_v3, run that file first.

library(clipr)
library(openxlsx)
lawS1Path <- file.path(rootDir, "data", "Law_et_al_2018",
                       "PP2018-RA-00062DR1_Supplemental_Table_1.xlsx")

df <- read.xlsx(xlsxFile=lawS1Path, sheet=1, skipEmptyRows=FALSE, startRow=2)
col <- colnames(df)
col <- gsub(".-.", "_", col)
col <- gsub("\\.", "_", col)
colnames(df) <- col
temp3 <- df

# # take data from clipboard and parse
# temp <- clipr::read_clip()
# temp2 <- stringr::str_split(temp, "\t")
# temp3 <- data.frame(do.call(rbind, temp2 ), stringsAsFactors=FALSE)
# cols <- temp3[1, ]
# cols <- gsub(" - ", "_", cols)
# cols <- gsub(" ", "_", cols)
# colnames(temp3) <- cols
# temp3 <- temp3[-1, ]
# # convert character columns to numeric
# temp3$T0 <- as.numeric(temp3$T0)
# temp3$DP_D3 <- as.numeric(temp3$DP_D3)
# temp3$DP_D6 <- as.numeric(temp3$DP_D6)
# temp3$IDL_D3 <- as.numeric(temp3$IDL_D3)
# temp3$IDL_D6 <- as.numeric(temp3$IDL_D6)

# calculate fold changes
temp3$idl_d3_vs_d0 <- log(temp3$IDL_D3/temp3$T0, 2)
temp3$idl_d6_vs_d0 <- log(temp3$IDL_D6/temp3$T0, 2)
temp3$dp_d3_vs_d0 <- log(temp3$DP_D3/temp3$T0, 2)
temp3$dp_d6_vs_d0 <- log(temp3$DP_D6/temp3$T0, 2)


temp4 <- allGenesWT_SvP
temp4$Gene_ID <- row.names(temp4)
temp4 <- dplyr::left_join(temp4, edinburghGeneInfo, by="Gene_ID")

commonGenes <- temp3[temp3$Locus_ID %in% temp4$A._thaliana_best_hit, "Locus_ID"]

temp5 <- dplyr::left_join(temp3, temp4[,c(1:8)], by=c("Locus_ID" = "A._thaliana_best_hit"))

sum(is.na(temp5$idl_d3_vs_d0))
temp5 <- temp5[!is.na(temp5$idl_d3_vs_d0), ]
temp5 <- temp5[!is.na(temp5$logFC), ]

# standard deviations of idl and dp columns
sdD3Idl <- sd(temp5$idl_d3_vs_d0)
sdD6Idl <- sd(temp5$idl_d6_vs_d0)
sdD3Dp <- sd(temp5$dp_d3_vs_d0)
sdD6Dp <- sd(temp5$dp_d6_vs_d0)

# Analysis below -------------------------------------------

ggplot(temp5, aes(x=idl_d3_vs_d0, y=logFC)) + geom_point(stroke=0, alpha=0.5)
ggplot(temp5, aes(x=idl_d6_vs_d0, y=logFC)) + geom_point(stroke=0, alpha=0.5)
ggplot(temp5, aes(x=idl_d3_vs_d0, y=idl_d6_vs_d0)) + geom_point(stroke=0, alpha=0.5)

cor(temp5[, c("idl_d3_vs_d0",
              "idl_d6_vs_d0",
              "dp_d3_vs_d0",
              "dp_d6_vs_d0",
              "logFC")], method="pearson")

cor(temp5[, c("idl_d3_vs_d0",
              "idl_d6_vs_d0",
              "dp_d3_vs_d0",
              "dp_d6_vs_d0",
              "logFC")], method="spearman")


ggplot(temp5, aes(x=dp_d3_vs_d0, y=logFC)) + geom_point(stroke=0, alpha=0.5)
ggplot(temp5, aes(x=dp_d6_vs_d0, y=logFC)) + geom_point(stroke=0, alpha=0.5)
ggplot(temp5, aes(x=dp_d3_vs_d0, y=dp_d6_vs_d0)) + geom_point(stroke=0, alpha=0.5)

ggplot(temp5, aes(x=idl_d6_vs_d0, y=dp_d6_vs_d0)) + geom_point(stroke=0, alpha=0.5)

ggplot(temp3, aes(x=idl_d3_vs_d0, y=idl_d6_vs_d0)) + geom_point()
plot(temp5$idl_d3_vs_d0, temp5$logFC)

nrow(subset(temp4, FDR<0.01 & abs(logFC)>1 & (A._thaliana_best_hit %in% commonGenes)
            , select=A._thaliana_best_hit))


sum(temp3$Locus_ID[1:1000] %in% subset(temp4, FDR<0.01&abs(logFC)>1)$A._thaliana_best_hit)

nrow(subset(temp5, FDR<0.01 & abs(logFC)>1))
nrow(subset(temp5, abs(idl_d3_vs_d0)>2))
nrow(subset(temp5, FDR<0.01 & abs(logFC)>1 & abs(idl_d3_vs_d0)>2))

sdD3Idl <- sd(temp5$idl_d3_vs_d0)
sdD6Idl <- sd(temp5$idl_d6_vs_d0)
sdD3Dp <- sd(temp5$dp_d3_vs_d0)
sdD6Dp <- sd(temp5$dp_d6_vs_d0)

# D3 vs D0 IDL only
nrow(subset(temp5, abs(idl_d3_vs_d0)>2 & abs(dp_d3_vs_d0)<2))
nrow(subset(temp5, abs(idl_d3_vs_d0)>2 & abs(dp_d3_vs_d0)<2 & FDR<0.01 & abs(logFC)>1))

# D3 vs D0 DP only
nrow(subset(temp5, abs(dp_d3_vs_d0)>2 & abs(idl_d3_vs_d0)<2))
nrow(subset(temp5, abs(dp_d3_vs_d0)>2 & abs(idl_d3_vs_d0)<2 & FDR<0.01 & abs(logFC)>1))

# D6 vs D0 IDL only
nrow(subset(temp5, abs(idl_d6_vs_d0)>2 & abs(dp_d6_vs_d0)<2))
nrow(subset(temp5, abs(idl_d6_vs_d0)>2 & abs(dp_d6_vs_d0)<2 & FDR<0.01 & abs(logFC)>1))

# D6 vs D0 DP only
nrow(subset(temp5, abs(dp_d6_vs_d0)>2 & abs(idl_d6_vs_d0)<2))
nrow(subset(temp5, abs(dp_d6_vs_d0)>2 & abs(idl_d6_vs_d0)<2 & FDR<0.01 & abs(logFC)>1))


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

temp6 <- subset(temp5, Locus_ID %in% df$AGI)

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



sum(temp6$dp_d6_vs_d0 * temp6$litSign > 0)
sum(temp6$idl_d6_vs_d0 * temp6$litSign > 0)



