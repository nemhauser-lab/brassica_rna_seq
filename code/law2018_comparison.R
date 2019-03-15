# comparison to Law et. al. 2018
# NOTE: this script references lists from analysis_v3, run that file first.

library(clipr)

temp <- clipr::read_clip()
temp2 <- stringr::str_split(temp, "\t")
temp3 <- data.frame(do.call(rbind, temp2 ), stringsAsFactors=FALSE)
cols <- temp3[1, ]
cols <- gsub(" - ", "_", cols)
cols <- gsub(" ", "_", cols)
colnames(temp3) <- cols
temp3 <- temp3[-1, ]

# convert character columns to numeric
temp3$T0 <- as.numeric(temp3$T0)
temp3$DP_D3 <- as.numeric(temp3$DP_D3)
temp3$DP_D6 <- as.numeric(temp3$DP_D6)
temp3$IDL_D3 <- as.numeric(temp3$IDL_D3)
temp3$IDL_D6 <- as.numeric(temp3$IDL_D6)

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
# temp5$idl_d3_vs_d0 <- as.numeric(temp5$idl_d3_vs_d0)
# temp5$idl_d6_vs_d0 <- as.numeric(temp5$idl_d6_vs_d0)

# temp5$DP_D3 <- as.numeric(temp5$DP_D3)
# temp5$DP_D6 <- as.numeric(temp5$DP_D6)
#temp5$T0 <- as.numeric(temp5$T0)
# temp5$dp_d3_vs_d0 <- log(temp5$DP_D3/temp5$T0, 2)
# temp5$dp_d6_vs_d0 <- log(temp5$DP_D6/temp5$T0, 2)


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






