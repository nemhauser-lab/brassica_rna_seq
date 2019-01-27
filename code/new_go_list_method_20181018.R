temp <- c(
  rownames(DEGsWT_RvS)[!(rownames(DEGsWT_RvS) %in% genesWithHomologs)],
  rownames(DEGsphyB_RvS)[!(rownames(DEGsphyB_RvS) %in% genesWithHomologs)],
  rownames(DEGsWT_SvP)[!(rownames(DEGsWT_SvP) %in% genesWithHomologs)],
  rownames(DEGsphyB_SvP)[!(rownames(DEGsphyB_SvP) %in% genesWithHomologs)]
)

# 17 genes missing homologs in RvS and SvP sets
unique(temp)

###############################################################################


setGroup <- function(row){
  group <- "none"
  if (row$FDR.x < 0.05 & abs(row$logFC.x) > 1.5) {
    group <- "WT_sig"
      if (row$logFC.x - row$logFC.y > 1.25) {
        group <- "WT_up"
      }
      if (row$logFC.x - row$logFC.y < -1.25) {
        group <- "WT_down"
      }
  }

  return(group)
}

group <- plyr::aaply(.data=allGenes_RvS, .margins=1, .fun=setGroup, .expand=FALSE)
allGenes_RvS$group <- group


# allGenes_RvS$group <- rep("Neither", nrow(allGenes_RvS))
# allGenes_RvS[row.names(RvS_WT_Only), ]$group <- "WT_Only"
# allGenes_RvS[row.names(RvS_phyB_Only), ]$group <- "phyB_Only"
# allGenes_RvS[RvS_CommonGenes$Gene_ID, ]$group <- "Both"
# allGenes_RvS[row.names(DEGs_RvS), ]$group <- "WT(R-S) vs. phyB(R-S) FDR<0.1"

cols <- c("none"="grey", "WT_sig"="blue", "WT_up"="red", "WT_down"="green")

ggplot(allGenes_RvS, aes(x=logFC.x, y=logFC.y, color=group)) +
  geom_point(data=subset(allGenes_RvS, group=="none"), alpha=0.2) +
  geom_point(data=subset(allGenes_RvS, group=="WT_sig"), alpha=0.5) +
  geom_point(data=subset(allGenes_RvS, group=="WT_up")) +
  geom_point(data=subset(allGenes_RvS, group=="WT_down")) +
  scale_color_manual(values=cols) +
  labs(x = "WT logFC", y="phyB logFC", title="logFC R vs. S")

ggplot(allGenes_RvS, aes(x=PValue.x, y=PValue.y, color=group)) +
  geom_point(data=subset(allGenes_RvS, group=="none"), alpha=0.2) +
  geom_point(data=subset(allGenes_RvS, group=="WT_sig"), alpha=0.5) +
  geom_point(data=subset(allGenes_RvS, group=="WT_up")) +
  geom_point(data=subset(allGenes_RvS, group=="WT_down")) +
  scale_color_manual(values=cols) +
  scale_x_log10() + scale_y_log10() +
  labs(x = "WT p-value", y="phyB p-value", title="p-values R vs. S")


# -----------------------------------------------------------------------------
# WT up more
WT_up <- rownames(allGenes_RvS[allGenes_RvS$FDR.x < 0.05 & allGenes_RvS$logFC.x > 1.5, ])
WT_up_more <- rownames(allGenes_RvS[(rownames(allGenes_RvS) %in% WT_up) & (allGenes_RvS$group == "WT_up"), ])

WT_up_At <- subset(allGenesEdin, Gene_ID %in% WT_up)$A._thaliana_best_hit
WT_up_more_At <- subset(allGenesEdin, Gene_ID %in% WT_up_more)$A._thaliana_best_hit

fname <- file.path(dirPath, "WT_up.csv")
write.table(WT_up_At, fname, quote=FALSE, row.names=FALSE, na="", col.names=FALSE, sep=",")

fname <- file.path(dirPath, "WT_up_more.csv")
write.table(WT_up_more_At, fname, quote=FALSE, row.names=FALSE, na="", col.names=FALSE, sep=",")

# -----------------------------------------------------------------------------
# WT down more
WT_down <- rownames(allGenes_RvS[allGenes_RvS$FDR.x < 0.05 & allGenes_RvS$logFC.x < -1.5, ])
WT_down_more <- rownames(allGenes_RvS[(rownames(allGenes_RvS) %in% WT_down) & (allGenes_RvS$group == "WT_down"), ])

WT_down_At <- subset(allGenesEdin, Gene_ID %in% WT_down)$A._thaliana_best_hit
WT_down_more_At <- subset(allGenesEdin, Gene_ID %in% WT_down_more)$A._thaliana_best_hit

fname <- file.path(dirPath, "WT_down.csv")
write.table(WT_down_At, fname, quote=FALSE, row.names=FALSE, na="", col.names=FALSE, sep=",")

fname <- file.path(dirPath, "WT_down_more.csv")
write.table(WT_down_more_At, fname, quote=FALSE, row.names=FALSE, na="", col.names=FALSE, sep=",")
# -----------------------------------------------------------------------------

################################################################################


setGroup <- function(row){
  group <- "none"
  if (row$FDR.y < 0.05 & abs(row$logFC.y) > 1.5) {
    group <- "phyB_sig"
    if (row$logFC.y - row$logFC.x > 1.25) {
      group <- "phyB_up"
    }
    if (row$logFC.y - row$logFC.x < -1.25) {
      group <- "phyB_down"
    }
  }

  return(group)
}

group <- plyr::aaply(.data=allGenes_RvS, .margins=1, .fun=setGroup, .expand=FALSE)
allGenes_RvS$group <- group


cols <- c("none"="grey", "phyB_sig"="blue", "phyB_up"="red", "phyB_down"="green")

ggplot(allGenes_RvS, aes(x=logFC.x, y=logFC.y, color=group)) +
  geom_point(data=subset(allGenes_RvS, group=="none"), alpha=0.2) +
  geom_point(data=subset(allGenes_RvS, group=="phyB_sig"), alpha=0.5) +
  geom_point(data=subset(allGenes_RvS, group=="phyB_up")) +
  geom_point(data=subset(allGenes_RvS, group=="phyB_down")) +
  scale_color_manual(values=cols) +
  labs(x = "WT logFC", y="phyB logFC", title="logFC R vs. S")

ggplot(allGenes_RvS, aes(x=PValue.x, y=PValue.y, color=group)) +
  geom_point(data=subset(allGenes_RvS, group=="none"), alpha=0.2) +
  geom_point(data=subset(allGenes_RvS, group=="phyB_sig"), alpha=0.5) +
  geom_point(data=subset(allGenes_RvS, group=="phyB_up")) +
  geom_point(data=subset(allGenes_RvS, group=="phyB_down")) +
  scale_color_manual(values=cols) +
  scale_x_log10() + scale_y_log10() +
  labs(x = "WT p-value", y="phyB p-value", title="p-values R vs. S")

# -----------------------------------------------------------------------------
# phyB up more
phyB_up <- rownames(allGenes_RvS[allGenes_RvS$FDR.y < 0.05 & allGenes_RvS$logFC.y > 1.5, ])
phyB_up_more <- rownames(allGenes_RvS[(rownames(allGenes_RvS) %in% phyB_up) & (allGenes_RvS$group == "phyB_up"), ])

phyB_up_At <- subset(allGenesEdin, Gene_ID %in% phyB_up)$A._thaliana_best_hit
phyB_up_more_At <- subset(allGenesEdin, Gene_ID %in% phyB_up_more)$A._thaliana_best_hit

fname <- file.path(dirPath, "phyB_up.csv")
write.table(phyB_up_At, fname, quote=FALSE, row.names=FALSE, na="", col.names=FALSE, sep=",")

fname <- file.path(dirPath, "phyB_up_more.csv")
write.table(phyB_up_more_At, fname, quote=FALSE, row.names=FALSE, na="", col.names=FALSE, sep=",")

# -----------------------------------------------------------------------------
# WT down more
phyB_down <- rownames(allGenes_RvS[allGenes_RvS$FDR.y < 0.05 & allGenes_RvS$logFC.y < -1.5, ])
phyB_down_more <- rownames(allGenes_RvS[(rownames(allGenes_RvS) %in% phyB_down) & (allGenes_RvS$group == "phyB_down"), ])

phyB_down_At <- subset(allGenesEdin, Gene_ID %in% phyB_down)$A._thaliana_best_hit
phyB_down_more_At <- subset(allGenesEdin, Gene_ID %in% phyB_down_more)$A._thaliana_best_hit

fname <- file.path(dirPath, "phyB_down.csv")
write.table(phyB_down_At, fname, quote=FALSE, row.names=FALSE, na="", col.names=FALSE, sep=",")

fname <- file.path(dirPath, "phyB_down_more.csv")
write.table(phyB_down_more_At, fname, quote=FALSE, row.names=FALSE, na="", col.names=FALSE, sep=",")
# -----------------------------------------------------------------------------

################################################################################


setGroup <- function(row){
  group <- "none"
  if ((row$FDR.x < 0.01 & row$logFC.x > 1) | (row$FDR.y < 0.01 & row$logFC.y > 1)) {
    group <- "sig"
    if (row$FDR.y > 0.01) {
      group <- "WT_sig"
    }else if (row$FDR.x > 0.01) {
      group <- "phyB_sig"
    } else {
      group <- "Both_sig"
    }

  }

  return(group)
}

group <- plyr::aaply(.data=allGenes_RvS, .margins=1, .fun=setGroup, .expand=FALSE)
allGenes_RvS$group <- group

sum(allGenes_RvS$group == "phyB_sig")


cols <- c("none"="grey50", "phyB_sig"="red", "WT_sig"="blue", "Both_sig"="purple")

ggplot(allGenes_RvS, aes(x=logFC.x, y=logFC.y, color=group)) +
  geom_point(data=subset(allGenes_RvS, group=="none"), alpha=0.2, stroke=0) +
  geom_point(data=subset(allGenes_RvS, group=="Both_sig"), alpha=0.5, stroke=0) +
  geom_point(data=subset(allGenes_RvS, group=="phyB_sig"), alpha=0.5, stroke=0) +
  geom_point(data=subset(allGenes_RvS, group=="WT_sig"), alpha=0.5, stroke=0) +
  scale_color_manual(values=cols) +
  labs(x = "WT logFC", y="phyB logFC", title="logFC R vs. S")


# Up
phyB_only_up <- rownames(subset(allGenes_RvS, group=="phyB_sig"))
WT_only_up <- rownames(subset(allGenes_RvS, group=="WT_sig"))
both_sig_up <- rownames(subset(allGenes_RvS, group=="Both_sig"))

phyB_up_At <- subset(allGenesEdin, Gene_ID %in% phyB_only_up)$A._thaliana_best_hit
WT_up_At <- subset(allGenesEdin, Gene_ID %in% WT_only_up)$A._thaliana_best_hit
All_up_At <- subset(allGenesEdin, Gene_ID %in% c(WT_only_up, phyB_only_up, both_sig_up))$A._thaliana_best_hit

fname <- file.path(dirPath, "phyB_only_sig_up.csv")
write.table(phyB_up_At, fname, quote=FALSE, row.names=FALSE, na="", col.names=FALSE, sep=",")

fname <- file.path(dirPath, "WT_only_sig_up.csv")
write.table(WT_up_At, fname, quote=FALSE, row.names=FALSE, na="", col.names=FALSE, sep=",")

fname <- file.path(dirPath, "either_sig_up.csv")
write.table(All_up_At, fname, quote=FALSE, row.names=FALSE, na="", col.names=FALSE, sep=",")


#-------------------------------------------------------------------------------
setGroup <- function(row){
  group <- "none"
  if ((row$FDR.x < 0.01 & row$logFC.x < -1) | (row$FDR.y < 0.01 & row$logFC.y < -1)) {
    group <- "sig"
    if (row$FDR.y > 0.01) {
      group <- "WT_sig"
    }else if (row$FDR.x > 0.01) {
      group <- "phyB_sig"
    } else {
      group <- "Both_sig"
    }
  }
  return(group)
}

group <- plyr::aaply(.data=allGenes_RvS, .margins=1, .fun=setGroup, .expand=FALSE)
allGenes_RvS$group <- group


# Down
phyB_only_down <- rownames(subset(allGenes_RvS, group=="phyB_sig"))
WT_only_down <- rownames(subset(allGenes_RvS, group=="WT_sig"))
both_sig_down <- rownames(subset(allGenes_RvS, group=="Both_sig"))

phyB_down_At <- subset(allGenesEdin, Gene_ID %in% phyB_only_down)$A._thaliana_best_hit
WT_down_At <- subset(allGenesEdin, Gene_ID %in% WT_only_down)$A._thaliana_best_hit
All_down_At <- subset(allGenesEdin, Gene_ID %in% c(WT_only_down, phyB_only_down, both_sig_down))$A._thaliana_best_hit

fname <- file.path(dirPath, "phyB_only_sig_down.csv")
write.table(phyB_down_At, fname, quote=FALSE, row.names=FALSE, na="", col.names=FALSE, sep=",")

fname <- file.path(dirPath, "WT_only_sig_down.csv")
write.table(WT_down_At, fname, quote=FALSE, row.names=FALSE, na="", col.names=FALSE, sep=",")

fname <- file.path(dirPath, "either_sig_down.csv")
write.table(All_down_At, fname, quote=FALSE, row.names=FALSE, na="", col.names=FALSE, sep=",")


################################################################################


temp <- geneHeatmap(c(phyB_only_down),
                    DGEdata[, !grepl("x|(24)", colnames(DGEdata))],
                    rsOnly=FALSE)

# _At_homologs_FDR01_LFC_01_20190125
