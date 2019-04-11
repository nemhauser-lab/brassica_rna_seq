# utility functions for RNA-seq analysis

#' generate folder structure for results,
#'
#' @param dirPath location to put the results folder in
#' @return a path to the top level results folder
#'
makeResultsFolder <- function(dirPath){

  date <- format(Sys.time(), "%Y%m%d")
  resultsPath <- file.path(dirPath, paste("analysis_output", date, sep="_"))
  if (file.exists(resultsPath)) {
    i=2
    while (file.exists(paste(resultsPath, "_", i, sep=""))) {
      i <- i + 1
    }
    resultsPath <- paste(resultsPath, "_", i, sep="")
  }
  dir.create(resultsPath)

  folders <- c("geno_x_timepoint",
               "geno_x_timepoint/GO_lists",
               "geno_x_timepoint/plots",
               "geno_x_timepoint/promoters",
               "phyB_v_WT_at_timepoint",
               "phyB_v_WT_at_timepoint/GO_lists",
               "phyB_v_WT_at_timepoint/plots",
               "phyB_v_WT_at_timepoint/promoters",
               "wound_control",
               "wound_control/GO_lists",
               "wound_control/plots",
               "wound_control/promoters",
               "timepoint_v_timepoint",
               "timepoint_v_timepoint/GO_lists",
               "timepoint_v_timepoint/plots",
               "timepoint_v_timepoint/promoters",
               "timepoint_v_timepoint/RvP_venn",
               "timepoint_v_timepoint/RvP_venn/GO_lists",
               "timepoint_v_timepoint/RvP_venn/promoters",
               "timepoint_v_timepoint/RvD_venn",
               "timepoint_v_timepoint/RvD_venn/GO_lists",
               "timepoint_v_timepoint/RvD_venn/promoters",
               "timepoint_v_timepoint/DvP_venn",
               "timepoint_v_timepoint/DvP_venn/GO_lists",
               "timepoint_v_timepoint/DvP_venn/promoters"
  )
  for (folder in folders) {
    dir.create(file.path(resultsPath, folder))
  }
  return(resultsPath)
}



#' generate a filename for an output file
#'
#' @param folder the folder where the file will be saved
#' @param cont the contrast used to generate the output
#' @param change "_UP" or "_DOWN" or default is ""
#' @param dtype description of the data in the file to be saved
#' @param FDRmax the FDR cutoff used to generate the results
#' @param logFCmin the log fold change cutoff used to generate the results
#' @return file name includign folder
#'
makeFname <- function(folder, cont, change="", dtype="", FDRmax, logFCmin,
                      dirPath){
  filename <- paste(folder, "/", cont, change, dtype, "_FDR",
                    str_pad(FDRmax*100, 2, pad="0"), "_LFC",
                    str_pad(logFCmin,2, pad=0), "_",
                    format(Sys.time(), "%Y%m%d"), sep="")
  if (dtype %in% c("", "_At_homologs")){
    filename <- paste(filename, ".csv", sep="")
  } else if (grepl("plot|hist|diag", dtype)){
    filename <- paste(filename, ".png", sep="")
  } else if (dtype == "_promoters") {
    filename <- paste(filename, ".fasta", sep="")
  }
  filename <- file.path(dirPath, filename)
  return(filename)
}

#' generate and saves various files associated with a list of differentially
#' expressed genes
#'
#'
save_DEGs <- function(DEGList, folder, cont, FDRmax, logFCmin, dir){
  DEGs <- DEGList
  DEGs$Gene_ID <- row.names(DEGs)
  DEGs <- dplyr::left_join(DEGs, edinburghGeneInfo, by="Gene_ID")
  row.names(DEGs) <- DEGs$Gene_ID
  upDEGs <- subset(DEGs, logFC > 0)
  downDEGs <- subset(DEGs, logFC < 0)
  # UP DEGs
  fname <- makeFname(folder=folder, cont=cont, change="_UP", FDRmax=FDRmax,
                     logFCmin=logFCmin, dirPath=dir)
  write.csv(upDEGs, fname, row.names=TRUE)
  # DOWN DEGs
  fname <- makeFname(folder=folder, cont=cont, change="_DOWN", FDRmax=FDRmax,
                     logFCmin=logFCmin, dirPath=dir)
  write.csv(downDEGs, fname, row.names=TRUE)
  # UP GO lists
  fname <- makeFname(folder=paste(folder, "/GO_lists", sep=""), cont=cont,
                     change="_UP", dtype="_At_homologs", FDRmax=FDRmax,
                     logFCmin=logFCmin, dirPath=dir)
  write.table(subset(upDEGs, (!is.na(Symbol)), select="A._thaliana_best_hit"),
              fname, quote=FALSE, row.names=FALSE, na="", col.names=FALSE,
              sep=",")
  # DOWN GO lists
  fname <- makeFname(folder=paste(folder, "/GO_lists", sep=""), cont=cont,
                     change="_DOWN", dtype="_At_homologs", FDRmax=FDRmax,
                     logFCmin=logFCmin, dirPath=dir)
  write.table(subset(downDEGs, (!is.na(Symbol)), select="A._thaliana_best_hit"),
              fname, quote=FALSE, row.names=FALSE, na="", col.names=FALSE,
              sep=",")
  # UP Promoter sequences
  fname <- makeFname(folder=paste(folder, "/promoters", sep=""), cont=cont,
                     change="_UP", dtype="_promoters", FDRmax=FDRmax,
                     logFCmin=logFCmin, dirPath=dir)
  writeXStringSet(allPromoters2[names(allPromoters) %in% upDEGs$Gene_ID], fname)
  # DOWN Promoter sequences
  fname <- makeFname(folder=paste(folder, "/promoters", sep=""), cont=cont,
                     change="_DOWN", dtype="_promoters", FDRmax=FDRmax,
                     logFCmin=logFCmin, dirPath=dir)
  writeXStringSet(allPromoters2[names(allPromoters) %in% downDEGs$Gene_ID],
                  fname)
}



save_Venn_Lists <- function(allGenesList, folder, cont, FDRmax, logFCmin,
                            change, dir){
  # AllGenesList should have edinInfo, and group columns
  # change should be of form "_UP", or "_DOWN"

  allGenesList$Gene_ID <- row.names(allGenesList)

  # WT significant:
  fname <- makeFname(folder=folder, cont=paste(cont, "WT_SIG", sep="_"),
                     change=change, FDRmax=FDRmax, logFCmin=logFCmin, dirPath=dir)
  write.csv(subset(allGenesList, group=="WT_sig"), fname, row.names=TRUE)

  fname <- makeFname(folder=paste(folder, "/GO_lists", sep=""),
                     cont=paste(cont, "WT_SIG", sep="_"),
                     change=change, dtype="_At_homologs", FDRmax=FDRmax, logFCmin=logFCmin, dirPath=dir)
  write.table(subset(allGenesList, (group=="WT_sig")&(!is.na(Symbol)),
                     select="A._thaliana_best_hit"),
              fname, quote=FALSE, row.names=FALSE, na="", col.names=FALSE, sep=",")

  fname <- makeFname(folder=paste(folder, "/promoters", sep=""),
                     cont=paste(cont, "WT_SIG", sep="_"),
                     change=change, dtype="_promoters", FDRmax=FDRmax, logFCmin=logFCmin, dirPath=dir)
  writeXStringSet(allPromoters2[names(allPromoters) %in% row.names(allGenesList[allGenesList$group == "WT_sig", ])], fname)

  # phyB significant:
  fname <- makeFname(folder=folder, cont=paste(cont, "phyB_SIG", sep="_"),
                     change=change, FDRmax=FDRmax, logFCmin=logFCmin, dirPath=dir)
  write.csv(subset(allGenesList, group=="phyB_sig"), fname, row.names=TRUE)

  fname <- makeFname(folder=paste(folder, "/GO_lists", sep=""),
                     cont=paste(cont, "phyB_SIG", sep="_"),
                     change=change, dtype="_At_homologs", FDRmax=FDRmax, logFCmin=logFCmin, dirPath=dir)
  write.table(subset(allGenesList, (group=="phyB_sig")&(!is.na(Symbol)),
                     select="A._thaliana_best_hit"),
              fname, quote=FALSE, row.names=FALSE, na="", col.names=FALSE, sep=",")

  fname <- makeFname(folder=paste(folder, "/promoters", sep=""),
                     cont=paste(cont, "phyB_SIG", sep="_"),
                     change=change, dtype="_promoters", FDRmax=FDRmax, logFCmin=logFCmin, dirPath=dir)
  writeXStringSet(allPromoters2[names(allPromoters) %in% row.names(allGenesList[allGenesList$group == "phyB_sig", ])], fname)

  # both significant:
  fname <- makeFname(folder=folder, cont=paste(cont, "BOTH_SIG", sep="_"),
                     change=change, FDRmax=FDRmax, logFCmin=logFCmin, dirPath=dir)
  write.csv(subset(allGenesList, group=="Both_sig"), fname, row.names=TRUE)

  fname <- makeFname(folder=paste(folder, "/GO_lists", sep=""),
                     cont=paste(cont, "BOTH_SIG", sep="_"),
                     change=change, dtype="_At_homologs", FDRmax=FDRmax, logFCmin=logFCmin, dirPath=dir)
  write.table(subset(allGenesList, (group=="Both_sig")&(!is.na(Symbol)),
                     select="A._thaliana_best_hit"),
              fname, quote=FALSE, row.names=FALSE, na="", col.names=FALSE, sep=",")

  fname <- makeFname(folder=paste(folder, "/promoters", sep=""),
                     cont=paste(cont, "BOTH_SIG", sep="_"),
                     change=change, dtype="_promoters", FDRmax=FDRmax, logFCmin=logFCmin, dirPath=dir)
  writeXStringSet(allPromoters2[names(allPromoters) %in% row.names(allGenesList[allGenesList$group == "Both_sig", ])], fname)
}



# adds a column labeled "woundCtrl" with TRUE/FALSE if gene is in the woundIDs set
addWoundCtrl <- function(DEGs, woundIDs){
  output <- DEGs
  output$woundCtrl <- row.names(DEGs) %in% woundIDs
  return(output)
}

#' returns a merged DEG list dataframe
mergeGenes <- function(genesA, genesB){
  genesA$Gene_ID <- row.names(genesA)
  genesB$Gene_ID <- row.names(genesB)
  commonGenes <- dplyr::inner_join(genesA, genesB, by="Gene_ID")
  row.names(commonGenes) <- commonGenes$Gene_ID
  return(commonGenes)
}

#' returns common gene IDs between DEGsA and DEGsB
commonGenes <- function(DEGsA, DEGsB){
  commonIDs <- row.names(DEGsA)[row.names(DEGsA) %in% row.names(DEGsB)]
  DEGsA[commonIDs, ]
}

# standard error
se <- function(x){
  sd(x)/sqrt(length(x))
}



# plot a single gene's average expression at all time points
plotDE <- function(geneId, normExpr, display=TRUE, save=FALSE, saveDir="", fileType="svg", title=NULL, ...){
  group2Line <- Vectorize(function(group){
    switch(group,
           phyB.24 = 1,
           phyB.P = 1,
           phyB.R = 2,
           phyB.D = 2,
           WT.24 = 3,
           WT.P = 3,
           WT.R = 4,
           WT.D = 4)
  })
  normExpr <- data.frame(normExpr)
  group <- gsub("_\\d_", ".", colnames(normExpr))
  geneData <- as.data.frame(t(normExpr[geneId, ]))
  geneData$group <- group
  geneData <- dplyr::group_by(geneData, group)
  dataSummary <- dplyr::summarise(geneData, mean=mean(!!rlang::sym(geneId)),
                                  SE=se(!!rlang::sym(geneId)))
  dataSummary <- tidyr::separate(dataSummary, group,
                                 into=c("Genotype", "Timepoint"), sep="[.]",
                                 remove=FALSE)
  dataSummary$Timepoint <- factor(dataSummary$Timepoint, levels=c("P", "24", "D", "R"))
  dataSummary$line <- group2Line(dataSummary$group)
  if (!is.null(title)){
    plotTitle <- title
  } else{
    plotTitle <- geneId
  }
  # make plot
  plot <- ggplot(dataSummary, aes(x=Timepoint, y=mean, color=Genotype)) + geom_point(size=2) +
    geom_errorbar(aes(ymin=mean-SE, ymax=mean+SE), width=0.2) + geom_line(aes(group=line)) +
    expand_limits(y=0) + labs(title=plotTitle, y="Relative Expression")
  if (display) {
    print(plot)
  }
  if (save) {
    if (!file.exists(saveDir)) {
      stop("the save directory does not exist")
    }
    saveFile <- file.path(saveDir, paste(geneId, " expression plot.", fileType, sep=""))
    if (file.exists(saveFile)) {
      warning(paste(saveFile, "is being overwritten"))
    }
    ggsave(saveFile, plot, device=fileType, ...)
  }
}


geneHeatmap <- function(genes, DGEdata, groupAv=TRUE, rsOnly=TRUE,
                        labels=NULL){

  DEGCpm <- cpm(DGEdata[genes, ])
  sampleOrder <- c(
    "WT_3_P",   "WT_4_P",   "WT_5_P",   "WT_6_P",
    "phyB_3_P", "phyB_4_P", "phyB_5_P",
    "WT_3_24",  "WT_4_24",  "WT_5_24",  "WT_6_24",
    "phyB_3_24","phyB_4_24","phyB_5_24",
    "WT_3_D",   "WT_4_D",   "WT_5_D",   "WT_6_D",
    "phyB_3_D", "phyB_4_D", "phyB_5_D",
    "WT_3_R",   "WT_4_R",   "WT_5_R",   "WT_6_R",
    "phyB_3_R", "phyB_4_R", "phyB_5_R",
    "WT.D_P",   "phyB.D_P", "WT.R_D",   "phyB.R_D")
  sampleOrder <- sampleOrder[sampleOrder %in% colnames(DEGCpm)]
  DEGCpm <- DEGCpm[, sampleOrder]
  group <- gsub("_\\d_", ".", colnames(DEGCpm))
  if (groupAv) {
    DEGCpm2 <- as.data.frame(t(DEGCpm))
    DEGCpm2 <- stats::aggregate(DEGCpm2, by=list(group), FUN=mean)
    row.names(DEGCpm2) <- DEGCpm2$Group.1
    DEGCpm2 <- DEGCpm2[, -1]
    DEGCpm <- t(DEGCpm2)
    groupOrder <- c("WT.P", "phyB.P", "WT.24", "phyB.24", "WT.D", "phyB.D", "WT.R", "phyB.R")
    groupOrder <- groupOrder[groupOrder %in% colnames(DEGCpm)]
    DEGCpm <- DEGCpm[, groupOrder]
  }
  if (rsOnly){
    DEGCpm <- DEGCpm[, !grepl("P|(24)", colnames(DEGCpm))]
  }
  if(!is.null(labels)) {
    hMapLabels <- labels
  } else{
    hMapLabels <- row.names(DEGCpm)
  }
  hr <- hclust(as.dist(1-cor(t(DEGCpm), method="pearson")), method="complete")
  fig <- heatmap.2(DEGCpm, dendrogram="row", Rowv=as.dendrogram(hr), trace="none",
                   col="bluered", scale="row", Colv=NA, margins=c(8,7),
                   labRow=hMapLabels)
  return(fig)
}


doubleVenn <- function(DEG1, DEG2, name1="Set1", name2="Set2"){
  overlap <- nrow(commonGenes(DEG1, DEG2))

  draw.pairwise.venn(area1=nrow(DEG1), area2=nrow(DEG2),
                     cross.area=overlap,
                     fill=c("red", "blue"), category=c(name1, name2))
}

tripVenn <- function(DEG1, DEG2, DEG3, names=c("Set1", "Set2", "Set3")){
  genes12 <- commonGenes(DEG1, DEG2)
  genes23 <- commonGenes(DEG2, DEG3)
  genes13 <- commonGenes(DEG1, DEG3)
  genes123 <- commonGenes(DEG1, genes23)

  draw.triple.venn(area1=nrow(DEG1), area2=nrow(DEG2), area3=nrow(DEG3),
                   n12=nrow(genes12), n23=nrow(genes23), n13=nrow(genes13),
                   n123=nrow(genes123),
                   fill=c("red", "blue", "green"), category=names, scaled=FALSE)
}

nVenn <- function(DEGsList){
  DEGVenn <- Venn(DEGsList)
  gp <- VennThemes(compute.Venn(DEGVenn))
  gp$SetText <- lapply(gp$SetText,function(x) {x$fontsize<-10; return(x)})
  plot(DEGVenn, doWeights=FALSE,gp=gp)
}

# set group names for ven diagram
setGroup <- function(row, change="up", sigLvl=0.01, LFCutoff=1){
  if (change=="up"){
    subsetMask <- (row$logFC_WT > LFCutoff | row$logFC_phyB > LFCutoff) &
      (row$FDR_WT < sigLvl | row$FDR_phyB < sigLvl)
  }
  if (change=="down"){
    subsetMask <- (row$FDR_WT < sigLvl | row$FDR_phyB < sigLvl) &
      (row$logFC_WT < -LFCutoff | row$logFC_phyB < -LFCutoff)
  }
  wtSigTest <-row$FDR_WT < sigLvl & row$FDR_phyB > sigLvl
  phyBSigTest <- row$FDR_phyB < sigLvl & row$FDR_WT > sigLvl
  group <- "none"
  if (subsetMask) {
    group <- "sig"
    if (wtSigTest) {
      group <- "WT_sig"
    }else if (phyBSigTest) {
      group <- "phyB_sig"
    } else {
      group <- "Both_sig"
    }
  }
  return(group)
}

#' add average counts per million to a gene list/data frame where row names are
#' gene ids
addAvgCpm <- function(geneList, DGEdata, groups=NULL){

  genes <- row.names(geneList)

  DEGCpm <- cpm(DGEdata[genes, ])
  sampleOrder <- c(
    "WT_3_P",   "WT_4_P",   "WT_5_P",   "WT_6_P",
    "phyB_3_P", "phyB_4_P", "phyB_5_P",
    "WT_3_24",  "WT_4_24",  "WT_5_24",  "WT_6_24",
    "phyB_3_24","phyB_4_24","phyB_5_24",
    "WT_3_D",   "WT_4_D",   "WT_5_D",   "WT_6_D",
    "phyB_3_D", "phyB_4_D", "phyB_5_D",
    "WT_3_R",   "WT_4_R",   "WT_5_R",   "WT_6_R",
    "phyB_3_R", "phyB_4_R", "phyB_5_R",
    "WT.D_P",   "phyB.D_P", "WT.R_D",   "phyB.R_D")
  sampleOrder <- sampleOrder[sampleOrder %in% colnames(DEGCpm)]
  DEGCpm <- DEGCpm[, sampleOrder]

  if (class(DEGCpm)=="matrix") {
    group <- gsub("_\\d_", ".", colnames(DEGCpm))
    DEGCpm2 <- as.data.frame(t(DEGCpm))
  }else {
    group <- gsub("_\\d_", ".", names(DEGCpm))
    DEGCpm2 <- as.data.frame(DEGCpm)
    colnames(DEGCpm2) <- genes
  }
  DEGCpm2 <- stats::aggregate(DEGCpm2, by=list(group), FUN=mean)
  row.names(DEGCpm2) <- DEGCpm2$Group.1
  DEGCpm <- data.frame(t(DEGCpm2))

  groupOrder <- c("WT.P", "phyB.P", "WT.24", "phyB.24", "WT.D", "phyB.D", "WT.R", "phyB.R")
  groupOrder <- groupOrder[groupOrder %in% colnames(DEGCpm)]
  DEGCpm <- DEGCpm[-1, groupOrder]
  if (!is.null(groups)){
    DEGCpm <- DEGCpm[, groups]
  }
  colnames(DEGCpm) <- paste(colnames(DEGCpm), "_Avg_CPM", sep="")
  return(cbind(geneList, DEGCpm))

}
