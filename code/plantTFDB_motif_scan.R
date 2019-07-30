# pull motifs from plantTFDB data, search for those motifs in brassica promoter regions.

library(Biostrings)
library(GenomicRanges)
library(rprojroot)
library(TFBSTools)
library(motifmatchr)
library(stringr)
library(tibble)
library(openxlsx)
library(dplyr)
library(plyr)

rootDir <- find_root(is_rstudio_project)

tfdbFolder <- file.path(rootDir, "data", "plantTFDB")
infoFile <- file.path(tfdbFolder, "Ath_TF_binding_motifs_information.txt")
memeFolder <- file.path(tfdbFolder, "Ath_TF_binding_motifs_individual")

resultsFoldeToRead <- file.path(rootDir, "results/local", "analysis_output_20190726")


# ------------------------------------------------------------------------------
Files <- c(
# RvD UP________________________________________________________________________
  "RvD_WT_up" = file.path(resultsFoldeToRead,
                          "timepoint_v_timepoint/RvD_venn/promoters",
                          "RvD_WT_SIG_UP_promoters_FDR01_LFC01_20190726.fasta"),
  "RvD_phyB_up" = file.path(resultsFoldeToRead,
                            "timepoint_v_timepoint/RvD_venn/promoters",
                            "RvD_phyB_SIG_UP_promoters_FDR01_LFC01_20190726.fasta"),
  "RvD_both_up" = file.path(resultsFoldeToRead,
                            "timepoint_v_timepoint/RvD_venn/promoters",
                            "RvD_BOTH_SIG_UP_promoters_FDR01_LFC01_20190726.fasta"),
# RvD DOWN______________________________________________________________________
  "RvD_WT_down" = file.path(resultsFoldeToRead,
                            "timepoint_v_timepoint/RvD_venn/promoters",
                            "RvD_WT_SIG_DOWN_promoters_FDR01_LFC01_20190726.fasta"),
  "RvD_phyB_down" = file.path(resultsFoldeToRead,
                              "timepoint_v_timepoint/RvD_venn/promoters",
                              "RvD_phyB_SIG_DOWN_promoters_FDR01_LFC01_20190726.fasta"),
  "RvD_both_dpwm" = file.path(resultsFoldeToRead,
                              "timepoint_v_timepoint/RvD_venn/promoters",
                              "RvD_BOTH_SIG_DOWN_promoters_FDR01_LFC01_20190726.fasta"),
# DvP UP________________________________________________________________________
  "DvP_WT_up" = file.path(resultsFoldeToRead,
                          "timepoint_v_timepoint/DvP_venn/promoters",
                          "DvP_WT_SIG_UP_promoters_FDR01_LFC01_20190726.fasta"),
  "DvP_phyB_up" = file.path(resultsFoldeToRead,
                            "timepoint_v_timepoint/DvP_venn/promoters",
                            "DvP_phyB_SIG_UP_promoters_FDR01_LFC01_20190726.fasta"),
  "DvP_both_up" = file.path(resultsFoldeToRead,
                            "timepoint_v_timepoint/DvP_venn/promoters",
                            "DvP_BOTH_SIG_UP_promoters_FDR01_LFC01_20190726.fasta"),
# DvP DOWN______________________________________________________________________
  "DvP_WT_down" = file.path(resultsFoldeToRead,
                            "timepoint_v_timepoint/DvP_venn/promoters",
                            "DvP_WT_SIG_DOWN_promoters_FDR01_LFC01_20190726.fasta"),
  "DvP_phyB_down" = file.path(resultsFoldeToRead,
                              "timepoint_v_timepoint/DvP_venn/promoters",
                              "DvP_phyB_SIG_DOWN_promoters_FDR01_LFC01_20190726.fasta"),
  "DvP_both_down" = file.path(resultsFoldeToRead,
                              "timepoint_v_timepoint/DvP_venn/promoters",
                              "DvP_BOTH_SIG_DOWN_promoters_FDR01_LFC01_20190726.fasta"),
# WT vs. PHYB at PRE____________________________________________________________
  "phyB.P_up" = file.path(resultsFoldeToRead,
                          "phyB_v_WT_at_timepoint/promoters",
                          "phyB.P_vs_WT.P_exact_UP_promoters_FDR05_LFC01_20190726.fasta"),
  "phyB.P_down" = file.path(resultsFoldeToRead,
                            "phyB_v_WT_at_timepoint/promoters",
                            "phyB.P_vs_WT.P_exact_DOWN_promoters_FDR05_LFC01_20190726.fasta"),
# WT vs. PHYB at Recovery_______________________________________________________
"phyB.R_up" = file.path(resultsFoldeToRead,
                        "phyB_v_WT_at_timepoint/promoters",
                        "phyB.R_vs_WT.R_exact_UP_promoters_FDR05_LFC01_20190726.fasta"),
"phyB.R_down" = file.path(resultsFoldeToRead,
                          "phyB_v_WT_at_timepoint/promoters",
                          "phyB.R_vs_WT.R_exact_DOWN_promoters_FDR05_LFC01_20190726.fasta")
)

# make a list of genes w/ chloroplast GO cellular component term
chloroplastFile <- file.path(rootDir, "data", "GO_terms", "chloroplast_genes.txt")
conn <- file(chloroplastFile, open="r")
temp <- readLines(conn)
close(conn)
temp <- sub("^.*?=", "", temp)
temp <- sub("\\|.*", "", temp)
temp <- sub("g", "G", temp)
temp <- sub("t", "T", temp)
chloroplastGenes <- temp

# plantTFDB info on all A.t motifs
info <- read.table(infoFile, header=TRUE, sep="\t", stringsAsFactors = FALSE)

# get gene Info from Blast results
AtBlastPath <- file.path(rootDir, "data", "blast_top_hits.xlsx")
df <- read.xlsx(xlsxFile=AtBlastPath, sheet=2, skipEmptyRows=FALSE)
colnames(df) <- gsub("(\\.)([a-z])", "_\\2", colnames(df))
colnames(df)[1] <- "Gene_ID"
AtBlastGeneInfo <- df
info <- left_join(info, unique(AtBlastGeneInfo[, 2:4]), by=c("Gene_id"="A._thaliana_best_hit"))

# get motif cluster information
fname <- file.path(rootDir, "data/plantTFDB",
                   "matrix-clustering_2019-07-05.192043_zxAuzQ",
                   "matrix-clustering_tables", "clusters_motif_names.tab")
conn <- file(fname, "r")
lines <- readLines(conn)
close(conn)
df <- str_split(lines, "\t", simplify=TRUE)
df <- adply(df, 1, .fun=function(row){
  outputChunk <- data.frame("id"=str_split(row[2], ",")[[1]],
                            "cluster"= row[1], stringsAsFactors=FALSE)
}, .id=NULL)
info <- left_join(info, df, by=c("Gene_id"="id") )

# ------------------------------------------------------------------------------
# FUNCTIONS

readMeme <- function(fname) {
  conn <- file(fname, "r")
  lines <- readLines(conn)
  close(conn)
  # background
  background <- lines[8]
  background <- as.numeric(str_split(gsub("[A-Z] ", "", background), " ")[[1]][1:4])
  names(background) <- c("A", "C", "G", "T")
  ID <- str_match(lines[10], " (MP\\d*)")[2]
  # matrix
  matInfo <- lines[12]
  matInfo <- sub("letter-probability matrix: ", "", matInfo)
  matRows <- as.numeric(str_match(matInfo, " w= (\\d*)")[2])
  mat <- read.table(text=lines[13:(12+matRows)], sep="\t",
                    stringsAsFactors=FALSE)
  mat <- mat[, 1:4]
  colnames(mat) <- c("A", "C", "G", "T")
  mat <- as.matrix(mat)
  #deal with rounding error, matchMotifs() is sensitive to PWM rows adding to 1
  remainder <- 1- rowSums(mat)
  mat <- mat + remainder/4
  tempPWM <- PWMatrix(ID=ID, name=ID, bg=background, profileMatrix=t(mat))
  return(tempPWM)
}

loadPlantTFDB <- function(memeFolder, infoFile) {
  info <- read.table(infoFile, header=TRUE, sep="\t", stringsAsFactors = FALSE)
  xList <- list()
  for (i in 1:nrow(info)){
    fname <- file.path(memeFolder, paste(info[i, "Gene_id"], ".meme", sep=""))
    tempPWM <- readMeme(fname)
    xList[info[i, "Matrix_id"]] <- tempPWM
  }
  PWMList <- do.call(PWMatrixList, xList)
  return(PWMList)
}

fisherEnrich <- function(setA, setB, PWMList, promoterSeqs, set="all", namesAB=NULL){
  # setA and setB are A.t gene ID lists
  aPromoters <- setA
  bPromoters <- setB

  if (set=="chloro"){
    promoterAnnos <- str_split(names(aPromoters), " \\| ", simplify=TRUE)
    aPromoters <- aPromoters[(promoterAnnos[,2] %in% chloroplastGenes)]
    promoterAnnos <- str_split(names(bPromoters), " \\| ", simplify=TRUE)
    bPromoters <- bPromoters[(promoterAnnos[,2] %in% chloroplastGenes)]
  } else if (set=="non-chloro") {
    promoterAnnos <- str_split(names(aPromoters), " \\| ", simplify=TRUE)
    aPromoters <- aPromoters[!(promoterAnnos[,2] %in% chloroplastGenes)]
    promoterAnnos <- str_split(names(bPromoters), " \\| ", simplify=TRUE)
    bPromoters <- bPromoters[!(promoterAnnos[,2] %in% chloroplastGenes)]
  }

  output <- data.frame("motif_ID"=names(PWMList), stringsAsFactors=FALSE)
  fArray <- array(dim=c(length(PWMList), 2, 2), dimnames=list(names(PWMList), c("A","B"), c("w.Motif", "wo.Motif")))

  # A promoters
  results <- matchMotifs(PWMList, aPromoters)
  resultsMatrix <- motifMatches(results)
  motifCounts <- colSums(resultsMatrix)
  fArray[,1,1] <- motifCounts
  fArray[,1,2] <- length(aPromoters) - motifCounts # genes w/o motif counts
  output$aPercent <- motifCounts/length(aPromoters)*100

  # B promoters
  results <- matchMotifs(PWMList, bPromoters)
  resultsMatrix <- motifMatches(results)
  motifCounts <- colSums(resultsMatrix)
  fArray[,2,1] <- motifCounts
  fArray[,2,2] <- length(bPromoters) - motifCounts
  output$bPercent <- motifCounts/length(bPromoters)*100

  if (!is.null(namesAB)){
      colnames(output)[2:3] <- paste0("%_", namesAB, "_w_motif")
  }

  output$p.val <- NA
  for (i in 1:length(PWMList)) {
    output[i, "p.val"] <- fisher.test(fArray[i,,])$p.value
  }

  output$FDR <- p.adjust(output$p.val, method="hochberg")
  output <- output[order(output$p.val), ]
  return(output)
}

# ------------------------------------------------------------------------------
# MAIN CODE

TfdbMotifs <- loadPlantTFDB(memeFolder, infoFile)


df <- tribble(~name,      ~A,            ~B,              ~namesAB,
              "RvD_up",   "RvD_WT_up",   "RvD_phyB_up",   c("WT", "phyB"),
              "RvD_down", "RvD_WT_down", "RvD_phyB_down", c("WT", "phyB"),
              "DvP_up",   "DvP_WT_up",   "DvP_phyB_up",   c("WT", "phyB"),
              "DvP_down", "DvP_WT_down", "DvP_phyB_down", c("WT", "phyB"),
              "phyB.P_v_WT", "phyB.P_up",  "phyB.P_down",   c("phyB_up", "phyB_down"),
              "phyB.R_v_WT", "phyB.R_up",  "phyB.R_down",   c("phyB_up", "phyB_down")
)

xlWorkbk <- createWorkbook()
for (i in 1:nrow(df)){
  row <- df[i, ]
  aPromoters <- readDNAStringSet(Files[row$A])
  bPromoters <- readDNAStringSet(Files[row$B])
  # Genes w/ Chloroplast GO term
  enrich <- fisherEnrich(aPromoters, bPromoters, TfdbMotifs, allPromoters, set="chloro", namesAB=row$namesAB[[1]])
  enrich <- left_join(enrich, info[, c(1:3,10,8,9)], by=c("motif_ID"="Matrix_id"))
  name <- paste0(row$name, "_chloro_GO_term")
  addWorksheet(xlWorkbk, name)
  writeData(xlWorkbk, name, enrich)
  # Genes W/o Chloroplast GO terms
  enrich <- fisherEnrich(aPromoters, bPromoters, TfdbMotifs, allPromoters, set="non-chloro", namesAB=row$namesAB[[1]])
  enrich <- left_join(enrich, info[, c(1:3,10,8,9)], by=c("motif_ID"="Matrix_id"))
  name <- paste0(row$name, "_NO_chloro_GO_term")
  addWorksheet(xlWorkbk, name)
  writeData(xlWorkbk, name, enrich)
}
fname <- file.path(rootDir, "results/local", "motif_enrichment.xlsx")
saveWorkbook(xlWorkbk, file=fname, overwrite=TRUE)

