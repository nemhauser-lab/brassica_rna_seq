# pull motifs from plantTFDB data, search for those motifs in brassica promoter regions.

library(rprojroot)
library(TFBSTools)
library(motifmatchr)
library(stringr)
library(tibble)
library(openxlsx)

rootDir <- find_root(is_rstudio_project)

tfdbFolder <- file.path(rootDir, "data", "plantTFDB")
infoFile <- file.path(tfdbFolder, "Ath_TF_binding_motifs_information.txt")
memeFolder <- file.path(tfdbFolder, "Ath_TF_binding_motifs_individual")

resultsFoldeToRead <- file.path(rootDir, "results", "analysis_output_20190510")


# ------------------------------------------------------------------------------
Files <- c(
# RvD UP
  "RvD_WT_up" = file.path(resultsFoldeToRead,
                          "timepoint_v_timepoint/RvD_venn/promoters",
                          "RvD_WT_SIG_UP_promoters_FDR01_LFC01_20190510.fasta"),
  "RvD_phyB_up" = file.path(resultsFoldeToRead,
                            "timepoint_v_timepoint/RvD_venn/promoters",
                            "RvD_phyB_SIG_UP_promoters_FDR01_LFC01_20190510.fasta"),
  "RvD_both_up" = file.path(resultsFoldeToRead,
                            "timepoint_v_timepoint/RvD_venn/promoters",
                            "RvD_BOTH_SIG_UP_promoters_FDR01_LFC01_20190510.fasta"),
# RvD DOWN
  "RvD_WT_down" = file.path(resultsFoldeToRead,
                            "timepoint_v_timepoint/RvD_venn/promoters",
                            "RvD_WT_SIG_DOWN_promoters_FDR01_LFC01_20190510.fasta"),
  "RvD_phyB_down" = file.path(resultsFoldeToRead,
                              "timepoint_v_timepoint/RvD_venn/promoters",
                              "RvD_phyB_SIG_DOWN_promoters_FDR01_LFC01_20190510.fasta"),
  "RvD_both_dpwm" = file.path(resultsFoldeToRead,
                              "timepoint_v_timepoint/RvD_venn/promoters",
                              "RvD_BOTH_SIG_DOWN_promoters_FDR01_LFC01_20190510.fasta"),
# DvP UP
  "DvP_WT_up" = file.path(resultsFoldeToRead,
                          "timepoint_v_timepoint/DvP_venn/promoters",
                          "DvP_WT_SIG_UP_promoters_FDR01_LFC01_20190510.fasta"),
  "DvP_phyB_up" = file.path(resultsFoldeToRead,
                            "timepoint_v_timepoint/DvP_venn/promoters",
                            "DvP_phyB_SIG_UP_promoters_FDR01_LFC01_20190510.fasta"),
  "DvP_both_up" = file.path(resultsFoldeToRead,
                            "timepoint_v_timepoint/DvP_venn/promoters",
                            "DvP_BOTH_SIG_UP_promoters_FDR01_LFC01_20190510.fasta"),
# DvP DOWN
  "DvP_WT_down" = file.path(resultsFoldeToRead,
                            "timepoint_v_timepoint/DvP_venn/promoters",
                            "DvP_WT_SIG_DOWN_promoters_FDR01_LFC01_20190510.fasta"),
  "DvP_phyB_down" = file.path(resultsFoldeToRead,
                              "timepoint_v_timepoint/DvP_venn/promoters",
                              "DvP_phyB_SIG_DOWN_promoters_FDR01_LFC01_20190510.fasta"),
  "DvP_both_down" = file.path(resultsFoldeToRead,
                              "timepoint_v_timepoint/DvP_venn/promoters",
                              "DvP_BOTH_SIG_DOWN_promoters_FDR01_LFC01_20190510.fasta"),
# WT vs. PHYB at PRE
  "phyB.P_up" = file.path(resultsFoldeToRead,
                          "phyB_v_WT_at_timepoint/promoters",
                          "phyB.P_vs_WT.P_exact_UP_promoters_FDR05_LFC01_20190510.fasta"),
  "phyB.P_down" = file.path(resultsFoldeToRead,
                            "phyB_v_WT_at_timepoint/promoters",
                            "phyB.P_vs_WT.P_exact_DOWN_promoters_FDR05_LFC01_20190510.fasta"),

# WT vs. PHYB at Recovery
"phyB.R_up" = file.path(resultsFoldeToRead,
                        "phyB_v_WT_at_timepoint/promoters",
                        "phyB.R_vs_WT.R_exact_UP_promoters_FDR05_LFC01_20190510.fasta"),
"phyB.R_down" = file.path(resultsFoldeToRead,
                          "phyB_v_WT_at_timepoint/promoters",
                          "phyB.R_vs_WT.R_exact_DOWN_promoters_FDR05_LFC01_20190510.fasta")
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

xList <- list()
for (i in 1:nrow(info)){
  fname <- file.path(memeFolder, paste(info[i, "Gene_id"], ".meme", sep=""))
  tempPWM <- readMeme(fname)
  xList[info[i, "Matrix_id"]] <- tempPWM
}
PWMList <- do.call(PWMatrixList, xList)

# results <- matchMotifs(PWMList, promoters)
# resultsMatrix <- motifMatches(results)
# motifCounts <- colSums(resultsMatrix)

# ------------------------------------------------------------------------------
# Fisher's test (could be made into function)
# inputs: A & B file names, PWMList, info

A <- "RvD_WT_down"
B <- "RvD_phyB_down"

A <- "DvP_WT_down"
B <- "DvP_phyB_down"

A <- "phyB.P_up"
B <- "phyB.P_down"

A <- "phyB.R_up"
B <- "phyB.R_down"

# set
set <- "chloro"
# set <- "non-chloro"
# ...

aPromoters <- readDNAStringSet(Files[A])
bPromoters <- readDNAStringSet(Files[B])

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


output <- data.frame("motif_ID"=names(PWMList))

fArray <- array(dim=c(length(PWMList), 2, 2), dimnames=list(names(PWMList), c(A,B), c("w.Motif", "wo.Motif")))

# A promoters
results <- matchMotifs(PWMList, aPromoters)
resultsMatrix <- motifMatches(results)
motifCounts <- colSums(resultsMatrix)
fArray[,1,1] <- motifCounts
fArray[,1,2] <- length(aPromoters) - motifCounts # genes w/o motif counts
output$aRatio <- motifCounts/length(aPromoters)

# B promoters
results <- matchMotifs(PWMList, bPromoters)
resultsMatrix <- motifMatches(results)
motifCounts <- colSums(resultsMatrix)
fArray[,2,1] <- motifCounts
fArray[,2,2] <- length(bPromoters) - motifCounts
output$bRatio <- motifCounts/length(bPromoters)

# Fisher test for each motif
output$p.val <- NA
for (i in 1:length(PWMList)) {
  output[i, "p.val"] <- fisher.test(fArray[i,,])$p.value
}

output$FDR <- p.adjust(output$p.val, method="hochberg")
output <- cbind(output, info)

# ------------------------------------------------------------------------------
# add At gene info
# get gene Info from Blast results
AtBlastPath <- file.path(rootDir, "data", "blast_top_hits.xlsx")
df <- read.xlsx(xlsxFile=AtBlastPath, sheet=2, skipEmptyRows=FALSE)
colnames(df) <- gsub("(\\.)([a-z])", "_\\2", colnames(df))
colnames(df)[1] <- "Gene_ID"
AtBlastGeneInfo <- df

output2 <- dplyr::left_join(output, unique(AtBlastGeneInfo[, 2:4]), by=c("Gene_id"="A._thaliana_best_hit"))
output2 <- output2[order(output2$p.val), ]

# clipr::write_clip(output2)
hist(output2$p.val, breaks=40)

# ==============================================================================
# clusteredMotifs:
fname <- file.path(rootDir, "data", "plantTFDB",
                   "matrix-clustering_2019-07-05.192043_zxAuzQ",
                   "matrix-clustering_cluster_root_motifs.tf")

curratedClusters <- c("cluster_2"= "node_8/cluster_2_node_8_trimmed_matrices.tf",
                      "cluster_35"= "node_2/cluster_35_node_2_trimmed_matrices.tf",
                      "cluster_1"= "node_24/cluster_1_node_24_trimmed_matrices.tf",
                      "cluster_1.2"= "node_22/cluster_1_node_22_trimmed_matrices.tf")

readTf <- function(fname){
  conn <- file(fname, "r")
  lines <- readLines(conn)
  close(conn)

  splits <- lines=="//"
  chunk <- rep(0, length(lines))
  currentChunk <- 1
  for (i in 1:length(lines)) {
    chunk[i] <- currentChunk
    if (lines[i] == "//"){
      currentChunk <- currentChunk + 1
    }
  }
  chunkedLines <- split(lines, chunk)


  clustList <- list()
  clust_IDs <- c()
  # process each root
  for (i in 1:length(chunkedLines)){
    clustLines <- chunkedLines[[i]]
    # get cluster name
    clustName <- sub("AC  ", "", clustLines[1])

    if (clustName %in% names(curratedClusters)) {
      fname <- file.path(rootDir, "data", "plantTFDB", "matrix-clustering_2019-07-05.192043_zxAuzQ",
                "matrix-clustering_clusters_information", clustName,
                "merged_consensuses", curratedClusters[clustName])
      conn <- file(fname, "r")
      clustLines <- readLines(conn)
      close(conn)
      clustName <- sub("AC  ", "", clustLines[1])
    }

    matrixLines <- clustLines[grepl("^\\d", clustLines)]
    mat <- str_split(matrixLines, " +", simplify = TRUE)[, 2:5]
    colnames(mat) <- c("A", "C", "G", "T")
    class(mat) <- "numeric"

    clustList[clustName] <- toPWM( PFMatrix(ID=clustName, name=clustName,
                                            bg=background, profileMatrix=t(mat)) )
    clust_ID <- sub("CC  merged_ID: ", "",
                     clustLines[grepl("CC  merged_ID", clustLines)] )
    clust_ID <- gsub("'", "", clust_IDs)
    clust_IDs[clustName] <- clust_ID
  }
  PwmClustList <- do.call(PWMatrixList, clustList)
}

PWMList <- readTf(fname)
# then re-run Fisher's test section

head(output[order(output$p.val), ], 56)

