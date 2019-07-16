# generate some tables for the paper.


rootDir <- find_root(is_rstudio_project)
resultsFoldeToRead <- file.path(rootDir, "results", "analysis_output_20190510")
outputDir <- file.path(rootDir, "results", "tables")
dir.create(outputDir)

#load utility functions
source(file.path(rootDir, "code", "utilities.R"))

# load counts file =============================================================
countsFile <- file.path(rootDir, "data", "raw_counts.csv")
counts <- read.csv(countsFile, row.names="X")
# rename S (for sensence) to D (for dark-treatment)
colnames(counts) <- sub("_S", "_D", colnames(counts))
group <- gsub("_\\d_", ".", colnames(counts))
# create the DGEList object edgeR uses
DGEdata <- DGEList(counts=counts, group=group)
# filter out rows with fewer than 3 samples that have more than 10 counts/million
keep <- rowSums(cpm(DGEdata) > 10) >= 3
DGEdata <- DGEdata[keep, , keep.lib.sizes=FALSE]
# ==============================================================================


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

# files for the outer wedges of the Venn Diagram
Files <- c(
  "RvD_WT_up" = file.path(resultsFoldeToRead,
                          "timepoint_v_timepoint/RvD_venn",
                          "RvD_WT_SIG_UP_FDR01_LFC01_20190510.csv"),
  "RvD_BOTH_up" = file.path(resultsFoldeToRead,
                            "timepoint_v_timepoint/RvD_venn",
                            "RvD_BOTH_SIG_UP_FDR01_LFC01_20190510.csv"),
  "RvD_phyB_up" = file.path(resultsFoldeToRead,
                            "timepoint_v_timepoint/RvD_venn",
                            "RvD_phyB_SIG_UP_FDR01_LFC01_20190510.csv"),

  "RvD_WT_down" = file.path(resultsFoldeToRead,
                            "timepoint_v_timepoint/RvD_venn",
                            "RvD_WT_SIG_DOWN_FDR01_LFC01_20190510.csv"),
  "RvD_BOTH_down" = file.path(resultsFoldeToRead,
                            "timepoint_v_timepoint/RvD_venn",
                            "RvD_BOTH_SIG_DOWN_FDR01_LFC01_20190510.csv"),
  "RvD_phyB_down" = file.path(resultsFoldeToRead,
                              "timepoint_v_timepoint/RvD_venn",
                              "RvD_phyB_SIG_DOWN_FDR01_LFC01_20190510.csv"),

  "DvP_WT_up" = file.path(resultsFoldeToRead,
                          "timepoint_v_timepoint/DvP_venn",
                          "DvP_WT_SIG_UP_FDR01_LFC01_20190510.csv"),
  "DvP_BOTH_up" = file.path(resultsFoldeToRead,
                          "timepoint_v_timepoint/DvP_venn",
                          "DvP_BOTH_SIG_UP_FDR01_LFC01_20190510.csv"),
  "DvP_phyB_up" = file.path(resultsFoldeToRead,
                            "timepoint_v_timepoint/DvP_venn",
                            "DvP_phyB_SIG_UP_FDR01_LFC01_20190510.csv"),

  "DvP_WT_down" = file.path(resultsFoldeToRead,
                            "timepoint_v_timepoint/DvP_venn",
                            "DvP_WT_SIG_DOWN_FDR01_LFC01_20190510.csv"),
  "DvP_WT_down" = file.path(resultsFoldeToRead,
                            "timepoint_v_timepoint/DvP_venn",
                            "DvP_BOTH_SIG_DOWN_FDR01_LFC01_20190510.csv"),
  "DvP_phyB_down" = file.path(resultsFoldeToRead,
                              "timepoint_v_timepoint/DvP_venn",
                              "DvP_phyB_SIG_DOWN_FDR01_LFC01_20190510.csv")
)

promoterFiles <- c(
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
                            "phyB.P_vs_WT.P_exact_DOWN_promoters_FDR05_LFC01_20190510.fasta")
)

# RvD_up_files <- Files[1:3]
# RvD_down_files <- Files[4:6]
# DvP_up_files <- Files[7:9]
# DvP_downFiles <- Files[10:12]

fileLists <- list("RvD_up"=Files[1:3], "RvD_down"=Files[4:6],
                  "DvP_up"=Files[7:9], "DvP_down"=Files[10:12])




for (i in 1:4){
  # make blank DF
  outputTable <- data.frame()
  for(j in 1:3){
    file <- fileLists[[i]][j]
    rawTable <- read.table(file, header=TRUE, sep=",", row.names=1)
    # subset to chloroplast genes
    tempTable <- rawTable
    tempTable$sig_conditions <- sub("^.*?_", "", names(fileLists[[i]][j]))
    tempTable <- subset(tempTable, A._thaliana_best_hit %in% chloroplastGenes,
                        select=c("Gene_ID",
                                 "logFC_WT", "logFC_phyB",
                                 "logCPM_WT", "logCPM_phyB",
                                 "PValue_WT", "PValue_phyB",
                                 "FDR_WT", "FDR_phyB",
                                 "sig_conditions",
                                 "A._thaliana_best_hit", "Symbol", "Description"))

    # append new stuff to DF
    outputTable <- rbind(outputTable, tempTable)
  }
  # add average cpm values for all time points
  outputTable <- addAvgCpm(outputTable, DGEdata)
  # save DF
  fname <- paste(names(fileLists[i]), "_chloroplast_related.csv", sep="")
  fileOutName <- file.path(outputDir, fname)
  write.table(outputTable, fileOutName, sep=",", row.names=FALSE)
}

# for (i in 1:length(Files)){
#   file <- Files[i]
#   rawTable <- read.table(file, header=TRUE, sep=",", row.names=1)
#   # subset to chloroplast genes
#   tempTable <- subset(rawTable, A._thaliana_best_hit %in% chloroplastGenes,
#                       select=c("Gene_ID",
#                                "logFC_WT", "logFC_phyB",
#                                "logCPM_WT", "logCPM_phyB",
#                                "PValue_WT", "PValue_phyB",
#                                "FDR_WT", "FDR_phyB",
#                                "A._thaliana_best_hit", "Symbol", "Description"))
#
#   # add average cpm values for all time points
#   tempTable <- addAvgCpm(tempTable, DGEdata)
#
#   # save table
#   fname <- sub("\\.csv", "_chloroplast_related.csv", sub(".*\\/", "", Files[i]))
#   fileOutName <- file.path(outputDir, fname)
#   write.table(tempTable, fileOutName, sep=",", row.names=FALSE)
#
#
# }
#






