# calculate stats on the presense of G-Box 6-mers (CACGTG) in promoter sequences

library(Biostrings)
library(GenomicRanges)
library(rprojroot)
library(plyr)
library(stringr)
library(tibble)

rootDir <- find_root(is_rstudio_project)
resultsFoldeToRead <- file.path(rootDir, "results", "analysis_output_20190510")

Files <- c(
  "phyB.P_up" = file.path(resultsFoldeToRead,
                          "phyB_v_WT_at_timepoint/promoters",
                          "phyB.P_vs_WT.P_exact_DOWN_promoters_FDR05_LFC01_20190510.fasta"),
  "phyB.P_down" = file.path(resultsFoldeToRead,
                            "phyB_v_WT_at_timepoint/promoters",
                            "phyB.P_vs_WT.P_exact_UP_promoters_FDR05_LFC01_20190510.fasta"),

  # R v D_______________________________________________________________________
  "RvD WT_up" = file.path(resultsFoldeToRead,
                          "timepoint_v_timepoint/RvD_venn/promoters",
                          "RvD_WT_SIG_UP_promoters_FDR01_LFC01_20190510.fasta"),
  "RvD phyB_up" = file.path(resultsFoldeToRead,
                            "timepoint_v_timepoint/RvD_venn/promoters",
                            "RvD_phyB_SIG_UP_promoters_FDR01_LFC01_20190510.fasta"),
  "RvD both_up" = file.path(resultsFoldeToRead,
                            "timepoint_v_timepoint/RvD_venn/promoters",
                            "RvD_BOTH_SIG_UP_promoters_FDR01_LFC01_20190510.fasta"),

  "RvD WT_down" = file.path(resultsFoldeToRead,
                          "timepoint_v_timepoint/RvD_venn/promoters",
                          "RvD_WT_SIG_DOWN_promoters_FDR01_LFC01_20190510.fasta"),
  "RvD phyB_down" = file.path(resultsFoldeToRead,
                            "timepoint_v_timepoint/RvD_venn/promoters",
                            "RvD_phyB_SIG_DOWN_promoters_FDR01_LFC01_20190510.fasta"),
  "RvD both_down" = file.path(resultsFoldeToRead,
                            "timepoint_v_timepoint/RvD_venn/promoters",
                            "RvD_BOTH_SIG_DOWN_promoters_FDR01_LFC01_20190510.fasta"),

  # D v P_______________________________________________________________________
  "DvP WT_up" = file.path(resultsFoldeToRead,
                          "timepoint_v_timepoint/DvP_venn/promoters",
                          "DvP_WT_SIG_UP_promoters_FDR01_LFC01_20190510.fasta"),
  "DvP phyB_up" = file.path(resultsFoldeToRead,
                            "timepoint_v_timepoint/DvP_venn/promoters",
                            "DvP_phyB_SIG_UP_promoters_FDR01_LFC01_20190510.fasta"),
  "DvP both_up" = file.path(resultsFoldeToRead,
                            "timepoint_v_timepoint/DvP_venn/promoters",
                            "DvP_BOTH_SIG_UP_promoters_FDR01_LFC01_20190510.fasta"),

  "DvP WT_down" = file.path(resultsFoldeToRead,
                            "timepoint_v_timepoint/DvP_venn/promoters",
                            "DvP_WT_SIG_DOWN_promoters_FDR01_LFC01_20190510.fasta"),
  "DvP phyB_down" = file.path(resultsFoldeToRead,
                              "timepoint_v_timepoint/DvP_venn/promoters",
                              "DvP_phyB_SIG_DOWN_promoters_FDR01_LFC01_20190510.fasta"),
  "DvP both_down" = file.path(resultsFoldeToRead,
                              "timepoint_v_timepoint/DvP_venn/promoters",
                              "DvP_BOTH_SIG_DOWN_promoters_FDR01_LFC01_20190510.fasta"),

  # R v P_______________________________________________________________________
  "RvP WT_up" = file.path(resultsFoldeToRead,
                          "timepoint_v_timepoint/RvP_venn/promoters",
                          "RvP_WT_SIG_UP_promoters_FDR01_LFC01_20190510.fasta"),
  "RvP phyB_up" = file.path(resultsFoldeToRead,
                            "timepoint_v_timepoint/RvP_venn/promoters",
                            "RvP_phyB_SIG_UP_promoters_FDR01_LFC01_20190510.fasta"),
  "RvP both_up" = file.path(resultsFoldeToRead,
                            "timepoint_v_timepoint/RvP_venn/promoters",
                            "RvP_BOTH_SIG_UP_promoters_FDR01_LFC01_20190510.fasta"),

  "RvP WT_down" = file.path(resultsFoldeToRead,
                            "timepoint_v_timepoint/RvP_venn/promoters",
                            "RvP_WT_SIG_DOWN_promoters_FDR01_LFC01_20190510.fasta"),
  "RvP phyB_down" = file.path(resultsFoldeToRead,
                              "timepoint_v_timepoint/RvP_venn/promoters",
                              "RvP_phyB_SIG_DOWN_promoters_FDR01_LFC01_20190510.fasta"),
  "RvP both_down" = file.path(resultsFoldeToRead,
                              "timepoint_v_timepoint/RvP_venn/promoters",
                              "RvP_BOTH_SIG_DOWN_promoters_FDR01_LFC01_20190510.fasta"),

  # WT vs phyB at R_____________________________________________________________
  "phyB.R_up" = file.path(resultsFoldeToRead,
                          "phyB_v_WT_at_timepoint/promoters",
                          "phyB.R_vs_WT.R_exact_DOWN_promoters_FDR05_LFC01_20190510.fasta"),
  "phyB.R_down" = file.path(resultsFoldeToRead,
                            "phyB_v_WT_at_timepoint/promoters",
                            "phyB.R_vs_WT.R_exact_UP_promoters_FDR05_LFC01_20190510.fasta")

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
# unique(substr(temp2, 1,4))  # check that all genes are of form AT_G
# should be 5248 chloroplast genes total

# try using different source fro chloroplast GO terms, these should be identical.
chloroplastFile <- file.path(rootDir, "data", "GO_terms", "amiGOterm0009507_genes_products.txt")
conn <- file(chloroplastFile, open="r")
temp <- readLines(conn)
close(conn)
temp2 <- str_split(temp, "\t", simplify=TRUE)
tempGenes <- str_extract(temp2[,2:4], "A(T|t).(G|g)\\d{5}")
tempGenes <- tempGenes[!is.na(tempGenes)]
tempGenes <- sub("g", "G", tempGenes)
tempGenes <- sub("t", "T", tempGenes)
tempGenes <- unique(tempGenes)
# unique(substr(tempGenes, 1, 5))
misfits <- temp2[!( grepl("A(T|t).(G|g)\\d{5}", temp2[,2]) |
                     grepl("A(T|t).(G|g)\\d{5}", temp2[,3]) |
                     grepl("A(T|t).(G|g)\\d{5}", temp2[,4]) ),  ]
nonMisfits <- temp2[( grepl("A(T|t).(G|g)\\d{5}", temp2[,2]) |
                         grepl("A(T|t).(G|g)\\d{5}", temp2[,3]) |
                         grepl("A(T|t).(G|g)\\d{5}", temp2[,4]) ),  ]
# conclusion: lists are not identical, but overlap by >96%




# # 189 photosynthesis genes total, 182 of them in the chloroplast set
# photosynFile <- file.path(rootDir, "data", "GO_terms", "photosynthesis_genes.txt")
# conn <- file(photosynFile, open="r")
# temp <- readLines(conn)
# close(conn)
# temp2 <- sub("^.*?=", "", temp)
# temp2 <- sub("\\|.*", "", temp2)
# temp2 <- sub("g", "G", temp2)
# temp2 <- sub("t", "T", temp2)
# photosynthesisGenes <- temp2

# ==============================================================================
# functions

makeResTable <- function(Files, chloroplastGenes, motif){
  # chloroplast genes
  resTable <- data.frame()
  for (i in 1:length(Files)){
    file <- Files[i]
    promoters <- readDNAStringSet(file)
    promoterAnnos <- str_split(names(promoters), " \\| ", simplify=TRUE)
    promoters2 <- promoters[(promoterAnnos[,2] %in% chloroplastGenes)]
    promoterStrs <- as.character(promoters2)
    counts <- str_count(promoterStrs, motif)
    row <- data.frame(
      "name" = names(Files)[i],
      "nGenes" = length(counts),
      "genes_w_motif" = sum(counts != 0),
      "Motifs" = sum(counts),
      "percent_w_motif" = sum(counts != 0)/length(counts),
      "ave_motif_per_gene" = sum(counts)/sum(counts != 0),
      "subset" = "Chloroplast"
    )
    resTable <- rbind(resTable, row)
  }
  # Non-Chloroplast
  for (i in 1:length(Files)){
    file <- Files[i]
    promoters <- readDNAStringSet(file)
    promoterAnnos <- str_split(names(promoters), " \\| ", simplify=TRUE)
    promoters2 <- promoters[!(promoterAnnos[,2] %in% chloroplastGenes)]
    promoterStrs <- as.character(promoters2)
    counts <- str_count(promoterStrs, "CACGTG")
    row <- data.frame(
      "name" = names(Files)[i],
      "nGenes" = length(counts),
      "genes_w_motif" = sum(counts != 0),
      "Motifs" = sum(counts),
      "percent_w_motif" = sum(counts != 0)/length(counts),
      "ave_motif_per_gene" = sum(counts)/sum(counts != 0),
      "subset" = "Non-Chloroplast"
    )
    resTable <- rbind(resTable, row)
  }
  return(resTable)
}

plotBars <- function(resTable, title){
  resTable$genes_wo_motif <- resTable$nGenes-resTable$genes_w_motif
  resTable <- tidyr::separate(resTable, name, into=c("category", "sign"), sep="_")
  # resTable$category <- sub("(_up|_down)", "", resTable$name)

  resTable[resTable$sign=="down", "percent_w_motif"] <- -resTable[resTable$sign=="down", "percent_w_motif"]
  # resTable[c(4:6, 10:12), c("genes_w_motif", "genes_wo_motif", "percent_w_motif")] <- -resTable[c(4:6, 10:12), c("genes_w_motif", "genes_wo_motif", "percent_w_motif")]

  ggplot(resTable, aes(fill=subset, x=category, y=percent_w_motif*100)) +
    geom_bar(stat="identity", position="dodge", colour="black") +
    geom_text(aes(x=category, y=(1-2*(sign=="down"))*abs(percent_w_motif*50),
                  label=paste0(abs(round(percent_w_motif*100)),
                               "%\n(", abs(genes_w_motif), " / ",
                               nGenes, " genes)") ),
              position=position_dodge(width=1)) +
    coord_flip() +  scale_fill_manual(values=c("#99e354", "#54e3dc")) +
    labs(title=title, y="% of genes with motif (negative = down expressed)", x="significance category")

}

rawCount <- function(file, motif){
  promoters <- readDNAStringSet(file)
  promoterStrs <- as.character(promoters)
  counts <- str_count(promoterStrs, motif)
  output <- data.frame("counts"=counts, "gene_info"=names(promoters) )
  output <- tidyr::separate(output, gene_info,
                            into=c("gene_ID", "At_homolog", "short_name", "desription"),
                            sep=" \\| ")
}


# ==============================================================================
# look for G-boxes





# make table for all genes G-box presence
resTable <- data.frame()
for (i in 3:5){
  file <- Files[i]
  promoters <- readDNAStringSet(file)
  promoterStrs <- as.character(promoters)
  counts <- str_count(promoterStrs, "CACGTG")
  row <- data.frame(
    "name" = names(Files)[i],
    "nGenes" = length(counts),
    "genes_w_Gbox" = sum(counts != 0),
    "Gboxes" = sum(counts),
    "percent_w_Gbox" = sum(counts != 0)/length(counts),
    "ave_Gbox_per_gene" = sum(counts)/sum(counts != 0)
  )
  resTable <- rbind(resTable, row)
}

# do analysis for chloroplast related genes.
# 5248 chloroplast genes total
resTable <- data.frame()
for (i in 3:8){
  file <- Files[i]
  promoters <- readDNAStringSet(file)
  promoterAnnos <- str_split(names(promoters), " \\| ", simplify=TRUE)
  promoters2 <- promoters[(promoterAnnos[,2] %in% chloroplastGenes)]
  promoterStrs <- as.character(promoters2)
  counts <- str_count(promoterStrs, "CACGTG")
  row <- data.frame(
    "name" = names(Files)[i],
    "nGenes" = length(counts),
    "genes_w_motif" = sum(counts != 0),
    "Motifs" = sum(counts),
    "percent_w_motif" = sum(counts != 0)/length(counts),
    "ave_motif_per_gene" = sum(counts)/sum(counts != 0),
    "subset" = "Chloroplast"
  )
  resTable <- rbind(resTable, row)
}
# Non-Chloroplast
for (i in 3:8){
  file <- Files[i]
  promoters <- readDNAStringSet(file)
  promoterAnnos <- str_split(names(promoters), " \\| ", simplify=TRUE)
  promoters2 <- promoters[!(promoterAnnos[,2] %in% chloroplastGenes)]
  promoterStrs <- as.character(promoters2)
  counts <- str_count(promoterStrs, "CACGTG")
  row <- data.frame(
    "name" = names(Files)[i],
    "nGenes" = length(counts),
    "genes_w_motif" = sum(counts != 0),
    "Motifs" = sum(counts),
    "percent_w_motif" = sum(counts != 0)/length(counts),
    "ave_motif_per_gene" = sum(counts)/sum(counts != 0),
    "subset" = "Non-Chloroplast"
  )
  resTable <- rbind(resTable, row)
}


#===============================================================================
# evening element: AAAATATCT

# Chloroplast genes
resTable2 <- data.frame()
for (i in 3:8){
  file <- Files[i]
  promoters <- readDNAStringSet(file)
  promoterAnnos <- str_split(names(promoters), " \\| ", simplify=TRUE)
  promoters2 <- promoters[(promoterAnnos[,2] %in% chloroplastGenes)]
  promoterStrs <- as.character(promoters2)
  counts <- str_count(promoterStrs, "AAAATATCT|AGATATTTT") # Evenig element
  #counts <- str_count(promoterStrs, "AAAAAATCT|AGATTTTTT") # CCA1/morning element
  row <- data.frame(
    "name" = names(Files)[i],
    "nGenes" = length(counts),
    "genes_w_motif" = sum(counts != 0),
    "Motifs" = sum(counts),
    "percent_w_motif" = sum(counts != 0)/length(counts),
    "ave_motif_per_gene" = sum(counts)/sum(counts != 0),
    "subset" = "Chloroplast"
  )
  resTable2 <- rbind(resTable2, row)
}

# NON-Chloroplast genes
for (i in 3:8){
  file <- Files[i]
  promoters <- readDNAStringSet(file)
  promoterAnnos <- str_split(names(promoters), " \\| ", simplify=TRUE)
  promoters2 <- promoters[!(promoterAnnos[,2] %in% chloroplastGenes)]
  promoterStrs <- as.character(promoters2)
  counts <- str_count(promoterStrs, "AAAATATCT|AGATATTTT") # Evenig element
  # counts <- str_count(promoterStrs, "AAAAAATCT|AGATTTTTT") # CCA1/morning element
  row <- data.frame(
    "name" = names(Files)[i],
    "nGenes" = length(counts),
    "genes_w_motif" = sum(counts != 0),
    "Motifs" = sum(counts),
    "percent_w_motif" = sum(counts != 0)/length(counts),
    "ave_motif_per_gene" = sum(counts)/sum(counts != 0),
    "subset" = "Non-Chloroplast"
  )
  resTable2 <- rbind(resTable2, row)
}




###############################################################################

# G-BOX
resTable <- makeResTable(Files[3:8], chloroplastGenes, "CACGTG")
plotBars(resTable, "G-Box: CACGTG")
resTable <- makeResTable(Files[9:14], chloroplastGenes, "CACGTG")
plotBars(resTable, "G-Box: CACGTG")
# get genes
output <- rawCount(Files["RvD WT_up"], "CACGTG")
temp1 <- subset(output, counts>0 & (At_homolog %in% chloroplastGenes))

# EVENING ELEMENT
resTable <- makeResTable(Files[3:8], chloroplastGenes, "AAAATATCT|AGATATTTT")
plotBars(resTable, "Evening Element AAAATATCT")
resTable <- makeResTable(Files[9:14], chloroplastGenes, "AAAATATCT|AGATATTTT")
plotBars(resTable, "Evening Element AAAATATCT")
# get genes
output <- rawCount(Files["RvD phyB_down"], "AAAATATCT|AGATATTTT")
temp1 <- subset(output, counts>0 & (At_homolog %in% chloroplastGenes))

# TELOBOX
resTable <- makeResTable(Files[3:8], chloroplastGenes, "AAACCCTAA|TTAGGGTTT")
plotBars(resTable, "Telobox motif AAACCCTAA")
resTable <- makeResTable(Files[9:14], chloroplastGenes, "AAACCCTAA|TTAGGGTTT")
plotBars(resTable, "telobox motif AAACCCTAA")
# get telobox genes
# WT down
output <- rawCount(Files["RvD WT_down"], "AAACCCTAA|TTAGGGTTT")
temp1 <- subset(output, counts>0 & (At_homolog %in% chloroplastGenes))
# phyB up
output <- rawCount(Files["RvD phyB_up"], "AAACCCTAA|TTAGGGTTT")
temp1 <- subset(output, counts>0 & (At_homolog %in% chloroplastGenes))


# Evening Element in phyB.R vs WT.R
resTable <- makeResTable(Files[c("phyB.P_up", "phyB.P_down", "phyB.R_up", "phyB.R_down")],
                         chloroplastGenes, "AAAATATCT|AGATATTTT")
plotBars(resTable, "Evening Element AAAATATCT")

# Evening Element in phyB.R vs WT.R
resTable <- makeResTable(Files[c("phyB.P_up", "phyB.P_down", "phyB.R_up", "phyB.R_down")],
                         chloroplastGenes, "(A|T)AAATATCT|AGATATTT(T|A)")
plotBars(resTable, "full evening AGATATTT(T|A)")

output <- rawCount(Files["phyB.R_down"], "AAATATCT|AGATATTT")
temp1 <- subset(output, counts>0 & (At_homolog %in% chloroplastGenes))

# temp3 <- temp2[!(temp2$gene_ID %in% temp1$gene_ID), ]

clipr::write_clip(temp3)




################################################################################
# GENE PLOTTING

folder <- file.path(rootDir, "results", "temp")
for (i in 1:nrow(temp1)){
  row <- temp1[i, ]
  stylePlotDE(row$gene_ID, cpm(DGEdata), At_hom=row$short_name, save=TRUE, saveDir=folder)
}



