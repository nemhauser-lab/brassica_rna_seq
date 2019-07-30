# calculate stats on the presense of G-Box 6-mers (CACGTG) in promoter sequences

library(Biostrings)
library(GenomicRanges)
library(rprojroot)
library(plyr)
library(stringr)
library(tibble)

rootDir <- find_root(is_rstudio_project)
resultsFoldeToRead <- file.path(rootDir, "results", "analysis_output_20190510")

# R v D_______________________________________________________________________
RvD_files <- c(
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
                            "RvD_BOTH_SIG_DOWN_promoters_FDR01_LFC01_20190510.fasta")
)

# D v P_______________________________________________________________________
DvP_files <- c(
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
                              "DvP_BOTH_SIG_DOWN_promoters_FDR01_LFC01_20190510.fasta")
)

# R v P_______________________________________________________________________
RvP_files <- c(
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
                              "RvP_BOTH_SIG_DOWN_promoters_FDR01_LFC01_20190510.fasta")
)

# WT vs phyB at P_____________________________________________________________
Pre_files <- c(
  "phyB.P_up" = file.path(resultsFoldeToRead,
                          "phyB_v_WT_at_timepoint/promoters",
                          "phyB.P_vs_WT.P_exact_DOWN_promoters_FDR05_LFC01_20190510.fasta"),
  "phyB.P_down" = file.path(resultsFoldeToRead,
                            "phyB_v_WT_at_timepoint/promoters",
                            "phyB.P_vs_WT.P_exact_UP_promoters_FDR05_LFC01_20190510.fasta")
)

# WT vs phyB at R_____________________________________________________________
Rec_files <- c(
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
# should be 5248 chloroplast genes total

# ==============================================================================
# functions

makeResTable <- function(Files, chloroplastGenes, motif){
  resTable <- data.frame()
  # chloroplast genes
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
  resTable[resTable$sign=="down", "percent_w_motif"] <- -resTable[resTable$sign=="down", "percent_w_motif"]
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
# Search for motifs

motifs <- tribble(~name, ~pattern,
                  "G-box", "CACGTG",
                  "Evening_element", "AAAATATCT|AGATATTTT",
                  "Telo-box", "AAACCCTAA|TTAGGGTTT"
                  )

for (i in 1:nrow(motifs)){
  motif <- motifs[i, ]
  resTable <- makeResTable(RvD_files, chloroplastGenes, motif$pattern)
  p <- plotBars(resTable, paste0(motif$name, " motif in RvD differentially regulated genes"))
  saveFile <- file.path(rootDir, "results/local/motif_enrichment", paste0(motif$name, "_RvD_barplot.svg"))
  ggsave(saveFile, p, device=svg)

  resTable <- makeResTable(DvP_files, chloroplastGenes, motif$pattern)
  p <- plotBars(resTable, paste0(motif$name, " motif in DvP differentially regulated genes"))
  saveFile <- file.path(rootDir, "results/local/motif_enrichment", paste0(motif$name, "_DvP_barplot.svg"))
  ggsave(saveFile, p, device=svg)
}

  saveFile <- file.path(rootDir, "results/local/motif_enrichment", paste0(motif$name, "_RvD_barplot.svg"))
  ggsave(saveFile, p, device=svg)
