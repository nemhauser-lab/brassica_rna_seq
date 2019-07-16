library(Biostrings)
library(GenomicRanges)
library(rprojroot)
library(plyr)
library(BSgenome)

rootDir <- find_root(is_rstudio_project)
dirPath <- file.path(rootDir, "data", "Brapa2.5")

genomeFile <- file.path(dirPath, "BrapaV2.5_Chr.fa")
BrapaGenome <- readDNAStringSet(genomeFile)
gffFile <- file.path(dirPath, "BrapaV2.5_Chr.gene.gff")
gffData <- ape::read.gff(gffFile)

# parse the attributes field
gffData$ID <- stringr::str_match(gffData$attributes,
                                       "ID=(.*?)(;.*)?$")[,2]
gffData$Parent <- stringr::str_match(gffData$attributes,
                                 "Parent=(.*?)(;.*)?$")[,2]
gffData$Name <- stringr::str_match(gffData$attributes,
                                     "Name=(.*?)(;.*)?$")[,2]
gffData$Gene_ID <- stringr::str_match(gffData$ID,
                                      "BraA\\d{8}|BraSca\\d{6}")

# process entries of the "gene" type
geneData <- gffData[gffData$type %in% "gene", ]
# geneData$gene_id <- stringr::str_match(geneData$attributes,
#                                        "ID=evm\\.TU\\.(.*?);")[,2]

geneIDs <- read.table(file.path(dirPath, "all_gene_ids.csv"), stringsAsFactors=FALSE)
# ourGeneData <- subset(geneData, Gene_ID %in% geneIDs$V1)
ourGeneData <- subset(geneData, !is.na(Gene_ID))

# gene sequences
ourGenesGR <- GRanges(seqnames=ourGeneData$seqid,
                      ranges=IRanges(start=ourGeneData$start, end=ourGeneData$end,
                                     names=ourGeneData$Gene_ID),
                      strand=ourGeneData$strand)

ourGenesDSS <- getSeq(BrapaGenome, ourGenesGR)


writeXStringSet(ourGenesDSS, file.path(dirPath, "all_genes.fasta"))


# Promoters
promoterLength <- 1000
df <- adply(.data=ourGeneData,
            .margins=1,
            .fun=function(row){
              strand <- as.character(row$strand)
              if (strand=="+"){
                promoterStart <- max(row$start - promoterLength, 1)
                promoterEnd <- max(row$start - 1, 1)
              } else if(strand=="-"){
                promoterStart <- min(row$end + 1, width(BrapaGenome[row$seqid]))
                promoterEnd <- min(row$end + promoterLength, width(BrapaGenome[row$seqid]))
              }
              output <- cbind(row, data.frame(promoterStart, promoterEnd))
            })

ourPromoterGR <- GRanges(seqnames=df$seqid,
                      ranges=IRanges(start=df$promoterStart, end=df$promoterEnd,
                                     names=df$Gene_ID),
                      strand=df$strand)

ourPromoterDSS <- getSeq(BrapaGenome, ourPromoterGR)
fname <- paste("all_genes_promoters_",
               as.character(round(promoterLength/1000, 0)),
               "k.fasta", sep="")
writeXStringSet(ourPromoterDSS, file.path(dirPath, fname))

# ==============================================================================
# Misc. stuff

# to read in a fasta:
dirPath <- file.path(rootDir, "data", "Brapa2.5")
Promoters <- readDNAStringSet(file.path(dirPath, "all_genes_promoters_1k.fasta"))

# temp df to test limits of string set indexing.
df <- data.frame("seqid"=rep("A01", 3),
                 "strand"=rep("+", 3),
                 "promoterStart"=c(1,1,33885991),
                 "promoterEnd"=c(1,1,33885992))


