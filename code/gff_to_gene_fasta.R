library(Biostrings)
library(GenomicRanges)

dirPath <- "C:/Users/Morgan/Desktop/B_Rapa2.5"


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
                                      "BraA\\d{8}")

# process entries of the "gene" type
geneData <- gffData[gffData$type %in% "gene", ]
# geneData$gene_id <- stringr::str_match(geneData$attributes,
#                                        "ID=evm\\.TU\\.(.*?);")[,2]

geneIDs <- read.table(file.path(dirPath, "all_gene_ids.csv"), stringsAsFactors=FALSE)
ourGeneData <- subset(geneData, Gene_ID %in% geneIDs$V1)


ourGenesGR <- GRanges(seqnames=ourGeneData$seqid,
                      ranges=IRanges(start=ourGeneData$start, end=ourGeneData$end,
                                     names=ourGeneData$Gene_ID),
                      strand=ourGeneData$strand)

ourGenesDSS <- BrapaGenome[ourGenesGR]
names(ourGenesDSS) <- names(ourGenesGR)

writeXStringSet(ourGenesDSS, file.path(dirPath, "all_genes.fasta"))

# --------------------------------------------------------------

# try to make coding sequences from CDS entries



