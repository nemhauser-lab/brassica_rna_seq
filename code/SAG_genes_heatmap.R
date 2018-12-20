library(rprojroot)

rootDir <- find_root(is_rstudio_project)

# this assumes you ran the first part of DGE_analysis_v3.Rmd or equiv.
DGEdata.no24 <- DGEdata[,-c(5:8,20:22)]

tairGenesFile <- file.path(rootDir, "data", "TAIR_SAG_Genes.csv")
tairGenes <- read.csv(tairGenesFile, stringsAsFactors=FALSE)


sagGenes <- subset(allGenesEdin, A._thaliana_best_hit %in% tairGenes$Gene.ID)

sagGenes <- dplyr::left_join(tairGenes, allGenesEdin, by=c("Gene.ID"="A._thaliana_best_hit"))
sagGenes <- subset(sagGenes, !is.na(Gene_ID))

temp <- geneHeatmap(sagGenes$Gene_ID, DGEdata.no24, groupAv=TRUE, rsOnly=FALSE,
                    labels=sagGenes$Symbol.y)

SAGs <- subset(sagGenes, Group=="SAG")
#SAGs <- subset(allGenesEdin, A._thaliana_best_hit %in% SAG_IDs$Gene.ID)
temp <- geneHeatmap(SAGs$Gene_ID, DGEdata.no24, groupAv=TRUE, rsOnly=FALSE, labels=SAGs$SAG.symbol)


SAG_IDs <- subset(tairGenes, Group %in% c( "SA_", "SRG"))
SAGs <- subset(allGenesEdin, A._thaliana_best_hit %in% SAG_IDs$Gene.ID)
temp <- geneHeatmap(SAGs$Gene_ID, DGEdata.no24, groupAv=TRUE, rsOnly=FALSE, labels=SAGs$Symbol)




