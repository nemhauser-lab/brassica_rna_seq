subset(allGenesEdin, A._thaliana_best_hit %in% c("AT1G66580",
                                                 "AT1G71190",
                                                 "AT2G29350",
                                                 "AT2G45210",
                                                 "AT3G10985",
                                                 "AT4G02380",
                                                 "AT4G35770",
                                                 "AT5G13170",
                                                 "AT5G14930",
                                                 "AT5G20230",
                                                 "AT5G45890",
                                                 "AT5G51070",
                                                 "AT5G59220",
                                                 "AT5G60360"))


DGEdata.no24 <- DGEdata[,-c(5:8,20:22)]


tairGenes <- read.csv("C:/Users/Morgan/Documents/Lab Stuff/RNA seq/TAIR_SAG_Genes.csv",
                      stringsAsFactors=FALSE)


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




