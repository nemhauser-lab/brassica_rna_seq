
# note: much of this code is taken directly from the Law 2018 supplemental methods

library(rprojroot)
library(pd.aragene.1.1.st)
library(oligo)
library(limma)
library(scatterplot3d)
library(vcd)
library(vsn)
library(DESeq)
library(LSD)
library(RColorBrewer)
library(VennDiagram)
cols <- brewer.pal(8,"Dark2")

rootDir <- find_root(is_rstudio_project)
law_data_path <- file.path(rootDir, "data", "Law_et_al_2018", "CEL_files")
celfiles <- list.celfiles(path=law_data_path, full.names=TRUE)
rawdata<-read.celfiles(files=celfiles)

pms <- pm(rawdata)
pms<-vsn::vsnMatrix(pms)
meanSdPlot(pms)
pm(rawdata) <- exprs(pms)
eset <- rma(rawdata,normalize=FALSE,background=FALSE)

# annotation
annot <- pData(getNetAffx(eset,'transcript'))
stopifnot(all(rownames(exprs(eset)) == rownames(annot)))
sel <- rowSums(is.na(annot))==18
annot <- annot[!sel,]
atg <- regexpr("at\\d{1}g\\d{5}",annot[,"mrnaassignment"],ignore.case=TRUE)
ID <- substr(annot[,"mrnaassignment"],atg,atg+attr(atg,"match.length")-1)
annot$ID <- ID

# "between array normalization"
n.mat <- 2^exprs(eset)[!sel,]
rownames(n.mat)[-which(is.na(ID) | ID =="")] <- ID[-which(is.na(ID) | ID =="")]
colnames(n.mat) <- sub("_Ara.*","",colnames(n.mat))
mat <- normalizeBetweenArrays(n.mat,method="cyclicloess")

# remove duplicated genes
sprintf("There are %s duplicated genes (i.e. 2 or more probes per gene)",sum(duplicated(rownames(mat))))
sel <- rownames(mat) %in% rownames(mat)[duplicated(rownames(mat))]
print("Number of probes per gene for these genes:")
table(table(rownames(mat)[sel]))

mat <- rbind(mat[!sel,],
             t(sapply(unique(rownames(mat)[duplicated(rownames(mat))]),
                      function(nam,mat){
                        lmat <- mat[rownames(mat) %in% nam,]
                        lmat[ rowMeans(lmat) == max(rowMeans(lmat)), ,drop=FALSE][1,]
                      },mat)))
mat <- mat[grepl("A[t,T]",rownames(mat)),]
mat <- mat[order(rownames(mat)),]

# quality control (just to see if I have the same data they do)
par(mar=c(6.1,4.1,4.1,2.1))
boxplot(mat,las=2,main="normalized")

# DEG analysis
# design
sample.type <- sub("Wt","IDL",sub("Wt-t0","t0",gsub("_Rep\\d","",colnames(mat))))
sample.type <- factor(sub("^.*_", "", sample.type)) # removes GSE..._
levels(sample.type) <- sub("-",".",levels(sample.type))
design <- model.matrix(~0+sample.type)
colnames(design) <- levels(sample.type)
fit <- lmFit(mat,design)

# model fit DP.t3 vs. IDL.t3 (should get 3819 DEGs)
count.matrix <- makeContrasts(T3=DP.t3-IDL.t3,levels=design)
fit2 <- contrasts.fit(fit,count.matrix)
fit2 <- eBayes(fit2)
res <- topTable(fit2,number=nrow(mat))
sum(res$adj.P.Val <= 0.01) # should get 3819 DEGs

#===============================================================================

# WT D vs P
DvP_DF <- read.table(file.path(rootDir, "results/local", "WT.D_vs_WT.P_all_genes_20190516.csv"),
                     row.names=1, header=TRUE, sep=",", stringsAsFactors=FALSE)

# load arabidopsis blast results for Brassica genes
AtBlastPath <- file.path(rootDir, "data", "blast_top_hits.xlsx")
df <- read.xlsx(xlsxFile=AtBlastPath, sheet=2, skipEmptyRows=FALSE)
colnames(df) <- gsub("(\\.)([a-z])", "_\\2", colnames(df))
colnames(df)[1] <- "Gene_ID"
AtBlastGeneInfo <- df

# DvP_DF <- allGenesWT_DvP
# DvP_DF <- WT_DvP_results$allGenes
DvP_DF$Gene_ID <- row.names(DvP_DF)
DvP_DF <- dplyr::left_join(DvP_DF, AtBlastGeneInfo, by="Gene_ID")

# reset design as Law's desing after running DGE_analysis...
# design <- model.matrix(~0+sample.type)
# colnames(design) <- levels(sample.type)

#===============================================================================
#Comparisons

groupKey <- data.frame(DvP = c(FALSE, FALSE, FALSE, FALSE, TRUE,  TRUE,  TRUE,  TRUE ),
                       DP =  c(FALSE, FALSE, TRUE,  TRUE,  FALSE, FALSE, TRUE,  TRUE),
                       IDL = c(FALSE, TRUE,  FALSE, TRUE,  FALSE, TRUE,  FALSE, TRUE),
                       group=c("none", "IDL", "DP", "IDL&DP", "DvP", "DvP&IDL", "DvP&DP", "DvP&DP&IDL"))



# ---------------------------------
# D3
count.matrix <- makeContrasts(DP_D3=DP.t3-t0,levels=design)
fit2 <- contrasts.fit(fit,count.matrix)
fit2 <- eBayes(fit2)
DP_res <- topTable(fit2,number=nrow(mat))
sum(DP_res$adj.P.Val <= 0.01) # should be 2826 + 3968 + 1001 + 899 = 8694 DEGs
DP_res$At_Gene_ID <- row.names(DP_res)

count.matrix <- makeContrasts(IDL_D3=IDL.t3-t0,levels=design)
fit2 <- contrasts.fit(fit,count.matrix)
fit2 <- eBayes(fit2)
IDL_res <- topTable(fit2,number=nrow(mat))
sum(IDL_res$adj.P.Val <= 0.01) # should be 2826 + 3968 + 1793 + 1104 = 9691 DEGs
IDL_res$At_Gene_ID <- row.names(IDL_res)

res <- dplyr::inner_join(DP_res, IDL_res, by="At_Gene_ID")
colnames(res) <- sub("\\.x", ".DP", colnames(res))
colnames(res) <- sub("\\.y", ".IDL", colnames(res))

# recreate T3  venn diagram from Law fig. 2B
resVennDF <- data.frame("DP.t3" = with(res, (adj.P.Val.DP <= 0.01)*sign(logFC.DP)),
                        "IDL.t3" = with(res, (adj.P.Val.IDL <= 0.01)*sign(logFC.IDL)),
                        row.names = res$At_Gene_ID)
vennDiagram(resVennDF, include=c("up","down"))
title(list("Recreated Law 2018 fig 2B", font=1, cex=2),  line=-1.5)
# the fact that the gene counts in our diagram match the Law 2018 figure exactly
# gives us confidence that we are re-creating their gene lists correctly.

# make a dataframe combining our WT DvP results with thier T3 results
commonDF <- dplyr::inner_join(res, DvP_DF[,c(1:8)], by=c("At_Gene_ID" = "A._thaliana_best_hit"))

# Venn-diagrams:
vennDF <- data.frame("DvP" = with(commonDF, (abs(logFC) >= 1 & FDR <= 0.01) * sign(logFC)),
                     "DP.t3" = with(commonDF, (adj.P.Val.DP <= 0.01) * sign(logFC.DP)),
                     "IDL.t3" = with(commonDF, (adj.P.Val.IDL <= 0.01) * sign(logFC.IDL)))
vennDiagram(vennDF, include=c("up","down"))
title(list("Non-unique A.t homologs of Brassica genes", font=1, cex=2),  line=-1.5)
# vennDiagram(vennDF, include="both")

# Add venn categories to commonDF and save
commonDF$t3.upVennGroup <- dplyr::left_join(data.frame(vennDF==1), groupKey,
                                            by=c("DvP"="DvP", "DP.t3"="DP",
                                                 "IDL.t3"="IDL"))$group
commonDF$t3.downVennGroup <- dplyr::left_join(data.frame(vennDF==-1), groupKey,
                                              by=c("DvP"="DvP", "DP.t3"="DP",
                                                   "IDL.t3"="IDL"))$group

# generate unique A.t Venn diagram lists and counts
t3.upUnique <- unique(commonDF[, c("At_Gene_ID", "t3.upVennGroup")])
table(t3.upUnique$t3.upVennGroup)

t3.downUnique <- unique(commonDF[, c("At_Gene_ID", "t3.downVennGroup")])
table(t3.downUnique$t3.downVennGroup)

# save files
fname <-

# ---------------------------------------
# D6

count.matrix <- makeContrasts(DP_D6=DP.t6-t0,levels=design)
fit2 <- contrasts.fit(fit,count.matrix)
fit2 <- eBayes(fit2)
DP_res <- topTable(fit2,number=nrow(mat))
sum(DP_res$adj.P.Val <= 0.01) # should be 2714 + 3963 + 1708 + 1422 = 9807 DEGs
DP_res$At_Gene_ID <- row.names(DP_res)

count.matrix <- makeContrasts(IDL_D6=IDL.t6-t0,levels=design)
fit2 <- contrasts.fit(fit,count.matrix)
fit2 <- eBayes(fit2)
IDL_res <- topTable(fit2,number=nrow(mat))
sum(IDL_res$adj.P.Val <= 0.01) # should be 2714 + 3963 + 2313 + 1270 = 10260 DEGs
IDL_res$At_Gene_ID <- row.names(IDL_res)

res <- dplyr::inner_join(DP_res, IDL_res, by="At_Gene_ID")
colnames(res) <- sub("\\.x", ".DP", colnames(res))
colnames(res) <- sub("\\.y", ".IDL", colnames(res))

# recreate T6 venn diagram from Law fig. 2B
resVennDF <- data.frame("DP.t6" = with(res, (adj.P.Val.DP <= 0.01)*sign(logFC.DP)),
                        "IDL.t6" = with(res, (adj.P.Val.IDL <= 0.01)*sign(logFC.IDL)),
                        row.names = res$At_Gene_ID)
colnames
vennDiagram(resVennDF, include=c("up","down"))

commonDF <- dplyr::inner_join(res, DvP_DF[,c(1:8)], by=c("At_Gene_ID" = "A._thaliana_best_hit"))


vennDF <- data.frame("DvP" = with(commonDF, (abs(logFC)>1 & FDR<0.01)*sign(logFC)),
                     "DP.t6" = with(commonDF, (adj.P.Val.DP <= 0.01)*sign(logFC.DP)),
                     "IDL.t6" = with(commonDF, (adj.P.Val.IDL <= 0.01)*sign(logFC.IDL)))
vennDiagram(vennDF, include=c("up","down"))
title(list("Non-unique A.t homologs of Brassica genes", font=1, cex=2),  line=-1.5)
# vennDiagram(vennDF, include="both")






################################################################################
# MISC.

temp <- aaply(t(vennDF), .margins=2,  paste, collapse="", .expand=FALSE)

temp <- vennDF$DvP!=-1 & vennDF$DP.t3!=-1 & vennDF$IDL.t3!=-1
sum(temp)
length(unique(commonDF$At_Gene_ID[temp]))

# coppy list of unique At_Gene_IDs to clipboard for specific condition
temp <- vennDF$DvP!=-1 & vennDF$DP.t6==-1 & vennDF$IDL.t6==-1
sum(temp)
clipr::write_clip(unique(commonDF$At_Gene_ID[temp]))

clipr::write_clip(commonDF)

groupKey <- data.frame(DvP = c(FALSE, FALSE, FALSE, FALSE, TRUE,  TRUE,  TRUE,  TRUE ),
                       DP =  c(FALSE, FALSE, TRUE,  TRUE,  FALSE, FALSE, TRUE,  TRUE),
                       IDL = c(FALSE, TRUE,  FALSE, TRUE,  FALSE, TRUE,  FALSE, TRUE),
                       group=c("none", "IDL", "DP", "IDL&DP", "DvP", "DvP&IDL", "DvP&DP", "DvP&DP&IDL"))

t6.upGroup <- dplyr::left_join(data.frame(vennDF==1), groupKey, by=c("DvP"="DvP", "DP.t6"="DP", "IDL.t6"="IDL"))$group
t6.downGroup <- dplyr::left_join(data.frame(vennDF==-1), groupKey, by=c("DvP"="DvP", "DP.t6"="DP", "IDL.t6"="IDL"))$group
temp <- cbind(vennDF, t6.upGroup, t6.downGroup)

