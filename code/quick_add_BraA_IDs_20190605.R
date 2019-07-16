
# select and copy portion of excel sheet to clipboard.
temp <- clipr::read_clip()
temp <- data.frame("temp"=temp)
temp2 <- tidyr::separate(temp, "temp", into=stringr::str_split(temp[1,], "\t")[[1]], sep="\t")
temp2 <- temp2[-1, ]

# column 14 needs a name other than ""
colnames(temp2)[14] <- "notes"

temp3 <- dplyr::left_join(temp2,
                          AtBlastGeneInfo[, c("Gene_ID", "A._thaliana_best_hit")],
                          by=c("AGI " = "A._thaliana_best_hit"))

clipr::write_clip(temp3)
