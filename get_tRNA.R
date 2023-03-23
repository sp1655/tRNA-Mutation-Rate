# Getting the tRNA data out of the gff3 file and writing it to its own file

rm(list = ls())

##library(GenomicRanges)
library(genomation)
library(dplyr)

setwd("C:/Users/Sarah/Desktop/Lab/RNAPIII/RNAPIII-mutation-rate-main/Data/Scerevisiae/")

# Isolating out the tRNA data
gff.list = gffToGRanges("DB/Saccharomyces_cerevisiae.R64-1-1.84.gff3")
tRNA.list = gff.list[which(gff.list$biotype == "tRNA" & gff.list$type == "tRNA_gene")]

# Converting to a dataframe
tRNA.tab = as.data.frame(tRNA.list)
tRNA.tab = tRNA.tab[,c("seqnames", "biotype", "start", "end", "strand")]

# Renaming columns - if this isn't working, restart R
tRNA.tab = rename(tRNA.tab, geneType = biotype)
tRNA.tab = rename(tRNA.tab, chromosome = seqnames)

# Checking to make sure we only got the tRNA data
unique(tRNA.tab$geneType)

# Changing the geneType data to match the other files
tRNA.tab$geneType = "tRNA_gene"

# Writing the table
write.table(tRNA.tab, file = "DB/RNAPIII-pos.tab", row.names = F, col.names = T, sep = '\t', quote = F)
