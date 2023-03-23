# Makes the tpmLow.tab and tpmHigh.tab files in the S. cerevisiae folder

rm(list = ls())

library(genomation)

setwd("C:/Users/Sarah/Desktop/Lab/RNAPIII/RNAPIII-mutation-rate-main/Data/Scerevisiae/DB/")

# Reading in the table and gff3 file
readVals.tab = read.table("readVals.tab", h = T)

gff.list = gffToGRanges("Saccharomyces_cerevisiae.R64-1-1.84.gff3")

# Removing any genes not found in the gff3 file
nrow(readVals.tab[which(readVals.tab$gene %in% gff.list$gene_id),])
readVals.tab = readVals.tab[which(readVals.tab$gene %in% gff.list$gene_id),]

# Isolating the tRNA genes from the gff3 file that are also in readVal.tab
tRNA.list = gff.list[which(gff.list$gene_id %in% readVals.tab$gene)]
tRNA.tab = as.data.frame(tRNA.list)

tCol = c("seqnames", "start", "end", "strand", "gene_id", "biotype")
tRNA.tab = tRNA.tab[,tCol]

colnames(tRNA.tab) = c("chromosome", "start", "end", "strand", "gene", "geneType")

# Merging tRNA position and read number
tRNAvals.tab = merge(readVals.tab, tRNA.tab)
tCol = c("chromosome", "geneType", "start", "end", "strand", "numReads")
tRNAvals.tab = tRNAvals.tab[,tCol]

hist(tRNAvals.tab$numReads)

# Getting high and low ends of read number
split = median(tRNAvals.tab$numReads)

low.tab = tRNAvals.tab[which(tRNAvals.tab$numReads < split),]
high.tab = tRNAvals.tab[which(tRNAvals.tab$numReads >= split),]

low.tab = low.tab[,-6]
high.tab = high.tab[,-6]

# Writing the tables
write.table(low.tab, file = "tpmLow.tab", quote = F, sep = '\t', row.names = F, col.names = T)
write.table(high.tab, file = "tpmHigh.tab", quote = F, sep = '\t', row.names = F, col.names = T)

############################################################################################################################################
# Graph
# sz = 1
# 
# tpmPlot <- ggplot(tRNAvals.tab, aes(x = numReads)) + 
#   geom_histogram(color = "#2E4052", fill = "#BDD9BF", alpha = 0.8, bins = 12) + # Making the density plot
#   geom_vline(aes(xintercept = split), color = "#A997DF", linetype = "solid", size = sz) + # Making the cutoff line
#   geom_text(aes(minVal, 30, label = "TPM Split Point", angle= 90, vjust = -0.5)) + # Labeling the cutoff line
#   labs(title = "S. cerevisiae tRNA TPM Values", x = "TPM Value", y = "Count") + # Labeling the graph
#   theme_classic()
# tpmPlot
