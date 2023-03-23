# Finds high & low expression genes by tpm value in dataset

rm(list = ls())

library(readxl)
library(genomation)
library(GenomicRanges)
library(BSgenome)
##library(tidyr)

setwd("C:/Users/Sarah/Desktop/Lab/RNAPIII/RNAPIII-mutation-rate-main/Data/Ecoli/")

############################################################################################################################################
# Finding the high and low expression tRNA genes

tpmData = read_xlsx("C:/Users/Sarah/Desktop/Lab/RNAPIII/RNAPIII-mutation-rate-main/Data/Ecoli/41598_2019_39369_MOESM2_ESM.xlsx", 
                    sheet = "Ecoli174")
tpmData = as.data.frame(tpmData)

tpmData = tpmData[,c("target_id", "tpm")]
tpmData = tpmData[which(is.na(tpmData$target_id) == F),]

target_id = tpmData$target_id

# Getting out just the lines with only one tRNA
onetRNA <- c()
for(i in 1:length(target_id)) {
  id = target_id[i]
  tmp = unlist(strsplit(id, split = "[.]"))
  
  if(length(tmp) == 2) {
    onetRNA <- append(onetRNA, id)
  }
}

# Looking at just the lines with only one tRNA
onetRNA.tab = tpmData[which(tpmData$target_id %in% onetRNA),]
names <- c()

# Getting the shortened names
for(i in 1:length(onetRNA)) {
  id = onetRNA[i]
  tmp = unlist(strsplit(id, split = "[_]"))[[4]]
  tmp = unlist(strsplit(tmp, split = "[-]"))[[1]]
  
  names <- append(names, tmp)
}

# Putting the names into the onetRNA.tab table
onetRNA.tab$name = names
onetRNA.tab = onetRNA.tab[,c(1,3,2)]

# Getting the start and end positions from the comp.tab table
##for(i in 1:nrow(onetRNA.tab)) {
##  id = onetRNA.tab$name[i]
  
##  if(id %in% comp.tab$name) {
##    startPos = comp.tab$gffStart[which(comp.tab$name == id)]
##    endPos = comp.tab$gffEnd[which(comp.tab$name == id)]
    
##    onetRNA.tab$start[i] = startPos
##    onetRNA.tab$end[i] = endPos
##  }
  
##  else {
##    onetRNA.tab$start[i] = NA
##    onetRNA.tab$end[i] = NA
##  }
##}

# Seeing if any of the tRNAs were not in the tmptRNAsite table
##onetRNA.tab[which(is.na(onetRNA.tab$start) == T),]

# Removing the tRNAs that were not in the tmptRNAsite table
##onetRNA.tab = onetRNA.tab[which(is.na(onetRNA.tab$start) == F),]

# Getting the names for the site data into onetRNA.tab
correct.tab = read.table("DB/correct-tRNA-pos.tab", h = T)
coords.tab = read.table("DB/correct-coords.tab", h = T)
onetRNA.tab = onetRNA.tab[which(onetRNA.tab$name %in% correct.tab$oldName),]

# Replacing the old name with the new names
for(i in 1:nrow(onetRNA.tab)) {
  id = onetRNA.tab$name[i]
  
  # Replacing the old id with the new one
  new_id = correct.tab$name[which(correct.tab$oldName == id)]
  
  # Adding in the v3 and v2 coordinates
  startv3 = coords.tab$startv3[which(coords.tab$name == new_id)]
  endv3 = coords.tab$endv3[which(coords.tab$name == new_id)]
  
  startv2 = coords.tab$startv2[which(coords.tab$name == new_id)]
  endv2 = coords.tab$endv2[which(coords.tab$name == new_id)]
  
  onetRNA.tab$name[i] = new_id
  onetRNA.tab$startv3[i] = startv3
  onetRNA.tab$endv3[i] = endv3
  onetRNA.tab$startv2[i] = startv2
  onetRNA.tab$endv2[i] = endv2
}

tRNAmatched  = onetRNA.tab$name

#########################################################################################################################
# Looking at tpm distribution to determine cutoffs for high/low expression

tpmRes = onetRNA.tab$tpm

##tpmDensity = density(tpmRes) ## This makes up a negative minimum value?? Only 2 below min value
##plot(tpmDensity)

hist(tpmRes)

midLow = onetRNA.tab[which(onetRNA.tab$tpm < 15000 & onetRNA.tab$tpm > 10000),]
midHigh = onetRNA.tab[which(onetRNA.tab$tpm < 20000 & onetRNA.tab$tpm > 15000),]

max(midLow$tpm)
min(midHigh$tpm)

min(midLow$tpm)
max(midHigh$tpm)

minVal = min(midLow$tpm)
maxVal = max(midHigh$tpm)

#########################################################################################################################
# Getting position tables for comparison

min.tab = onetRNA.tab[which(onetRNA.tab$tpm <= minVal),]
max.tab = onetRNA.tab[which(onetRNA.tab$tpm >= maxVal),]

# Getting the names of tRNAs in the low and high tables
lowNames = min.tab$name
highNames = max.tab$name

# Keeping only the startv2 and endv2 columns so they can be used as tRNA-pos files
min.tab = min.tab[,c(6,7)]
max.tab = max.tab[,c(6,7)]

colnames(min.tab) = c("start", "end")
colnames(max.tab) = c("start", "end")

# Adding and reordering columns
min.tab$chromosome = "U00096.2"
min.tab$geneType = "tRNA_gene"
min.tab$strand = "+"

max.tab$chromosome = "U00096.2"
max.tab$geneType = "tRNA_gene"
max.tab$strand = "+"

min.tab = min.tab[,c(3,4,1,2,5)]
max.tab = max.tab[,c(3,4,1,2,5)]

#########################################################################################################################
# Writing the tables

write.table(onetRNA.tab, file = "DB/onetRNA.tab", quote = F, sep = '\t', row.names = F, col.names = T)

write.table(min.tab, file = "DB/tpmLowFoster.tab", quote = F, sep = '\t', row.names = F, col.names = T)
write.table(max.tab, file = "DB/tpmHighFoster.tab", quote = F, sep = '\t', row.names = F, col.names = T)

# Writing the matched-tRNAs.tab file
write.table(tRNAmatched, file = "DB/matched-tRNAs.tab", quote = F, sep = '\t', row.names = F, col.names = F)

# Writing the matched-tRNAs-low.tab and matched-tRNAs-high.tab files
write.table(lowNames, file = "DB/matched-tRNAs-low.tab", quote = F, sep = '\t', row.names = F, col.names = F)
write.table(highNames, file = "DB/matched-tRNAs-high.tab", quote = F, sep = '\t', row.names = F, col.names = F)

# Writing the MG16655-pos.tab table - need to have the v2 start/end positions
MG1655.tab = onetRNA.tab[,c(2,6,7)]
colnames(MG1655.tab) = c("name", "start", "end")

MG1655.tab$chromosome = "U00096.2"
MG1655.tab = MG1655.tab[,c(4,1,2,3)]

write.table(MG1655.tab, file = "DB/MG1655-pos.tab", quote = F, sep = '\t', row.names = F, col.names = T)

#########################################################################################################################
# Graph

# library(ggplot2)
# 
# tpm.df = as.data.frame(tpmRes)
# sz = 1
# 
# tpmPlot <- ggplot(tpm.df, aes(x = tpmRes)) + 
#   geom_histogram(color = "#2E4052", fill = "#BDD9BF", alpha = 0.8, bins = 8) + # Making the density plot
#   geom_vline(aes(xintercept = minVal), color = "#A997DF", linetype = "solid", size = sz) + # Making the cutoff lines
#   geom_vline(aes(xintercept = maxVal), color = "#A997DF", linetype = "solid", size = sz) +
#   geom_text(aes(minVal, 5.5, label = "Low TPM Cutoff", angle= 90, vjust = 1)) + # Labeling the cutoff lines
#   geom_text(aes(maxVal, 5.5, label = "High TPM Cutoff", angle= 90, vjust = 1)) +
#   labs(title = "E. coli tRNA TPM Values", x = "TPM Value", y = "Count") + # Labeling the graph
#   theme_classic()
# tpmPlot

############################################################################################################################################

# Getting list of tRNAs with just the first of each listed
tRNA_list <- c()
for(i in 1:length(target_id)) {
  id = target_id[i]
  tmp = unlist(strsplit(id, split = "_Escherichia"))
  
  tRNA_id = tmp[1]
  
  tRNA_list <- append(tRNA_list, tRNA_id)
}

# Fixing the names
correct.tab = read.table("DB/correct-tRNA-pos.tab", h = T)


write.table(tRNA_list, file = "DB/first-tRNA-list.tab", sep = "\t", quote = F, row.name = F, col.names = F)
