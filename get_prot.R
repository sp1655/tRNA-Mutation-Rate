# Code for getting the protein-coding gene data from the gff3 file

rm(list = ls())

library(genomation)
library(dplyr)

setwd("C:/Users/Sarah/Desktop/Lab/RNAPIII/RNAPIII-mutation-rate-main/Data/Scerevisiae/DB/")

############################################################################################################################################
# Getting the protein-coding genes

gff.list = gffToGRanges("Saccharomyces_cerevisiae.R64-1-1.84.gff3")
prot.list = gff.list[which(gff.list$biotype == "protein_coding" & gff.list$type == "gene")]

# Converting to dataframe
prot.tab = as.data.frame(prot.list)
prot.tab = prot.tab[,c("seqnames", "biotype", "start", "end", "strand", "gene_id")]

# Renaming columns - if this isn't working, restart R
prot.tab = rename(prot.tab, geneType = biotype)
prot.tab = rename(prot.tab, chromosome = seqnames)

############################################################################################################################################
# Separating by expression level

low.tab = read.table("lowProt.tab", h = T)
high.tab = read.table("highProt.tab", h = T)
med.tab = read.table("medProt.tab", h = T)

lowProt.tab = prot.tab[which(prot.tab$gene_id %in% low.tab$Sistematic.name),]
highProt.tab = prot.tab[which(prot.tab$gene_id %in% high.tab$Sistematic.name),]
medProt.tab = prot.tab[which(prot.tab$gene_id %in% med.tab$Sistematic.name),]

# Removing the gene names
prot.tab = prot.tab[,-6]
lowProt.tab = lowProt.tab[,-6]
highProt.tab = highProt.tab[,-6]
medProt.tab = medProt.tab[,-6]

############################################################################################################################################
# Writing the tables
write.table(prot.tab, file = "prot-pos.tab", row.names = F, col.names = T, sep = '\t', quote = F)
write.table(lowProt.tab, file = "low-prot-pos.tab", row.names = F, col.names = T, sep = '\t', quote = F)
write.table(highProt.tab, file = "high-prot-pos.tab", row.names = F, col.names = T, sep = '\t', quote = F)
write.table(medProt.tab, file = "med-prot-pos.tab", row.names = F, col.names = T, sep = '\t', quote = F)




############################################################################################################################################
# Fixing the E. coli protein coding to have strand data

# MG1655 (Foster)

setwd("C:/Users/Sarah/Desktop/Lab/RNAPIII/RNAPIII-mutation-rate-main/Data/Ecoli/DB/")

gff.list = gffToGRanges("U00096.3.gff3")
gff.tab = as.data.frame(gff.list)

low.tab = read.table("mg1655-low-prot-pos.tab", h = T)
med.tab = read.table("mg1655-med-prot-pos.tab", h = T)
high.tab = read.table("mg1655-high-prot-pos.tab", h = T)

# Fixing low expression table
for(i in 1:nrow(low.tab)) {
  start = low.tab$start[i]
  end = low.tab$end[i]
  
  str = as.character(gff.tab$strand[which(gff.tab$start == start & gff.tab$end == end)])
  
  if(length(str == 1)) {
    low.tab$strand[i] = str
  }
  
  else {
    str = as.character(str[1])
    low.tab$strand[i] = str[1]
  }
}

# Fixing medium expression table
for(i in 1:nrow(med.tab)) {
  start = med.tab$start[i]
  end = med.tab$end[i]
  
  str = as.character(gff.tab$strand[which(gff.tab$start == start & gff.tab$end == end)])
  
  if(length(str == 1)) {
    med.tab$strand[i] = str
  }
  
  else {
    str = as.character(str[1])
    med.tab$strand[i] = str[1]
  }
}

# Fixing high expression table
for(i in 1:nrow(high.tab)) {
  start = high.tab$start[i]
  end = high.tab$end[i]
  
  str = as.character(gff.tab$strand[which(gff.tab$start == start & gff.tab$end == end)])
  
  if(length(str == 1)) {
    high.tab$strand[i] = str
  }
  
  else {
    str = as.character(str[1])
    high.tab$strand[i] = str[1]
  }
}

# Writing the fixed tables over the old
write.table(low.tab, file = "mg1655-low-prot-pos.tab", quote = F, sep = '\t', row.names = F, col.names = T)
write.table(med.tab, file = "mg1655-med-prot-pos.tab", quote = F, sep = '\t', row.names = F, col.names = T)
write.table(high.tab, file = "mg1655-high-prot-pos.tab", quote = F, sep = '\t', row.names = F, col.names = T)

########################################################################################
# ATCC 8739 (Zhang)

gff.list = gffToGRanges("CP000946.1.gff3")
gff.tab = as.data.frame(gff.list)

low.tab = read.table("atcc-low-prot-pos.tab", h = T)
med.tab = read.table("atcc-med-prot-pos.tab", h = T)
high.tab = read.table("atcc-high-prot-pos.tab", h = T)

# Fixing low expression table
for(i in 1:nrow(low.tab)) {
  start = low.tab$start[i]
  end = low.tab$end[i]
  
  str = as.character(gff.tab$strand[which(gff.tab$start == start & gff.tab$end == end)])
  
  if(length(str == 1)) {
    low.tab$strand[i] = str
  }
  
  else {
    str = as.character(str[1])
    low.tab$strand[i] = str[1]
  }
}

# Fixing medium expression table
for(i in 1:nrow(med.tab)) {
  start = med.tab$start[i]
  end = med.tab$end[i]
  
  str = as.character(gff.tab$strand[which(gff.tab$start == start & gff.tab$end == end)])
  
  if(length(str == 1)) {
    med.tab$strand[i] = str
  }
  
  else {
    str = as.character(str[1])
    med.tab$strand[i] = str[1]
  }
}

# Fixing high expression table
for(i in 1:nrow(high.tab)) {
  start = high.tab$start[i]
  end = high.tab$end[i]
  
  str = as.character(gff.tab$strand[which(gff.tab$start == start & gff.tab$end == end)])
  
  if(length(str == 1)) {
    high.tab$strand[i] = str
  }
  
  else {
    str = as.character(str[1])
    high.tab$strand[i] = str[1]
  }
}

# Writing the fixed tables over the old
write.table(low.tab, file = "atcc-low-prot-pos.tab", quote = F, sep = '\t', row.names = F, col.names = T)
write.table(med.tab, file = "atcc-med-prot-pos.tab", quote = F, sep = '\t', row.names = F, col.names = T)
write.table(high.tab, file = "atcc-high-prot-pos.tab", quote = F, sep = '\t', row.names = F, col.names = T)
