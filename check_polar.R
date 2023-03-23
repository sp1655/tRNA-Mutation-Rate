# Code to check polarization of data - is it already polarized or does it need to be?

rm(list = ls())

library(GenomicRanges)
library(Biostrings)
library(BSgenome)
library(genomation)

setwd("C:/Users/Sarah/Desktop/Lab/RNAPIII/RNAPIII-mutation-rate-main/Data/Celegans/")

fMut = "Denver2009/snm.tab"
genomeFastaFile = "DB/Celegans_genome.fasta"
gffFile = "DB/genomic.gff"

mut.tab = read.table(fMut, h = T, sep = '\t')
genome = readDNAStringSet(genomeFastaFile)
gff.list = gffToGRanges(gffFile)

############################################################################################################################################

# Fixing genome names
names = names(genome)
for(i in 1:length(names)) {
  tmp = strsplit(names[i], " ")
  tmp = unlist(tmp)
  chr = tmp[1]
  
  names[i] = chr
}

names(genome) = names

# Getting nucleotide from reference added to table
mut.tab = mut.tab[which(mut.tab$type == "BASE_SUB"),]
for(i in 1:nrow(mut.tab)) {
  chr = mut.tab$chromosome[i]
  pos = mut.tab$gPos[i]
  
  refNuc = as.character(substr(genome[chr], pos, pos))
  mut.tab$refNuc[i] = refNuc
}

# Comparing reference nucleotide to "from" category
# If this command returns 0, means results are not polarized
nrow(mut.tab[which(mut.tab$from != mut.tab$refNuc),])

############################################################################################################################################
# Assigning refNuc based on strand - run if results not polarized

mut.tab$strand = "*"
mut.tab$strand = factor(mut.tab$strand, levels = c("+", "-", "*"))

# Removing unnecessary lines
#gene_list = c("gene", "ncRNA_gene", "tRNA_gene", "snoRNA_gene", "snRNA_gene", "rRNA_gene")
gene_list = c("gene")
gff.list.gene = gff.list[which(gff.list$type %in% gene_list)]

chr_list = c("NC_003279.8", "NC_003280.10", "NC_003281.10", "NC_003282.8", "NC_003283.11", "NC_003284.9")
names(chr_list) = c("I", "II", "III", "IV", "V", "X")

# Assigning strand
for(i in 1:nrow(mut.tab)) {
  chr = mut.tab$chromosome[i]
  chr = as.character(chr_list[chr])
  
  pos = mut.tab$gPos[i]
  
  tmp.tab = as.data.frame(gff.list.gene[which(gff.list.gene@seqnames == chr)])
  
  str = tmp.tab$strand[which(tmp.tab$start <= pos & tmp.tab$end >= pos)]
  
  if(isEmpty(str) == F) {
    mut.tab$strand[i] = str
  }
}

############################################################################################################################################
# Run code below if polarized data is needed - data must NOT be polarized for use in source_me.R

# # Changing mutations to correctly reflect strand
# for(i in 1:nrow(mut.tab)) {
#   if(mut.tab$strand[i] == "+" | mut.tab$strand[i] == "*") {
#     mut.tab$refBase[i] = mut.tab$from[i]
#     mut.tab$mutBase[i] = mut.tab$to[i]
#   }
#   
#   else {
#     tmpFrom = mut.tab$from[i]
#     tmpTo = mut.tab$to[i]
#     
#     if(tmpFrom == "A") {
#       mut.tab$refBase[i] = "T"
#     }
#     if(tmpTo == "A") {
#       mut.tab$mutBase[i] = "T"
#     }
#     
#     if(tmpFrom == "T") {
#       mut.tab$refBase[i] = "A"
#     }
#     if(tmpTo == "T") {
#       mut.tab$mutBase[i] = "A"
#     }
#     
#     if(tmpFrom == "G") {
#       mut.tab$refBase[i] = "C"
#     }
#     if(tmpTo == "G") {
#       mut.tab$mutBase[i] = "C"
#     }
#     
#     if(tmpFrom == "C") {
#       mut.tab$refBase[i] = "G"
#     }
#     if(tmpTo == "C") {
#       mut.tab$mutBase[i] = "G"
#     }
#   }
# }
# 
# mut_clean.tab = mut.tab[,c("RefSeqID", "sample", "genomic_position", "refBase", "mutBase", "gPos", "chromosome", "type")]
# colnames(mut_clean.tab) = c("RefSeqID", "sample", "genomic_positions", "from", "to", "gPos", "chromosome", "type")
# 
# ############################################################################################################################################
# # Writing polarized tables
# 
# write.table(mut_clean.tab, file = "Liu2021/mutations-second-no-mito-polar.tab", quote = F, sep = '\t', row.names = F, col.names = T)
