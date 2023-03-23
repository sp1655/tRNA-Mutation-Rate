# Code for checking that the dataset and the assembly line up properly

rm(list=ls())

library(GenomicRanges)
library(Biostrings)
library(BSgenome)
library(genomation)

setwd("C:/Users/Sarah/Desktop/Lab/RNAPIII/RNAPIII-mutation-rate-main/Data/Ecoli/")

fMut = "Foster2018/mutAll.tab"
genomeFastaFile = "DB/NC_000913.2.fasta"

vBases = c("A", "T", "C", "G")

# Reading the table of mutations from the MA line ("from" column contains the value for the non-mutated allele)
tm = read.table(fMut, h=T, as.is=T, sep="\t")
gm = makeGRangesFromDataFrame(tm, keep.extra.columns = TRUE, start.field = "gPos", end.field = "gPos")

# Loading the genome in R
genome = readDNAStringSet(genomeFastaFile)

# Making genome name match seqlevels in gm (U00096.2 and NC000913.2 are the same)
names(genome) = "U00096.2"

# Had to do this to get it to work
##names(genome) = "I"
for(i in (1:length(names(genome)))){
  fullID = names(genome)[i]
  chrID = unlist(strsplit(fullID, " "))[1]
  names(genome)[i] = chrID
}

# Making sure that the base calls in the MA line match that of the reference genome used here:
for(base in vBases){
  gmb = gm[which(gm$from==base)]
  seq = getSeq(genome, gmb)
  if( length(which(seq==base)) != length(seq) ){
    cat("WARNING: out of ", length(seq), " ", base, " in the MA line data, only ", length(which(seq==base)), 
        " match in the genome.\n")
    cat("Check if mutations are given in zero-based or 1-based coordinates + check assembly version.\n")
  }
}

# Use if above gave multiple warnings
for(i in 1:nrow(tm)) {
  tm$gPos[i] = tm$gPos[i] - 1
}
