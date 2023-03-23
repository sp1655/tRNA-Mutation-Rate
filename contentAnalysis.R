################################################################################################
#                                                                                              #
#  Program description:                                                                        #
#                                                                                              #
#  Checking mutations to make sure nothing is strange with the data that would affect results  #
#  Checks number of mutations found in genes on the +/- strands and compares frequencies       #
#  Checks frequency of G->T mutations in dataset                                               #
#  Frequency output can be checked for mutation distribution abnormalities                     #
#                                                                                              #
################################################################################################

rm(list=ls())

BD = "C:/Users/Sarah/Desktop/Lab/RNAPIII/RNAPIII-mutation-rate-main" # <- This should be changed for your local computer.
# The next few lines will attempt to automatically find the correct location:

se <- Sys.getenv()
if( as.numeric(se["RSTUDIO"]) == 1 ){
  # Running in RStudio
  BD = paste(dirname(rstudioapi::getActiveDocumentContext()$path), "..", sep="/")
} else {
  # Running outside RStudio
  BD = paste(getSrcDirectory()[1], "..", sep="/")
}

setwd(BD)

library(GenomicRanges)
library(Biostrings)
library(BSgenome)
library(genomation)

source( paste(BD, "Rcode/contentAnalysisFunctions.R", sep="/") )

############################################################################################################################################
# Reading in the species information

fsp = paste(BD, "Data/species-list-tRNA.tab", sep="/") # <- File containing information about species
tsp = read.table(fsp, h=T, as.is=T, sep="\t")

############################################################################################################################################
# Running E. coli data to check G->T content

fx = paste(BD, "Data/MAecoliOnly.tab", sep="/")
tx = read.table(fx, h=T, as.is=T, sep="\t")

# If dplyr is loaded, need to detach before running the loop
if("dplyr" %in% (.packages())){
  detach("package:dplyr", unload=TRUE) 
}

mg1655.tab = data.frame(matrix(ncol = 5, nrow = 0))
atcc.tab = data.frame(matrix(ncol = 5, nrow = 0))

for(i in (1:nrow(tx))){
  sp = tx$species[i]
  dataset = tx$dataset[i]
  fMut = paste(BD, "Data", tx$file[i], sep="/")
  
  wsp = which(tsp$species==sp)
  genomeFastaFile = paste(BD, "Data", tsp[wsp, "genomeFile"], sep="/")
  fChrOK = paste(BD, "Data", tsp[wsp, "chrOK_file"], sep="/")
  fPos = paste(BD, "Data", tsp[wsp, "rnap3_pos_file"], sep="/")
  
  cat(sp, "\t-\t", dataset, ":\n")
  vRes = contentAnalysis(genomeFastaFile, fChrOK, fPos, fMut)
  cat("\n\n-------------------------\n\n")
  
  if(vRes["nbMut"] != 0) {
    lsTmp = c(dataset, vRes["GTfreqPolar"], vRes["GTfreqNonpolar"], vRes["minusFreq"], vRes["nbMut"])
    
    if(sp == "Escherichia_coli_mg1655") {
      mg1655.tab = rbind(mg1655.tab, lsTmp)
    }
  
    if(sp == "Escherichia_coli_atcc") {
      atcc.tab = rbind(atcc.tab, lsTmp)
    }
  }
}

colnames(mg1655.tab) = c("dataset", "GTfreqPolar", "GTfreqNonpolar", "minusFreq", "nbMut")
colnames(atcc.tab) = c("dataset", "GTfreqPolar", "GTfreqNonpolar", "minusFreq", "nbMut")

mg1655.tab$GTfreqPolar = as.numeric(mg1655.tab$GTfreqPolar)
mg1655.tab$GTfreqNonpolar = as.numeric(mg1655.tab$GTfreqNonpolar)
mg1655.tab$minusFreq = as.numeric(mg1655.tab$minusFreq)
mg1655.tab$nbMut = as.numeric(mg1655.tab$nbMut)

atcc.tab$GTfreqPolar = as.numeric(atcc.tab$GTfreqPolar)
atcc.tab$GTfreqNonpolar = as.numeric(atcc.tab$GTfreqNonpolar)
atcc.tab$minusFreq = as.numeric(atcc.tab$minusFreq)
atcc.tab$nbMut = as.numeric(atcc.tab$nbMut)

# Getting averages of each column
sum(mg1655.tab$GTfreqPolar)/nrow(mg1655.tab)
sum(mg1655.tab$GTfreqNonpolar)/nrow(mg1655.tab)
sum(mg1655.tab$minusFreq)/nrow(mg1655.tab)
sum(mg1655.tab$nbMut)

sum(atcc.tab$GTfreqPolar)/nrow(atcc.tab)
sum(atcc.tab$GTfreqNonpolar)/nrow(atcc.tab)
sum(atcc.tab$minusFreq)/nrow(atcc.tab)
sum(atcc.tab$nbMut)

############################################################################################################################################
# Running S. cerevisiae data to check G->T content

fx = paste(BD, "Data/MAcerevisiaeOnly.tab", sep="/")
tx = read.table(fx, h=T, as.is=T, sep="\t")

# If dplyr is loaded, need to detach before running the loop
if("dplyr" %in% (.packages())){
  detach("package:dplyr", unload=TRUE) 
}

scerevisiae.tab = data.frame(matrix(ncol = 5, nrow = 0))

for(i in (1:nrow(tx))){
  sp = tx$species[i]
  dataset = tx$dataset[i]
  fMut = paste(BD, "Data", tx$file[i], sep="/")
  
  wsp = which(tsp$species==sp)
  genomeFastaFile = paste(BD, "Data", tsp[wsp, "genomeFile"], sep="/")
  fChrOK = paste(BD, "Data", tsp[wsp, "chrOK_file"], sep="/")
  fPos = paste(BD, "Data", tsp[wsp, "rnap3_pos_file"], sep="/")
  
  cat(sp, "\t-\t", dataset, ":\n")
  vRes = contentAnalysis(genomeFastaFile, fChrOK, fPos, fMut)
  cat("\n\n-------------------------\n\n")
  
  if(vRes["nbMut"] != 0) {
    lsTmp = c(dataset, vRes["GTfreqPolar"], vRes["GTfreqNonpolar"], vRes["minusFreq"], vRes["nbMut"])
    scerevisiae.tab = rbind(scerevisiae.tab, lsTmp)
  }
}

colnames(scerevisiae.tab) = c("dataset", "GTfreqPolar", "GTfreqNonpolar", "minusFreq", "nbMut")

scerevisiae.tab$GTfreqPolar = as.numeric(scerevisiae.tab$GTfreqPolar)
scerevisiae.tab$GTfreqNonpolar = as.numeric(scerevisiae.tab$GTfreqNonpolar)
scerevisiae.tab$minusFreq = as.numeric(scerevisiae.tab$minusFreq)
scerevisiae.tab$nbMut = as.numeric(scerevisiae.tab$nbMut)

# Getting averages of each column
sum(scerevisiae.tab$GTfreqPolar)/nrow(scerevisiae.tab)
sum(scerevisiae.tab$GTfreqNonpolar)/nrow(scerevisiae.tab)
sum(scerevisiae.tab$minusFreq)/nrow(scerevisiae.tab)
sum(scerevisiae.tab$nbMut)
