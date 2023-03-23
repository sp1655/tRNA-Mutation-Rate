# MUST BE RUN ONLY ON NON-POLARIZED DATASETS

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
library(foreach)
library(doParallel)
#library(ggplot2)
#library(tidyr)

source( paste(BD, "Rcode/RNAPIII-mutations-functions.R", sep="/") )

############################################################################################################################################
# Reading in the species information

fsp = paste(BD, "Data/species-list-tRNA.tab", sep="/") # <- File containing information about species
tsp = read.table(fsp, h=T, as.is=T, sep="\t")

############################################################################################################################################
# Running ALL datasets, including subsets

fx = paste(BD, "Data/MA_lines_tRNA.tab", sep="/") # <- File containing the list of MA lines to analyze
tx = read.table(fx, h=T, as.is=T, sep="\t")

# If dplyr is loaded, need to detach before running the loop
if("dplyr" %in% (.packages())){
  detach("package:dplyr", unload=TRUE) 
}

# Calling function for each dataset to get mutation data
for(i in (1:nrow(tx))){
  sp = tx$species[i]
  dataset = tx$dataset[i]
  fMut = paste(BD, "Data", tx$file[i], sep="/")
  
  wsp = which(tsp$species==sp)
  genomeFastaFile = paste(BD, "Data", tsp[wsp, "genomeFile"], sep="/")
  fChrOK = paste(BD, "Data", tsp[wsp, "chrOK_file"], sep="/")
  fPos = paste(BD, "Data", tsp[wsp, "rnap3_pos_file"], sep="/")
  
  cat(sp, "\t-\t", dataset, ":\n")
  vRes = RNAPIII_mutations_analysis(genomeFastaFile, fChrOK, fPos, fMut)
  cat("\n\n-------------------------\n\n")
}

############################################################################################################################################
# Running ONLY full datasets
# Note: need to make a separate MA line list containing only the full datasets

fx = paste(BD, "Data/MA_lines_tRNA_fullOnly.tab", sep="/") # <- File containing the list of MA lines to analyze
tx = read.table(fx, h=T, as.is=T, sep="\t")

speciesMatrix = matrix(nrow = nrow(tsp), ncol = 4)
colnames(speciesMatrix) = c("species", "observedTot", "expectedTotGC", "expectedTotNoGC")
speciesMatrix[,1] = tsp$species

species_tRNA.df = as.data.frame(speciesMatrix)

species_tRNA.df$observedTot = 0
species_tRNA.df$expectedTotGC = 0
species_tRNA.df$expectedTotNoGC = 0
species_tRNA.df$basesCovered = 0

speciesVec <- rep("", times = nrow(tsp))
spNum = length(speciesVec)

# Tables to store results from function
tRNA_mut.df = data.frame(species = character(spNum), AT = numeric(spNum), AG = numeric(spNum), AC = numeric(spNum), TA = numeric(spNum),
                         TG = numeric(spNum), TC = numeric(spNum), GA = numeric(spNum), GT = numeric(spNum), GC = numeric(spNum),
                         CA = numeric(spNum), CT = numeric(spNum), CG = numeric(spNum))

rest_mut.df = data.frame(species = character(spNum), AT = numeric(spNum), AG = numeric(spNum), AC = numeric(spNum), TA = numeric(spNum),
                         TG = numeric(spNum), TC = numeric(spNum), GA = numeric(spNum), GT = numeric(spNum), GC = numeric(spNum),
                         CA = numeric(spNum), CT = numeric(spNum), CG = numeric(spNum))

tRNA_bps.df = data.frame(species = character(spNum), AT_GC = numeric(spNum), GC_AT = numeric(spNum), AT_TA = numeric(spNum), 
                         AT_CG = numeric(spNum), GC_TA = numeric(spNum), GC_CG = numeric(spNum))

rest_bps.df = data.frame(species = character(spNum), AT_GC = numeric(spNum), GC_AT = numeric(spNum), AT_TA = numeric(spNum), 
                         AT_CG = numeric(spNum), GC_TA = numeric(spNum), GC_CG = numeric(spNum))

nFreqIII.df = data.frame(species = character(spNum), freq_A = numeric(spNum), freq_T = numeric(spNum), freq_G = numeric(spNum), 
                      freq_C = numeric(spNum))

nFreqRest.df = data.frame(species = character(spNum), freq_A = numeric(spNum), freq_T = numeric(spNum), freq_G = numeric(spNum), 
                         freq_C = numeric(spNum))

tRNA_mut.df$species = tsp$species
rest_mut.df$species = tsp$species

tRNA_bps.df$species = tsp$species
rest_bps.df$species = tsp$species

nFreqIII.df$species = tsp$species
nFreqRest.df$species = tsp$species

# Calling function for each dataset to get mutation data
for(i in (1:nrow(tx))){
  sp = tx$species[i]
  dataset = tx$dataset[i]
  fMut = paste(BD, "Data", tx$file[i], sep="/")
  
  wsp = which(tsp$species==sp)
  genomeFastaFile = paste(BD, "Data", tsp[wsp, "genomeFile"], sep="/")
  fChrOK = paste(BD, "Data", tsp[wsp, "chrOK_file"], sep="/")
  fPos = paste(BD, "Data", tsp[wsp, "rnap3_pos_file"], sep="/")
  
  cat(sp, "\t-\t", dataset, ":\n")
  vRes = RNAPIII_mutations_analysis(genomeFastaFile, fChrOK, fPos, fMut)
  cat("\n\n-------------------------\n\n")
  
  speciesRow = as.numeric(row.names(species_tRNA.df[which(species_tRNA.df$species == sp),]))
  species_tRNA.df$observedTot[speciesRow] = species_tRNA.df$observedTot[speciesRow] + vRes[1]
  species_tRNA.df$expectedTotGC[speciesRow] = species_tRNA.df$expectedTotGC[speciesRow] + vRes[2]
  species_tRNA.df$expectedTotNoGC[speciesRow] = species_tRNA.df$expectedTotNoGC[speciesRow] + vRes[3]
  species_tRNA.df$basesCovered[speciesRow] = species_tRNA.df$basesCovered[speciesRow] + vRes[4]
  
  tRNA_mut.df$AT[speciesRow] = tRNA_mut.df$AT[speciesRow] + vRes["AT_intr"]
  tRNA_mut.df$AG[speciesRow] = tRNA_mut.df$AG[speciesRow] + vRes["AG_intr"]
  tRNA_mut.df$AC[speciesRow] = tRNA_mut.df$AC[speciesRow] + vRes["AC_intr"]
  tRNA_mut.df$TA[speciesRow] = tRNA_mut.df$TA[speciesRow] + vRes["TA_intr"]
  tRNA_mut.df$TG[speciesRow] = tRNA_mut.df$TG[speciesRow] + vRes["TG_intr"]
  tRNA_mut.df$TC[speciesRow] = tRNA_mut.df$TC[speciesRow] + vRes["TC_intr"]
  tRNA_mut.df$GA[speciesRow] = tRNA_mut.df$GA[speciesRow] + vRes["GA_intr"]
  tRNA_mut.df$GT[speciesRow] = tRNA_mut.df$GT[speciesRow] + vRes["GT_intr"]
  tRNA_mut.df$GC[speciesRow] = tRNA_mut.df$GC[speciesRow] + vRes["GC_intr"]
  tRNA_mut.df$CA[speciesRow] = tRNA_mut.df$CA[speciesRow] + vRes["CA_intr"]
  tRNA_mut.df$CT[speciesRow] = tRNA_mut.df$CT[speciesRow] + vRes["CT_intr"]
  tRNA_mut.df$CG[speciesRow] = tRNA_mut.df$CG[speciesRow] + vRes["CG_intr"]
  
  rest_mut.df$AT[speciesRow] = rest_mut.df$AT[speciesRow] + vRes["AT_rest"]
  rest_mut.df$AG[speciesRow] = rest_mut.df$AG[speciesRow] + vRes["AG_rest"]
  rest_mut.df$AC[speciesRow] = rest_mut.df$AC[speciesRow] + vRes["AC_rest"]
  rest_mut.df$TA[speciesRow] = rest_mut.df$TA[speciesRow] + vRes["TA_rest"]
  rest_mut.df$TG[speciesRow] = rest_mut.df$TG[speciesRow] + vRes["TG_rest"]
  rest_mut.df$TC[speciesRow] = rest_mut.df$TC[speciesRow] + vRes["TC_rest"]
  rest_mut.df$GA[speciesRow] = rest_mut.df$GA[speciesRow] + vRes["GA_rest"]
  rest_mut.df$GT[speciesRow] = rest_mut.df$GT[speciesRow] + vRes["GT_rest"]
  rest_mut.df$GC[speciesRow] = rest_mut.df$GC[speciesRow] + vRes["GC_rest"]
  rest_mut.df$CA[speciesRow] = rest_mut.df$CA[speciesRow] + vRes["CA_rest"]
  rest_mut.df$CT[speciesRow] = rest_mut.df$CT[speciesRow] + vRes["CT_rest"]
  rest_mut.df$CG[speciesRow] = rest_mut.df$CG[speciesRow] + vRes["CG_rest"]
  
  tRNA_bps.df$AT_GC[speciesRow] = tRNA_bps.df$AT_GC[speciesRow] + vRes["AT_GC_intr"]
  tRNA_bps.df$GC_AT[speciesRow] = tRNA_bps.df$GC_AT[speciesRow] + vRes["GC_AT_intr"]
  tRNA_bps.df$AT_TA[speciesRow] = tRNA_bps.df$AT_TA[speciesRow] + vRes["AT_TA_intr"]
  tRNA_bps.df$AT_CG[speciesRow] = tRNA_bps.df$AT_CG[speciesRow] + vRes["AT_CG_intr"]
  tRNA_bps.df$GC_TA[speciesRow] = tRNA_bps.df$GC_TA[speciesRow] + vRes["GC_TA_intr"]
  tRNA_bps.df$GC_CG[speciesRow] = tRNA_bps.df$GC_CG[speciesRow] + vRes["GC_CG_intr"]
  
  rest_bps.df$AT_GC[speciesRow] = rest_bps.df$AT_GC[speciesRow] + vRes["AT_GC_rest"]
  rest_bps.df$GC_AT[speciesRow] = rest_bps.df$GC_AT[speciesRow] + vRes["GC_AT_rest"]
  rest_bps.df$AT_TA[speciesRow] = rest_bps.df$AT_TA[speciesRow] + vRes["AT_TA_rest"]
  rest_bps.df$AT_CG[speciesRow] = rest_bps.df$AT_CG[speciesRow] + vRes["AT_CG_rest"]
  rest_bps.df$GC_TA[speciesRow] = rest_bps.df$GC_TA[speciesRow] + vRes["GC_TA_rest"]
  rest_bps.df$GC_CG[speciesRow] = rest_bps.df$GC_CG[speciesRow] + vRes["GC_CG_rest"]
  
  nFreqIII.df$freq_A[speciesRow] = vRes["aFreq"]
  nFreqIII.df$freq_T[speciesRow] = vRes["tFreq"]
  nFreqIII.df$freq_G[speciesRow] = vRes["gFreq"]
  nFreqIII.df$freq_C[speciesRow] = vRes["cFreq"]
  
  nFreqRest.df$freq_A[speciesRow] = vRes["aFreqRest"]
  nFreqRest.df$freq_T[speciesRow] = vRes["tFreqRest"]
  nFreqRest.df$freq_G[speciesRow] = vRes["gFreqRest"]
  nFreqRest.df$freq_C[speciesRow] = vRes["cFreqRest"]
}

colnames(tRNA_mut.df) = c("species", "A>T", "A>G", "A>C", "T>A", "T>G", "T>C", "G>A", "G>T", "G>C", "C>A", "C>T", "C>G")
colnames(rest_mut.df) = c("species", "A>T", "A>G", "A>C", "T>A", "T>G", "T>C", "G>A", "G>T", "G>C", "C>A", "C>T", "C>G")

colnames(tRNA_bps.df) = c("species", "A:T>G:C", "G:C>A:T", "A:T>T:A", "A:T>C:G", "G:C>T:A", "G:C>C:G")
colnames(rest_bps.df) = c("species", "A:T>G:C", "G:C>A:T", "A:T>T:A", "A:T>C:G", "G:C>T:A", "G:C>C:G")

species_tRNA.df$observedTot = as.numeric(species_tRNA.df$observedTot)
species_tRNA.df$expectedTotGC = as.numeric(species_tRNA.df$expectedTotGC)
species_tRNA.df$expectedTotNoGC = as.numeric(species_tRNA.df$expectedTotNoGC)

#species_tRNA.df$ratioGC = species_tRNA.df$observedTot/species_tRNA.df$expectedTotGC
species_tRNA.df$mutRate = species_tRNA.df$observedTot/species_tRNA.df$basesCovered

species_tRNA.df$type = "tRNA"
species_tRNA.df$expLevel = "all"

species_tRNA.df = species_tRNA.df[,c(1,7,8,2,3,4,5,6)]
species_tRNA.df

############################################################################################################################################
# E. coli and S. cerevisiae tRNA by expression level (tpm)

fsp = paste(BD, "Data/species-list-tRNA-tpm.tab", sep="/") # <- File containing information about species
tsp = read.table(fsp, h=T, as.is=T, sep="\t")

fx = paste(BD, "Data/MA_lines_tRNA_tpm.tab", sep="/") # <- File containing the list of MA lines to analyze
tx = read.table(fx, h=T, as.is=T, sep="\t")

speciesMatrix = matrix(nrow = nrow(tsp), ncol = 4)
colnames(speciesMatrix) = c("species", "observedTot", "expectedTotGC", "expectedTotNoGC")
speciesMatrix[,1] = tsp$species

species_tRNA_tpm.df = as.data.frame(speciesMatrix)

species_tRNA_tpm.df$observedTot = 0
species_tRNA_tpm.df$expectedTotGC = 0
species_tRNA_tpm.df$expectedTotNoGC = 0
species_tRNA_tpm.df$basesCovered = 0

speciesVec <- rep("", times = nrow(tsp))

# Calling function for each dataset to get mutation data
for(i in (1:nrow(tx))){
  sp = tx$species[i]
  dataset = tx$dataset[i]
  fMut = paste(BD, "Data", tx$file[i], sep="/")
  
  wsp = which(tsp$species==sp)
  genomeFastaFile = paste(BD, "Data", tsp[wsp, "genomeFile"], sep="/")
  fChrOK = paste(BD, "Data", tsp[wsp, "chrOK_file"], sep="/")
  fPos = paste(BD, "Data", tsp[wsp, "rnap3_pos_file"], sep="/")
  
  cat(sp, "\t-\t", dataset, ":\n")
  vRes = RNAPIII_mutations_analysis(genomeFastaFile, fChrOK, fPos, fMut)
  cat("\n\n-------------------------\n\n")
  
  speciesRow = as.numeric(row.names(species_tRNA_tpm.df[which(species_tRNA_tpm.df$species == sp),]))
  species_tRNA_tpm.df$observedTot[speciesRow] = species_tRNA_tpm.df$observedTot[speciesRow] + vRes[1]
  species_tRNA_tpm.df$expectedTotGC[speciesRow] = species_tRNA_tpm.df$expectedTotGC[speciesRow] + vRes[2]
  species_tRNA_tpm.df$expectedTotNoGC[speciesRow] = species_tRNA_tpm.df$expectedTotNoGC[speciesRow] + vRes[3]
  species_tRNA_tpm.df$basesCovered[speciesRow] = species_tRNA_tpm.df$basesCovered[speciesRow] + vRes[4]
}

species_tRNA_tpm.df$observedTot = as.numeric(species_tRNA_tpm.df$observedTot)
species_tRNA_tpm.df$expectedTotGC = as.numeric(species_tRNA_tpm.df$expectedTotGC)
species_tRNA_tpm.df$expectedTotNoGC = as.numeric(species_tRNA_tpm.df$expectedTotNoGC)

#species_tRNA.df$ratioGC = species_tRNA.df$observedTot/species_tRNA.df$expectedTotGC
species_tRNA_tpm.df$mutRate = species_tRNA_tpm.df$observedTot/species_tRNA_tpm.df$basesCovered

species_tRNA_tpm.df$type = "tRNA"
species_tRNA_tpm.df$expLevel = "X"

for(i in 1:nrow(species_tRNA_tpm.df)) {
  tmp = unlist(strsplit(species_tRNA_tpm.df$species[i], "_tRNA_"))
  
  species_tRNA_tpm.df$species[i] = tmp[1]
  species_tRNA_tpm.df$expLevel[i] = tmp[2]
}

species_tRNA_tpm.df = species_tRNA_tpm.df[,c(1,7,8,2,3,4,5,6)]
species_tRNA_tpm.df

############################################################################################################################################
# Getting protein-coding and intergenic information for relevant species

########################################################################################
# Protein-coding

fsp = paste(BD, "Data/species-list-prot.tab", sep="/") # <- File containing information about species
tsp = read.table(fsp, h=T, as.is=T, sep="\t")

fx = paste(BD, "Data/MA_lines_prot.tab", sep="/") # <- File containing the list of MA lines to analyze
tx = read.table(fx, h=T, as.is=T, sep="\t")

speciesMatrix = matrix(nrow = nrow(tsp), ncol = 4)
colnames(speciesMatrix) = c("species", "observedTot", "expectedTotGC", "expectedTotNoGC")
speciesMatrix[,1] = tsp$species

species_prot.df = as.data.frame(speciesMatrix)

species_prot.df$observedTot = 0
species_prot.df$expectedTotGC = 0
species_prot.df$expectedTotNoGC = 0
species_prot.df$basesCovered = 0

speciesVec <- rep("", times = nrow(tsp))
spNum = length(speciesVec)

# Tables to store results from function
prot_mut.df = data.frame(species = character(spNum), AT = numeric(spNum), AG = numeric(spNum), AC = numeric(spNum), TA = numeric(spNum),
                         TG = numeric(spNum), TC = numeric(spNum), GA = numeric(spNum), GT = numeric(spNum), GC = numeric(spNum),
                         CA = numeric(spNum), CT = numeric(spNum), CG = numeric(spNum))

rest_prot_mut.df = data.frame(species = character(spNum), AT = numeric(spNum), AG = numeric(spNum), AC = numeric(spNum), TA = numeric(spNum),
                         TG = numeric(spNum), TC = numeric(spNum), GA = numeric(spNum), GT = numeric(spNum), GC = numeric(spNum),
                         CA = numeric(spNum), CT = numeric(spNum), CG = numeric(spNum))

prot_bps.df = data.frame(species = character(spNum), AT_GC = numeric(spNum), GC_AT = numeric(spNum), AT_TA = numeric(spNum), 
                         AT_CG = numeric(spNum), GC_TA = numeric(spNum), GC_CG = numeric(spNum))

rest_prot_mut_bps.df = data.frame(species = character(spNum), AT_GC = numeric(spNum), GC_AT = numeric(spNum), AT_TA = numeric(spNum), 
                         AT_CG = numeric(spNum), GC_TA = numeric(spNum), GC_CG = numeric(spNum))

nFreqProt.df = data.frame(species = character(spNum), freq_A = numeric(spNum), freq_T = numeric(spNum), freq_G = numeric(spNum), 
                         freq_C = numeric(spNum))

nFreqRestProt.df = data.frame(species = character(spNum), freq_A = numeric(spNum), freq_T = numeric(spNum), freq_G = numeric(spNum), 
                          freq_C = numeric(spNum))

prot_mut.df$species = tsp$species
rest_prot_mut.df$species = tsp$species

prot_bps.df$species = tsp$species
rest_prot_mut_bps.df$species = tsp$species

nFreqProt.df$species = tsp$species
nFreqRestProt.df$species = tsp$species

# Calling function for each dataset to get mutation data
for(i in (1:nrow(tx))){
  sp = tx$species[i]
  dataset = tx$dataset[i]
  fMut = paste(BD, "Data", tx$file[i], sep="/")
  
  wsp = which(tsp$species==sp)
  genomeFastaFile = paste(BD, "Data", tsp[wsp, "genomeFile"], sep="/")
  fChrOK = paste(BD, "Data", tsp[wsp, "chrOK_file"], sep="/")
  fPos = paste(BD, "Data", tsp[wsp, "rnap3_pos_file"], sep="/")
  
  cat(sp, "\t-\t", dataset, ":\n")
  vRes = RNAPIII_mutations_analysis(genomeFastaFile, fChrOK, fPos, fMut)
  cat("\n\n-------------------------\n\n")
  
  speciesRow = as.numeric(row.names(species_prot.df[which(species_prot.df$species == sp),]))
  species_prot.df$observedTot[speciesRow] = species_prot.df$observedTot[speciesRow] + vRes[1]
  species_prot.df$expectedTotGC[speciesRow] = species_prot.df$expectedTotGC[speciesRow] + vRes[2]
  species_prot.df$expectedTotNoGC[speciesRow] = species_prot.df$expectedTotNoGC[speciesRow] + vRes[3]
  species_prot.df$basesCovered[speciesRow] = species_prot.df$basesCovered[speciesRow] + vRes[4]
  
  prot_mut.df$AT[speciesRow] = prot_mut.df$AT[speciesRow] + vRes["AT_intr"]
  prot_mut.df$AG[speciesRow] = prot_mut.df$AG[speciesRow] + vRes["AG_intr"]
  prot_mut.df$AC[speciesRow] = prot_mut.df$AC[speciesRow] + vRes["AC_intr"]
  prot_mut.df$TA[speciesRow] = prot_mut.df$TA[speciesRow] + vRes["TA_intr"]
  prot_mut.df$TG[speciesRow] = prot_mut.df$TG[speciesRow] + vRes["TG_intr"]
  prot_mut.df$TC[speciesRow] = prot_mut.df$TC[speciesRow] + vRes["TC_intr"]
  prot_mut.df$GA[speciesRow] = prot_mut.df$GA[speciesRow] + vRes["GA_intr"]
  prot_mut.df$GT[speciesRow] = prot_mut.df$GT[speciesRow] + vRes["GT_intr"]
  prot_mut.df$GC[speciesRow] = prot_mut.df$GC[speciesRow] + vRes["GC_intr"]
  prot_mut.df$CA[speciesRow] = prot_mut.df$CA[speciesRow] + vRes["CA_intr"]
  prot_mut.df$CT[speciesRow] = prot_mut.df$CT[speciesRow] + vRes["CT_intr"]
  prot_mut.df$CG[speciesRow] = prot_mut.df$CG[speciesRow] + vRes["CG_intr"]
  
  rest_prot_mut.df$AT[speciesRow] = rest_prot_mut.df$AT[speciesRow] + vRes["AT_rest"]
  rest_prot_mut.df$AG[speciesRow] = rest_prot_mut.df$AG[speciesRow] + vRes["AG_rest"]
  rest_prot_mut.df$AC[speciesRow] = rest_prot_mut.df$AC[speciesRow] + vRes["AC_rest"]
  rest_prot_mut.df$TA[speciesRow] = rest_prot_mut.df$TA[speciesRow] + vRes["TA_rest"]
  rest_prot_mut.df$TG[speciesRow] = rest_prot_mut.df$TG[speciesRow] + vRes["TG_rest"]
  rest_prot_mut.df$TC[speciesRow] = rest_prot_mut.df$TC[speciesRow] + vRes["TC_rest"]
  rest_prot_mut.df$GA[speciesRow] = rest_prot_mut.df$GA[speciesRow] + vRes["GA_rest"]
  rest_prot_mut.df$GT[speciesRow] = rest_prot_mut.df$GT[speciesRow] + vRes["GT_rest"]
  rest_prot_mut.df$GC[speciesRow] = rest_prot_mut.df$GC[speciesRow] + vRes["GC_rest"]
  rest_prot_mut.df$CA[speciesRow] = rest_prot_mut.df$CA[speciesRow] + vRes["CA_rest"]
  rest_prot_mut.df$CT[speciesRow] = rest_prot_mut.df$CT[speciesRow] + vRes["CT_rest"]
  rest_prot_mut.df$CG[speciesRow] = rest_prot_mut.df$CG[speciesRow] + vRes["CG_rest"]
  
  prot_bps.df$AT_GC[speciesRow] = prot_bps.df$AT_GC[speciesRow] + vRes["AT_GC_intr"]
  prot_bps.df$GC_AT[speciesRow] = prot_bps.df$GC_AT[speciesRow] + vRes["GC_AT_intr"]
  prot_bps.df$AT_TA[speciesRow] = prot_bps.df$AT_TA[speciesRow] + vRes["AT_TA_intr"]
  prot_bps.df$AT_CG[speciesRow] = prot_bps.df$AT_CG[speciesRow] + vRes["AT_CG_intr"]
  prot_bps.df$GC_TA[speciesRow] = prot_bps.df$GC_TA[speciesRow] + vRes["GC_TA_intr"]
  prot_bps.df$GC_CG[speciesRow] = prot_bps.df$GC_CG[speciesRow] + vRes["GC_CG_intr"]
  
  rest_prot_mut_bps.df$AT_GC[speciesRow] = rest_prot_mut_bps.df$AT_GC[speciesRow] + vRes["AT_GC_rest"]
  rest_prot_mut_bps.df$GC_AT[speciesRow] = rest_prot_mut_bps.df$GC_AT[speciesRow] + vRes["GC_AT_rest"]
  rest_prot_mut_bps.df$AT_TA[speciesRow] = rest_prot_mut_bps.df$AT_TA[speciesRow] + vRes["AT_TA_rest"]
  rest_prot_mut_bps.df$AT_CG[speciesRow] = rest_prot_mut_bps.df$AT_CG[speciesRow] + vRes["AT_CG_rest"]
  rest_prot_mut_bps.df$GC_TA[speciesRow] = rest_prot_mut_bps.df$GC_TA[speciesRow] + vRes["GC_TA_rest"]
  rest_prot_mut_bps.df$GC_CG[speciesRow] = rest_prot_mut_bps.df$GC_CG[speciesRow] + vRes["GC_CG_rest"]
  
  nFreqProt.df$freq_A[speciesRow] = vRes["aFreq"]
  nFreqProt.df$freq_T[speciesRow] = vRes["tFreq"]
  nFreqProt.df$freq_G[speciesRow] = vRes["gFreq"]
  nFreqProt.df$freq_C[speciesRow] = vRes["cFreq"]
  
  nFreqRestProt.df$freq_A[speciesRow] = vRes["aFreqRest"]
  nFreqRestProt.df$freq_T[speciesRow] = vRes["tFreqRest"]
  nFreqRestProt.df$freq_G[speciesRow] = vRes["gFreqRest"]
  nFreqRestProt.df$freq_C[speciesRow] = vRes["cFreqRest"]
}

species_prot.df$observedTot = as.numeric(species_prot.df$observedTot)
species_prot.df$expectedTotGC = as.numeric(species_prot.df$expectedTotGC)
species_prot.df$expectedTotNoGC = as.numeric(species_prot.df$expectedTotNoGC)

#species_prot.df$ratioGC = species_prot.df$observedTot/species_prot.df$expectedTotGC
species_prot.df$mutRate = species_prot.df$observedTot/species_prot.df$basesCovered

library(stringr)

# Splitting species column into species name and expression level
species_prot.df$expLevel = "X"
species_prot.df$type = "protein_coding"
for(i in 1:nrow(species_prot.df)) {
  tmp = str_split(species_prot.df$species, "_prot_")[i]
  tmp = unlist(tmp)
  
  species_prot.df$species[i] = tmp[1]
  species_prot.df$expLevel[i] = tmp[2]
}

species_prot.df = species_prot.df[,c(1,8,7,2,3,4,5,6)]
species_prot.df

########################################################################################
# Intergenic

fsp = paste(BD, "Data/species-list-inter.tab", sep="/") # <- File containing information about species
tsp = read.table(fsp, h=T, as.is=T, sep="\t")

fx = paste(BD, "Data/MA_lines_inter.tab", sep="/") # <- File containing the list of MA lines to analyze
tx = read.table(fx, h=T, as.is=T, sep="\t")

speciesMatrix = matrix(nrow = nrow(tsp), ncol = 4)
colnames(speciesMatrix) = c("species", "observedTot", "expectedTotGC", "expectedTotNoGC")
speciesMatrix[,1] = tsp$species

species_inter.df = as.data.frame(speciesMatrix)

species_inter.df$observedTot = 0
species_inter.df$expectedTotGC = 0
species_inter.df$expectedTotNoGC = 0
species_inter.df$basesCovered = 0

speciesVec <- rep("", times = nrow(tsp))

# Calling function for each dataset to get mutation data
for(i in (1:nrow(tx))){
  sp = tx$species[i]
  dataset = tx$dataset[i]
  fMut = paste(BD, "Data", tx$file[i], sep="/")
  
  wsp = which(tsp$species==sp)
  genomeFastaFile = paste(BD, "Data", tsp[wsp, "genomeFile"], sep="/")
  fChrOK = paste(BD, "Data", tsp[wsp, "chrOK_file"], sep="/")
  fPos = paste(BD, "Data", tsp[wsp, "rnap3_pos_file"], sep="/")
  
  cat(sp, "\t-\t", dataset, ":\n")
  vRes = RNAPIII_mutations_analysis(genomeFastaFile, fChrOK, fPos, fMut)
  cat("\n\n-------------------------\n\n")
  
  speciesRow = as.numeric(row.names(species_inter.df[which(species_inter.df$species == sp),]))
  species_inter.df$observedTot[speciesRow] = species_inter.df$observedTot[speciesRow] + vRes[1]
  species_inter.df$expectedTotGC[speciesRow] = species_inter.df$expectedTotGC[speciesRow] + vRes[2]
  species_inter.df$expectedTotNoGC[speciesRow] = species_inter.df$expectedTotNoGC[speciesRow] + vRes[3]
  species_inter.df$basesCovered[speciesRow] = species_inter.df$basesCovered[speciesRow] + vRes[4]
}

species_inter.df$observedTot = as.numeric(species_inter.df$observedTot)
species_inter.df$expectedTotGC = as.numeric(species_inter.df$expectedTotGC)
species_inter.df$expectedTotNoGC = as.numeric(species_inter.df$expectedTotNoGC)

#species_inter.df$ratioGC = species_inter.df$observedTot/species_inter.df$expectedTotGC
species_inter.df$mutRate = species_inter.df$observedTot/species_inter.df$basesCovered

species_inter.df$type = "intergenic"
species_inter.df$expLevel = "all"
species_inter.df$strain = "X"
for(i in 1:nrow(species_inter.df)) {
  tmp = strsplit(species_inter.df$species, "_inter")[i]
  tmp = unlist(tmp)
  
  species_inter.df$species[i] = tmp[1]
  
  if(length(tmp) > 1) {
    tmp2 = unlist(strsplit(tmp[2], "_"))
    species_inter.df$strain[i] = tmp2[2]
  }
  
  else {
    species_inter.df$strain[i] = "X"
  }
}

species_inter.df = species_inter.df[,c(1,7,8,2,3,4,5,6,9)]
species_inter.df

############################################################################################################################################
# Writing the tables

write.table(species_tRNA.df, file = "Data/species_tRNA.tab", sep = '\t', quote = F, row.names = F, col.names = T)
write.table(tRNA_mut.df, file = "Data/mutSpec_tRNA.tab", sep = '\t', quote = F, row.names = F, col.names = T)
write.table(rest_mut.df, file = "Data/mutSpec_rest_tRNA.tab", sep = '\t', quote = F, row.names = F, col.names = T)
write.table(tRNA_bps.df, file = "Data/BPS_tRNA.tab", sep = '\t', quote = F, row.names = F, col.names = T)
write.table(rest_bps.df, file = "Data/BPS_rest_tRNA.tab", sep = '\t', quote = F, row.names = F, col.names = T)
write.table(nFreqIII.df, file = "Data/ntFreq_tRNA.tab", sep = '\t', quote = F, row.names = F, col.names = T)
write.table(nFreqRest.df, file = "Data/ntFreq_rest_tRNA.tab", sep = '\t', quote = F, row.names = F, col.names = T)

write.table(species_tRNA_tpm.df, file = "Data/species_tRNA_tpm.tab", sep = '\t', quote = F, row.names = F, col.names = T)

write.table(species_prot.df, file = "Data/species_prot.tab", sep = '\t', quote = F, row.names = F, col.names = T)
write.table(prot_mut.df, file = "Data/mutSpec_prot.tab", sep = '\t', quote = F, row.names = F, col.names = T)
write.table(rest_prot_mut.df, file = "Data/mutSpec_rest_prot.tab", sep = '\t', quote = F, row.names = F, col.names = T)
write.table(prot_bps.df, file = "Data/BPS_prot.tab", sep = '\t', quote = F, row.names = F, col.names = T)
write.table(rest_prot_mut_bps.df, file = "Data/BPS_rest_prot.tab", sep = '\t', quote = F, row.names = F, col.names = T)
write.table(nFreqProt.df, file = "Data/ntFreq_prot.tab", sep = '\t', quote = F, row.names = F, col.names = T)
write.table(nFreqRestProt.df, file = "Data/ntFreq_rest_prot.tab", sep = '\t', quote = F, row.names = F, col.names = T)

write.table(species_inter.df, file = "Data/species_inter.tab", sep = '\t', quote = F, row.names = F, col.names = T)



############################################################################################################################################
# Running ALL datasets, including subsets

# library(GenomicRanges)
# library(Biostrings)
# library(BSgenome)
# library(genomation)
# 
# rm(list = ls())
# 
# BD = "C:/Users/Sarah/Desktop/Lab/RNAPIII/RNAPIII-mutation-rate-main"
# 
# source( paste(BD, "Rcode/RNAPIII-mutations-functions.R", sep="/") )
# 
# fsp = paste(BD, "Data/species-tmp.tab", sep="/") # <- File containing information about species
# tsp = read.table(fsp, h=T, as.is=T, sep="\t")
# 
# fx = paste(BD, "Data/MA_tmp.tab", sep="/") # <- File containing the list of MA lines to analyze
# tx = read.table(fx, h=T, as.is=T, sep="\t")
# 
# # If dplyr is loaded, need to detach before running the loop
# if("dplyr" %in% (.packages())){
#   detach("package:dplyr", unload=TRUE) 
# }
# 
# for(i in (1:nrow(tx))){
#   sp = tx$species[i]
#   dataset = tx$dataset[i]
#   fMut = paste(BD, "Data", tx$file[i], sep="/")
#   
#   wsp = which(tsp$species==sp)
#   genomeFastaFile = paste(BD, "Data", tsp[wsp, "genomeFile"], sep="/")
#   fChrOK = paste(BD, "Data", tsp[wsp, "chrOK_file"], sep="/")
#   fPos = paste(BD, "Data", tsp[wsp, "rnap3_pos_file"], sep="/")
#   
#   cat(sp, "\t-\t", dataset, ":\n")
#   vRes = RNAPIII_mutations_analysis(genomeFastaFile, fChrOK, fPos, fMut)
#   cat("\n\n-------------------------\n\n")
# }
