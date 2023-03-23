# Loads and splits Vibrio data by species into two separate file

############################################################################################################################################
# Loading and editing the Vibrio data

rm(list = ls())

setwd("c:/Users/Sarah/Desktop/Lab/RNAPIII/RNAPIII-mutation-rate-main/Data/")

library(readxl)

cwt = read_excel("Vcholerae/Dillon2016/Supplementary_Dataset_new.xlsx", sheet = "V. cholerae WT")
cmut = read_excel("Vcholerae/Dillon2016/Supplementary_Dataset_new.xlsx", sheet = "V. cholerae mut")

fwt = read_excel("Vfischeri/Dillon2016/Supplementary_Dataset_new.xlsx", sheet = "V. fischeri WT")
fmut = read_excel("Vfischeri/Dillon2016/Supplementary_Dataset_new.xlsx", sheet = "V. fischeri mut")

cwt = as.data.frame(cwt)
cmut = as.data.frame(cmut)
fwt = as.data.frame(fwt)
fmut = as.data.frame(fmut)

############################################################################################################################################
# Splitting by species to avoid confusion with genome names

# V. cholerae
genomeFileC = "Vcholerae/DB/GCF_001683415.1_ASM168341v1_genomic.fna.gz"
genomeC = readDNAStringSet(genomeFileC)

for(i in (1:length(names(genomeC)))){
  fullID = names(genomeC)[i]
  chrID = unlist(strsplit(fullID, " "))[1]
  names(genomeC)[i] = chrID
}

chr1 = names(genomeC)[1]
chr2 = names(genomeC)[2]

# Replacing chromosome names in both V. cholerae files
for(i in 1:nrow(cwt)) {
  if(cwt$chromosome[i] == "1") {
    cwt$chromosome[i] = chr1
  }
  
  else if(cwt$chromosome[i] == "2") {
    cwt$chromosome[i] = chr2
  }
}

for(i in 1:nrow(cmut)) {
  if(cmut$chromosome[i] == "1") {
    cmut$chromosome[i] = chr1
  }
  
  else if(cmut$chromosome[i] == "2") {
    cmut$chromosome[i] = chr2
  }
}

#########################################################################################################################
# V. fischeri
genomeFileF = "Vfischeri/DB/GCF_000011805.1_ASM1180v1_genomic.fna.gz"
genomeF = readDNAStringSet(genomeFileF)

for(i in (1:length(names(genomeF)))){
  fullID = names(genomeF)[i]
  chrID = unlist(strsplit(fullID, " "))[1]
  names(genomeF)[i] = chrID
}

chr1 = names(genomeF)[1]
chr2 = names(genomeF)[2]
chrp = names(genomeF)[3]

# Replacing chromosome names in both V. fischeri files
for(i in 1:nrow(fwt)) {
  if(fwt$chromosome[i] == "1") {
    fwt$chromosome[i] = chr1
  }
  
  else if(fwt$chromosome[i] == "2") {
    fwt$chromosome[i] = chr2
  }
  
  else if(fwt$chromosome[i] == "P") {
    fwt$chromosome[i] = chrp
  }
}


for(i in 1:nrow(fmut)) {
  if(fmut$chromosome[i] == "1") {
    fmut$chromosome[i] = chr1
  }
  
  else if(fmut$chromosome[i] == "2") {
    fmut$chromosome[i] = chr2
  }
  
  else if(fmut$chromosome[i] == "P") {
    fmut$chromosome[i] = chrp
  }
}

############################################################################################################################################
# Saving the tables
cwt$type = "BASE_SUB"
cmut$type = "BASE_SUB"
fwt$type = "BASE_SUB"
fmut$type = "BASE_SUB"

write.table(cwt, file = "Vcholerae/Dillon2016/wt_mut.tab", sep = '\t', quote = F, row.names = F, col.names = T)
write.table(cmut, file = "Vcholerae/Dillon2016/mmr_mut.tab", sep = '\t', quote = F, row.names = F, col.names = T)
write.table(fwt, file = "Vfischeri/Dillon2016/wt_mut.tab", sep = '\t', quote = F, row.names = F, col.names = T)
write.table(fmut, file = "Vfischeri/Dillon2016/mmr_mut.tab", sep = '\t', quote = F, row.names = F, col.names = T)
