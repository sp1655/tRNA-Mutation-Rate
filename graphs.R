# Making graphs from the data table outputs of source_me.R

rm(list = ls())

library(ggplot2)
library(ggpattern)
library(tidyr)
library(stringr)
library(BSgenome)
library(genomation)

setwd("C:/Users/Sarah/Desktop/Lab/RNAPIII/RNAPIII-mutation-rate-main/Data/")

species_tRNA.df = read.table("species_tRNA.tab", h = T)
species_tRNA_tpm.df = read.table("species_tRNA_tpm.tab", h = T)
species_prot.df = read.table("species_prot.tab", h = T)
species_inter.df = read.table("species_inter.tab", h = T)

# tRNA data
tRNA_mut.df = read.table("mutSpec_tRNA.tab", h = T)
rest_tRNA_mut.df = read.table("mutSpec_rest_tRNA.tab", h = T)

colnames(tRNA_mut.df) = c("species", "A>T", "A>G", "A>C", "T>A", "T>G", "T>C", "G>A", "G>T", "G>C", "C>A", "C>T", "C>G")
colnames(rest_tRNA_mut.df) = c("species", "A>T", "A>G", "A>C", "T>A", "T>G", "T>C", "G>A", "G>T", "G>C", "C>A", "C>T", "C>G")

tRNA_bps.df = read.table("BPS_tRNA.tab", h = T)
rest_tRNA_bps.df = read.table("BPS_rest_tRNA.tab", h = T)

colnames(tRNA_bps.df) = c("species", "A:T>G:C", "G:C>A:T", "A:T>T:A", "A:T>C:G", "G:C>T:A", "G:C>C:G")
colnames(rest_tRNA_bps.df) = c("species", "A:T>G:C", "G:C>A:T", "A:T>T:A", "A:T>C:G", "G:C>T:A", "G:C>C:G")

# Protein-coding data
prot_mut.df = read.table("mutSpec_prot.tab", h = T)
rest_prot_mut.df = read.table("mutSpec_rest_prot.tab", h = T)

colnames(prot_mut.df) = c("species", "A>T", "A>G", "A>C", "T>A", "T>G", "T>C", "G>A", "G>T", "G>C", "C>A", "C>T", "C>G")
colnames(rest_prot_mut.df) = c("species", "A>T", "A>G", "A>C", "T>A", "T>G", "T>C", "G>A", "G>T", "G>C", "C>A", "C>T", "C>G")

prot_bps.df = read.table("BPS_prot.tab", h = T)
rest_prot_bps.df = read.table("BPS_rest_prot.tab", h = T)

colnames(prot_bps.df) = c("species", "A:T>G:C", "G:C>A:T", "A:T>T:A", "A:T>C:G", "G:C>T:A", "G:C>C:G")
colnames(rest_prot_bps.df) = c("species", "A:T>G:C", "G:C>A:T", "A:T>T:A", "A:T>C:G", "G:C>T:A", "G:C>C:G")

# Color palettes
cbPalette <- c("#BDD9BF", "#A997DF", "#FFC857", "#E5323B")
bpsPalette <- c("#A997DF", "#E5323B")

############################################################################################################################################
# Making the graphs for species totals (S. cerevisiae, E. coli, S. pombe)

# Combining all the E. coli data
er1 = as.numeric(rownames(species_tRNA.df[which(species_tRNA.df$species == "Escherichia_coli_mg1655"),]))
er2 = as.numeric(rownames(species_tRNA.df[which(species_tRNA.df$species == "Escherichia_coli_atcc"),]))

ecoli_combined <- c("Escherichia_coli", species_tRNA.df$observedTot[er1] + species_tRNA.df$observedTot[er2], 
                    species_tRNA.df$expectedTotGC[er1] + species_tRNA.df$expectedTotGC[er2], 
                    species_tRNA.df$expectedTotNoGC[er1] + species_tRNA.df$expectedTotNoGC[er2])

ecoli_combined[5] = as.numeric(ecoli_combined[2])/as.numeric(ecoli_combined[3])

species_tRNA.df = rbind(species_tRNA.df, ecoli_combined)

# Keeping only the species with sufficient data
speciesKept = c("Saccharomyces_cerevisiae", "Saccharomyces_pombe", "Escherichia_coli")
species_tRNA_trimmed.df = species_tRNA.df[which(species_tRNA.df$species %in% speciesKept),]

species_tRNA_trimmed.df$observedTot = as.numeric(species_tRNA_trimmed.df$observedTot)
species_tRNA_trimmed.df$expectedTotGC = as.numeric(species_tRNA_trimmed.df$expectedTotGC)

df <- gather(species_tRNA_trimmed.df, event, total, observedTot:expectedTotGC)

# Making the plot
plot <- ggplot(df, aes(species, total, fill = event)) +
  geom_bar(stat = "identity", position = 'dodge') +
  ggtitle("Observed vs. Expected Mutations in tRNA") + 
  labs(x = "Species",y = "Count")
plot

############################################################################################################################################
# S. cerevisiae tRNA data only

scRow = as.numeric(rownames(species_tRNA.df[which(species_tRNA.df$species == "Saccharomyces_cerevisiae"),]))
sc_data = species_tRNA.df[scRow,]

# Making sure the columns contain numeric data
sc_data$observedTot = as.numeric(sc_data$observedTot)
sc_data$expectedTotGC = as.numeric(sc_data$expectedTotGC)
sc_data$expectedTotNoGC = as.numeric(sc_data$expectedTotNoGC)

vNum <- c(numObs = 0, expGC = 0, expNoGC = 0)

# Adding up totals to put in vNum
totNumObs = 0
totExpGC = 0
totExpNoGC = 0
for(i in 1:nrow(sc_data)) {
  totNumObs = totNumObs + sc_data$observedTot[i]
  totExpGC = totExpGC + sc_data$expectedTotGC[i]
  totExpNoGC = totExpNoGC + sc_data$expectedTotNoGC[i]
}

vNum["numObs"] = totNumObs
vNum["expGC"] = totExpGC
vNum["expNoGC"] = totExpNoGC

cnt <- labels(vNum)

sc.df = as.data.frame(cbind(cnt, vNum))

# Getting vNum as it needs to be
sc.df$vNum = as.numeric(sc.df$vNum)
sc.df$vNum = round(sc.df$vNum, digits = 1)

# Getting columns in the right order
sc.df$cnt = factor(sc.df$cnt, levels = c("numObs", "expGC", "expNoGC"))

# Making the plot
plot_labs = c("Observed", "Expected (No GC Correction)", "Expected (GC Correction)")

sc_oe <- ggplot(sc.df, aes(x = cnt, y = vNum, fill = cnt)) +
  geom_bar(stat = "identity", position = 'dodge') +
  ggtitle("S. cerevisiae Observed vs. Expected Mutations in tRNA") +
  #geom_text(aes(label = vNum), position = position_dodge(width = 1), vjust = -0.25) +
  labs(x = "",y = "Number of Mutations") +
  scale_x_discrete(labels = plot_labs) +
  scale_fill_manual(values = cbPalette, guide = "none") +
  theme_classic()
sc_oe

############################################################################################################################################
# S. cerevisiae protein-coding, intergenic, and tRNA graph

scRow = as.numeric(rownames(species_tRNA.df[which(species_tRNA.df$species == "Saccharomyces_cerevisiae"),]))
sc_data = species_tRNA.df[scRow,]
sc_data$expLevel = "all"

sc_prot = species_prot.df[which(species_prot.df$species == "Saccharomyces_cerevisiae"),]
sc_inter = species_inter.df[which(species_inter.df$species == "Saccharomyces_cerevisiae"),]

# Combining protein-coding for Figure 1
protObs = sum(sc_prot$observedTot)
protExpGC = sum(sc_prot$expectedTotGC)
protExp = sum(sc_prot$expectedTotNoGC)
protBase = sum(sc_prot$basesCovered)
protMut = protObs/protBase

vProt = c("Saccharomyces_cerevisiae", "protein_coding", "all", protObs, protExpGC, protExp, protBase, protMut)

sc_prot = rbind(sc_prot, vProt)

# Removing strain column
sc_inter = sc_inter[,-9]

# Combining protein-coding and tRNA data
sc_data = rbind(sc_prot, sc_data)
sc_data = rbind(sc_data, sc_inter)

sc_all = sc_data
sc_all = sc_all[which(sc_all$expLevel == "all"),]

sc_data = sc_data[-4,]

# Combining expression bin and gene type into one column
sc_data$expressLvl = "X"
for(i in 1:nrow(sc_data)) {
  sc_data$expressLvl[i] = paste(sc_data$type[i], sc_data$expLevel[i], sep = "_")
}

# Getting the columns in the ideal order
scLevels <- c("tRNA_all", "intergenic_all", "protein_coding_high", "protein_coding_med", "protein_coding_low")
sc_data$expressLvl <- factor(sc_data$expressLvl, levels = scLevels)

sc_type = c("tRNA", "intergenic", "protein_coding")
sc_data$type = factor(sc_data$type, levels = sc_type)

# Making sure mutRate is a numeric value
sc_data$mutRate = as.numeric(sc_data$mutRate)
sc_data$observedTot = as.numeric(sc_data$observedTot)
sc_data$basesCovered = as.numeric(sc_data$basesCovered)

# Getting the standard deviation
sc_data$mutRateSD = 0
for(i in 1:nrow(sc_data)) {
  tmp = prop.test(sc_data$observedTot[i], sc_data$basesCovered[i])
  sd = sc_data$mutRate[i] - tmp$conf.int[1]
  sc_data$mutRateSD[i] = sd
}

# Making the plot
plot_labs = c("All", "All", "High", "Med", "Low")

sc_plot <- ggplot(sc_data, aes(expressLvl, mutRate, fill = type)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = mutRate-mutRateSD, ymax = mutRate+mutRateSD), width = 0.15) +
  labs(title = "S. cerevisiae Mutation Rates by Gene Type", x = "Expression Level", y = "Mutation Rate") +
  scale_x_discrete(labels = plot_labs) +
  scale_fill_manual(values = cbPalette, name = "Gene Type", labels = c("tRNA", "Intergenic", "Protein Coding")) +
  theme_classic()
sc_plot

# Getting the standard deviation
sc_all$mutRate = as.numeric(sc_all$mutRate)
sc_all$observedTot = as.numeric(sc_all$observedTot)
sc_all$basesCovered = as.numeric(sc_all$basesCovered)

sc_all$mutRateSD = 0
for(i in 1:nrow(sc_all)) {
  tmp = prop.test(sc_all$observedTot[i], sc_all$basesCovered[i])
  sd = sc_all$mutRate[i] - tmp$conf.int[1]
  sc_all$mutRateSD[i] = sd
}

sc_all$type <- factor(sc_all$type, levels = sc_type)

# Figure 1
plot_labs = c("tRNA", "Intergenic", "Protein Coding")

sc_plot <- ggplot(sc_all, aes(type, mutRate, fill = type)) +
  geom_bar(stat = "identity", position = "dodge", show.legend = F) +
  geom_errorbar(aes(ymin = mutRate-mutRateSD, ymax = mutRate+mutRateSD), width = 0.15) +
  labs(title = "S. cerevisiae Mutation Rates by Gene Type", x = "Gene Type", y = "Mutation Rate") +
  scale_x_discrete(labels = plot_labs) +
  scale_fill_manual(values = c("#BDD9BF", "#A997DF", "#FFC857", "#E5323B"), name = "Gene Type", 
                    labels = c("tRNA", "Intergenic", "Protein Coding")) +
  theme_classic()
sc_plot

############################################################################################################################################
# S. cerevisiae graph with TPM

# Adding in the tpm split
sc_tpm = species_tRNA_tpm.df[which(species_tRNA_tpm.df$species == "Saccharomyces_cerevisiae"),]

# Getting expressLvl
sc_tpm$expressLvl = "X"
for(i in 1:nrow(sc_tpm)) {
  sc_tpm$expressLvl[i] = paste(sc_tpm$type[i], sc_tpm$expLevel[i], sep = "_")
}

# Getting the standard deviation
sc_tpm$mutRateSD = 0
for(i in 1:nrow(sc_tpm)) {
  tmp = prop.test(sc_tpm$observedTot[i], sc_tpm$basesCovered[i])
  sd = sc_tpm$mutRate[i] - tmp$conf.int[1]
  sc_tpm$mutRateSD[i] = sd
}

# Combining with sc_data
sc_data = rbind(sc_data, sc_tpm)

# Removing tRNA_all data
sc_data = sc_data[which(sc_data$expressLvl != "tRNA_all"),]

# Getting the columns in the ideal order
scLevels <- c("tRNA_high", "tRNA_low", "intergenic_all", "protein_coding_high", "protein_coding_med", "protein_coding_low")
sc_data$expressLvl <- factor(sc_data$expressLvl, levels = scLevels)

# Making the plots
plot_labs = c("High", "Low", "All", "High", "Med", "Low")

sc_plot <- ggplot(sc_data, aes(expressLvl, mutRate, fill = type)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = mutRate-mutRateSD, ymax = mutRate+mutRateSD), width = 0.15) +
  labs(title = "S. cerevisiae Mutation Rates by Gene Type", x = "Expression Level", y = "Mutation Rate") +
  scale_x_discrete(labels = plot_labs) +
  scale_fill_manual(values = cbPalette, name = "Gene Type", labels = c("tRNA", "Intergenic", "Protein Coding")) +
  theme_classic()
sc_plot

############################################################################################################################################
# E. coli protein-coding, intergenic, and tRNA graph

# Getting table with just E. coli data
ecRow = as.numeric(rownames(species_tRNA.df[which(species_tRNA.df$species == "Escherichia_coli_mg1655" | 
                                                    species_tRNA.df$species == "Escherichia_coli_atcc"),]))
ec_data = species_tRNA.df[ecRow,]
ec_data$type = "tRNA"
ec_data$expLevel = "all"

# Separating strain into its own column
ec_data$strain = "X"
for(i in 1:nrow(ec_data)) {
  tmp = str_split(ec_data$species, "_")[i]
  tmp = unlist(tmp)
  
  ec_data$species[i] = "Escherichia_coli"
  ec_data$strain[i] = tmp[3]
}

# Getting tables with just E. coli data
ec_prot = species_prot.df[which(species_prot.df$species == "Escherichia_coli"),]
ec_inter = species_inter.df[which(species_inter.df$species == "Escherichia_coli"),]

# Separating strain into its own column
ec_prot$strain = "X"
for(i in 1:nrow(ec_prot)) {
  tmp = str_split(ec_prot$expLevel, "_")[i]
  tmp = unlist(tmp)
  
  ec_prot$expLevel [i] = tmp[1]
  ec_prot$strain[i] = tmp[2]
}

# Combining protein-coding and tRNA data
ec_data = rbind(ec_prot, ec_data)
ec_data = rbind(ec_data, ec_inter)

ec_levels = c("tRNA", "intergenic", "protein_coding")
ec_data$type = factor(ec_data$type, levels = ec_levels)

# Combining expression bin and gene type into one column
ec_data$expressLvl = "X"
for(i in 1:nrow(ec_data)) {
  ec_data$expressLvl[i] = paste(ec_data$type[i], ec_data$expLevel[i], sep = "_")
}

# Getting the columns in the ideal order
ecLevels <- c("tRNA_all", "intergenic_all", "protein_coding_low", "protein_coding_med", "protein_coding_high")
ec_data$expressLvl <- factor(ec_data$expressLvl, levels = ecLevels)

# Making sure mutRate is a numeric value
ec_data$mutRate = as.numeric(ec_data$mutRate)
ec_data$observedTot = as.numeric(ec_data$observedTot)
ec_data$basesCovered = as.numeric(ec_data$basesCovered)

# Getting the standard deviation
ec_data$mutRateSD = 0
for(i in 1:nrow(ec_data)) {
  tmp = prop.test(ec_data$observedTot[i], ec_data$basesCovered[i])
  sd = ec_data$mutRate[i] - tmp$conf.int[1]
  ec_data$mutRateSD[i] = sd
}

ec_mg1655 = ec_data[which(ec_data$strain == "mg1655"),]
ec_atcc = ec_data[which(ec_data$strain == "atcc"),]

# Making the plots
plot_labs = c("All", "All", "Low", "Med", "High")

mg_plot <- ggplot(ec_mg1655, aes(expressLvl, mutRate, fill = type)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = mutRate-mutRateSD, ymax = mutRate+mutRateSD), width = 0.15) +
  labs(title = "E. coli Mutation Rates by Gene Type", subtitle = "MG1655 Strain", x = "Expression Level", y = "Mutation Rate") +
  scale_x_discrete(labels = plot_labs) +
  scale_fill_manual(values = cbPalette, name = "Gene Type", labels = c("tRNA", "Intergenic", "Protein Coding")) +
  theme_classic()
mg_plot

atcc_plot <- ggplot(ec_atcc, aes(expressLvl, mutRate, fill = type)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = mutRate-mutRateSD, ymax = mutRate+mutRateSD), width = 0.15) +
  labs(title = "E. coli Mutation Rates by Gene Type", subtitle = "ATCC 8739 Strain", x = "Expression Level", y = "Mutation Rate") +
  scale_x_discrete(labels = plot_labs) +
  scale_fill_manual(values = cbPalette, name = "Gene Type", labels = c("tRNA", "Intergenic", "Protein Coding")) +
  theme_classic()
atcc_plot

############################################################################################################################################
# E. coli graphs with TPM

# Adding in the tpm split
ec_tpm = species_tRNA_tpm.df[which(species_tRNA_tpm.df$species == "Escherichia_coli"),]

# Separating expLevel and strain
ec_tpm$strain = "X"
for(i in 1:nrow(ec_tpm)) {
  tmp = str_split(ec_tpm$expLevel, "_")[i]
  tmp = unlist(tmp)
  
  ec_tpm$expLevel[i] = tmp[1]
  ec_tpm$strain[i] = tmp[2]
}

# Getting expressLvl
ec_tpm$expressLvl = "X"
for(i in 1:nrow(ec_tpm)) {
  ec_tpm$expressLvl[i] = paste(ec_tpm$type[i], ec_tpm$expLevel[i], sep = "_")
}

# Getting the standard deviation
ec_tpm$mutRateSD = 0
for(i in 1:nrow(ec_tpm)) {
  tmp = prop.test(ec_tpm$observedTot[i], ec_tpm$basesCovered[i])
  sd = ec_tpm$mutRate[i] - tmp$conf.int[1]
  ec_tpm$mutRateSD[i] = sd
}

# Combining with ec_data
ec_data = rbind(ec_data, ec_tpm)

# Removing tRNA_all data
ec_data = ec_data[which(ec_data$expressLvl != "tRNA_all"),]

# Getting the columns in the ideal order
ecLevels <- c("tRNA_low", "tRNA_high", "intergenic_all", "protein_coding_low", "protein_coding_med", "protein_coding_high")
ec_data$expressLvl <- factor(ec_data$expressLvl, levels = ecLevels)

# Separating strains
ec_mg1655 = ec_data[which(ec_data$strain == "mg1655"),]
ec_atcc = ec_data[which(ec_data$strain == "atcc"),]

# Making the plots
plot_labs = c("Low", "High", "All", "Low", "Med", "High")

mg_plot <- ggplot(ec_mg1655, aes(expressLvl, mutRate, fill = type)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = mutRate-mutRateSD, ymax = mutRate+mutRateSD), width = 0.15) +
  labs(title = "E. coli Mutation Rates by Gene Type", subtitle = "MG1655 Strain", x = "Expression Level", y = "Mutation Rate") +
  scale_x_discrete(labels = plot_labs) +
  scale_fill_manual(values = cbPalette, name = "Gene Type", labels = c("tRNA", "Intergenic", "Protein Coding")) +
  theme_classic()
mg_plot

atcc_plot <- ggplot(ec_atcc, aes(expressLvl, mutRate, fill = type)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = mutRate-mutRateSD, ymax = mutRate+mutRateSD), width = 0.15) +
  labs(title = "E. coli Mutation Rates by Gene Type", subtitle = "ATCC 8739 Strain", x = "Expression Level", y = "Mutation Rate") +
  scale_x_discrete(labels = plot_labs) +
  scale_fill_manual(values = cbPalette, name = "Gene Type", labels = c("tRNA", "Intergenic", "Protein Coding")) +
  theme_classic()
atcc_plot

#################################################################################

# Combined dataset
ec_full = ec_data[,-c(8,9,11)] # Removing mutation rate, error bars, and strain

ec_full$observedTot = as.numeric(ec_full$observedTot)
ec_full$expectedTotGC = as.numeric(ec_full$expectedTotGC)
ec_full$expectedTotNoGC = as.numeric(ec_full$expectedTotNoGC)
ec_full$basesCovered = as.numeric(ec_full$basesCovered)

vExp = unique(ec_full$expressLvl)

fullMatrix = matrix(nrow = 0, ncol = ncol(ec_full))
colnames(fullMatrix) = colnames(ec_full)

for(lvl in vExp) {
  vNew = c()
  tmp.df = ec_full[which(ec_full$expressLvl == lvl),]
  
  totObs = 0
  totExpGC = 0
  totExpNoGC = 0
  totBases = 0
  for(i in 1:nrow(tmp.df)) {
    totObs = totObs + tmp.df$observedTot[i]
    totExpGC = totExpGC + tmp.df$expectedTotGC[i]
    totExpNoGC = totExpNoGC + tmp.df$expectedTotNoGC[i]
    totBases = totBases + tmp.df$basesCovered[i]
  }
  
  vNew[1] <- tmp.df$species[1]
  vNew[2] <- tmp.df$type[1]
  vNew[3] <- tmp.df$expLevel[1]
  vNew[4] <- totObs
  vNew[5] <- totExpGC
  vNew[6] <- totExpNoGC
  vNew[7] <- totBases
  vNew[8] <- lvl
  
  fullMatrix = rbind(fullMatrix, vNew)
}

ec_full = as.data.frame(fullMatrix, row.names = c(1:nrow(fullMatrix)))

ec_full$observedTot = as.numeric(ec_full$observedTot)
ec_full$basesCovered = as.numeric(ec_full$basesCovered)
ec_full$mutRate = ec_full$observedTot/ec_full$basesCovered

# Getting the standard deviation
ec_full$mutRateSD = 0
for(i in 1:nrow(ec_full)) {
  tmp = prop.test(ec_full$observedTot[i], ec_full$basesCovered[i])
  sd = ec_full$mutRate[i] - tmp$conf.int[1]
  ec_full$mutRateSD[i] = sd
}

ecLevels <- c("tRNA_low", "tRNA_high", "intergenic_all", "protein_coding_low", "protein_coding_med", "protein_coding_high")
ec_full$expressLvl <- factor(ec_full$expressLvl, levels = ecLevels)

# Combined plot
ec_plot <- ggplot(ec_full, aes(expressLvl, mutRate, fill = type)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = mutRate-mutRateSD, ymax = mutRate+mutRateSD), width = 0.15) +
  labs(title = "E. coli Mutation Rates by Gene Type", subtitle = "MG1655 & ATCC 8739 Strains", x = "Expression Level", y = "Mutation Rate") +
  scale_x_discrete(labels = plot_labs) +
  scale_fill_manual(values = cbPalette, name = "Gene Type", labels = c("tRNA", "Intergenic", "Protein Coding")) +
  theme_classic()
ec_plot


############################################################################################################################################
# BPS spectra plot
tRNA_sp = c("Saccharomyces_cerevisiae", "Escherichia_coli_mg1655", "Escherichia_coli_atcc")

# Reducing to only the species we're interested in
tRNA_bpsRed.df = tRNA_bps.df[which(tRNA_bps.df$species %in% tRNA_sp),]
rest_bpsRed.df = rest_tRNA_bps.df[which(rest_tRNA_bps.df$species %in% tRNA_sp),]

ec_AT_GC_tRNA = 0
ec_GC_AT_tRNA = 0
ec_AT_TA_tRNA = 0
ec_AT_CG_tRNA = 0
ec_GC_TA_tRNA = 0
ec_GC_CG_tRNA = 0

ec_AT_GC_rest = 0
ec_GC_AT_rest = 0
ec_AT_TA_rest = 0
ec_AT_CG_rest = 0
ec_GC_TA_rest = 0
ec_GC_CG_rest = 0

ec_strain = c("Escherichia_coli_mg1655", "Escherichia_coli_atcc")

# Combining tRNA and rest data from both strains
for(strain in ec_strain) {
  ec_AT_GC_tRNA = ec_AT_GC_tRNA + tRNA_bpsRed.df$`A:T>G:C`[which(tRNA_bpsRed.df$species == strain)]
  ec_GC_AT_tRNA = ec_GC_AT_tRNA + tRNA_bpsRed.df$`G:C>A:T`[which(tRNA_bpsRed.df$species == strain)]
  ec_AT_TA_tRNA = ec_AT_TA_tRNA + tRNA_bpsRed.df$`A:T>T:A`[which(tRNA_bpsRed.df$species == strain)]
  ec_AT_CG_tRNA = ec_AT_CG_tRNA + tRNA_bpsRed.df$`A:T>C:G`[which(tRNA_bpsRed.df$species == strain)]
  ec_GC_TA_tRNA = ec_GC_TA_tRNA + tRNA_bpsRed.df$`G:C>T:A`[which(tRNA_bpsRed.df$species == strain)]
  ec_GC_CG_tRNA = ec_GC_CG_tRNA + tRNA_bpsRed.df$`G:C>C:G`[which(tRNA_bpsRed.df$species == strain)]
  
  ec_AT_GC_rest = ec_AT_GC_rest + rest_bpsRed.df$`A:T>G:C`[which(rest_bpsRed.df$species == strain)]
  ec_GC_AT_rest = ec_GC_AT_rest + rest_bpsRed.df$`G:C>A:T`[which(rest_bpsRed.df$species == strain)]
  ec_AT_TA_rest = ec_AT_TA_rest + rest_bpsRed.df$`A:T>T:A`[which(rest_bpsRed.df$species == strain)]
  ec_AT_CG_rest = ec_AT_CG_rest + rest_bpsRed.df$`A:T>C:G`[which(rest_bpsRed.df$species == strain)]
  ec_GC_TA_rest = ec_GC_TA_rest + rest_bpsRed.df$`G:C>T:A`[which(rest_bpsRed.df$species == strain)]
  ec_GC_CG_rest = ec_GC_CG_rest + rest_bpsRed.df$`G:C>C:G`[which(rest_bpsRed.df$species == strain)]
}

ec_comb_tRNA = c("Escherichia_coli", ec_AT_GC_tRNA, ec_GC_AT_tRNA, ec_AT_TA_tRNA, ec_AT_CG_tRNA, ec_GC_TA_tRNA, ec_GC_CG_tRNA)
ec_comb_rest = c("Escherichia_coli", ec_AT_GC_rest, ec_GC_AT_rest, ec_AT_TA_rest, ec_AT_CG_rest, ec_GC_TA_rest, ec_GC_CG_rest)

# Adding combined data to the table
tRNA_bpsRed.df = rbind(tRNA_bpsRed.df, ec_comb_tRNA)
rest_bpsRed.df = rbind(rest_bpsRed.df, ec_comb_rest)

# Adding type indicator
tRNA_bpsRed.df$type = "tRNA"
rest_bpsRed.df$type = "rest"

# Reducing to just the strain-combined data
spNew = c("Saccharomyces_cerevisiae", "Escherichia_coli")
tRNA_bpsRed.df = tRNA_bpsRed.df[which(tRNA_bpsRed.df$species %in% spNew),]
rest_bpsRed.df = rest_bpsRed.df[which(rest_bpsRed.df$species %in% spNew),]

# Making a table to hold the combined data
BPS_row = nrow(tRNA_bpsRed.df) * 2
BPS_full.df = data.frame(species = character(BPS_row), v1 = numeric(BPS_row), v2 = numeric(BPS_row), v3 = numeric(BPS_row), 
                         v4 = numeric(BPS_row), v5 = numeric(BPS_row), v6 = numeric(BPS_row), type = character(BPS_row))

colnames(BPS_full.df) = c("species", "A:T>G:C", "G:C>A:T", "A:T>T:A", "A:T>C:G", "G:C>T:A", "G:C>C:G", "type")

# Combining tRNA table with rest table
rowCnt = 1
for(i in 1:BPS_row) {
  if(BPS_row/i >= 2) {
    BPS_full.df[i,] = tRNA_bpsRed.df[i,]
  }
  
  else{
    BPS_full.df[i,] = rest_bpsRed.df[rowCnt,]
    rowCnt = rowCnt + 1
  }
}

df <- gather(BPS_full.df, BPS, mutNum, 2:7)
df$mutNum = as.numeric(df$mutNum)
df$fullType = paste(df$species, df$type, sep = "_")

typeList = unique(df$fullType)
typeNum = length(typeList)

# Adding the "total" column
df$total = 0
for(i in 1:typeNum) {
  tmpType = typeList[i]
  tmpSum = sum(df$mutNum[which(df$fullType == tmpType)])
  df$total[which(df$fullType == tmpType)] = tmpSum
}

# Changing type to something more appropriate for the graph
df$type[which(df$type == "rest")] = "Protein-Coding Genes & Intergenic Regions"

# Getting the columns to display in the right order
typeLvls = c("tRNA", "Protein-Coding Genes & Intergenic Regions")
df$type = factor(df$type, levels = typeLvls)

# Getting the BPSs to display in the right order
BPSlvls = c("A:T>G:C", "G:C>A:T", "A:T>T:A", "A:T>C:G", "G:C>T:A", "G:C>C:G")
df$BPS = factor(df$BPS, levels = BPSlvls)

# Getting the standard deviation
df$mutProp = df$mutNum / df$total
df$mutPropSD = 0
for(i in 1:nrow(df)) {
  tmp = prop.test(df$mutNum[i], df$total[i])
  sd = df$mutProp[i] - tmp$conf.int[1]
  df$mutPropSD[i] = sd
}

# Making the plot
BPS_plot <- ggplot(data = df, aes(x = BPS, y = mutNum/total, fill = species, pattern = type)) +
  geom_bar_pattern(stat = "identity", position = position_dodge(preserve = "single"),
                   color = "black", 
                   pattern_fill = "black",
                   pattern_angle = 45,
                   pattern_density = 0.1,
                   pattern_spacing = 0.015,
                   pattern_key_scale_factor = 0.6) + 
  geom_errorbar(aes(ymin = mutProp-mutPropSD, ymax = mutProp+mutPropSD), position = "dodge") +
  #stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.8, position = position_dodge(0.9)) +
  scale_fill_manual(values = colorRampPalette(bpsPalette)(2), name = "Species", labels = c("E. coli", "S. cerevisiae")) +
  scale_pattern_manual(values = c(tRNA = "stripe", `Protein-Coding Genes & Intergenic Regions` = "none")) +
  labs(title = "Base Pair Substitution Frequencies by Species", x = "BPS", y = "Frequency", pattern = "Type") + 
  guides(pattern = guide_legend(override.aes = list(fill = "white")),
         fill = guide_legend(override.aes = list(pattern = "none"))) +
  theme(legend.title = element_text(size=20), legend.text = element_text(size=20)) +
  theme_classic()

BPS_plot

############################################################################################################################################
# Mutation spectra plot - tRNA
tRNA_sp = c("Saccharomyces_cerevisiae", "Escherichia_coli_mg1655", "Escherichia_coli_atcc")

# Reducing to only the species we're interested in
tRNA_mutRed.df = tRNA_mut.df[which(tRNA_mut.df$species %in% tRNA_sp),]
rest_mutRed.df = rest_tRNA_mut.df[which(rest_tRNA_mut.df$species %in% tRNA_sp),]

ec_AT_tRNA = 0
ec_AG_tRNA = 0
ec_AC_tRNA = 0
ec_TA_tRNA = 0
ec_TG_tRNA = 0
ec_TC_tRNA = 0
ec_GA_tRNA = 0
ec_GT_tRNA = 0
ec_GC_tRNA = 0
ec_CA_tRNA = 0
ec_CT_tRNA = 0
ec_CG_tRNA = 0

ec_AT_rest = 0
ec_AG_rest = 0
ec_AC_rest = 0
ec_TA_rest = 0
ec_TG_rest = 0
ec_TC_rest = 0
ec_GA_rest = 0
ec_GT_rest = 0
ec_GC_rest = 0
ec_CA_rest = 0
ec_CT_rest = 0
ec_CG_rest = 0

ec_strain = c("Escherichia_coli_mg1655", "Escherichia_coli_atcc")

# Combining tRNA and rest data from both strains
for(strain in ec_strain) {
  ec_AT_tRNA = ec_AT_tRNA + tRNA_mutRed.df$`A>T`[which(tRNA_mutRed.df$species == strain)]
  ec_AG_tRNA = ec_AG_tRNA + tRNA_mutRed.df$`A>G`[which(tRNA_mutRed.df$species == strain)]
  ec_AC_tRNA = ec_AC_tRNA + tRNA_mutRed.df$`A>C`[which(tRNA_mutRed.df$species == strain)]
  ec_TA_tRNA = ec_TA_tRNA + tRNA_mutRed.df$`T>A`[which(tRNA_mutRed.df$species == strain)]
  ec_TG_tRNA = ec_TG_tRNA + tRNA_mutRed.df$`T>G`[which(tRNA_mutRed.df$species == strain)]
  ec_TC_tRNA = ec_TC_tRNA + tRNA_mutRed.df$`T>C`[which(tRNA_mutRed.df$species == strain)]
  ec_GA_tRNA = ec_GA_tRNA + tRNA_mutRed.df$`G>A`[which(tRNA_mutRed.df$species == strain)]
  ec_GT_tRNA = ec_GT_tRNA + tRNA_mutRed.df$`G>T`[which(tRNA_mutRed.df$species == strain)]
  ec_GC_tRNA = ec_GC_tRNA + tRNA_mutRed.df$`G>C`[which(tRNA_mutRed.df$species == strain)]
  ec_CA_tRNA = ec_CA_tRNA + tRNA_mutRed.df$`C>A`[which(tRNA_mutRed.df$species == strain)]
  ec_CT_tRNA = ec_CT_tRNA + tRNA_mutRed.df$`C>T`[which(tRNA_mutRed.df$species == strain)]
  ec_CG_tRNA = ec_CG_tRNA + tRNA_mutRed.df$`C>G`[which(tRNA_mutRed.df$species == strain)]
  
  ec_AT_rest = ec_AT_rest + rest_mutRed.df$`A>T`[which(rest_mutRed.df$species == strain)]
  ec_AG_rest = ec_AG_rest + rest_mutRed.df$`A>G`[which(rest_mutRed.df$species == strain)]
  ec_AC_rest = ec_AC_rest + rest_mutRed.df$`A>C`[which(rest_mutRed.df$species == strain)]
  ec_TA_rest = ec_TA_rest + rest_mutRed.df$`T>A`[which(rest_mutRed.df$species == strain)]
  ec_TG_rest = ec_TG_rest + rest_mutRed.df$`T>G`[which(rest_mutRed.df$species == strain)]
  ec_TC_rest = ec_TC_rest + rest_mutRed.df$`T>C`[which(rest_mutRed.df$species == strain)]
  ec_GA_rest = ec_GA_rest + rest_mutRed.df$`G>A`[which(rest_mutRed.df$species == strain)]
  ec_GT_rest = ec_GT_rest + rest_mutRed.df$`G>T`[which(rest_mutRed.df$species == strain)]
  ec_GC_rest = ec_GC_rest + rest_mutRed.df$`G>C`[which(rest_mutRed.df$species == strain)]
  ec_CA_rest = ec_CA_rest + rest_mutRed.df$`C>A`[which(rest_mutRed.df$species == strain)]
  ec_CT_rest = ec_CT_rest + rest_mutRed.df$`C>T`[which(rest_mutRed.df$species == strain)]
  ec_CG_rest = ec_CG_rest + rest_mutRed.df$`C>G`[which(rest_mutRed.df$species == strain)]
}

ec_comb_tRNA = c("Escherichia_coli", ec_AT_tRNA, ec_AG_tRNA, ec_AC_tRNA, ec_TA_tRNA, ec_TG_tRNA, ec_TC_tRNA, ec_GA_tRNA, ec_GT_tRNA, 
                 ec_GC_tRNA, ec_CA_tRNA, ec_CT_tRNA, ec_CG_tRNA)
ec_comb_rest = c("Escherichia_coli", ec_AT_rest, ec_AG_rest, ec_AC_rest, ec_TA_rest, ec_TG_rest, ec_TC_rest, ec_GA_rest, ec_GT_rest, 
                 ec_GC_rest, ec_CA_rest, ec_CT_rest, ec_CG_rest)

# Adding combined data to the tables
tRNA_mutRed.df = rbind(tRNA_mutRed.df, ec_comb_tRNA)
rest_mutRed.df = rbind(rest_mutRed.df, ec_comb_rest)

# Reading in the table of nucleotide numbers in tRNA genes and the rest of the genome
nt_III.df = read.table("ntFreq_tRNA.tab", h = T)
nt_III.df = nt_III.df[which(nt_III.df$species %in% tRNA_sp),]

nt_rest.df = read.table("ntFreq_rest_tRNA.tab", h = T)
nt_rest.df = nt_rest.df[which(nt_rest.df$species %in% tRNA_sp),]

ec_A_III = 0
ec_T_III = 0
ec_G_III = 0
ec_C_III = 0

ec_A_rest = 0
ec_T_rest = 0
ec_G_rest = 0
ec_C_rest = 0

# Making sure numbers are numerics in table
nt_III.df$freq_A = as.numeric(nt_III.df$freq_A)
nt_III.df$freq_T = as.numeric(nt_III.df$freq_T)
nt_III.df$freq_G = as.numeric(nt_III.df$freq_G)
nt_III.df$freq_C = as.numeric(nt_III.df$freq_C)

nt_rest.df$freq_A = as.numeric(nt_rest.df$freq_A)
nt_rest.df$freq_T = as.numeric(nt_rest.df$freq_T)
nt_rest.df$freq_G = as.numeric(nt_rest.df$freq_G)
nt_rest.df$freq_C = as.numeric(nt_rest.df$freq_C)

for(strain in ec_strain) {
  ec_A_III = ec_A_III + nt_III.df$freq_A[which(nt_III.df$species == strain)]
  ec_T_III = ec_T_III + nt_III.df$freq_T[which(nt_III.df$species == strain)]
  ec_G_III = ec_G_III + nt_III.df$freq_G[which(nt_III.df$species == strain)]
  ec_C_III = ec_C_III + nt_III.df$freq_C[which(nt_III.df$species == strain)]
  
  ec_A_rest = ec_A_rest + nt_rest.df$freq_A[which(nt_rest.df$species == strain)]
  ec_T_rest = ec_T_rest + nt_rest.df$freq_T[which(nt_rest.df$species == strain)]
  ec_G_rest = ec_G_rest + nt_rest.df$freq_G[which(nt_rest.df$species == strain)]
  ec_C_rest = ec_C_rest + nt_rest.df$freq_C[which(nt_rest.df$species == strain)]
}

ec_nt_III = c("Escherichia_coli", ec_A_III, ec_T_III, ec_G_III, ec_C_III)
ec_nt_rest = c("Escherichia_coli", ec_A_rest, ec_T_rest, ec_G_rest, ec_C_rest)

nt_III.df = rbind(nt_III.df, ec_nt_III)
nt_rest.df = rbind(nt_rest.df, ec_nt_rest)

# Reducing to just combined data
spNew = c("Saccharomyces_cerevisiae", "Escherichia_coli")

tRNA_mutRed.df = tRNA_mutRed.df[which(tRNA_mutRed.df$species %in% spNew),]
rest_mutRed.df = rest_mutRed.df[which(rest_mutRed.df$species %in% spNew),]
nt_III.df = nt_III.df[which(nt_III.df$species %in% spNew),]
nt_rest.df = nt_rest.df[which(nt_rest.df$species %in% spNew),]

##########################################################################################
# Getting mutations/total # nucleotide

# Making copies of the tables for ntMutPlot graphs
tRNA_nt.df = tRNA_mutRed.df
rest_nt.df = rest_mutRed.df

# Columns by ref base
Amut = c("A>T", "A>G", "A>C")
Tmut = c("T>A", "T>G", "T>C")
Gmut = c("G>A", "G>T", "G>C")
Cmut = c("C>A", "C>T", "C>G")

# Dividing number of mutations by number of reference base occurrences in genes of interest
# Replacing number of mutations in table with new number calculated
for(sp in spNew) {
  fA_III = as.numeric(nt_III.df$freq_A[which(nt_III.df$species == sp)])
  fT_III = as.numeric(nt_III.df$freq_T[which(nt_III.df$species == sp)])
  fG_III = as.numeric(nt_III.df$freq_G[which(nt_III.df$species == sp)])
  fC_III = as.numeric(nt_III.df$freq_C[which(nt_III.df$species == sp)])
  
  fA_rest = as.numeric(nt_rest.df$freq_A[which(nt_rest.df$species == sp)])
  fT_rest = as.numeric(nt_rest.df$freq_T[which(nt_rest.df$species == sp)])
  fG_rest = as.numeric(nt_rest.df$freq_G[which(nt_rest.df$species == sp)])
  fC_rest = as.numeric(nt_rest.df$freq_C[which(nt_rest.df$species == sp)])
  
  for(a in Amut) {
    tmp = as.numeric(tRNA_nt.df[which(tRNA_nt.df$species == sp), a])
    newA = tmp/fA_III
    tRNA_nt.df[which(tRNA_nt.df$species == sp), a] = newA
    
    tmp = as.numeric(rest_nt.df[which(rest_nt.df$species == sp), a])
    newA = tmp/fA_rest
    rest_nt.df[which(rest_nt.df$species == sp), a] = newA
  }
  
  for(t in Tmut) {
    tmp = as.numeric(tRNA_nt.df[which(tRNA_nt.df$species == sp), t])
    newT = tmp/fT_III
    tRNA_nt.df[which(tRNA_nt.df$species == sp), t] = newT
    
    tmp = as.numeric(rest_nt.df[which(rest_nt.df$species == sp), t])
    newT = tmp/fT_rest
    rest_nt.df[which(rest_nt.df$species == sp), t] = newT
  }
  
  for(g in Gmut) {
    tmp = as.numeric(tRNA_nt.df[which(tRNA_nt.df$species == sp), g])
    newG = tmp/fG_III
    tRNA_nt.df[which(tRNA_nt.df$species == sp), g] = newG
    
    tmp = as.numeric(rest_nt.df[which(rest_nt.df$species == sp), g])
    newG = tmp/fG_rest
    rest_nt.df[which(rest_nt.df$species == sp), g] = newG
  }
  
  for(c in Cmut) {
    tmp = as.numeric(tRNA_nt.df[which(tRNA_nt.df$species == sp), c])
    newC = tmp/fC_III
    tRNA_nt.df[which(tRNA_nt.df$species == sp), c] = newC
    
    tmp = as.numeric(rest_nt.df[which(rest_nt.df$species == sp), c])
    newC = tmp/fC_rest
    rest_nt.df[which(rest_nt.df$species == sp), c] = newC
  }
}


# Making a table to hold the combined data
tRNA_mutRed.df$type = "tRNA"
rest_mutRed.df$type = "rest"

tRNA_nt.df$type = "tRNA"
rest_nt.df$type = "rest"

mut_row = nrow(tRNA_mutRed.df) * 2
mut_full.df = data.frame(species = character(mut_row), v1 = numeric(mut_row), v2 = numeric(mut_row), v3 = numeric(mut_row), 
                         v4 = numeric(mut_row), v5 = numeric(mut_row), v6 = numeric(mut_row), v7 = numeric(mut_row), v8 = numeric(mut_row),
                         v9 = numeric(mut_row), v10 = numeric(mut_row), v11 = numeric(mut_row), v12 = numeric(mut_row),
                         type = character(mut_row))

mut_full_nt.df = data.frame(species = character(mut_row), v1 = numeric(mut_row), v2 = numeric(mut_row), v3 = numeric(mut_row), 
                         v4 = numeric(mut_row), v5 = numeric(mut_row), v6 = numeric(mut_row), v7 = numeric(mut_row), v8 = numeric(mut_row),
                         v9 = numeric(mut_row), v10 = numeric(mut_row), v11 = numeric(mut_row), v12 = numeric(mut_row),
                         type = character(mut_row))

colnames(mut_full.df) = c("species", "A>T", "A>G", "A>C", "T>A", "T>G", "T>C", "G>A", "G>T", "G>C", "C>A", "C>T", "C>G", "type")
colnames(mut_full_nt.df) = c("species", "A>T", "A>G", "A>C", "T>A", "T>G", "T>C", "G>A", "G>T", "G>C", "C>A", "C>T", "C>G", "type")

# Combining tRNA table with rest table
rowCnt = 1
for(i in 1:mut_row) {
  if(mut_row/i >= 2) {
    mut_full.df[i,] = tRNA_mutRed.df[i,]
    mut_full_nt.df[i,] = tRNA_nt.df[i,]
  }
  
  else{
    mut_full.df[i,] = rest_mutRed.df[rowCnt,]
    mut_full_nt.df[i,] = rest_nt.df[rowCnt,]
    rowCnt = rowCnt + 1
  }
}

df <- gather(mut_full.df, mut, mutNum, 2:13)
df$mutNum = as.numeric(df$mutNum)
df$fullType = paste(df$species, df$type, sep = "_")

df2 <- gather(mut_full_nt.df, mut, mutNum, 2:13)
df2$mutNum = as.numeric(df2$mutNum)
df2$fullType = paste(df2$species, df2$type, sep = "_")

typeList = unique(df$fullType)
typeNum = length(typeList)

# Adding the "total" column
df$total = 0
df2$total = 0
for(i in 1:typeNum) {
  tmpType = typeList[i]
  tmpSum = sum(df$mutNum[which(df$fullType == tmpType)])
  df$total[which(df$fullType == tmpType)] = tmpSum
  
  tmpSum2 = sum(df2$mutNum[which(df2$fullType == tmpType)])
  df2$total[which(df2$fullType == tmpType)] = tmpSum2
}

# Changing type to something more appropriate for the graph
df$type[which(df$type == "rest")] = "Protein-Coding Genes & Intergenic Regions"
df2$type[which(df2$type == "rest")] = "Protein-Coding Genes & Intergenic Regions"

# Getting the columns to display in the right order
typeLvls = c("tRNA", "Protein-Coding Genes & Intergenic Regions")
df$type = factor(df$type, levels = typeLvls)
df2$type = factor(df2$type, levels = typeLvls)

# Getting the BPSs to display in the right order
mutlvls = c("A>T", "A>G", "A>C", "T>A", "T>G", "T>C", "G>A", "G>T", "G>C", "C>A", "C>T", "C>G")
df$BPS = factor(df$mut, levels = mutlvls)
df2$BPS = factor(df2$mut, levels = mutlvls)

# Getting the standard deviation
df$mutProp = df$mutNum / df$total
df$mutPropSD = 0

df2$mutProp = df2$mutNum / df2$total
df2$mutPropSD = 0

for(i in 1:nrow(df)) {
  tmp = prop.test(df$mutNum[i], df$total[i])
  sd = df$mutProp[i] - tmp$conf.int[1]
  df$mutPropSD[i] = sd
  
  tmp2 = prop.test(df2$mutNum[i], df2$total[i])
  sd2 = df2$mutProp[i] - tmp2$conf.int[1]
  df2$mutPropSD[i] = sd2
}

################################################################################################
# Plots

# Number of specific mutation divided by total number of mutations
tRNAMutPlot <- ggplot(data = df, aes(x = mut, y = mutNum/total, fill = species, pattern = type)) +
  geom_bar_pattern(stat = "identity", position = position_dodge(preserve = "single"),
                   color = "black", 
                   pattern_fill = "black",
                   pattern_angle = 45,
                   pattern_density = 0.1,
                   pattern_spacing = 0.015,
                   pattern_key_scale_factor = 0.6) + 
  geom_errorbar(aes(ymin = mutProp-mutPropSD, ymax = mutProp+mutPropSD), position = "dodge") +
  #stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.8, position = position_dodge(0.9)) +
  scale_fill_manual(values = colorRampPalette(bpsPalette)(2), name = "Species", labels = c("E. coli", "S. cerevisiae")) +
  scale_pattern_manual(values = c(tRNA = "stripe", `Protein-Coding Genes & Intergenic Regions` = "none")) +
  labs(title = "Mutation Frequencies by Species", x = "Mutation", y = "Frequency", pattern = "Type") + 
  guides(pattern = guide_legend(override.aes = list(fill = "white")),
         fill = guide_legend(override.aes = list(pattern = "none"))) +
  theme(legend.title = element_text(size=20), legend.text = element_text(size=20)) +
  theme_classic()

tRNAMutPlot    # <-- Saved as mutFreqTot.png locally


# Plot of number of mutated nucleotide divided by total number of same nucleotide (# A muts / # As)
# Error bars are excluded because they are not informative
tRNAntMutPlot <- ggplot(data = df2, aes(x = mut, y = mutNum, fill = species, pattern = type)) +
  geom_bar_pattern(stat = "identity", position = position_dodge(preserve = "single"),
                   color = "black", 
                   pattern_fill = "black",
                   pattern_angle = 45,
                   pattern_density = 0.1,
                   pattern_spacing = 0.015,
                   pattern_key_scale_factor = 0.6) + 
  #geom_errorbar(aes(ymin = mutProp-mutPropSD, ymax = mutProp+mutPropSD), position = "dodge") +
  #stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.8, position = position_dodge(0.9)) +
  scale_fill_manual(values = colorRampPalette(bpsPalette)(2), name = "Species", labels = c("E. coli", "S. cerevisiae")) +
  scale_pattern_manual(values = c(tRNA = "stripe", `Protein-Coding Genes & Intergenic Regions` = "none")) +
  labs(title = "Mutation Frequencies by Species", x = "Mutation", y = "Frequency", pattern = "Type") + 
  guides(pattern = guide_legend(override.aes = list(fill = "white")),
         fill = guide_legend(override.aes = list(pattern = "none"))) +
  theme(legend.title = element_text(size=20), legend.text = element_text(size=20)) +
  theme_classic()

tRNAntMutPlot


# ntMutPlot but only E. coli data
df_ec = df2[which(df2$species == "Escherichia_coli"),]
tRNAecMutPlot <- ggplot(data = df_ec, aes(x = mut, y = mutNum, fill = species, pattern = type)) +
  geom_bar_pattern(stat = "identity", position = position_dodge(preserve = "single"),
                   color = "black", 
                   pattern_fill = "black",
                   pattern_angle = 45,
                   pattern_density = 0.1,
                   pattern_spacing = 0.015,
                   pattern_key_scale_factor = 0.6) + 
  #geom_errorbar(aes(ymin = mutProp-mutPropSD, ymax = mutProp+mutPropSD), position = "dodge") +
  #stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.8, position = position_dodge(0.9)) +
  scale_fill_manual(values = colorRampPalette(bpsPalette)(2), name = "Species", labels = c("E. coli", "S. cerevisiae")) +
  scale_pattern_manual(values = c(tRNA = "stripe", `Protein-Coding Genes & Intergenic Regions` = "none")) +
  labs(title = "Mutation Frequencies by Species", x = "Mutation", y = "Frequency", pattern = "Type") + 
  guides(pattern = guide_legend(override.aes = list(fill = "white")),
         fill = guide_legend(override.aes = list(pattern = "none"))) +
  theme(legend.title = element_text(size=20), legend.text = element_text(size=20)) +
  theme_classic()

tRNAecMutPlot


# ntMutPlot but only S. cerevisiae data
df_sc = df2[which(df2$species == "Saccharomyces_cerevisiae"),]
tRNAscMutPlot <- ggplot(data = df_sc, aes(x = mut, y = mutNum, fill = species, pattern = type)) +
  geom_bar_pattern(stat = "identity", position = position_dodge(preserve = "single"),
                   color = "black", 
                   pattern_fill = "black",
                   pattern_angle = 45,
                   pattern_density = 0.1,
                   pattern_spacing = 0.015,
                   pattern_key_scale_factor = 0.6) + 
  #geom_errorbar(aes(ymin = mutProp-mutPropSD, ymax = mutProp+mutPropSD), position = "dodge") +
  #stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.8, position = position_dodge(0.9)) +
  scale_fill_manual(values = c("#E5323B"), name = "Species", labels = c("S. cerevisiae")) +
  scale_pattern_manual(values = c(tRNA = "stripe", `Protein-Coding Genes & Intergenic Regions` = "none")) +
  labs(title = "Mutation Frequencies by Species", x = "Mutation", y = "Frequency", pattern = "Type") + 
  guides(pattern = guide_legend(override.aes = list(fill = "white")),
         fill = guide_legend(override.aes = list(pattern = "none"))) +
  theme(legend.title = element_text(size=20), legend.text = element_text(size=20)) +
  theme_classic()

tRNAscMutPlot


# ntMutPlot but normalized so that frequencies = 1
tRNAntMutPlot2 <- ggplot(data = df, aes(x = mut, y = mutNum/total, fill = species, pattern = type)) +
  geom_bar_pattern(stat = "identity", position = position_dodge(preserve = "single"),
                   color = "black", 
                   pattern_fill = "black",
                   pattern_angle = 45,
                   pattern_density = 0.1,
                   pattern_spacing = 0.015,
                   pattern_key_scale_factor = 0.6) + 
  #geom_errorbar(aes(ymin = mutProp-mutPropSD, ymax = mutProp+mutPropSD), position = "dodge") +
  #stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.8, position = position_dodge(0.9)) +
  scale_fill_manual(values = colorRampPalette(bpsPalette)(2), name = "Species", labels = c("E. coli", "S. cerevisiae")) +
  scale_pattern_manual(values = c(tRNA = "stripe", `Protein-Coding Genes & Intergenic Regions` = "none")) +
  labs(title = "Mutation Frequencies by Species", x = "Mutation", y = "Frequency", pattern = "Type") + 
  guides(pattern = guide_legend(override.aes = list(fill = "white")),
         fill = guide_legend(override.aes = list(pattern = "none"))) +
  theme(legend.title = element_text(size=20), legend.text = element_text(size=20)) +
  theme_classic()

tRNAntMutPlot2





############################################################################################################################################
# Mutation spectra plot - protein-coding
prot_sp = c("Saccharomyces_cerevisiae_low", "Saccharomyces_cerevisiae_med", "Saccharomyces_cerevisiae_high",
            "Escherichia_coli_low_mg1655", "Escherichia_coli_med_mg1655", "Escherichia_coli_high_mg1655",
            "Escherichia_coli_low_atcc", "Escherichia_coli_med_atcc", "Escherichia_coli_high_atcc")

# Removing the _prot string from species names in prot_mut.df & rest_prot_mut.df
for(i in 1:nrow(prot_mut.df)) {
  tmp = str_remove(prot_mut.df$species[i], "_prot")
  prot_mut.df$species[i] = tmp
}

for(i in 1:nrow(rest_prot_mut.df)) {
  tmp = str_remove(rest_prot_mut.df$species[i], "_prot")
  rest_prot_mut.df$species[i] = tmp
}

# Reducing to only the species we're interested in
prot_mutRed.df = prot_mut.df[which(prot_mut.df$species %in% prot_sp),]
rest_prot_mutRed.df = rest_prot_mut.df[which(rest_prot_mut.df$species %in% prot_sp),]

ec_low_AT_prot = 0
ec_low_AG_prot = 0
ec_low_AC_prot = 0
ec_low_TA_prot = 0
ec_low_TG_prot = 0
ec_low_TC_prot = 0
ec_low_GA_prot = 0
ec_low_GT_prot = 0
ec_low_GC_prot = 0
ec_low_CA_prot = 0
ec_low_CT_prot = 0
ec_low_CG_prot = 0

ec_med_AT_prot = 0
ec_med_AG_prot = 0
ec_med_AC_prot = 0
ec_med_TA_prot = 0
ec_med_TG_prot = 0
ec_med_TC_prot = 0
ec_med_GA_prot = 0
ec_med_GT_prot = 0
ec_med_GC_prot = 0
ec_med_CA_prot = 0
ec_med_CT_prot = 0
ec_med_CG_prot = 0

ec_high_AT_prot = 0
ec_high_AG_prot = 0
ec_high_AC_prot = 0
ec_high_TA_prot = 0
ec_high_TG_prot = 0
ec_high_TC_prot = 0
ec_high_GA_prot = 0
ec_high_GT_prot = 0
ec_high_GC_prot = 0
ec_high_CA_prot = 0
ec_high_CT_prot = 0
ec_high_CG_prot = 0

ec_low_AT_rest_prot = 0
ec_low_AG_rest_prot = 0
ec_low_AC_rest_prot = 0
ec_low_TA_rest_prot = 0
ec_low_TG_rest_prot = 0
ec_low_TC_rest_prot = 0
ec_low_GA_rest_prot = 0
ec_low_GT_rest_prot = 0
ec_low_GC_rest_prot = 0
ec_low_CA_rest_prot = 0
ec_low_CT_rest_prot = 0
ec_low_CG_rest_prot = 0

ec_med_AT_rest_prot = 0
ec_med_AG_rest_prot = 0
ec_med_AC_rest_prot = 0
ec_med_TA_rest_prot = 0
ec_med_TG_rest_prot = 0
ec_med_TC_rest_prot = 0
ec_med_GA_rest_prot = 0
ec_med_GT_rest_prot = 0
ec_med_GC_rest_prot = 0
ec_med_CA_rest_prot = 0
ec_med_CT_rest_prot = 0
ec_med_CG_rest_prot = 0

ec_high_AT_rest_prot = 0
ec_high_AG_rest_prot = 0
ec_high_AC_rest_prot = 0
ec_high_TA_rest_prot = 0
ec_high_TG_rest_prot = 0
ec_high_TC_rest_prot = 0
ec_high_GA_rest_prot = 0
ec_high_GT_rest_prot = 0
ec_high_GC_rest_prot = 0
ec_high_CA_rest_prot = 0
ec_high_CT_rest_prot = 0
ec_high_CG_rest_prot = 0

ec_strain = c("Escherichia_coli_low_mg1655", "Escherichia_coli_med_mg1655", "Escherichia_coli_high_mg1655",
              "Escherichia_coli_low_atcc", "Escherichia_coli_med_atcc", "Escherichia_coli_high_atcc")

# Making all transition/transversion columns numeric
for(i in 2:ncol(prot_mutRed.df)) {
  prot_mutRed.df[,i] = as.numeric(prot_mutRed.df[,i])
  rest_prot_mutRed.df[,i] = as.numeric(rest_prot_mutRed.df[,i])
}

# Combining protein-coding and rest data from both strains
for(strain in ec_strain) {
  expLvl = unlist(strsplit(strain, "_"))[3]
  
  if(expLvl == "low") {
    ec_low_AT_prot = ec_low_AT_prot + prot_mutRed.df$`A>T`[which(prot_mutRed.df$species == strain)]
    ec_low_AG_prot = ec_low_AG_prot + prot_mutRed.df$`A>G`[which(prot_mutRed.df$species == strain)]
    ec_low_AC_prot = ec_low_AC_prot + prot_mutRed.df$`A>C`[which(prot_mutRed.df$species == strain)]
    ec_low_TA_prot = ec_low_TA_prot + prot_mutRed.df$`T>A`[which(prot_mutRed.df$species == strain)]
    ec_low_TG_prot = ec_low_TG_prot + prot_mutRed.df$`T>G`[which(prot_mutRed.df$species == strain)]
    ec_low_TC_prot = ec_low_TC_prot + prot_mutRed.df$`T>C`[which(prot_mutRed.df$species == strain)]
    ec_low_GA_prot = ec_low_GA_prot + prot_mutRed.df$`G>A`[which(prot_mutRed.df$species == strain)]
    ec_low_GT_prot = ec_low_GT_prot + prot_mutRed.df$`G>T`[which(prot_mutRed.df$species == strain)]
    ec_low_GC_prot = ec_low_GC_prot + prot_mutRed.df$`G>C`[which(prot_mutRed.df$species == strain)]
    ec_low_CA_prot = ec_low_CA_prot + prot_mutRed.df$`C>A`[which(prot_mutRed.df$species == strain)]
    ec_low_CT_prot = ec_low_CT_prot + prot_mutRed.df$`C>T`[which(prot_mutRed.df$species == strain)]
    ec_low_CG_prot = ec_low_CG_prot + prot_mutRed.df$`C>G`[which(prot_mutRed.df$species == strain)]
    
    ec_low_AT_rest_prot = ec_low_AT_rest_prot + rest_prot_mutRed.df$`A>T`[which(rest_prot_mutRed.df$species == strain)]
    ec_low_AG_rest_prot = ec_low_AG_rest_prot + rest_prot_mutRed.df$`A>G`[which(rest_prot_mutRed.df$species == strain)]
    ec_low_AC_rest_prot = ec_low_AC_rest_prot + rest_prot_mutRed.df$`A>C`[which(rest_prot_mutRed.df$species == strain)]
    ec_low_TA_rest_prot = ec_low_TA_rest_prot + rest_prot_mutRed.df$`T>A`[which(rest_prot_mutRed.df$species == strain)]
    ec_low_TG_rest_prot = ec_low_TG_rest_prot + rest_prot_mutRed.df$`T>G`[which(rest_prot_mutRed.df$species == strain)]
    ec_low_TC_rest_prot = ec_low_TC_rest_prot + rest_prot_mutRed.df$`T>C`[which(rest_prot_mutRed.df$species == strain)]
    ec_low_GA_rest_prot = ec_low_GA_rest_prot + rest_prot_mutRed.df$`G>A`[which(rest_prot_mutRed.df$species == strain)]
    ec_low_GT_rest_prot = ec_low_GT_rest_prot + rest_prot_mutRed.df$`G>T`[which(rest_prot_mutRed.df$species == strain)]
    ec_low_GC_rest_prot = ec_low_GC_rest_prot + rest_prot_mutRed.df$`G>C`[which(rest_prot_mutRed.df$species == strain)]
    ec_low_CA_rest_prot = ec_low_CA_rest_prot + rest_prot_mutRed.df$`C>A`[which(rest_prot_mutRed.df$species == strain)]
    ec_low_CT_rest_prot = ec_low_CT_rest_prot + rest_prot_mutRed.df$`C>T`[which(rest_prot_mutRed.df$species == strain)]
    ec_low_CG_rest_prot = ec_low_CG_rest_prot + rest_prot_mutRed.df$`C>G`[which(rest_prot_mutRed.df$species == strain)]
  }
  
  if(expLvl == "med") {
    ec_med_AT_prot = ec_med_AT_prot + prot_mutRed.df$`A>T`[which(prot_mutRed.df$species == strain)]
    ec_med_AG_prot = ec_med_AG_prot + prot_mutRed.df$`A>G`[which(prot_mutRed.df$species == strain)]
    ec_med_AC_prot = ec_med_AC_prot + prot_mutRed.df$`A>C`[which(prot_mutRed.df$species == strain)]
    ec_med_TA_prot = ec_med_TA_prot + prot_mutRed.df$`T>A`[which(prot_mutRed.df$species == strain)]
    ec_med_TG_prot = ec_med_TG_prot + prot_mutRed.df$`T>G`[which(prot_mutRed.df$species == strain)]
    ec_med_TC_prot = ec_med_TC_prot + prot_mutRed.df$`T>C`[which(prot_mutRed.df$species == strain)]
    ec_med_GA_prot = ec_med_GA_prot + prot_mutRed.df$`G>A`[which(prot_mutRed.df$species == strain)]
    ec_med_GT_prot = ec_med_GT_prot + prot_mutRed.df$`G>T`[which(prot_mutRed.df$species == strain)]
    ec_med_GC_prot = ec_med_GC_prot + prot_mutRed.df$`G>C`[which(prot_mutRed.df$species == strain)]
    ec_med_CA_prot = ec_med_CA_prot + prot_mutRed.df$`C>A`[which(prot_mutRed.df$species == strain)]
    ec_med_CT_prot = ec_med_CT_prot + prot_mutRed.df$`C>T`[which(prot_mutRed.df$species == strain)]
    ec_med_CG_prot = ec_med_CG_prot + prot_mutRed.df$`C>G`[which(prot_mutRed.df$species == strain)]
    
    ec_med_AT_rest_prot = ec_med_AT_rest_prot + rest_prot_mutRed.df$`A>T`[which(rest_prot_mutRed.df$species == strain)]
    ec_med_AG_rest_prot = ec_med_AG_rest_prot + rest_prot_mutRed.df$`A>G`[which(rest_prot_mutRed.df$species == strain)]
    ec_med_AC_rest_prot = ec_med_AC_rest_prot + rest_prot_mutRed.df$`A>C`[which(rest_prot_mutRed.df$species == strain)]
    ec_med_TA_rest_prot = ec_med_TA_rest_prot + rest_prot_mutRed.df$`T>A`[which(rest_prot_mutRed.df$species == strain)]
    ec_med_TG_rest_prot = ec_med_TG_rest_prot + rest_prot_mutRed.df$`T>G`[which(rest_prot_mutRed.df$species == strain)]
    ec_med_TC_rest_prot = ec_med_TC_rest_prot + rest_prot_mutRed.df$`T>C`[which(rest_prot_mutRed.df$species == strain)]
    ec_med_GA_rest_prot = ec_med_GA_rest_prot + rest_prot_mutRed.df$`G>A`[which(rest_prot_mutRed.df$species == strain)]
    ec_med_GT_rest_prot = ec_med_GT_rest_prot + rest_prot_mutRed.df$`G>T`[which(rest_prot_mutRed.df$species == strain)]
    ec_med_GC_rest_prot = ec_med_GC_rest_prot + rest_prot_mutRed.df$`G>C`[which(rest_prot_mutRed.df$species == strain)]
    ec_med_CA_rest_prot = ec_med_CA_rest_prot + rest_prot_mutRed.df$`C>A`[which(rest_prot_mutRed.df$species == strain)]
    ec_med_CT_rest_prot = ec_med_CT_rest_prot + rest_prot_mutRed.df$`C>T`[which(rest_prot_mutRed.df$species == strain)]
    ec_med_CG_rest_prot = ec_med_CG_rest_prot + rest_prot_mutRed.df$`C>G`[which(rest_prot_mutRed.df$species == strain)]
  }
  
  if(expLvl == "high") {
    ec_high_AT_prot = ec_high_AT_prot + prot_mutRed.df$`A>T`[which(prot_mutRed.df$species == strain)]
    ec_high_AG_prot = ec_high_AG_prot + prot_mutRed.df$`A>G`[which(prot_mutRed.df$species == strain)]
    ec_high_AC_prot = ec_high_AC_prot + prot_mutRed.df$`A>C`[which(prot_mutRed.df$species == strain)]
    ec_high_TA_prot = ec_high_TA_prot + prot_mutRed.df$`T>A`[which(prot_mutRed.df$species == strain)]
    ec_high_TG_prot = ec_high_TG_prot + prot_mutRed.df$`T>G`[which(prot_mutRed.df$species == strain)]
    ec_high_TC_prot = ec_high_TC_prot + prot_mutRed.df$`T>C`[which(prot_mutRed.df$species == strain)]
    ec_high_GA_prot = ec_high_GA_prot + prot_mutRed.df$`G>A`[which(prot_mutRed.df$species == strain)]
    ec_high_GT_prot = ec_high_GT_prot + prot_mutRed.df$`G>T`[which(prot_mutRed.df$species == strain)]
    ec_high_GC_prot = ec_high_GC_prot + prot_mutRed.df$`G>C`[which(prot_mutRed.df$species == strain)]
    ec_high_CA_prot = ec_high_CA_prot + prot_mutRed.df$`C>A`[which(prot_mutRed.df$species == strain)]
    ec_high_CT_prot = ec_high_CT_prot + prot_mutRed.df$`C>T`[which(prot_mutRed.df$species == strain)]
    ec_high_CG_prot = ec_high_CG_prot + prot_mutRed.df$`C>G`[which(prot_mutRed.df$species == strain)]
    
    ec_high_AT_rest_prot = ec_high_AT_rest_prot + rest_prot_mutRed.df$`A>T`[which(rest_prot_mutRed.df$species == strain)]
    ec_high_AG_rest_prot = ec_high_AG_rest_prot + rest_prot_mutRed.df$`A>G`[which(rest_prot_mutRed.df$species == strain)]
    ec_high_AC_rest_prot = ec_high_AC_rest_prot + rest_prot_mutRed.df$`A>C`[which(rest_prot_mutRed.df$species == strain)]
    ec_high_TA_rest_prot = ec_high_TA_rest_prot + rest_prot_mutRed.df$`T>A`[which(rest_prot_mutRed.df$species == strain)]
    ec_high_TG_rest_prot = ec_high_TG_rest_prot + rest_prot_mutRed.df$`T>G`[which(rest_prot_mutRed.df$species == strain)]
    ec_high_TC_rest_prot = ec_high_TC_rest_prot + rest_prot_mutRed.df$`T>C`[which(rest_prot_mutRed.df$species == strain)]
    ec_high_GA_rest_prot = ec_high_GA_rest_prot + rest_prot_mutRed.df$`G>A`[which(rest_prot_mutRed.df$species == strain)]
    ec_high_GT_rest_prot = ec_high_GT_rest_prot + rest_prot_mutRed.df$`G>T`[which(rest_prot_mutRed.df$species == strain)]
    ec_high_GC_rest_prot = ec_high_GC_rest_prot + rest_prot_mutRed.df$`G>C`[which(rest_prot_mutRed.df$species == strain)]
    ec_high_CA_rest_prot = ec_high_CA_rest_prot + rest_prot_mutRed.df$`C>A`[which(rest_prot_mutRed.df$species == strain)]
    ec_high_CT_rest_prot = ec_high_CT_rest_prot + rest_prot_mutRed.df$`C>T`[which(rest_prot_mutRed.df$species == strain)]
    ec_high_CG_rest_prot = ec_high_CG_rest_prot + rest_prot_mutRed.df$`C>G`[which(rest_prot_mutRed.df$species == strain)]
  }
}

ec_low_comb_prot = c("Escherichia_coli_low", ec_low_AT_prot, ec_low_AG_prot, ec_low_AC_prot, ec_low_TA_prot, ec_low_TG_prot, ec_low_TC_prot, 
                     ec_low_GA_prot, ec_low_GT_prot, ec_low_GC_prot, ec_low_CA_prot, ec_low_CT_prot, ec_low_CG_prot)
ec_low_comb_rest_prot = c("Escherichia_coli_low", ec_low_AT_rest_prot, ec_low_AG_rest_prot, ec_low_AC_rest_prot, ec_low_TA_rest_prot, 
                          ec_low_TG_rest_prot, ec_low_TC_rest_prot, ec_low_GA_rest_prot, ec_low_GT_rest_prot, ec_low_GC_rest_prot, 
                          ec_low_CA_rest_prot, ec_low_CT_rest_prot, ec_low_CG_rest_prot)

ec_med_comb_prot = c("Escherichia_coli_med", ec_med_AT_prot, ec_med_AG_prot, ec_med_AC_prot, ec_med_TA_prot, ec_med_TG_prot, ec_med_TC_prot, 
                     ec_med_GA_prot, ec_med_GT_prot, ec_med_GC_prot, ec_med_CA_prot, ec_med_CT_prot, ec_med_CG_prot)
ec_med_comb_rest_prot = c("Escherichia_coli_med", ec_med_AT_rest_prot, ec_med_AG_rest_prot, ec_med_AC_rest_prot, ec_med_TA_rest_prot, 
                          ec_med_TG_rest_prot, ec_med_TC_rest_prot, ec_med_GA_rest_prot, ec_med_GT_rest_prot, ec_med_GC_rest_prot, 
                          ec_med_CA_rest_prot, ec_med_CT_rest_prot, ec_med_CG_rest_prot)

ec_high_comb_prot = c("Escherichia_coli_high", ec_high_AT_prot, ec_high_AG_prot, ec_high_AC_prot, ec_high_TA_prot, ec_high_TG_prot, ec_high_TC_prot, 
                     ec_high_GA_prot, ec_high_GT_prot, ec_high_GC_prot, ec_high_CA_prot, ec_high_CT_prot, ec_high_CG_prot)
ec_high_comb_rest_prot = c("Escherichia_coli_high", ec_high_AT_rest_prot, ec_high_AG_rest_prot, ec_high_AC_rest_prot, ec_high_TA_rest_prot, 
                          ec_high_TG_rest_prot, ec_high_TC_rest_prot, ec_high_GA_rest_prot, ec_high_GT_rest_prot, ec_high_GC_rest_prot, 
                          ec_high_CA_rest_prot, ec_high_CT_rest_prot, ec_high_CG_rest_prot)

# Adding combined data to the tables
prot_mutRed.df = rbind(prot_mutRed.df, ec_low_comb_prot)
prot_mutRed.df = rbind(prot_mutRed.df, ec_high_comb_prot)
prot_mutRed.df = rbind(prot_mutRed.df, ec_med_comb_prot)

rest_prot_mutRed.df = rbind(rest_prot_mutRed.df, ec_low_comb_rest_prot)
rest_prot_mutRed.df = rbind(rest_prot_mutRed.df, ec_high_comb_rest_prot)
rest_prot_mutRed.df = rbind(rest_prot_mutRed.df, ec_med_comb_rest_prot)

#########################################################################################
# Reading in the table of nucleotide numbers in tRNA genes and the rest of the genome
nt_prot.df = read.table("ntFreq_prot.tab", h = T)
nt_rest_prot.df = read.table("ntFreq_rest_prot.tab", h = T)

# Removing the _prot string from species names in nt_prot.df & nt_rest_prot.df
for(i in 1:nrow(nt_prot.df)) {
  tmp = str_remove(nt_prot.df$species[i], "_prot")
  nt_prot.df$species[i] = tmp
}

for(i in 1:nrow(nt_rest_prot.df)) {
  tmp = str_remove(nt_rest_prot.df$species[i], "_prot")
  nt_rest_prot.df$species[i] = tmp
}

# Making sure the tables only have the species we're interested in
nt_prot.df = nt_prot.df[which(nt_prot.df$species %in% prot_sp),]
nt_rest_prot.df = nt_rest_prot.df[which(nt_rest_prot.df$species %in% prot_sp),]

# Reducing nt freqs to only non-protein-coding areas (removing high/med freqs from low, etc)
sc_genome = readDNAStringSet("Scerevisiae/DB/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz")
ec_atcc_genome = readDNAStringSet("Ecoli/DB/CP000946.1.fasta")
ec_mg_genome = readDNAStringSet("Ecoli/DB/U00096.2.fasta")

############################################################################################################
# Fixing the S. cerevisiae nt frequencies
vSC = c(freqA = 0, freqT = 0, freqG = 0, freqC = 0)
sc_list = c("Saccharomyces_cerevisiae_low", "Saccharomyces_cerevisiae_high", "Saccharomyces_cerevisiae_med")

for(strain in sc_list) {
  vSC["freqA"] = vSC["freqA"] + nt_prot.df$freq_A[which(nt_prot.df$species == strain)]
  vSC["freqT"] = vSC["freqT"] + nt_prot.df$freq_T[which(nt_prot.df$species == strain)]
  vSC["freqG"] = vSC["freqG"] + nt_prot.df$freq_G[which(nt_prot.df$species == strain)]
  vSC["freqC"] = vSC["freqC"] + nt_prot.df$freq_C[which(nt_prot.df$species == strain)]
}

restSC = c(restA = 0, restT = 0, restG = 0, restC = 0)
restSC["restA"] = sum(letterFrequency(sc_genome, "A")) - vSC["freqA"]
restSC["restT"] = sum(letterFrequency(sc_genome, "T")) - vSC["freqT"]
restSC["restG"] = sum(letterFrequency(sc_genome, "G")) - vSC["freqG"]
restSC["restC"] = sum(letterFrequency(sc_genome, "C")) - vSC["freqC"]

nt_rest_prot.df$freq_A[which(nt_rest_prot.df$species %in% sc_list)] = restSC["restA"]
nt_rest_prot.df$freq_T[which(nt_rest_prot.df$species %in% sc_list)] = restSC["restT"]
nt_rest_prot.df$freq_G[which(nt_rest_prot.df$species %in% sc_list)] = restSC["restG"]
nt_rest_prot.df$freq_C[which(nt_rest_prot.df$species %in% sc_list)] = restSC["restC"]

############################################################################################################
# Fixing the E. coli nt frequencies
vECMG = c(freqA = 0, freqT = 0, freqG = 0, freqC = 0)
vECAT = c(freqA = 0, freqT = 0, freqG = 0, freqC = 0)

# MG1655 strain
ecmg_list = c("Escherichia_coli_low_mg1655", "Escherichia_coli_high_mg1655", "Escherichia_coli_med_mg1655")

for(strain in ecmg_list) {
  vECMG["freqA"] = vECMG["freqA"] + nt_prot.df$freq_A[which(nt_prot.df$species == strain)]
  vECMG["freqT"] = vECMG["freqT"] + nt_prot.df$freq_T[which(nt_prot.df$species == strain)]
  vECMG["freqG"] = vECMG["freqG"] + nt_prot.df$freq_G[which(nt_prot.df$species == strain)]
  vECMG["freqC"] = vECMG["freqC"] + nt_prot.df$freq_C[which(nt_prot.df$species == strain)]
}

restECMG = c(restA = 0, restT = 0, restG = 0, restC = 0)
restECMG["restA"] = sum(letterFrequency(ec_mg_genome, "A")) - vECMG["freqA"]
restECMG["restT"] = sum(letterFrequency(ec_mg_genome, "T")) - vECMG["freqT"]
restECMG["restG"] = sum(letterFrequency(ec_mg_genome, "G")) - vECMG["freqG"]
restECMG["restC"] = sum(letterFrequency(ec_mg_genome, "C")) - vECMG["freqC"]

# ATCC 8739 strain
ecat_list = c("Escherichia_coli_low_atcc", "Escherichia_coli_high_atcc", "Escherichia_coli_med_atcc")

for(strain in ecat_list) {
  vECAT["freqA"] = vECAT["freqA"] + nt_prot.df$freq_A[which(nt_prot.df$species == strain)]
  vECAT["freqT"] = vECAT["freqT"] + nt_prot.df$freq_T[which(nt_prot.df$species == strain)]
  vECAT["freqG"] = vECAT["freqG"] + nt_prot.df$freq_G[which(nt_prot.df$species == strain)]
  vECAT["freqC"] = vECAT["freqC"] + nt_prot.df$freq_C[which(nt_prot.df$species == strain)]
}

restECAT = c(restA = 0, restT = 0, restG = 0, restC = 0)
restECAT["restA"] = sum(letterFrequency(ec_atcc_genome, "A")) - vECAT["freqA"]
restECAT["restT"] = sum(letterFrequency(ec_atcc_genome, "T")) - vECAT["freqT"]
restECAT["restG"] = sum(letterFrequency(ec_atcc_genome, "G")) - vECAT["freqG"]
restECAT["restC"] = sum(letterFrequency(ec_atcc_genome, "C")) - vECAT["freqC"]

# Combining the strains
restEC = c(restA = 0, restT = 0, restG = 0, restC = 0)
restEC["restA"] = restECMG["restA"] + restECAT["restA"]
restEC["restT"] = restECMG["restT"] + restECAT["restT"]
restEC["restG"] = restECMG["restG"] + restECAT["restG"]
restEC["restC"] = restECMG["restC"] + restECAT["restC"]

ec_low_nt_rest_comb = c("Escherichia_coli_low", restEC["restA"], restEC["restT"], restEC["restG"], restEC["restC"])
ec_high_nt_rest_comb = c("Escherichia_coli_high", restEC["restA"], restEC["restT"], restEC["restG"], restEC["restC"])
ec_med_nt_rest_comb = c("Escherichia_coli_med", restEC["restA"], restEC["restT"], restEC["restG"], restEC["restC"])

nt_rest_prot.df = rbind(nt_rest_prot.df, ec_low_nt_rest_comb)
nt_rest_prot.df = rbind(nt_rest_prot.df, ec_high_nt_rest_comb)
nt_rest_prot.df = rbind(nt_rest_prot.df, ec_med_nt_rest_comb)

############################################################################################################

freqA_low = 0
freqT_low = 0
freqG_low = 0
freqC_low = 0

freqA_med = 0
freqT_med = 0
freqG_med = 0
freqC_med = 0

freqA_high = 0
freqT_high = 0
freqG_high = 0
freqC_high = 0

ec_list = append(ecat_list, ecmg_list)

for(strain in ec_list) {
  expLvl = unlist(strsplit(strain, "_"))[3]
  
  if(expLvl == "low") {
    freqA_low = freqA_low + nt_prot.df$freq_A[which(nt_prot.df$species == strain)]
    freqT_low = freqT_low + nt_prot.df$freq_T[which(nt_prot.df$species == strain)]
    freqG_low = freqG_low + nt_prot.df$freq_G[which(nt_prot.df$species == strain)]
    freqC_low = freqC_low + nt_prot.df$freq_C[which(nt_prot.df$species == strain)]
  }
  
  if(expLvl == "med") {
    freqA_med = freqA_med + nt_prot.df$freq_A[which(nt_prot.df$species == strain)]
    freqT_med = freqT_med + nt_prot.df$freq_T[which(nt_prot.df$species == strain)]
    freqG_med = freqG_med + nt_prot.df$freq_G[which(nt_prot.df$species == strain)]
    freqC_med = freqC_med + nt_prot.df$freq_C[which(nt_prot.df$species == strain)]
  }
  
  if(expLvl == "high") {
    freqA_high = freqA_high + nt_prot.df$freq_A[which(nt_prot.df$species == strain)]
    freqT_high = freqT_high + nt_prot.df$freq_T[which(nt_prot.df$species == strain)]
    freqG_high = freqG_high + nt_prot.df$freq_G[which(nt_prot.df$species == strain)]
    freqC_high = freqC_high + nt_prot.df$freq_C[which(nt_prot.df$species == strain)]
  }
}

ec_low = c("Escherichia_coli_low", freqA_low, freqT_low, freqG_low, freqC_low)
ec_med = c("Escherichia_coli_med", freqA_med, freqT_med, freqG_med, freqC_med)
ec_high = c("Escherichia_coli_high", freqA_high, freqT_high, freqG_high, freqC_high)

nt_prot.df = rbind(nt_prot.df, ec_low, ec_high, ec_med)

# Reducing to just combined data
spNew = c("Saccharomyces_cerevisiae_low", "Saccharomyces_cerevisiae_med", "Saccharomyces_cerevisiae_high", "Escherichia_coli_low", 
          "Escherichia_coli_med", "Escherichia_coli_high")

prot_mutRed.df = prot_mutRed.df[which(prot_mutRed.df$species %in% spNew),]
rest_prot_mutRed.df = rest_prot_mutRed.df[which(rest_prot_mutRed.df$species %in% spNew),]
nt_prot.df = nt_prot.df[which(nt_prot.df$species %in% spNew),]
nt_rest_prot.df = nt_rest_prot.df[which(nt_rest_prot.df$species %in% spNew),]

##########################################################################################
# Combining all expression levels to get just protein-coding vs. rest - no bins
ec_list = c("Escherichia_coli_low", "Escherichia_coli_high", "Escherichia_coli_med")

# Making sure all numeric columns are actually numeric
prot_mutRed.df[,2:13] = sapply(prot_mutRed.df[,2:13], as.numeric)
rest_prot_mutRed.df[,2:13] = sapply(rest_prot_mutRed.df[,2:13], as.numeric)
nt_prot.df[,2:5] = sapply(nt_prot.df[,2:5], as.numeric)

# Combining for nt_prot.df
scA = sum(nt_prot.df$freq_A[which(nt_prot.df$species %in% sc_list)])
scT = sum(nt_prot.df$freq_T[which(nt_prot.df$species %in% sc_list)])
scG = sum(nt_prot.df$freq_G[which(nt_prot.df$species %in% sc_list)])
scC = sum(nt_prot.df$freq_C[which(nt_prot.df$species %in% sc_list)])

ecA = sum(nt_prot.df$freq_A[which(nt_prot.df$species %in% ec_list)])
ecT = sum(nt_prot.df$freq_T[which(nt_prot.df$species %in% ec_list)])
ecG = sum(nt_prot.df$freq_G[which(nt_prot.df$species %in% ec_list)])
ecC = sum(nt_prot.df$freq_C[which(nt_prot.df$species %in% ec_list)])

scTot = c("Saccharomyces_cerevisiae", scA, scT, scG, scC)
ecTot = c("Escherichia_coli", ecA, ecT, ecG, ecC)

nt_prot_tot.df = rbind(scTot, ecTot)
colnames(nt_prot_tot.df) = colnames(nt_prot.df)

nt_prot_tot.df = as.data.frame(nt_prot_tot.df)
nt_prot_tot.df[,2:5] = sapply(nt_prot_tot.df[,2:5], as.numeric)

# Getting an equivalent version of nt_rest_prot.df
nt_rest_prot_tot.df = nt_rest_prot.df[c(1,4),]
nt_rest_prot_tot.df$species[1] = "Saccharomyces_cerevisiae"
nt_rest_prot_tot.df$species[2] = "Escherichia_coli"

# Combining for prot_mutRed.df
prot_mutRed.df[7,1] = "Saccharomyces_cerevisiae"
prot_mutRed.df[8,1] = "Escherichia_coli"

for(i in 2:13) {
  prot_mutRed.df[7,i] = sum(prot_mutRed.df[which(prot_mutRed.df$species %in% sc_list),i])
  prot_mutRed.df[8,i] = sum(prot_mutRed.df[which(prot_mutRed.df$species %in% ec_list),i])
}

prot_mutRed_tot.df = prot_mutRed.df[7:8,]

# Combining for rest_prot_mutRed.df
rest_prot_mutRed.df[7,1] = "Saccharomyces_cerevisiae"
rest_prot_mutRed.df[8,1] = "Escherichia_coli"

for(i in 2:13) {
  rest_prot_mutRed.df[7,i] = sum(rest_prot_mutRed.df[which(rest_prot_mutRed.df$species %in% sc_list),i])
  rest_prot_mutRed.df[8,i] = sum(rest_prot_mutRed.df[which(rest_prot_mutRed.df$species %in% ec_list),i])
  
  # Need to subtract out mutations in other protein-coding genes
  # Each would've been counted twice - once for each other expression level
  rest_prot_mutRed.df[7,i] = rest_prot_mutRed.df[7,i] - (2*prot_mutRed.df[7,i])
}

rest_prot_mutRed_tot.df = rest_prot_mutRed.df[7:8,]

# Removing new rows from original tables in case they are needed again
prot_mutRed.df = prot_mutRed.df[-c(7:8),]
rest_prot_mutRed.df = rest_prot_mutRed.df[-c(7:8),]

##########################################################################################
# Getting mutations/total # nucleotide

# Making copies of the tables for ntMutPlot graphs
prot_nt.df = prot_mutRed_tot.df
rest_prot_nt.df = rest_prot_mutRed_tot.df

# Columns by ref base
Amut = c("A>T", "A>G", "A>C")
Tmut = c("T>A", "T>G", "T>C")
Gmut = c("G>A", "G>T", "G>C")
Cmut = c("C>A", "C>T", "C>G")

# Dividing number of mutations by number of reference base occurrences in genes of interest
# Replacing number of mutations in table with new number calculated
spNew = c("Saccharomyces_cerevisiae", "Escherichia_coli")
for(sp in spNew) {
  fA_prot = as.numeric(nt_prot_tot.df$freq_A[which(nt_prot_tot.df$species == sp)])
  fT_prot = as.numeric(nt_prot_tot.df$freq_T[which(nt_prot_tot.df$species == sp)])
  fG_prot = as.numeric(nt_prot_tot.df$freq_G[which(nt_prot_tot.df$species == sp)])
  fC_prot = as.numeric(nt_prot_tot.df$freq_C[which(nt_prot_tot.df$species == sp)])
  
  fA_rest = as.numeric(nt_rest_prot_tot.df$freq_A[which(nt_rest_prot_tot.df$species == sp)])
  fT_rest = as.numeric(nt_rest_prot_tot.df$freq_T[which(nt_rest_prot_tot.df$species == sp)])
  fG_rest = as.numeric(nt_rest_prot_tot.df$freq_G[which(nt_rest_prot_tot.df$species == sp)])
  fC_rest = as.numeric(nt_rest_prot_tot.df$freq_C[which(nt_rest_prot_tot.df$species == sp)])
  
  for(a in Amut) {
    tmp = as.numeric(prot_nt.df[which(prot_nt.df$species == sp), a])
    newA = tmp/fA_prot
    prot_nt.df[which(prot_nt.df$species == sp), a] = newA
    
    tmp = as.numeric(rest_prot_nt.df[which(rest_prot_nt.df$species == sp), a])
    newA = tmp/fA_rest
    rest_prot_nt.df[which(rest_prot_nt.df$species == sp), a] = newA
  }
  
  for(t in Tmut) {
    tmp = as.numeric(prot_nt.df[which(prot_nt.df$species == sp), t])
    newT = tmp/fT_prot
    prot_nt.df[which(prot_nt.df$species == sp), t] = newT
    
    tmp = as.numeric(rest_prot_nt.df[which(rest_prot_nt.df$species == sp), t])
    newT = tmp/fT_rest
    rest_prot_nt.df[which(rest_prot_nt.df$species == sp), t] = newT
  }
  
  for(g in Gmut) {
    tmp = as.numeric(prot_nt.df[which(prot_nt.df$species == sp), g])
    newG = tmp/fG_prot
    prot_nt.df[which(prot_nt.df$species == sp), g] = newG
    
    tmp = as.numeric(rest_prot_nt.df[which(rest_prot_nt.df$species == sp), g])
    newG = tmp/fG_rest
    rest_prot_nt.df[which(rest_prot_nt.df$species == sp), g] = newG
  }
  
  for(c in Cmut) {
    tmp = as.numeric(prot_nt.df[which(prot_nt.df$species == sp), c])
    newC = tmp/fC_prot
    prot_nt.df[which(prot_nt.df$species == sp), c] = newC
    
    tmp = as.numeric(rest_prot_nt.df[which(rest_prot_nt.df$species == sp), c])
    newC = tmp/fC_rest
    rest_prot_nt.df[which(rest_prot_nt.df$species == sp), c] = newC
  }
}


# Making a table to hold the combined data
prot_mutRed_tot.df$type = "Protein Coding"
rest_prot_mutRed_tot.df$type = "rest"

prot_nt.df$type = "Protein Coding"
rest_prot_nt.df$type = "rest"

mut_row = nrow(prot_mutRed_tot.df) * 2
mut_full.df = data.frame(species = character(mut_row), v1 = numeric(mut_row), v2 = numeric(mut_row), v3 = numeric(mut_row), 
                         v4 = numeric(mut_row), v5 = numeric(mut_row), v6 = numeric(mut_row), v7 = numeric(mut_row), v8 = numeric(mut_row),
                         v9 = numeric(mut_row), v10 = numeric(mut_row), v11 = numeric(mut_row), v12 = numeric(mut_row),
                         type = character(mut_row))

mut_full_nt.df = data.frame(species = character(mut_row), v1 = numeric(mut_row), v2 = numeric(mut_row), v3 = numeric(mut_row), 
                            v4 = numeric(mut_row), v5 = numeric(mut_row), v6 = numeric(mut_row), v7 = numeric(mut_row), v8 = numeric(mut_row),
                            v9 = numeric(mut_row), v10 = numeric(mut_row), v11 = numeric(mut_row), v12 = numeric(mut_row),
                            type = character(mut_row))

colnames(mut_full.df) = c("species", "A>T", "A>G", "A>C", "T>A", "T>G", "T>C", "G>A", "G>T", "G>C", "C>A", "C>T", "C>G", "type")
colnames(mut_full_nt.df) = c("species", "A>T", "A>G", "A>C", "T>A", "T>G", "T>C", "G>A", "G>T", "G>C", "C>A", "C>T", "C>G", "type")

# Combining prot table with rest table
rowCnt = 1
for(i in 1:mut_row) {
  if(mut_row/i >= 2) {
    mut_full.df[i,] = prot_mutRed_tot.df[i,]
    mut_full_nt.df[i,] = prot_nt.df[i,]
  }
  
  else{
    mut_full.df[i,] = rest_prot_mutRed_tot.df[rowCnt,]
    mut_full_nt.df[i,] = rest_prot_nt.df[rowCnt,]
    rowCnt = rowCnt + 1
  }
}

df <- gather(mut_full.df, mut, mutNum, 2:13)
df$mutNum = as.numeric(df$mutNum)
df$fullType = paste(df$species, df$type, sep = "_")

df2 <- gather(mut_full_nt.df, mut, mutNum, 2:13)
df2$mutNum = as.numeric(df2$mutNum)
df2$fullType = paste(df2$species, df2$type, sep = "_")

typeList = unique(df$fullType)
typeNum = length(typeList)

# Adding the "total" column
df$total = 0
df2$total = 0
for(i in 1:typeNum) {
  tmpType = typeList[i]
  tmpSum = sum(df$mutNum[which(df$fullType == tmpType)])
  df$total[which(df$fullType == tmpType)] = tmpSum
  
  tmpSum2 = sum(df2$mutNum[which(df2$fullType == tmpType)])
  df2$total[which(df2$fullType == tmpType)] = tmpSum2
}

# Changing type to something more appropriate for the graph
df$type[which(df$type == "rest")] = "tRNA Genes & Intergenic Regions"
df2$type[which(df2$type == "rest")] = "tRNA Genes & Intergenic Regions"

# Getting the columns to display in the right order
typeLvls = c("Protein Coding", "tRNA Genes & Intergenic Regions")
df$type = factor(df$type, levels = typeLvls)
df2$type = factor(df2$type, levels = typeLvls)

# Getting the BPSs to display in the right order
mutlvls = c("A>T", "A>G", "A>C", "T>A", "T>G", "T>C", "G>A", "G>T", "G>C", "C>A", "C>T", "C>G")
df$BPS = factor(df$mut, levels = mutlvls)
df2$BPS = factor(df2$mut, levels = mutlvls)

# Getting the standard deviation
df$mutProp = df$mutNum / df$total
df$mutPropSD = 0

df2$mutProp = df2$mutNum / df2$total
df2$mutPropSD = 0

for(i in 1:nrow(df)) {
  tmp = prop.test(df$mutNum[i], df$total[i])
  sd = df$mutProp[i] - tmp$conf.int[1]
  df$mutPropSD[i] = sd
  
  tmp2 = prop.test(df2$mutNum[i], df2$total[i])
  sd2 = df2$mutProp[i] - tmp2$conf.int[1]
  df2$mutPropSD[i] = sd2
}

################################################################################################
# Plots

# Number of specific mutation divided by total number of mutations
protMutPlot <- ggplot(data = df, aes(x = mut, y = mutNum/total, fill = species, pattern = type)) +
  geom_bar_pattern(stat = "identity", position = position_dodge(preserve = "single"),
                   color = "black", 
                   pattern_fill = "black",
                   pattern_angle = 45,
                   pattern_density = 0.1,
                   pattern_spacing = 0.015,
                   pattern_key_scale_factor = 0.6) + 
  geom_errorbar(aes(ymin = mutProp-mutPropSD, ymax = mutProp+mutPropSD), position = "dodge") +
  #stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.8, position = position_dodge(0.9)) +
  scale_fill_manual(values = colorRampPalette(bpsPalette)(2), name = "Species", labels = c("E. coli", "S. cerevisiae")) +
  scale_pattern_manual(values = c(`Protein Coding` = "stripe", `tRNA Genes & Intergenic Regions` = "none")) +
  labs(title = "Mutation Frequencies by Species", x = "Mutation", y = "Frequency", pattern = "Type") + 
  guides(pattern = guide_legend(override.aes = list(fill = "white")),
         fill = guide_legend(override.aes = list(pattern = "none"))) +
  theme(legend.title = element_text(size=20), legend.text = element_text(size=20)) +
  theme_classic()

protMutPlot    # <-- Saved as mutFreqTot_prot_2.png locally


# Plot of number of mutated nucleotide divided by total number of same nucleotide (# A muts / # As)
# Error bars are excluded because they are not informative
protNTMutPlot <- ggplot(data = df2, aes(x = mut, y = mutNum, fill = species, pattern = type)) +
  geom_bar_pattern(stat = "identity", position = position_dodge(preserve = "single"),
                   color = "black", 
                   pattern_fill = "black",
                   pattern_angle = 45,
                   pattern_density = 0.1,
                   pattern_spacing = 0.015,
                   pattern_key_scale_factor = 0.6) + 
  #geom_errorbar(aes(ymin = mutProp-mutPropSD, ymax = mutProp+mutPropSD), position = "dodge") +
  #stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.8, position = position_dodge(0.9)) +
  scale_fill_manual(values = colorRampPalette(bpsPalette)(2), name = "Species", labels = c("E. coli", "S. cerevisiae")) +
  scale_pattern_manual(values = c(`Protein Coding` = "stripe", `tRNA Genes & Intergenic Regions` = "none")) +
  labs(title = "Mutation Frequencies by Species", x = "Mutation", y = "Frequency", pattern = "Type") + 
  guides(pattern = guide_legend(override.aes = list(fill = "white")),
         fill = guide_legend(override.aes = list(pattern = "none"))) +
  theme(legend.title = element_text(size=20), legend.text = element_text(size=20)) +
  theme_classic()

protNTMutPlot


# ntMutPlot but only E. coli data
df_ec = df2[which(df2$species == "Escherichia_coli"),]
protECMutPlot <- ggplot(data = df_ec, aes(x = mut, y = mutNum, fill = species, pattern = type)) +
  geom_bar_pattern(stat = "identity", position = position_dodge(preserve = "single"),
                   color = "black", 
                   pattern_fill = "black",
                   pattern_angle = 45,
                   pattern_density = 0.1,
                   pattern_spacing = 0.015,
                   pattern_key_scale_factor = 0.6) + 
  #geom_errorbar(aes(ymin = mutProp-mutPropSD, ymax = mutProp+mutPropSD), position = "dodge") +
  #stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.8, position = position_dodge(0.9)) +
  scale_fill_manual(values = colorRampPalette(bpsPalette)(2), name = "Species", labels = c("E. coli", "S. cerevisiae")) +
  scale_pattern_manual(values = c(`Protein Coding` = "stripe", `tRNA Genes & Intergenic Regions` = "none")) +
  labs(title = "Mutation Frequencies by Species", x = "Mutation", y = "Frequency", pattern = "Type") + 
  guides(pattern = guide_legend(override.aes = list(fill = "white")),
         fill = guide_legend(override.aes = list(pattern = "none"))) +
  theme(legend.title = element_text(size=20), legend.text = element_text(size=20)) +
  theme_classic()

protECMutPlot


# ntMutPlot but only S. cerevisiae data
df_sc = df2[which(df2$species == "Saccharomyces_cerevisiae"),]
tRNAscMutPlot <- ggplot(data = df_sc, aes(x = mut, y = mutNum, fill = species, pattern = type)) +
  geom_bar_pattern(stat = "identity", position = position_dodge(preserve = "single"),
                   color = "black", 
                   pattern_fill = "black",
                   pattern_angle = 45,
                   pattern_density = 0.1,
                   pattern_spacing = 0.015,
                   pattern_key_scale_factor = 0.6) + 
  #geom_errorbar(aes(ymin = mutProp-mutPropSD, ymax = mutProp+mutPropSD), position = "dodge") +
  #stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.8, position = position_dodge(0.9)) +
  scale_fill_manual(values = c("#E5323B"), name = "Species", labels = c("S. cerevisiae")) +
  scale_pattern_manual(values = c(`Protein Coding` = "stripe", `tRNA Genes & Intergenic Regions` = "none")) +
  labs(title = "Mutation Frequencies by Species", x = "Mutation", y = "Frequency", pattern = "Type") + 
  guides(pattern = guide_legend(override.aes = list(fill = "white")),
         fill = guide_legend(override.aes = list(pattern = "none"))) +
  theme(legend.title = element_text(size=20), legend.text = element_text(size=20)) +
  theme_classic()

tRNAscMutPlot


# ntMutPlot but normalized so that frequencies = 1
tRNAntMutPlot2 <- ggplot(data = df, aes(x = mut, y = mutNum/total, fill = species, pattern = type)) +
  geom_bar_pattern(stat = "identity", position = position_dodge(preserve = "single"),
                   color = "black", 
                   pattern_fill = "black",
                   pattern_angle = 45,
                   pattern_density = 0.1,
                   pattern_spacing = 0.015,
                   pattern_key_scale_factor = 0.6) + 
  #geom_errorbar(aes(ymin = mutProp-mutPropSD, ymax = mutProp+mutPropSD), position = "dodge") +
  #stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.8, position = position_dodge(0.9)) +
  scale_fill_manual(values = colorRampPalette(bpsPalette)(2), name = "Species", labels = c("E. coli", "S. cerevisiae")) +
  scale_pattern_manual(values = c(`Protein Coding` = "stripe", `tRNA Genes & Intergenic Regions` = "none")) +
  labs(title = "Mutation Frequencies by Species", x = "Mutation", y = "Frequency", pattern = "Type") + 
  guides(pattern = guide_legend(override.aes = list(fill = "white")),
         fill = guide_legend(override.aes = list(pattern = "none"))) +
  theme(legend.title = element_text(size=20), legend.text = element_text(size=20)) +
  theme_classic()

tRNAntMutPlot2







############################################################################################################################################
# Mutation spectra plot - tRNA & protein-coding only
tRNA_sp = c("Saccharomyces_cerevisiae", "Escherichia_coli_mg1655", "Escherichia_coli_atcc")

prot_sp = c("Saccharomyces_cerevisiae_low", "Saccharomyces_cerevisiae_med", "Saccharomyces_cerevisiae_high",
            "Escherichia_coli_low_mg1655", "Escherichia_coli_med_mg1655", "Escherichia_coli_high_mg1655",
            "Escherichia_coli_low_atcc", "Escherichia_coli_med_atcc", "Escherichia_coli_high_atcc")

# Removing the _prot string from species names in prot_mut.df & rest_prot_mut.df
for(i in 1:nrow(prot_mut.df)) {
  tmp = str_remove(prot_mut.df$species[i], "_prot")
  prot_mut.df$species[i] = tmp
}

# Reducing to only the species we're interested in
tRNA_mutRed.df = tRNA_mut.df[which(tRNA_mut.df$species %in% tRNA_sp),]
prot_mutRed.df = prot_mut.df[which(prot_mut.df$species %in% prot_sp),]

##########################################################################################################
# Combining data from both strains

ec_low_AT_prot = 0
ec_low_AG_prot = 0
ec_low_AC_prot = 0
ec_low_TA_prot = 0
ec_low_TG_prot = 0
ec_low_TC_prot = 0
ec_low_GA_prot = 0
ec_low_GT_prot = 0
ec_low_GC_prot = 0
ec_low_CA_prot = 0
ec_low_CT_prot = 0
ec_low_CG_prot = 0

ec_med_AT_prot = 0
ec_med_AG_prot = 0
ec_med_AC_prot = 0
ec_med_TA_prot = 0
ec_med_TG_prot = 0
ec_med_TC_prot = 0
ec_med_GA_prot = 0
ec_med_GT_prot = 0
ec_med_GC_prot = 0
ec_med_CA_prot = 0
ec_med_CT_prot = 0
ec_med_CG_prot = 0

ec_high_AT_prot = 0
ec_high_AG_prot = 0
ec_high_AC_prot = 0
ec_high_TA_prot = 0
ec_high_TG_prot = 0
ec_high_TC_prot = 0
ec_high_GA_prot = 0
ec_high_GT_prot = 0
ec_high_GC_prot = 0
ec_high_CA_prot = 0
ec_high_CT_prot = 0
ec_high_CG_prot = 0

ec_strain = c("Escherichia_coli_low_mg1655", "Escherichia_coli_med_mg1655", "Escherichia_coli_high_mg1655",
              "Escherichia_coli_low_atcc", "Escherichia_coli_med_atcc", "Escherichia_coli_high_atcc")

# Making all transition/transversion columns numeric
for(i in 2:ncol(prot_mutRed.df)) {
  prot_mutRed.df[,i] = as.numeric(prot_mutRed.df[,i])
}

for(strain in ec_strain) {
  expLvl = unlist(strsplit(strain, "_"))[3]
  
  if(expLvl == "low") {
    ec_low_AT_prot = ec_low_AT_prot + prot_mutRed.df$`A>T`[which(prot_mutRed.df$species == strain)]
    ec_low_AG_prot = ec_low_AG_prot + prot_mutRed.df$`A>G`[which(prot_mutRed.df$species == strain)]
    ec_low_AC_prot = ec_low_AC_prot + prot_mutRed.df$`A>C`[which(prot_mutRed.df$species == strain)]
    ec_low_TA_prot = ec_low_TA_prot + prot_mutRed.df$`T>A`[which(prot_mutRed.df$species == strain)]
    ec_low_TG_prot = ec_low_TG_prot + prot_mutRed.df$`T>G`[which(prot_mutRed.df$species == strain)]
    ec_low_TC_prot = ec_low_TC_prot + prot_mutRed.df$`T>C`[which(prot_mutRed.df$species == strain)]
    ec_low_GA_prot = ec_low_GA_prot + prot_mutRed.df$`G>A`[which(prot_mutRed.df$species == strain)]
    ec_low_GT_prot = ec_low_GT_prot + prot_mutRed.df$`G>T`[which(prot_mutRed.df$species == strain)]
    ec_low_GC_prot = ec_low_GC_prot + prot_mutRed.df$`G>C`[which(prot_mutRed.df$species == strain)]
    ec_low_CA_prot = ec_low_CA_prot + prot_mutRed.df$`C>A`[which(prot_mutRed.df$species == strain)]
    ec_low_CT_prot = ec_low_CT_prot + prot_mutRed.df$`C>T`[which(prot_mutRed.df$species == strain)]
    ec_low_CG_prot = ec_low_CG_prot + prot_mutRed.df$`C>G`[which(prot_mutRed.df$species == strain)]
  }
  
  if(expLvl == "med") {
    ec_med_AT_prot = ec_med_AT_prot + prot_mutRed.df$`A>T`[which(prot_mutRed.df$species == strain)]
    ec_med_AG_prot = ec_med_AG_prot + prot_mutRed.df$`A>G`[which(prot_mutRed.df$species == strain)]
    ec_med_AC_prot = ec_med_AC_prot + prot_mutRed.df$`A>C`[which(prot_mutRed.df$species == strain)]
    ec_med_TA_prot = ec_med_TA_prot + prot_mutRed.df$`T>A`[which(prot_mutRed.df$species == strain)]
    ec_med_TG_prot = ec_med_TG_prot + prot_mutRed.df$`T>G`[which(prot_mutRed.df$species == strain)]
    ec_med_TC_prot = ec_med_TC_prot + prot_mutRed.df$`T>C`[which(prot_mutRed.df$species == strain)]
    ec_med_GA_prot = ec_med_GA_prot + prot_mutRed.df$`G>A`[which(prot_mutRed.df$species == strain)]
    ec_med_GT_prot = ec_med_GT_prot + prot_mutRed.df$`G>T`[which(prot_mutRed.df$species == strain)]
    ec_med_GC_prot = ec_med_GC_prot + prot_mutRed.df$`G>C`[which(prot_mutRed.df$species == strain)]
    ec_med_CA_prot = ec_med_CA_prot + prot_mutRed.df$`C>A`[which(prot_mutRed.df$species == strain)]
    ec_med_CT_prot = ec_med_CT_prot + prot_mutRed.df$`C>T`[which(prot_mutRed.df$species == strain)]
    ec_med_CG_prot = ec_med_CG_prot + prot_mutRed.df$`C>G`[which(prot_mutRed.df$species == strain)]
  }
  
  if(expLvl == "high") {
    ec_high_AT_prot = ec_high_AT_prot + prot_mutRed.df$`A>T`[which(prot_mutRed.df$species == strain)]
    ec_high_AG_prot = ec_high_AG_prot + prot_mutRed.df$`A>G`[which(prot_mutRed.df$species == strain)]
    ec_high_AC_prot = ec_high_AC_prot + prot_mutRed.df$`A>C`[which(prot_mutRed.df$species == strain)]
    ec_high_TA_prot = ec_high_TA_prot + prot_mutRed.df$`T>A`[which(prot_mutRed.df$species == strain)]
    ec_high_TG_prot = ec_high_TG_prot + prot_mutRed.df$`T>G`[which(prot_mutRed.df$species == strain)]
    ec_high_TC_prot = ec_high_TC_prot + prot_mutRed.df$`T>C`[which(prot_mutRed.df$species == strain)]
    ec_high_GA_prot = ec_high_GA_prot + prot_mutRed.df$`G>A`[which(prot_mutRed.df$species == strain)]
    ec_high_GT_prot = ec_high_GT_prot + prot_mutRed.df$`G>T`[which(prot_mutRed.df$species == strain)]
    ec_high_GC_prot = ec_high_GC_prot + prot_mutRed.df$`G>C`[which(prot_mutRed.df$species == strain)]
    ec_high_CA_prot = ec_high_CA_prot + prot_mutRed.df$`C>A`[which(prot_mutRed.df$species == strain)]
    ec_high_CT_prot = ec_high_CT_prot + prot_mutRed.df$`C>T`[which(prot_mutRed.df$species == strain)]
    ec_high_CG_prot = ec_high_CG_prot + prot_mutRed.df$`C>G`[which(prot_mutRed.df$species == strain)]
  }
}

ec_low_comb_prot = c("Escherichia_coli_low", ec_low_AT_prot, ec_low_AG_prot, ec_low_AC_prot, ec_low_TA_prot, ec_low_TG_prot, ec_low_TC_prot, 
                     ec_low_GA_prot, ec_low_GT_prot, ec_low_GC_prot, ec_low_CA_prot, ec_low_CT_prot, ec_low_CG_prot)

ec_med_comb_prot = c("Escherichia_coli_med", ec_med_AT_prot, ec_med_AG_prot, ec_med_AC_prot, ec_med_TA_prot, ec_med_TG_prot, ec_med_TC_prot, 
                     ec_med_GA_prot, ec_med_GT_prot, ec_med_GC_prot, ec_med_CA_prot, ec_med_CT_prot, ec_med_CG_prot)

ec_high_comb_prot = c("Escherichia_coli_high", ec_high_AT_prot, ec_high_AG_prot, ec_high_AC_prot, ec_high_TA_prot, ec_high_TG_prot, ec_high_TC_prot, 
                      ec_high_GA_prot, ec_high_GT_prot, ec_high_GC_prot, ec_high_CA_prot, ec_high_CT_prot, ec_high_CG_prot)

# Adding combined data to the tables
prot_mutRed.df = rbind(prot_mutRed.df, ec_low_comb_prot)
prot_mutRed.df = rbind(prot_mutRed.df, ec_high_comb_prot)
prot_mutRed.df = rbind(prot_mutRed.df, ec_med_comb_prot)


ec_AT_tRNA = 0
ec_AG_tRNA = 0
ec_AC_tRNA = 0
ec_TA_tRNA = 0
ec_TG_tRNA = 0
ec_TC_tRNA = 0
ec_GA_tRNA = 0
ec_GT_tRNA = 0
ec_GC_tRNA = 0
ec_CA_tRNA = 0
ec_CT_tRNA = 0
ec_CG_tRNA = 0

ec_strain = c("Escherichia_coli_mg1655", "Escherichia_coli_atcc")

for(strain in ec_strain) {
  ec_AT_tRNA = ec_AT_tRNA + tRNA_mutRed.df$`A>T`[which(tRNA_mutRed.df$species == strain)]
  ec_AG_tRNA = ec_AG_tRNA + tRNA_mutRed.df$`A>G`[which(tRNA_mutRed.df$species == strain)]
  ec_AC_tRNA = ec_AC_tRNA + tRNA_mutRed.df$`A>C`[which(tRNA_mutRed.df$species == strain)]
  ec_TA_tRNA = ec_TA_tRNA + tRNA_mutRed.df$`T>A`[which(tRNA_mutRed.df$species == strain)]
  ec_TG_tRNA = ec_TG_tRNA + tRNA_mutRed.df$`T>G`[which(tRNA_mutRed.df$species == strain)]
  ec_TC_tRNA = ec_TC_tRNA + tRNA_mutRed.df$`T>C`[which(tRNA_mutRed.df$species == strain)]
  ec_GA_tRNA = ec_GA_tRNA + tRNA_mutRed.df$`G>A`[which(tRNA_mutRed.df$species == strain)]
  ec_GT_tRNA = ec_GT_tRNA + tRNA_mutRed.df$`G>T`[which(tRNA_mutRed.df$species == strain)]
  ec_GC_tRNA = ec_GC_tRNA + tRNA_mutRed.df$`G>C`[which(tRNA_mutRed.df$species == strain)]
  ec_CA_tRNA = ec_CA_tRNA + tRNA_mutRed.df$`C>A`[which(tRNA_mutRed.df$species == strain)]
  ec_CT_tRNA = ec_CT_tRNA + tRNA_mutRed.df$`C>T`[which(tRNA_mutRed.df$species == strain)]
  ec_CG_tRNA = ec_CG_tRNA + tRNA_mutRed.df$`C>G`[which(tRNA_mutRed.df$species == strain)]
}

ec_comb_tRNA = c("Escherichia_coli", ec_AT_tRNA, ec_AG_tRNA, ec_AC_tRNA, ec_TA_tRNA, ec_TG_tRNA, ec_TC_tRNA, ec_GA_tRNA, ec_GT_tRNA, 
                 ec_GC_tRNA, ec_CA_tRNA, ec_CT_tRNA, ec_CG_tRNA)

# Adding combined data to the tables
tRNA_mutRed.df = rbind(tRNA_mutRed.df, ec_comb_tRNA)

# Reading in the table of nucleotide numbers in tRNA genes and the rest of the genome
nt_III.df = read.table("ntFreq_tRNA.tab", h = T)
nt_prot.df = read.table("ntFreq_prot.tab", h = T)

# Removing the _prot string from species names in nt_prot.df & nt_rest_prot.df
for(i in 1:nrow(nt_prot.df)) {
  tmp = str_remove(nt_prot.df$species[i], "_prot")
  nt_prot.df$species[i] = tmp
}

# Making sure the tables only have the species we're interested in
nt_prot.df = nt_prot.df[which(nt_prot.df$species %in% prot_sp),]

freqA_low = 0
freqT_low = 0
freqG_low = 0
freqC_low = 0

freqA_med = 0
freqT_med = 0
freqG_med = 0
freqC_med = 0

freqA_high = 0
freqT_high = 0
freqG_high = 0
freqC_high = 0

ecat_list = c("Escherichia_coli_low_atcc", "Escherichia_coli_high_atcc", "Escherichia_coli_med_atcc")
ecmg_list = c("Escherichia_coli_low_mg1655", "Escherichia_coli_high_mg1655", "Escherichia_coli_med_mg1655")

ec_list = append(ecat_list, ecmg_list)

for(strain in ec_list) {
  expLvl = unlist(strsplit(strain, "_"))[3]
  
  if(expLvl == "low") {
    freqA_low = freqA_low + nt_prot.df$freq_A[which(nt_prot.df$species == strain)]
    freqT_low = freqT_low + nt_prot.df$freq_T[which(nt_prot.df$species == strain)]
    freqG_low = freqG_low + nt_prot.df$freq_G[which(nt_prot.df$species == strain)]
    freqC_low = freqC_low + nt_prot.df$freq_C[which(nt_prot.df$species == strain)]
  }
  
  if(expLvl == "med") {
    freqA_med = freqA_med + nt_prot.df$freq_A[which(nt_prot.df$species == strain)]
    freqT_med = freqT_med + nt_prot.df$freq_T[which(nt_prot.df$species == strain)]
    freqG_med = freqG_med + nt_prot.df$freq_G[which(nt_prot.df$species == strain)]
    freqC_med = freqC_med + nt_prot.df$freq_C[which(nt_prot.df$species == strain)]
  }
  
  if(expLvl == "high") {
    freqA_high = freqA_high + nt_prot.df$freq_A[which(nt_prot.df$species == strain)]
    freqT_high = freqT_high + nt_prot.df$freq_T[which(nt_prot.df$species == strain)]
    freqG_high = freqG_high + nt_prot.df$freq_G[which(nt_prot.df$species == strain)]
    freqC_high = freqC_high + nt_prot.df$freq_C[which(nt_prot.df$species == strain)]
  }
}

ec_low = c("Escherichia_coli_low", freqA_low, freqT_low, freqG_low, freqC_low)
ec_med = c("Escherichia_coli_med", freqA_med, freqT_med, freqG_med, freqC_med)
ec_high = c("Escherichia_coli_high", freqA_high, freqT_high, freqG_high, freqC_high)

nt_prot.df = rbind(nt_prot.df, ec_low, ec_high, ec_med)

# Reducing to just combined data
spNew = c("Saccharomyces_cerevisiae_low", "Saccharomyces_cerevisiae_med", "Saccharomyces_cerevisiae_high", "Escherichia_coli_low", 
          "Escherichia_coli_med", "Escherichia_coli_high")

prot_mutRed.df = prot_mutRed.df[which(prot_mutRed.df$species %in% spNew),]
nt_prot.df = nt_prot.df[which(nt_prot.df$species %in% spNew),]

##########################################################################################################

ec_A_III = 0
ec_T_III = 0
ec_G_III = 0
ec_C_III = 0

# Making sure numbers are numerics in table
nt_III.df$freq_A = as.numeric(nt_III.df$freq_A)
nt_III.df$freq_T = as.numeric(nt_III.df$freq_T)
nt_III.df$freq_G = as.numeric(nt_III.df$freq_G)
nt_III.df$freq_C = as.numeric(nt_III.df$freq_C)

for(strain in ec_strain) {
  ec_A_III = ec_A_III + nt_III.df$freq_A[which(nt_III.df$species == strain)]
  ec_T_III = ec_T_III + nt_III.df$freq_T[which(nt_III.df$species == strain)]
  ec_G_III = ec_G_III + nt_III.df$freq_G[which(nt_III.df$species == strain)]
  ec_C_III = ec_C_III + nt_III.df$freq_C[which(nt_III.df$species == strain)]
}

ec_nt_III = c("Escherichia_coli", ec_A_III, ec_T_III, ec_G_III, ec_C_III)

nt_III.df = rbind(nt_III.df, ec_nt_III)

# Reducing to just combined data
spNew = c("Saccharomyces_cerevisiae", "Escherichia_coli")

tRNA_mutRed.df = tRNA_mutRed.df[which(tRNA_mutRed.df$species %in% spNew),]
nt_III.df = nt_III.df[which(nt_III.df$species %in% spNew),]

##########################################################################################
# Combining all expression levels to get just protein-coding - no bins
ec_list = c("Escherichia_coli_low", "Escherichia_coli_high", "Escherichia_coli_med")
sc_list = c("Saccharomyces_cerevisiae_low", "Saccharomyces_cerevisiae_high", "Saccharomyces_cerevisiae_med")

# Making sure all numeric columns are actually numeric
prot_mutRed.df[,2:13] = sapply(prot_mutRed.df[,2:13], as.numeric)
nt_prot.df[,2:5] = sapply(nt_prot.df[,2:5], as.numeric)

# Combining for nt_prot.df
scA = sum(nt_prot.df$freq_A[which(nt_prot.df$species %in% sc_list)])
scT = sum(nt_prot.df$freq_T[which(nt_prot.df$species %in% sc_list)])
scG = sum(nt_prot.df$freq_G[which(nt_prot.df$species %in% sc_list)])
scC = sum(nt_prot.df$freq_C[which(nt_prot.df$species %in% sc_list)])

ecA = sum(nt_prot.df$freq_A[which(nt_prot.df$species %in% ec_list)])
ecT = sum(nt_prot.df$freq_T[which(nt_prot.df$species %in% ec_list)])
ecG = sum(nt_prot.df$freq_G[which(nt_prot.df$species %in% ec_list)])
ecC = sum(nt_prot.df$freq_C[which(nt_prot.df$species %in% ec_list)])

scTot = c("Saccharomyces_cerevisiae", scA, scT, scG, scC)
ecTot = c("Escherichia_coli", ecA, ecT, ecG, ecC)

nt_prot_tot.df = rbind(scTot, ecTot)
colnames(nt_prot_tot.df) = colnames(nt_prot.df)

nt_prot_tot.df = as.data.frame(nt_prot_tot.df)
nt_prot_tot.df[,2:5] = sapply(nt_prot_tot.df[,2:5], as.numeric)

# Combining for prot_mutRed.df
prot_mutRed.df[7,1] = "Saccharomyces_cerevisiae"
prot_mutRed.df[8,1] = "Escherichia_coli"

for(i in 2:13) {
  prot_mutRed.df[7,i] = sum(prot_mutRed.df[which(prot_mutRed.df$species %in% sc_list),i])
  prot_mutRed.df[8,i] = sum(prot_mutRed.df[which(prot_mutRed.df$species %in% ec_list),i])
}

prot_mutRed_tot.df = prot_mutRed.df[7:8,]

# Removing new rows from original table in case it is needed again
prot_mutRed.df = prot_mutRed.df[-c(7:8),]

##########################################################################################
# Getting mutations/total # nucleotide

# Making copies of the tables for ntMutPlot graphs
prot_nt.df = prot_mutRed_tot.df

# Columns by ref base
Amut = c("A>T", "A>G", "A>C")
Tmut = c("T>A", "T>G", "T>C")
Gmut = c("G>A", "G>T", "G>C")
Cmut = c("C>A", "C>T", "C>G")

# Dividing number of mutations by number of reference base occurrences in genes of interest
# Replacing number of mutations in table with new number calculated
spNew = c("Saccharomyces_cerevisiae", "Escherichia_coli")
for(sp in spNew) {
  fA_prot = as.numeric(nt_prot_tot.df$freq_A[which(nt_prot_tot.df$species == sp)])
  fT_prot = as.numeric(nt_prot_tot.df$freq_T[which(nt_prot_tot.df$species == sp)])
  fG_prot = as.numeric(nt_prot_tot.df$freq_G[which(nt_prot_tot.df$species == sp)])
  fC_prot = as.numeric(nt_prot_tot.df$freq_C[which(nt_prot_tot.df$species == sp)])
  
  for(a in Amut) {
    tmp = as.numeric(prot_nt.df[which(prot_nt.df$species == sp), a])
    newA = tmp/fA_prot
    prot_nt.df[which(prot_nt.df$species == sp), a] = newA
  }
  
  for(t in Tmut) {
    tmp = as.numeric(prot_nt.df[which(prot_nt.df$species == sp), t])
    newT = tmp/fT_prot
    prot_nt.df[which(prot_nt.df$species == sp), t] = newT
  }
  
  for(g in Gmut) {
    tmp = as.numeric(prot_nt.df[which(prot_nt.df$species == sp), g])
    newG = tmp/fG_prot
    prot_nt.df[which(prot_nt.df$species == sp), g] = newG
  }
  
  for(c in Cmut) {
    tmp = as.numeric(prot_nt.df[which(prot_nt.df$species == sp), c])
    newC = tmp/fC_prot
    prot_nt.df[which(prot_nt.df$species == sp), c] = newC
  }
}


# Making copies of the tables for ntMutPlot graphs
tRNA_nt.df = tRNA_mutRed.df

# Columns by ref base
Amut = c("A>T", "A>G", "A>C")
Tmut = c("T>A", "T>G", "T>C")
Gmut = c("G>A", "G>T", "G>C")
Cmut = c("C>A", "C>T", "C>G")

# Dividing number of mutations by number of reference base occurrences in genes of interest
# Replacing number of mutations in table with new number calculated
for(sp in spNew) {
  fA_III = as.numeric(nt_III.df$freq_A[which(nt_III.df$species == sp)])
  fT_III = as.numeric(nt_III.df$freq_T[which(nt_III.df$species == sp)])
  fG_III = as.numeric(nt_III.df$freq_G[which(nt_III.df$species == sp)])
  fC_III = as.numeric(nt_III.df$freq_C[which(nt_III.df$species == sp)])
  
  for(a in Amut) {
    tmp = as.numeric(tRNA_nt.df[which(tRNA_nt.df$species == sp), a])
    newA = tmp/fA_III
    tRNA_nt.df[which(tRNA_nt.df$species == sp), a] = newA
  }
  
  for(t in Tmut) {
    tmp = as.numeric(tRNA_nt.df[which(tRNA_nt.df$species == sp), t])
    newT = tmp/fT_III
    tRNA_nt.df[which(tRNA_nt.df$species == sp), t] = newT
  }
  
  for(g in Gmut) {
    tmp = as.numeric(tRNA_nt.df[which(tRNA_nt.df$species == sp), g])
    newG = tmp/fG_III
    tRNA_nt.df[which(tRNA_nt.df$species == sp), g] = newG
  }
  
  for(c in Cmut) {
    tmp = as.numeric(tRNA_nt.df[which(tRNA_nt.df$species == sp), c])
    newC = tmp/fC_III
    tRNA_nt.df[which(tRNA_nt.df$species == sp), c] = newC
  }
}

##########################################################################################
# Making a table to hold the combined data

tRNA_mutRed.df$type = "tRNA"
tRNA_nt.df$type = "tRNA"

prot_mutRed_tot.df$type = "Protein Coding"
prot_nt.df$type = "Protein Coding"

mut_row = nrow(tRNA_mutRed.df) * 2
mut_full.df = data.frame(species = character(mut_row), v1 = numeric(mut_row), v2 = numeric(mut_row), v3 = numeric(mut_row), 
                         v4 = numeric(mut_row), v5 = numeric(mut_row), v6 = numeric(mut_row), v7 = numeric(mut_row), v8 = numeric(mut_row),
                         v9 = numeric(mut_row), v10 = numeric(mut_row), v11 = numeric(mut_row), v12 = numeric(mut_row),
                         type = character(mut_row))

mut_full_nt.df = data.frame(species = character(mut_row), v1 = numeric(mut_row), v2 = numeric(mut_row), v3 = numeric(mut_row), 
                            v4 = numeric(mut_row), v5 = numeric(mut_row), v6 = numeric(mut_row), v7 = numeric(mut_row), v8 = numeric(mut_row),
                            v9 = numeric(mut_row), v10 = numeric(mut_row), v11 = numeric(mut_row), v12 = numeric(mut_row),
                            type = character(mut_row))

colnames(mut_full.df) = c("species", "A>T", "A>G", "A>C", "T>A", "T>G", "T>C", "G>A", "G>T", "G>C", "C>A", "C>T", "C>G", "type")
colnames(mut_full_nt.df) = c("species", "A>T", "A>G", "A>C", "T>A", "T>G", "T>C", "G>A", "G>T", "G>C", "C>A", "C>T", "C>G", "type")

# Combining tRNA table with protein-coding table
rowCnt = 1
for(i in 1:mut_row) {
  if(mut_row/i >= 2) {
    mut_full.df[i,] = tRNA_mutRed.df[i,]
    mut_full_nt.df[i,] = tRNA_nt.df[i,]
  }
  
  else{
    mut_full.df[i,] = prot_mutRed_tot.df[rowCnt,]
    mut_full_nt.df[i,] = prot_nt.df[rowCnt,]
    rowCnt = rowCnt + 1
  }
}

df <- gather(mut_full.df, mut, mutNum, 2:13)
df$mutNum = as.numeric(df$mutNum)
df$fullType = paste(df$species, df$type, sep = "_")

df2 <- gather(mut_full_nt.df, mut, mutNum, 2:13)
df2$mutNum = as.numeric(df2$mutNum)
df2$fullType = paste(df2$species, df2$type, sep = "_")

typeList = unique(df$fullType)
typeNum = length(typeList)

# Adding the "total" column
df$total = 0
df2$total = 0
for(i in 1:typeNum) {
  tmpType = typeList[i]
  tmpSum = sum(df$mutNum[which(df$fullType == tmpType)])
  df$total[which(df$fullType == tmpType)] = tmpSum
  
  tmpSum2 = sum(df2$mutNum[which(df2$fullType == tmpType)])
  df2$total[which(df2$fullType == tmpType)] = tmpSum2
}

# Changing type to something more appropriate for the graph
df$type[which(df$type == "Protein Coding")] = "Protein-Coding"
df2$type[which(df2$type == "Protein Coding")] = "Protein-Coding"

# Getting the columns to display in the right order
typeLvls = c("tRNA", "Protein-Coding")
df$type = factor(df$type, levels = typeLvls)
df2$type = factor(df2$type, levels = typeLvls)

# Getting the BPSs to display in the right order
mutlvls = c("A>T", "A>G", "A>C", "T>A", "T>G", "T>C", "G>A", "G>T", "G>C", "C>A", "C>T", "C>G")
df$BPS = factor(df$mut, levels = mutlvls)
df2$BPS = factor(df2$mut, levels = mutlvls)

# Getting the standard deviation
df$mutProp = df$mutNum / df$total
df$mutPropSD = 0

df2$mutProp = df2$mutNum / df2$total
df2$mutPropSD = 0

for(i in 1:nrow(df)) {
  tmp = prop.test(df$mutNum[i], df$total[i])
  sd = df$mutProp[i] - tmp$conf.int[1]
  df$mutPropSD[i] = sd
  
  tmp2 = prop.test(df2$mutNum[i], df2$total[i])
  sd2 = df2$mutProp[i] - tmp2$conf.int[1]
  df2$mutPropSD[i] = sd2
}

################################################################################################
# Plots

# Number of specific mutation divided by total number of mutations
tRNA_prot_MutPlot <- ggplot(data = df, aes(x = mut, y = mutNum/total, fill = species, pattern = type)) +
  geom_bar_pattern(stat = "identity", position = position_dodge(preserve = "single"),
                   color = "black", 
                   pattern_fill = "black",
                   pattern_angle = 45,
                   pattern_density = 0.1,
                   pattern_spacing = 0.015,
                   pattern_key_scale_factor = 0.6) + 
  geom_errorbar(aes(ymin = mutProp-mutPropSD, ymax = mutProp+mutPropSD), position = "dodge") +
  #stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.8, position = position_dodge(0.9)) +
  scale_fill_manual(values = colorRampPalette(bpsPalette)(2), name = "Species", labels = c("E. coli", "S. cerevisiae")) +
  scale_pattern_manual(values = c(tRNA = "stripe", `Protein-Coding` = "none")) +
  labs(title = "Mutation Frequencies by Species", x = "Mutation", y = "Frequency", pattern = "Gene Type") + 
  guides(pattern = guide_legend(override.aes = list(fill = "white")),
         fill = guide_legend(override.aes = list(pattern = "none"))) +
  theme(legend.title = element_text(size=20), legend.text = element_text(size=20)) +
  theme_classic()

tRNA_prot_MutPlot    # <-- Saved as mutFreqTot.png locally
