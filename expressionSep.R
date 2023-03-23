# Separate protein-coding genes in S. cerevisiae into three bins by expression level

rm(list = ls())

library(readxl)

setwd("C:/Users/Sarah/Desktop/Lab/RNAPIII/RNAPIII-mutation-rate-main/Data/Scerevisiae/DB/")

# Data source: https://doi.org/10.1371/journal.pone.0015442
data.xl <- read_excel("pone.0015442.s002.xls")

data.tab = data.xl[which(is.na(data.xl$`TR (mol/min) dilution corrected`) == F),]
data.tab = as.data.frame(data.tab)

# Getting only the relevant data
tmp.tab = data.tab[,c(1,2,8)]
tmp.tab$logTR = log(tmp.tab$`TR (mol/min) dilution corrected`)

# Ordering by logTR (ascending order)
new.tab = tmp.tab[order(tmp.tab$logTR),]
names(new.tab) = make.names(names(new.tab), unique = T)

# Getting 10% of rows for high and low split
rowTot = nrow(new.tab)
per10 = rowTot * 0.1
per10 = round(per10, 0)

per10Low = per10
per10High = rowTot - per10 + 1

low.tab = new.tab[1:per10Low,]
high.tab = new.tab[per10High:rowTot,]

# Getting everything in the middle
medLow = per10Low + 1
medHigh = per10High - 1

med.tab = new.tab[medLow:medHigh,]

# Writing the tables
write.table(low.tab, file = "lowProt.tab", quote = F, sep = '\t', row.names = F, col.names = T)
write.table(high.tab, file = "highProt.tab", quote = F, sep = '\t', row.names = F, col.names = T)
write.table(med.tab, file = "medProt.tab", quote = F, sep = '\t', row.names = F, col.names = T)
