# Finding the tRNA coordinates from the MG1655 strain (Foster2015/2018) in the ATCC 8739 strain (Zhang2018)

# Only needs to be done once
##old_path = Sys.getenv()
##Sys.setenv(PATH = paste(old_path, "C:/Users/Sarah/Desktop/", sep = .Platform$path.sep))

rm(list = ls())

setwd("C:/Users/Sarah/Desktop/Lab/RNAPIII/RNAPIII-mutation-rate-main/Data/Ecoli/")

library(GenomicRanges)
library(BSgenome)
library(rBLAST)

fGenome = readDNAStringSet("DB/U00096.2.fasta")
zGenome = readDNAStringSet("DB/CP000946.1.fasta")

names(fGenome) = c("U00096.2")

fPos.tab = read.table("DB/tRNApos_low_NC.2.tab", h = T)
fPos.gr = makeGRangesFromDataFrame(fPos.tab)

seq = getSeq(fGenome, fPos.gr)
names(seq) = fPos.tab$name

# Making the BLAST database (only needs to be done once)
##makeblastdb(file = "DB/NC_000913.3.fasta", dbtype = "nucl")

# Loading the BLAST database
dbf = blast(db = "./DB/U00096.2.fasta")
dbz = blast(db = "./DB/CP000946.1.fasta")
##dbn = blast(db = "./DB/U00096.3.fasta")

# Blasting the sequences against the database
fRes = predict(dbf, seq)
fCleanRes = fRes[which(fRes$Perc.Ident == 100),]

fCleanRes$nb = 1
fvnb = by(fCleanRes$nb, fCleanRes$QueryID, sum)
length(which(fvnb == 1))

zRes = predict(dbz, seq)
zCleanRes = zRes[which(zRes$Perc.Ident == 100),]

zCleanRes$nb = 1
zvnb = by(zCleanRes$nb, zCleanRes$QueryID, sum)
length(which(zvnb == 1))

### 45 and 16 are getting cut from the Zhang dataset - should we do anything?

# Making a table of all the positions found through the search
tRNAATCCpos.tab = zCleanRes[,c(2,1,9,10)]
colnames(tRNAATCCpos.tab) = c("chromosome", "name", "start", "end")

#########################################################################################################################
# Separating high/low tpm from MG1655

# Reading in the names of the tRNAs we are using right now
tRNAmatched = readLines("DB/matched-tRNAs.tab")
tRNAlow = readLines("DB/matched-tRNAs-low.tab")
tRNAhigh = readLines("DB/matched-tRNAs-high.tab")

# Pulling only the results which correspond to a tRNA listed in tRNAmatched
newRes = zCleanRes[which(zCleanRes$QueryID %in% tRNAmatched),]

for(i in 1:nrow(newRes)) {
  if(newRes$S.start[i] > newRes$S.end[i]) {
    newStart = newRes$S.end[i]
    newEnd = newRes$S.start[i]
    
    newRes$S.start[i] = newStart
    newRes$S.end[i] = newEnd
  }
}

lowRes = newRes[which(newRes$QueryID %in% tRNAlow),]
highRes = newRes[which(newRes$QueryID %in% tRNAhigh),]

# Making the high/low tables
tpmLow = lowRes[,c(2,1,9,10)]
tpmHigh = highRes[,c(2,1,9,10)]

names(tpmLow) = c("chromosome", "geneType", "start", "end")
names(tpmHigh) = c("chromosome", "geneType", "start", "end")

tpmLow$geneType = "tRNA_gene"
tpmHigh$geneType = "tRNA_gene"

tpmLow$strand = "+"
tpmHigh$strand = "+"

#########################################################################################################################
# Writing the tables

write.table(tpmLow, file = "DB/tpmLowZhang.tab", quote = F, sep = '\t', row.names = F, col.names = T)
write.table(tpmHigh, file = "DB/tpmHighZhang.tab", quote = F, sep = '\t', row.names = F, col.names = T)

# Writing the tRNA-ATCC8739-pos.tab table
write.table(tRNAATCCpos.tab, file = "DB/ATCC8739-pos.tab", quote = F, sep = '\t', row.names = F, col.names = T)





#####
##matched = read.table("DB/matched-tRNAs.tab")

##nGenome = readDNAStringSet("DB/U00096.3.fasta")
##names(nGenome) = "U00096.3"

##nPos.tab = read.table("DB/comp.tab", h = T)
##nPos.tab = nPos.tab[,c(1,3)]
##colnames(nPos.tab) = c("start", "end")

##nPos.tab$chromosome = "U00096.3"
##nPos.tab$name = fPos.tab$name

##nPos.gr = makeGRangesFromDataFrame(nPos.tab)

##nSeq = getSeq(nGenome, nPos.gr)
##names(nSeq) = nPos.tab$name

##nRes = predict(dbn, nSeq)
##nCleanRes = nRes[which(nRes$Perc.Ident == 100),]

##nCleanRes$nb = 1
##nvnb = by(nCleanRes$nb, nCleanRes$QueryID, sum)
##length(which(nvnb == 1))


##nrow(fCleanRes[which(fCleanRes$QueryID %in% matched$V1),])
##nrow(nCleanRes[which(nCleanRes$QueryID %in% matched$V1),])
#####

#########################################################################################################################
# Comparing copy numbers of tRNA genes in MG1655 and ATCC 8739 strains

# Getting copy number of each tRNA gene in each database
##cleanNames = unique(fCleanRes$QueryID)

##fCopyNum = data.frame(matrix(ncol = 2, nrow = length(cleanNames)))
##colnames(fCopyNum) = c("name", "copyNum")
##for(i in 1:length(cleanNames)) {
##  id = cleanNames[i]
##  rnum = nrow(fCleanRes[which(fCleanRes$QueryID == id),])

##  fCopyNum$name[i] = id
##  fCopyNum$copyNum[i] = rnum
##}

##zCopyNum = data.frame(matrix(ncol = 2, nrow = length(cleanNames)))
##colnames(zCopyNum) = c("name", "copyNum")
##for(i in 1:length(cleanNames)) {
##  id = cleanNames[i]
##  rnum = nrow(zCleanRes[which(zCleanRes$QueryID == id),])

##  zCopyNum$name[i] = id
##  zCopyNum$copyNum[i] = rnum
##}

# Comparing copy numbers between the two
##copyComp = merge(fCopyNum, zCopyNum, by = "name")
##colnames(copyComp) = c("name", "fCopy", "zCopy")
##copyComp$diff = copyComp$fCopy - copyComp$zCopy

# Looking at copy numbers in lowRes and highRes
##lowNames = unique(lowRes$QueryID)
##highNames = unique(highRes$QueryID)

##lowComp = copyComp[which(copyComp$name %in% lowNames),]
##highComp = copyComp[which(copyComp$name %in% highNames),]

# Writing the copyNum.tab table
##write.table(copyComp, file = "DB/copyNum.tab", quote = F, sep = '\t', row.names = F, col.names = T)
