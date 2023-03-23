# Getting the E. coli expression data for protein-coding genes

rm(list = ls())

library(genomation)

setwd("C:/Users/Sarah/Desktop/Lab/RNAPIII/RNAPIII-mutation-rate-main/Data/Ecoli/")

tmp = read.table("DB/E_coli_v4_Build_5_affy_cdf/avg_E_coli_v4_Build_5_affy_cdf_exps380probes4345.tab", h = T)

vMG1655 = tmp$MG1655_wt_24hr_biofilm
vGeneNames = tmp$E_coli_v4_Build_5_affy_cdf.genes

genes = cbind(vGeneNames, vMG1655)
genes.tab = as.data.frame(genes)

genes.tab$vGeneNames = as.character(genes.tab$vGeneNames)
genes.tab$vMG1655 = as.numeric(genes.tab$vMG1655)

# Fixing the gene names
for(i in 1:nrow(genes.tab)) {
  tmpL = unlist(strsplit(genes.tab$vGeneNames[i], "_"))
  
  if(length(tmpL[which(substr(tmpL,1,1) == "b")]) > 1) {
    tmpB = tmpL[which(substr(tmpL,1,1) == "b")]
    newName = tmpB[2]
  }
  
  else {
    newName = tmpL[which(substr(tmpL,1,1) == "b")]
  }
  
  if(length(newName) == 0) {
    newName = "NA"
  }
  
  genes.tab$vGeneNames[i] = newName
}

# Removing "NA" gene name lines
if(nrow(genes.tab[which(genes.tab$vGeneNames == "NA"),]) > 0) {
  genes.tab = genes.tab[which(genes.tab$vGeneNames != "NA"),]
}

# Removing "NA" gene name lines
genes.tab = genes.tab[which(genes.tab$vGeneNames != "NA"),]

# GFF files
gff.list = gffToGRanges("DB/NC_000913.2.gff3")
gffATCC.list = gffToGRanges("DB/CP000946.1.gff3")

gff_cut.list = gff.list[which(gff.list$locus_tag %in% genes.tab$vGeneNames & gff.list$type == "gene" & 
                                    gff.list$gene_biotype == "protein_coding")]

# Removing genes not in gff_cut.list
genes.tab = genes.tab[which(genes.tab$vGeneNames %in% gff_cut.list$locus_tag),]

# Getting start and end coordinates for genes
gff_cut.tab = as.data.frame(gff_cut.list)

genes.tab$start = 0
genes.tab$end = 0
for(i in 1:nrow(genes.tab)) {
  genes.tab$start[i] = gff_cut.tab$start[which(gff_cut.list$locus_tag == genes.tab$vGeneNames[i])]
  genes.tab$end[i] = gff_cut.tab$end[which(gff_cut.list$locus_tag == genes.tab$vGeneNames[i])]
}

genes.tab$chromosome = gff_cut.tab$seqnames[1]
genes.tab$strand = "*"

colOrder = c("chromosome", "")

############################################################################################################################################
# Sorting genes by expression level

genes.tab$logVal = log(genes.tab$vMG1655)
genes.tab = genes.tab[order(genes.tab$logVal),]

# Classifying genes as high/medium/low expression
per10 = round(nrow(genes.tab) * 0.1)

rowHigh = nrow(genes.tab) - per10 + 1

prot_low.tab = genes.tab[1:per10,]
prot_high.tab = genes.tab[rowHigh:nrow(genes.tab),]
prot_med.tab = genes.tab[(per10 + 1):(rowHigh - 1),]
