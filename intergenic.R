#################################################################################
#                                                                               #
# This code generates GenomicRanges objects that contain the location of all    #
# transcribed and untranscribed regions in S. cerevisiae and E. coli            #
#                                                                               #
# For S. cerevisiae:                                                            #
# The transcribed regions are obtained by combining the location occupied by    #
# all type of genes in the ensembl gff3 annotation with a list of UTR locations #
# obtained from yeastmine (https://yeastmine.yeastgenome.org/)                  #
#                                                                               #
# For E. coli:                                                                  #
# Presence of UTRs is ignored, as they are not significant in this species for  #
# the study being conducted                                                     #
#                                                                               #
#################################################################################

rm(list = ls())

library(GenomicRanges)
library(Biostrings)
library(BSgenome)
library(genomation)

############################################################################################################################################
# S. cerevisiae intergenic data

setwd("C:/Users/Sarah/Desktop/Lab/RNAPIII/RNAPIII-mutation-rate-main/Data/Scerevisiae/DB/")

gff_file <- "Saccharomyces_cerevisiae.R64-1-1.84.gff3"
utr_file <- "UTRs.tsv.gz"

gff <- gffToGRanges(gff.file = gff_file)

# Excluding the mitochondria:
gff <- gff[which(seqnames(gff)!="Mito")]

# Selecting all types of genes:
genes <- gff[which( is.element(gff$type, c("gene", "ncRNA_gene", "tRNA_gene", "snoRNA_gene", "pseudogene", "snRNA_gene", "rRNA_gene"))==T )]
strand(genes) = "*"
# This list of chromosomes contains all positions in the genome:
chromosomes <- gff[which(gff$type == "chromosome")]

#################################################################################
# Loading the UTR annotation:
t_UTR <- read.table(utr_file, h=F, sep="\t")
colnames(t_UTR) <- c("utr_id", "chromosome", "end", "start", "strand")

# Chromosome names need to be modified to remove the leading "chr"
t_UTR$chromosome <- substr(t_UTR$chromosome, 4, nchar(t_UTR$chromosome))
# We are not going to use the strand information for this analyzis
t_UTR$strand = "*"
utr <- makeGRangesFromDataFrame(t_UTR, keep.extra.columns = F)

#################################################################################
# Combining the genic regions from the original gff3 annotation with the UTRs:

transcribed_regions <- reduce(c(genes, utr), ignore.strand = T)
untranscribed_regions <- reduce(setdiff(chromosomes, transcribed_regions, ignore.strand = T))

# Sanity check: the total size of transcribed + untranscribed regions should equal the total genome size
size_untranscribed <- sum(width(untranscribed_regions))
size_transcribed <- sum(width(transcribed_regions))
size_genome <- sum(width(chromosomes))

if( (size_transcribed + size_untranscribed) == size_genome ){
  cat("Sanity check passed!\n")
} else {
  cat("PROBLEM, the sizes do not match:\n")
  cat("Transcribed regions: ", size_transcribed, "\n")
  cat("Unranscribed regions: ", size_untranscribed, "\n")
  cat("Entire genome: ", size_genome, "\n")
}

#################################################################################
# Getting in RNAPIII-pos.tab format

intergenic.tab = as.data.frame(untranscribed_regions)
intergenic.tab = intergenic.tab[,-c(4)]
colnames(intergenic.tab) = c("chromosome", "start", "end", "strand")

intergenic.tab$geneType = "intergenic"
intergenic.tab = intergenic.tab[,c(1,5,2,3,4)]

#################################################################################
# Writing the table
write.table(intergenic.tab, file = "intergenic-pos.tab", quote = F, sep = '\t', row.names = F, col.names = T)

############################################################################################################################################
# E. coli intergenic regions

setwd("C:/Users/Sarah/Desktop/Lab/RNAPIII/RNAPIII-mutation-rate-main/Data/Ecoli/DB/")

#################################################################################
# ATCC 8739 (Zhang)

gff_atcc <- gffToGRanges("CP000946.1.gff3")

# Selecting all types of genes:
atcc_genes <- gff_atcc[which( is.element(gff_atcc$type, c("gene", "pseudogene"))==T )]
strand(atcc_genes) = "*"

# Getting all positions in the genome
atcc_chromosomes <- gff_atcc[which(gff_atcc$type == "region")]

# Transcribed regions
transcribed_regions_atcc <- reduce(atcc_genes, ignore.strand = T)

# Getting the untranscribed regions
untranscribed_regions_atcc <- reduce(setdiff(atcc_chromosomes, transcribed_regions_atcc, ignore.strand = T))

#################################################################################

# Sanity check: the total size of transcribed + untranscribed regions should equal the total genome size
size_untranscribed_atcc <- sum(width(untranscribed_regions_atcc))
size_transcribed_atcc <- sum(width(transcribed_regions_atcc))
size_genome_atcc <- sum(width(atcc_chromosomes))

if( (size_transcribed_atcc + size_untranscribed_atcc) == size_genome_atcc ){
  cat("Sanity check passed!\n")
} else {
  cat("PROBLEM, the sizes do not match:\n")
  cat("Transcribed regions: ", size_transcribed_atcc, "\n")
  cat("Unranscribed regions: ", size_untranscribed_atcc, "\n")
  cat("Entire genome: ", size_genome_atcc, "\n")
}

#################################################################################
# MG1655 (Foster)

gff_mg1655 <- gffToGRanges("U00096.2.gff3")

# Selecting all types of genes:
mg1655_genes <- gff_mg1655[which( is.element(gff_mg1655$type, c("gene", "pseudogene"))==T )]
strand(mg1655_genes) = "*"

# Getting all positions in the genome
mg1655_chromosomes <- gff_mg1655[which(gff_mg1655$type == "region")]

# Transcribed regions
transcribed_regions_mg1655 <- reduce(mg1655_genes, ignore.strand = T)

# Getting the untranscribed regions
untranscribed_regions_mg1655 <- reduce(setdiff(mg1655_chromosomes, transcribed_regions_mg1655, ignore.strand = T))

#################################################################################

# Sanity check: the total size of transcribed + untranscribed regions should equal the total genome size
size_untranscribed_mg1655 <- sum(width(untranscribed_regions_mg1655))
size_transcribed_mg1655 <- sum(width(transcribed_regions_mg1655))
size_genome_mg1655 <- sum(width(mg1655_chromosomes))

if( (size_transcribed_mg1655 + size_untranscribed_mg1655) == size_genome_mg1655 ){
  cat("Sanity check passed!\n")
} else {
  cat("PROBLEM, the sizes do not match:\n")
  cat("Transcribed regions: ", size_transcribed_mg1655, "\n")
  cat("Unranscribed regions: ", size_untranscribed_mg1655, "\n")
  cat("Entire genome: ", size_genome_mg1655, "\n")
}

#################################################################################
# Getting both in RNAPIII-pos.tab format

# ATCC 8739 (Zhang)
intergenic_atcc.tab = as.data.frame(untranscribed_regions_atcc)
intergenic_atcc.tab = intergenic_atcc.tab[,-c(4)]
colnames(intergenic_atcc.tab) = c("chromosome", "start", "end", "strand")

intergenic_atcc.tab$geneType = "intergenic"
intergenic_atcc.tab = intergenic_atcc.tab[,c(1,5,2,3,4)]

# MG1655 (Foster)
intergenic_mg1655.tab = as.data.frame(untranscribed_regions_mg1655)
intergenic_mg1655.tab = intergenic_mg1655.tab[,-c(4)]
colnames(intergenic_mg1655.tab) = c("chromosome", "start", "end", "strand")

intergenic_mg1655.tab$geneType = "intergenic"
intergenic_mg1655.tab = intergenic_mg1655.tab[,c(1,5,2,3,4)]

#################################################################################
# Writing tables

write.table(intergenic_atcc.tab, file = "atcc-intergenic-pos.tab", quote = F, sep = '\t', row.names = F, col.names = T)
write.table(intergenic_mg1655.tab, file = "mg1655-intergenic-pos.tab", quote = F, sep = '\t', row.names = F, col.names = T)
