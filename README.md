File information

check_lineup.R
==============
code for checking that the dataset and the assembly line up properly  
reads in table of mutations for MA line and genome  
makes sure base calls in MA line match those of ref genome  


get_tRNA.R
==========
gets only the tRNA data out of the gff3 file and writes to new file  


RNAPIII-mutations-functions.R
=============================
functions for analysis of datasets  


source_me.R
===========
runs analysis of datasets  


tpm.R
=====
finds highly & lowly expressed genes by tpm value in dataset  
only works for singular tRNAs - cannot distinguish tpm for multiple  
puts results in two files (high & low)  


tpm_abundance.R
===============
finds high/low expression cutoffs with tpm and tRNA abundance  
didn't use  


tpm_strain_conversion.R
=======================
finds coordinates of tRNAs in one strain from another  
didn't use after Zhang data was found  


vibrioFiles.R
=============
splits the two vibrio species into two separate datasets  
writes mutation files for each species and each condition  
