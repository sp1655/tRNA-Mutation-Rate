# Finding high/low expression cutoffs with tpm AND tRNA abundance

rm(list = ls())

library(readxl)

setwd("C:/Users/Sarah/Desktop/Lab/RNAPIII/RNAPIII-mutation-rate-main/Data/Ecoli/")

# Reading in the tpm data with the tRNA names
tpm.tab = read.table("DB/onetRNA.tab", h = T)

#########################################################################################################################
# Dong et al. data

dAbund.tab = read_xlsx("41598_2019_39369_MOESM2_ESM.xlsx", sheet = "Ecoli174.Dong et al.1996", range = "A3:I40")
dAbund.tab = as.data.frame(dAbund.tab)

# Checking that the plot is the same as in the excel sheet - it is!
plot(x = dAbund.tab$`average tRNA abundance`, y = dAbund.tab$`total tpm (Ecoli174)`)

#########################################################################################################################
# Ikemura data - we haven't done the tpm stuff for this yet

##iAbund.tab = read_xlsx("41598_2019_39369_MOESM2_ESM.xlsx", sheet = "Ecoli174.Ikrmura1985", range = "A3:G37")

# Removing the rows for which there is no abundance data
##abund.tab = abund.tab[which(is.na(abund.tab$`Ikemura 1985 adapted`) == F),]

# Checking that the plot is the same as in the excel sheet - it is!
##plot(x = abund.tab$`Ikemura 1985 adapted`, y = abund.tab$`tpm (Ecoli 174)`)
