rm(list=ls())
library(INLA)
library(bigDM)

################################
## Load cancer mortality data ##
################################
data("Data_LungCancer")
str(Data_LungCancer)


########################################################################################
## Load final model (1st-order neighbourhood + Type IV interaction) fitted using INLA ##
########################################################################################
