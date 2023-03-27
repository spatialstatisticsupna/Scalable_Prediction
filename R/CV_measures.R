rm(list=ls())
library(bigDM)
library(INLA)


#########################################
## Load the data and cartography files ##
#########################################
data("Data_LungCancer")
str(Data_LungCancer)

data("Carto_SpainMUN")
str(Carto_SpainMUN)


#########################################################
## Select the data corresponding to a certain province ##
#########################################################
Carto_SpainMUN$ID.prov <- substr(Carto_SpainMUN$ID,1,2)
unique(Carto_SpainMUN$ID.prov)

carto <- Carto_SpainMUN[Carto_SpainMUN$ID.prov=="01", ]
data <- Data_LungCancer[substr(Data_LungCancer$ID,1,2)=="01", ]


#####################################################
## Set to 'NA' the count data for prediction years ##
#####################################################
data$obs[data$year %in% 2013:2015] <- NA


####################################################
## Fit the models using the STCAR_INLA() function ##
####################################################
Model <- STCAR_INLA(carto=carto, data=data,
                    ID.area="ID", ID.year="year", O="obs", E="pop",
                    spatial="BYM2", temporal="rw1", interaction="TypeIV",
                    model="global",  strategy="gaussian")
summary(Model)


#######################################
## Compute cross validation measures ##
#######################################
help("inla.group.cv")

sdunif="expression:
          logdens=-log_precision/2;
          return(logdens)"

LOOCV <- inla.group.cv(result=Model, num.level.sets=-1)
LGOCV <- inla.group.cv(result=Model, num.level.sets=3)

c(LOOCV=-sum(log(LOOCV$cv), na.rm=T),
  LGOCV=-sum(log(LGOCV$cv), na.rm=T))
