rm(list=ls())
library(bigDM)
library(INLA)
library(tmap)


#########################################
## Load the data and cartography files ##
#########################################
data("Data_LungCancer")
str(Data_LungCancer)

data("Carto_SpainMUN")
Carto_SpainMUN$ID.prov <- substr(Carto_SpainMUN$ID,1,2)
Carto.prov <- aggregate(Carto_SpainMUN[,"geometry"], list(ID.prov=Carto_SpainMUN$ID.prov), head)

Map <- tm_shape(Carto_SpainMUN) + tm_polygons() + 
  tm_shape(Carto.prov) + tm_borders(col="blue", lwd=2)
print(Map)


####################################################
## Fit the models using the STCAR_INLA() function ##
####################################################
help("STCAR_INLA")

## Set to 'NA' the count data for prediction years ##
data <- Data_LungCancer
data$obs[data$year %in% 2013:2015] <- NA


## Classical models (not run) 
##############################
global_typeI <- STCAR_INLA(carto=Carto_SpainMUN, data=data,
                           ID.area="ID", ID.year="year", O="obs", E="pop",
                           spatial="BYM2", temporal="rw1", interaction="TypeI",
                           model="global", strategy="gaussian")

global_typeIII <- STCAR_INLA(carto=Carto_SpainMUN, data=data,
                             ID.area="ID", ID.year="year", O="obs", E="pop",
                             spatial="BYM2", temporal="rw1", interaction="TypeIII",
                             model="global", strategy="gaussian")

## Disjoint models 
###################
workers <- future::availableWorkers()[-1]

k0_typeI <- STCAR_INLA(carto=Carto_SpainMUN, data=data, ID.group="ID.prov",
                       ID.area="ID", ID.year="year", O="obs", E="pop",
                       spatial="BYM2", temporal="rw1", interaction="TypeI",
                       model="partition", k=0, strategy="gaussian",
                       plan="cluster", workers=workers, inla.mode="compact")

k0_typeII <- STCAR_INLA(carto=Carto_SpainMUN, data=data, ID.group="ID.prov",
                        ID.area="ID", ID.year="year", O="obs", E="pop",
                        spatial="BYM2", temporal="rw1", interaction="TypeII",
                        model="partition", k=0, strategy="gaussian",
                        plan="cluster", workers=workers, inla.mode="compact")

k0_typeIII <- STCAR_INLA(carto=Carto_SpainMUN, data=data, ID.group="ID.prov",
                         ID.area="ID", ID.year="year", O="obs", E="pop",
                         spatial="BYM2", temporal="rw1", interaction="TypeIII",
                         model="partition", k=0, strategy="gaussian",
                         plan="cluster", workers=workers, inla.mode="compact")

k0_typeIV <- STCAR_INLA(carto=Carto_SpainMUN, data=data, ID.group="ID.prov",
                         ID.area="ID", ID.year="year", O="obs", E="pop",
                         spatial="BYM2", temporal="rw1", interaction="TypeIV",
                         model="partition", k=0, strategy="gaussian",
                         plan="cluster", workers=workers, inla.mode="compact")


## 1st-order neighbourhood models 
##################################
workers <- future::availableWorkers()-1

k1_typeI <- STCAR_INLA(carto=Carto_SpainMUN, data=data, ID.group="ID.prov",
                       ID.area="ID", ID.year="year", O="obs", E="pop",
                       spatial="BYM2", temporal="rw1", interaction="TypeI",
                       model="partition", k=1, strategy="gaussian",
                       plan="cluster", workers=workers, inla.mode="compact")

k1_typeII <- STCAR_INLA(carto=Carto_SpainMUN, data=data, ID.group="ID.prov",
                        ID.area="ID", ID.year="year", O="obs", E="pop",
                        spatial="BYM2", temporal="rw1", interaction="TypeII",
                        model="partition", k=1, strategy="gaussian",
                        plan="cluster", workers=workers, inla.mode="compact")

k1_typeIII <- STCAR_INLA(carto=Carto_SpainMUN, data=data, ID.group="ID.prov",
                         ID.area="ID", ID.year="year", O="obs", E="pop",
                         spatial="BYM2", temporal="rw1", interaction="TypeIII",
                         model="partition", k=1, strategy="gaussian",
                         plan="cluster", workers=workers, inla.mode="compact")

k1_typeIV <- STCAR_INLA(carto=Carto_SpainMUN, data=data, ID.group="ID.prov",
                        ID.area="ID", ID.year="year", O="obs", E="pop",
                        spatial="BYM2", temporal="rw1", interaction="TypeIV",
                        model="partition", k=1, strategy="gaussian",
                        plan="cluster", workers=workers, inla.mode="compact")
