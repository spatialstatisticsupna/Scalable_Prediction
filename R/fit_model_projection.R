rm(list = ls())
######################
## Packages loading ##
######################

# install.packages("bigDM")
# install.packages("INLA", repos=c(getOption("repos"), INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
library(bigDM)
library(INLA)

################################################################################
## data preparation for predicting from 2013 to 2015 with data from 1991-2012 ##
################################################################################

data <- read.table(file = "data/Data_LungCancer.txt", header = TRUE)
# data <- read.table(file = "data/Data_GlobalCancer.txt", header = TRUE)
Carto_SpainMUN$ID.prov <- substr(Carto_SpainMUN$ID, 1, 2)
pred <- 2013:2015
data$obs[data$year %in% pred] <- NA

####################
## fitting models ##
####################

## classical models 

global_typeI <- STCAR_INLA(carto = Carto_SpainMUN, data = data, ID.area = "ID",
                           ID.year = "year", O = "obs", E = "pob", spatial = "BYM2",
                           temporal = "rw1", interaction = "TypeI", model = "global")

global_typeIII <- STCAR_INLA(carto = Carto_SpainMUN, data = data, ID.area = "ID",
                             ID.year = "year", O = "obs", E = "pob", spatial = "BYM2",
                             temporal = "rw1", interaction = "TypeIII", model = "global")

## disjoint models

nodes <- rep("localhost", 2)

k0_typeI <- STCAR_INLA(carto=Carto_SpainMUN, data = data, ID.area="ID",
                       ID.year="year", O="obs", E="pob", ID.group="ID.prov",
                       spatial="BYM2", temporal="rw1", interaction="TypeI",
                       model="partition", k=0, merge.strategy = "original", 
                       plan = "cluster", workers=nodes, compute.DIC = FALSE, 
                       save.models = TRUE)

k0_typeII <- STCAR_INLA(carto=Carto_SpainMUN, data = data, ID.area="ID",
                        ID.year="year", O="obs", E="pob", ID.group="ID.prov",
                        spatial="BYM2", temporal="rw1", interaction="TypeII",
                        model="partition", k=0, merge.strategy = "original", 
                        plan = "cluster", workers=nodes, compute.DIC = FALSE, 
                        save.models = TRUE)

k0_typeIII <- STCAR_INLA(carto=Carto_SpainMUN, data = data, ID.area="ID",
                         ID.year="year", O="obs", E="pob", ID.group="ID.prov",
                         spatial="BYM2", temporal="rw1", interaction="TypeIII",
                         model="partition", k=0, merge.strategy = "original", 
                         plan = "cluster", workers=nodes, compute.DIC = FALSE, 
                         save.models = TRUE)

k0_typeIV <- STCAR_INLA(carto=Carto_SpainMUN, data = data, ID.area="ID",
                        ID.year="year", O="obs", E="pob", ID.group="ID.prov",
                        spatial="BYM2", temporal="rw1", interaction="TypeIV",
                        model="partition", k=0, merge.strategy = "original", 
                        plan = "cluster", workers=nodes, compute.DIC = FALSE, 
                        save.models = TRUE)

## 1st order neighbourhood models 

k1_typeI <- STCAR_INLA(carto=Carto_SpainMUN, data = data, ID.area="ID",
                       ID.year="year", O="obs", E="pob", ID.group="ID.prov",
                       spatial="BYM2", temporal="rw1", interaction="TypeI",
                       model="partition", k=1, merge.strategy = "original", 
                       plan = "cluster", workers=nodes, compute.DIC = FALSE, 
                       save.models = TRUE)

k1_typeII <- STCAR_INLA(carto=Carto_SpainMUN, data = data, ID.area="ID",
                        ID.year="year", O="obs", E="pob", ID.group="ID.prov",
                        spatial="BYM2", temporal="rw1", interaction="TypeII",
                        model="partition", k=1, merge.strategy = "original", 
                        plan = "cluster", workers=nodes, compute.DIC = FALSE, 
                        save.models = TRUE)

k1_typeIII <- STCAR_INLA(carto=Carto_SpainMUN, data = data, ID.area="ID",
                         ID.year="year", O="obs", E="pob", ID.group="ID.prov",
                         spatial="BYM2", temporal="rw1", interaction="TypeIII",
                         model="partition", k=1, merge.strategy = "original", 
                         plan = "cluster", workers=nodes, compute.DIC = FALSE, 
                         save.models = TRUE)

k1_typeIV <- STCAR_INLA(carto=Carto_SpainMUN, data = data, ID.area="ID",
                        ID.year="year", O="obs", E="pob", ID.group="ID.prov",
                        spatial="BYM2", temporal="rw1", interaction="TypeIV",
                        model="partition", k=1, merge.strategy = "original", 
                        plan = "cluster", workers=nodes, compute.DIC = FALSE, 
                        save.models = TRUE)

fitted_MODELS <- list(global_typeI, global_typeIII,
                      k0_typeI, k0_typeII, k0_typeIII, k0_typeIV,
                      k1_typeI, k1_typeII, k1_typeIII, k1_typeIV)

save(list=c("fitted_MODELS"), file = "data/fitted_MODELS.Rdata")

