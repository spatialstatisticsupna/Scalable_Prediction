rm(list=ls())
library(bigDM)
library(INLA)
library(tmap)


###################################################################
## Load male lung and overall cancer (all sites) mortality data  ##
###################################################################
data("Data_LungCancer")
str(Data_LungCancer)

load(url("https://emi-sstcdapp.unavarra.es/bigDM/inst/Rdata/Data_OverallCancer.Rdata"))
str(Data_OverallCancer)


#####################################################
## Load cartography file of Spanish municipalities ##
#####################################################
data("Carto_SpainMUN")
head(Carto_SpainMUN)

Carto_SpainMUN$ID.prov <- substr(Carto_SpainMUN$ID,1,2)
Carto.prov <- aggregate(Carto_SpainMUN[,"geometry"], list(ID.prov=Carto_SpainMUN$ID.prov), head)

Map <- tm_shape(Carto_SpainMUN) + tm_polygons() + 
  tm_shape(Carto.prov) + tm_borders(col="blue", lwd=2)
print(Map)


####################################################
## Fit the models using the STCAR_INLA() function ##
####################################################
help("STCAR_INLA")


## Select the data ##
cancer <- "Lung"  # cancer <- "Overall"

if(cancer=="Lung") data <- Data_LungCancer
if(cancer=="Overall") data <- Data_OverallCancer


## Set to 'NA' the count data for prediction years ##
data$obs[data$year %in% 2013:2015] <- NA


## Set the aproximation strategy for INLA ##
strategy <- "gaussian"
inla.mode <- "classic"   # inla.mode <- "compact" 


## Classical models (not run) 
##############################
global_typeI <- STCAR_INLA(carto=Carto_SpainMUN, data=data,
                           ID.area="ID", ID.year="year", O="obs", E="pop",
                           spatial="BYM2", temporal="rw1", interaction="TypeI",
                           model="global", strategy=strategy, inla.mode=inla.mode)

save(global_typeI, file="INLAmodels_global_typeI.Rdata")

global_typeIII <- STCAR_INLA(carto=Carto_SpainMUN, data=data,
                             ID.area="ID", ID.year="year", O="obs", E="pop",
                             spatial="BYM2", temporal="rw1", interaction="TypeIII",
                             model="global", strategy=strategy, inla.mode=inla.mode)

save(global_typeIII, file="INLAmodels_global_typeIII.Rdata")


## Disjoint models 
###################
workers <- future::availableWorkers()[-1]

type <- list(TypeI="TypeI", TypeII="TypeII", TypeIII="TypeIII", TypeIV="TypeIV")

MODELS.k0 <- lapply(type, function(x){
  STCAR_INLA(carto=Carto_SpainMUN, data=data, ID.group="ID.prov",
             ID.area="ID", ID.year="year", O="obs", E="pop",
             spatial="BYM2", temporal="rw1", interaction=x,
             model="partition", k=0, compute.fitted.values=TRUE,
             strategy=strategy, inla.mode=inla.mode,
             plan="cluster", workers=workers)
})

save(MODELS.k0, file="INLAmodels_k0.Rdata")


## 1st-order neighbourhood models 
##################################
workers <- future::availableWorkers()[-1]

type <- list(TypeI="TypeI", TypeII="TypeII", TypeIII="TypeIII", TypeIV="TypeIV")

MODELS.k1 <- lapply(type, function(x){
  STCAR_INLA(carto=Carto_SpainMUN, data=data, ID.group="ID.prov",
             ID.area="ID", ID.year="year", O="obs", E="pop",
             spatial="BYM2", temporal="rw1", interaction=x,
             model="partition", k=0, compute.fitted.values=TRUE,
             strategy=strategy, inla.mode=inla.mode,
             plan="cluster", workers=workers)
})

save(MODELS.k1, file="INLAmodels_k1.Rdata")


#####################################################################################################
## Table 3: Logarithmic score using both LOOCV and LGOCV techniques, model selection criteria and  ##
##          computational time (in minutes) for cancer mortality data with models fitted with INLA ##
#####################################################################################################
sdunif="expression:
          logdens=-log_precision/2;
          return(logdens)"

CV <- function(x, m=3){
  x$.args$inla.call <- NULL
  
  LOOCV <- inla.group.cv(result=x, num.level.sets=-1)
  LGOCV <- inla.group.cv(result=x, num.level.sets=m)
  
  data <- x$.args$data
  data$ID <- paste(data$Year, data$Area, sep=".")
  
  return(list(LOOCV=LOOCV, LGOCV=LGOCV, data=data))
}

compute.DIC <- function(x){
  data.frame(DIC=x$dic$dic,
             WAIC=x$waic$waic,
             Time=x$cpu.used[4]/60)
}

######################
## Classical models ##
######################
aux <- lapply(MODELS.k0, CV)
Table3a <- lapply(aux, function(x){
  data.frame(LOOCV=-sum(log(x$LOOCV$cv),na.rm=T),
             LGOCV=-sum(log(x$LGOCV$cv),na.rm=T))
})
Table3a <- do.call(rbind,Table3a)

Table3a$DIC <- unlist(lapply(MODELS, function(x) x$dic$dic))
Table3a$WAIC <- unlist(lapply(MODELS, function(x) x$waic$waic))
Table3a$Time <- unlist(lapply(MODELS, function(x) x$cpu.used[4]/60))
Table3a


#####################
## Disjoint models ##
#####################
Table3b <- data.frame(LOOCV=rep(NA,4), LGOCV=rep(NA,4),
                      row.names=paste("k0_type",c("I","II","III","IV"),sep=""))

## k0 - TypeI ##
load("Lung_Cancer/k0/temp/INLAsubmodels_202302192312.Rdata")
aux <- lapply(inla.models, CV)
rm(inla.models)
gc()

k0_typeI <- lapply(aux, function(x){
  data.frame(LOOCV=-sum(log(x$LOOCV$cv),na.rm=T),
             LGOCV=-sum(log(x$LGOCV$cv),na.rm=T))
})
Table3b["k0_typeI",] <- apply(do.call(rbind,k0_typeI),2,sum)


## k0 - TypeII ##
load("Lung_Cancer/k0/temp/INLAsubmodels_202302200424.Rdata")
aux <- lapply(inla.models, CV)
rm(inla.models)
gc()

k0_typeII <- lapply(aux, function(x){
  data.frame(LOOCV=-sum(log(x$LOOCV$cv),na.rm=T),
             LGOCV=-sum(log(x$LGOCV$cv),na.rm=T))
})
Table3b["k0_typeII",] <- apply(do.call(rbind,k0_typeII),2,sum)


## k0 - TypeIII ##
load("Lung_Cancer/k0/temp/INLAsubmodels_202302200524.Rdata")
aux <- lapply(inla.models, CV)
rm(inla.models)
gc()

k0_typeIII <- lapply(aux, function(x){
  data.frame(LOOCV=-sum(log(x$LOOCV$cv),na.rm=T),
             LGOCV=-sum(log(x$LGOCV$cv),na.rm=T))
})
Table3b["k0_typeIII",] <- apply(do.call(rbind,k0_typeIII),2,sum)


## k0 - TypeIV ##
load("Lung_Cancer/k0/temp/INLAsubmodels_202302201202.Rdata")
aux <- lapply(inla.models, CV)
rm(inla.models)
gc()

k0_typeIV <- lapply(aux, function(x){
  data.frame(LOOCV=-sum(log(x$LOOCV$cv),na.rm=T),
             LGOCV=-sum(log(x$LGOCV$cv),na.rm=T))
})
Table3b["k0_typeIV",] <- apply(do.call(rbind,k0_typeIV),2,sum)

