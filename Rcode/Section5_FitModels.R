###################################################################################
## R script to reproduce the illustration section of Orozco-Acosta et al. (2023) ##
###################################################################################

## We recommend to install the latest version of 'bigDM' package
# devtools::install_github("spatialstatisticsupna/bigDM")
library(bigDM)
library(INLA)
library(sf)
library(tmap)

setwd(dirname(rstudioapi::getSourceEditorContext()$path))


######################################################################
## 1) Load male lung and overall cancer (all sites) mortality data  ##
######################################################################
load("../data/Data_LungCancer.Rdata")
str(Data_LungCancer)

load("../data/Data_OverallCancer.Rdata")
str(Data_OverallCancer)


## Set to 'NA' the count data for prediction years ##
Data_LungCancer$obs[Data_LungCancer$year %in% 2013:2015] <- NA
Data_OverallCancer$obs[Data_OverallCancer$year %in% 2013:2015] <- NA


########################################################
## 2) Load cartography file of Spanish municipalities ##
########################################################
data("Carto_SpainMUN")
head(Carto_SpainMUN)

Carto_SpainMUN$ID.prov <- substr(Carto_SpainMUN$ID,1,2)
Carto.prov <- aggregate(Carto_SpainMUN[,"geometry"], list(ID.prov=Carto_SpainMUN$ID.prov), head)

Map <- tm_shape(Carto_SpainMUN) + tm_polygons() + 
  tm_shape(Carto.prov) + tm_borders(col="blue", lwd=2)
print(Map)


#######################################################
## 3) Fit the models using the STCAR_INLA() function ##
#######################################################
help("STCAR_INLA")

## Select the data ##
cancer <- "Lung"  # cancer <- "Overall"

if(cancer=="Lung") data <- Data_LungCancer
if(cancer=="Overall") data <- Data_OverallCancer


## Additional INLA parameters (faster computations) ##
## To use the same strategy of the paper set: inla.mode <- "classic"
inla.mode <- "compact"


## CAUTION: These computations are very time consuming! ##

## Classical models 
#####################
Classical_typeI <- bigDM::STCAR_INLA(carto=Carto_SpainMUN, data=data,
                                     ID.area="ID", ID.year="year", O="obs", E="pop",
                                     spatial="BYM2", temporal="rw1", interaction="TypeI",
                                     model="global", compute.fitted.values=TRUE, inla.mode=inla.mode)

Classical_typeIII <- bigDM::STCAR_INLA(carto=Carto_SpainMUN, data=data,
                                       ID.area="ID", ID.year="year", O="obs", E="pop",
                                       spatial="BYM2", temporal="rw1", interaction="TypeIII",
                                       model="global", compute.fitted.values=TRUE, inla.mode=inla.mode)


## Disjoint models 
###################
workers <- future::availableWorkers()[-1]

type <- list(TypeI="TypeI", TypeII="TypeII", TypeIII="TypeIII", TypeIV="TypeIV")

MODELS.k0 <- lapply(type, function(x){
  res <- STCAR_INLA(carto=Carto_SpainMUN, data=data, ID.group="ID.prov",
                    ID.area="ID", ID.year="year", O="obs", E="pop",
                    spatial="BYM2", temporal="rw1", interaction=x,
                    model="partition", k=0, compute.fitted.values=TRUE,
                    inla.mode=inla.mode, plan="cluster", workers=workers,
                    save.models=TRUE)
  
  aux <- file.info(list.files(paste0(getwd(),"/temp"), full.names=T))
  filename.old <- rownames(aux)[order(aux$mtime, decreasing=T)[1]]
  filename.new <- paste0(getwd(),"/temp/INLAsubmodels_k0_",x,".Rdata")
  file.rename(filename.old, filename.new)
  
  return(res)
})


## 1st-order neighbourhood models 
##################################
workers <- future::availableWorkers()[-1]

type <- list(TypeI="TypeI", TypeII="TypeII", TypeIII="TypeIII", TypeIV="TypeIV")

MODELS.k1 <- lapply(type, function(x){
  res <- STCAR_INLA(carto=Carto_SpainMUN, data=data, ID.group="ID.prov",
                    ID.area="ID", ID.year="year", O="obs", E="pop",
                    spatial="BYM2", temporal="rw1", interaction=x,
                    model="partition", k=1, compute.fitted.values=TRUE,
                    inla.mode=inla.mode, plan="cluster", workers=workers,
                    save.models=TRUE)
  
  aux <- file.info(list.files(paste0(getwd(),"/temp"), full.names=T))
  filename.old <- rownames(aux)[order(aux$mtime, decreasing=T)[1]]
  filename.new <- paste0(getwd(),"/temp/INLAsubmodels_k1_",x,".Rdata")
  file.rename(filename.old, filename.new)
  
  return(res)
})


## Save the final model (1st-order neighbourhood + Type IV interaction) ##
k1_TypeIV <- MODELS.k1$TypeIV
save(k1_TypeIV, file=paste("../results/INLAmodels_",cancer,"Cancer.Rdata",sep=""))


########################
## 4) Compute Table 3 ##
########################
source("Auxiliary_functions.R")


## CAUTION: These computations are very time consuming! ##

## Classical models 
#####################
MODELS <- list(TypeI=Classical_typeI, TypeIII=Classical_typeIII)

aux <- lapply(MODELS, CV)
Table3a <- lapply(aux, function(x){
  data.frame(LOOCV=-sum(log(x$LOOCV$cv),na.rm=T),
             LGOCV=-sum(log(x$LGOCV$cv),na.rm=T))
})
Table3a <- do.call(rbind,Table3a)

Table3a$DIC <- unlist(lapply(MODELS, function(x) x$dic$dic))
Table3a$WAIC <- unlist(lapply(MODELS, function(x) x$waic$waic))
Table3a$Time <- unlist(lapply(MODELS, function(x) x$cpu.used[4]/60))
Table3a


## Disjoint models 
###################
Table3b <- data.frame(LOOCV=rep(NA,4), LGOCV=rep(NA,4),
                      row.names=paste("k0_Type",c("I","II","III","IV"),sep=""))

## k0 - TypeI ##
load("temp/INLAsubmodels_k0_TypeI.Rdata")
aux <- lapply(inla.models, CV)
rm(inla.models)
gc()

k0_TypeI <- lapply(aux, function(x){
  data.frame(LOOCV=-sum(log(x$LOOCV$cv),na.rm=T),
             LGOCV=-sum(log(x$LGOCV$cv),na.rm=T))
})
Table3b["k0_TypeI",] <- apply(do.call(rbind,k0_TypeI),2,sum)



## k0 - TypeII ##
load("temp/INLAsubmodels_k0_TypeII.Rdata")
aux <- lapply(inla.models, CV)
rm(inla.models)
gc()

k0_TypeII <- lapply(aux, function(x){
  data.frame(LOOCV=-sum(log(x$LOOCV$cv),na.rm=T),
             LGOCV=-sum(log(x$LGOCV$cv),na.rm=T))
})
Table3b["k0_TypeII",] <- apply(do.call(rbind,k0_TypeII),2,sum)


## k0 - TypeIII ##
load("temp/INLAsubmodels_k0_TypeIII.Rdata")
aux <- lapply(inla.models, CV)
rm(inla.models)
gc()

k0_TypeIII <- lapply(aux, function(x){
  data.frame(LOOCV=-sum(log(x$LOOCV$cv),na.rm=T),
             LGOCV=-sum(log(x$LGOCV$cv),na.rm=T))
})
Table3b["k0_TypeIII",] <- apply(do.call(rbind,k0_TypeIII),2,sum)


## k0 - TypeIV ##
load("temp/INLAsubmodels_k0_TypeIV.Rdata")
aux <- lapply(inla.models, CV)
rm(inla.models)
gc()

k0_TypeIV <- lapply(aux, function(x){
  data.frame(LOOCV=-sum(log(x$LOOCV$cv),na.rm=T),
             LGOCV=-sum(log(x$LGOCV$cv),na.rm=T))
})
Table3b["k0_TypeIV",] <- apply(do.call(rbind,k0_TypeIV),2,sum)


## DIC/WAIC values ##
Table3b$DIC <- unlist(lapply(MODELS.k0, function(x) x$dic$dic))
Table3b$WAIC <- unlist(lapply(MODELS.k0, function(x) x$waic$waic))
Table3b$Time <- unlist(lapply(MODELS.k0, function(x) x$cpu.used[3]/60))
Table3b


## 1st-order neighbourhood models 
##################################
Table3c<- data.frame(LOOCV=rep(NA,4), LGOCV=rep(NA,4),
                     row.names=paste("k1_Type",c("I","II","III","IV"),sep=""))

data("Carto_SpainMUN")
Carto_SpainMUN$ID.prov <- substr(Carto_SpainMUN$ID, 1, 2)
ID.list <- lapply(unique(Carto_SpainMUN$ID.prov), function(x) Carto_SpainMUN[Carto_SpainMUN$ID.prov==x,]$ID)


## k1 - TypeI ##
load("temp/INLAsubmodels_k1_TypeI.Rdata")
aux <- lapply(inla.models, CV)
rm(inla.models)
gc()

k1_TypeI <- mapply(function(x,y){
  pos <- which(x$data$Area %in% y)
  data.frame(LOOCV=-sum(log(x$LOOCV$cv[pos]),na.rm=T),
             LGOCV=-sum(log(x$LGOCV$cv[pos]),na.rm=T))
}, x=aux, y=ID.list, SIMPLIFY=F)

Table3c["k1_TypeI",] <- apply(do.call(rbind,k1_TypeI),2,sum)


## k1 - TypeII ##
load("temp/INLAsubmodels_k1_TypeII.Rdata")
aux <- lapply(inla.models, CV)
rm(inla.models)
gc()

k1_TypeII <- mapply(function(x,y){
  pos <- which(x$data$Area %in% y)
  data.frame(LOOCV=-sum(log(x$LOOCV$cv[pos]),na.rm=T),
             LGOCV=-sum(log(x$LGOCV$cv[pos]),na.rm=T))
}, x=aux, y=ID.list, SIMPLIFY=F)

Table3c["k1_TypeII",] <- apply(do.call(rbind,k1_TypeII),2,sum)


## k1 - TypeIII ##
load("temp/INLAsubmodels_k1_TypeIII.Rdata")
aux <- lapply(inla.models, CV)
rm(inla.models)
gc()

k1_TypeIII <- mapply(function(x,y){
  pos <- which(x$data$Area %in% y)
  data.frame(LOOCV=-sum(log(x$LOOCV$cv[pos]),na.rm=T),
             LGOCV=-sum(log(x$LGOCV$cv[pos]),na.rm=T))
}, x=aux, y=ID.list, SIMPLIFY=F)

Table3c["k1_TypeIII",] <- apply(do.call(rbind,k1_TypeIII),2,sum)


## k1 - TypeIV ##
load("temp/INLAsubmodels_k1_TypeIV.Rdata")
aux <- lapply(inla.models, CV)
rm(inla.models)
gc()

k1_TypeIV <- mapply(function(x,y){
  pos <- which(x$data$Area %in% y)
  data.frame(LOOCV=-sum(log(x$LOOCV$cv[pos]),na.rm=T),
             LGOCV=-sum(log(x$LGOCV$cv[pos]),na.rm=T))
}, x=aux, y=ID.list, SIMPLIFY=F)

Table3c["k1_TypeIV",] <- apply(do.call(rbind,k1_TypeIV),2,sum)


## DIC/WAIC values ##
Table3c$DIC <- unlist(lapply(MODELS.k1, function(x) x$dic$dic))
Table3c$WAIC <- unlist(lapply(MODELS.k1, function(x) x$waic$waic))
Table3c$Time <- unlist(lapply(MODELS.k1, function(x) x$cpu.used[3]/60))
Table3c


#################
## Final table ##
#################
Table3 <- rbind(Table3a,Table3b,Table3c)
Table3$LOOCV <- round(Table3$LOOCV-min(Table3$LOOCV))
Table3$LGOCV <- round(Table3$LGOCV-min(Table3$LGOCV))
Table3$DIC <- round(Table3$DIC-min(Table3$DIC))
Table3$WAIC <- round(Table3$WAIC-min(Table3$WAIC))
Table3$Time <- round(Table3$Time)
print(Table3)

save(Table3, file="../results/Illustration_Table3.Rdata")
