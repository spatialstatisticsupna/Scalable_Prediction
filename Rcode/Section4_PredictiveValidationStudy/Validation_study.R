####################################################################################################
## R script to reproduce the predictive validation study describen in Orozco-Acosta et al., 2023) ##
####################################################################################################

## We recommend to install the latest version of 'bigDM' package
# devtools::install_github("spatialstatisticsupna/bigDM")
library(bigDM)
library(INLA)
library(parallel)
library(sf)

setwd(dirname(rstudioapi::getSourceEditorContext()$path))


#################################################
## 1) Generate the 8 validation configurations ##
#################################################
data(Data_LungCancer, package="bigDM")

S <- length(unique(Data_LungCancer$ID))
T <- length(unique(Data_LungCancer$year))

t.min <- min(Data_LungCancer$year)
t.max <- max(Data_LungCancer$year)
t.length <- 15
t.periods <- 3

Data_pred <- lapply(1:8, function(x){
  cc <- seq(t.min+(x-1), t.min+t.length+(x+1))
  data <- Data_LungCancer[Data_LungCancer$year %in% cc,]
  data$obs.true <- data$obs
  data$obs[data$year %in% tail(unique(data$year), n=t.periods)] <- NA
  
  return(data)
})

names(Data_pred) <- paste("config",1:8,sep=".")
str(Data_pred,2)


#################################
## 2) Fit the models with INLA ##
#################################
data(Carto_SpainMUN)
head(Carto_SpainMUN)

Carto_SpainMUN$ID.prov <- substr(Carto_SpainMUN$ID,1,2)


## Create a list to save the inla objects ##
aux <- vector("list",8)
names(aux) <- paste("config",1:8,sep=".")

models.INLA <- list("TypeI"=aux, "TypeII"=aux, "TypeIII"=aux, "TypeIV"=aux)
str(models.INLA)


## Select the model to be fitted ##
model <- "Disjoint"      # "Classical" or "1st-order nb"
interaction <- "TypeIV"  # "TypeI", "TypeI" or "TypeIII"


## Additional INLA parameters (faster computations) ##
## To use the same strategy of the paper set: inla.mode <- "classic"
inla.mode <- "compact"

## CAUTION: These computations are very time consuming! ##
for(i in seq(1,8)){
  if(model=="Classical"){
    res <- bigDM::STCAR_INLA(carto=Carto_SpainMUN, data=Data_pred[[i]],
                             ID.area="ID", ID.year="year", O="obs", E="pop", 
                             spatial="intrinsic", temporal="rw1", interaction=interaction, 
                             model="global", compute.fitted.values=TRUE,
                             inla.mode=inla.mode)
  }
  if(model=="Disjoint"){
    res <- bigDM::STCAR_INLA(carto=Carto_SpainMUN, data=Data_pred[[i]], ID.group="ID.prov",
                             ID.area="ID", ID.year="year", O="obs", E="pop", 
                             spatial="intrinsic", temporal="rw1", interaction=interaction, 
                             model="partition", k=0, compute.fitted.values=TRUE,
                             inla.mode=inla.mode)
  }
  if(model=="1st-order nb"){
    res <- bigDM::STCAR_INLA(carto=Carto_SpainMUN, data=Data_pred[[i]], ID.group="ID.prov",
                             ID.area="ID", ID.year="year", O="obs", E="pop", 
                             spatial="intrinsic", temporal="rw1", interaction=interaction, 
                             model="partition", k=1, compute.fitted.values=TRUE,
                             inla.mode=inla.mode)
  }
  
  res$.args$data$obs_true <- Data_pred[[i]]$obs.true
  
  models.INLA[[interaction]][[i]] <- list(data=res$.args$data,
                                          summary.fitted.values=res$summary.fitted.values,
                                          marginals.fitted.values=res$marginals.fitted.values,
                                          cpu.used=res$cpu.used)
}


##########################################
## 3) Compute model assessment criteria ##
##########################################
source("Auxiliary_functions.R")

## CAUTION: These computations are very time consuming in Windows OS!
count.pred <- lapply(models.INLA[[interaction]], function(x) compute.pred(x, ID.area=ID.area))

## Compute which areas have proportion of zero observed cases during the study period less or equal than 0.2##
cases <- matrix(Data_LungCancer$obs, nrow=S, ncol=T, byrow=F)
zero.prop <- apply(cases, 1, function(x) sum(x==0)/T)
loc.prop <- which(zero.prop <= 0.2)
  
## Compute which areas lies at the boundary between two or more provinces with at least 2 observed cases per 100000 inhabitants ##
rates <- matrix(1e+5*Data_LungCancer$obs/Data_LungCancer$pop, nrow=S, ncol=T, byrow=F)
rates.loc <- which(apply(rates,1,function(x) all(x>2)))

carto.k1 <- divide_carto(Carto_SpainMUN, ID.group="ID.prov", k=1)
ID <- unlist(lapply(carto.k1, function(x) x$ID))
ID.border <- names(which(table(ID)>1))
border.loc <- intersect(rates.loc, which(Carto_SpainMUN$ID %in% ID.border))

## Compute average values of Table 2 / Table A1 / Table A2 ##
table.loc <- list("Table 2"=1:S, "Table A1"=loc.prop, "Table A2"=border.loc)

Table <- lapply(table.loc, function(i){
  
  ## Interval score ##
  aux <- lapply(count.pred, function(x){
    matrix((x$quant0.975-x$quant0.025) + 2/0.05*(x$quant0.025-x$obs_true)*ifelse(x$obs_true<x$quant0.025,1,0) + 2/0.05*(x$obs_true-x$quant0.975)*ifelse(x$obs_true>x$quant0.975,1,0),ncol=3,byrow=F)[i,]
  })
  IS <- apply(Reduce("+",aux)/8, 2, mean)
  names(IS) <- c("1-year ahead","2-year ahead","3-year ahead")
  
  ## Mean absolute error ##
  aux <- lapply(count.pred, function(x){
    matrix(abs(x$obs_true-x$obs_pred),ncol=3,byrow=F)[i,]
  })
  MAE <- apply(Reduce("+",aux)/8, 2, mean)
  names(MAE) <- c("1-year ahead","2-year ahead","3-year ahead")
  
  ## Root mean square error ##
  aux <- lapply(count.pred, function(x){
    matrix((x$obs_true-x$obs_pred)^2,ncol=3,byrow=F)[i,]
  })
  RMSE <- apply(sqrt(Reduce("+",aux)/8), 2, mean)
  names(RMSE) <- c("1-year ahead","2-year ahead","3-year ahead")

  
  ## Computational time (minutes) ##
  Time <- mean(unlist(lapply(models.INLA[[interaction]], function(x) x$cpu.used[2]/60)))
  
  return(list(IS=IS, MAE=MAE, RMSE=RMSE, Time=Time))
})

lapply(Table, function(x) round(do.call(rbind,x[1:3]),2))


#########################
## 4) Compute Figure 2 ##
#########################