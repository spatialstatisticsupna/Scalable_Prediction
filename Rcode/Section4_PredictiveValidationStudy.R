####################################################################################################
## R script to reproduce the predictive validation study described in Orozco-Acosta et al. (2023) ##
####################################################################################################

## We recommend to install the latest version of 'bigDM' package
# devtools::install_github("spatialstatisticsupna/bigDM")
library(bigDM)
library(ggplot2)
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


#################################################################
## 2) Fit the models with INLA using the STCAR_INLA() function ##
#################################################################
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
                             spatial="BYM2", temporal="rw1", interaction=interaction, 
                             model="global", compute.fitted.values=TRUE,
                             inla.mode=inla.mode)
  }
  if(model=="Disjoint"){
    res <- bigDM::STCAR_INLA(carto=Carto_SpainMUN, data=Data_pred[[i]], ID.group="ID.prov",
                             ID.area="ID", ID.year="year", O="obs", E="pop", 
                             spatial="BYM2", temporal="rw1", interaction=interaction, 
                             model="partition", k=0, compute.fitted.values=TRUE,
                             inla.mode=inla.mode)
  }
  if(model=="1st-order nb"){
    res <- bigDM::STCAR_INLA(carto=Carto_SpainMUN, data=Data_pred[[i]], ID.group="ID.prov",
                             ID.area="ID", ID.year="year", O="obs", E="pop", 
                             spatial="BYM2", temporal="rw1", interaction=interaction, 
                             model="partition", k=1, compute.fitted.values=TRUE,
                             inla.mode=inla.mode)
  }
  
  res$.args$data$obs_true <- Data_pred[[i]]$obs.true
  
  models.INLA[[interaction]][[i]] <- list(data=res$.args$data,
                                          summary.fitted.values=res$summary.fitted.values,
                                          marginals.fitted.values=res$marginals.fitted.values,
                                          cpu.used=res$cpu.used)
}

# save(models.INLA, file=paste0("ValidationStudy_",model,".Rdata"))


##########################################
## 3) Compute model assessment criteria ##
##########################################
# load(paste0("ValidationStudy_",model,".Rdata"))
source("../Auxiliary_functions.R")

## CAUTION: These computations are very time consuming in Windows OS!
count.pred <- lapply(models.INLA[[interaction]], function(x) compute.pred(x))

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


## Create and save tables ##
aux <- lapply(Table, function(x) round(do.call(cbind,x),2))
table.columns <- list("1-year ahead"="1-year ahead","2-year ahead"="2-year ahead","3-year ahead"="3-year ahead")

Table2 <- lapply(table.columns, function(x) aux$`Table 2`[x,])
print(Table2)

TableA1 <- lapply(table.columns, function(x) aux$`Table A1`[x,])
print(TableA1)

TableA2 <- lapply(table.columns, function(x) aux$`Table A2`[x,])
print(TableA2)

save(list=c("Table2","TableA1","TableA2"), file="../../results/ValidationStudy_Tables.Rdata")


#########################
## 4) Compute Figure 2 ##
#########################
count.pred <- lapply(models.INLA[[interaction]], function(x) compute.pred(x, ID.area=c("28079","34120","05019")))

rate.pred <- lapply(names(count.pred), function(i){
  x <- Data_pred[[i]]
  y <- count.pred[[i]]
  
  ID.x <- paste(x$ID,x$year,sep=".")
  ID.y <- paste(y$Area,y$Year,sep=".")
  
  aux <- x[ID.x %in% ID.y,]
  y$rate_true <- y$obs_true/aux$pop*1e+5
  y$rate_pred <- y$obs_pred/aux$pop*1e+5
  y$`rate_q0.025` <- y$quant0.025/aux$pop*1e+5
  y$`rate_q0.975` <- y$quant0.975/aux$pop*1e+5
  
  return(y)
})
rate.pred <- do.call(rbind,rate.pred)
rate.pred <- split(rate.pred, rate.pred$Area)


## MADRID (28079) ##
Figure2a <- plot.Figure2(rate.pred$`28079`, title=paste(model,"model -",interaction,"- Madrid"))

## PALENCIA (34120) ##
Figure2b <- plot.Figure2(rate.pred$`34120`, title=paste(model,"model -",interaction,"- Palencia"))

## AVILA (05019) ##
Figure2c <- plot.Figure2(rate.pred$`05019`, title=paste(model,"model -",interaction,"- Palencia"))

## Save the plots ##
ggarrange(Figure2a, Figure2b, Figure2c, nrow=3)
ggsave(filename="../../results/Figure2.pdf", width=7.4, height=9)