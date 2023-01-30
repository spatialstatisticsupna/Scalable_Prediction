rm(list = ls())
######################
## Packages loading ##
######################

library(bigDM)
library(INLA)

#########################
## Load simulated data ##
#########################

data <- read.table(file = "data/Data_LungCancer.txt", header = TRUE)
# data <- read.table(file = "data/Data_GlobalCancer.txt", header = TRUE)
S <- 7907
T <- 25

Carto_SpainMUN$ID.prov <- substr(Carto_SpainMUN$ID, 1, 2)
ID.list <- lapply(unique(Carto_SpainMUN$ID.prov), function(x) Carto_SpainMUN[Carto_SpainMUN$ID.prov==x,]$ID)

carto.k0 <- divide_carto(carto=Carto_SpainMUN, ID.group="ID.prov", k=0)
carto.k1 <- divide_carto(carto=Carto_SpainMUN, ID.group="ID.prov", k=1)
ID.loc <- mapply(function(x,y) which(y$ID %in% x$ID), x=carto.k0, y=carto.k1, SIMPLIFY=F)

sdunif="expression:
          logdens=-log_precision/2;
          return(logdens)"

############################
## Compute LGOCV measures ##
############################

compute <- TRUE

CV <- function(x, m=3){
  x$.args$inla.call <- NULL
  n <- nrow(x$.args$data)
  LOOCV <- inla.group.cv(result=x, num.level.sets=-1)
  LGOCV <- inla.group.cv(result=x, num.level.sets=m)
  
  return(list(LOOCV=LOOCV, LGOCV=LGOCV))
}

if(compute){
  
  load("data/fitted_MODELS.Rdata")
  global_typeI <- CV(global_typeI)
  global_typeIII <- CV(global_typeIII)
  
  ##################################################################
  ## previously it is recommended to change the name by an alias  ##
  ## to the submodels obtained in the file fit_model_projection.R ##
  ##################################################################
  
  load("temp/INLAsubmodels_k0_typeI.Rdata")
  k0_typeI <- lapply(inla.models, function(x) CV(x))
  load("temp/INLAsubmodels_k0_typeII.Rdata")
  k0_typeII <- lapply(inla.models, function(x) CV(x))
  load("temp/INLAsubmodels_k0_typeIII.Rdata")
  k0_typeIII <- lapply(inla.models, function(x) CV(x))
  load("temp/INLAsubmodels_k0_typeIV.Rdata")
  k0_typeIV <- lapply(inla.models, function(x) CV(x))
  
  load("temp/INLAsubmodels_k1_typeI.Rdata")
  k1_typeI <- lapply(inla.models, function(x) CV(x))
  load("temp/INLAsubmodels_k1_typeII.Rdata")
  k1_typeII <- lapply(inla.models, function(x) CV(x))
  load("temp/INLAsubmodels_k1_typeIII.Rdata")
  k1_typeIII <- lapply(inla.models, function(x) CV(x))
  load("temp/INLAsubmodels_k1_typeIV.Rdata")
  k1_typeIV <- lapply(inla.models, function(x) CV(x))
  
  ## Save values ##
  MODELS <- list(global_typeI=global_typeI, global_typeIII=global_typeIII,
                 k0_typeI=k0_typeI, k0_typeII=k0_typeII, k0_typeIII=k0_typeIII, k0_typeIV=k0_typeIV,
                 k1_typeI=k1_typeI, k1_typeII=k1_typeII, k1_typeIII=k1_typeIII, k1_typeIV=k1_typeIV)
  
  save(list=c("MODELS"), file="data/CV_measures.Rdata")
}

#########################
## Compute LS and MSPE ##
#########################
load("data/CV_measures.Rdata")
cases <- matrix(data$obs, nrow = S, ncol = T, byrow = FALSE)
zero.prop <- apply(cases, 1, function(x) sum(x==0)/T)

cpo.LOO <- list(global_typeI=matrix(MODELS$global_typeI$LOOCV$cv, ncol=T, byrow=F),
                global_typeIII=matrix(MODELS$global_typeIII$LOOCV$cv, ncol=T, byrow=F),
                k0_typeI=do.call(rbind,lapply(MODELS$k0_typeI, function(x) matrix(x$LOOCV$cv,ncol=T,byrow=F))),
                k0_typeII=do.call(rbind,lapply(MODELS$k0_typeII, function(x) matrix(x$LOOCV$cv,ncol=T,byrow=F))),
                k0_typeIII=do.call(rbind,lapply(MODELS$k0_typeIII, function(x) matrix(x$LOOCV$cv,ncol=T,byrow=F))),
                k0_typeIV=do.call(rbind,lapply(MODELS$k0_typeIV, function(x) matrix(x$LOOCV$cv,ncol=T,byrow=F))),
                k1_typeI=do.call(rbind,mapply(function(x,y) matrix(x$LOOCV$cv, ncol=T, byrow=F)[y,], x=MODELS$k1_typeI, y=ID.loc, SIMPLIFY=F)),
                k1_typeII=do.call(rbind,mapply(function(x,y) matrix(x$LOOCV$cv, ncol=T, byrow=F)[y,], x=MODELS$k1_typeII, y=ID.loc, SIMPLIFY=F)),
                k1_typeIII=do.call(rbind,mapply(function(x,y) matrix(x$LOOCV$cv, ncol=T, byrow=F)[y,], x=MODELS$k1_typeIII, y=ID.loc, SIMPLIFY=F)),
                k1_typeIV=do.call(rbind,mapply(function(x,y) matrix(x$LOOCV$cv, ncol=T, byrow=F)[y,], x=MODELS$k1_typeIV, y=ID.loc, SIMPLIFY=F)))

cpo.LGO <- list(global_typeI=matrix(MODELS$global_typeI$LGOCV$cv, ncol=T, byrow=F),
                global_typeIII=matrix(MODELS$global_typeIII$LGOCV$cv, ncol=T, byrow=F),
                k0_typeI=do.call(rbind,lapply(MODELS$k0_typeI, function(x) matrix(x$LGOCV$cv,ncol=T,byrow=F))),
                k0_typeII=do.call(rbind,lapply(MODELS$k0_typeII, function(x) matrix(x$LGOCV$cv,ncol=T,byrow=F))),
                k0_typeIII=do.call(rbind,lapply(MODELS$k0_typeIII, function(x) matrix(x$LGOCV$cv,ncol=T,byrow=F))),
                k0_typeIV=do.call(rbind,lapply(MODELS$k0_typeIV, function(x) matrix(x$LGOCV$cv,ncol=T,byrow=F))),
                k1_typeI=do.call(rbind,mapply(function(x,y) matrix(x$LGOCV$cv, ncol=T, byrow=F)[y,], x=MODELS$k1_typeI, y=ID.loc, SIMPLIFY=F)),
                k1_typeII=do.call(rbind,mapply(function(x,y) matrix(x$LGOCV$cv, ncol=T, byrow=F)[y,], x=MODELS$k1_typeII, y=ID.loc, SIMPLIFY=F)),
                k1_typeIII=do.call(rbind,mapply(function(x,y) matrix(x$LGOCV$cv, ncol=T, byrow=F)[y,], x=MODELS$k1_typeIII, y=ID.loc, SIMPLIFY=F)),
                k1_typeIV=do.call(rbind,mapply(function(x,y) matrix(x$LGOCV$cv, ncol=T, byrow=F)[y,], x=MODELS$k1_typeIV, y=ID.loc, SIMPLIFY=F)))

LS.LOO <- data.frame(do.call(rbind,lapply(cpo.LOO, function(x) -sum(log(x), na.rm=T))),
                     do.call(rbind,lapply(cpo.LOO, function(x) -sum(log(x)[which(zero.prop<=0.2),], na.rm=T))),
                     do.call(rbind,lapply(cpo.LOO, function(x) -sum(log(x)[which(zero.prop>0.2 & zero.prop<=0.6),], na.rm=T))),
                     do.call(rbind,lapply(cpo.LOO, function(x) -sum(log(x)[which(zero.prop>0.6),], na.rm=T))))
colnames(LS.LOO) <- c("All areas", "Po<=0.2", "0.2<Po<=0.6", "Po>0.6")
round(apply(LS.LOO, 2, function(x) x-min(x)))

LS.LGO <- data.frame(do.call(rbind,lapply(cpo.LGO, function(x) -sum(log(x), na.rm=T))),
                     do.call(rbind,lapply(cpo.LGO, function(x) -sum(log(x)[which(zero.prop<=0.2),], na.rm=T))),
                     do.call(rbind,lapply(cpo.LGO, function(x) -sum(log(x)[which(zero.prop>0.2 & zero.prop<=0.6),], na.rm=T))),
                     do.call(rbind,lapply(cpo.LGO, function(x) -sum(log(x)[which(zero.prop>0.6),], na.rm=T))))
colnames(LS.LGO) <- c("All areas", "Po<=0.2", "0.2<Po<=0.6", "Po>0.6")
round(apply(LS.LGO, 2, function(x) x-min(x)))
