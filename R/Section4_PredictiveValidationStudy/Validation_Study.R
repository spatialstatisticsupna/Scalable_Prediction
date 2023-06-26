rm(list=ls())
library(bigDM)
library(ggplot2)
library(ggpubr)
library(INLA)



#########################
## Auxiliary functions ##
#########################
compute.pred <- function(x, ns=5000){
  
  loc_pred <- which(is.na(x$summary.fitted$O))
  aux <- x$summary.fitted[loc_pred,c("Area","Year")]
  aux$period <- rep(c("1-year ahead","2-year ahead","3-year ahead"), length(ID.area))
  
  set.seed(1234)
  for(i in loc_pred){
    s <- inla.rmarginal(ns, x$marginals.fitted[[i]])
    my.x <- s*x$summary.fitted$E[i]
    samples <- rpois(ns, my.x)
    quant <- quantile(samples/x$summary.fitted$E[i]*1e+5, probs=c(0.5, 0.025, 0.975))
    
    aux[as.character(i),"rate_pred"] <- quant[1]
    aux[as.character(i),"quant0.025"] <- quant[2]
    aux[as.character(i),"quant0.975"] <- quant[3]
  }
  
  return(aux)
}

plot.Figure2 <- function(aux, title=NULL){
  
  xx <- ggplot(data=aux,
               mapping=aes(x=Year, y=rate_pred, ymin=quant0.025, ymax=quant0.975,
                           group=period, color=period))
  
  resul <- xx + geom_pointrange(position=position_dodge(width=0.7), shape=21, fatten=0.5, size=0.2) +
    scale_x_discrete(breaks=as.numeric(unique(aux$Year))) +
    theme(legend.position="bottom", plot.title=element_text(size=8, hjust=0.5),
          axis.title=element_text(size=7),
          axis.text=element_text(size=7),
          axis.text.x=element_text(angle=90),
          legend.title=element_text(size=7),
          legend.text=element_text(size=7),
          legend.box.spacing=unit(0,"pt")) +
    labs(x="Year", y="Predicted counts (per 100,000 inhabitants)", color="") +
    ylim(0,200) +
    ggtitle(title)
  
  return(resul)
}


#########################################################################
## Generate data frames to replicate the validation study of Section 4 ##
#########################################################################
data("Data_LungCancer")
head(Data_LungCancer)

t.min <- min(Data_LungCancer$year)
t.max <- max(Data_LungCancer$year)

t.length <- 15
t.periods <- 3

Data_pred <- lapply(1:8, function(x){
  cc <- seq(t.min+(x-1), t.min+t.length+(x+1))
  data <- Data_LungCancer[Data_LungCancer$year %in% cc,]
  data$obs[data$year %in% tail(unique(data$year), n=t.periods)] <- NA

  return(data)
})

names(Data_pred) <- paste("config",1:8,sep=".")
str(Data_pred,2)


##################################################################################################
## Figure 2: One, two and three-year ahead predictions for the municipalities of Madrid (top),  ##
##           Palencia (middle) and Ávila (bottom) using the disjoint model (left column) and    ##
##           1st-order neighbourhood model (right column) with Type IV interactions.            ##
##################################################################################################
load("ValidationStudy_Figure2.Rdata")

ID.area <- unique(ValidationStudy_k0_typeIV[[1]]$summary.fitted$Area)

data.pred.k0 <- do.call(rbind,lapply(ValidationStudy_k0_typeIV, compute.pred))
data.pred.k1 <- do.call(rbind,lapply(ValidationStudy_k1_typeIV, compute.pred))


####################
## MADRID (28079) ##
####################
Figure2a <- plot.Figure2(data.pred.k0[data.pred.k0$Area=="28079",], title="Disjoint model - Type IV - Madrid")
Figure2b <- plot.Figure2(data.pred.k1[data.pred.k1$Area=="28079",], title="1st order neighbourhood model - Type IV - Madrid")
Madrid <- ggarrange(Figure2a, Figure2b, ncol=2)

######################
## PALENCIA (34120) ##
######################
Figure2a <- plot.Figure2(data.pred.k0[data.pred.k0$Area=="34120",], title="Disjoint model - Type IV - Palencia")
Figure2b <- plot.Figure2(data.pred.k1[data.pred.k1$Area=="34120",], title="1st order neighbourhood model - Type IV - Palencia")
Palencia <- ggarrange(Figure2a, Figure2b, ncol=2)

###################
## AVILA (05019) ##
###################
Figure2a <- plot.Figure2(data.pred.k0[data.pred.k0$Area=="05019",], title="Disjoint model - Type IV - Ávila")
Figure2b <- plot.Figure2(data.pred.k1[data.pred.k1$Area=="05019",], title="1st order neighbourhood model - Type IV - Avila")
Avila <- ggarrange(Figure2a, Figure2b, ncol=2)

ggarrange(Madrid, Palencia, Avila, nrow=3)
ggsave(filename="Figure2.pdf", width=7.4, height=9)
