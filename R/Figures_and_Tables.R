rm(list=ls())
library(bigDM)
library(Hmisc)
library(INLA)
library(RColorBrewer)
library(tmap)


############################################################################
## Load the cartography with the 7907 municipalities of continental Spain ##
############################################################################
data("Carto_SpainMUN")
head(Carto_SpainMUN)

Carto_SpainMUN$prov <- substr(Carto_SpainMUN$ID,1,2)
carto.prov <- aggregate(Carto_SpainMUN[,"geometry"], list(ID.group=Carto_SpainMUN$prov), head)


#####################################################################################################
## Load the final model (1st-order neighbourhood + Type IV interaction) fitted using bigDM,        ##
## which can be downloaded from: https://emi-sstcdapp.unavarra.es/bigDM/k1_typeIV_simulated.Rdata  ##
#####################################################################################################
load("k1_typeIV_simulated.Rdata")

Model <- k1_typeIV
summary(Model)


#################################################################################################
## Table 3: Posterior median estimates of the predicted lung cancer mortality rates per        ##
##          100 000 males and its corresponding 95% credible intervals for years 2013 and 2015 ##
##          for the 47 municipalities that form the provincial capitals.                       ##
#################################################################################################
loc.pred <- is.na(Model$.args$data$O)
E <- Model$.args$data$E[loc.pred]
marginals <- Model$marginals.fitted.values[loc.pred]

N <- length(marginals)
ns <-  5000
my.samples <- matrix(NA, nrow=ns, ncol=N)

set.seed(1234)
for(i in 1:N){
  x <- inla.rmarginal(ns, marginals[[i]])
  my.x <- x*E[i]
  my.samples[,i] <- rpois(ns, my.x)
}

quant0025 <- matrix(apply(my.samples, 2, quantile, 0.025), ncol=3, nrow=7907, byrow=FALSE)
quant0025 <- 10^5*quant0025/matrix(E,ncol=3,nrow=7907,byrow=FALSE)
quant0500 <- matrix(apply(my.samples, 2, quantile, 0.50), ncol=3, nrow=7907, byrow=FALSE)
quant0500 <- 10^5*quant0500/matrix(E,ncol=3,nrow=7907,byrow=FALSE)
quant0975 <- matrix(apply(my.samples, 2, quantile, 0.975), ncol=3, nrow=7907, byrow=FALSE)
quant0975 <- 10^5*quant0975/matrix(E,ncol=3,nrow=7907,byrow=FALSE)

## Provincial capitals ##
ID.capitals.provinces <- c("01059","02003","03014","04013","05019","06015",
                           "08019","09059","10037","11012","12040","13034",
                           "14021","15030","16078","17079","18087","19130",
                           "20069","21041","22125","23050","24089","25120",
                           "26089","27028","28079","29067","30030","31201",
                           "32054","33044","34120","36038","37274","39075",
                           "40194","41091","42173","43148","44216","45168",
                           "46250","47186","48020","49275","50297")

loc.capital <- which(Carto_SpainMUN$ID %in% ID.capitals.provinces)
names.capital <- Carto_SpainMUN$name[Carto_SpainMUN$ID %in% ID.capitals.provinces]

quant0025 <- quant0025[loc.capital,]
quant0500 <- quant0500[loc.capital,]
quant0975 <- quant0975[loc.capital,]

IC95_2013 <- paste("(",round(quant0025[,1],1),",",round(quant0975[,1],1),")",sep = "")
IC95_2015 <- paste("(",round(quant0025[,3],1),",",round(quant0975[,3],1),")",sep = "")
tab <- data.frame(Obs2013 = quant0500[,1], 
                  IC95_2013=IC95_2013, 
                  Obs2015 = quant0500[,3], 
                  IC95_2015=IC95_2015)
rownames(tab) <- names.capital
tab <- tab[order(tab$Obs2013),]
print(tab)


###################################################################################
## Figure 3: Posterior median estimates of estimated lung cancer mortality rates ##
##           per 100000 males for the 7907 municipalities of continental during  ##
##           the period 1991-2015. The years 2013 to 2015 were predicted.        ##
####################################################################################
rates <- matrix(Model$summary.fitted.values$`0.5quant`, ncol=T, nrow=S, byrow = FALSE)*100000
colnames(rates) <- paste("rates",seq(t.from,t.to),sep=".")

select.years <- round(seq(1,T,length.out=6))
carto <- cbind(Carto_SpainMUN, rates[,select.years])

n.color <- 7
paleta <- brewer.pal(n.color,"RdYlGn")[n.color:1]
values <- c(-Inf,60,70,80,85,95,105,Inf)

Map.Rates <- tm_shape(carto) +
  tm_polygons(col=paste("rates",round(seq(t.from,t.to,length.out=6)),sep="."),
              palette=paleta, border.alpha=0, title="Est. rates", legend.show=T, legend.reverse=T,
              style="fixed", breaks=values, interval.closure="left") +
  tm_shape(carto.prov) + tm_borders(col="black", lwd=1.5, alpha=0.25) +
  tm_layout(main.title="Lung cancer mortality data", main.title.position=0.2, panel.label.size=1.5,
            legend.outside=T, legend.outside.position="right", legend.frame=F,
            panel.labels=paste("Year",round(seq(t.from,t.to,length.out=6))),
            legend.outside.size=0.3, outer.margins=c(0.0,0.03,0.01,0)) + 
  tm_facets(nrow=2, ncol=3)

tmap_mode("plot")
tmap_save(Map.Rates, file="Figure3.pdf")


##################################################################################################
## Figure 4: Posterior predictive median and its corresponding 95% CI of lung cancer (top) and  ##
##           overall cancer (bottom) deaths during the period 1991-2015 for the municipalities  ##
##           of Guadalajara, Madrid and Bilbao. Crude rates are shown as a filled circle.       ##
##           The vertical dotted line indicates where the prediction starts.                    ##
##################################################################################################
obs.real <- matrix(data$O, ncol=T, nrow=S)

small.area <- c("17079","28079","48020")
name.area <- c("Gerona","Madrid","Bilbao")

graphics.off()
pdf(file="Figure4.pdf", width=7.5, height=2.5)

par(mfrow=c(1,3))
par(mar=c(2,2,1,0.2))

for(s in 1:3){
  loc <- which(Model$.args$data$Area==small.area[s])
  obs <- obs.real[loc[1],]
  E <- Model$.args$data$E[loc]
  marginals <- Model$marginals.fitted.values[loc]
  
  set.seed(1234)
  ns <- 5000
  my.samples <- matrix(NA, nrow=ns, ncol=T)
  for(i in 1:T){
    x = inla.rmarginal(ns, marginals[[i]])
    my.x = x*E[i]
    my.samples[,i]= rpois(ns, my.x)
  }
  
  lower <- 1e+5/E*apply(my.samples, 2, quantile, 0.025)
  upper <- 1e+5/E*apply(my.samples, 2, quantile, 0.975)
  median <- 1e+5/E*apply(my.samples, 2, quantile, 0.5)
  
  time <- t.from:t.to
  plot(time, obs, type="n", ylab="Predicted rates", xlab="", yaxt="n", xaxt="n", ylim=c(0,150), xlim=range(time), cex.main = 0.8)
  polygon(c(time, rev(time)), c(lower, rev(upper)), col ="gray70", border = NA)
  lines(time, median, col="black", lwd=1, lty=1)
  points(time, 1e+5*obs/E, pch=20, cex=0.5)     
  abline(v=2012, col="black", lty=3, lwd=1)
  axis(1, at = c(1991,2002,2012), cex.axis=0.8)
  minor.tick(nx=5, ny=5, tick.ratio=0.3)
  axis(2, at=seq(0,150,by=25), cex.axis=0.8)
  title(main=name.area[s], adj=0.5, line=0.3, cex.main=1)
  legend("topleft", inset=0.01, legend=c("Crude rates","Spatio-temporal prediction","95% credible interval"), 
         col=c(NA,"black",NA), lty=c(0,1,0), cex=0.9, box.lty=0)
  legend("topleft", inset=0.01, legend=c("","",""), col=c("black",NA,NA), pch=20, cex=0.9, box.lty=0, bty="n")
  legend("topleft", inset=0.01, legend=c("","","","",""), fill=c(NA,NA,"gray70"), cex=0.9, box.lty=0, bty="n",border = NA)
}
dev.off()


##############################################################################################
## Figure 5: Posterior median of the percentage change of mortality rates from 2013 to 2015 ##
##############################################################################################

## PENDIENTE ##