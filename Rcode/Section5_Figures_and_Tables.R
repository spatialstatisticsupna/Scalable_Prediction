###################################################################################
## R script to reproduce the illustration section of Orozco-Acosta et al. (2023) ##
###################################################################################
library(bigDM)
library(Hmisc)
library(INLA)
library(RColorBrewer)
library(tmap)
options(timeout=600)

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source("Auxiliary_functions.R")


######################################################
######    Section 5.1: Lung cancer mortality    ######
######################################################

## Load mortality data and its associated cartography ##
load("../data/Data_LungCancer.Rdata")
str(Data_LungCancer)

data("Carto_SpainMUN")
head(Carto_SpainMUN)


## Load final model (1st-order neighbourhood + Type IV interaction) fitted using INLA ##
## Option 1: Load the previously fitted model using the 'Section5_FitModel.R' script ##
# load("../results/INLAmodel_LungCancer.Rdata")

## Option 2: Download the fitted model from fighshare ##
load(url("https://figshare.com/ndownloader/files/42527605"))

model <- k1_TypeIV
summary(model)


## Figure 3a: Posterior median estimates of LUNG CANCER mortality rates per 100,000 males ##
##            for the 7907 municipalities of continental during the period 1991-2015.     ##
##            The years 2013 to 2015 were predicted.                                      ##
S <- length(unique(Data_LungCancer$ID))
T <- length(unique(Data_LungCancer$year))
t.from <- min(unique(Data_LungCancer$year))
t.to <- max(unique(Data_LungCancer$year))

rates <- matrix(model$summary.fitted.values$`0.5quant`, ncol=T, nrow=S, byrow=FALSE)*100000
colnames(rates) <- paste("rates",seq(t.from,t.to),sep=".")

select.years <- round(seq(1,T,length.out=6))
carto <- cbind(Carto_SpainMUN, rates[,select.years])
carto$prov <- substr(carto$ID,1,2)
carto.prov <- aggregate(Carto_SpainMUN[,"geometry"], list(ID.group = carto$prov), head)

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
tmap_save(Map.Rates, file="../results/Illustration_Figure3a.jpg")
# tmap_save(Map.Rates, file="../results/Illustration_Figure3a.pdf")


## Figure 4a: Posterior predictive median estimates of LUNG CANCER mortality rates and its corresponding 95% CIs        ##
##            per 100,000 males during the period 1991-2015 for the municipalities of Gerona, Madrid and Bilbao.        ##
##            Crude rates are shown as a filled circle. The vertical dotted line indicates where the prediction starts. ##
obs_real <- matrix(Data_LungCancer$obs, nrow=S, ncol=T, byrow=F)
time <- t.from:t.to
timerev <- t.to:t.from

small.area <- c("17079","28079","48020")
name.area <- c("Gerona","Madrid","Bilbao")

graphics.off()
pdf(file="../results/Illustration_Figure4a.pdf", width=7.5, height=2.5)

par(mfrow=c(1,3))
par(mar = c(2,2,1,0.2))
for(s in 1:3){
  loc <- which(model$.args$data$Area==small.area[s])
  obs <- obs_real[loc[1],]
  E <- model$.args$data$E[loc]
  marginals <- model$marginals.fitted.values[loc]
  
  set.seed(1234)
  ns <- 5000
  N <- 25
  
  my.samples = matrix(NA, nrow=ns, ncol=N)
  for(i in 1:N){
    x = inla.rmarginal(ns, marginals[[i]])
    my.x = x*E[i]
    my.samples[,i]= rpois(ns, my.x)
  }
  
  lower <- 1e+5/E*apply(my.samples, 2, quantile, 0.025)
  upper <- 1e+5/E*apply(my.samples, 2, quantile, 0.975)
  median <- 1e+5/E*apply(my.samples, 2, quantile, 0.5)
  
  plot(time, obs, type="n", ylab="Predicted rates", xlab="", yaxt="n", xaxt="n", ylim=c(0,150), xlim=c(min(time),max(time)), cex.main = 0.8)
  polygon(c(time, timerev), c(lower, rev(upper)), col ="gray70", border = NA)
  lines(time, median, col="black", lwd=1, lty = 1)
  points(time, 1e+5*obs/E, pch=20, cex=0.5)
  abline(v=2012, col="black", lty=3, lwd=1)
  axis(1, at = c(1991,2002,2012), cex.axis = 0.8)
  minor.tick(nx = 5, ny = 5, tick.ratio = 0.3)
  axis(2, at=seq(0,150,by=25), cex.axis = 0.8)
  title(main=name.area[s], adj=0.5, line=0.3, cex.main=1)
  legend("topleft", inset = 0.01, legend=c("Crude rates","Spatio-temporal prediction","95% credible interval"),
         col=c(NA,"black", NA), lty=c(0,1,0), cex=0.9, box.lty = 0)
  legend("topleft", inset = 0.01, legend=c("","",""), col=c("black",NA,NA), pch = 20, cex = 0.9, box.lty=0, bty = "n")
  legend("topleft", inset = 0.01, legend=c("","","","",""),fill = c(NA,NA,"gray70"), cex = 0.9, box.lty=0, bty = "n",border = NA)
}
dev.off()


## Table 4: Posterior median estimates of the predicted LUNG CANCER mortality rates per 100,000 males, ##
##          its corresponding 95% credible intervals (CI) and width of the CIS for years 2013 and 2015 ##
##          for the 47 municipalities that form the provincial capitals (sorted by increasing order)   ##
ID.capitals.provinces <- c("01059","02003","03014","04013","05019","06015",
                           "08019","09059","10037","11012","12040","13034",
                           "14021","15030","16078","17079","18087","19130",
                           "20069","21041","22125","23050","24089","25120",
                           "26089","27028","28079","29067","30030","31201",
                           "32054","33044","34120","36038","37274","39075",
                           "40194","41091","42173","43148","44216","45168",
                           "46250","47186","48020","49275","50297")

names.capital <- c("Vitoria","Albacete","Alicante","Almeria","Avila","Badajoz",
                   "Barcelona","Burgos","Caceres","Cadiz","Castellon","Ciudad Real",
                   "Cordoba","A Coruna","Cuenca","Girona","Granada","Guadalajara",
                   "San Sebastian","Huelva","Huesca","Jaen","Leon","Lerida",
                   "Logrono","Lugo","Madrid","Malaga","Murcia","Pamplona","Ourense",
                   "Oviedo","Palencia","Pontevedra","Salamanca","Santander",
                   "Segovia","Sevilla","Soria","Tarragona","Teruel","Toledo",
                   "Valencia","Valladolid","Bilbao","Zamora","Zaragoza")

res <- list(data=model$.args$data,
            summary.fitted.values=model$summary.fitted.values,
            marginals.fitted.values=model$marginals.fitted.values)
res$data$obs_true <- Data_LungCancer$obs

aux <- compute.pred(res, ns=5000, ID.area=ID.capitals.provinces)
aux$rate_pred <- aux$obs_pred/aux$E*1e+5
aux$`rate_q0.025` <- aux$quant0.025/aux$E*1e+5
aux$`rate_q0.975` <- aux$quant0.975/aux$E*1e+5

Table4 <- data.frame(rate2013 = round(aux$rate_pred[which(aux$Year==2013)],1),
                     IC95_2013=paste("(",round(aux$rate_q0.025[which(aux$Year==2013)],1),",",round(aux$rate_q0.975[which(aux$Year==2013)],1),")",sep = ""),
                     Width2013=round(aux$rate_q0.975[which(aux$Year==2013)],1)-round(aux$rate_q0.025[which(aux$Year==2013)],1),
                     rate2015 = round(aux$rate_pred[which(aux$Year==2015)],1),
                     IC95_2015=paste("(",round(aux$rate_q0.025[which(aux$Year==2015)],1),",",round(aux$rate_q0.975[which(aux$Year==2015)],1),")",sep = ""),
                     Width2015=round(aux$rate_q0.975[which(aux$Year==2015)],1)-round(aux$rate_q0.025[which(aux$Year==2015)],1))
rownames(Table4) <- names.capital
Table4 <- Table4[order(Table4$rate2013),]
Table4 <- Table4[,c("rate2013","IC95_2013","Width2013","rate2015","IC95_2015","Width2015")]
print(Table4)          

save(Table4, file="../results/Illustration_Table4.Rdata")


#########################################################
######    Section 5.2: Overall cancer mortality    ######
#########################################################

## Load mortality data and its associated cartography ##
load("../data/Data_OverallCancer.Rdata")
str(Data_LungCancer)

data("Carto_SpainMUN")
head(Carto_SpainMUN)


## Load final model (1st-order neighbourhood + Type IV interaction) fitted using INLA ##
## Option 1: Load the previously fitted model using the 'Section5_FitModel.R' script ##
# load("../results/INLAmodel_OverallCancer.Rdata")

## Option 2: Download the fitted model from fighshare ##
load(url("https://figshare.com/ndownloader/files/42527620"))

model <- k1_TypeIV
summary(model)


## Figure 3b: Posterior median estimates of OVERALL CANCER mortality rates per 100,000 males ##
##            for the 7907 municipalities of continental during the period 1991-2015.        ##
##            The years 2013 to 2015 were predicted.                                         ##
S <- length(unique(Data_OverallCancer$ID))
T <- length(unique(Data_OverallCancer$year))
t.from <- min(unique(Data_OverallCancer$year))
t.to <- max(unique(Data_OverallCancer$year))

rates <- matrix(model$summary.fitted.values$`0.5quant`, ncol=T, nrow=S, byrow=FALSE)*100000
colnames(rates) <- paste("rates",seq(t.from,t.to),sep=".")

select.years <- round(seq(1,T,length.out=6))
carto <- cbind(Carto_SpainMUN, rates[,select.years])
carto$prov <- substr(carto$ID,1,2)
carto.prov <- aggregate(Carto_SpainMUN[,"geometry"], list(ID.group = carto$prov), head)

n.color <- 7
paleta <- brewer.pal(n.color,"RdYlGn")[n.color:1]
values <- c(-Inf,250,300,350,400,450,500,Inf)

Map.Rates <- tm_shape(carto) +
  tm_polygons(col=paste("rates",round(seq(t.from,t.to,length.out=6)),sep="."),
              palette=paleta, border.alpha=0, title="Est. rates", legend.show=T, legend.reverse=T,
              style="fixed", breaks=values, interval.closure="left") +
  tm_shape(carto.prov) + tm_borders(col="black", lwd=1.5, alpha=0.25) +
  tm_layout(main.title="Overall cancer mortality data", main.title.position=0.2, panel.label.size=1.5,
            legend.outside=T, legend.outside.position="right", legend.frame=F,
            panel.labels=paste("Year",round(seq(t.from,t.to,length.out=6))),
            legend.outside.size=0.3, outer.margins=c(0.0,0.03,0.01,0)) +
  tm_facets(nrow=2, ncol=3)

tmap_mode("plot")
tmap_save(Map.Rates, file="../results/Illustration_Figure3b.jpg")
# tmap_save(Map.Rates, file="../results/Illustration_Figure3b.pdf")


##########################################################################################################################
## Figure 4b: Posterior predictive median estimates of OVERALL CANCER mortality rates and its corresponding 95% CIs     ##
##            per 100,000 males during the period 1991-2015 for the municipalities of Gerona, Madrid and Bilbao.        ##
##            Crude rates are shown as a filled circle. The vertical dotted line indicates where the prediction starts. ##
##########################################################################################################################
obs_real <- matrix(Data_OverallCancer$obs, nrow=S, ncol=T, byrow=F)
time <- t.from:t.to
timerev <- t.to:t.from

small.area <- c("17079","28079","48020")
name.area <- c("Gerona","Madrid","Bilbao")

graphics.off()
pdf(file="../results/Illustration_Figure4b.pdf", width=7.5, height=2.5)

par(mfrow=c(1,3))
par(mar = c(2,2,1,0.2))
for(s in 1:3){
  loc <- which(model$.args$data$Area==small.area[s])
  obs <- obs_real[loc[1],]
  E <- model$.args$data$E[loc]
  marginals <- model$marginals.fitted.values[loc]
  
  set.seed(1234)
  ns <- 5000
  N <- 25
  
  my.samples = matrix(NA, nrow=ns, ncol=N)
  for(i in 1:N){
    x = inla.rmarginal(ns, marginals[[i]])
    my.x = x*E[i]
    my.samples[,i]= rpois(ns, my.x)
  }
  lower <- 1e+5/E*apply(my.samples, 2, quantile, 0.025)
  upper <- 1e+5/E*apply(my.samples, 2, quantile, 0.975)
  median <- 1e+5/E*apply(my.samples, 2, quantile, 0.5)
  
  plot(time, obs, type="n", ylab="Predicted rates", xlab="", yaxt="n", xaxt="n", ylim=c(0,600), xlim=c(min(time),max(time)), cex.main = 0.8)
  polygon(c(time, timerev), c(lower, rev(upper)), col ="gray70", border = NA)
  lines(time, median, col="black", lwd=1, lty = 1)
  points(time, 1e+5*obs/E, pch=20, cex=0.5)
  abline(v=2012, col="black", lty=3, lwd=1)
  axis(1, at = c(1991,2002,2012), cex.axis = 0.8)
  minor.tick(nx = 5, ny = 5, tick.ratio = 0.3)
  axis(2, at=seq(0,600,by=100), cex.axis = 0.8)
  title(main=name.area[s], adj=0.5, line=0.3, cex.main=1)
  legend("topleft", inset = 0.01, legend=c("Crude rates","Spatio-temporal prediction","95% credible interval"),
         col=c(NA,"black", NA), lty=c(0,1,0), cex=0.9, box.lty = 0)
  legend("topleft", inset = 0.01, legend=c("","",""), col=c("black",NA,NA), pch = 20, cex = 0.9, box.lty=0, bty = "n")
  legend("topleft", inset = 0.01, legend=c("","","","",""),fill = c(NA,NA,"gray70"), cex = 0.9, box.lty=0, bty = "n",border = NA)
}
dev.off()


## Table 5: Posterior median estimates of the predicted OVERALL CANCER mortality rates per 100,000 males, ##
##          its corresponding 95% credible intervals (CI) and width of the CIS for years 2013 and 2015    ##
##          for the 47 municipalities that form the provincial capitals (sorted by increasing order)      ##

ID.capitals.provinces <- c("01059","02003","03014","04013","05019","06015",
                           "08019","09059","10037","11012","12040","13034",
                           "14021","15030","16078","17079","18087","19130",
                           "20069","21041","22125","23050","24089","25120",
                           "26089","27028","28079","29067","30030","31201",
                           "32054","33044","34120","36038","37274","39075",
                           "40194","41091","42173","43148","44216","45168",
                           "46250","47186","48020","49275","50297")

names.capital <- c("Vitoria","Albacete","Alicante","Almeria","Avila","Badajoz",
                   "Barcelona","Burgos","Caceres","Cadiz","Castellon","Ciudad Real",
                   "Cordoba","A Coruna","Cuenca","Girona","Granada","Guadalajara",
                   "San Sebastian","Huelva","Huesca","Jaen","Leon","Lerida",
                   "Logrono","Lugo","Madrid","Malaga","Murcia","Pamplona","Ourense",
                   "Oviedo","Palencia","Pontevedra","Salamanca","Santander",
                   "Segovia","Sevilla","Soria","Tarragona","Teruel","Toledo",
                   "Valencia","Valladolid","Bilbao","Zamora","Zaragoza")

res <- list(data=model$.args$data,
            summary.fitted.values=model$summary.fitted.values,
            marginals.fitted.values=model$marginals.fitted.values)
res$data$obs_true <- Data_OverallCancer$obs

aux <- compute.pred(res, ns=5000, ID.area=ID.capitals.provinces)
aux$rate_pred <- aux$obs_pred/aux$E*1e+5
aux$`rate_q0.025` <- aux$quant0.025/aux$E*1e+5
aux$`rate_q0.975` <- aux$quant0.975/aux$E*1e+5

Table5 <- data.frame(rate2013 = round(aux$rate_pred[which(aux$Year==2013)],1),
                     IC95_2013=paste("(",round(aux$rate_q0.025[which(aux$Year==2013)],1),",",round(aux$rate_q0.975[which(aux$Year==2013)],1),")",sep = ""),
                     Width2013=round(aux$rate_q0.975[which(aux$Year==2013)],1)-round(aux$rate_q0.025[which(aux$Year==2013)],1),
                     rate2015 = round(aux$rate_pred[which(aux$Year==2015)],1),
                     IC95_2015=paste("(",round(aux$rate_q0.025[which(aux$Year==2015)],1),",",round(aux$rate_q0.975[which(aux$Year==2015)],1),")",sep = ""),
                     Width2015=round(aux$rate_q0.975[which(aux$Year==2015)],1)-round(aux$rate_q0.025[which(aux$Year==2015)],1))
rownames(Table5) <- names.capital
Table5 <- Table5[order(Table5$rate2013),]
Table5 <- Table5[,c("rate2013","IC95_2013","Width2013","rate2015","IC95_2015","Width2015")]
print(Table5)

save(Table5, file="../results/Illustration_Table5.Rdata")
