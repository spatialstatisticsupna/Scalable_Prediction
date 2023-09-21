#' Compute summary statistics (expected value and 95% credible interval) for the posterior predictive counts
#'
#' @description The function takes an inla object and computes the expected values and quantiles for the posterior predictive counts by sampling from the posterior marginal estimates of the mortality rates.
#'
#' @param res a list containing some elements (original data, marginals of the fitted values and computational time) of the \code{inla} model.
#' @param ns number of samples (default to 5000). See \help(inla.rmarginal) for further details.
#' @param ID.area optional argument (default \code{NULL}) to specify for which areas the posterior predictive counts should be computed.
#' 
#' @return This function returns a \code{data.frame} object with the following variables:
#' \itemize{
#'   \item \code{Area}: character vector of geographic identifiers
#'   \item \code{Year}: numeric vector of years identifiers
#'   \item \code{obs_true}: numeric vector with the real observed number of cases
#'   \item \code{period}: character vector indicating one, two or three-year ahead predictions
#'   \item \code{obs_pred}: numeric vector with the expected values of the posterior predictive counts
#'   \item \code{quant0.025}: numeric vector with the 2.5 quantile of the posterior predictive counts
#'   \item \code{quant0.975}: numeric vector with the 97.5 quantile of the posterior predictive counts
#' }
#'
compute.pred <- function(res, ns=5000, ID.area=NULL){
  
  if(is.null(ID.area)) ID.area <- unique(res$data$Area)
  
  loc_pred <- which(is.na(res$data$O) & (res$data$Area %in% ID.area))
  aux <- res$data[loc_pred,c("Area","Year","E","obs_true")]
  aux$period <- rep(c("1-year ahead","2-year ahead","3-year ahead"), each=length(unique(aux$Area)))
  
  if(Sys.info()[1]=="Windows"){
    set.seed(1234)
    
    pred <- do.call(rbind,lapply(loc_pred, function(i) {
      s <- inla.rmarginal(ns, res$marginals.fitted.values[[i]])
      my.x <- s*res$data$E[i]
      samples <- rpois(ns, my.x)
      quant <- quantile(samples, probs=c(0.025, 0.975))
      
      data.frame(obs_pred=res$summary.fitted.values$mean[[i]]*res$data$E[i],
                 `quant0.025`=quant[1], `quant0.975`=quant[2])
    }))
  }else{
    pred <- do.call(rbind,mclapply(loc_pred, function(i) {
      set.seed(1234)
      
      s <- inla.rmarginal(ns, res$marginals.fitted.values[[i]])
      my.x <- s*res$data$E[i]
      samples <- rpois(ns, my.x)
      quant <- quantile(samples, probs=c(0.025, 0.975))
      
      data.frame(obs_pred=res$summary.fitted.values$mean[[i]]*res$data$E[i],
                 `quant0.025`=quant[1], `quant0.975`=quant[2])
    },mc.cores=detectCores()))
  }
  
  aux <- cbind(aux,pred)
  
  return(aux)
}


#' Auxiliary function to plot Figure 2
#'
#' @description The function takes a data frame with the expected values and quantiles for the posterior predictive rates (cases per 100,000 inhabitants) obtained from the differente configurations of the validation study.
#'
#' @param aux \code{data.frame} object with the expected values and quantiles for the posterior predictive rates (cases per 100,000 inhabitants)
#' @param title character string with the title of the plot
#' 
#' @return A \code{ggplot} object
#'
plot.Figure2 <- function(aux, title=NULL){
  
  xx <- ggplot(data=aux,
               mapping=aes(x=Year, y=rate_pred, ymin=rate_q0.025, ymax=rate_q0.975,
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
    ggtitle(title) + 
    geom_point(shape="asterisk", aes(x=Year, y=rate_true), colour="black", size = 0.5)
}


#' From a fitted model, compute \code{group.cv} values
#'
#' @description The function takes an inla object and computes measures for leave-one-out cross-validation (LOOCV) and leave-group-out cross-validation (LGOCV).
#'
#' @param x an object of class \code{inla}.
#' @param m number of level sets to use. See \help(inla.group.cv) for further details.
#' 
#' @return This function returns a \code{list} object with the following elements:
#' \itemize{
#'   \item \code{LOOCV}: a list related to leave-one-out cross-validation measures
#'   \item \code{LGOCV}: a list related to leave-group-out cross-validation measures
#'   \item \code{data}: data frame extracted from the inla model
#' }
#'
CV <- function(x, m=3){
  
  sdunif="expression:
          logdens=-log_precision/2;
          return(logdens)"
  
  x$.args$inla.call <- NULL
  
  LOOCV <- inla.group.cv(result=x, num.level.sets=-1)
  LGOCV <- inla.group.cv(result=x, num.level.sets=m)
  
  data <- x$.args$data
  data$ID <- paste(data$Year, data$Area, sep=".")
  
  return(list(LOOCV=LOOCV, LGOCV=LGOCV, data=data))
}