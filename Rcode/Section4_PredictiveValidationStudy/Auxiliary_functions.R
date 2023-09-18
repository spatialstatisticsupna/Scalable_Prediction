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
  aux <- res$data[loc_pred,c("Area","Year","obs_true")]
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