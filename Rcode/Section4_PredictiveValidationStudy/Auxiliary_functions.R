#' Compute summary statistics (expected value and 95% credible interval) for the posterior predictive counts
#'
#' @description The function takes an inla object and computes the expected values and quantiles for the posterior predictive counts by sampling from the posterior marginal estimates of the mortality rates.
#'
#' @param res a list containing some elements (original data, marginals of the fitted values and computational time) of the \code{inla} model.
#' @param ns number of samples (default to 5000). See \help(inla.rmarginal) for further details.
compute.pred <- function(res, ns=5000){
  
  loc_pred <- which(is.na(res$data$O))
  aux <- res$data[loc_pred,c("Area","Year","obs_true")]
  aux$period <- rep(c("1-year ahead","2-year ahead","3-year ahead"), each=length(unique(aux$Area)))
  
  set.seed(1234)
  for(i in loc_pred){
    s <- inla.rmarginal(ns, res$marginals.fitted.values[[i]])
    my.x <- s*res$data$E[i]
    samples <- rpois(ns, my.x)
    quant <- quantile(samples, probs=c(0.5, 0.025, 0.975))
    
    aux[as.character(i),"obs_pred"] <- quant[1]
    aux[as.character(i),"quant0.025"] <- quant[2]
    aux[as.character(i),"quant0.975"] <- quant[3]
  }
  
  return(aux)
}