#'@title Class for fit to evolutionary sequence (time-series) models
#'
#'@description A function that combines useful information summarizing model fit.
#'
#'@param converge info on model convergence 
#'
#'@param modelName description of the model
#'
#'@param logL log-likelihood of model
#'
#'@param ancestral.values maximum-likelihood estimates of the ancestral trait values
#'
#'@param SE.anc standard errors of the estimated ancestral states
#'
#'@param R maximum-likelihood estimates of the parameters in the R matrix
#'
#'@param SE.R standard errors of the parameters in the R matrix
#'
#'@param shift.point the sample in the time series that represents the first sample in the second segment.
#'
#'@param method the parameterization used: Joint
#'
#'@param K number of parameters in the model
#'
#'@param n sample size
#'
#'@param iter the number of times the optimization method is run from different starting points. Default is NULL, meaning the optimization is run once.
#'
#'@details This function is used by the model-fitting routines for the multivariate Unbiased Random Walk models to create standardized output
#'
#'@note This function is not likely to be called directly by the user.
#'
#'@author Kjetil Lysne Voje


as.evoTS.multi.URW.shift.fit<-function (converge, modelName, logL, ancestral.values, SE.anc, R, SE.R, shift.point, method, K, n, iter)
{
  ic <- paleoTS::IC(logL = logL, K = K, n = n, method = "AICc")
  y <- list(converge = converge, modelName = modelName, logL = logL, AICc = ic, ancestral.values = ancestral.values, SE.anc = SE.anc, R = R, SE.R = SE.R,
            shift.point = shift.point, method = method, K = K, n = n, iter = iter)
  class(y) <- "evoTSmvFit"
  return(y)
}
