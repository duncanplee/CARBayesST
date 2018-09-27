ST.CARclustrends <- function(formula, family, data=NULL, trials=NULL, W, burnin, n.sample, thin=1, trends=NULL, changepoint=NULL, knots=NULL, prior.mean.beta=NULL, prior.var.beta=NULL, prior.mean.gamma=NULL, prior.var.gamma=NULL, prior.lambda=NULL, prior.tau2=NULL, Nchains=4, verbose=TRUE)
{
  ## This is a wrapper function for the following two functions.
  ## binomial.CARclustrends
  ## poisson.CARclustrends
  if(is.null(family)) stop("the family argument is missing", call.=FALSE)
  
  #### Run the appropriate model according to the family argument
  if(family=="binomial")
  {
    if(is.null(trials)) stop("a binomial model was specified but the trials arugment was not specified", call.=FALSE)
    model <- binomial.CARclustrends(formula=formula, data=data, trials=trials, W=W, burnin=burnin, n.sample=n.sample, thin=thin, trends=trends, changepoint=changepoint, knots=knots,
                                    prior.mean.beta=prior.mean.beta, prior.var.beta=prior.var.beta, prior.mean.gamma=prior.mean.gamma, prior.var.gamma=prior.var.gamma,
                                    prior.lambda=prior.lambda, prior.tau2=prior.tau2, Nchains=Nchains, verbose=TRUE)
  }else if(family=="poisson")
  {
    if(!is.null(trials)) stop("you do not need a trials arugment as a binomial model was not specified", call.=FALSE)
    model <- poisson.CARclustrends(formula=formula, data=data, W=W, burnin=burnin, n.sample=n.sample, thin=thin, trends=trends, changepoint=changepoint, knots=knots,
                                   prior.mean.beta=prior.mean.beta, prior.var.beta=prior.var.beta, prior.mean.gamma=prior.mean.gamma, prior.var.gamma=prior.var.gamma,
                                   prior.lambda=prior.lambda, prior.tau2=prior.tau2, Nchains=Nchains, verbose=TRUE)    
  }else
  {
    stop("the family argument is not one of `binomial' or `poisson'.", call.=FALSE)
  }
  return(model)     
}
