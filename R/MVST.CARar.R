MVST.CARar <- function(formula, family, data=NULL,  trials=NULL, W, burnin, n.sample, thin=1, prior.mean.beta=NULL, prior.var.beta=NULL, prior.nu2=NULL, prior.Sigma.df=NULL, prior.Sigma.scale=NULL, AR=NULL, rho.S=NULL, rho.T=NULL, MALA=FALSE, verbose=TRUE)
{
    ## This is a wrapper function for the following six functions.
    ## binomial.MVCARar1
    ## gaussian.MVCARar1
    ## poisson.MVCARar1
    ## binomial.MVCARar2
    ## gaussian.MVCARar2
    ## poisson.MVCARar2
    if(is.null(family)) stop("the 'family' argument is missing.", call.=FALSE)
    if(is.null(AR)) stop("the 'AR' argument is missing, please specify 1 for an AR(1) model and 2 for an AR(2) model.", call.=FALSE)
    
    
    #### Run the appropriate model according to the family argument and ar argument
    if(family=="binomial" & AR==1)
    {
        if(is.null(trials)) stop("a binomial model was specified but the trials arugment was not specified.", call.=FALSE)
    model <- binomial.MVCARar1(formula=formula, data=data,  trials=trials, W=W, burnin=burnin, n.sample=n.sample, thin=thin, prior.mean.beta=prior.mean.beta, prior.var.beta=prior.var.beta, prior.Sigma.df=prior.Sigma.df, prior.Sigma.scale=prior.Sigma.scale, rho.S=rho.S, rho.T=rho.T, MALA=MALA, verbose=verbose)
    }else if(family=="binomial" & AR==2)
    {
        if(is.null(trials)) stop("a binomial model was specified but the trials arugment was not specified.", call.=FALSE)
    model <- binomial.MVCARar2(formula=formula, data=data,  trials=trials, W=W, burnin=burnin, n.sample=n.sample, thin=thin, prior.mean.beta=prior.mean.beta, prior.var.beta=prior.var.beta, prior.Sigma.df=prior.Sigma.df, prior.Sigma.scale=prior.Sigma.scale, rho.S=rho.S, rho.T=rho.T, MALA=MALA, verbose=verbose)
    }else if(family=="gaussian" & AR==1)
    {
    model <- gaussian.MVCARar1(formula=formula, data=data,  W=W, burnin=burnin, n.sample=n.sample, thin=thin, prior.mean.beta=prior.mean.beta, prior.var.beta=prior.var.beta, prior.nu2=prior.nu2, prior.Sigma.df=prior.Sigma.df, prior.Sigma.scale=prior.Sigma.scale,  rho.S=rho.S, rho.T=rho.T, verbose=verbose)          
    }else if(family=="gaussian" & AR==2)
    {
    model <- gaussian.MVCARar2(formula=formula, data=data,  W=W, burnin=burnin, n.sample=n.sample, thin=thin, prior.mean.beta=prior.mean.beta, prior.var.beta=prior.var.beta, prior.nu2=prior.nu2, prior.Sigma.df=prior.Sigma.df, prior.Sigma.scale=prior.Sigma.scale,  rho.S=rho.S, rho.T=rho.T, verbose=verbose)          
    }else if(family=="poisson" & AR==1)
    {
    model <- poisson.MVCARar1(formula=formula, data=data,  W=W, burnin=burnin, n.sample=n.sample, thin=thin, prior.mean.beta=prior.mean.beta, prior.var.beta=prior.var.beta, prior.Sigma.df=prior.Sigma.df, prior.Sigma.scale=prior.Sigma.scale,  rho.S=rho.S, rho.T=rho.T, MALA=MALA, verbose=verbose)          
    }else if(family=="poisson" & AR==2)
    {
    model <- poisson.MVCARar2(formula=formula, data=data,  W=W, burnin=burnin, n.sample=n.sample, thin=thin, prior.mean.beta=prior.mean.beta, prior.var.beta=prior.var.beta, prior.Sigma.df=prior.Sigma.df, prior.Sigma.scale=prior.Sigma.scale,  rho.S=rho.S, rho.T=rho.T, MALA=MALA, verbose=verbose)          
    }else
    {
        stop("the 'family' arugment is not one of `binomial', `gaussian' or `poisson' or the 'AR' argument is not '1' or '2'.", call.=FALSE)     
    }  
    
    
    return(model)     
}