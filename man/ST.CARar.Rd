\name{ST.CARar}
\alias{ST.CARar}


\title{
    Fit a spatio-temporal generalised linear mixed model to data, with a spatio-temporal
    autoregressive process.
}

\description{
    Fit a spatio-temporal generalised linear mixed model to areal unit data, where the response
    variable can be binomial, Gaussian or Poisson. The linear predictor is modelled by known 
    covariates and a vector of random effects. The latter follows a multivariate first order 
    autoregressive time series process, where spatial autocorrelation is modelled via the
    precision matrix. This precision matrix comes from the conditional autoregressive (CAR) 
    prior proposed by Leroux et al. (2000), and the full model was proposed by Rushworth et al.
    (2014).  Missing values are allowed in the response in this model. Further details are 
    given in the vignette accompanying this package. Inference is conducted in a Bayesian 
    setting using Markov chain Monte Carlo (MCMC) simulation. 
}



\usage{
    ST.CARar(formula, family, data=NULL,  trials=NULL, W, burnin, n.sample, thin=1, 
    prior.mean.beta=NULL, prior.var.beta=NULL, prior.nu2=NULL, prior.tau2=NULL,
    fix.rho.S=FALSE, rho.S=NULL, fix.rho.T=FALSE, rho.T=NULL, MALA=TRUE, verbose=TRUE)
}


%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{formula}{
        A formula for the covariate part of the model using the syntax of the
        lm() function. Offsets can be included here using the offset() function. 
        The response and each covariate should be vectors of length (KN)*1, where
        K is the number of spatial units and N is the number of time periods. Each
        vector should be ordered so that the first K data points are the set of all
        K spatial locations at time 1, the next K are the set of spatial locations for
        time 2 and so on. Missing values are allowed in the response.
    }
    \item{family}{
        One of either `binomial', `gaussian' or `poisson', which respectively specify a 
        binomial likelihood model with a logistic link function, a Gaussian likelihood 
        model with an identity link function, or a Poisson likelihood model with a 
        log link function. 
    }
    \item{data}{
        An optional data.frame containing the  variables in the formula.
    }
    \item{trials}{
        A vector the same length and in the same order as the response containing the 
        total number of trials for each area and time period. Only used if family=`binomial'. 
    }
        \item{W}{
        A non-negative K by K neighbourhood matrix (where K is the number of spatial units). 
        Typically a binary specification is used, where the jkth element equals one
        if areas (j, k) are spatially close (e.g. share a common border) and is zero 
        otherwise. The matrix can be non-binary, but each row must contain at least 
        one non-zero entry. Note, if fix.rho.S=TRUE and rho.S=0 then the W matrix is 
        not actually used by the model so a K by K matrix of zeros could be entered
        instead if a real W matrix cannot be obtained.
        }  
        \item{burnin}{
        The number of MCMC samples to discard as the burn-in period.
        }
        \item{n.sample}{
        The number of MCMC samples to generate.
        }
        \item{thin}{
        The level of thinning to apply to the MCMC samples to reduce their temporal 
        autocorrelation. Defaults to 1 (no thinning).
        }
        \item{prior.mean.beta}{
        A vector of prior means for the regression parameters beta (Gaussian priors are 
        assumed). Defaults to a vector of zeros.
        }
        \item{prior.var.beta}{
        A vector of prior variances for the regression parameters beta (Gaussian priors 
        are assumed). Defaults to a vector with values 1000.
        }  
        \item{prior.nu2}{
        The prior shape and scale in the form of c(shape, scale) for an Inverse-Gamma(shape, scale) 
        prior for the Gaussian error variance nu2. Defaults to c(1, 0.01) and only used if
        family=`Gaussian'.   
        }
        \item{prior.tau2}{
        The prior shape and scale in the form of c(shape, scale) for an Inverse-Gamma(shape, scale) 
        prior for the random effect variance tau2. Defaults to c(1, 0.01).  
        }
        \item{fix.rho.S}{
        Logical, should the spatial dependence parameter rho.S be fixed when fitting the model.  
        }
        \item{rho.S}{
        If fix.rho.S=TRUE, then the value in the interval [0,1] that rho.S should be fixed to.
        Not needed if fix.rho.S=FALSE.
        }
        \item{fix.rho.T}{
        Logical, should the temporal dependence parameter rho.T be fixed when fitting the model.  
        }
        \item{rho.T}{
        If fix.rho.T=TRUE, then the value in the interval [0,1] that rho.T should be fixed to.
        Not needed if fix.rho.T=FALSE.
        }
        \item{MALA}{
        Logical, should the function use MALA (TRUE, default) or simple random walk (FALSE) 
        updates for the random effects. Not applicable if family=`gaussian'.   
        }
        \item{verbose}{
        Logical, should the function update the user on its progress.  
        }
    }


\value{
    \item{summary.results }{A summary table of the parameters.}
    \item{samples }{A list containing the MCMC samples from the model.}
    \item{fitted.values }{A vector of fitted values for each area and time period.}
    \item{residuals }{A matrix with 3 columns where each column is a type of 
        residuals for each area and time period. The types are "response" (raw), 
        "pearson", and "deviance".}
        \item{modelfit }{Model fit criteria including the Deviance Information Criterion 
(DIC) and its corresponding estimated effective number of parameters (p.d), the Log 
Marginal Predictive Likelihood (LMPL), the Watanabe-Akaike Information Criterion 
(WAIC) and its corresponding estimated number of effective parameters (p.w), and the
estimated loglikelihood.}
    \item{localised.structure }{NULL, for compatability with the other models.}
    \item{formula }{The formula for the covariate and offset part of the model.}
    \item{model }{A text string describing the model fit.}
    \item{X }{The design matrix of covariates.}
}

\references{
    Leroux, B., X. Lei, and N. Breslow (2000). Estimation of disease rates in small 
    areas: A new mixed model for spatial dependence, Chapter Statistical Models in 
    Epidemiology, the Environment and Clinical Trials, Halloran, M and Berry, D (eds),
    pp. 135-178. Springer-Verlag, New York.
    
    Rushworth, A., D. Lee, and R. Mitchell (2014). A spatio-temporal model for estimating
    the long-term effects of air pollution on respiratory hospital admissions in Greater
    London. Spatial and Spatio-temporal Epidemiology 10, 29-38.

}

\author{
    Duncan Lee
}




\examples{
##################################################
#### Run the model on simulated data on a lattice
##################################################
#### set up the regular lattice    
x.easting <- 1:10
x.northing <- 1:10
Grid <- expand.grid(x.easting, x.northing)
K <- nrow(Grid)
N <- 10
N.all <- N * K

    
#### set up spatial neighbourhood matrix W
distance <- as.matrix(dist(Grid))
W <-array(0, c(K,K))
W[distance==1] <-1 	


#### Simulate the elements in the linear predictor and the data
gamma <- rnorm(n=N.all, mean=0, sd=0.001)
x <- rnorm(n=N.all, mean=0, sd=1)
beta <- 0.1
    
Q.W <- 0.99 * (diag(apply(W, 2, sum)) - W) + 0.01 * diag(rep(1,K))
Q.W.inv <- solve(Q.W)
phi.temp <- mvrnorm(n=1, mu=rep(0,K), Sigma=(0.1 * Q.W.inv))
phi <- phi.temp
    for(i in 2:N)
    {
    phi.temp2 <- mvrnorm(n=1, mu=(0.8 * phi.temp), Sigma=(0.1 * Q.W.inv))
    phi.temp <- phi.temp2
    phi <- c(phi, phi.temp)
    }
    
LP <- 3 + x * beta  + phi
mean <- exp(LP)
Y <- rpois(n=N.all, lambda=mean)
    

#### Run the model
\dontrun{model <- ST.CARar(formula=Y~x, family="poisson",  W=W, burnin=10000,
    n.sample=50000, fix.rho.T=TRUE, rho.T=0)}
}
