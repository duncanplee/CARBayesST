\name{ST.CARar}
\alias{ST.CARar}


\title{
Fit a spatio-temporal generalised linear mixed model to data, with a spatio-temporal
autoregressive process.
}

\description{
Fit a spatio-temporal generalised linear mixed model to areal unit data, where the 
response variable can be binomial, Gaussian or Poisson. The linear predictor is 
modelled by known covariates and a vector of random effects. The latter follows either 
a multivariate first order autoregressive time series process or a multivariate
second order autoregressive time series process. In both cases the spatial 
autocorrelation is modelled via the precision matrix corresponding to the conditional 
autoregressive (CAR) prior proposed by Leroux et al. (2000), and the initial 
AR(1) model was proposed by Rushworth et al. (2014).  Missing values are 
allowed in the response in this model, and are sampled from in the model using 
data augmentation. Further details are given in the vignette accompanying this 
package. Inference is conducted in a Bayesian setting using Markov chain Monte 
Carlo (MCMC) simulation. 
}



\usage{
ST.CARar(formula, family, data=NULL,  trials=NULL, W, burnin, n.sample, thin=1, 
prior.mean.beta=NULL, prior.var.beta=NULL, prior.nu2=NULL, prior.tau2=NULL,
AR=NULL, rho.S=NULL, rho.T=NULL, MALA=FALSE, verbose=TRUE)
}


%- maybe also 'usage' for other objects documented here.
\arguments{
\item{formula}{
A formula for the covariate part of the model using the syntax of the lm() 
function. Offsets can be included here using the offset() function. The response 
and each covariate should be vectors of length (KN)*1, where K is the number of 
spatial units and N is the number of time periods. Each vector should be ordered 
so that the first K data points are the set of all K spatial locations at time 1, 
the next K are the set of spatial locations for time 2 and so on. The response can 
contain missing (NA) values.
}
\item{family}{
One of either "binomial", "gaussian" or "poisson", which respectively specify a 
binomial likelihood model with a logistic link function, a Gaussian likelihood 
model with an identity link function, or a Poisson likelihood model with a log 
link function.
}
\item{data}{
An optional data.frame containing the  variables in the formula.
}
\item{trials}{
A vector the same length and in the same order as the response containing the 
total number of trials for each area and time period. Only used if family="binomial". 
}
\item{W}{A non-negative K by K neighbourhood matrix (where K is the number of 
spatial units). Typically a binary specification is used, where the jkth 
element equals one if areas (j, k) are spatially close (e.g. share a common 
border) and is zero otherwise. The matrix can be non-binary, but each row must 
contain at least one non-zero entry.
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
are assumed). Defaults to a vector with values 100000.
}   
\item{prior.nu2}{
The prior shape and scale in the form of c(shape, scale) for an Inverse-Gamma(shape, scale) 
prior for nu2. Defaults to c(1, 0.01) and only used if family="Gaussian".   
}
\item{prior.tau2}{
The prior shape and scale in the form of c(shape, scale) for an Inverse-Gamma(shape, scale) 
prior for tau2. Defaults to c(1, 0.01).  
}
\item{AR}{
The order of the autoregressive time series process that must be either 1 or 2. 
}
\item{rho.S}{
The value in the interval [0, 1] that the spatial dependence parameter rho.S is 
fixed at if it should not be estimated. If this arugment is NULL then rho.S is
estimated in the model. 
}
\item{rho.T}{
Whether to fix or estimate the temporal dependence parameter(s) rho.T in the model. 
If this arugment is NULL then they are estimated in the model. If you want to fix
them and AR=1 then it must be a value in the interval [0, 1]. If  AR=2 then it must be
a vector of length two with the first and second order autoregressive coefficients.
}
\item{MALA}{
Logical, should the function use Metropolis adjusted Langevin algorithm (MALA)
updates (TRUE) or simple random walk (FALSE, default) updates for the regression 
parameters. Not applicable if family="gaussian".    
}
\item{verbose}{
Logical, should the function update the user on its progress.  
}
}


\value{
\item{summary.results }{A summary table of the parameters.}
\item{samples }{A list containing the MCMC samples from the model.}
\item{fitted.values }{A vector of fitted values for each area and time period.}
\item{residuals }{A matrix with 2 columns where each column is a type of 
residual and each row relates to an area and time period. The types are "response" 
(raw), and "pearson".}
\item{modelfit }{Model fit criteria including the Deviance Information Criterion 
(DIC) and its corresponding estimated effective number of parameters (p.d), the Log 
Marginal Predictive Likelihood (LMPL), the Watanabe-Akaike Information Criterion 
(WAIC) and its corresponding estimated number of effective parameters (p.w), and
the loglikelihood.}
\item{accept }{The acceptance probabilities for the parameters.}
\item{localised.structure }{NULL, for compatability with the other models.}
\item{formula }{The formula (as a text string) for the response, covariate and 
offset parts of the model.}
\item{model }{A text string describing the model fit.}
\item{X }{The design matrix of covariates.}
}

\references{
Rushworth, A., D. Lee, and R. Mitchell (2014). A spatio-temporal model for estimating
the long-term effects of air pollution on respiratory hospital admissions in Greater
London. Spatial and Spatio-temporal Epidemiology 10, 29-38.
}

\author{
Duncan Lee
}




\examples{
#################################################
#### Run the model on simulated data on a lattice
#################################################
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
    n.sample=50000, AR=1)}
    
    
#### Toy example for checking  
model <- ST.CARar(formula=Y~x, family="poisson",  W=W, burnin=10,
    n.sample=50, AR=1)
}
