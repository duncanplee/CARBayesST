\name{ST.CARadaptive}
\alias{ST.CARadaptive}


\title{
Fit a spatio-temporal generalised linear mixed model to data, with a spatio-temporal
autoregressive process that has an adaptive autocorrelation stucture.
}

\description{
Fit a spatio-temporal generalised linear mixed model to areal unit data, where 
the response variable can be binomial, Gaussian or Poisson. The linear predictor 
is modelled by known covariates and a vector of random effects. The latter follows 
a multivariate first order autoregressive time series process, where spatial 
autocorrelation is modelled via the precision matrix, which is based on a CAR type
structure and a neighbourhood (adjacency) matrix W. The non-zero elements of W 
associated with each pair of geographically adjacent areal units are treated as 
random variables with ranges in the unit interval, which allows step changes to 
be identified in the random effects surface between geographically adjacent 
regions. The model is similar to that proposed by Rushworth et al. (2017). Further 
details are given in the vignette accompanying this package. Inference is conducted 
in a Bayesian setting using Markov chain Monte Carlo (MCMC) simulation. 
}


\usage{
ST.CARadaptive(formula, family, data=NULL, trials=NULL, W, burnin, n.sample, thin=1,  
prior.mean.beta=NULL, prior.var.beta=NULL, prior.nu2=NULL, prior.tau2=NULL, 
rho = NULL, epsilon = 0, MALA=FALSE, verbose=TRUE) 
}


\arguments{
\item{formula}{
A formula for the covariate part of the model using the syntax of the lm() 
function. Offsets can be included here using the offset() function. The response 
and each covariate should be vectors of length (KN)*1, where K is the number of 
spatial units and N is the number of time periods. Each vector should be ordered 
so that the first K data points are the set of all K spatial locations at time 1, 
the next K are the set of spatial locations for time 2 and so on. The response 
must NOT contain missing (NA) values.
}
\item{family}{
One of either "binomial", "gaussian", or "poisson", which respectively specify a 
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
\item{W}{
A non-negative K by K neighbourhood matrix (where K is the number of spatial units). 
Typically a binary specification is used, where the jkth element equals one if 
areas (j, k) are spatially close (e.g. share a common border) and is zero otherwise. 
For this model the matrix must be binary.
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
\item{rho}{
The value in the interval [0, 1] that the spatial dependence parameter rho is 
fixed at if it should not be estimated. If this arugment is NULL then rho is
estimated in the model. Setting rho=1, reduces the random effects prior to the 
intrinsic CAR model but does require epsilon>0.
}
\item{epsilon}{
Diagonal ridge parameter to add to the random effects prior precision matrix, only 
required when rho = 1, and the prior precision is improper.  Defaults to 0.
}
\item{MALA}{
Logical, should the function use Metropolis adjusted Langevin algorithm (MALA)
updates (TRUE) or simple random walk (FALSE, default) updates for the random 
effects and the regression parameters. Not applicable if family="gaussian".   
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
\item{localised.structure }{A list with 2 K*K matrices, Wmedian and W99 summarising 
the estimated adjacency relationships. Wmedian contains the posterior median for 
each w_ij element estimated in the model for adjacent areal units, while W99 
contains binary indicator variables for whether Prob(w_ij < 0.5|data)>0.99. 
For both matrices, elements corresponding to non-adjacent pairs of areas have 
NA values.}
\item{formula }{The formula (as a text string) for the response, covariate and 
offset parts of the model.}
\item{model }{A text string describing the model fit.}
\item{X }{The design matrix of covariates.}
}




\references{
Rushworth, A., Lee, D., and Sarran, C (2017).  An adaptive spatio-temporal smoothing 
model for estimating trends and step changes in disease risk. Journal of the Royal 
Statistical Society Series C, 66, 141-157.
}

\author{
Alastair Rushworth
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
jump <- rep(c(rep(2, 70), rep(4, 30)),N)
LP <- jump + phi
fitted <- exp(LP)
Y <- rpois(n=N.all, lambda=fitted)


#### Run the model     
\dontrun{model <- ST.CARadaptive(formula=Y~1, family="poisson", W=W, burnin=10000,
n.sample=50000)}


#### Toy example for checking    
model <- ST.CARadaptive(formula=Y~1, family="poisson", W=W, burnin=10,
n.sample=50)
}