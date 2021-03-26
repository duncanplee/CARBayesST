\name{ST.CARsepspatial}
\alias{ST.CARsepspatial}


\title{
Fit a spatio-temporal generalised linear mixed model to data, with a common 
temporal main effect and separate spatial surfaces with individual variances.
}

\description{
Fit a spatio-temporal generalised linear mixed model to areal unit data, where the 
response variable can be binomial or Poisson. The linear predictor is modelled by
known covariates and two sets of random effects. These include a common temporal 
main effect, and separate time period specific spatial effects with a common 
spatial dependence parameter but separate variance parameters. Each component is 
modelled by the conditional autoregressive (CAR) prior proposed by Leroux et al. 
(2000). Further details are given in Napier et al. (2016) and in the vignette 
accompanying this package. Inference is conducted in a Bayesian setting using 
Markov chain Monte Carlo (MCMC) simulation. 
}

\usage{
ST.CARsepspatial(formula, family, data=NULL,  trials=NULL, W, burnin, n.sample,
thin=1, prior.mean.beta=NULL, prior.var.beta=NULL, prior.tau2=NULL, rho.S=NULL, 
rho.T=NULL, MALA=FALSE, verbose=TRUE)
}



%- maybe also 'usage' for other objects documented here.
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
One of either "binomial" or "poisson", which respectively specify a binomial 
likelihood model with a logistic link function, or a Poisson likelihood model with
a log link function.
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
\item{prior.tau2}{
The prior shape and scale in the form of c(shape, scale) for an Inverse-Gamma(shape, scale) 
prior for tau2. Defaults to c(1, 0.01).  
}
\item{rho.S}{
The value in the interval [0, 1] that the spatial dependence parameter rho.S is 
fixed at if it should not be estimated. If this arugment is NULL then rho.S is
estimated in the model. 
}
\item{rho.T}{
The value in the interval [0, 1] that the temporal dependence parameter rho.T is 
fixed at if it should not be estimated. If this arugment is NULL then rho.T is
estimated in the model. 
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
Napier, G, D. Lee, C. Robertson, A. Lawson, and K. Pollock (2016). A model to 
estimate the impact of changes in MMR vaccination uptake on inequalities in 
measles susceptibility in Scotland, Statistical Methods in Medical Research, 25, 
1185-1200.
}
        
\author{
Gary Napier
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
N <- 5
N.all <- N * K
  
        
#### set up spatial neighbourhood matrix W
distance <- as.matrix(dist(Grid))
W <-array(0, c(K,K))
W[distance==1] <-1 	


#### Create the spatial covariance matrix
Q.W <- 0.99 * (diag(apply(W, 2, sum)) - W) + 0.01 * diag(rep(1,K))
Q.W.inv <- solve(Q.W)
  
           
#### Simulate the elements in the linear predictor and the data
x <- rnorm(n=N.all, mean=0, sd=1)
beta <- 0.1

phi1 <- mvrnorm(n=1, mu=rep(0,K), Sigma=(0.01 * Q.W.inv))
phi2 <- mvrnorm(n=1, mu=rep(0,K), Sigma=(0.01 * Q.W.inv))
phi3 <- mvrnorm(n=1, mu=rep(0,K), Sigma=(0.01 * Q.W.inv))
phi4 <- mvrnorm(n=1, mu=rep(0,K), Sigma=(0.05 * Q.W.inv))
phi5 <- mvrnorm(n=1, mu=rep(0,K), Sigma=(0.05 * Q.W.inv))
  
delta <- c(0, 0.5, 0, 0.5, 0)
phi.long <- c(phi1, phi2, phi3, phi4, phi5)
delta.long <- kronecker(delta, rep(1,K))
LP <- 4 +  x * beta + phi.long +  delta.long
mean <- exp(LP)
Y <- rpois(n=N.all, lambda=mean)
  
                
#### Run the model
\dontrun{model <- ST.CARsepspatial(formula=Y~x, family="poisson", W=W, burnin=10000, 
n.sample=50000)}


#### Toy example for checking
model <- ST.CARsepspatial(formula=Y~x, family="poisson", W=W, burnin=10, 
n.sample=50)
}
        