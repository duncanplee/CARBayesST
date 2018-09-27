\name{ST.CARclustrends}
\alias{ST.CARclustrends}


\title{
Fit a spatio-temporal generalised linear mixed model to data, with a clustering 
of temporal trend functions and a temporally common spatial surface.
}

\description{
Fit a spatio-temporal generalised linear mixed model to areal unit data, where the 
response variable can be binomial or Poisson. The linear predictor is modelled by 
known covariates, a temoporally common spatial surface, and a mixture of temporal 
trend functions. The spatial component is modelled by the conditional autoregressive 
(CAR) prior proposed by Leroux et al. (2000). The temporal trend functions are 
user-specified and are fixed parametric forms (e.g. linear, step-change) or 
constrained shapes (e.g. monotonically increasing). Further details are given 
in Napier et al. (2018) and in the vignette accompanying this package. Inference 
is conducted in a Bayesian setting using Metropolis coupled Markov chain Monte 
Carlo (MCMCMC) simulation. 
}

\usage{
    ST.CARclustrends(formula, family, data=NULL, trials=NULL, W, burnin, n.sample,
    thin=1, trends=NULL, changepoint=NULL, knots=NULL, prior.mean.beta=NULL,
    prior.var.beta=NULL, prior.mean.gamma=NULL, prior.var.gamma=NULL,
    prior.lambda=NULL, prior.tau2=NULL, Nchains=4, verbose=TRUE)
}


%- maybe also 'usage' for other objects documented here.
\arguments{
\item{formula}{
A formula using the syntax of the lm() function. However, due to identifiability 
issues covariates are not allowed. So the only elements allowed on the right side
of the formula are an intercept term and an offset term using the offset() 
function. The response variable and the offset (if specified) should be vectors 
of length (KN)*1, where K is the number of spatial units and N is the number of 
time periods. Each vector should be ordered so that the first K data points are 
the set of all K spatial locations at time 1, the next K are the set of spatial 
locations for time 2 and so on.
}
\item{family}{
One of either "binomial" or"poisson", which respectively specify a binomial 
likelihood model with a logistic link function, or a Poisson likelihood model 
with a log link function.
}
\item{data}{
An optional data.frame containing the variables in the formula.
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
\item{trends}{
A vector containing the temporal trend functions to include in the model, which 
include: constant ("Constant""); linear decreasing ("LD"); linear increasing 
("LI"); Known change point, where the trend can increase towards the change point 
before subsequently decreasing ("CP"); or decrease towards the change point before 
subsequently increasing ("CT"); and monotonic cubic splines which are decreasing 
("MD") or increasing ("MI"). At least two trends have to be selected, with the 
constant trend always included. To avoid identifiability problems only one of "LI" 
or "MI" can be included at a given time (similarily for "LD" and "MD").
}
\item{changepoint}{
A scalar indicating the position of the change point should one of the change point 
trend functions be included in the trends vector, i.e. if "CP" or "CT" is specified.
}
\item{knots}{
A scalar indicating the number of knots to use should one of the monotonic cubic 
splines trend functions be included in the trends vector, i.e. if "MD" or "MI" is 
specified.
}
\item{prior.mean.beta}{
A vector of prior means for the regression parameters beta (Gaussian priors are 
assumed). Defaults to a vector of zeros.
}
\item{prior.var.beta}{
A vector of prior variances for the regression parameters beta (Gaussian priors 
are assumed). Defaults to a vector with values 100000.
}  
\item{prior.mean.gamma}{
A vector of prior means for the temporal trend parameters (Gaussian priors are 
assumed). Defaults to a vector of zeros.
}
\item{prior.var.gamma}{
A vector of prior variances for the temporal trend parameters (Gaussian priors 
are assumed). Defaults to a vector with values 100000.
}  
\item{prior.lambda}{
A vector of prior samples sizes for the Dirichlet prior controlling the 
probabilities that each trend function is chosen. The vector should be the same 
length as the trends vector and defaults to a vector of ones.
}
\item{prior.tau2}{
The prior shape and scale in the form of c(shape, scale) for an Inverse-Gamma(shape, scale) 
prior for the random effect variances tau2. Defaults to c(1, 0.01).  
}
\item{Nchains}{
The number of parallel Markov chains to be used in the Metropolis coupled Markov 
chain Monte Carlo (MCMCMC) simulations. Defaults to 4.
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
\item{localised.structure }{A list containing two elements. The first is "trends",
which is a vector the same length and in the same order as the number of areas. 
The kth element specifies which trend area k has been allocated to based on the 
posterior mode. The second element is "trend.probs", which is a matrix containing 
the probabilities associated with each trend for each areal unit.}
\item{formula }{The formula for the covariate and offset parts of the model.}
\item{model }{A text string describing the model fit.}
\item{X }{The design matrix of covariates.}
}





        
\references{
Leroux, B., Lei, X., and Breslow, N. (2000). Estimation of disease rates in small 
areas: A new mixed model for spatial dependence, Chapter Statistical Models in 
Epidemiology, the Environment and Clinical Trials, Halloran, M and Berry, D (eds),
pp. 135-178. Springer-Verlag, New York.

Napier, G., Lee, D., Robertson, C., and Lawson, A. (2018). A Bayesian space-time 
model for clustering areal units based on their disease trends, Biostatistics,
to appear, 10.1093/biostatistics/kxy024.
}
        
\author{
Gary Napier
}
        
        
        
        
\examples{
##################################################
#### Run the model on simulated data on a lattice
##################################################
#### Load the libraries required
library(truncdist)
library(gtools)


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


#### Create the spatial covariance matrix
Q.W <- 0.99 * (diag(apply(W, 2, sum)) - W) + 0.01 * diag(rep(1,K))
Q.W.inv <- solve(Q.W)
  
           
#### Simulate the elements in the linear predictor and the data
beta <- 0.01
gamma <- 0.7
phi <- mvrnorm(n=1, mu=rep(0,K), Sigma=(0.01 * Q.W.inv))

lambda <- rep(1/2, 2)
w <- t(rmultinom(K, 1, lambda))

Y <- matrix(NA, nrow = K, ncol = N)
for (i in 1:N)
{
  LP <- beta + w[, 2] * (gamma * i) + phi
  mean <- exp(LP)
  Y[, i] <- rpois(n=K, lambda=mean)
 }
Y <- as.numeric(Y)


#### Run the model
\dontrun{model <- ST.CARclustrends(formula=Y~1, family="poisson", W=W, burnin=10000, 
n.sample=50000, trends=c("Constant", "LI"))}


#### Toy example for checking
model <- ST.CARclustrends(formula=Y~1, family="poisson", W=W, burnin=10, 
n.sample=50, trends=c("Constant", "LI"))
}