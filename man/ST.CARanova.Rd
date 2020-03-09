\name{ST.CARanova}
\alias{ST.CARanova}


\title{
Fit a spatio-temporal generalised linear mixed model to data, with spatial and  temporal main
effects and a spatio-temporal interaction.
}

\description{
Fit a spatio-temporal generalised linear mixed model to areal unit data, where the 
response variable can be binomial, Gaussian or Poisson. The linear predictor is 
modelled by known covariates and three vectors of random effects. The latter 
include spatial and temporal main effects and a spatio-temporal interaction. 
The spatial and temporal main effects are modelled by the conditional 
autoregressive (CAR) prior proposed by Leroux et al. (2000), while the 
spatio-temporal interaction random effects are independent. Due to the lack of 
identifiability between the interactions and the Gaussian errors, only main 
effects  are allowed in the Gaussian model. Missing values are allowed in the 
response in this model, and are sampled from in the model using data augmentation.
The model is similar to that proposed by Knorr-Held (2000), and further details 
are given in the vignette accompanying this package. Inference is conducted in a 
Bayesian setting using Markov chain Monte Carlo (MCMC) simulation. 
}

\usage{
ST.CARanova(formula, family, data=NULL,  trials=NULL, W, interaction=TRUE, burnin, 
n.sample, thin=1, prior.mean.beta=NULL, prior.var.beta=NULL, prior.nu2=NULL, 
prior.tau2=NULL, rho.S=NULL, rho.T=NULL, MALA=FALSE, verbose=TRUE)
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
\item{interaction}{
TRUE or FALSE indicating whether the spatio-temporal interaction random effects
should be included. Defaults to TRUE unless family="gaussian" in which case 
interactions are not allowed.
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
\item{localised.structure }{NULL, for compatability with the other models.}
\item{formula }{The formula (as a text string) for the response, covariate and 
offset parts of the model.}
\item{model }{A text string describing the model fit.}
\item{X }{The design matrix of covariates.}
}



\references{
Knorr-Held, L. (2000). Bayesian modelling of inseparable space-time variation in 
disease risk. Statistics in Medicine, 19, 2555-2567.

Leroux, B., X. Lei, and N. Breslow (2000). Estimation of disease rates in small 
areas: A new mixed model for spatial dependence, Chapter Statistical Models in 
Epidemiology, the Environment and Clinical Trials, Halloran, M and Berry, D (eds),
pp. 135-178. Springer-Verlag, New York.
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


#### set up spatial (W) and temporal (D) neighbourhood matrices
distance <- as.matrix(dist(Grid))
W <-array(0, c(K,K))
W[distance==1] <-1 	
	
D <-array(0, c(N,N))
for(i in 1:N)
{
    for(j in 1:N)
    {
        if(abs((i-j))==1)  D[i,j] <- 1 
    }	
}


#### Simulate the elements in the linear predictor and the data
gamma <- rnorm(n=N.all, mean=0, sd=0.001)
x <- rnorm(n=N.all, mean=0, sd=1)
beta <- 0.1

Q.W <- 0.99 * (diag(apply(W, 2, sum)) - W) + 0.01 * diag(rep(1,K))
Q.W.inv <- solve(Q.W)
phi <- mvrnorm(n=1, mu=rep(0,K), Sigma=(0.01 * Q.W.inv))

Q.D <- 0.99 * (diag(apply(D, 2, sum)) - D) + 0.01 * diag(rep(1,N))
Q.D.inv <- solve(Q.D)
delta <- mvrnorm(n=1, mu=rep(0,N), Sigma=(0.01 * Q.D.inv))

phi.long <- rep(phi, N)
delta.long <- kronecker(delta, rep(1,K))
LP <- 4 +  x * beta + phi.long +  delta.long + gamma
mean <- exp(LP)
Y <- rpois(n=N.all, lambda=mean)


#### Run the model
\dontrun{model <- ST.CARanova(formula=Y~x, family="poisson", interaction=TRUE, 
W=W, burnin=10000, n.sample=50000)}


#### Toy example for checking
model <- ST.CARanova(formula=Y~x, family="poisson", interaction=TRUE, 
W=W, burnin=10, n.sample=50)
}
