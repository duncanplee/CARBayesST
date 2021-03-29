\name{MVST.CARar}
\alias{MVST.CARar}


\title{
Fit a multivariate spatio-temporal generalised linear mixed model to data, with a multivariate spatio-temporal autoregressive process.
}


\description{
Fit a multivariate spatio-temporal generalised linear mixed model to 
multivariate areal unit data, where the response variable can be binomial, 
Gaussian or Poisson. The linear predictor is modelled by known covariates and a 
vector of random effects. The latter allows for correlations over: (i) K areal
units; (ii) N time periods; and (iii) J outcomes. These random effects are 
modelled by either a multivariate first order autoregressive time series 
process or a multivariate second order autoregressive time series process. 
In both cases the spatial and between outcome correlation is modelled via the 
precision matrix, and the spatial correlation is represented by the 
conditional autoregressive (CAR) prior proposed by Leroux et al. (2000). In 
contrast, the between outcome correlation structure is estimated from the data,
and no prior form is assumed. Missing values are allowed in the response in this 
model, and are sampled from in the model using data augmentation. Further details are 
given in the vignette accompanying this package. Inference is conducted in a 
Bayesian setting using Markov chain Monte Carlo (MCMC) simulation.
}


\usage{
MVST.CARar(formula, family, data=NULL,  trials=NULL, W, burnin, n.sample, thin=1, 
prior.mean.beta=NULL, prior.var.beta=NULL, prior.nu2=NULL, prior.Sigma.df=NULL, 
prior.Sigma.scale=NULL, AR=NULL, rho.S=NULL, rho.T=NULL, MALA=FALSE, verbose=TRUE)
}


\arguments{
\item{formula}{
A formula for the covariate part of the model using the syntax of the lm() 
function. Offsets can be included here using the offset() function. The response 
and the offset (if included) should be matrices of dimension (KN)*J, where K is 
the number of spatial units, N is the number of time periods and J is the number 
of different variables. Each column of the response and offset matrices relates to
a different outcome variable. The values in each column of these matrices should be 
ordered so that the first K data points are the set of all K spatial locations at 
time 1, the next K are the set of spatial locations for time 2 and so on. The 
covariates should each be a (KN)*1 vector, and different regression parameters are 
estimated for each of the J variables. The response can contain missing (NA) values. 
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
A (KN)*J matrix of the same dimension as the response. Only used if 
family="binomial". 
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
\item{prior.Sigma.df}{
    The prior degrees of freedom for the Inverse-Wishart prior for Sigma. 
    Defaults to J+1.   
}
\item{prior.Sigma.scale}{
    The prior J times J scale matrix for the Inverse-Wishart prior for Sigma. 
    Defaults to the identity matrix divided by 1000.
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
them and AR=1 then it must be a single value. If  AR=2 then it must be
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
\item{fitted.values }{A matrix of fitted values for each area, time period and response 
variable in the same order as the response variable.}
\item{residuals }{A list with 2 elements, where each element is a matrix of a type  
of residual. Each row of a matrix relates to an area and time period and each column to a 
response (category). The types of residual are "response" (raw), and "pearson".
}
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


\author{
Duncan Lee
}


\examples{
#################################################
#### Run the model on simulated data on a lattice
#################################################
#### Set up a square lattice region
x.easting <- 1:8
x.northing <- 1:8
Grid <- expand.grid(x.easting, x.northing)


#### Set up the coordinate dimensions
K <- nrow(Grid)
N <- 15
J <- 2
N.all <- N * K * J


#### set up distance and neighbourhood (W, based on sharing a common border) matrices
distance <- as.matrix(dist(Grid))
W <-array(0, c(K,K))
W[distance==1] <-1 	


#### Set up the spatial covariance matrix
Q.W <- 0.8 * (diag(apply(W, 2, sum)) - W) + 0.2 * diag(rep(1,K))
Q.W.inv <- solve(Q.W)


#### Set up the multivariate outcome covariance matrix
Sigma <- 0.01 * array(c(1, 1, 1, 2), c(2,2))
Sigma.inv <- solve(Sigma)


#### Spatial and between outcome covariance
QSig.prec <- kronecker(Q.W, Sigma.inv)
QSig.var <-solve(QSig.prec)


#### Generate the covariate
x1 <- rnorm(n=N * K, mean=0, sd=1)
lp.regression.mat <- cbind(0.1 + 0.1 * x1, 0.1 - 0.1*x1)
lp.regression <- as.numeric(t(lp.regression.mat))


#### Spatio-temporal random effects
phi.temp <- mvrnorm(n=1, mu=rep(0,(J*K)), Sigma=QSig.var)
phi <- phi.temp
    for(i in 2:N)
    {
    phi.temp2 <- mvrnorm(n=1, mu=(0.8 * phi.temp), Sigma=QSig.var)
    phi.temp <- phi.temp2
    phi <- c(phi, phi.temp)
    }
phi <- phi - mean(phi)
phi.true <- matrix(phi, ncol=2, byrow=TRUE)


#### Generate the binomial counts
lp <- lp.regression + phi
p <- exp(lp) / (1+exp(lp))
trials <- rpois(N.all, lambda=100)
Y <- rbinom(n=N.all, size=trials, prob=p)
Y.mat <- matrix(Y, nrow=(K*N), ncol=J, byrow=TRUE)
trials.mat <- matrix(trials, nrow=(K*N), ncol=J, byrow=TRUE)
formula <- Y.mat~x1

#### Run the model
formula <- Y.mat ~ x1
\dontrun{mod <- MVST.CARar(formula=formula, family="binomial", trials=trials.mat, W=W, 
burnin=10000, n.sample=50000, prior.Sigma.scale=matrix(c(0.01, 0, 0, 0.01), nrow=2), 
AR=1, MALA=FALSE)}

#### Toy example for checking
mod <- MVST.CARar(formula=formula, family="binomial", trials=trials.mat, W=W, 
burnin=10, n.sample=50, prior.Sigma.scale=matrix(c(0.01, 0, 0, 0.01), nrow=2), 
AR=1, MALA=FALSE)
}
