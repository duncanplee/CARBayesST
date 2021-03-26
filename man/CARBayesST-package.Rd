\name{CARBayesST-package}
\alias{CARBayesST-package}
\alias{CARBayesST}
\docType{package}

\title{
Spatio-Temporal Generalised Linear Mixed Models For Areal Unit Data
}

\description{
Implements a class of sunivariate and multivariate patio-temporal generalised 
linear mixed models for areal unit data, with inference in a Bayesian setting 
using Markov chain Monte Carlo (MCMC) simulation. The response variable can be 
binomial, Gaussian or Poisson, but for some models only the binomial and Poisson 
data likelihoods are available. The spatio-temporal autocorrelation is modelled 
by random effects,  which are assigned conditional autoregressive (CAR) style 
prior distributions. A number of different random effects structures are available, 
and full details are given in the vignette accompanying this package and the 
references below. The creation of this package was supported by the Engineering and 
Physical Sciences Research Council  (EPSRC) grant EP/J017442/1 and the Medical Research 
Council (MRC) grant MR/L022184/1.
}

\details{
\tabular{ll}{
Package: \tab CARBayesST\cr
Type: \tab Package\cr
Version: \tab 3.2\cr
Date: \tab 2021-03-30\cr
License: \tab GPL (>= 2)\cr
}
}

\author{
Author: Duncan Lee, Alastair Rushworth, Gary Napier and William Pettersson 

Maintainer: Duncan Lee <Duncan.Lee@glasgow.ac.uk>
}

\references{
Bernardinelli, L., D. Clayton, C.Pascuto, C.Montomoli, M.Ghislandi, and M. Songini
(1995). Bayesian analysis of space-time variation in disease risk. Statistics in 
Medicine, 14, 2433-2443.

Knorr-Held, L. (2000). Bayesian modelling of inseparable space-time variation in disease risk. 
Statistics in Medicine, 19, 2555-2567.

Lee, D and Lawson, C  (2016). Quantifying the spatial inequality and temporal trends in
maternal smoking rates in Glasgow, Annals of Applied Statistics, 10, 1427-1446.

Lee, D and Rushworth, A and Napier, G (2018). Spatio-Temporal Areal Unit Modeling in 
R with Conditional Autoregressive Priors Using the CARBayesST Package, 
Journal of Statistical Software, 84, 9, 1-39.

Lee, D and Meeks, K (2020). Improved inference for areal unit count data using graph-based 
optimisation. arXiv:2010.10893.

Napier, G, D. Lee, C. Robertson, A. Lawson, and K. Pollock (2016). A model to estimate the 
impact of changes in MMR vaccination uptake on inequalities in measles susceptibility in 
Scotland, Statistical Methods in Medical Research, 25, 1185-1200.

Napier, G., Lee, D., Robertson, C., and Lawson, A. (2019). A Bayesian space-time 
model for clustering areal units based on their disease trends, Biostatistics,
20, 681-697.

Rushworth, A., D. Lee, and R. Mitchell (2014). A spatio-temporal model for estimating the 
long-term effects of air pollution on respiratory hospital admissions in Greater London. 
Spatial and Spatio-temporal Epidemiology 10, 29-38.

Rushworth, A., Lee, D., and Sarran, C (2017).  An adaptive spatio-temporal smoothing model for 
estimating trends and step changes in disease risk. Journal of the Royal Statistical Society 
Series C, 66, 141-157.
}


\examples{
## See the examples in the function specific help files and in the vignette
## accompanying this package.
}
