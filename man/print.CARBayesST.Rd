\name{print.CARBayesST}
\alias{print.CARBayesST}


\title{
Print a summary of the fitted model to the screen.
}
\description{
This function takes a CARBayesST object and returns a summary of the fitted model.
The summary includes, for selected parameters, posterior medians and 95 percent 
credible intervals, the number of samples, the acceptance rate, the effective number 
of independent samples and the Geweke convergence diagnostic in the form of a Z-score.
}

\usage{
\method{print}{CARBayesST}(x, \dots)
}



\arguments{
  \item{x}{
A CARBayesST fitted model object.
}
  \item{\dots}{
Ignored.
}
}
\author{
Duncan Lee
}
