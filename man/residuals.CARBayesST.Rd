\name{residuals.CARBayesST}
\alias{residuals.CARBayesST}


\title{
    Extract the residuals from a model.
}
\description{
    This function takes a CARBayesST object and returns a set of residuals. 
    The allowable types of residual are "response" (raw), and "pearson" (the 
    default). In each case the fitted values are based on posterior means.
}

\usage{
    \method{residuals}{CARBayesST}(object, type, \dots)
}



\arguments{
    \item{object}{
        A CARBayesST fitted model object.
    }
    \item{type}{
        A text string and one of c("response", "pearson"). If this
        argument is omitted the default is "pearson".
    }
    \item{\dots}{
        Ignored.
    }
}
\author{
    Duncan Lee
}
