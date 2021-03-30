coef.CARBayesST <- function(object,...)
{
    #### Return the estimated regression coefficient
    if(is.null(nrow(object$samples$beta)))
    {
        return(NULL)  
    }else
    {
    beta <- apply(object$samples$beta, 2, median)
    names(beta) <- rownames(object$summary.results)[1:length(beta)]
    return(beta)
    }
}