residuals.CARBayesST <- function(object, type="deviance", ...)
{
    #### Return one of three types of residuals
    if(type=="response")
    {
        return(object$residuals$response)
    }else if(type=="pearson")
    {
        return(object$residuals$pearson)
    }else if(type=="deviance")
    {
        return(object$residuals$deviance)
    }else
    {
        return("Error. That is not one of the allowable residual types.")   
    }
}