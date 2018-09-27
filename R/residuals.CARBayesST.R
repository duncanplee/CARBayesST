residuals.CARBayesST <- function(object, type="pearson", ...)
{
    #### Return one of two types of residuals
    if(type=="response")
    {
        return(object$residuals$response)
    }else if(type=="pearson")
    {
        return(object$residuals$pearson)
    }else
    {
        return("Error. That is not one of the allowable residual types.")   
    }
}