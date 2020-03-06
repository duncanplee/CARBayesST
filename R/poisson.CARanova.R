poisson.CARanova <- function(formula, data=NULL, W, interaction=TRUE, burnin, n.sample, thin=1,  prior.mean.beta=NULL, prior.var.beta=NULL, prior.tau2=NULL, rho.S=NULL, rho.T=NULL, MALA=FALSE, verbose=TRUE)
{
##############################################
#### Format the arguments and check for errors
##############################################
#### Verbose
a <- common.verbose(verbose)  
     
     
#### Frame object
frame.results <- common.frame(formula, data, "poisson")
N.all <- frame.results$n
p <- frame.results$p
X <- frame.results$X
X.standardised <- frame.results$X.standardised
X.sd <- frame.results$X.sd
X.mean <- frame.results$X.mean
X.indicator <- frame.results$X.indicator 
offset <- frame.results$offset
Y <- frame.results$Y
which.miss <- frame.results$which.miss
n.miss <- frame.results$n.miss  
Y.DA <- Y  


#### Check on MALA argument
    if(length(MALA)!=1) stop("MALA is not length 1.", call.=FALSE)
    if(!is.logical(MALA)) stop("MALA is not logical.", call.=FALSE) 


#### Check on the rho arguments
    if(is.null(rho.S))
    {
    rho <- runif(1)
    fix.rho.S <- FALSE   
    }else
    {
    rho <- rho.S
    fix.rho.S <- TRUE
    }
    if(!is.numeric(rho)) stop("rho.S is fixed but is not numeric.", call.=FALSE)  
    if(rho<0 ) stop("rho.S is outside the range [0, 1].", call.=FALSE)  
    if(rho>1 ) stop("rho.S is outside the range [0, 1].", call.=FALSE)    

    if(is.null(rho.T))
    {
    lambda <- runif(1)
    fix.rho.T <- FALSE   
    }else
    {
    lambda <- rho.T
    fix.rho.T <- TRUE
    }
    if(!is.numeric(lambda)) stop("rho.T is fixed but is not numeric.", call.=FALSE)  
    if(lambda<0 ) stop("rho.T is outside the range [0, 1].", call.=FALSE)  
    if(lambda>1 ) stop("rho.T is outside the range [0, 1].", call.=FALSE)  


#### CAR quantities
W.quants <- common.Wcheckformat.leroux(W)
K <- W.quants$n
N <- N.all / K
W <- W.quants$W
W.triplet <- W.quants$W.triplet
W.n.triplet <- W.quants$n.triplet
W.triplet.sum <- W.quants$W.triplet.sum
n.neighbours <- W.quants$n.neighbours 
W.begfin <- W.quants$W.begfin



#### Checks on the interaction flag
    if(sum(interaction==c(TRUE, FALSE)) != 1) stop("interaction must be either TRUE or FALSE.", call.=FALSE)
    if(length(interaction) != 1) stop("interaction must be of length 1.", call.=FALSE)


#### Priors
    if(is.null(prior.mean.beta)) prior.mean.beta <- rep(0, p)
    if(is.null(prior.var.beta)) prior.var.beta <- rep(100000, p)
    if(is.null(prior.tau2)) prior.tau2 <- c(1, 0.01)
prior.beta.check(prior.mean.beta, prior.var.beta, p)
prior.var.check(prior.tau2)


## Compute the blocking structure for beta     
block.temp <- common.betablock(p)
beta.beg  <- block.temp[[1]]
beta.fin <- block.temp[[2]]
n.beta.block <- block.temp[[3]]
list.block <- as.list(rep(NA, n.beta.block*2))
    for(r in 1:n.beta.block)
    {
    list.block[[r]] <- beta.beg[r]:beta.fin[r]-1
    list.block[[r+n.beta.block]] <- length(list.block[[r]])
    }


#### MCMC quantities - burnin, n.sample, thin
common.burnin.nsample.thin.check(burnin, n.sample, thin)



#############################
#### Initial parameter values
#############################
mod.glm <- glm(Y~X.standardised-1, offset=offset, family="quasipoisson")
beta.mean <- mod.glm$coefficients
beta.sd <- sqrt(diag(summary(mod.glm)$cov.scaled))
beta <- rnorm(n=length(beta.mean), mean=beta.mean, sd=beta.sd)

log.Y <- log(Y)
log.Y[Y==0] <- -0.1  
res.temp <- log.Y - X.standardised %*% beta - offset
res.sd <- sd(res.temp, na.rm=TRUE)/5
phi <- rnorm(n=K, mean=0, sd = res.sd)
delta <- rnorm(n=N, mean=0, sd = res.sd)
tau2.phi <- var(phi)/10
tau2.delta <- var(delta)/10
    if(interaction)
    {
    gamma <- rnorm(n=N.all, mean=0, sd = res.sd)
    tau2.gamma <- var(gamma)/10
    }else
    {}


#### Matrix versions
offset.mat <- matrix(offset, nrow=K, ncol=N, byrow=FALSE) 
regression.mat <- matrix(X.standardised %*% beta, nrow=K, ncol=N, byrow=FALSE)   
phi.mat <- matrix(rep(phi, N), byrow=F, nrow=K)
delta.mat <- matrix(rep(delta, K), byrow=T, nrow=K)

    if(interaction)
    {
    gamma.mat <- matrix(gamma, byrow=F, nrow=K)
    }else
    {
    gamma.mat <- matrix(rep(0, N.all), byrow=F, nrow=K)
    }
fitted <- exp(as.numeric(offset.mat + regression.mat + phi.mat  + delta.mat + gamma.mat))



###############################    
#### Set up the MCMC quantities    
###############################
#### Matrices to store samples
n.keep <- floor((n.sample - burnin)/thin)
samples.beta <- array(NA, c(n.keep, p))
samples.phi <- array(NA, c(n.keep, K))
samples.delta <- array(NA, c(n.keep, N))
    if(!fix.rho.S) samples.rho <- array(NA, c(n.keep, 1))
    if(!fix.rho.T) samples.lambda <- array(NA, c(n.keep, 1))
samples.fitted <- array(NA, c(n.keep, N.all))
samples.loglike <- array(NA, c(n.keep, N.all))
    if(n.miss>0) samples.Y <- array(NA, c(n.keep, n.miss))
    if(interaction)
    {
    samples.gamma <- array(NA, c(n.keep, N.all))
    samples.tau2 <- array(NA, c(n.keep, 3))
    colnames(samples.tau2) <- c("tau2.phi", "tau2.delta", "tau2.gamma")
    }else
    {
    samples.tau2 <- array(NA, c(n.keep, 2))
    colnames(samples.tau2) <- c("tau2.phi", "tau2.delta")    
    }


#### Specify the Metropolis quantities
accept.all <- rep(0,12)
accept <- accept.all
proposal.sd.beta <- 0.01
proposal.sd.phi <- 0.1
proposal.sd.delta <- 0.1
proposal.sd.rho <- 0.02
proposal.sd.lambda <- 0.02
tau2.phi.shape <- prior.tau2[1] + K/2
tau2.delta.shape <- prior.tau2[1] + N/2
    if(interaction)
    {
    proposal.sd.gamma <- 0.1
    tau2.gamma.shape <- prior.tau2[1] + N*K/2
    }else
    {
    }



##########################################
#### Specify spatial and temporal elements
##########################################
#### Spatial determinant
    if(!fix.rho.S) 
    {
    Wstar <- diag(apply(W,1,sum)) - W
    Wstar.eigen <- eigen(Wstar)
    Wstar.val <- Wstar.eigen$values
    det.Q.W <-  0.5 * sum(log((rho * Wstar.val + (1-rho))))     
    }else
    {}
  
  
#### Temporal neighbourhood matrix
D <-array(0, c(N,N))
    for(i in 1:N)
    {
        for(j in 1:N)
        {
            if(abs((i-j))==1)  D[i,j] <- 1 
        }    
    }


#### Temporal triplet object
D.triplet <- c(NA, NA, NA)
    for(i in 1:N)
    {
        for(j in 1:N)
        {
            if(D[i,j]>0)
            {
            D.triplet <- rbind(D.triplet, c(i,j, D[i,j]))     
            }else{}
        }
    }
D.triplet <- D.triplet[-1, ]     
D.n.triplet <- nrow(D.triplet) 
D.triplet.sum <- tapply(D.triplet[ ,3], D.triplet[ ,1], sum)
D.neighbours <- tapply(D.triplet[ ,3], D.triplet[ ,1], length)


#### Temporal begfin argument
D.begfin <- array(NA, c(N, 2))     
temp <- 1
    for(i in 1:N)
    {
    D.begfin[i, ] <- c(temp, (temp + D.neighbours[i]-1))
    temp <- temp + D.neighbours[i]
    }


#### Temporal determinant     
    if(!fix.rho.T) 
    {
    Dstar <- diag(apply(D,1,sum)) - D
    Dstar.eigen <- eigen(Dstar)
    Dstar.val <- Dstar.eigen$values
    det.Q.D <-  0.5 * sum(log((lambda * Dstar.val + (1-lambda))))    
    }else
    {}



##########################################
#### Specify quantities that do not change
##########################################
#### Check for islands
W.list<- mat2listw(W)
W.nb <- W.list$neighbours
W.islands <- n.comp.nb(W.nb)
islands <- W.islands$comp.id
n.islands <- max(W.islands$nc)
    if(rho==1) tau2.phi.shape <- prior.tau2[1] + 0.5 * (K-n.islands)   
    if(lambda==1) tau2.delta.shape <- prior.tau2[1] + 0.5 * (N-1)   



###########################
#### Run the Bayesian model
###########################
#### Start timer
     if(verbose)
     {
     cat("Generating", n.keep, "post burnin and thinned (if requested) samples.\n", sep = " ")
     progressBar <- txtProgressBar(style = 3)
     percentage.points<-round((1:100/100)*n.sample)
     }else
     {
     percentage.points<-round((1:100/100)*n.sample)     
     }
     
  
#### Create the MCMC samples
    for(j in 1:n.sample)
    {
    ####################################
    ## Sample from Y - data augmentation
    ####################################
        if(n.miss>0)
        {
        Y.DA[which.miss==0] <- rpois(n=n.miss, lambda=fitted[which.miss==0])    
        }else
        {}
    Y.DA.mat <- matrix(Y.DA, nrow=K, ncol=N, byrow=FALSE)
        

    
    ###################
    ## Sample from beta
    ###################
    offset.temp <- offset + as.numeric(phi.mat) + as.numeric(delta.mat) + as.numeric(gamma.mat)   
        if(MALA)
        {
        temp <- poissonbetaupdateMALA(X.standardised, N.all, p, beta, offset.temp, Y.DA, prior.mean.beta, prior.var.beta, n.beta.block, proposal.sd.beta, list.block)
        }else
        {
        temp <- poissonbetaupdateRW(X.standardised, N.all, p, beta, offset.temp, Y.DA, prior.mean.beta, prior.var.beta, n.beta.block, proposal.sd.beta, list.block)
        }
    beta <- temp[[1]]
    accept[1] <- accept[1] + temp[[2]]
    accept[2] <- accept[2] + n.beta.block  
    regression.mat <- matrix(X.standardised %*% beta, nrow=K, ncol=N, byrow=FALSE)  
 
      
      
    ####################
    ## Sample from phi
    ####################
    phi.offset <- offset.mat + regression.mat + delta.mat + gamma.mat
        if(MALA)
        {
        temp1 <- poissoncarupdateMALA(W.triplet, W.begfin, W.triplet.sum, K, phi, tau2.phi, Y.DA.mat, proposal.sd.phi, rho, phi.offset, N, rep(1,N))
        }else
        {
        temp1 <- poissoncarupdateRW(W.triplet, W.begfin, W.triplet.sum, K, phi, tau2.phi, Y.DA.mat, proposal.sd.phi, rho, phi.offset, N, rep(1,N))
        }
    phi <- temp1[[1]]
        if(rho<1)
        {
        phi <- phi - mean(phi)
        }else
        {
        phi[which(islands==1)] <- phi[which(islands==1)] - mean(phi[which(islands==1)])   
        }
    phi.mat <- matrix(rep(phi, N), byrow=F, nrow=K)    
    accept[3] <- accept[3] + temp1[[2]]
    accept[4] <- accept[4] + K  
    
    
    
    ####################
    ## Sample from delta
    ####################
    delta.offset <- t(offset.mat + regression.mat + phi.mat + gamma.mat)
        if(MALA)
        {
        temp2 <- poissoncarupdateMALA(D.triplet, D.begfin, D.triplet.sum, N, delta, tau2.delta, t(Y.DA.mat), proposal.sd.delta, lambda, delta.offset, K, rep(1,K))
        }else
        {
        temp2 <- poissoncarupdateRW(D.triplet, D.begfin, D.triplet.sum, N, delta, tau2.delta, t(Y.DA.mat), proposal.sd.delta, lambda, delta.offset, K, rep(1,K))
        }
    delta <- temp2[[1]]
    delta <- delta - mean(delta)
    delta.mat <- matrix(rep(delta, K), byrow=T, nrow=K)
    accept[5] <- accept[5] + temp2[[2]]
    accept[6] <- accept[6] + N      
    
    

        if(interaction)
        {
        ####################
        ## Sample from gamma
        ####################
        gamma.offset <- offset.mat + regression.mat + phi.mat +  delta.mat
        gamma.offset.vec <- as.numeric(gamma.offset)
            if(MALA)
            {
            temp5 <- poissonindepupdateMALA(N.all, gamma, tau2.gamma, Y.DA, proposal.sd.gamma, gamma.offset.vec)
            }else
            {
            temp5 <- poissonindepupdateRW(N.all, gamma, tau2.gamma, Y.DA, proposal.sd.gamma, gamma.offset.vec)
            }
        gamma <- temp5[[1]]
        gamma <- gamma - mean(gamma)
        gamma.mat <- matrix(gamma, byrow=F, nrow=K)
        accept[7] <- accept[7] + temp5[[2]]
        accept[8] <- accept[8] + N * K
    
    
        
        #########################
        ## Sample from tau2.gamma
        #########################
        tau2.gamma.scale <- prior.tau2[2]  + sum(gamma.mat^2)/2
        tau2.gamma <- 1 / rgamma(1, tau2.gamma.shape, scale=(1/tau2.gamma.scale)) 
        }else
        {}
        
    
    
    #######################
    ## Sample from tau2.phi
    #######################
    temp2.phi <- quadform(W.triplet, W.triplet.sum, W.n.triplet, K, phi, phi, rho)
    tau2.phi.scale <- temp2.phi + prior.tau2[2] 
    tau2.phi <- 1 / rgamma(1, tau2.phi.shape, scale=(1/tau2.phi.scale))
    
    
    
    #########################
    ## Sample from tau2.delta
    #########################
    temp2.delta <- quadform(D.triplet, D.triplet.sum, D.n.triplet, N, delta, delta, lambda)
    tau2.delta.scale <- temp2.delta + prior.tau2[2] 
    tau2.delta <- 1 / rgamma(1, tau2.delta.shape, scale=(1/tau2.delta.scale))

    

    ##################
    ## Sample from rho
    ##################
        if(!fix.rho.S)
        {
        proposal.rho <- rtruncnorm(n=1, a=0, b=1, mean=rho, sd=proposal.sd.rho)   
        temp3 <- quadform(W.triplet, W.triplet.sum, W.n.triplet, K, phi, phi, proposal.rho)
        det.Q.proposal <- 0.5 * sum(log((proposal.rho * Wstar.val + (1-proposal.rho))))              
        logprob.current <- det.Q.W - temp2.phi / tau2.phi
        logprob.proposal <- det.Q.proposal - temp3 / tau2.phi
        hastings <- log(dtruncnorm(x=rho, a=0, b=1, mean=proposal.rho, sd=proposal.sd.rho)) - log(dtruncnorm(x=proposal.rho, a=0, b=1, mean=rho, sd=proposal.sd.rho)) 
        prob <- exp(logprob.proposal - logprob.current + hastings)
        
        #### Accept or reject the proposal
            if(prob > runif(1))
            {
            rho <- proposal.rho
            det.Q.W <- det.Q.proposal
            accept[9] <- accept[9] + 1           
            }else
            {
            }              
        accept[10] <- accept[10] + 1           
        }else
        {}


    
    #####################
    ## Sample from lambda
    #####################
        if(!fix.rho.T)
        {
        proposal.lambda <- rtruncnorm(n=1, a=0, b=1, mean=lambda, sd=proposal.sd.lambda)   
        temp3 <- quadform(D.triplet, D.triplet.sum, D.n.triplet, N, delta, delta, proposal.lambda)
        det.Q.proposal <- 0.5 * sum(log((proposal.lambda * Dstar.val + (1-proposal.lambda))))              
        logprob.current <- det.Q.D - temp2.delta / tau2.delta
        logprob.proposal <- det.Q.proposal - temp3 / tau2.delta
        hastings <- log(dtruncnorm(x=lambda, a=0, b=1, mean=proposal.lambda, sd=proposal.sd.lambda)) - log(dtruncnorm(x=proposal.lambda, a=0, b=1, mean=lambda, sd=proposal.sd.lambda)) 
        prob <- exp(logprob.proposal - logprob.current + hastings)

        #### Accept or reject the proposal
            if(prob > runif(1))
            {
            lambda <- proposal.lambda
            det.Q.D <- det.Q.proposal
            accept[11] <- accept[11] + 1           
            }else
            {
            }              
        accept[12] <- accept[12] + 1           
        }else
        {}


    
    #########################
    ## Calculate the deviance
    #########################
    lp <- as.numeric(offset.mat + regression.mat + phi.mat  + delta.mat + gamma.mat)
    fitted <- exp(lp)
    loglike <- dpois(x=as.numeric(Y), lambda=fitted, log=TRUE)

    
    
    ###################
    ## Save the results
    ###################
        if(j > burnin & (j-burnin)%%thin==0)
        {
        ele <- (j - burnin) / thin
        samples.beta[ele, ] <- beta
        samples.phi[ele, ] <- phi
        samples.delta[ele, ] <- delta
            if(!fix.rho.S) samples.rho[ele, ] <- rho
            if(!fix.rho.T) samples.lambda[ele, ] <- lambda
        samples.fitted[ele, ] <- fitted
        samples.loglike[ele, ] <- loglike
            if(n.miss>0) samples.Y[ele, ] <- Y.DA[which.miss==0]
            if(interaction)
            {
            samples.gamma[ele, ] <- gamma
            samples.tau2[ele, ] <- c(tau2.phi, tau2.delta, tau2.gamma)        
            }else
            {
            samples.tau2[ele, ] <- c(tau2.phi, tau2.delta)
            }
        }else
        {}
    
    
    
    ########################################
    ## Self tune the acceptance probabilties
    ########################################
    k <- j/100
        if(ceiling(k)==floor(k))
        {
        #### Update the proposal sds
            if(p>2)
            {
            proposal.sd.beta <- common.accceptrates1(accept[1:2], proposal.sd.beta, 40, 50)
            }else
            {
            proposal.sd.beta <- common.accceptrates1(accept[1:2], proposal.sd.beta, 30, 40)    
            }
        proposal.sd.phi <- common.accceptrates1(accept[3:4], proposal.sd.phi, 40, 50)
        proposal.sd.delta <- common.accceptrates1(accept[5:6], proposal.sd.delta, 40, 50)
            if(interaction) proposal.sd.gamma <- common.accceptrates1(accept[7:8], proposal.sd.gamma, 40, 50)
            if(!fix.rho.S) proposal.sd.rho <- common.accceptrates2(accept[9:10], proposal.sd.rho, 40, 50, 0.5)
            if(!fix.rho.T) proposal.sd.lambda <- common.accceptrates2(accept[11:12], proposal.sd.lambda, 40, 50, 0.5)
        accept.all <- accept.all + accept
        accept <- rep(0,12)
        }else
        {}

    
        
    ################################       
    ## print progress to the console
    ################################
        if(j %in% percentage.points & verbose)
        {
         setTxtProgressBar(progressBar, j/n.sample)
         }
    }

    
#### end timer
     if(verbose)
     {
     cat("\nSummarising results.")
     close(progressBar)
     }else
     {}
  
 

###################################
#### Summarise and save the results 
###################################
#### Compute the acceptance rates
accept.beta <- 100 * accept.all[1] / accept.all[2]
accept.phi <- 100 * accept.all[3] / accept.all[4]
accept.delta <- 100 * accept.all[5] / accept.all[6]
    if(!fix.rho.S)
    {
    accept.rho <- 100 * accept.all[9] / accept.all[10]
    }else
    {
    accept.rho <- NA    
    }

    if(!fix.rho.T)
    {
    accept.lambda <- 100 * accept.all[11] / accept.all[12]
    }else
    {
    accept.lambda <- NA    
    }

    if(interaction)
    {
    accept.gamma <- 100 * accept.all[7] / accept.all[8]
    accept.final <- c(accept.beta, accept.phi, accept.delta, accept.gamma, accept.rho, accept.lambda)
    names(accept.final) <- c("beta", "phi", "delta", "gamma", "rho.S", "rho.T")
    }else
    {
    accept.final <- c(accept.beta, accept.phi, accept.delta, accept.rho, accept.lambda)
    names(accept.final) <- c("beta", "phi", "delta", "rho.S", "rho.T")        
    }


#### Compute the fitted deviance
mean.phi <- apply(samples.phi, 2, mean)
mean.delta <- apply(samples.delta, 2, mean)  
mean.phi.mat <- matrix(rep(mean.phi, N), byrow=F, nrow=K)
mean.delta.mat <- matrix(rep(mean.delta, K), byrow=T, nrow=K)
mean.beta <- apply(samples.beta,2,mean)
regression.mat <- matrix(X.standardised %*% mean.beta, nrow=K, ncol=N, byrow=FALSE)   
    if(interaction)
    {
    mean.gamma <- apply(samples.gamma, 2,mean)
    mean.gamma.mat <-  matrix(mean.gamma, byrow=F, nrow=K)   
    fitted.mean <- as.numeric(exp(offset.mat + regression.mat + mean.phi.mat + mean.delta.mat + mean.gamma.mat))
    }else
    {
    fitted.mean <- as.numeric(exp(offset.mat + regression.mat + mean.phi.mat + mean.delta.mat))    
    }

deviance.fitted <- -2 * sum(dpois(x=as.numeric(Y), lambda=fitted.mean, log=TRUE), na.rm=TRUE)


#### Model fit criteria
modelfit <- common.modelfit(samples.loglike, deviance.fitted)


#### Create the fitted values and residuals
fitted.values <- apply(samples.fitted, 2, mean)
response.residuals <- as.numeric(Y) - fitted.values
pearson.residuals <- response.residuals /sqrt(fitted.values)
residuals <- data.frame(response=response.residuals, pearson=pearson.residuals)
 
 
#### Transform the parameters back to the origianl covariate scale.
samples.beta.orig <- common.betatransform(samples.beta, X.indicator, X.mean, X.sd, p, FALSE)


#### Create a summary object
samples.beta.orig <- mcmc(samples.beta.orig)
summary.beta <- t(apply(samples.beta.orig, 2, quantile, c(0.5, 0.025, 0.975))) 
summary.beta <- cbind(summary.beta, rep(n.keep, p), rep(accept.beta,p), effectiveSize(samples.beta.orig), geweke.diag(samples.beta.orig)$z)
rownames(summary.beta) <- colnames(X)
colnames(summary.beta) <- c("Median", "2.5%", "97.5%", "n.sample", "% accept", "n.effective", "Geweke.diag")


    if(interaction)
    {
    summary.hyper <- array(NA, c(5, 7))     
    summary.hyper[1,1:3] <- quantile(samples.tau2[ ,1], c(0.5, 0.025, 0.975))
    summary.hyper[2,1:3] <- quantile(samples.tau2[ ,2], c(0.5, 0.025, 0.975))
    summary.hyper[3,1:3] <- quantile(samples.tau2[ ,3], c(0.5, 0.025, 0.975))
    rownames(summary.hyper) <- c("tau2.S", "tau2.T", "tau2.I", "rho.S", "rho.T")     
    summary.hyper[1, 4:7] <- c(n.keep, 100, effectiveSize(mcmc(samples.tau2[ ,1])), geweke.diag(mcmc(samples.tau2[ ,1]))$z)     
    summary.hyper[2, 4:7] <- c(n.keep, 100, effectiveSize(mcmc(samples.tau2[ ,2])), geweke.diag(mcmc(samples.tau2[ ,2]))$z)   
    summary.hyper[3, 4:7] <- c(n.keep, 100, effectiveSize(mcmc(samples.tau2[ ,3])), geweke.diag(mcmc(samples.tau2[ ,3]))$z)   
    
        if(!fix.rho.S)
        {
        summary.hyper[4, 1:3] <- quantile(samples.rho, c(0.5, 0.025, 0.975))
        summary.hyper[4, 4:7] <- c(n.keep, accept.rho, effectiveSize(samples.rho), geweke.diag(samples.rho)$z)
        }else
        {
        summary.hyper[4, 1:3] <- c(rho, rho, rho)
        summary.hyper[4, 4:7] <- rep(NA, 4)
        }
    
        if(!fix.rho.T)
        {
        summary.hyper[5, 1:3] <- quantile(samples.lambda, c(0.5, 0.025, 0.975))
        summary.hyper[5, 4:7] <- c(n.keep, accept.lambda, effectiveSize(samples.lambda), geweke.diag(samples.lambda)$z)
        }else
        {
        summary.hyper[5, 1:3] <- c(lambda, lambda, lambda)
        summary.hyper[5, 4:7] <- rep(NA, 4)
        }   
    }else
    {
    summary.hyper <- array(NA, c(4, 7))     
    summary.hyper[1,1:3] <- quantile(samples.tau2[ ,1], c(0.5, 0.025, 0.975))
    summary.hyper[2,1:3] <- quantile(samples.tau2[ ,2], c(0.5, 0.025, 0.975))
    rownames(summary.hyper) <- c("tau2.S", "tau2.T",  "rho.S", "rho.T")     
    summary.hyper[1, 4:7] <- c(n.keep, 100, effectiveSize(mcmc(samples.tau2[ ,1])), geweke.diag(mcmc(samples.tau2[ ,1]))$z)     
    summary.hyper[2, 4:7] <- c(n.keep, 100, effectiveSize(mcmc(samples.tau2[ ,2])), geweke.diag(mcmc(samples.tau2[ ,2]))$z)   
    
        if(!fix.rho.S)
        {
        summary.hyper[3, 1:3] <- quantile(samples.rho, c(0.5, 0.025, 0.975))
        summary.hyper[3, 4:7] <- c(n.keep, accept.rho, effectiveSize(samples.rho), geweke.diag(samples.rho)$z)
        }else
        {
        summary.hyper[3, 1:3] <- c(rho, rho, rho)
        summary.hyper[3, 4:7] <- rep(NA, 4)
        }
    
        if(!fix.rho.T)
        {
        summary.hyper[4, 1:3] <- quantile(samples.lambda, c(0.5, 0.025, 0.975))
        summary.hyper[4, 4:7] <- c(n.keep, accept.lambda, effectiveSize(samples.lambda), geweke.diag(samples.lambda)$z)
        }else
        {
        summary.hyper[4, 1:3] <- c(lambda, lambda, lambda)
        summary.hyper[4, 4:7] <- rep(NA, 4)
        }   
    }
summary.results <- rbind(summary.beta, summary.hyper)
summary.results[ , 1:3] <- round(summary.results[ , 1:3], 4)
summary.results[ , 4:7] <- round(summary.results[ , 4:7], 1)


#### Compile and return the results
#### Harmonise samples in case of them not being generated
    if(fix.rho.S & fix.rho.T)
    {
    samples.rhoext <- NA
    }else if(fix.rho.S & !fix.rho.T)
    {
    samples.rhoext <- samples.lambda
    names(samples.rhoext) <- "rho.T"
    }else if(!fix.rho.S & fix.rho.T)
    {
    samples.rhoext <- samples.rho  
    names(samples.rhoext) <- "rho.S"
    }else
    {
    samples.rhoext <- cbind(samples.rho, samples.lambda)
    colnames(samples.rhoext) <- c("rho.S", "rho.T")
    }
    if(n.miss==0) samples.Y = NA
    if(!interaction) samples.gamma = NA


    if(interaction)
    {
    model.string <- c("Likelihood model - Poisson (log link function)", "\nLatent structure model - spatial and temporal main effects and an interaction\n")
    colnames(samples.tau2) <- c("tau2.S", "tau2.T", "tau2.I")
    }else
    {
    model.string <- c("Likelihood model - Poisson (log link function)", "\nLatent structure model - spatial and temporal main effects\n")
    colnames(samples.tau2) <- c("tau2.S", "tau2.T")
    }

samples <- list(beta=mcmc(samples.beta.orig), phi=mcmc(samples.phi),  delta=mcmc(samples.delta), gamma=mcmc(samples.gamma), tau2=mcmc(samples.tau2), rho=mcmc(samples.rhoext), fitted=mcmc(samples.fitted), Y=mcmc(samples.Y))        
results <- list(summary.results=summary.results, samples=samples, fitted.values=fitted.values, residuals=residuals, modelfit=modelfit, accept=accept.final, localised.structure=NULL, formula=formula, model=model.string,  X=X)
class(results) <- "CARBayesST"


#### Finish by stating the time taken 
     if(verbose)
     {
     b<-proc.time()
     cat("Finished in ", round(b[3]-a[3], 1), "seconds.\n")
     }else
     {}
  return(results)
}
