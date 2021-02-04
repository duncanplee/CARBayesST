binomial.CARsepspatial <- function(formula, data=NULL, trials, W, burnin, n.sample, thin=1,  prior.mean.beta=NULL, prior.var.beta=NULL, prior.tau2=NULL, rho.S=NULL, rho.T=NULL, MALA=FALSE, verbose=TRUE)
{
  ##############################################
  #### Format the arguments and check for errors
  ##############################################
  #### Verbose
  a <- common.verbose(verbose)  
  
  #### Frame object
  frame.results <- common.frame(formula, data, "binomial")
  N.all <- frame.results$n
  p <- frame.results$p
  X <- frame.results$X
  X.standardised <- frame.results$X.standardised
  X.sd <- frame.results$X.sd
  X.mean <- frame.results$X.mean
  X.indicator <- frame.results$X.indicator 
  offset <- frame.results$offset
  Y <- frame.results$Y
  failures <- trials - Y
  n.miss <- frame.results$n.miss  
    if(n.miss>0) stop("the response has missing 'NA' values.", call.=FALSE)
  #### Check on MALA argument
    if(length(MALA)!=1) stop("MALA is not length 1.", call.=FALSE)
    if(!is.logical(MALA)) stop("MALA is not logical.", call.=FALSE) 
  
  ## Check for errors
  if(sum(is.na(trials))>0) stop("the numbers of trials has missing 'NA' values.", call.=FALSE)
  if(!is.numeric(trials)) stop("the numbers of trials has non-numeric values.", call.=FALSE)
  if(min(trials)<=0) stop("the numbers of trials has zero or negative values.", call.=FALSE)
  if(sum(Y>trials)>0) stop("the response variable has larger values that the numbers of trials.", call.=FALSE)
  
  
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
  

  
  #### Specify the initial parameter values
  dat <- cbind(Y, failures)
  mod.glm <- glm(dat~X.standardised-1, offset=offset, family="quasibinomial")
  beta.mean <- mod.glm$coefficients
  beta.sd <- sqrt(diag(summary(mod.glm)$cov.scaled))
  beta <- rnorm(n=length(beta.mean), mean=beta.mean, sd=beta.sd)
  
  theta.hat <- Y / trials
  theta.hat[theta.hat==0] <- 0.01
  theta.hat[theta.hat==1] <- 0.99
  res.temp <- log(theta.hat / (1 - theta.hat)) - X.standardised %*% beta - offset
  res.sd <- sd(res.temp, na.rm=TRUE)/5
  phi <- rnorm(n=N.all, mean=0, sd = res.sd)
  phi.mat <- matrix(phi, nrow=K, ncol=N, byrow=FALSE)  
  delta <- rnorm(n=N, mean=0, sd = res.sd)
  tau2 <- apply(phi.mat, 2, var) / 10
  sig2 <- var(delta)/10

  
   #### Check and specify the priors
  ## Put in default priors
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
  
  #### Set up matrices to store samples
  n.keep <- floor((n.sample - burnin)/thin)
  samples.beta <- array(NA, c(n.keep, p))
  samples.phi <- array(NA, c(n.keep, N.all))
  samples.tau2 <- array(NA, c(n.keep, N))
  samples.sig2 <- array(NA, c(n.keep, 1))
  if(!fix.rho.S) samples.rho <- array(NA, c(n.keep, 1))
  if(!fix.rho.T) samples.lambda <- array(NA, c(n.keep, 1))
  samples.delta <- array(NA, c(n.keep, N))     
  samples.fitted <- array(NA, c(n.keep, N.all))
  samples.loglike <- array(NA, c(n.keep, N.all))
  
  #### Specify the Metropolis quantities
  accept <- rep(0,10)
  proposal.sd.phi <- 0.1
  proposal.sd.rho <- 0.05
  proposal.sd.beta <- 0.01
  proposal.sd.delta <- 0.05
  proposal.sd.lambda <- 0.02
  proposal.corr.beta <- solve(t(X.standardised) %*% X.standardised)
  chol.proposal.corr.beta <- chol(proposal.corr.beta)     
  tau2.shape <- prior.tau2[1] + K/2
  sig2.shape <- prior.tau2[1] + N/2
  
  #### Spatial quantities
  ## Create the determinant     
  if(!fix.rho.S) 
  {
    Wstar <- diag(apply(W,1,sum)) - W
    Wstar.eigen <- eigen(Wstar)
    Wstar.val <- Wstar.eigen$values
    det.Q.W <-  0.5 * sum(log((rho * Wstar.val + (1-rho))))     
  }else
  {}
  
  #### .T quantities
  ## .T neighbourhood matrix
  D <-array(0, c(N,N))
  for(i in 1:N)
  {
    for(j in 1:N)
    {
      if(abs((i-j))==1)  D[i,j] <- 1 
    }    
  }
  
  ## Create the triplet object
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
  
  ## Create the start and finish points for W updating
  D.begfin <- array(NA, c(N, 2))     
  temp <- 1
  for(i in 1:N)
  {
    D.begfin[i, ] <- c(temp, (temp + D.neighbours[i]-1))
    temp <- temp + D.neighbours[i]
  }
  
  ## Create the determinant     
  if(!fix.rho.T) 
  {
    Dstar <- diag(apply(D,1,sum)) - D
    Dstar.eigen <- eigen(Dstar)
    Dstar.val <- Dstar.eigen$values
    det.Q.D <-  0.5 * sum(log((lambda * Dstar.val + (1-lambda))))    
  }else
  {} 
  
  #### Specify quantities that do not change
  offset.mat <- matrix(offset, nrow=K, ncol=N, byrow=FALSE) 
  regression.mat <- matrix(X.standardised %*% beta, nrow=K, ncol=N, byrow=FALSE)
  Y.mat <- matrix(Y, nrow=K, ncol=N, byrow=FALSE)
  trials.mat <- matrix(trials, nrow=K, ncol=N, byrow=FALSE)
  failures.mat <- matrix(failures, nrow=K, ncol=N, byrow=FALSE)
  delta.mat <- matrix(delta, nrow=K, ncol=N, byrow=TRUE)
  
  ## Check for islands
  W.list<- mat2listw(W)
  W.nb <- W.list$neighbours
  W.islands <- n.comp.nb(W.nb)
  islands <- W.islands$comp.id
  n.islands <- max(W.islands$nc)
  n.island1 <- length(which(islands==1))
  if(rho==1) tau2.shape <- prior.tau2[1] + 0.5 * (K-n.islands)   
  if(lambda==1) sig2.shape <- prior.tau2[1] + 0.5 * (N-1)   
  
  ###########################
  #### Run the Bayesian model
  ###########################
  ## Start timer
  if(verbose)
  {
    cat("Generating", n.keep, "post burnin and thinned (if requested) samples.\n", sep = " ")
    progressBar <- txtProgressBar(style = 3)
    percentage.points<-round((1:100/100)*n.sample)
  }else
  {
    percentage.points<-round((1:100/100)*n.sample)     
  }
  
  for(j in 1:n.sample)
  {
    ###################
    ## Sample from beta
    ###################
    offset.temp <- as.numeric(offset.mat + phi.mat + delta.mat)    
    if(MALA)
    {
      temp <- binomialbetaupdateMALA(X.standardised, N.all, p, beta, offset.temp, Y, failures, trials, prior.mean.beta, prior.var.beta, n.beta.block, proposal.sd.beta, list.block)
    }else
    {
      temp <- binomialbetaupdateRW(X.standardised, N.all, p, beta, offset.temp, Y, failures, prior.mean.beta, prior.var.beta, n.beta.block, proposal.sd.beta, list.block)
    }
    beta <- temp[[1]]
    accept[1] <- accept[1] + temp[[2]]
    accept[2] <- accept[2] + n.beta.block  
    regression.mat <- matrix(X.standardised %*% beta, nrow=K, ncol=N, byrow=FALSE)  

    ####################
    ## Sample from phi
    ####################
    phi.offset <- offset.mat + regression.mat + delta.mat
    den.offset <- rho * W.triplet.sum + 1 - rho
    if(MALA==TRUE)
    {
      temp1 <- binomialsrecarupdateMALA(W.triplet, W.begfin, W.triplet.sum, K, N, phi.mat, rho, Y.mat, failures.mat, trials.mat, proposal.sd.phi, phi.offset, den.offset, tau2)
    }else
    {
      temp1 <- binomialsrecarupdateRW(W.triplet, W.begfin, W.triplet.sum, K, N, phi.mat, rho, Y.mat, failures.mat, proposal.sd.phi, phi.offset, den.offset, tau2)
    }
    phi.temp <- temp1[[1]]
    phi.mean <- apply(phi.temp,2,mean)
    if(rho<1)
    {
      phi <- as.numeric(phi.temp) - kronecker(phi.mean, rep(1,K))
    }else
    {
      phi.temp[which(islands==1), ] <- phi.temp[which(islands==1), ] - matrix(kronecker(phi.mean, rep(1,n.island1)), ncol=N, byrow=F) 
      phi <- as.numeric(phi.temp)
    }
    phi.mat <- matrix(phi, nrow=K, ncol=N, byrow=FALSE)
    
    accept[3] <- accept[3] + temp1[[2]]
    accept[4] <- accept[4] + N.all
    
    #####################
    ## Samples from delta
    #####################
    delta.offset <- t(offset.mat + phi.mat + regression.mat)
    if(MALA==TRUE)
    {
      temp2 <- binomialcarupdateMALA(D.triplet, D.begfin, D.triplet.sum, N, delta, sig2, t(Y.mat), t(failures.mat), t(trials.mat), proposal.sd.delta, lambda, delta.offset, K, rep(1,K))
    }else
    {
      temp2 <- binomialcarupdateRW(D.triplet, D.begfin, D.triplet.sum, N, delta, sig2, t(Y.mat), t(failures.mat), proposal.sd.delta, lambda, delta.offset, K, rep(1,K))
    }
    delta <- temp2[[1]]
    delta <- delta - mean(delta)
    delta.mat <- matrix(rep(delta, K), byrow=T, nrow=K)
    accept[7] <- accept[7] + temp2[[2]]
    accept[8] <- accept[8] + N      

    ####################
    ## Samples from tau2
    ####################
    tau2.temp <- tauquadformcompute2(W.triplet, W.triplet.sum, W.n.triplet, K, N, phi.mat, rho)
    tau2 <- tau2compute(tau2, tau2.temp, tau2.shape, prior.tau2[2], N)
    
    ####################
    ## Samples from sig2
    ####################
    temp2.delta <- quadform(D.triplet, D.triplet.sum, D.n.triplet, N, delta, delta, lambda)
    sig2.scale <- temp2.delta + prior.tau2[2] 
    sig2 <- 1 / rgamma(1, sig2.shape, scale=(1/sig2.scale))
    
    ##################
    ## Sample from rho
    ##################
    if(!fix.rho.S)
    {
      temp3 <- rhoquadformcompute(W.triplet, W.triplet.sum, W.n.triplet, K, N, phi.mat, rho, tau2)
      proposal.rho <- rtruncnorm(n=1, a=0, b=1, mean=rho, sd=proposal.sd.rho)
      temp4 <- rhoquadformcompute(W.triplet, W.triplet.sum, W.n.triplet,  K, N, phi.mat, proposal.rho, tau2)
      det.Q.W.proposal <- 0.5 * sum(log((proposal.rho * Wstar.val + (1-proposal.rho))))
      logprob.current <- N * det.Q.W - temp3
      logprob.proposal <- N * det.Q.W.proposal - temp4
      hastings <- log(dtruncnorm(x=rho, a=0, b=1, mean=proposal.rho, sd=proposal.sd.rho)) - log(dtruncnorm(x=proposal.rho, a=0, b=1, mean=rho, sd=proposal.sd.rho)) 
      prob <- exp(logprob.proposal - logprob.current + hastings)
      if(prob > runif(1))
      {
        rho <- proposal.rho
        det.Q.W <- det.Q.W.proposal
        accept[5] <- accept[5] + 1
      }else
      {
      }
      accept[6] <- accept[6] + 1
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
      logprob.current <- det.Q.D - temp2.delta / sig2
      logprob.proposal <- det.Q.proposal - temp3 / sig2
      hastings <- log(dtruncnorm(x=lambda, a=0, b=1, mean=proposal.lambda, sd=proposal.sd.lambda)) - log(dtruncnorm(x=proposal.lambda, a=0, b=1, mean=lambda, sd=proposal.sd.lambda)) 
      prob <- exp(logprob.proposal - logprob.current + hastings)
      
      #### Accept or reject the proposal
      if(prob > runif(1))
      {
        lambda <- proposal.lambda
        det.Q.D <- det.Q.proposal
        accept[9] <- accept[9] + 1           
      }else
      {
      }              
      accept[10] <- accept[10] + 1           
    }else
    {}
    
    #########################
    ## Calculate the deviance
    #########################
    lp <- as.numeric(offset.mat + regression.mat + phi.mat + delta.mat)
    prob <- exp(lp) / (1+exp(lp))
    fitted <- trials * prob
    loglike <- dbinom(x=Y, size=trials, prob=prob, log=TRUE)


    
    ###################
    ## Save the results
    ###################
    if(j > burnin & (j-burnin)%%thin==0)
    {
      ele <- (j - burnin) / thin
      samples.beta[ele, ] <- beta
      samples.phi[ele, ] <- as.numeric(phi)
      if(!fix.rho.S) samples.rho[ele, ] <- rho
      if(!fix.rho.T) samples.lambda[ele, ] <- lambda
      samples.tau2[ele, ] <- tau2
      samples.sig2[ele, ] <- sig2
      samples.delta[ele, ] <- delta
      samples.fitted[ele, ] <- fitted
      samples.loglike[ele, ] <- loglike
    }else
    {
    }
    
    ########################################
    ## Self tune the acceptance probabilties
    ########################################
    if(ceiling(j/100)==floor(j/100) & j < burnin)
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
      proposal.sd.delta <- common.accceptrates1(accept[7:8], proposal.sd.delta, 40, 50)
      if(!fix.rho.S) proposal.sd.rho <- common.accceptrates2(accept[5:6], proposal.sd.rho, 40, 50, 0.5)
      if(!fix.rho.T) proposal.sd.lambda <- common.accceptrates2(accept[9:10], proposal.sd.lambda, 40, 50, 0.5)
      accept <- rep(0,10)
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
  
  # end timer
  if(verbose)
  {
      cat("\nSummarising results.")
    close(progressBar)
  }else
  {}
  
  
###################################
#### Summarise and save the results 
###################################
## Compute the acceptance rates
accept.beta <- 100 * accept[1] / accept[2]
accept.phi <- 100 * accept[3] / accept[4]
accept.delta <- 100 * accept[7] / accept[8]
   if(!fix.rho.S)
   {
   accept.rho <- 100 * accept[5] / accept[6]
   }else
   {
   accept.rho <- NA    
   }
  
   if(!fix.rho.T)
   {
   accept.lambda <- 100 * accept[9] / accept[10]
   }else
   {
   accept.lambda <- NA    
   }
accept.final <- c(accept.beta, accept.phi, accept.delta, accept.rho, accept.lambda)
names(accept.final) <- c("beta", "phi", "delta", "rho.S", "rho.T")
  
  
#### Compute the fitted deviance
mean.beta <- apply(samples.beta,2,mean)
regression.mat <- matrix(X.standardised %*% mean.beta, nrow=K, ncol=N, byrow=FALSE)   
mean.phi <- matrix(apply(samples.phi, 2, mean), nrow=K, ncol=N)
mean.delta <- apply(samples.delta,2,mean)
delta.mat <- matrix(mean.delta, nrow=K, ncol=N, byrow=TRUE)
lp.mean <- as.numeric(offset.mat + mean.phi + regression.mat + delta.mat)   
mean.prob <- exp(lp.mean)  / (1 + exp(lp.mean))
deviance.fitted <- -2 * sum(dbinom(x=Y, size=trials, prob=mean.prob, log=TRUE), na.rm = T)

    
#### Model fit criteria
modelfit <- common.modelfit(samples.loglike, deviance.fitted)
  
  
#### Create the fitted values and residuals
fitted.values <- apply(samples.fitted, 2, mean)
response.residuals <- as.numeric(Y) - fitted.values
pearson.residuals <- response.residuals /sqrt(fitted.values * (1 - mean.prob))
residuals <- data.frame(response=response.residuals, pearson=pearson.residuals)

  
  #### Transform the parameters back to the origianl covariate scale.
  samples.beta.orig <- common.betatransform(samples.beta, X.indicator, X.mean, X.sd, p, FALSE)
  
  #### Create a summary object
  samples.beta.orig <- mcmc(samples.beta.orig)
  summary.beta <- t(apply(samples.beta.orig, 2, quantile, c(0.5, 0.025, 0.975))) 
  summary.beta <- cbind(summary.beta, rep(n.keep, p), rep(accept.beta,p), effectiveSize(samples.beta.orig), geweke.diag(samples.beta.orig)$z)
  rownames(summary.beta) <- colnames(X)
  colnames(summary.beta) <- c("Median", "2.5%", "97.5%", "n.sample", "% accept", "n.effective", "Geweke.diag")
  
  summary.hyper <- array(NA, c(3 + N, 7))    
  for (tt in  1:N) {
    summary.hyper[tt,1:3] <- quantile(samples.tau2[, tt], c(0.5, 0.025, 0.975))
    summary.hyper[tt, 4:7] <- c(n.keep, 100, effectiveSize(mcmc(samples.tau2[, tt])), geweke.diag(mcmc(samples.tau2[, tt]))$z) 
  }
  summary.hyper[N+1,1:3] <- quantile(samples.sig2, c(0.5, 0.025, 0.975))
  summary.hyper[N+1, 4:7] <- c(n.keep, 100, effectiveSize(mcmc(samples.sig2)), geweke.diag(mcmc(samples.sig2))$z)  
  
  if(!fix.rho.S)
  {
    summary.hyper[N+2, 1:3] <- quantile(samples.rho, c(0.5, 0.025, 0.975))
    summary.hyper[N+2, 4:7] <- c(n.keep, accept.rho, effectiveSize(samples.rho), geweke.diag(samples.rho)$z)
  }else
  {
    summary.hyper[N+2, 1:3] <- c(rho, rho, rho)
    summary.hyper[N+2, 4:7] <- rep(NA, 4)
  }
  
  if(!fix.rho.T)
  {
    summary.hyper[N+3, 1:3] <- quantile(samples.lambda, c(0.5, 0.025, 0.975))
    summary.hyper[N+3, 4:7] <- c(n.keep, accept.lambda, effectiveSize(samples.lambda), geweke.diag(samples.lambda)$z)
  }else
  {
    summary.hyper[N+3, 1:3] <- c(lambda, lambda, lambda)
    summary.hyper[N+3, 4:7] <- rep(NA, 4)
  }    
  
  rownames(summary.hyper) <- c(paste("tau2.", c(1:N), sep = ""), "tau2.T", "rho.S","rho.T")  
  summary.results <- rbind(summary.beta, summary.hyper)
  summary.results[ , 1:3] <- round(summary.results[ , 1:3], 4)
  summary.results[ , 4:7] <- round(summary.results[ , 4:7], 1)
  
  
#### Compile and return the results
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
  
  samples <- list(beta=mcmc(samples.beta.orig), phi=mcmc(samples.phi),  rho=mcmc(samples.rhoext), tau2=mcmc(samples.tau2), tau2.T=mcmc(samples.sig2),
                  delta=mcmc(samples.delta), fitted=mcmc(samples.fitted))
  model.string <- c("Likelihood model - binomial (logit link function)", "\nLatent structure model - An overall time trend with temporal specific spatial effects\n")
  results <- list(summary.results=summary.results, samples=samples, fitted.values=fitted.values, residuals=residuals, modelfit=modelfit, accept=accept.final, localised.structure=NULL, formula=formula, model=model.string,  X=X)
  class(results) <- "CARBayesST"
  if(verbose)
  {
    b<-proc.time()
    cat("Finished in ", round(b[3]-a[3], 1), "seconds.\n")
  }else
  {}
  return(results)
}