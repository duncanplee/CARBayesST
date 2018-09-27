binomial.CARlocalised <- function(formula, data=NULL, G, trials,  W, burnin, n.sample, thin=1,  prior.mean.beta=NULL, prior.var.beta=NULL, prior.delta=NULL, prior.tau2=NULL, MALA=TRUE, verbose=TRUE)
{
##############################################
#### Format the arguments and check for errors
##############################################
#### Verbose
a <- common.verbose(verbose)  
    
    
#### Frame object
frame.results <- common.frame.localised(formula, data, "binomial")
N.all <- frame.results$n
p <- frame.results$p
X <- frame.results$X
X.standardised <- frame.results$X.standardised
X.sd <- frame.results$X.sd
X.mean <- frame.results$X.mean
X.indicator <- frame.results$X.indicator 
offset <- frame.results$offset
Y <- frame.results$Y
which.miss <- as.numeric(!is.na(Y))
n.miss <- N.all - sum(which.miss)
    if(n.miss>0) stop("the response has missing 'NA' values.", call.=FALSE)


#### Check on MALA argument
    if(length(MALA)!=1) stop("MALA is not length 1.", call.=FALSE)
    if(!is.logical(MALA)) stop("MALA is not logical.", call.=FALSE)  


#### Check the trials argument
    if(sum(is.na(trials))>0) stop("the numbers of trials has missing 'NA' values.", call.=FALSE)
    if(!is.numeric(trials)) stop("the numbers of trials has non-numeric values.", call.=FALSE)
int.check <- N.all-sum(ceiling(trials)==floor(trials))
    if(int.check > 0) stop("the numbers of trials has non-integer values.", call.=FALSE)
    if(min(trials)<=0) stop("the numbers of trials has zero or negative values.", call.=FALSE)
    if(sum(Y>trials, na.rm=TRUE)>0) stop("the response variable has larger values that the numbers of trials.", call.=FALSE)
failures <- trials - Y


#### Compute a starting value for beta
    if(!is.null(X))
    {
    dat <- cbind(Y, failures)
    mod.glm <- glm(dat~X.standardised-1, offset=offset, family="quasibinomial")
    beta.mean <- mod.glm$coefficients
    beta.sd <- sqrt(diag(summary(mod.glm)$cov.scaled))
    beta <- rnorm(n=length(beta.mean), mean=beta.mean, sd=beta.sd)
    regression.vec <- X.standardised %*% beta
    }else
    {
    regression.vec <- rep(0, N.all)    
    }


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


#### Format and check the number of clusters G     
    if(length(G)!=1) stop("G is the wrong length.", call.=FALSE)    
    if(!is.numeric(G)) stop("G is not numeric.", call.=FALSE)    
    if(G<=1) stop("G is less than 2.", call.=FALSE)    
    if(G!=round(G)) stop("G is not an integer.", call.=FALSE) 
    if(floor(G/2)==ceiling(G/2))
    {
    Gstar <- G/2    
    }else
    {
    Gstar <- (G+1)/2          
    }


#### Priors
    if(!is.null(X))
    {
    if(is.null(prior.mean.beta)) prior.mean.beta <- rep(0, p)
    if(is.null(prior.var.beta)) prior.var.beta <- rep(100000, p)
    prior.beta.check(prior.mean.beta, prior.var.beta, p)
    }else
    {}
    if(is.null(prior.tau2)) prior.tau2 <- c(1, 0.01)
prior.var.check(prior.tau2)
    if(is.null(prior.delta)) prior.delta <- 10
    if(length(prior.delta)!=1) stop("the prior value for delta is the wrong length.", call.=FALSE)    
    if(!is.numeric(prior.delta)) stop("the prior value for delta is not numeric.", call.=FALSE)    
    if(sum(is.na(prior.delta))!=0) stop("the prior value for delta has missing values.", call.=FALSE)    
    if(prior.delta<=0) stop("the prior value for delta is not positive.", call.=FALSE)    


#### Compute the blocking structure for beta     
    if(!is.null(X))
    {
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
    }else
    {}


#### MCMC quantities - burnin, n.sample, thin
common.burnin.nsample.thin.check(burnin, n.sample, thin)  



#############################
#### Initial parameter values
#############################
theta.hat <- Y / trials
theta.hat[theta.hat==0] <- 0.01
theta.hat[theta.hat==1] <- 0.99
res.temp <- log(theta.hat / (1 - theta.hat)) - regression.vec - offset
clust <- kmeans(res.temp,G)
lambda <- clust$centers[order(clust$centers)]
lambda.mat <- matrix(rep(lambda, N), nrow=N, byrow=TRUE)
Z <- rep(1, N.all)
    for(j in 2:G)
    {
    Z[clust$cluster==order(clust$centers)[j]] <- j    
    }
Z.mat <- matrix(Z, nrow=K, ncol=N, byrow=FALSE)
mu <- matrix(lambda[Z], nrow=K, ncol=N, byrow=FALSE)
res.sd <- sd(res.temp, na.rm=TRUE)/5
phi.mat <- matrix(rnorm(n=N.all, mean=0, sd = res.sd), nrow=K, byrow=FALSE)
phi <- as.numeric(phi.mat)
tau2 <- var(phi)/10
gamma <- runif(1)
delta <- runif(1,1, min(2, prior.delta))


    
###############################    
#### Set up the MCMC quantities    
###############################
#### Matrices to store samples
n.keep <- floor((n.sample - burnin)/thin)
samples.Z <- array(NA, c(n.keep, N.all))
samples.lambda <- array(NA, c(n.keep, G))
samples.delta <- array(NA, c(n.keep, 1))
samples.tau2 <- array(NA, c(n.keep, 1))
samples.gamma <- array(NA, c(n.keep, 1))
samples.phi <- array(NA, c(n.keep, N.all))
samples.fitted <- array(NA, c(n.keep, N.all))
samples.loglike <- array(NA, c(n.keep, N.all))
    

#### Specify the Metropolis quantities
    if(!is.null(X))
    {
    samples.beta <- array(NA, c(n.keep, p))
    accept.all <- rep(0,8)
    proposal.corr.beta <- solve(t(X.standardised) %*% X.standardised)
    chol.proposal.corr.beta <- chol(proposal.corr.beta) 
    proposal.sd.beta <- 0.01
    }else
    {
    accept.all <- rep(0,6)    
    }

accept <- accept.all
proposal.sd.lambda <- 0.1
proposal.sd.delta <- 0.1
proposal.sd.phi <- 0.1
Y.extend <- matrix(rep(Y, G), byrow=F, ncol=G)
delta.update <- matrix(rep(1:G, N.all-K), ncol=G, byrow=T)
tau2.posterior.shape <- prior.tau2[1] + N * (K-1) /2
    
    
    
##########################################
#### Specify quantities that do not change
##########################################
which.miss.mat <- matrix(which.miss, nrow=K, ncol=N, byrow=FALSE)
Y.mat <- matrix(Y, nrow=K, ncol=N, byrow=FALSE) 
offset.mat <- matrix(offset, nrow=K, ncol=N, byrow=FALSE)    
regression.mat <- matrix(regression.vec, nrow=K, ncol=N, byrow=FALSE)
failures.mat <- matrix(failures, nrow=K, ncol=N, byrow=FALSE)  
trials.mat <- matrix(trials, nrow=K, ncol=N, byrow=FALSE)

    
    
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
    ####################
    ## Sample from beta
    ####################
        if(!is.null(X))
        {
        offset.temp <- offset + as.numeric(mu) + as.numeric(phi.mat)   
            if(p>2)
            {
            temp <- binomialbetaupdateMALA(X.standardised, N.all, p, beta, offset.temp, Y, failures, trials, prior.mean.beta, prior.var.beta, n.beta.block, proposal.sd.beta, list.block)
            }else
            {
            temp <- binomialbetaupdateRW(X.standardised, N.all, p, beta, offset.temp, Y, failures, prior.mean.beta, prior.var.beta, proposal.sd.beta)
            }
        beta <- temp[[1]]
        accept[7] <- accept[7] + temp[[2]]
        accept[8] <- accept[8] + n.beta.block  
        regression.vec <- X.standardised %*% beta
        regression.mat <- matrix(regression.vec, nrow=K, ncol=N, byrow=FALSE)    
        }else
        {}
        
        
        
    #######################     
    #### Sample from lambda
    #######################
    #### Propose a new value
    proposal.extend <- c(-100, lambda, 100) 
        for(r in 1:G)
        {
        proposal.extend[(r+1)] <- rtruncnorm(n=1, a=proposal.extend[r], b=proposal.extend[(r+2)], mean=proposal.extend[(r+1)], sd=proposal.sd.lambda)
        }
    proposal <- proposal.extend[-c(1, (G+2))]
        
    #### Compute the data likelihood
    lp.current <- lambda[Z] + offset + as.numeric(regression.mat) + as.numeric(phi.mat)
    lp.proposal <- proposal[Z] + offset + as.numeric(regression.mat) + as.numeric(phi.mat)
    p.current <- exp(lp.current) / (1 + exp(lp.current))
    p.proposal <- exp(lp.proposal) / (1 + exp(lp.proposal))
    like.current <- Y * log(p.current) + failures * log(1 - p.current)
    like.proposal <- Y * log(p.proposal) + failures * log(1 - p.proposal)
    prob <- exp(sum(like.proposal - like.current))
        if(prob > runif(1))
        {
        lambda <- proposal
        lambda.mat <- matrix(rep(lambda, N), nrow=N, byrow=TRUE)
        mu <- matrix(lambda[Z], nrow=K, ncol=N, byrow=FALSE)
        accept[1] <- accept[1] + 1  
        }else
        {}
        accept[2] <- accept[2] + 1           
        
        
        
    ##################     
    #### Sample from Z
    ##################
    prior.offset <- rep(NA, G)
        for(r in 1:G)
        {
        prior.offset[r] <-  log(sum(exp(-delta * ((1:G - r)^2 + (1:G - Gstar)^2)))) 
        }
    mu.offset <- offset.mat + regression.mat + phi.mat
    test <- Zupdatesqbin(Z=Z.mat, Offset=mu.offset, Y=Y.mat, delta=delta, lambda=lambda, nsites=K, ntime=N, G=G, SS=1:G, prioroffset=prior.offset, Gstar=Gstar, failures=failures.mat)          
    Z.mat <- test
    Z <- as.numeric(Z.mat)
    mu <- matrix(lambda[Z], nrow=K, ncol=N, byrow=FALSE)
        
        
        
    ######################
    #### Sample from delta
    ######################
    proposal.delta <-  rtruncnorm(n=1, a=1, b=prior.delta, mean=delta, sd=proposal.sd.delta)
    sum.delta1 <- sum((Z - Gstar)^2)
    sum.delta2 <- sum((Z.mat[ ,-1] - Z.mat[ ,-N])^2)
    current.fc1 <- -delta * (sum.delta1 + sum.delta2) - K *  log(sum(exp(-delta * (1:G - Gstar)^2))) 
    proposal.fc1 <- -proposal.delta * (sum.delta1 + sum.delta2) - K *  log(sum(exp(-proposal.delta * (1:G - Gstar)^2)))                 
    Z.temp <- matrix(rep(as.numeric(Z.mat[ ,-N]),G), ncol=G, byrow=FALSE)
    Z.temp2 <- (delta.update - Z.temp)^2 + (delta.update - Gstar)^2
    current.fc <- current.fc1 - sum(log(apply(exp(-delta * Z.temp2),1,sum)))
    proposal.fc <- proposal.fc1 - sum(log(apply(exp(-proposal.delta * Z.temp2),1,sum)))        
    prob <- exp(proposal.fc - current.fc)       
        if(prob > runif(1))
        {
        delta <- proposal.delta
        accept[3] <- accept[3] + 1  
        }else
        {}
    accept[4] <- accept[4] + 1   
        
        
        
    ####################
    #### Sample from phi
    ####################
    phi.offset <- mu + offset.mat + regression.mat
        if(MALA)
        {
        temp1 <- binomialarcarupdateMALA(W.triplet, W.begfin, W.triplet.sum,  K, N, phi.mat, tau2, gamma, 1, Y.mat, failures.mat, trials.mat, proposal.sd.phi, phi.offset, W.triplet.sum)      
        }else
        {
        temp1 <- binomialarcarupdateRW(W.triplet, W.begfin, W.triplet.sum,  K, N, phi.mat, tau2, gamma, 1, Y.mat, failures.mat, proposal.sd.phi, phi.offset, W.triplet.sum)      
        }
    phi.temp <- temp1[[1]]
    phi <- as.numeric(phi.temp)
            for(i in 1:G)
            {
            phi[which(Z==i)] <- phi[which(Z==i)] - mean(phi[which(Z==i)])
            }
    phi.mat <- matrix(phi, nrow=K, ncol=N, byrow=FALSE)
    accept[5] <- accept[5] + temp1[[2]]
    accept[6] <- accept[6] + K*N    
  
        
        
    ####################
    ## Sample from gamma
    ####################
    temp2 <- gammaquadformcompute(W.triplet, W.triplet.sum, W.n.triplet,  K, N, phi.mat, 1)
    mean.gamma <- temp2[[1]] / temp2[[2]]
    sd.gamma <- sqrt(tau2 / temp2[[2]]) 
    gamma <- rtruncnorm(n=1, a=0, b=1, mean=mean.gamma, sd=sd.gamma)   

        
        
    ####################
    ## Samples from tau2
    ####################
    temp3 <- tauquadformcompute(W.triplet, W.triplet.sum, W.n.triplet,  K, N, phi.mat, 1, gamma)
    tau2.posterior.scale <- temp3 + prior.tau2[2] 
    tau2 <- 1 / rgamma(1, tau2.posterior.shape, scale=(1/tau2.posterior.scale))          
        
        
        
    #########################
    ## Calculate the deviance
    #########################
    lp <- as.numeric(mu + offset.mat + regression.mat + phi.mat)
    prob <- exp(lp) / (1+exp(lp))
    fitted <- trials * prob
    loglike <- dbinom(x=Y, size=trials, prob=prob, log=TRUE)
        
        
    ###################
    ## Save the results
    ###################
        if(j > burnin & (j-burnin)%%thin==0)
        {
        ele <- (j - burnin) / thin
        samples.delta[ele, ] <- delta
        samples.lambda[ele, ] <- lambda
        samples.Z[ele, ] <- Z
        samples.phi[ele, ] <- as.numeric(phi.mat)
        samples.tau2[ele, ] <- tau2
        samples.gamma[ele, ] <- gamma
        samples.fitted[ele, ] <- fitted
        samples.loglike[ele, ] <- loglike
            if(!is.null(X)) samples.beta[ele, ] <- beta        
        }else
        {}
        

            
    ########################################
    ## Self tune the acceptance probabilties
    ########################################
    k <- j/100
        if(ceiling(k)==floor(k))
        {
            if(!is.null(X))
            {
                if(p>2)
                {
                proposal.sd.beta <- common.accceptrates1(accept[7:8], proposal.sd.beta, 40, 50)
                }else
                {
                proposal.sd.beta <- common.accceptrates1(accept[7:8], proposal.sd.beta, 30, 40)    
                }
            proposal.sd.phi <- common.accceptrates1(accept[5:6], proposal.sd.phi, 40, 50)
            proposal.sd.lambda <- common.accceptrates2(accept[1:2], proposal.sd.lambda, 20, 40, 10)
            proposal.sd.delta <- common.accceptrates2(accept[3:4], proposal.sd.delta, 40, 50, prior.delta/6)
            accept.all <- accept.all + accept
            accept <- rep(0,8)  
            }else
            {
            proposal.sd.phi <- common.accceptrates1(accept[5:6], proposal.sd.phi, 40, 50)
            proposal.sd.lambda <- common.accceptrates2(accept[1:2], proposal.sd.lambda, 20, 40, 10)
            proposal.sd.delta <- common.accceptrates2(accept[3:4], proposal.sd.delta, 40, 50, prior.delta/6)
            accept.all <- accept.all + accept
            accept <- rep(0,6)     
            }
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
accept.lambda <- 100 * accept.all[1] / accept.all[2]
accept.delta <- 100 * accept.all[3] / accept.all[4]
accept.phi <- 100 * accept.all[5] / accept.all[6]
accept.gamma <- 100
    if(!is.null(X))
    {
    accept.beta <- 100 * accept.all[7] / accept.all[8]   
    accept.final <- c(accept.beta, accept.lambda, accept.delta, accept.phi, accept.gamma)
    names(accept.final) <- c("beta", "lambda", "delta", "phi", "rho.T")   
    }else
    {
    accept.final <- c(accept.lambda,  accept.delta, accept.phi, accept.gamma)
    names(accept.final) <- c("lambda", "delta", "phi", "rho.T")   
    }
    
    
    
#### Compute the fitted deviance
mean.Z <- round(apply(samples.Z,2,mean), 0)       
mean.lambda <- apply(samples.lambda, 2, mean)
mean.mu <- matrix(mean.lambda[mean.Z], nrow=K, ncol=N, byrow=FALSE)
    if(!is.null(X))
    {
    mean.beta <- apply(samples.beta,2,mean)
    regression.mat <- matrix(X.standardised %*% mean.beta, nrow=K, ncol=N, byrow=FALSE)     
    }else
    {}
    
mean.phi <- matrix(apply(samples.phi, 2, mean), nrow=K, byrow=FALSE)
lp.mean <- as.numeric(mean.mu + offset.mat + mean.phi + regression.mat)   
mean.prob <- exp(lp.mean)  / (1 + exp(lp.mean))
fitted.mean <- trials * mean.prob
deviance.fitted <- -2 * sum(dbinom(x=Y, size=trials, prob=mean.prob, log=TRUE))

    
#### Model fit criteria
modelfit <- common.modelfit(samples.loglike, deviance.fitted)


#### Create the fitted values and residuals
fitted.values <- apply(samples.fitted, 2, mean)
response.residuals <- as.numeric(Y) - fitted.values
pearson.residuals <- response.residuals /sqrt(fitted.values * (1 - mean.prob))
residuals <- data.frame(response=response.residuals, pearson=pearson.residuals)

    
#### Transform the parameters back to the original covariate scale
    if(!is.null(X))
    {    
    samples.beta.orig <- common.betatransform(samples.beta, X.indicator, X.mean, X.sd, p, FALSE)
    }else
    {}    
    
    
#### Create a summary object
summary.hyper <- array(NA, c(3, 7))     
summary.hyper[1,1:3] <- quantile(samples.delta, c(0.5, 0.025, 0.975))
summary.hyper[2,1:3] <- quantile(samples.tau2, c(0.5, 0.025, 0.975))
summary.hyper[3,1:3] <- quantile(samples.gamma, c(0.5, 0.025, 0.975))
rownames(summary.hyper) <- c("delta", "tau2", "rho.T")      
summary.hyper[1, 4:7] <- c(n.keep, accept.delta, effectiveSize(mcmc(samples.delta)), geweke.diag(mcmc(samples.delta))$z)   
summary.hyper[2, 4:7] <- c(n.keep, 100, effectiveSize(mcmc(samples.tau2)), geweke.diag(mcmc(samples.tau2))$z)   
summary.hyper[3, 4:7] <- c(n.keep, 100, effectiveSize(mcmc(samples.gamma)), geweke.diag(mcmc(samples.gamma))$z)   
    
summary.lambda <- array(NA, c(G,1))
summary.lambda <- t(apply(samples.lambda, 2, quantile, c(0.5, 0.025, 0.975)))
summary.lambda <- cbind(summary.lambda, rep(n.keep, G), rep(accept.lambda, G), effectiveSize(mcmc(samples.lambda)), geweke.diag(mcmc(samples.lambda))$z)
summary.lambda <- matrix(summary.lambda, ncol=7)
rownames(summary.lambda) <- paste("lambda", 1:G, sep="")
    
    if(!is.null(X))
    {
    samples.beta.orig <- mcmc(samples.beta.orig)
    summary.beta <- t(apply(samples.beta.orig, 2, quantile, c(0.5, 0.025, 0.975))) 
    summary.beta <- cbind(summary.beta, rep(n.keep, p), rep(accept.beta,p), effectiveSize(samples.beta.orig), geweke.diag(samples.beta.orig)$z)
    rownames(summary.beta) <- colnames(X)
    colnames(summary.beta) <- c("Median", "2.5%", "97.5%", "n.sample", "% accept", "n.effective", "Geweke.diag")
    summary.results <- rbind(summary.beta, summary.lambda, summary.hyper)    
    }else
    {
    summary.results <- rbind(summary.lambda, summary.hyper)    
    }
    
summary.results[ , 1:3] <- round(summary.results[ , 1:3], 4)
summary.results[ , 4:7] <- round(summary.results[ , 4:7], 1)
colnames(summary.results) <- c("Median", "2.5%", "97.5%", "n.sample", "% accept", "n.effective", "Geweke.diag")    
    
    
    
#### Compile and return the results
#### Harmonise samples in case of them not being generated
    if(is.null(X)) samples.beta.orig = NA
    
samples <- list(beta=mcmc(samples.beta.orig), lambda=mcmc(samples.lambda),  Z=mcmc(samples.Z), delta=mcmc(samples.delta), phi = mcmc(samples.phi), tau2=mcmc(samples.tau2), rho.T=mcmc(samples.gamma), fitted=mcmc(samples.fitted))
model.string <- c("Likelihood model - Binomial (logit link function)", "\nLatent structure model - Localised autoregressive CAR model\n")
results <- list(summary.results=summary.results, samples=samples, fitted.values=fitted.values, residuals=residuals, modelfit=modelfit, accept=accept.final, localised.structure=mean.Z, formula=formula, model=model.string,  X=X)
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
