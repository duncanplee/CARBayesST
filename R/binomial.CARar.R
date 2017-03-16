binomial.CARar <- function(formula, data=NULL, trials, W, burnin, n.sample, thin=1,  prior.mean.beta=NULL, prior.var.beta=NULL, prior.tau2=NULL, fix.rho.S=FALSE, rho.S=NULL, fix.rho.T=FALSE, rho.T=NULL, MALA=TRUE, verbose=TRUE)
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
Y.miss <- frame.results$Y.miss
which.miss <- frame.results$which.miss
n.miss <- frame.results$n.miss  

    
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
failures.miss <- failures
failures.miss[which.miss==0] <- median(failures, na.rm=TRUE)


#### Check on the rho arguments
    if(!is.logical(fix.rho.S)) stop("fix.rho.S is not logical.", call.=FALSE)   
    if(fix.rho.S & is.null(rho.S)) stop("rho.S is fixed but an initial value was not set.", call.=FALSE)   
    if(fix.rho.S & !is.numeric(rho.S) ) stop("rho.S is not numeric.", call.=FALSE)  
    if(!is.logical(fix.rho.T)) stop("fix.rho.T is not logical.", call.=FALSE)   
    if(fix.rho.T & is.null(rho.T)) stop("rho.T is fixed but an initial value was not set.", call.=FALSE)   
    if(fix.rho.T & !is.numeric(rho.T) ) stop("rho.T is not numeric.", call.=FALSE)  

    if(fix.rho.S)
    {
    rho <- rho.S
    }else
    {
    rho <- runif(1)       
    }
    if(fix.rho.T)
    {
    gamma <- rho.T
    }else
    {
    gamma <- runif(1)       
    }   
    if(rho<0 ) stop("rho.S is outside the range [0, 1].", call.=FALSE)  
    if(rho>1 ) stop("rho.S is outside the range [0, 1].", call.=FALSE)  
    if(gamma<0 ) stop("rho.T is outside the range [0, 1].", call.=FALSE)  
    if(gamma>1 ) stop("rho.T is outside the range [0, 1].", call.=FALSE)  



#### CAR quantities
W.quants <- common.Wcheckformat.leroux(W, fix.rho.S, rho)
K <- W.quants$n
N <- N.all / K
W <- W.quants$W
W.triplet <- W.quants$W.triplet
W.n.triplet <- W.quants$n.triplet
W.triplet.sum <- W.quants$W.triplet.sum
n.neighbours <- W.quants$n.neighbours 
W.begfin <- W.quants$W.begfin


#### Priors
    if(is.null(prior.mean.beta)) prior.mean.beta <- rep(0, p)
    if(is.null(prior.var.beta)) prior.var.beta <- rep(1000, p)
    if(is.null(prior.tau2)) prior.tau2 <- c(1, 0.01)
prior.beta.check(prior.mean.beta, prior.var.beta, p)
prior.var.check(prior.tau2)


#### Compute the blocking structure for beta     
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
tau2 <- var(phi)/10



###############################    
#### Set up the MCMC quantities    
###############################
#### Matrices to store samples 
n.keep <- floor((n.sample - burnin)/thin)
samples.beta <- array(NA, c(n.keep, p))
samples.phi <- array(NA, c(n.keep, N.all))
samples.tau2 <- array(NA, c(n.keep, 1))
    if(!fix.rho.S) samples.rho <- array(NA, c(n.keep, 1))
    if(!fix.rho.T) samples.gamma <- array(NA, c(n.keep, 1))
samples.fitted <- array(NA, c(n.keep, N.all))
samples.like <- array(NA, c(n.keep, N.all))
samples.deviance <- array(NA, c(n.keep, 1))
    if(n.miss>0) samples.Y <- array(NA, c(n.keep, n.miss))
    
    
#### Specify the Metropolis quantities
accept.all <- rep(0,6)
accept <- accept.all
proposal.sd.phi <- 0.1
proposal.sd.rho <- 0.05
proposal.sd.beta <- 0.01
proposal.corr.beta <- solve(t(X.standardised) %*% X.standardised)
chol.proposal.corr.beta <- chol(proposal.corr.beta)     
tau2.shape <- prior.tau2[1] + N.all/2
   

 
#############################
#### Specify spatial elements
#############################
#### Spatial determinant
    if(!fix.rho.S) 
    {
    Wstar <- diag(apply(W,1,sum)) - W
    Wstar.eigen <- eigen(Wstar)
    Wstar.val <- Wstar.eigen$values
    det.Q.W <-  0.5 * sum(log((rho * Wstar.val + (1-rho))))     
    }else
    {}


#### Specify quantities that do not change
offset.mat <- matrix(offset, nrow=K, ncol=N, byrow=FALSE) 
regression.mat <- matrix(X.standardised %*% beta, nrow=K, ncol=N, byrow=FALSE)
Y.mat <- matrix(Y, nrow=K, ncol=N, byrow=FALSE)
trials.mat <- matrix(trials, nrow=K, ncol=N, byrow=FALSE)
Y.mat.miss <- matrix(Y.miss, nrow=K, ncol=N, byrow=FALSE)
which.miss.mat <- matrix(which.miss, nrow=K, ncol=N, byrow=FALSE)
failures.mat <- matrix(failures, nrow=K, ncol=N, byrow=FALSE)
failures.mat.miss <- matrix(failures.miss, nrow=K, ncol=N, byrow=FALSE)
phi.mat <- matrix(phi, nrow=K, ncol=N, byrow=FALSE)   

 
#### Check for islands
W.list<- mat2listw(W)
W.nb <- W.list$neighbours
W.islands <- n.comp.nb(W.nb)
islands <- W.islands$comp.id
n.islands <- max(W.islands$nc)
    if(rho==1 & gamma==1) 
    {
    tau2.phi.shape <- prior.tau2[1] + prior.tau2[1] + ((N-1) * (K-1))/2
    }else if(rho==1)
    {
    tau2.phi.shape <- prior.tau2[1] + prior.tau2[1] + (N * (K-1))/2        
    }else if(gamma==1)
    {
    tau2.phi.shape <- prior.tau2[1] + prior.tau2[1] + ((N-1) * K)/2          
    }else
    {}



###########################
#### Run the Bayesian model
###########################
#### Start timer
    if(verbose)
    {
    cat("Generating", n.keep, "post burnin and thinned (if requested) samples\n", sep = " ")
    progressBar <- txtProgressBar(style = 3)
    percentage.points<-round((1:100/100)*n.sample)
    }else
    {
    percentage.points<-round((1:100/100)*n.sample)     
    }
    
    
#### Matrices to store samples
    for(j in 1:n.sample)
    {
    ####################
    ## Sample from beta
    ####################
    offset.temp <- as.numeric(offset.mat + phi.mat)     
        if(p>2)
        {
        temp <- binomialbetaupdateMALA(X.standardised, N.all, p, beta, offset.temp, Y.miss, failures.miss, trials, prior.mean.beta, prior.var.beta, which.miss, n.beta.block, proposal.sd.beta, list.block)
        }else
        {
        temp <- binomialbetaupdateRW(X.standardised, N.all, p, beta, offset.temp, Y.miss, failures.miss, prior.mean.beta, prior.var.beta, which.miss, proposal.sd.beta)
        }
    beta <- temp[[1]]
    accept[1] <- accept[1] + temp[[2]]
    accept[2] <- accept[2] + n.beta.block  
    regression.mat <- matrix(X.standardised %*% beta, nrow=K, ncol=N, byrow=FALSE)           


    
    ####################
    ## Sample from phi
    ####################
    phi.offset <- offset.mat + regression.mat
    den.offset <- rho * W.triplet.sum + 1 - rho
        if(MALA)
        {
        temp1 <- binomialarcarupdateMALA(W.triplet, W.begfin, W.triplet.sum,  K, N, phi.mat, tau2, gamma, rho, Y.mat.miss, failures.mat.miss, trials.mat, proposal.sd.phi, phi.offset, den.offset, which.miss.mat)      
        }else
        {
        temp1 <- binomialarcarupdateRW(W.triplet, W.begfin, W.triplet.sum,  K, N, phi.mat, tau2, gamma, rho, Y.mat.miss, failures.mat.miss, proposal.sd.phi, phi.offset, den.offset, which.miss.mat)      
         }
    phi.temp <- temp1[[1]]
    phi <- as.numeric(phi.temp)  - mean(as.numeric(phi.temp))
    phi.mat <- matrix(phi, nrow=K, ncol=N, byrow=FALSE)
    accept[3] <- accept[3] + temp1[[2]]
    accept[4] <- accept[4] + K*N    
        
        

    ####################
    ## Sample from gamma
    ####################
        if(!fix.rho.T)
        {
        temp2 <- gammaquadformcompute(W.triplet, W.triplet.sum, W.n.triplet,  K, N, phi.mat, rho)
        mean.gamma <- temp2[[1]] / temp2[[2]]
        sd.gamma <- sqrt(tau2 / temp2[[2]])
        gamma <- rtruncnorm(n=1, a=0, b=1, mean=mean.gamma, sd=sd.gamma)  
        }else
        {}
        

                
    ####################
    ## Samples from tau2
    ####################
    temp3 <- tauquadformcompute(W.triplet, W.triplet.sum, W.n.triplet,  K, N, phi.mat, rho, gamma)
    tau2.scale <- temp3 + prior.tau2[2] 
    tau2 <- 1 / rgamma(1, tau2.shape, scale=(1/tau2.scale)) 
        
        

    ##################
    ## Sample from rho
    ##################
        if(!fix.rho.S)
        {
        proposal.rho <- rtruncnorm(n=1, a=0, b=1, mean=rho, sd=proposal.sd.rho)   
        temp4 <- tauquadformcompute(W.triplet, W.triplet.sum, W.n.triplet,  K, N, phi.mat, proposal.rho, gamma)
        det.Q.W.proposal <- 0.5 * sum(log((proposal.rho * Wstar.val + (1-proposal.rho))))
        logprob.current <- N * det.Q.W - temp3 / tau2
        logprob.proposal <- N * det.Q.W.proposal - temp4 / tau2
        prob <- exp(logprob.proposal - logprob.current)
            if(prob > runif(1))
            {
            rho <- proposal.rho
            det.Q.W <- det.Q.W.proposal
            accept[5] <- accept[5] + 1           
            }else
            {}              
        accept[6] <- accept[6] + 1       
        }else
        {}
        
        
    
    #########################
    ## Calculate the deviance
    #########################
    lp <- as.numeric(offset.mat + regression.mat + phi.mat)
    prob <- exp(lp) / (1+exp(lp))
    fitted <- trials * prob
    deviance.all <- dbinom(x=Y, size=trials, prob=prob, log=TRUE)
    like <- exp(deviance.all)
    deviance <- -2 * sum(deviance.all, na.rm=TRUE)             
   

    
    ###################
    ## Save the results
    ###################
        if(j > burnin & (j-burnin)%%thin==0)
        {
        ele <- (j - burnin) / thin
        samples.beta[ele, ] <- beta
        samples.phi[ele, ] <- as.numeric(phi)
            if(!fix.rho.S) samples.rho[ele, ] <- rho
            if(!fix.rho.T) samples.gamma[ele, ] <- gamma
        samples.tau2[ele, ] <- tau2
        samples.deviance[ele, ] <- deviance
        samples.fitted[ele, ] <- fitted
        samples.like[ele, ] <- like
            if(n.miss>0) samples.Y[ele, ] <- rbinom(n=n.miss, size=trials[which.miss==0], prob=prob[which.miss==0])
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
            if(!fix.rho.S) proposal.sd.rho <- common.accceptrates2(accept[5:6], proposal.sd.rho, 40, 50, 0.5) 
        accept.all <- accept.all + accept
        accept <- rep(0,6)  
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
    cat("\nSummarising results")
    close(progressBar)
    }else
    {}
    
    

###################################
#### Summarise and save the results 
###################################
#### Compute the acceptance rates
accept.beta <- 100 * accept.all[1] / accept.all[2]
accept.phi <- 100 * accept.all[3] / accept.all[4]
    if(!fix.rho.S)
    {
    accept.rho <- 100 * accept.all[5] / accept.all[6]
    }else
    {
    accept.rho <- NA    
    }
accept.final <- c(accept.beta, accept.phi, accept.rho, 100)
names(accept.final) <- c("beta", "phi", "rho.S", "rho.T")

    
#### Deviance information criterion (DIC)
median.beta <- apply(samples.beta,2,median)
regression.mat <- matrix(X.standardised %*% median.beta, nrow=K, ncol=N, byrow=FALSE)   
median.phi <- matrix(apply(samples.phi, 2, median), nrow=K, ncol=N)
lp.median <- as.numeric(offset.mat + median.phi + regression.mat)   
median.prob <- exp(lp.median)  / (1 + exp(lp.median))
fitted.median <- trials * median.prob
deviance.fitted <- -2 * sum(dbinom(x=Y, size=trials, prob=median.prob, log=TRUE), na.rm=TRUE)
p.d <- median(samples.deviance) - deviance.fitted
DIC <- 2 * median(samples.deviance) - deviance.fitted     


#### Watanabe-Akaike Information Criterion (WAIC)
LPPD <- sum(log(apply(samples.like,2,mean)), na.rm=TRUE)
p.w <- sum(apply(log(samples.like),2,var), na.rm=TRUE)
WAIC <- -2 * (LPPD - p.w)


#### Compute the LMPL
CPO <- rep(NA, N.all)
    for(j in 1:N.all)
    {
    CPO[j] <- 1/median((1 / dbinom(x=Y[j], size=trials[j], prob=(samples.fitted[ ,j] / trials[j]))))    
    }
LMPL <- sum(log(CPO), na.rm=TRUE)  
    

#### Create the fitted values and residuals
fitted.values <- apply(samples.fitted, 2, median)
response.residuals <- as.numeric(Y) - fitted.values
pearson.residuals <- response.residuals /sqrt(fitted.values * (1 - median.prob))
deviance.residuals <- sign(response.residuals) * sqrt(2 * (Y * log(Y/fitted.values) + (trials-Y) * log((trials-Y)/(trials - fitted.values))))
residuals <- data.frame(response=response.residuals, pearson=pearson.residuals, deviance=deviance.residuals)

    
#### Transform the parameters back to the origianl covariate scale.
samples.beta.orig <- common.betatransform(samples.beta, X.indicator, X.mean, X.sd, p, FALSE)


#### Create a summary object
samples.beta.orig <- mcmc(samples.beta.orig)
summary.beta <- t(apply(samples.beta.orig, 2, quantile, c(0.5, 0.025, 0.975))) 
summary.beta <- cbind(summary.beta, rep(n.keep, p), rep(accept.beta,p), effectiveSize(samples.beta.orig), geweke.diag(samples.beta.orig)$z)
rownames(summary.beta) <- colnames(X)
colnames(summary.beta) <- c("Median", "2.5%", "97.5%", "n.sample", "% accept", "n.effective", "Geweke.diag")
    
summary.hyper <- array(NA, c(3, 7))    
rownames(summary.hyper) <- c("tau2", "rho.S", "rho.T")     
summary.hyper[1,1:3] <- quantile(samples.tau2, c(0.5, 0.025, 0.975))
summary.hyper[1, 4:7] <- c(n.keep, 100, effectiveSize(mcmc(samples.tau2)), geweke.diag(mcmc(samples.tau2))$z)     
    if(!fix.rho.S)
    {
    summary.hyper[2, 1:3] <- quantile(samples.rho, c(0.5, 0.025, 0.975))
    summary.hyper[2, 4:7] <- c(n.keep, accept.rho, effectiveSize(samples.rho), geweke.diag(samples.rho)$z)
    }else
    {
    summary.hyper[2, 1:3] <- c(rho, rho, rho)
    summary.hyper[2, 4:7] <- rep(NA, 4)
    }
    if(!fix.rho.T)
    {
    summary.hyper[3, 1:3] <- quantile(samples.gamma, c(0.5, 0.025, 0.975))
    summary.hyper[3, 4:7] <- c(n.keep, 100, effectiveSize(mcmc(samples.gamma)), geweke.diag(mcmc(samples.gamma))$z)       
    }else
    {
    summary.hyper[3, 1:3] <- c(gamma, gamma, gamma)
    summary.hyper[3, 4:7] <- rep(NA, 4)
    }   

summary.results <- rbind(summary.beta, summary.hyper)
summary.results[ , 1:3] <- round(summary.results[ , 1:3], 4)
summary.results[ , 4:7] <- round(summary.results[ , 4:7], 1)

    
#### Compile and return the results
loglike <- (-0.5 * deviance.fitted)
modelfit <- c(DIC, p.d, WAIC, p.w, LMPL, loglike)
names(modelfit) <- c("DIC", "p.d", "WAIC", "p.w", "LMPL", "loglikelihood")


#### Harmonise samples in case of them not being generated
    if(fix.rho.S & fix.rho.T)
    {
    samples.rhoext <- NA
    }else if(fix.rho.S & !fix.rho.T)
    {
    samples.rhoext <- samples.gamma
    names(samples.rhoext) <- "rho.T"
    }else if(!fix.rho.S & fix.rho.T)
    {
    samples.rhoext <- samples.rho  
    names(samples.rhoext) <- "rho.S"
    }else
    {
    samples.rhoext <- cbind(samples.rho, samples.gamma)
    colnames(samples.rhoext) <- c("rho.S", "rho.T")
    }
    if(n.miss==0) samples.Y = NA

samples <- list(beta=mcmc(samples.beta.orig), phi=mcmc(samples.phi),  rho=mcmc(samples.rhoext), tau2=mcmc(samples.tau2), fitted=mcmc(samples.fitted), Y=mcmc(samples.Y))
model.string <- c("Likelihood model - binomial (logit link function)", "\nLatent structure model - Autoregressive CAR model\n")
results <- list(summary.results=summary.results, samples=samples, fitted.values=fitted.values, residuals=residuals, modelfit=modelfit, accept=accept.final, localised.structure=NULL, formula=formula, model=model.string,  X=X)
class(results) <- "CARBayesST"


#### Finish by stating the time taken 
    if(verbose)
    {
    b<-proc.time()
    cat(" finished in ", round(b[3]-a[3], 1), "seconds")
    }else
    {}
return(results)
}
