poisson.CARadaptive <- function(formula, data = NULL, W, burnin, n.sample, thin = 1, prior.mean.beta = NULL, prior.var.beta = NULL, prior.tau2 = NULL, rho = NULL, epsilon = 0, MALA=FALSE, verbose = TRUE)
{ 
#### Verbose
a <- common.verbose(verbose)  

    
  blocksize.beta <- 5
  blocksize.v    <- 10
  z              <- which(W > 0, arr.ind = T)
  locs           <- z[which(z[,1] < z[,2]), ]
  char.locs      <- paste(locs[,1], ".", locs[,2], sep = "")
  n.edges        <- nrow(locs)
  
  # convert the supplied adjacency matrix into a spam matrix, if required.
  if(!is.symmetric.matrix(W)) stop("W is not symmetric.", call.=FALSE)
  #if(class(W) == "matrix") W <- as.spam(W)
  #if(!class(W) %in% c("matrix", "spam")) stop("W must be an object with class \"matrix\" or \"spam\"", call.=FALSE)  

  
  logit     <- function(p) log(p/(1-p))
  inv_logit <- function(v) 1/(1+exp(-v))
  
  
  #### Check on MALA argument
  if(length(MALA)!=1) stop("MALA is not length 1.", call.=FALSE)
  if(!is.logical(MALA)) stop("MALA is not logical.", call.=FALSE) 
  
  
  # interpret the formula
  frame <- try(suppressWarnings(model.frame(formula, data = data, na.action=na.pass)), silent=TRUE)
  if(class(frame)=="try-error") stop("the formula inputted contains an error, e.g the variables may be different lengths.", call.=FALSE)
  X <- try(suppressWarnings(model.matrix(object=attr(frame, "terms"), data=frame)), silent=TRUE)
  #if(class(X)=="try-error") stop("the covariate matrix contains inappropriate values.", call.=FALSE)
  if(sum(is.na(X))>0) stop("the covariate matrix contains missing 'NA' values.", call.=FALSE)
  
  # get summaries of the model matrix
  p       <- ncol(X)
  y       <- model.response(frame)
  which.miss <- as.numeric(!is.na(y))
  n.sites <- as.integer(nrow(W))
  n.time  <- as.integer(length(y)/n.sites)
  k       <- as.integer(round(n.sites*n.time, 0))
  
  #### Check and specify the priors
  if(is.null(prior.mean.beta)) prior.mean.beta <- rep(0, p)
  if(is.null(prior.var.beta)) prior.var.beta <- rep(100000, p)
  if(is.null(prior.tau2)) prior.tau2 <- c(1, 0.01)
  prior.beta.check(prior.mean.beta, prior.var.beta, p)
  prior.var.check(prior.tau2)
  
  
  
  # identify and error check the offset term, if it exists.
  offset <- try(model.offset(frame), silent=TRUE)
  #if(class(offset)=="try-error")   stop("the offset is not numeric.", call.=FALSE)
  if(is.null(offset))              offset <- rep(0,(n.time * n.sites))
  if(sum(is.na(offset))>0)         stop("the offset has missing 'NA' values.", call.=FALSE)
  if(!is.numeric(offset))          stop("the offset variable has non-numeric values.", call.=FALSE) 
  
  #### Format and check the MCMC quantities
  common.burnin.nsample.thin.check(burnin, n.sample, thin)


  ## Standardise the model matrix, 
  X.standardised   <- X
  X.sd             <- apply(X, 2, sd)
  X.mean           <- apply(X, 2, mean)
  X.indicator      <- rep(NA, p)       # To determine which parameter estimates to transform back
  for(j in 1:p){
    if(length(table(X[ ,j])) > 2){
      X.indicator[j] <- 1
      X.standardised[ ,j] <- (X[ ,j] - mean(X[ ,j])) / sd(X[ ,j])
    }else if(length(table(X[ ,j]))==1){
      X.indicator[j] <- 2
    }else{
      X.indicator[j] <- 0
    }
  } 
  
  # based on the blocksize.v provided create lists with relevent bits for untransformed edge parameter update
  if(is.numeric(blocksize.v)){
    ## Compute the blocking structure for v
    fromto     <- seq(0, n.edges, by = blocksize.v)
    fromto[1]  <- 0
    if(!n.edges %in% fromto) fromto <- c(fromto, n.edges)
    n.blocks   <- length(fromto) - 1
    blockinds  <- vector("list", length = n.blocks)
    for(i in 1:n.blocks)  blockinds[[i]]    <- (fromto[i] + 1):fromto[i + 1]
  } 
  
  # propose starting values for the adjacency elements (very close to 1)
  # current temporary version of the adacency is W_current
  v                                <- logit(rtruncnorm(n.edges, mean = 0.999, sd = 0.001, a = 0, b=1))
  v_15                             <- v - 15
  vqform_current                   <- sum(v_15^2)
  W_current                        <- W
  W_current[locs][1:n.edges]       <- inv_logit(v)
  W_current[locs[,2:1]][1:n.edges] <- inv_logit(v)
  
  # given the temporary adjacency, construct a temporary (Q.space) and proposal (Q.space.prop)
  # for the prior ICAR precision for phi.  Associated with these is the triplet form tripList.
  # get the cholesky of Q.space, and its determinant
  # if rho is not fixed, then ridge must be fixed
  rhofix <- rho
  rho                     <- ifelse(!is.null(rhofix), rhofix, 0.99)
  fixedridge              <- epsilon
  if(rho==1) fixedridge <- 0.0001
  if(!is.numeric(rho) ) stop("rho is fixed but is not numeric.", call.=FALSE)  
  if(rho<0 ) stop("rho is outside the range [0, 1].", call.=FALSE)  
  if(rho>1 ) stop("rho is outside the range [0, 1].", call.=FALSE)    
  
  tripList                <- vector("list", length = 2)
  tripList[[1]]           <- cbind(1:nrow(W_current), 1:nrow(W_current), rowSums(W_current) + fixedridge)
  tripList[[2]]           <- cbind(rbind(locs, locs[,2:1]), -rep(inv_logit(v), 2))
  Q.space.trip            <- rbind(tripList[[1]], tripList[[2]])
  Q.space.trip            <- updatetriplets_rho(trips = Q.space.trip, nsites = n.sites, rho_old = 1, rho_new = rho, fixedridge = fixedridge)  
  Q.space <- Q.space.prop <- spam(list(i = Q.space.trip[,1], j = Q.space.trip[,2], Q.space.trip[,3]))
  chol.Q.space            <- chol.spam(Q.space)
  Q.space.det.old         <- n.time*2*determinant(chol.Q.space, logarithm = T)$modulus
  
  # propose an initial value for alpha, the temporal correlation parameter
  # using alpha, create initial temporal precision matrices Q.time
  alpha <- 1
  if(n.time > 1){
    # this bit constructs Q.time, temporal precision, its determinant, and triplet form
    Q.block      <- as.spam(crossprod(diff(diag(n.time))))
    Q.block[1,1] <- Q.block[1,1] + 1
    Dg           <- diag.spam(diag.spam(Q.block))
    R            <- Q.block - Dg
    Dtime        <- diag.spam( c(rep(1,nrow(Q.block)-1), 0))
    Dg           <- Dg - Dtime
    Q.time       <- Dg + Dtime*alpha^2+ R*alpha
    Q.time[n.time,n.time] <- 1
    Q.det        <- determinant(Q.time, logarithm = T)
    detTime      <- as.numeric(0.5*n.sites*(Q.det$m)*(Q.det$s))
    Q.time.trip  <- Reduce("cbind", triplet(Q.time)) 
  }  else {
    # if n.time == 1, then Q is Q.space, and detTime is just replaced with 1.
    Q.time       <- 1
    detTime      <- 1
    Q.time.trip  <- matrix(rep(1, 3), ncol = 3) 
  }
  
  # MCMC parameter starting values
  phi_tune      <- 0.5
  W.tune        <- 1
  rho.tune      <- 0.1
  tau_v         <- 200
  prior.max.tau <- 1000
  increment     <- 0
  glm_mod       <- glm(y ~-1+X.standardised, family = "quasipoisson", offset = offset)
  beta.mean <- glm_mod$coefficients
  beta.sd <- sqrt(diag(summary(glm_mod)$cov.scaled))
  beta_par      <- rnorm(n=length(beta.mean), mean=beta.mean, sd=beta.sd)
  
  log.y <- log(y)
  log.y[y==0] <- -0.1  
  res.temp <- log.y - X.standardised %*% beta_par - offset
  res.sd <- sd(res.temp, na.rm=TRUE)/5
  phi           <- rnorm(k, mean=0, sd=res.sd)
  tau           <- var(phi)/10
  phiQphi       <- qform_ST(Qspace = Q.space.trip, Qtime = Q.time.trip, phi = phi, nsites = n.sites) 
  XB            <- X.standardised %*% beta_par
  tau_v.shape   <- (n.edges/2) +  prior.tau2[1]
  tau_phi_shape <- (n.sites*n.time/2) + prior.tau2[1]
  # general MCMC housekeeping
  n.save        <- ifelse(thin == 1, (n.sample - burnin), (n.sample - burnin) / thin)
  accept    <- rep(0, 8)

  # storage of parameters in the MCMC, 
  samples.beta  <- array(NA, c(n.save, p))
  samples.phi   <- array(NA, c(n.save, n.sites * n.time))
  samples.tau2  <- samples.vtau2 <- samples.alpha <- samples.rho <- matrix(0, n.save, 1)
  samples.v     <- matrix(0, ncol = n.edges, nrow = c(n.save, n.sites*n.time))
  samples.fit   <- array(NA, c(n.save, n.sites * n.time))  
  samples.loglike <- array(NA, c(n.save, n.sites*n.time))
  
  # turn off spam check options to speed things up (a bit)
  options(spam.cholsymmetrycheck = FALSE)
  options(spam.cholpivotcheck = FALSE)
  options(spam.safemode = c(F, F, F))
  
  
  ## Compute the blocking structure for beta     
  #if(blocksize.beta >= p){
  #  n.beta.block <- 1
  #  beta.beg <- 1
  #  beta.fin <- p
  #} else {
  #  n.standard <- 1 + floor((p-blocksize.beta) / blocksize.beta)
  #  remainder <- p - n.standard * blocksize.beta 
  #  if(remainder==0){
  #    beta.beg <- c(1,seq((blocksize.beta+1), p, blocksize.beta))
  #    beta.fin <- seq(blocksize.beta, p, blocksize.beta)
  #    n.beta.block <- length(beta.beg)
  #  } else {
  #    beta.beg <- c(1, seq((blocksize.beta+1), p, blocksize.beta))
  #    beta.fin <- c(seq((blocksize.beta), p, blocksize.beta), p)
  #    n.beta.block <- length(beta.beg)
  #  }
  #}    
  
  
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
  
  
  
  proposal.sd.beta        <- 0.01
  proposal.corr.beta      <- solve(t(X.standardised) %*% X.standardised)
  chol.proposal.corr.beta <- chol(proposal.corr.beta)     
  
  # the perm ordering is used to map the @entries slot ordering to the ordering used when 'triplet' is called
  perm           <- order(Q.space.trip[,1], Q.space.trip[,2])
  diags.space    <- which(Q.space.trip[perm,1] == Q.space.trip[perm,2])
  if(n.time > 1) diag.time      <- Reduce("cbind", triplet(diag.spam(n.time - 1)))
  time.last.diag <- which((Q.time.trip[,1] == Q.time.trip[,2]) & (Q.time.trip[,1] == n.time))
  lastblock      <- (k - n.sites+1):k
  firstblock     <- 1:n.sites
  
  ## Start timer
  n.keep <- floor((n.sample - burnin)/thin)
  if(verbose){
      cat("Generating", n.keep, "post burnin and thinned (if requested) samples.\n", sep = " ")
      progressBar            <- txtProgressBar(style = 3)
    percentage.points      <- round((1:100/100)*n.sample)
  } else percentage.points <- round((1:100/100)*n.sample)     
  
  
  #   -------------------------------------------------------------------------------------------
  #   START THE MCMC SAMPLING
  #   -------------------------------------------------------------------------------------------
  
  for(j in 1:n.sample){ 
    # START ITERATING, ONLY SAVE thin^th ITERATION
    save.iter <- j > burnin && ((j %% thin == 0) | thin == 0)
    if(save.iter) increment <- increment+1
    
    
    # update ALPHA
    if(n.time > 1){
      phifirst         <- phi[-firstblock]
      philast          <- phi[-lastblock]
      philastQphilast  <- qform_ST(Qspace = Q.space.trip, Qtime = diag.time, phi = philast, nsites = n.sites)   
      phifirstQphilast <- qform_ST_asym(Qspace = Q.space.trip, Qtime = diag.time, phi1 = phifirst, phi2 = philast, nsites = n.sites) 
      mu_alpha         <- phifirstQphilast/philastQphilast
      mu_sigmasq       <- tau/philastQphilast
      alpha            <- rtruncnorm(n=1, a=10^-5, b=1 - 10^-5,  mean=mu_alpha, sd = sqrt(mu_sigmasq))
      Q.time.trip      <- update_Qtime(Q.time.trip, alpha, time.last.diag - 1)
      phiQphi          <- qform_ST(Qspace = Q.space.trip, Qtime = Q.time.trip, phi = phi, nsites = n.sites)   
      detTime          <- determinant(Q.time, logarithm = TRUE)
      detTime          <- (detTime$m)*(detTime$s)
    }
    
    # Gibbs update of tau_v
    tau_scale  <- vqform_current/2 + prior.tau2[2]
    tau_v      <- 1/rtrunc(n=1, spec="gamma", a=0.000001, b=Inf, shape=tau_v.shape, scale=(1/tau_scale))
    v.proposal <- rtruncnorm(n = n.edges, a=-15, b=15,  mean = v, sd = W.tune)
    for(i in 1:n.blocks){
      # propose new v for the i^th block
      vnew                                 <- v
      block_inds                           <- blockinds[[i]]
      vnew[block_inds]                     <- v.proposal[block_inds] 
      # update the spatial precision matrix using c++ loop.  
      # This is efficient because changes are only made where vnew differs from v
      # combine the result back into triplet matrix (Q.space.trip.prop), and spam matrix (Q.space.prop)
      tripUpdate                           <- updatetripList2(Q.space.trip, vold = v, vnew = vnew, nedges = n.edges, 
                                                              nsites = n.sites, block = block_inds, 
                                                              block_length = length(block_inds), fixedridge = fixedridge, rho = rho) 
      Q.space.trip.prop                    <- tripUpdate[[1]]
      Q.space.trip.diff                    <- tripUpdate[[2]]
      # combine the result back into triplet matrix (Q.space.trip.prop), and spam matrix (Q.space.prop)      
      Q.space.prop@entries                 <- Q.space.trip.prop[perm,3]
      # acceptance ratio requires calculation of phi'Q_prop phi - phi'Q phi.
      # do this quickly by taking the difference between old and new triplets and working out the 
      # difference directly.  Much faster than working out quadratic forms seperately.
      Q.space.trip.diff[, 3]<- Q.space.trip[, 3] - Q.space.trip.prop[,3]
      phiQphi_phiQphiNew   <- qform_difference_ST(Qtrip = Q.space.trip.diff, Qtime = Q.time.trip, phi = phi, nsites = n.sites)
      # update the cholesky of the precision matrix & calculate the determinant
      chol.Q.space.prop    <- update(chol.Q.space, x = Q.space.prop) 
      detSpace             <- 2*determinant(chol.Q.space.prop, logarithm = T)$modulus
      Q.space.det.prop     <- n.sites*detTime + n.time*detSpace
      v_15_prop            <- vnew - 15
      vqform_prop          <- sum(v_15_prop^2)
      acceptance           <- exp(0.5*(Q.space.det.prop - Q.space.det.old) + (1/(2*tau))*(phiQphi_phiQphiNew) 
                                  + 0.5*(1/tau_v)*(vqform_current - vqform_prop))
      accept[8]      <- accept[8] + (1/n.blocks)
      if(runif(1)  <= acceptance){
        vqform_current   <- vqform_prop
        v                <- vnew
        accept[7]        <- accept[7] + (1/n.blocks)
        Q.space.det.old  <- Q.space.det.prop
        Q.space.trip     <- Q.space.trip.prop
        chol.Q.space     <- chol.Q.space.prop
        Q.space          <- Q.space.prop
      }
    }
    
    # update BETA
    offset.temp   <- offset + as.numeric(phi)         
        if(MALA)
        {
        temp <- poissonbetaupdateMALA(X.standardised, k, p, beta_par, offset.temp, y, prior.mean.beta, prior.var.beta, n.beta.block, proposal.sd.beta, list.block)
        }else
        {
        temp <- poissonbetaupdateRW(X.standardised, k, p, beta_par, offset.temp, y, prior.mean.beta, prior.var.beta, n.beta.block, proposal.sd.beta, list.block)
        }
    beta_par <- temp[[1]]
    accept[1] <- accept[1] + temp[[2]]
    accept[2] <- accept[2] + n.beta.block  
    XB <- X.standardised %*% beta_par     
    
    
    
    # update PHI using one at a time M-H sampling
    nneighbours   <- diag.spam(Q.space)
    W_current     <- diag(nneighbours) - as.matrix(Q.space)
    phi_update    <- SPTICARphiVarb(W = W_current, nsites = n.sites, ntimes = n.time, phiVarb = phi, 
                                            nneighbours = nneighbours, tau = tau, y = y,  E = offset, 
                                            phiVarb_tune = phi_tune, alpha = alpha, XB = XB)    
    phi       <- phi_update[[2]]
    phi       <- phi - mean(phi)
    accept[3] <- accept[3] + phi_update[[1]][2]
    accept[4] <- accept[4] + k
    

    # update rho, the spatial leroux parameter
    if(!is.null(rhofix)){
      proposal.rho <- rhofix
    } else {
      proposal.rho           <- rtruncnorm(n = 1, a=0, b=1, mean = rho, sd = rho.tune) 
    }
    Q.space.trip.prop      <- updatetriplets_rho(trips = Q.space.trip, nsites = n.sites, rho_old = rho, rho_new = proposal.rho, fixedridge = fixedridge)   
    Q.space.prop@entries   <- Q.space.trip.prop[perm,3]
    Q.space.trip.diff[, 3] <- Q.space.trip[, 3] - Q.space.trip.prop[,3]
    phiQphi_phiQphiNew     <- qform_difference_ST(Qtrip = Q.space.trip.diff, Qtime = Q.time.trip, phi = phi, nsites = n.sites)
    # update the cholesky of the precision matrix & calculate the determinant
    chol.Q.space.prop      <- update(chol.Q.space, x = Q.space.prop) 
    detSpace               <- 2*determinant(chol.Q.space.prop, logarithm = T)$modulus
    Q.space.det.prop       <- n.sites*detTime + n.time*detSpace
    hastings <- log(dtruncnorm(x=rho, a=0, b=1, mean=proposal.rho, sd=rho.tune)) - log(dtruncnorm(x=proposal.rho, a=0, b=1, mean=rho, sd=rho.tune)) 
    acceptance             <- exp(0.5*(Q.space.det.prop - Q.space.det.old) + (1/(2*tau))*(phiQphi_phiQphiNew) + hastings)
    
    accept[6]              <- accept[6] + 1
    if(runif(1)  <= acceptance){
      accept[5]        <- accept[5] + 1
      Q.space.det.old  <- Q.space.det.prop
      Q.space.trip     <- Q.space.trip.prop
      chol.Q.space     <- chol.Q.space.prop
      Q.space          <- Q.space.prop
      rho              <- proposal.rho
    }
  
    # Gibbs update TAU using the gamma distribution
    phiQphi    <- qform_ST(Qspace = Q.space.trip, Qtime = Q.time.trip, phi = phi, nsites = n.sites)     
    tau_scale  <- phiQphi/2 + prior.tau2[2]
    tau        <- 1/rtrunc(n=1, spec="gamma", a=0.000001, b=Inf, shape=tau_phi_shape, scale=(1/tau_scale))
    
    # calculate the deviance and fitted values
    fitted     <- exp(as.vector(XB) + phi + offset)
    loglike <- dpois(x=as.numeric(y), lambda=fitted, log=TRUE)

    
    # save samples if past burnin 
    if(save.iter){
      samples.beta[increment,]      <- beta_par
      samples.phi[increment,]       <- phi
      samples.fit[increment, ]      <- fitted
      samples.tau2[increment,]      <- tau
      samples.vtau2[increment,]     <- tau_v
      samples.v[increment,]         <- v
      samples.alpha[increment,]     <- alpha
      samples.rho[increment,]       <- rho
      samples.loglike[increment, ] <- loglike
    }
    
    # adjust the acceptance rate if required
    if(j %% 100 == 0 & j < burnin){
        accept.beta <- 100 * accept[1] / accept[2]
        accept.phi <- 100 * accept[3] / accept[4]
        accept.w <- 100 * accept[7] / accept[8]
        if(is.null(rhofix))
        {
            accept.rho <- 100 * accept[5] / accept[6]     
        }else
        {
            accept.rho <- 45 
        }
        
        #### beta tuning parameter
        if(accept.beta > 50)
        {
            proposal.sd.beta <- proposal.sd.beta + 0.1 * proposal.sd.beta
        }else if(accept.beta < 40)              
        {
            proposal.sd.beta <- proposal.sd.beta - 0.1 * proposal.sd.beta
        }else
        {
        }
        
        #### phi tuning parameter
        if(accept.phi > 50)
        {
            phi_tune <- phi_tune + 0.1 * phi_tune
        }else if(accept.phi < 40)              
        {
            phi_tune <- phi_tune - 0.1 * phi_tune
        }else
        {
        }       
        
        #### w tuning parameter
        if(accept.w > 40)
        {
            W.tune <- W.tune + 0.1 * W.tune
        }else if(accept.w < 20)              
        {
            W.tune <- W.tune - 0.1 * W.tune
        }else
        {
        }   
        
        #### rho tuning parameter
        if(accept.rho > 50)
        {
            rho.tune <- min(rho.tune + 0.1 * rho.tune, 0.5)
        }else if(accept.rho < 40)              
        {
            rho.tune <- rho.tune - 0.1 * rho.tune
        }else
        {
        }  
        accept             <- accept*0
    }else
    {}
    
    
    # print progress to the console
    if(j %in% percentage.points & verbose) setTxtProgressBar(progressBar, j/n.sample)
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
accept.beta  <- 100 * accept[1] / accept[2]
accept.phi   <- 100 * accept[3] / accept[4]
accept.rho <- 100 * accept[5] / accept[6]
accept.w     <- 100 * accept[7] / accept[8]
accept.alpha <- 100
  if(!is.null(rhofix))
  {
  accept.final <- c(accept.beta, accept.phi, accept.w)
  names(accept.final) <- c("beta", "phi", "w")  
  }else
  {
  accept.final <- c(accept.beta, accept.phi, accept.rho,accept.w)
  names(accept.final) <- c("beta", "phi", "rho", "w")  
  }
  
  
#### Compute the fitted deviance  
mean.beta        <- apply(samples.beta, 2, mean)
regression.mat     <- matrix(X.standardised %*% mean.beta, nrow = n.sites, ncol = n.time, byrow=FALSE)   
mean.phi         <- matrix(apply(samples.phi, 2, mean), nrow = n.sites, ncol = n.time)
offset.mat         <- matrix(offset, nrow = n.sites, ncol = n.time, byrow=FALSE) 
fitted.mean      <- as.numeric(exp(offset.mat + mean.phi + regression.mat))
deviance.fitted    <- -2 * sum(dpois(x=as.numeric(y), lambda=fitted.mean, log=TRUE))

  
#### Model fit criteria
modelfit <- common.modelfit(samples.loglike, deviance.fitted)
  
  
#### Create the fitted values and residuals
fitted.values <- apply(samples.fit, 2, mean)
response.residuals <- as.numeric(y) - fitted.values
pearson.residuals <- response.residuals /sqrt(fitted.values)
residuals <- data.frame(response=response.residuals, pearson=pearson.residuals)

  
#### transform the parameters back to the origianl covariate scale.
samples.beta.orig <- common.betatransform(samples.beta, X.indicator, X.mean, X.sd, p, FALSE)

  
  #### Create a summary object
  samples.beta.orig       <- mcmc(samples.beta.orig)
  summary.beta            <- t(apply(samples.beta.orig, 2, quantile, c(0.5, 0.025, 0.975))) 
  summary.beta            <- cbind(summary.beta, rep(n.save, p), rep(accept.beta,p), effectiveSize(samples.beta.orig), geweke.diag(samples.beta.orig)$z)
  rownames(summary.beta)  <- colnames(X)
  colnames(summary.beta)  <- c("Median", "2.5%", "97.5%", "n.sample", "% accept", "n.effective", "Geweke.diag")
  
  
  summary.hyper           <- array(NA, c(4, 7))     
  summary.hyper[1,1:3]    <- quantile(samples.tau2, c(0.5, 0.025, 0.975))
  summary.hyper[2,1:3]    <- quantile(samples.rho, c(0.5, 0.025, 0.975))
  summary.hyper[3,1:3]    <- quantile(samples.alpha, c(0.5, 0.025, 0.975))
  summary.hyper[4,1:3]    <- quantile(samples.vtau2, c(0.5, 0.025, 0.975))
  rownames(summary.hyper) <- c("tau2", "rho.S", "rho.T", "tau2.w")     
  summary.hyper[1, 4:7]   <- c(n.save, 100, effectiveSize(mcmc(samples.tau2)), geweke.diag(mcmc(samples.tau2))$z)     
  summary.hyper[2, 4:7]   <- c(n.save, accept.rho, effectiveSize(mcmc(samples.rho)), geweke.diag(mcmc(samples.rho))$z)   
  summary.hyper[3, 4:7]   <- c(n.save, accept.alpha, effectiveSize(mcmc(samples.alpha)), geweke.diag(mcmc(samples.alpha))$z)   
  summary.hyper[4, 4:7]   <- c(n.save, 100, effectiveSize(mcmc(samples.vtau2)), geweke.diag(mcmc(samples.vtau2))$z)    

  if(!is.null(rhofix))
  {
  summary.hyper[2, ] <- c(rep(rhofix, 3),rep(NA, 4))    
  }
  summary.results         <- rbind(summary.beta, summary.hyper)
  summary.results[ , 1:3] <- round(summary.results[ , 1:3], 4)
  summary.results[ , 4:7] <- round(summary.results[ , 4:7], 1)
  
  # convert v back to w, summarise and create a 'fitted' adjacency matrix
  samples.w <- inv_logit(samples.v)
  colnames(samples.w) <- char.locs
  get_prop_thresh <- function(v, thresh) as.numeric(!((sum(v < thresh)/length(v)) < 0.99))
  bdry99          <- apply(samples.w, 2, get_prop_thresh, thresh = 0.5)
  bdryMN          <- apply(samples.w, 2, mean)
  Wmn <- W99      <- matrix(NA, nrow = n.sites, ncol = n.sites)
  W99[locs]       <- bdry99
  W99[locs[ ,c(2,1)]] <- bdry99
  Wmn[locs]       <- bdryMN
  Wmn[locs[ ,c(2,1)]] <- bdryMN  
  
  
#### Compile and return the results
model.string    <- c("Likelihood model - Poisson (log link function)", 
                       "\nLatent structure model - Adaptive autoregressive order 1 CAR model\n")
  samples.tau2all <- cbind(samples.tau2, samples.vtau2)
  colnames(samples.tau2all) <- c("tau2", "tau2.w")
  if(is.null(rhofix))
  {
      samples.rhoext <- cbind(samples.rho, samples.alpha)
      colnames(samples.rhoext) <- c("rho.S", "rho.T")
  }else
  {
      samples.rhoext <- cbind(samples.alpha)
      names(samples.rhoext) <- c("rho.T")
  }
  
      
  samples         <- list(beta = mcmc(samples.beta.orig), phi = mcmc(samples.phi), rho = mcmc(samples.rhoext), 
                          tau2 = mcmc(samples.tau2all), w = mcmc(samples.w), fitted = mcmc(samples.fit))  
  localised.structure <- list(Wmedian = Wmn, W99 = W99)
  results <- list(summary.results=summary.results, samples=samples, fitted.values=fitted.values, residuals=residuals, modelfit=modelfit, accept=accept.final, localised.structure=localised.structure, formula=formula, model=model.string,  X=X)
  
  class(results) <- "CARBayesST"
  if(verbose)
  {
    b<-proc.time()
    cat("Finished in ", round(b[3]-a[3], 1), "seconds.\n")
  }else
  {}
  return(results)
}







