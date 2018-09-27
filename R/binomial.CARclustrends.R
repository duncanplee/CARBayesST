#--------------------------------------------------------------------------------------------------------------------------------------------------------
# Bayesian hierarchical mixed-effects model for clustering areas based on disease risk trends (Binomial)
#--------------------------------------------------------------------------------------------------------------------------------------------------------
binomial.CARclustrends <- function(formula, data=NULL, trials, W, burnin, n.sample, thin=1, trends=NULL, changepoint=NULL, knots=NULL,
                                   prior.mean.beta=NULL, prior.var.beta=NULL, prior.mean.gamma=NULL, prior.var.gamma=NULL,
                                   prior.lambda=NULL, prior.tau2=NULL, Nchains=4, verbose=TRUE)
{
  #------------------------------------------------------------------------------------------------------------------------------------------------------
  # Check on the verbose option
  #------------------------------------------------------------------------------------------------------------------------------------------------------
  a <- common.verbose(verbose) 
  #------------------------------------------------------------------------------------------------------------------------------------------------------
  # Check trends vector
  #------------------------------------------------------------------------------------------------------------------------------------------------------
  All.Trends <- c("Constant", "LD", "LI", "CP", "CT", "MD", "MI")
  Total.trends <- length(All.Trends) - 2 # minus 2 as can't include both LD/MD or LI/MI
  #------------------------------------------------------------------------------------------------------------------------------------------------------
  # Check that a trend vector has been given
  #------------------------------------------------------------------------------------------------------------------------------------------------------
  if(is.null(trends)) stop("At least two trends, with one being the constant trend, have to be given.", call.=FALSE)
  trends <- unique(trends)
  #------------------------------------------------------------------------------------------------------------------------------------------------------
  # Check that the constant trend is selected
  #------------------------------------------------------------------------------------------------------------------------------------------------------
  if((All.Trends[1] %in% trends) & length(trends) == 1 | !(All.Trends[1] %in% trends))
  {
    stop("The constant trend has to be selected alongside at least one other trend.", call.=FALSE)
  }
  #------------------------------------------------------------------------------------------------------------------------------------------------------
  # Check to see if correct trends inputted
  #------------------------------------------------------------------------------------------------------------------------------------------------------
  if(!all(trends %in% All.Trends)) stop("Incorrect trend selected.", call.=FALSE)
  #------------------------------------------------------------------------------------------------------------------------------------------------------
  # Check that both LI and MI are both not included
  #------------------------------------------------------------------------------------------------------------------------------------------------------
  if(all(c("LI", "MI") %in% trends)) stop("Select only one of LI or MI as the increasing trend.", call.=FALSE)
  #------------------------------------------------------------------------------------------------------------------------------------------------------
  # Check that both LD and MD are both not included
  #------------------------------------------------------------------------------------------------------------------------------------------------------
  if(all(c("LD", "MD") %in% trends)) stop("Select only one of LD or MD as the decreasing trend.", call.=FALSE)
  #------------------------------------------------------------------------------------------------------------------------------------------------------
  # Check that the changepoint is included and within the given time period
  #------------------------------------------------------------------------------------------------------------------------------------------------------
  if(any(c("CP", "CT") %in% trends) & is.null(changepoint)) stop("A changepoint needs to be included for the changepoint trends (CP, CT).", call.=FALSE)
  if(any(c("CP", "CT") %in% trends) & length(changepoint) != 1) stop("The changepoint should be a scalar.", call.=FALSE)
  if(any(c("CP", "CT") %in% trends) & !is.null(changepoint))
  {
    if(changepoint < 1) stop("The changepoint should be positive.", call.=FALSE)
  }
  #------------------------------------------------------------------------------------------------------------------------------------------------------
  # Check the number of knots for the monotonic trends
  #------------------------------------------------------------------------------------------------------------------------------------------------------
  if(any(c("MD", "MI") %in% trends) & is.null(knots)) stop("The number of knots has to be chosen for the monotonic trends (MD, MI).", call.=FALSE)
  if(any(c("MD", "MI") %in% trends) & length(knots) != 1) stop("The number of knots should be a scalar.", call.=FALSE)
  if(any(c("MD", "MI") %in% trends) & !is.null(knots))
  {
    if(knots < 1) stop("The number of knots should be positive.", call.=FALSE)
  }
  #------------------------------------------------------------------------------------------------------------------------------------------------------
  # The constant trend does not need to be included within the trends vector
  #------------------------------------------------------------------------------------------------------------------------------------------------------
  N.trends <- length(trends)
  #------------------------------------------------------------------------------------------------------------------------------------------------------
  # Set number of knots to 0 if monontonic trends not included
  #------------------------------------------------------------------------------------------------------------------------------------------------------
  if(!any(c("MD", "MI") %in% trends)) knots <- 0
  #------------------------------------------------------------------------------------------------------------------------------------------------------
  # Track positions of each of the possible trends
  #------------------------------------------------------------------------------------------------------------------------------------------------------
  Trends.pos <- c("Constant", "LD", "LI", rep("CP", 2), rep("CT", 2), rep("MD", knots + 1), rep("MI", knots + 1))
  Trends.pos.numeric <- c(1, 2, 3, rep(4, 2), rep(5, 2), rep(6, knots + 1), rep(7, knots + 1))
  Trends.pos.num <- length(Trends.pos)
  #------------------------------------------------------------------------------------------------------------------------------------------------------
  # Track positions of the chosen trends
  #------------------------------------------------------------------------------------------------------------------------------------------------------
  Trends.id <- which(Trends.pos %in% trends)
  Trends.sel <- length(Trends.pos[Trends.id])
  Trends.id <- Trends.id[Trends.id != 1]
  #------------------------------------------------------------------------------------------------------------------------------------------------------
  # Vector for the number of gamma parameters associated with each of the trends
  #------------------------------------------------------------------------------------------------------------------------------------------------------
  params.trends <- c(0, 1, 1, rep(1, 2), rep(1, 2), rep(1, knots + 1), rep(1, knots + 1))
  Total.params.trends <- sum(params.trends)
  params.selected <- sum(params.trends[Trends.id])
  #------------------------------------------------------------------------------------------------------------------------------------------------------
  # Matrix containing tracking positions of associated gamma parameters
  #------------------------------------------------------------------------------------------------------------------------------------------------------
  Trend.pairs <- matrix(c(1, 0,
                          2, 0,
                          3, 0,
                          4, 5,
                          6, 7), ncol = 2, byrow = TRUE)
  rownames(Trend.pairs) <- c("Constant", "LD", "LI", "CP", "CT")
  #------------------------------------------------------------------------------------------------------------------------------------------------------
  # Include corresponding information for the monotonic trends
  #------------------------------------------------------------------------------------------------------------------------------------------------------
  col.knots1 <- seq(from = max(Trend.pairs) + 1, to = max(Trend.pairs) + (2 * (knots + 1)), by = 1)
  col.knots2 <- c(0, rep(-1, knots), 0, rep(-1, knots))
  col.knots2[which(col.knots2 == 0)] <- col.knots1[which(col.knots2 == 0)]
  row.knots <- matrix(c(col.knots1, col.knots2), ncol = 2)
  rownames(row.knots) <- c(rep("MD", knots + 1), rep("MI", knots + 1))
  row.knots <- row.knots[which(rownames(row.knots) %in% trends), ]
  Trend.pairs <- rbind(Trend.pairs, row.knots)
  Trend.pairs <- Trend.pairs[which(rownames(Trend.pairs) %in% trends), ]
  n.sel <- nrow(Trend.pairs)
  #------------------------------------------------------------------------------------------------------------------------------------------------------
  # Update tracking positions for the selected gamma parameters
  #------------------------------------------------------------------------------------------------------------------------------------------------------
  Trend.pairs.update <- Trend.pairs
  CP.check <- 1
  for(i in 1:n.sel)
  {
    if(Trend.pairs[i, 2] == 0)
    {
      Trend.pairs.update[i, 1] <- i
    }else if(Trend.pairs[i, 2] > 0)
    {
      if(rownames(Trend.pairs)[i] %in% c("CP", "CT"))
      {
        Trend.pairs.update[i, 1] <- Trend.pairs.update[(i-1), 1] + CP.check
        Trend.pairs.update[i, 2] <- Trend.pairs.update[i, 1] + 1
        CP.check <- CP.check + 1
      }else if(rownames(Trend.pairs)[i] %in% c("MD", "MI"))
      {
        if(Trend.pairs.update[(i-1), 2] > 0)
        {
          Mono.check <- 2
        }else
        {
          Mono.check <- 1
        }
        Trend.pairs.update[i, ] <- Trend.pairs.update[(i-1), 1] + Mono.check
      }
    }else if(Trend.pairs[i, 2] < 0)
    {
      Trend.pairs.update[i, 1] <- Trend.pairs.update[(i-1), 1] + 1
    }
  }
  #------------------------------------------------------------------------------------------------------------------------------------------------------
  # Track positions of the gamma parameters selected by the given trends
  #------------------------------------------------------------------------------------------------------------------------------------------------------
  Trend.names <- rownames(Trend.pairs.update)
  gamma.pos <- rep(0, Trends.pos.num)
  pos.gamma <- unique(Trend.pairs.update[which(Trend.pairs %in% Trends.id)])
  gamma.pos[Trends.id] <- pos.gamma[order(pos.gamma)]
  #------------------------------------------------------------------------------------------------------------------------------------------------------
  # Check the number of MCMC chains is >= 2
  #------------------------------------------------------------------------------------------------------------------------------------------------------
  if(Nchains <= 1) stop("the number of chains has to be greater than or equal 2.", call.=FALSE)
  if(Nchains %% 1 != 0) stop("the number of chains needs to be an integer.", call.=FALSE)
  #------------------------------------------------------------------------------------------------------------------------------------------------------
  # Format the arguments and check for errors
  #------------------------------------------------------------------------------------------------------------------------------------------------------
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
  if(p>1) stop("No covariates are allowed in this model due to identifiability issues.", call.=FALSE)
  
  #------------------------------------------------------------------------------------------------------------------------------------------------------
  # Check that the changepoint is included and within time period
  #------------------------------------------------------------------------------------------------------------------------------------------------------
  if(any(c("CP", "CT") %in% trends) & !is.null(changepoint))
  {
    if(!(changepoint >= 1 & changepoint <=N)) stop("The changepoint needs to be within the time period.", call.=FALSE)
  }
  #------------------------------------------------------------------------------------------------------------------------------------------------------
  # Check that the number of knots is not greater than the number of time points
  #------------------------------------------------------------------------------------------------------------------------------------------------------
  if(any(c("MD", "MI") %in% trends) & !is.null(knots))
  {
    if(knots > N) stop("The number of knots cannot be greater than the number of time points.", call.=FALSE)
  }
  #------------------------------------------------------------------------------------------------------------------------------------------------------
  # Spatial quantities
  #------------------------------------------------------------------------------------------------------------------------------------------------------
  W.quants <- common.Wcheckformat.leroux(W)
  K <- W.quants$n
  N <- N.all / K
  W <- W.quants$W
  W.triplet <- W.quants$W.triplet
  W.n.triplet <- W.quants$n.triplet
  W.triplet.sum <- W.quants$W.triplet.sum
  W.neighbours <- W.quants$n.neighbours 
  W.begfin <- W.quants$W.begfin
  #------------------------------------------------------------------------------------------------------------------------------------------------------
  # Check and specify the priors
  #------------------------------------------------------------------------------------------------------------------------------------------------------
  if(is.null(prior.mean.beta)) prior.mean.beta <- rep(0, p)
  if(is.null(prior.var.beta)) prior.var.beta <- rep(100000, p)
  if(is.null(prior.mean.gamma)) prior.mean.gamma <- rep(0, params.selected)
  if(is.null(prior.var.gamma)) prior.var.gamma <- rep(100000, params.selected)
  prior.mean.trends <- rep(0, Trends.pos.num)
  prior.mean.trends[Trends.id] <- prior.mean.gamma
  prior.var.trends <- rep(1000,  Trends.pos.num)
  prior.var.trends[Trends.id] <- prior.var.gamma
  if(is.null(prior.lambda)) prior.lambda <- rep(1, N.trends)
  if(is.null(prior.tau2)) prior.tau2 <- c(1, 0.01)
  prior.beta.check(prior.mean.beta, prior.var.beta, p)
  prior.var.check(prior.tau2)
  if(length(prior.mean.gamma) != params.selected) stop("the prior mean for gamma is the wrong length.", call.=FALSE)
  if(!is.numeric(prior.mean.gamma)) stop("the prior mean for gamma is not numeric.", call.=FALSE)    
  if(sum(is.na(prior.mean.gamma)) != 0) stop("the prior mean for gamma is missing.", call.=FALSE)  
  if(prior.mean.trends[2] > 0) stop("the prior mean for the LD trend should be non-positive.", call.=FALSE)
  if(prior.mean.trends[3] < 0) stop("the prior mean for the LI trend should be non-negative.", call.=FALSE)
  if(prior.mean.trends[4] < 0) stop("the prior mean for the increase in CP trend should be non-negative.", call.=FALSE)
  if(prior.mean.trends[5] > 0) stop("the prior mean for the decrease in CP trend should be non-positive.", call.=FALSE)
  if(prior.mean.trends[6] > 0) stop("the prior mean for the decrease in CT trend should be non-positive.", call.=FALSE)
  if(prior.mean.trends[7] < 0) stop("the prior mean for the increase in CT trend should be non-negative.", call.=FALSE)
  if(any(prior.mean.trends[8:(8 + knots)] > 0)) stop("the prior mean for the MD trend should be non-positive.", call.=FALSE)
  if(any(prior.mean.trends[(8 + knots + 1):(8 + knots + 1) + knots] < 0)) stop("the prior mean for the MI trend should be non-negative.", call.=FALSE)
  if(length(prior.var.gamma)!= params.selected) stop("the prior variance for gamma is the wrong length.", call.=FALSE)    
  if(!is.numeric(prior.var.gamma)) stop("the prior variance for gamma is not numeric.", call.=FALSE)    
  if(sum(is.na(prior.var.gamma))!= 0) stop("the prior variance for gamma is missing.", call.=FALSE)    
  if(min(prior.var.gamma) <= 0) stop("the prior variance for gamma is less than zero", call.=FALSE)
  if(length(prior.lambda) != N.trends) stop("the prior value for lambda is the wrong length.", call.=FALSE)    
  if(!is.numeric(prior.lambda)) stop("the prior value for lambda is not numeric.", call.=FALSE)    
  if(sum(is.na(prior.lambda)) != 0) stop("the prior value for lambda has missing values.", call.=FALSE)  
  #------------------------------------------------------------------------------------------------------------------------------------------------------
  # Specify the initial parameter values
  #------------------------------------------------------------------------------------------------------------------------------------------------------
  dat <- cbind(Y, failures)
  beta <- glm(dat~X.standardised-1, offset=offset, family=binomial(link="logit"))$coefficients
  beta <- matrix(beta, nrow = p, ncol = Nchains)
  proposal.corr.beta <- solve(t(X.standardised) %*% X.standardised)
  chol.proposal.corr.beta <- chol(proposal.corr.beta) 
  #------------------------------------------------------------------------------------------------------------------------------------------------------
  # Different initial beta values for each chain
  #------------------------------------------------------------------------------------------------------------------------------------------------------
  proposal.sd.beta <- rep(0.1, Nchains)
  beta <- matcomp(chol.proposal.corr.beta, beta, proposal.sd.beta, p, Nchains)
  gamma <- array(0, c(Trends.sel, Nchains))
  for (i in Trends.id)
  {
    if(i == 2 | i == 5 | i == 6 | (i %in% 8:(8 + knots)))
    {
      gamma[gamma.pos[i], ] <- rtrunc(Nchains, spec = "norm", b = 0, mean = 0, sd = 0.1)
    }else if (i == 3 | i == 4 | i == 7 | (i %in% (8 + knots + 1):(8 + knots + 1 + knots)))
    {
      gamma[gamma.pos[i], ] <- rtrunc(Nchains, spec = "norm", a = 0, mean = 0, sd = 0.1)
    }
  }
  gamma.mat <- array(0, c(N.all, Trends.sel, Nchains))
  for (i in Trends.id) 
  {
    gamma.mat[,gamma.pos[i],] <- matN(gamma[gamma.pos[i], ], N.all, Nchains)
  }
  tau2 <- runif(Nchains, 0, 1)
  rho <- runif(Nchains, 0, 1)
  lambda <- t(rdirichlet(Nchains, prior.lambda))
  w <- array(NA, c(K, N.trends, Nchains))
  phi <- array(NA, c(K, Nchains))
  for (i in 1:Nchains)
  {
    w[, , i] <- t(rmultinom(K, 1, lambda[, i]))
    phi[, i] <- rnorm(K, mean = 0, sd = 0.01)
  }
  phi <- array(NA, c(K, Nchains))
  for (i in 1:Nchains)
  {
    phi[, i] <- rnorm(K, mean = 0, sd = 0.01)
  }
  kronN <- rep(1, N)
  phimat <- kronecker(kronN, phi)
  wmat <- kronecker(kronN, w)
  w.chain.mat <- matrix(aperm(w, c(1, 3, 2)), nrow = K * Nchains, ncol = N.trends) 
  #------------------------------------------------------------------------------------------------------------------------------------------------------
  # Compute the blocking structure for covariate beta's
  #------------------------------------------------------------------------------------------------------------------------------------------------------
  block.temp <- common.betablock(p)
  beta.beg  <- block.temp[[1]]
  beta.fin <- block.temp[[2]]
  n.beta.block <- block.temp[[3]]
  #------------------------------------------------------------------------------------------------------------------------------------------------------
  # MCMC quantities - burnin, n.sample, thin
  #------------------------------------------------------------------------------------------------------------------------------------------------------
  common.burnin.nsample.thin.check(burnin, n.sample, thin)  
  #------------------------------------------------------------------------------------------------------------------------------------------------------
  # Set up matrices to store samples
  #------------------------------------------------------------------------------------------------------------------------------------------------------
  n.keep <- floor((n.sample - burnin)/thin)
  samples.beta <- array(NA, c(n.keep, p, Nchains))
  samples.gamma <- array(NA, c(Trends.sel, n.keep, 1, Nchains))
  samples.w <- array(NA, c(n.keep, K, N.trends, Nchains))
  samples.lambda <- array(NA, c(n.keep, N.trends, Nchains))
  samples.tau2 <- array(NA, c(n.keep, 1, Nchains))
  samples.rho <- array(NA, c(n.keep, 1, Nchains))
  samples.phi <- array(NA, c(n.keep, K, Nchains))
  samples.fitted <- array(NA, c(n.keep, N.all, Nchains))
  samples.like <- array(NA, c(n.keep, N.all, Nchains))
  samples.deviance <- array(NA, c(n.keep, 1, Nchains))
  #------------------------------------------------------------------------------------------------------------------------------------------------------
  # Specify the Metropolis quantities
  #------------------------------------------------------------------------------------------------------------------------------------------------------
  accept.all <- rep(0, 2 * (Trends.sel + 1) * Nchains)
  accept <- accept.all
  begin.accept <- seq(from = 1, to = length(accept), by = 2)
  end.accept <- begin.accept + 1
  accept.blocks.num <- array(begin.accept, c(Nchains, 2))
  accept.blocks.den <- array(end.accept, c(Nchains, 2))
  accept.weight <- matrix(0, nrow = K, ncol = 2 * Nchains)
  accept.w.all <- matrix(0, nrow = K, ncol = 2 * Nchains)
  accept.phis <- matrix(0, nrow = K, ncol = 2 * Nchains)
  accept.phis.all <- matrix(0, nrow = K, ncol = 2 * Nchains)
  accept.gammas <- matrix(0, nrow = Trends.sel, ncol = 2 * Nchains)
  accept.gammas.all <- matrix(0, nrow = Trends.sel, ncol = 2 * Nchains)
  accept.couple <- rep(0, 2)
  couples <- accept.couple
  tau2.shape <- prior.tau2[1] + K/2
  #------------------------------------------------------------------------------------------------------------------------------------------------------
  # Create the determinant 
  #------------------------------------------------------------------------------------------------------------------------------------------------------
  Wstar <- diag(apply(W, 1, sum)) - W
  Wstar.eigen <- eigen(Wstar)
  Wstar.val <- Wstar.eigen$values
  det.Q.W <- Qdet(Nchains, rho, Wstar.val)
  #------------------------------------------------------------------------------------------------------------------------------------------------------
  # Specify quantities that do not change
  #------------------------------------------------------------------------------------------------------------------------------------------------------
  Y.mat <- matrix(Y, nrow=K, ncol=N, byrow=FALSE)
  failures.mat <- matrix(failures, nrow=K, ncol=N, byrow=FALSE)
  offset.mat <- matrix(offset, nrow=N.all, ncol=Nchains)
  tp <- rep(1:N, each=K)
  tp.mat <- matrix(tp, nrow=K, ncol=N)
  tp.mat.trends <- array(tp.mat, c(K, N, Trends.sel))
  tp.mat.trends <- aperm(tp.mat.trends, c(1, 3, 2))
  tpmat <- array(tp, c(N.all, Trends.sel, Nchains))
  #------------------------------------------------------------------------------------------------------------------------------------------------------
  # Update the matrix corresponding to time given trends CP/CT or MD/MI
  #------------------------------------------------------------------------------------------------------------------------------------------------------
  any.CP <- any(Trend.pairs[, 2] != 0)
  track.pos <- which(Trend.pairs[, 2] != 0)
  neg.pos <- which(Trend.pairs.update < 0, arr.ind = TRUE)
  Trend.pairs.update[neg.pos[, 1], 2] <- Trend.pairs.update[neg.pos[, 1], 1]
  track.add <- which(Trend.pairs[, 2] > 0)
  if(any(names(track.pos) %in% c("MD", "MI")))
  {
    track.0 <- track.pos[which(names(track.pos) %in% c("MD", "MI"))]
    track.0 <- track.0[which(track.0 %in% track.add)]
    track.add <- track.pos[-which(track.pos %in% track.0)]
    track.add <- Trend.pairs.update[track.add, 2]
  }else
  {
    track.add <- Trend.pairs.update[track.add, 2]
  }
  if(any.CP)
  {
    tp.pos <- Trend.pairs.update[track.pos, 2]
    if(any(names(tp.pos) %in% c("CP", "CT")))
    {
      tp.CP <- tp.pos[names(tp.pos) %in% c("CP", "CT")]
      tpmat[, tp.CP, ] <- tpmat[, tp.CP, ] - changepoint
      tpmat[tpmat < 0] <- 0
      tp.mat.trends[, tp.CP, ] <- tp.mat.trends[, tp.CP, ] - changepoint
      tp.mat.trends[tp.mat.trends < 0] <- 0
    }
    if(any(names(tp.pos) %in% c("MD", "MI")))
    {
      tp.CP <- tp.pos[names(tp.pos) %in% c("MD", "MI")]
      k.space <- seq(from = 1, to = N, length = knots + 2)
      k.space <- round(k.space[-c(1, (knots+2))], digits = 1)
      if(all(c("MD", "MI") %in% names(tp.pos)))
      {
        kmat.col <- 2 * knots
      }else
      {
        kmat.col <- knots
      }
      kmat <- matrix(k.space, nrow = N.all, ncol = kmat.col, byrow = TRUE)
      kmat.Nchains <- array(kmat, dim = c(N.all, kmat.col, Nchains))
      kmat.N <- matrix(k.space, nrow = K, ncol = kmat.col, byrow = TRUE)
      kmat.N <- array(kmat.N, dim = c(K, kmat.col, N))
      if(all(c("MD", "MI") %in% trends))
      {
        tp.pos.0 <- rep(NA, 2)
        tp.pos.0[1] <- which(names(tp.CP) == "MD")[1]
        tp.pos.0[2] <- which(names(tp.CP) == "MI")[1]
      } else if("MD" %in% trends & !("MI" %in% trends))
      {
        tp.pos.0 <- which(names(tp.CP) == "MD")[1]
      } else if("MI" %in% trends & !("MD" %in% trends))
      {
        tp.pos.0 <- which(names(tp.CP) == "MI")[1]
      }
      tp.pos.row <- tp.CP[tp.pos.0]
      tpmat[, tp.pos.row, ] <- tpmat[, tp.pos.row, ] / N
      tpmat[, tp.CP[-tp.pos.0], ] <- ((tpmat[, tp.CP[-tp.pos.0], ] - kmat.Nchains)^3) / N^3
      tpmat[tpmat < 0] <- 0
      kmax <- apply(tpmat[, tp.CP[-tp.pos.0], ], 2, max)
      kmax.mat <- matrix(kmax, nrow = N.all, ncol = kmat.col, byrow = TRUE)
      kmax.Nchains <- array(kmax.mat, dim = c(N.all, kmat.col, Nchains))
      kmax.N <- matrix(kmax, nrow = K, ncol = kmat.col, byrow = TRUE)
      kmax.N <- array(kmax.N, dim = c(K, kmat.col, N))
      tpmat[, tp.CP[-tp.pos.0], ] <- tpmat[, tp.CP[-tp.pos.0], ] / kmax.Nchains
      tp.mat.trends[, tp.pos.row, ] <- tp.mat.trends[, tp.pos.row, ] / N
      tp.mat.trends[, tp.CP[-tp.pos.0], ] <- ((tp.mat.trends[, tp.CP[-tp.pos.0], ] - kmat.N)^3) / N^3
      tp.mat.trends[tp.mat.trends < 0] <- 0
      tp.mat.trends[, tp.CP[-tp.pos.0], ] <- tp.mat.trends[, tp.CP[-tp.pos.0], ] / kmax.N
    }
  }
  #------------------------------------------------------------------------------------------------------------------------------------------------------
  # Keep track of the additional positions of the selected gamma parameters of the CP/CT and MD/MI trends
  #------------------------------------------------------------------------------------------------------------------------------------------------------
  Trends.chosen.names <- c("Constant", unique(Trends.pos[Trends.id]))
  New.trend.pos <- rep(NA, length(track.add))
    if(length(track.add) != 0)
    {
      for(i in 1:length(track.add))
      {
          New.trend.pos[i] <- which(Trends.chosen.names %in% names(track.add)[i])
      }
    }
  #------------------------------------------------------------------------------------------------------------------------------------------------------
  # tempering temperatures
  #------------------------------------------------------------------------------------------------------------------------------------------------------
  d.t <- Nchains / (Nchains + 4)
  temps <- tempupdate(Nchains, d.t) 
  #------------------------------------------------------------------------------------------------------------------------------------------------------
  # proposal standard deviations for M-H moves 
  #------------------------------------------------------------------------------------------------------------------------------------------------------
  proposal.sd.gamma <- matrix(0.1, nrow = Trends.sel, ncol = Nchains)
  proposal.sd.phi <- matrix(0.1, nrow = K, ncol = Nchains)
  proposal.sd.rho <- rep(0.01, Nchains)
  max.proposal.sd.rho <- 0.1
  min.proposal.sd.rho <- 0.001
  #------------------------------------------------------------------------------------------------------------------------------------------------------
  # begin/end of chains for use in c++ functions due to using arrays
  #------------------------------------------------------------------------------------------------------------------------------------------------------
  begin.chain <- seq(from = 1, to = K * Nchains, by = K)
  begin.chainN <- seq(from = 1, to = N.all * Nchains, by = N.all)
  beg.reg.chain <- seq(from = 1, to = N.all, by = K)
  log1 <- log(1)
  N.all.trends <- N.all * Trends.sel
  if(any.CP)
  {
    wmat.extend <- array(0, c(N.all, Trends.sel, Nchains))
    wmat.extend[, -track.add, ] <- wmat
    wmat.extend[, track.add, ] <- wmat[, New.trend.pos, ]
  }else
  {
    wmat.extend <- wmat
  }
  beg.trend <- seq(from = 1, to = N.all.trends, by = N.all)
  wmat.ar <- matrix(wmat.extend, nrow = N.all.trends, ncol = Nchains)
  gamma.mat.ar <- matrix(gamma.mat, nrow = N.all.trends, ncol = Nchains)
  tpmat.ar <- matrix(tpmat, nrow = N.all.trends, ncol = Nchains)
  trends.part <- offsetcompute(wmat.ar, gamma.mat.ar, tpmat.ar, Nchains, N.all, Trends.sel, beg.trend)
  #------------------------------------------------------------------------------------------------------------------------------------------------------
  # Run the Bayesian model
  #------------------------------------------------------------------------------------------------------------------------------------------------------
  # Start timer
  #------------------------------------------------------------------------------------------------------------------------------------------------------
  if(verbose)
  {
    cat("Generating", n.keep, "post burnin and thinned (if requested) samples\n", sep = " ")
    progressBar <- txtProgressBar(style = 3)
    percentage.points<-round((1:100/100)*n.sample)
  }else
  {
    percentage.points<-round((1:100/100)*n.sample)     
  }
  for(j in 1:n.sample)
  {
    #------------------------------------------------------------------------------------------------------------------------------------------------------
    # Sample from beta
    #------------------------------------------------------------------------------------------------------------------------------------------------------
    proposal <- matcomp(chol.proposal.corr.beta, beta, proposal.sd.beta, p, Nchains)
    proposal.beta <- beta
    offset.temp <- offset.mat + trends.part + phimat
    for(r in 1:n.beta.block)
    {
      proposal.beta[beta.beg[r]:beta.fin[r], ] <- proposal[beta.beg[r]:beta.fin[r], ]
      beta.linpred <- linpredcomputeNchains(X.standardised, N.all, p, beta, Nchains)
      proposal.linpred <- linpredcomputeNchains(X.standardised, N.all, p, proposal.beta, Nchains)
      prob <- binomialbetablockupdate(N.all, beta, proposal.beta, beta.linpred, proposal.linpred, offset.temp, Y, failures, prior.mean.beta,
                                      prior.var.beta, Nchains, temps, p)
      accept.beta.chain <- prob > runif(Nchains)
      beta[beta.beg[r]:beta.fin[r], accept.beta.chain] <- proposal.beta[beta.beg[r]:beta.fin[r], accept.beta.chain]
      accept[accept.blocks.num[, 1]] <- accept[accept.blocks.num[, 1]] + as.numeric(accept.beta.chain)
      proposal.beta[beta.beg[r]:beta.fin[r], !accept.beta.chain] <- beta[beta.beg[r]:beta.fin[r], !accept.beta.chain]
    }
    accept[accept.blocks.den[, 1]] <- accept[accept.blocks.den[, 1]] + n.beta.block
    regression.mat <- linpredcomputeNchains(X.standardised, N.all, p, beta, Nchains)
    #------------------------------------------------------------------------------------------------------------------------------------------------------
    # Sample trend gamma's
    #------------------------------------------------------------------------------------------------------------------------------------------------------
    W.areas <- apply(w, c(2, 3), sum)
    offset.temp <- offset.mat + regression.mat + trends.part + phimat
    for (i in Trends.id)
    {
      gamma.proposal <- gammaproposal(Nchains, gamma[gamma.pos[i], ], proposal.sd.gamma[gamma.pos[i], ], prior.var.trends[i],
                                      W.areas[which(Trends.chosen.names %in% Trends.pos[i])[1], ], i, knots)
      gamma.mat.proposal <- gamma.mat
      gamma.mat.proposal[, gamma.pos[i], ] <- matN(gamma.proposal, N.all, Nchains)  
      gamma.mat.proposal.ar <- matrix(gamma.mat.proposal, nrow = N.all.trends, ncol = Nchains)
      trends.proposal <- offsetcompute(wmat.ar, gamma.mat.proposal.ar, tpmat.ar, Nchains, N.all, Trends.sel, beg.trend)
      offset.proposal <- offset.mat + regression.mat + trends.proposal + phimat
      gamma.list <- binomialgammaupdate(N.all, gamma[gamma.pos[i], ], gamma.proposal, offset.temp, offset.proposal, Y, failures,
                                        prior.mean.trends[i], prior.var.trends[i], Nchains, temps)
      if(!all(gamma.list[[2]] == 0))
      {
        gamma[gamma.pos[i], ] <- gamma.list[[1]]
        gamma.mat[, gamma.pos[i], ] <- matN(gamma[gamma.pos[i], ], N.all, Nchains) 
        gamma.mat.ar <- matrix(gamma.mat, nrow = N.all.trends, ncol = Nchains)
        trends.part <- offsetcompute(wmat.ar, gamma.mat.ar, tpmat.ar, Nchains, N.all, Trends.sel, beg.trend)
        offset.temp <- offset.mat + regression.mat + trends.part + phimat
        accept.gammas[gamma.pos[i], accept.blocks.num[, 1]] <- accept.gammas[gamma.pos[i], accept.blocks.num[, 1]] + gamma.list[[2]]
      }
      accept.gammas[gamma.pos[i], accept.blocks.den[, 1]] <- accept.gammas[gamma.pos[i], accept.blocks.den[, 1]] + 1 
    }
    #------------------------------------------------------------------------------------------------------------------------------------------------------
    # Sample from w
    #------------------------------------------------------------------------------------------------------------------------------------------------------
    w.perm <- matrix(aperm(w, c(1, 3, 2)), nrow = K*Nchains, ncol = N.trends)
    w.props <- sample(N.trends)
    while (all(w.props == 1:N.trends))
    {
      w.props <- sample(N.trends)
    }
    w.proposal <- w.perm[, w.props]
    w.proposal.array <- array(w.proposal, c(K, Nchains, N.trends))
    w.proposal.array <- aperm(w.proposal.array, c(1, 3, 2))
    w.proposal.array <- kronecker(kronN, w.proposal.array)
    if(any.CP)
    {
      wmat.extend.proposal <- array(0, c(N.all, Trends.sel, Nchains))
      wmat.extend.proposal[, -track.add, ] <- w.proposal.array
      wmat.extend.proposal[, track.add, ] <- w.proposal.array[, New.trend.pos, ]
    }else
    {
      wmat.extend.proposal <- w.proposal.array
    }
    w.proposal.ar <- matrix(wmat.extend.proposal, nrow = N.all.trends, ncol = Nchains)
    trends.proposal <- offsetcompute(w.proposal.ar, gamma.mat.ar, tpmat.ar, Nchains, N.all, Trends.sel, beg.trend)
    offset.proposal <- offset.mat + regression.mat + trends.proposal + phimat
    w.list <- binomialwupdate(K, N, w.perm, offset.temp, offset.proposal, w.proposal, Y.mat, failures.mat, lambda, Nchains, temps, begin.chain,
                              beg.reg.chain, N.trends)
    if(!all(w.list[[2]] == 0))
    {
      w <- w.list[[1]]
      w.array <- array(w, c(K, Nchains, N.trends))
      w <- aperm(w.array, c(1, 3, 2))
      wmat <- kronecker(kronN, w)
      if(any.CP)
      {
        wmat.extend <- array(0, c(N.all, Trends.sel, Nchains))
        wmat.extend[, -track.add, ] <- wmat
        wmat.extend[, track.add, ] <- wmat[, New.trend.pos, ]
      }else
      {
        wmat.extend <- wmat
      }
      wmat.ar <- matrix(wmat.extend, nrow = N.all.trends, ncol = Nchains) 
      trends.part <- offsetcompute(wmat.ar, gamma.mat.ar, tpmat.ar, Nchains, N.all, Trends.sel, beg.trend)
      w.chain.mat <- matrix(aperm(w, c(1, 3, 2)), nrow = K * Nchains, ncol = N.trends)
      accept.weight[, accept.blocks.num[, 1]] <- accept.weight[, accept.blocks.num[, 1]] + w.list[[2]]
    }
    accept.weight[, accept.blocks.den[, 1]] <- accept.weight[, accept.blocks.den[, 1]] + 1 
    #------------------------------------------------------------------------------------------------------------------------------------------------------
    # Sample from lambda
    #------------------------------------------------------------------------------------------------------------------------------------------------------
    lambda.temp <- prior.lambda + apply(w, c(2, 3), sum)
    lambda <- lambdaupdate(Nchains, lambda.temp)
    #------------------------------------------------------------------------------------------------------------------------------------------------------
    # Sample from phi
    #------------------------------------------------------------------------------------------------------------------------------------------------------
    offset.temp <- offset.mat + regression.mat + trends.part
    phi.list <- binomialphiupdate(W.triplet, W.begfin, W.triplet.sum, K, N, phi, offset.temp, Y.mat, failures.mat, tau2, rho, Nchains,
                                  temps, proposal.sd.phi, beg.reg.chain)
    if(!all(phi.list[[2]] == 0))
    {
      phi.means <- apply(phi.list[[1]], 2, mean)
      phi <- phi.list[[1]] - matrix(phi.means, nrow = K, ncol = Nchains, byrow = TRUE)
      phimat <- kronecker(kronN, phi)
      accept.phis[, accept.blocks.num[, 1]] <- accept.phis[, accept.blocks.num[, 1]] + phi.list[[2]] 
    }
    accept.phis[, accept.blocks.den[, 1]] <- accept.phis[, accept.blocks.den[, 1]] + 1 
    #------------------------------------------------------------------------------------------------------------------------------------------------------
    # Samples from tau2
    #------------------------------------------------------------------------------------------------------------------------------------------------------
    tau2.temp <- tau2quadform(W.triplet, W.triplet.sum, W.n.triplet, K, phi, phi, rho, Nchains)
    tau2 <- tau2computeNchains(tau2.temp, tau2.shape, prior.tau2[2], Nchains)
    #------------------------------------------------------------------------------------------------------------------------------------------------------
    # Samples from rho
    #------------------------------------------------------------------------------------------------------------------------------------------------------
    rho.temp1 <- rhoquadformcomputeNchains(W.triplet, W.triplet.sum, W.n.triplet, K, Nchains, phi, rho, tau2)
    proposal.rho <- suppressWarnings(rtrunc(n = Nchains, spec = "norm", a = 0, b = 0.99, mean = rho, sd = proposal.sd.rho))
    rho.temp2 <- rhoquadformcomputeNchains(W.triplet, W.triplet.sum, W.n.triplet, K, Nchains, phi, proposal.rho, tau2)
    det.Q.W.proposal <- Qdet(Nchains, proposal.rho, Wstar.val)
    logprob.current <- det.Q.W - rho.temp1
    logprob.proposal <- det.Q.W.proposal - rho.temp2
    prob <- exp((logprob.proposal - logprob.current) * temps) # raised to temperature levels of each chain
    accept.rho.chain <- prob > runif(Nchains)
    rho[accept.rho.chain] <- proposal.rho[accept.rho.chain]
    det.Q.W[accept.rho.chain] <- det.Q.W.proposal[accept.rho.chain]
    accept[accept.blocks.num[, 2]] <- accept[accept.blocks.num[, 2]] + as.numeric(accept.rho.chain)
    accept[accept.blocks.den[, 2]] <- accept[accept.blocks.den[, 2]] + 1
    #------------------------------------------------------------------------------------------------------------------------------------------------------
    # Metropolis coupling
    #------------------------------------------------------------------------------------------------------------------------------------------------------
    swap <- sample(1:Nchains, 2)
    offset.temp <- offset.mat + regression.mat + trends.part + phimat
    accept.swap <- binomialcouplingAllupdate(N.all, K, p, w.chain.mat, offset.temp, beta, gamma, lambda, phi, rho, tau2, W.triplet.sum, W.triplet,
                                            W.begfin, Y, failures, prior.mean.beta, prior.var.beta, prior.mean.trends, prior.var.trends,
                                            prior.lambda, prior.tau2, swap, temps, begin.chain, N.trends, Trends.sel)
    if(accept.swap == 1)
    {
      rev.swap <- rev(swap)
      beta[rev.swap] <- beta[swap]
      regression.mat[, rev.swap] <- regression.mat[, swap]
      proposal.sd.beta[rev.swap] <- proposal.sd.beta[swap]
      gamma[, rev.swap] <- gamma[, swap]
      gamma.mat[, , rev.swap] <- gamma.mat[, , swap]
      gamma.mat.ar <- matrix(gamma.mat, nrow = N.all.trends, ncol = Nchains)
      proposal.sd.gamma[, rev.swap] <- proposal.sd.gamma[, swap]
      lambda[, rev.swap] <- lambda[, swap]
      w[, , rev.swap] <- w[, , swap]
      wmat[, , rev.swap] <- wmat[, , swap]
      if(any.CP)
      {
        wmat.extend <- array(0, c(N.all, Trends.sel, Nchains))
        wmat.extend[, -track.add, ] <- wmat
        wmat.extend[, track.add, ] <- wmat[, New.trend.pos, ]
      }else
      {
        wmat.extend <- wmat
      }
      w.chain.mat <- matrix(aperm(w, c(1, 3, 2)), nrow = K * Nchains, ncol = N.trends)
      wmat.ar <- matrix(wmat.extend, nrow = N.all.trends, ncol = Nchains)
      phi[, rev.swap] <- phi[, swap]
      proposal.sd.phi[, rev.swap] <- proposal.sd.phi[, swap]
      phimat[, rev.swap] <- phimat[, swap]
      tau2[rev.swap] <- tau2[swap]
      rho[rev.swap] <- rho[swap]
      proposal.sd.rho[rev.swap] <- proposal.sd.rho[swap]
      det.Q.W[rev.swap] <- det.Q.W[swap]
      trends.part <- offsetcompute(wmat.ar, gamma.mat.ar, tpmat.ar, Nchains, N.all, Trends.sel, beg.trend)
      offset.temp <- offset.mat + regression.mat + trends.part + phimat
    }else
    {}
    accept.couple[1] <- accept.couple[1] + accept.swap
    accept.couple[2] <- accept.couple[2] + 1
    #------------------------------------------------------------------------------------------------------------------------------------------------------
    # Update temperatures
    #------------------------------------------------------------------------------------------------------------------------------------------------------
    if(j%%10==0)
    {
      MC3.accept <- 100 * accept.couple[1] / accept.couple[2]
      if(MC3.accept > 30)
      {
        d.t <- max(runif(1, d.t * 0.8, d.t), 0.1)
        temps <- tempupdate(Nchains, d.t)
      }else if(MC3.accept < 20)
      {
        d.t <- min(runif(1, d.t, d.t * 1.2), 0.99)
        temps <- tempupdate(Nchains, d.t)
      }else
      {}
    }else
    {}
    #------------------------------------------------------------------------------------------------------------------------------------------------------
    # Calculate the deviance
    #------------------------------------------------------------------------------------------------------------------------------------------------------
    bin.probs <- exp(offset.temp) / (1 + exp(offset.temp))
    fitted <- trials * bin.probs
    dev.like <- binomialdevfit(Y, trials, bin.probs, N.all, Nchains)
    deviance <- dev.like[[1]]
    like <- dev.like[[2]]
    #------------------------------------------------------------------------------------------------------------------------------------------------------
    # Save the results
    #------------------------------------------------------------------------------------------------------------------------------------------------------
    if(j > burnin & (j-burnin)%%thin==0)
    {
      ele <- (j - burnin) / thin
      samples.beta[ele,,] <- beta
      samples.gamma[,ele,,] <- gamma
      samples.w[ele,,,] <- w
      samples.lambda[ele,,] <- lambda
      samples.tau2[ele,,] <- tau2
      samples.rho[ele,,] <- rho
      samples.phi[ele,,] <- phi
      samples.deviance[ele,,] <- deviance
      samples.fitted[ele,,] <- fitted
      samples.like[ele,,] <- like
    }else
    {
    }
    #------------------------------------------------------------------------------------------------------------------------------------------------------
    # Self tune the acceptance probabilties
    #------------------------------------------------------------------------------------------------------------------------------------------------------
    k <- j/100
    if(ceiling(k)==floor(k))
    {
      #------------------------------------------------------------------------------------------------------------------------------------------------------
      # Determine the acceptance probabilities
      #------------------------------------------------------------------------------------------------------------------------------------------------------
      accept.beta <- 100 * accept[accept.blocks.num[,1]] / accept[accept.blocks.den[,1]]
      accept.gamma <- 100 * accept.gammas[,accept.blocks.num[,1]] / accept.gammas[,accept.blocks.den[,1]]
      accept.gamma[1,] <- 0
      accept.gammas.all <- accept.gammas.all + accept.gammas
      accept.rho <- 100 * accept[accept.blocks.num[,2]] / accept[accept.blocks.den[,2]]
      accept.w <- 100 * accept.weight[,accept.blocks.num[,1]] / accept.weight[,accept.blocks.den[,1]]
      accept.w.all <- accept.w.all + accept.weight
      accept.phi <- 100 * accept.phis[,accept.blocks.num[,1]] / accept.phis[,accept.blocks.den[,1]]
      accept.phis.all <- accept.phis.all + accept.phis
      accept.all <- accept.all + accept
      accept <- rep(0, 2 * (Trends.sel + 1) * Nchains)
      accept.weight <- matrix(0, nrow=K, ncol=2*Nchains)
      accept.phis <- matrix(0, nrow=K, ncol=2*Nchains)
      accept.gammas <- matrix(0, nrow=Trends.sel, ncol=2*Nchains)
      #------------------------------------------------------------------------------------------------------------------------------------------------------
      # beta tuning parameter
      #------------------------------------------------------------------------------------------------------------------------------------------------------
      if(any(accept.beta > 50))
      {
        proposal.sd.beta[which(accept.beta > 50)] <- 2 * proposal.sd.beta[which(accept.beta > 50)]
      }else if(any(accept.beta < 40))         
      {
        proposal.sd.beta[which(accept.beta < 40)] <- 0.5 * proposal.sd.beta[which(accept.beta < 40)]
      }else
      {
      }
      #------------------------------------------------------------------------------------------------------------------------------------------------------
      # gamma tuning parameter
      #------------------------------------------------------------------------------------------------------------------------------------------------------
      if(any(accept.gamma > 50))
      {
        proposal.sd.gamma[which(accept.gamma > 50)] <- 2 * proposal.sd.gamma[which(accept.gamma > 50)]
      }else if(any(accept.gamma < 40))
      {
        proposal.sd.gamma[which(accept.gamma < 40)] <- 0.5 * proposal.sd.gamma[which(accept.gamma < 40)]
      }else
      {
      }
      #------------------------------------------------------------------------------------------------------------------------------------------------------
      # rho tuning parameter
      #------------------------------------------------------------------------------------------------------------------------------------------------------
      if(any(accept.rho > 50))
      {
        proposal.sd.rho[which(accept.rho > 50)] <- 2 * proposal.sd.rho[which(accept.rho > 50)]
        if(any(proposal.sd.rho > max.proposal.sd.rho))
        {
          proposal.sd.rho[which(proposal.sd.rho > max.proposal.sd.rho)] <- max.proposal.sd.rho
        }else
        {
        }
      }else if(any(accept.rho < 40))              
      {
        proposal.sd.rho[which(accept.rho < 40)] <- 0.5 * proposal.sd.rho[which(accept.rho < 40)]
        if(any(proposal.sd.rho < min.proposal.sd.rho))
        {
          proposal.sd.rho[which(proposal.sd.rho < min.proposal.sd.rho)] <- min.proposal.sd.rho
        }else
        {
        }
      }else
      {
      }
      #------------------------------------------------------------------------------------------------------------------------------------------------------
      # phi tuning parameter
      #------------------------------------------------------------------------------------------------------------------------------------------------------
      if(any(accept.phi > 50))
      {
        proposal.sd.phi[which(accept.phi > 50)] <- 2 * proposal.sd.phi[which(accept.phi > 50)]
      }else if(any(accept.phi < 40))              
      {
        proposal.sd.phi[which(accept.phi < 40)] <- 0.5 * proposal.sd.phi[which(accept.phi < 40)]
      }else
      {
      }
    }else
    {}
    #------------------------------------------------------------------------------------------------------------------------------------------------------    
    # Print progress to the console
    #------------------------------------------------------------------------------------------------------------------------------------------------------
    if(j %in% percentage.points)
    {
      setTxtProgressBar(progressBar, j/n.sample)
    }
  }
  #------------------------------------------------------------------------------------------------------------------------------------------------------
  # End timer
  #------------------------------------------------------------------------------------------------------------------------------------------------------
  if(verbose)
  {
    cat("\nSummarising results")
    close(progressBar)
  }else
  {}
  #------------------------------------------------------------------------------------------------------------------------------------------------------
  # Summarise and save the results 
  #------------------------------------------------------------------------------------------------------------------------------------------------------
  # Select untempered chain for inference
  #------------------------------------------------------------------------------------------------------------------------------------------------------
  chain.sel <- 1
  p.d <- DIC <- LMPL <- NA
  fitted.values <- residuals <- rep(NA, N.all)
  #------------------------------------------------------------------------------------------------------------------------------------------------------
  # Watanabe-Akaike Information Criterion (WAIC)
  #------------------------------------------------------------------------------------------------------------------------------------------------------
  LPPD <- sum(log(apply(samples.like[,,chain.sel],2,mean)), na.rm=TRUE)
  p.w <- sum(apply(log(samples.like[,,chain.sel]),2,var), na.rm=TRUE)
  WAIC <- -2 * (LPPD - p.w)
  #------------------------------------------------------------------------------------------------------------------------------------------------------
  # Compute information criterion (DIC, DIC3, WAIC)
  #------------------------------------------------------------------------------------------------------------------------------------------------------
  mode.w <- matrix(0, nrow=K, ncol=N.trends)
  Wsum <- apply(samples.w[,,,chain.sel],c(2,3),sum)
  Wtrend <- which(Wsum == rowMaxs(Wsum), arr.ind=TRUE)
  for (i in 1:K) {
    mode.w[Wtrend[i,1], Wtrend[i,2]] <- 1
  }
  mode.w <- array(mode.w, c(K,N.trends,N))
  mode.beta <- rep(NA, p)
  if(p == 1)
  {
    mode.beta <- density(samples.beta[,,chain.sel])
    mode.beta <- mean(mode.beta$x[which(mode.beta$y==max(mode.beta$y))])
  }else
  {
    for(i in 1:p)
    {
      betamode <- density(samples.beta[,i,chain.sel])
      mode.beta[i] <- mean(betamode$x[which(betamode$y==max(betamode$y))])
    }
  }
  reg.mat <- matrix(X.standardised %*% mode.beta, nrow=K, ncol=N, byrow=FALSE)
  gamma.mat <- array(0, c(K,Trends.sel,N))
  for(i in Trends.id)
  {
    gamma.dens <- density(samples.gamma[gamma.pos[i],,,chain.sel])
    gamma.mean <- mean(gamma.dens$x[which(gamma.dens$y==max(gamma.dens$y))])
    gamma.mat[,gamma.pos[i],] <- matN(rep(gamma.mean, N),K,N)
  }
  mode.phi <- rep(NA,K)
  for(i in 1:K)
  {
    phimode <- density(samples.phi[,i,chain.sel])
    mode.phi[i] <- mean(phimode$x[which(phimode$y==max(phimode$y))])
  }
  phi.mat <- matN(rep(mode.phi,N),K,N)
  wmat.extend <- array(0, c(K,Trends.sel,N))
  wmat.extend[,-track.add,] <- mode.w
  wmat.extend[,track.add,] <- mode.w[,New.trend.pos,]
  trends.part <- apply(wmat.extend * (gamma.mat * tp.mat.trends),c(1,3),sum)
  offset.temp <- as.numeric(offset.mat[,chain.sel] + reg.mat + trends.part + phi.mat)
  mode.prob <- exp(offset.temp) / (1 + exp(offset.temp))
  fitted.mode <- trials * mode.prob
  deviance.fitted <- -2 * sum(dbinom(x=Y, size=trials, prob=mode.prob, log=TRUE))
  p.d <- median(samples.deviance[,,chain.sel]) - deviance.fitted
  DIC <- 2 * median(samples.deviance[,,chain.sel]) - deviance.fitted     
  #------------------------------------------------------------------------------------------------------------------------------------------------------
  # Compute the LMPL
  #------------------------------------------------------------------------------------------------------------------------------------------------------
  CPO <- rep(NA, N.all)
  for(j in 1:N.all)
  {
    CPO[j] <- 1/median((1 / dbinom(x=Y[j], size=trials[j], prob=(samples.fitted[,j,chain.sel] / trials[j]))))  
  }
  LMPL <- sum(log(CPO))  
  #------------------------------------------------------------------------------------------------------------------------------------------------------
  # Create the Fitted values
  #------------------------------------------------------------------------------------------------------------------------------------------------------
  fitted.values <- apply(samples.fitted[,,chain.sel],2,mean)
  response.residuals <- as.numeric(Y) - fitted.values
  pearson.residuals <- response.residuals /sqrt(fitted.values * (1 - mode.prob))
  residuals <- data.frame(response=response.residuals, pearson=pearson.residuals)
  
  #------------------------------------------------------------------------------------------------------------------------------------------------------
  # Transform the parameters back to the original covariate scale
  #------------------------------------------------------------------------------------------------------------------------------------------------------
  samples.beta.orig <- common.betatransform(samples.beta[,,chain.sel], X.indicator, X.mean, X.sd, p, FALSE)
  #------------------------------------------------------------------------------------------------------------------------------------------------------
  # Compute the acceptance rates
  #------------------------------------------------------------------------------------------------------------------------------------------------------
  accept.beta <- 100 * accept.all[accept.blocks.num[,1]] / accept.all[accept.blocks.den[,1]]
  accept.beta <- accept.beta[chain.sel]
  accept.gammas <- 100 * accept.gammas.all[,accept.blocks.num[,1]] / accept.gammas.all[,accept.blocks.den[,1]]
  accept.gammas <- accept.gammas[,chain.sel]
  accept.rho <- 100 * accept.all[accept.blocks.num[,2]] / accept.all[accept.blocks.den[,2]]
  accept.rho <- accept.rho[chain.sel]
  accept.phis <- 100 * accept.phis.all[, accept.blocks.num[,1]] / accept.phis.all[,accept.blocks.den[,1]]
  accept.phis <- accept.phis[,chain.sel]
  coupled <- 100*accept.couple[1]/accept.couple[2]
  accept.final <- c(accept.beta, accept.gammas[-1], accept.rho, mean(accept.phis), coupled)
  names(accept.final) <- c("beta", paste("gamma.", Trends.pos[Trends.id], sep=""), "rho", "phi", "coupled")
  #------------------------------------------------------------------------------------------------------------------------------------------------------
  # Create a summary object
  #------------------------------------------------------------------------------------------------------------------------------------------------------
  samples.beta.orig <- mcmc(samples.beta.orig)
  mode.beta.orig <- rep(NA, p)
  HPD.beta.orig <- matrix(NA, nrow=2, ncol=p)
  if(p == 1)
  {
    mode.beta.orig <- density(samples.beta.orig)
    mode.beta.orig <- mean(mode.beta.orig$x[which(mode.beta.orig$y==max(mode.beta.orig$y))])
    HPD.beta.orig[1,] <- HPDinterval(samples.beta.orig, prob=0.95)[1]
    HPD.beta.orig[2,] <- HPDinterval(samples.beta.orig, prob=0.95)[2]
    summary.beta <- t(c(mode.beta.orig, HPD.beta.orig[1,], HPD.beta.orig[2,]))
  }else
  {
    summary.beta <- matrix(NA, nrow=p, ncol=3)
    for(i in 1:p)
    {
      origbetamode <- density(samples.beta.orig[,i])
      mode.beta.orig[i] <- mean(origbetamode$x[which(origbetamode$y==max(origbetamode$y))])
      HPD.beta.orig[1,i] <- HPDinterval(samples.beta.orig[,i], prob=0.95)[1]
      HPD.beta.orig[2,i] <- HPDinterval(samples.beta.orig[,i], prob=0.95)[2]
      summary.beta[i,1] <- mode.beta.orig[i]
      summary.beta[i,2] <- HPD.beta.orig[1,i]
      summary.beta[i,3] <- HPD.beta.orig[2,i]
    }
  }
  summary.beta <- cbind(summary.beta, rep(n.keep,p), rep(accept.beta,p), effectiveSize(samples.beta.orig), geweke.diag(samples.beta.orig)$z)
  rownames(summary.beta) <- colnames(X)
  colnames(summary.beta) <- c("Mode", "2.5%", "97.5%", "n.sample", "% accept", "n.effective", "Geweke.diag")
  summary.gamma <- matrix(NA, nrow=Trends.sel, ncol=7)
  for(i in Trends.id)
  {
    summary.gamma[gamma.pos[i],1] <- unique(as.numeric(gamma.mat[, gamma.pos[i],]))
    summary.gamma[gamma.pos[i],2:3] <- HPDinterval(mcmc(samples.gamma[gamma.pos[i],,,chain.sel]), prob=0.95)
    summary.gamma[gamma.pos[i],4] <- rep(n.keep,1)
    summary.gamma[gamma.pos[i],5] <- accept.gammas[gamma.pos[i]]
    summary.gamma[gamma.pos[i],6] <- effectiveSize(samples.gamma[gamma.pos[i],,,chain.sel])
    summary.gamma[gamma.pos[i],7] <- geweke.diag(samples.gamma[gamma.pos[i],,,chain.sel])$z
  }
  colnames(summary.gamma) <- c("Mode", "2.5%", "97.5%", "n.sample", "% accept", "n.effective", "Geweke.diag")
  rownames(summary.gamma) <- c("gamma.constant", paste("gamma.", Trends.pos[Trends.id], sep=""))
  summary.gamma <- summary.gamma[-1,]
  if(Trends.sel==2)
  {
    summary.gamma <- matrix(summary.gamma, nrow=1, ncol=7)
    colnames(summary.gamma) <- c("Mode", "2.5%", "97.5%", "n.sample", "% accept", "n.effective", "Geweke.diag")
    rownames(summary.gamma) <- paste("gamma.", Trends.pos[Trends.id], sep = "")
  }
  summary.lambda <- matrix(NA, nrow=N.trends, ncol=7)
  for(i in 1:N.trends)
  {
    lambda.dens <- density(samples.lambda[,i,chain.sel])
    lambda.mean <- mean(lambda.dens$x[which(lambda.dens$y==max(lambda.dens$y))])
    summary.lambda[i,1] <- lambda.mean
    summary.lambda[i,2:3] <- HPDinterval(mcmc(samples.lambda[,i,chain.sel]), prob=0.95)
    summary.lambda[i,4] <- rep(n.keep,1)
    summary.lambda[i,5] <- rep(100,1)
    summary.lambda[i,6] <- effectiveSize(samples.lambda[,i,chain.sel])
    summary.lambda[i,7] <- geweke.diag(samples.lambda[,i,chain.sel])$z
  }
  colnames(summary.lambda) <- c("Mode", "2.5%", "97.5%", "n.sample", "% accept", "n.effective", "Geweke.diag")
  rownames(summary.lambda) <- paste("lambda.", All.Trends[which(All.Trends %in% trends)], sep = "")
  mode.tau2 <- density(samples.tau2[,,chain.sel])
  mode.tau2 <- mean(mode.tau2$x[which(mode.tau2$y==max(mode.tau2$y))])
  summary.tau2 <- t(c(mode.tau2, HPDinterval(mcmc(samples.tau2[,,chain.sel]), prob=0.95)[1],
                      HPDinterval(mcmc(samples.tau2[,,chain.sel]), prob=0.95)[2]))
  summary.tau2 <- cbind(summary.tau2, rep(n.keep, 1), rep(100,1), effectiveSize(samples.tau2[,,chain.sel]),
                        geweke.diag(samples.tau2[,,chain.sel])$z)
  colnames(summary.tau2) <- c("Median", "2.5%", "97.5%", "n.sample", "% accept", "n.effective", "Geweke.diag")
  rownames(summary.tau2) <- c("tau2")
  mode.rho <- density(samples.rho[,,chain.sel])
  mode.rho <- mean(mode.rho$x[which(mode.rho$y==max(mode.rho$y))])
  summary.rho <- t(c(mode.rho, HPDinterval(mcmc(samples.rho[,,chain.sel]), prob=0.95)[1],
                     HPDinterval(mcmc(samples.rho[,,chain.sel]), prob=0.95)[2]))
  summary.rho <- cbind(summary.rho, rep(n.keep, 1), rep(accept.rho,1), effectiveSize(samples.rho[,,chain.sel]),
                       geweke.diag(samples.rho[,,chain.sel])$z)
  colnames(summary.rho) <- c("Median", "2.5%", "97.5%", "n.sample", "% accept", "n.effective", "Geweke.diag")
  rownames(summary.rho) <- c("rho")
  summary.results <- rbind(summary.beta, summary.gamma, summary.lambda, summary.tau2, summary.rho)
  summary.results[,1:3] <- round(summary.results[,1:3],4)
  summary.results[,4:7] <- round(summary.results[,4:7],1)
  #------------------------------------------------------------------------------------------------------------------------------------------------------
  # Allocated trends for each location
  #------------------------------------------------------------------------------------------------------------------------------------------------------
  trends <- apply(samples.w[,,,chain.sel], c(2,3), sum)
  trend.probs <- trends / n.keep
  trends <- which(trends==rowMaxs(trends), arr.ind=TRUE)
  trends <- trends[order(trends[,1]),]
  trends[ ,2] <- Trends.chosen.names[trends[ ,2]]
  #------------------------------------------------------------------------------------------------------------------------------------------------------
  # Compile and return the results
  #------------------------------------------------------------------------------------------------------------------------------------------------------
  loglike <- -0.5 * deviance.fitted
  modelfit <- c(DIC[chain.sel], p.d[chain.sel], WAIC[chain.sel], p.w[chain.sel], LMPL[chain.sel], loglike)
  names(modelfit) <- c("DIC", "p.d", "WAIC", "p.w", "LMPL", "loglikelihood")
  samples <- list(beta=mcmc(t(matrix(samples.beta.orig, ncol=n.keep))), gamma=mcmc(t(matrix(samples.gamma[-1,,,chain.sel], ncol=n.keep))), lambda=mcmc(samples.lambda[,,chain.sel]),
                  tau2=mcmc(as.matrix(samples.tau2[,,chain.sel])), rho=mcmc(as.matrix(samples.rho[,,chain.sel])), w=samples.w[,,,chain.sel],
                  phi=mcmc(samples.phi[,,chain.sel]), fitted=mcmc(samples.fitted[,,chain.sel]))
  model.string <- c("Likelihood model - binomial (logit link function)", "\nLatent structure model - spatial main effects and an area clustered trend\n")

  
  results <- list(summary.results=summary.results, samples=samples, fitted.values=fitted.values, residuals=residuals, modelfit=modelfit, 
                  accept=accept.final, localised.structure=list(trends=trends[ ,-1], trend.probs=trend.probs), formula=formula, model=model.string, X=X)
  class(results) <- "CARBayesST"
  if(verbose)
  {
    b<-proc.time()
    cat(" finished in ", round(b[3]-a[3], 1), "seconds")
  }else
  {}
  return(results)
}
#------------------------------------------------------------------------------------------------------------------------------------------------------