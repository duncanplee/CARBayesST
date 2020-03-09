context("check outputs of Gaussian CARBayesST models")


do_lattice_example <- function(fn, n.sample, burnin){
  x.easting <- 1:10
  x.northing <- 1:10
  Grid <- expand.grid(x.easting, x.northing)
  K <- nrow(Grid)
  N <- 10
  N.all <- K * N
  
  #### set up spatial neighbourhood matrix W
  W <-array(0, c(K,K))
  for(i in 1:K)
  {
    for(j in 1:K)
    {
      temp <- (Grid[i,1] - Grid[j,1])^2 + (Grid[i,2] - Grid[j,2])^2
      if(temp==1)  W[i,j] <- 1 
    }	
  }
  
  #### Simulate the elements in the linear predictor and the data
  x <- rnorm(n=N.all, mean=0, sd=1)
  beta <- 0.1
  Q.W <- 0.99 * (diag(apply(W, 2, sum)) - W) + 0.01 * diag(rep(1,K))
  Q.W.inv <- solve(Q.W)
  phi <- mvrnorm(n=1, mu=rep(0,K), Sigma=(0.1 * Q.W.inv))
  delta <- mvrnorm(n=1, mu=rep(0,K), Sigma=(0.1 * Q.W.inv))
  trend <- array(NA, c(K, N))
  time <-(1:N - mean(1:N))/N
  for(i in 1:K)
  {
    trend[i, ] <- phi[i] + delta[i] * time        
  }
  trend.vec <- as.numeric(trend)
  LP <-  x * beta + trend.vec
  Y  <- rnorm(n = length(LP), mean = LP, sd = 10^-1)
  if("G" %in% names(as.list(args(fn)))){
    try(do.call(fn, args = list(formula = Y~x, family="gaussian", W=W, burnin=burnin, n.sample = n.sample, G = 3)))
  } else {
    try(do.call(fn, args = list(formula = Y~x, family="gaussian", W=W, burnin=burnin, n.sample = n.sample)))
  }
}

# function to perform tests
test_outputs <- function(fn){
  # skip_on_cran()
  z <- do_lattice_example(fn, 200, 100)
  # check class of outputs
  expect_s3_class(z, "CARBayesST")
  expect_equal(class(z$summary.results), "matrix")
  expect_equal(class(z$samples), "list")
  expect_equal(class(z$fitted.values), "numeric")
  expect_equal(class(z$modelfit), "numeric")
  expect_equal(class(z$accept), "numeric")
  # different output behaviour expected for adaptive and localised models
  if(as.character(substitute(fn)) == c("ST.CARadaptive")){
    expect_equal(class(z$localised.structure), "list")
  } else if(as.character(substitute(fn)) == "ST.CARlocalised"){
    expect_equal(class(z$localised.structure), "numeric")
  } else {
    expect_equal(class(z$localised.structure), "NULL")
  }
  expect_equal(class(z$formula), "formula")
  expect_equal(class(z$X), "matrix")
  # check extractor functions work
  expect_is(coef(z), "numeric")
  expect_is(fitted(z), "numeric")
  expect_is(logLik(z), "numeric")
  expect_is(model.matrix(z), "matrix")
  expect_is(residuals(z, type = "response"), "numeric")
  expect_is(residuals(z, type = "pearson"), "numeric")
}


test_that("Gaussian ST.CARlinear model output is correct", {
  skip_on_cran()
  test_outputs(ST.CARlinear)
})

test_that("Gaussian ST.CARanova model output is correct", {
  skip_on_cran()
  test_outputs(ST.CARanova)
})

# GAUSSIAN NOT ALLOWED FOR THIS MODEL
# test_that("Poisson ST.CARsepspatial model output is correct", {
#   test_outputs(ST.CARsepspatial)
# })

test_that("Gaussian ST.CARar model output is correct", {
  skip_on_cran()
  test_outputs(ST.CARar)
})

# GAUSSIAN NOT ALLOWED FOR THIS MODEL
# test_that("Poisson ST.CARlocalised model output is correct", {
  # test_outputs(fn = ST.CARlocalised)
# })

test_that("Gaussian ST.CARadaptive model output is correct", {
  skip_on_cran()
  test_outputs(ST.CARadaptive)
})



