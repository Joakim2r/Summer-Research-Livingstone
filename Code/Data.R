set.seed(1234)
library(posteriordb)
library(bridgestan)
library(rmcmc)
library(mvtnorm)
library(coda)
library(bayesplot)
library(cmdstanr)
library(xtable)
library(posterior)
library(shiny)
library(bslib)

############################## DON'T TOUCH #####################################
Sys.setenv(
  CC      = "/usr/bin/clang",
  CXX     = "/usr/bin/clang++",
  SDKROOT = "/Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX.sdk"
)

############################## POSTERIOR DATABASE #####################################
pdb <- pdb_local("/Users/joakimderambures/posteriordb")
pos_names <- posterior_names(pdb)
########################## DATA/MODEL FUNC #####################################
posterior_draws <- function(name, refs=pos_names, seed = 1234) {
  p   <- posterior(name, pdb)
  rpd <- NULL
  rpd <- try(reference_posterior_draws(p), silent = TRUE)
  if (inherits(rpd, "try-error")) rpd <- NULL
  if (!is.null(rpd)) {
    X <- as_draws_matrix(rpd)
    truth_mean <- colMeans(X)
    truth_var  <- apply(X, 2, var)
  } else {
    truth_mean <- NULL
    truth_var  <- NULL
  }

  model <- StanModel$new(stan_code_file_path(p), data_file_path(p), seed = seed)
  target <- target_distribution_from_stan_model(model)
  d     <- model$param_unc_num()
  x0 <- rep(0,d)

  list(
    model = model,
    d = d,
    x0= x0,
    rpd = rpd,
    target=target,
    truth_mean=truth_mean,
    truth_var=truth_var
  )
}
