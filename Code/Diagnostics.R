################################################################################
############################### DIAG FUNC TESTING###############################
basic_diagnostics <- function(draws, probs = c(0.05, 0.5, 0.95), digits = 3) {

  tab <- summarise_draws(
    draws,
    mean,
    mcse_mean,
    sd,
    mad,
    ~posterior::quantile2(.x, probs = probs),  # adds q5, q50, q95
    ess_bulk,
    ess_tail,
    rhat
  )

  return(tab)
}

######################### DIAG #############################################

diagnostics <- function(chains, ess, truth_mean = NULL, truth_var = NULL) {

  compute_mcse_mean <- function(chains, ess) {
    pooled <- do.call(rbind,chains)
    apply(pooled, 2, sd) / sqrt(ess)
  }

  compute_esjd <- function(chains) {
    sq_list <- lapply(chains, function(ch) rowSums(diff(ch)^2)) # this takes the squared difference between rows, then computes the sum across the columns for each row. sq_list should be a n*1 after operations.
    esjd_per_chain <- sapply(sq_list, mean) # takes the mean of the squared differences between the rows. This should return a scalar : 1*1
    running_per_chain <- lapply(sq_list, function(v) cumsum(v) / seq_along(v))  #  returns a list of numeric vectors, Each vector being the running mean for each chain
    running_avg <- Reduce(`+`, running_per_chain) / length(running_per_chain) # reduce adds vectors element wise. So you combine the running means for all the chains, and then divide by the number of chains you have
                                                                              # length(running_per_chain) = number of chains

    list(
      per_chain_mean = esjd_per_chain,    # numeric vector (one per chain)
      running_per_chain  = running_per_chain, # list of curves
      running_avg = running_avg        # averaged curve
    )
  }

  mse_first_moment_avg_over_chains <- function(chains, truth_mean) {
    n <- nrow(chains[[1]]) # gives you the number of iterations
    truth_mat <- matrix(rep(truth_mean, each = n), nrow = n) # repeats truth mean on every row
    per_chain <- lapply(chains, function(M) {
      rs <- apply(M, 2, cumsum)  # gets the running sum for each of the chains
      rm <- sweep(rs, 1, seq_len(n), "/")  # sweep basically goes row by row, and apply the function '/' which means dividing, and it divides everything by seq_len(n), so each sum is divided by its own length
      rowMeans((rm - truth_mat)^2) # length n
    })
    Reduce(`+`, per_chain) / length(per_chain) # average over chains
  }

  compute_dt_curve <- function(chains, truth_var) {
    n <- nrow(chains[[1]]) #Selects the first chain
    n_seq <- seq_len(n)
    log_tv <- log(truth_var)

    chain_dts <- lapply(chains, function(M) {
      run_sum  <- apply(M,  2, cumsum)
      run_sum2 <- apply(M^2, 2, cumsum)
      run_var  <- run_sum2 / n_seq - (run_sum / n_seq)^2
      log_vars <- log(run_var)
      truth_log <- matrix(rep(log_tv, each = n), nrow = n)
      sqrt(rowMeans((log_vars - truth_log)^2))
    })
    Reduce(`+`, chain_dts) / length(chain_dts)
  }

  mse_curve <- dt_curve <- NULL
  if (!is.null(truth_mean) && !is.null(truth_var)) {
    mse_curve <- mse_first_moment_avg_over_chains(chains, truth_mean)
    dt_curve  <- compute_dt_curve(chains, truth_var)
    esjd <- compute_esjd(chains)
  }

  list(
    mcse_mean = compute_mcse_mean,
    esjd_per_chain = esjd$per_chain_mean,
    esjd_running   = esjd$running_avg,
    curves = list(dt = dt_curve, mse = mse_curve, esjd = esjd$running_avg)
  )
}


