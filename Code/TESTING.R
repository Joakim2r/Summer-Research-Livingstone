################################################################################
################################# TESTING FUNCTION #############################
evaluate_samplers <- function(posterior_name, adapter = "adapter1", samplers = list("bimodal", "barker", "amala", "amh", 'hmc', "hmc_nuts"),
                              chains = 5, main_iter = 20000, warmup_iter = 2000,
                              kappa = 0.6, iteration_offset = 10,
                              gamma = 0.05,thin = 1, scale_init=NULL, plot=FALSE, n_step=10, sample_auxiliary=TRUE) {

  # Samplers to test (with their label for output)

  info <- posterior_draws(name = posterior_name)
  truth_mean <- info$truth_mean
  truth_var  <- info$truth_var
  d <- info$d


  # Storage for diagnostics
  basic_stats_list <- list()
  extra_stats_list <- list()
  chains_store     <- list()

  # Loop through each sampler
  for (sampler_name in samplers) {
    method <- sampler_name

    if (method == "hmc_nuts") {
      # HMC-NUTS (no adapters choice, uses built-in adaptation)
      cat("Running", method, "...\n")
      result <- run_NUTS(posterior_name, warm_up = warmup_iter, main = main_iter,
                         chains = chains, seed = 1234, d=d)

      # Convert to coda mcmc.list for diagnostics compatibility
      draws <- result$mcl
      chains <- result$chains
      basic_d <- basic_diagnostics(draws)
      basic_stats_list[[method]] <- basic_d
      chains_store[[method]] <- result$chains
      extra_stats_list[[method]] <- diagnostics(
        chains,
           ess = min(basic_d$ess_bulk, na.rm = TRUE),
           truth_mean = truth_mean, truth_var = truth_var
           )
      print(basic_stats_list[[method]])

    } else {
      # Other samplers – loop over the two adapter schemes
      cat("Running", method, "...\n")
      result <- run_MCMC(method = method,
                           name = posterior_name,
                           warm_up = warmup_iter, main = main_iter,
                           adapters = adapter, scale = scale_init,
                           chains = chains, kappa = kappa, iteration_offset = iteration_offset,
                           gamma = gamma,
                           thin = thin, n_step = n_step, sample_auxiliary = sample_auxiliary)

      basic_stats_list[[method]] <- result$diag_basic
      if (!is.null(result$diag_extra)) extra_stats_list[[method]] <- result$diag_extra
      if (!is.null(result$chains)) chains_store[[method]] <- result$chains
      }
  }

  ## Basic diagnostics table:
  # For brevity, we might summarize each method by min ESS and max Rhat across parameters.
  basic_summary <- do.call(rbind, lapply(names(basic_stats_list), function(key) {
    stats <- basic_stats_list[[key]]
    data.frame(
      Sampler = key,
      Rhat_max = max(stats$rhat, na.rm=TRUE),
      Rhat_min = max(stats$rhat, na.rm=TRUE),
      ESS_min  = min(stats$ess_bulk, na.rm=TRUE),
      ESS_tail_min = min(stats$ess_tail, na.rm=TRUE)
    )
  }))
  print(basic_summary)  # print or save as needed (use xtable for LaTeX)
  #print(extra_stats_list)
  #print(basic_stats_list)
  #print(chains_store)

  write.csv(basic_summary,
            file = file.path(getwd(), paste0(posterior_name, "_diagnostics.csv")),
            row.names = FALSE)

  mse_curves <- sapply(extra_stats_list, function(diag) diag$curves$mse)  # matrix of [iterations x methods]
  dt_curves <- sapply(extra_stats_list, function(diag) diag$curves$dt)
  esjd_curves <- sapply(extra_stats_list, function(diag) diag$curves$esjd)

  if (plot) {
    if (!is.null(mse_curves)) {
      pdf(paste0(posterior_name, "_runningMSE.pdf"))
      matplot(mse_curves, type = "l", lty = 1, xlab = "Iterations", ylab = "Running MSE of Mean",
              main = paste("Running MSE -", posterior_name))
      legend("topright", legend = colnames(mse_curves), lty = 1, col = seq_len(ncol(mse_curves)), cex = 0.6)
      dev.off()
    } else{
      print("MSE CURVES IS NULL")
    }

    if (!is.null(dt_curves)) {
      pdf(paste0(posterior_name, "_runningDT.pdf"))
      matplot(dt_curves, type = "l", lty = 1, xlab = "Iterations", ylab = "Running DT (log-Var distance)",
              main = paste("Running Variance Distance -", posterior_name))
      legend("topright", legend = colnames(dt_curves), lty = 1, col = seq_len(ncol(dt_curves)), cex = 0.6)
      dev.off()
    } else{
      print("DT CURVES IS NULL")
    }

    if (!is.null(esjd_curves)) {
      pdf(paste0(posterior_name, "_runningESJD.pdf"))
      matplot(esjd_curves, type = "l", lty = 1, xlab = "Iterations", ylab = "Running ESJD",
              main = paste("Running ESJD -", posterior_name))
      legend("topright", legend = colnames(esjd_curves), lty = 1, col = seq_len(ncol(esjd_curves)), cex = 0.6)
      dev.off()
    } else{
      print("ESJD CURVES IS NULL")
    }
  }
    # Return summary stats (invisibly) for further use if needed
    invisible(list(basic=basic_summary, basic_all=basic_stats_list, extra=extra_stats_list, mse_curves = mse_curves , dt_curves=dt_curves, esjd_curves=esjd_curves))
  }


eval <- evaluate_samplers(posterior_name = "eight_schools-eight_schools_noncentered",
                           chains = 2, main_iter = 2000, adapter = "adapter2", warmup_iter = 200, plot=TRUE)

################################################################################


