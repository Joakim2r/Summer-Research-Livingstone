evaluate_samplers <- function(posterior_name,
                              adapter = "adapter1",                   # NEW-1: can be c("adapter1","adapter2")
                              samplers = list("bimodal", "barker", "amala", "amh", "hmc", "hmc_nuts"),
                              chains = 5, main_iter = 20000, warmup_iter = 2000,
                              kappa = 0.6, iteration_offset = 10,
                              gamma = 0.05, thin = 1, scale_init = NULL, plot = FALSE,
                              n_step = 10,                             # NEW-2: can be c(10, NA) or c(10, 0) etc.
                              sample_auxiliary = TRUE,                 # NEW-3: can be c(TRUE, FALSE)
                              sample_n_step = NULL, seed = 1234) {

  info <- posterior_draws(name = posterior_name)
  truth_mean <- info$truth_mean
  truth_var  <- info$truth_var
  d          <- info$d

  basic_stats_list <- list()
  extra_stats_list <- list()
  chains_store     <- list()

  # Helper to add a row-bound curves matrix with padding
  bind_curves <- function(lst) {
    if (!length(lst)) return(NULL)
    L <- max(vapply(lst, length, 0L))
    M <- sapply(lst, function(v) { length(v) <- L; v }, simplify = "matrix")
    colnames(M) <- names(lst)
    M
  }

  # Loop over adapters (NEW-1)
  adapters_to_try <- as.character(adapter)
  if (length(adapters_to_try) == 0L) adapters_to_try <- "adapter1"

  for (adpt in adapters_to_try) {

    for (sampler_name in samplers) {
      method <- sampler_name

      if (method == "hmc_nuts") {
        # Run NUTS once per call (ignores adapter choice)
        if (identical(adpt, adapters_to_try[1])) {
          cat("Running", method, "...\n")
          result <- run_NUTS(posterior_name,
                             warm_up = warmup_iter, main = main_iter,
                             chains = chains, seed = seed, d = d)

          draws <- result$mcl
          chains_arr <- result$chains
          basic_d <- basic_diagnostics(draws)
          lab <- paste0("hmc_nuts")
          basic_stats_list[[lab]] <- basic_d
          chains_store[[lab]]     <- chains_arr
          extra_stats_list[[lab]] <- diagnostics(
            chains     = chains_arr,
            ess        = min(basic_d$ess_bulk, na.rm = TRUE),
            truth_mean = truth_mean, truth_var = truth_var
          )
          print(basic_stats_list[[lab]])
        }
        next
      }

      if (method == "hmc") {
        fixed_steps <- unique(as.integer(n_step[is.finite(n_step) & n_step > 0]))
        include_no_step <- any(!is.finite(n_step) | is.na(n_step) | n_step <= 0) || length(n_step) == 0L
        aux_flags <- if (length(sample_auxiliary)) as.logical(sample_auxiliary) else TRUE

        run_specs <- list()
        # fixed-step runs
        for (s in fixed_steps) {
          for (aux in aux_flags) {
            run_specs[[length(run_specs) + 1L]] <- list(
              n_step       = s,
              sample_aux   = aux,
              label_suffix = paste0("n", s, "_aux", aux)
            )
          }
        }
        # no-fixed-step runs (if requested via NA/NULL/≤0)
        if (include_no_step) {
          for (aux in aux_flags) {
            run_specs[[length(run_specs) + 1L]] <- list(
              n_step       = NULL,
              sample_aux   = aux,
              label_suffix = paste0("no_step_aux", aux)
            )
          }
        }
        if (!length(run_specs)) {
          run_specs <- list(list(n_step = NULL, sample_aux = TRUE, label_suffix = "no_step"))
        }

        # --- Execute the HMC variants ---
        for (spec in run_specs) {
          # If user didn’t supply a sampler and we’re on no-fixed-step, let HMC use its default
          snstep <- sample_n_step
          if (is.null(spec$n_step) && is.null(snstep)) snstep <- TRUE

          cat("Running hmc ... adapter =", adpt,
              ", n_step =", if (is.null(spec$n_step)) "NULL" else spec$n_step,
              ", sample_auxiliary =", spec$sample_aux, "\n")

          result <- run_MCMC(
            method = "hmc",
            name   = posterior_name,
            warm_up = warmup_iter, main = main_iter,
            adapters = adpt, scale = scale_init,
            chains = chains, kappa = kappa, iteration_offset = iteration_offset,
            gamma = gamma, thin = thin,
            n_step = spec$n_step,
            sample_n_step = snstep,
            sample_auxiliary = spec$sample_aux
          )

          lab <- paste0("hmc_", adpt, "_", spec$label_suffix)
          basic_stats_list[[lab]] <- result$diag_basic
          if (!is.null(result$diag_extra)) extra_stats_list[[lab]] <- result$diag_extra
          if (!is.null(result$chains))     chains_store[[lab]]     <- result$chains
          print(basic_stats_list[[lab]])
        }

      } else {
        # All other samplers: run once per adapter (NEW-1 in effect here)
        cat("Running", method, "... adapter =", adpt, "\n")
        result <- run_MCMC(
          method = method,
          name = posterior_name,
          warm_up = warmup_iter, main = main_iter,
          adapters = adpt, scale = scale_init,
          chains = chains, kappa = kappa, iteration_offset = iteration_offset,
          gamma = gamma, thin = thin,
          n_step = if (length(n_step)) n_step[1] else NULL,
          sample_n_step = sample_n_step,
          sample_auxiliary = if (length(sample_auxiliary)) sample_auxiliary[1] else TRUE
        )

        lab <- paste0(method, "_", adpt)
        basic_stats_list[[lab]] <- result$diag_basic
        if (!is.null(result$diag_extra)) extra_stats_list[[lab]] <- result$diag_extra
        if (!is.null(result$chains))     chains_store[[lab]]     <- result$chains
        print(basic_stats_list[[lab]])
      }
    }
  }

  # ------------ summary tables & plots (unchanged except a small fix) ------------
  basic_summary <- do.call(rbind, lapply(names(basic_stats_list), function(key) {
    stats <- basic_stats_list[[key]]
    data.frame(
      Sampler       = key,
      Rhat_max      = max(stats$rhat, na.rm = TRUE),
      Rhat_min      = min(stats$rhat, na.rm = TRUE),   # fixed (was max)
      ESS_min       = min(stats$ess_bulk, na.rm = TRUE),
      ESS_tail_min  = min(stats$ess_tail, na.rm = TRUE),
      row.names = NULL
    )
  }))
  print(basic_summary)

  utils::write.csv(
    basic_summary,
    file = file.path(getwd(), paste0(posterior_name, "_diagnostics.csv")),
    row.names = FALSE
  )

  mse_curves  <- bind_curves(lapply(extra_stats_list, function(diag) diag$curves$mse))
  dt_curves   <- bind_curves(lapply(extra_stats_list, function(diag) diag$curves$dt))
  esjd_curves <- bind_curves(lapply(extra_stats_list, function(diag) diag$curves$esjd))

  if (plot) {
    if (!is.null(mse_curves)) {
      pdf(paste0(posterior_name, "_runningMSE.pdf"))
      matplot(mse_curves, type = "l", lty = 1, xlab = "Iterations", ylab = "Running MSE of Mean",
              main = paste("Running MSE -", posterior_name))
      legend("topright", legend = colnames(mse_curves), lty = 1, col = seq_len(ncol(mse_curves)), cex = 0.6)
      dev.off()
    } else message("MSE CURVES IS NULL")

    if (!is.null(dt_curves)) {
      pdf(paste0(posterior_name, "_runningDT.pdf"))
      matplot(dt_curves, type = "l", lty = 1, xlab = "Iterations", ylab = "Running DT (log-Var distance)",
              main = paste("Running Variance Distance -", posterior_name))
      legend("topright", legend = colnames(dt_curves), lty = 1, col = seq_len(ncol(dt_curves)), cex = 0.6)
      dev.off()
    } else message("DT CURVES IS NULL")

    if (!is.null(esjd_curves)) {
      pdf(paste0(posterior_name, "_runningESJD.pdf"))
      matplot(esjd_curves, type = "l", lty = 1, xlab = "Iterations", ylab = "Running ESJD",
              main = paste("Running ESJD -", posterior_name))
      legend("topright", legend = colnames(esjd_curves), lty = 1, col = seq_len(ncol(esjd_curves)), cex = 0.6)
      dev.off()
    } else message("ESJD CURVES IS NULL")
  }

  invisible(list(basic = basic_summary,
                 basic_all = basic_stats_list,
                 extra = extra_stats_list,
                 mse_curves = mse_curves,
                 dt_curves = dt_curves,
                 esjd_curves = esjd_curves,
                 chains = chains_store))
}

evaluate_samplers(
  posterior_name = "eight_schools-eight_schools_noncentered",
  samplers = "hmc",
  n_step = c(3, 7, 15),
  chains = 2, main_iter = 2000, warmup_iter = 200
)

evaluate_samplers(
  posterior_name = "eight_schools-eight_schools_noncentered",
  samplers = "hmc",
  adapter = c("adapter1"),
  n_step = c(3, 7, 15),
  sample_auxiliary = c(TRUE, NULL),
  chains = 2, main_iter = 2000, warmup_iter = 200
)


evaluate_samplers(
  posterior_name = "eight_schools-eight_schools_noncentered",
  samplers = c("hmc"),
  n_step = c(5, 50),           # length > 1 triggers the second "no_step" run
  chains = 2, main_iter = 2000, warmup_iter = 200
)

evaluate_samplers(
  posterior_name = "eight_schools-eight_schools_noncentered",
  samplers = c("hmc"),
  sample_auxiliary = c(TRUE, FALSE),   # length > 1 triggers a second run
  n_step = 10,
  chains = 2, main_iter = 2000, warmup_iter = 200
)



