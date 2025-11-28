################################################################################
#################################### CHAIN FUNCTION ############################
sample_mcmc <- function(method = c("bimodal", "barker", "amala", "amh", 'hmc'),
                     name,
                     warm_up = 100,
                     main = 1000,
                     adapters = NULL,
                     scale = 0.8,
                     thin = 1,
                     chains = 2, kappa = 0.6,
                     iteration_offset=10,
                     gamma=0.05, tau_star=0.40,
                     n_step, sample_n_step=NULL, sample_auxiliary, plot = TRUE, save_plot=TRUE) {

  method <- match.arg(method)

  # --- model + truths once ---
  info   <- posterior_draws(name)
  bsmod  <- info$model
  target <- info$target
  d     <- info$d
  x0    <- info$x0
  shape <- diag(d)

  # --- scale
  default_scale <- switch(method,
                          amh     = 2.38/sqrt(d),   # RWM optimal sd
                          amala   = 0.8,
                          barker  = 0.6,
                          bimodal = 0.6,
                          hmc = 0.8
  )

  scale0 <- if (is.null(scale)) default_scale else scale


  # --- starting points around x0 ---
  starts <- lapply(seq_len(chains), function(k) x0 + rnorm(d, 0, 0.5) * (k - 1))

  # --- adapter presets if requested by name ---
  adapers <- NULL
  if (is.character(adapters)) {
    adapter1 <- list(
      stochastic_approximation_scale_adapter(kappa = kappa, target_accept_prob = tau_star),
      variance_shape_adapter(kappa = kappa)
    )
    adapter2 <- list(
      dual_averaging_scale_adapter(kappa = kappa, gamma = gamma,
                                   iteration_offset = iteration_offset,
                                   target_accept_prob = tau_star),
      variance_shape_adapter(kappa = kappa)
    )
    adapters <- match.arg(adapters, c("adapter1", "adapter2"))
    adapters <- if (adapters == "adapter1") adapter1 else adapter2
  } else {
    adapters <- NULL
  }

  # --- pick sampler
  sampler <- switch(
    method,
    bimodal = function(x)
      bimodal_barker(target = target, init_state = x,
             warmup_iter = warm_up, main_iter = main,
             scale = scale, shape = shape, adapters = adapters),
    barker = function(x)
      barker(target = target, init_state = x,
             warmup_iter = warm_up, main_iter = main,
             scale = scale, shape = shape, adapters = adapters),
    amala  = function(x)
      AMALA(target = target, init_state = x,
             warmup_iter = warm_up, main_iter = main,
             adapters = adapters, scale=scale, shape=shape),
    amh  = function(x)
      AdaptiveMH(target = target, init_state = x,
             warmup_iter = warm_up, main_iter = main,
             scale = scale, shape = shape, adapters = adapters),
    hmc = function(x)
      HMC(target = target, init_state = x,
          warmup_iter = warm_up, main_iter = main,
          scale = scale, shape = shape, n_step = n_step, sample_n_step = sample_n_step,
          sample_auxiliary = sample_auxiliary, adapters=adapters)
  )

  # --- run chains & convert to mcmc ---
  chains_raw <- lapply(starts, function(x) sampler(x)$traces)

  # 2) Keep only the first d parameters, seq_len(d) is the equivalent of 0:d
  chains_mat <- lapply(chains_raw, function(ch) as.matrix(ch)[, seq_len(d), drop = FALSE])

  # 3) Convert to coda objects and build the mcmc.list
  mcl <- mcmc.list(lapply(chains_mat, function(m) mcmc(m, start = 1)))

  # 4) Drop warmup
  mcl <- window(mcl, start = warm_up + 1)


  # --- basic diag results ---
  diag_basic <- basic_diagnostics(mcl)
  print(diag_basic)

  #   # --- advanced diagnostics (only if truths available) ---
  diag_extra <- NULL
  if (!is.null(info$truth_mean) && !is.null(info$truth_var)) {
    diag_extra <- diagnostics(
    chains = mcl,
    ess = diag_basic$ess_bulk,
    truth_mean     = info$truth_mean,
    truth_var = info$truth_var
    )}

  mse_curve  <- diag_extra$curves$mse
  dt_curve   <- diag_extra$curves$dt
  esjd_curve <- diag_extra$curves$esjd

  # small draw util that skips legend when only one series
  graph <- function(mat, main, ylab) {
    if (is.null(mat)) return(invisible())
    matplot(mat, type = "l", lty = 1, xlab = "Iterations", ylab = ylab,
            main = main)
  }

  if (isTRUE(plot)) {
    if (!is.null(mse_curve))  {
      graph(mse_curve,  paste("Running MSE -", name),  "Running MSE of Mean")
      } else message("MSE CURVES IS NULL")
    if (!is.null(dt_curve)) {
      graph(dt_curve,   paste("Running Variance Distance -", name), "Running DT (log-Var distance)")
    } else message("DT CURVES IS NULL")
    if (!is.null(esjd_curve)) {
      graph(esjd_curve, paste("Running ESJD -", name), "Running ESJD")
    } else message("ESJD CURVES IS NULL")
  }

  if (isTRUE(save_plot)) {
    if (!is.null(mse_curve)) {
      pdf(paste0(name, "_runningMSE.pdf"))
      graph(mse_curve, paste("Running MSE -", name), "Running MSE of Mean")
      dev.off()
    } else message("MSE CURVES IS NULL")

    if (!is.null(dt_curve)) {
      pdf(paste0(name, "_runningDT.pdf"))
      graph(dt_curve, paste("Running Variance Distance -", name), "Running DT (log-Var distance)")
      dev.off()
    } else message("DT CURVES IS NULL")

    if (!is.null(esjd_curve)) {
      pdf(paste0(name, "_runningESJD.pdf"))
      graph(esjd_curve, paste("Running ESJD -", name), "Running ESJD")
      dev.off()
    } else message("ESJD CURVES IS NULL")
  }



  list(mcl= mcl, diag_basic = diag_basic, diag_extra=diag_extra)
}


