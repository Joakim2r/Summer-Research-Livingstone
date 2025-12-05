############################## SAMPLERS #############################

############################## BIMODAL BARKER FUNCTION #############################
bimodal_barker <- function(
    target,
    init_state,
    warmup_iter,
    main_iter,
    scale,
    shape,
    adapters,
    sigma) {

  prop <- bimodal_barker_proposal(
    sigma = sigma,
    scale = scale,
    shape = shape
  )

  fit_bimod_barker <- sample_chain(
    target_distribution = target,
    initial_state       = init_state,
    n_warm_up_iteration = warmup_iter,
    n_main_iteration    = main_iter,
    proposal            = prop,
    adapters            = adapters,
    show_progress       = TRUE
  )

  return(fit_bimod_barker)
}


############################## BARKER FUNCTION #############################
barker<- function(target,
                  init_state,
                  warmup_iter,
                  main_iter,
                  scale,
                  tau_star,
                  shape,
                  adapters) {

  prop <- barker_proposal(
    scale = scale,
    shape = shape
  )

  fit_barker <- sample_chain(
    target_distribution = target,
    initial_state       = init_state,
    n_warm_up_iteration = warmup_iter,
    n_main_iteration    = main_iter,
    proposal            = prop,
    adapters            = adapters,
    show_progress       = TRUE
  )

  return(fit_barker)
}


################################ AMALA ################################
AMALA <- function(
    target,
    init_state,
    warmup_iter,
    main_iter,
    adapters,  # or a custom list
    tau_star,
    kappa,
    gamma,
    scale,
    shape,
    iteration_offset
) {

  fit_amala <- sample_chain(
    target_distribution = target,
    initial_state       = init_state,
    n_warm_up_iteration = warmup_iter,
    n_main_iteration    = main_iter,
    proposal            = langevin_proposal(scale=scale, shape=shape),
    adapters            = adapters,
    show_progress       = TRUE
  )

  return(fit_amala)
}
############################## ADAPTIVE MH #####################################
AdaptiveMH <- function(target,
                       init_state,
                       warmup_iter,
                       main_iter,
                       shape,
                       adapters,
                       tau_star,
                       scale,
                       show_progress) {

  prop <- random_walk_proposal(scale = scale, shape = shape)

  fit <- sample_chain(
    target_distribution = target,
    initial_state       = init_state,
    n_warm_up_iteration = warmup_iter,
    n_main_iteration    = main_iter,
    proposal            = prop,
    adapters            = adapters,
    show_progress       = TRUE
  )

  return(fit)
}
############################## HMC #############################
HMC <- function(target,
                init_state,
                warmup_iter,
                main_iter,
                scale,
                shape,
                adapters,
                n_step,
                sample_n_step,
                sample_auxiliary) {


  partial_momentum_update <- function(state, phi = pi / 4) {
    momentum <- state$momentum()
    if (is.null(momentum)) {
      stats::rnorm(state$dimension())
    } else {
      cos(phi) * momentum + sin(phi) * stats::rnorm(length(momentum))
    }
  }

  aux_fun <- NULL
  if (isTRUE(sample_auxiliary)) {
    aux_fun <- partial_momentum_update
  } else if (isFALSE(sample_auxiliary) || is.null(sample_auxiliary)) {
    aux_fun <- NULL
  } else if (is.function(sample_auxiliary)) {
    aux_fun <- sample_auxiliary
  } else {
    stop("`sample_auxiliary` must be TRUE, FALSE/NULL, or a function(state, ...).")
  }

  adaptive_sample_state <- function(min = 5, max = 100) {
    # take first element if vectors / handle length-0, NA, non-finite
    mn <- suppressWarnings(as.numeric(min[1]))
    mx <- suppressWarnings(as.numeric(max[1]))

    if (length(mn) == 0L || !is.finite(mn)) mn <- 5
    if (length(mx) == 0L || !is.finite(mx)) mx <- 100

    # order & coerce to integers
    if (mx < mn) {
      tmp <- mn
      mn <- mx
      mx <- tmp
    }

    mn <- as.integer(round(mn))
    mx <- as.integer(round(mx))

    # ensure at least one admissible value
    n <- mx - mn + 1L
    if (n < 1L) {
      n <- 1L
      mx <- mn
      }

    mn + sample.int(n, 1L) - 1L
  }


  if (!is.null(n_step) && is.finite(n_step) && n_step > 0) {
    step_count <- as.integer(n_step)
    step_sampler <- NULL
  } else {
    step_count <- NULL
    if (isTRUE(sample_n_step) || is.null(sample_n_step)) {
      # default random step count in [5, 100]
      step_sampler <- adaptive_sample_state
    } else if (is.function(sample_n_step)) {
      step_sampler <- sample_n_step
    } else {
      stop("`sample_n_step` must be NULL or a function(min, max).")
    }
  }

  if (is.null(aux_fun)) {                                           # [FIX] Use default momentum if aux_fun is NULL
    prop <- hamiltonian_proposal(
      scale         = scale,
      shape         = shape,
      n_step        = step_count,
      sample_n_step = step_sampler
      # sample_auxiliary omitted: uses default (standard Gaussian momentum)       [FIX]
    )
  } else {
    prop <- hamiltonian_proposal(
      scale            = scale,
      shape            = shape,
      n_step           = step_count,
      sample_n_step    = step_sampler,
      sample_auxiliary = aux_fun
    )
  }

  # -------- run chain --------
  fit <- sample_chain(
    target_distribution = target,
    initial_state       = init_state,
    n_warm_up_iteration = warmup_iter,
    n_main_iteration    = main_iter,
    proposal            = prop,
    adapters            = adapters,
    show_progress       = TRUE
  )

  return(fit)
}


############################## HMC NUTS FUNCTION #############################
run_NUTS <- function(name,
                     init,
                     warm_up,
                     main,
                     adapt_delta,
                     max_treedepth,
                     chains,
                     seed,
                     thin) {

  p <- posterior(name, pdb)
  dfp <- data_file_path(p)
  cmdstan_model <- cmdstan_model(stan_code_file_path(p), quiet=TRUE)

  # 3. run Stan (NUTS is default)
  fit <- cmdstan_model$sample(
    data = dfp,
    init = init,
    chains  = chains,
    parallel_chains = chains/2,
    iter_warmup = warm_up,
    iter_sampling = main,
    adapt_delta = adapt_delta,
    max_treedepth = max_treedepth,
    seed = seed,
    thin=thin,
    show_messages = TRUE,
    save_warmup = FALSE,
    refresh= 0
  )


  mcl <- as_mcmc.list(fit)
  chains_mat <- lapply(mcl, function(ch) as.matrix(ch)[, seq_len(d), drop = FALSE])
  mcl <- mcmc.list(lapply(chains_mat, mcmc))

  list(mcl=mcl, chains=chains_mat)
}

