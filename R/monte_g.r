#' Collect mediator pool from a reference simulation
#'
#' @description
#' Runs the standard parametric g-formula forward simulation under
#' \code{intervention_ref} (the reference exposure, a*) and collects the
#' simulated mediator values at each time step. The returned pool is used by
#' the Phi10 arm for interventional mediation (Lin et al. 2017) to
#' implement the two-trajectory identification formula (Eq. 4): the mediator
#' distribution must be conditioned on l' (confounders evolving under a*), not
#' on the l trajectory used in the outcome model.
#'
#' For survival analyses, each time-step pool element is a list with
#' \code{vals} (mediator values) and \code{weights} (the cumulative survival
#' \eqn{Sc} under a* for each individual). The Phi10 arm then samples from
#' \code{vals} with probability proportional to \code{weights}, implementing
#' the \eqn{S(1{:}t-1)=1} condition from Eq. 4.
#' For non-survival analyses, \code{weights} is \code{NULL} (uniform sampling).
#'
#' @keywords internal
.collect_med_pool <- function(data, models, exposure, time_var, time_seq,
                               intervention_ref,
                               init_recode = NULL, in_recode = NULL,
                               out_recode = NULL) {

  med_flag <- sapply(models, function(m) m$mod_type == "mediator")
  if (!any(med_flag)) return(NULL)
  med_idx  <- which(med_flag)
  med_var  <- models[[med_idx]]$rsp_vars
  med_mod  <- models[[med_idx]]

  time_seq   <- sort(time_seq)
  min_time   <- min(time_seq)
  interv_vec <- rep(intervention_ref, length(time_seq))
  m_pool     <- list()

  for (indx in seq_along(time_seq)) {
    t_index <- time_seq[indx]
    set(data, j = time_var, value = t_index)

    if (t_index == min_time && !is.null(init_recode))
      apply_recodes(data, init_recode)
    if (!is.null(in_recode) && t_index != min_time)
      apply_recodes(data, in_recode)

    # Standard g-formula pass under a*: generates the l' trajectory and M values.
    # mediation_type = NA so no cross-world mediator logic is applied here.
    data <- simulate_data(
      data           = data,
      exposure       = exposure,
      models         = models,
      intervention   = interv_vec[indx],
      mediation_type = NA
    )

    # Collect mediator values for eligible individuals.
    # Start with the user-specified subset for the mediator model (if any).
    if (!is.null(med_mod$subset)) {
      cond <- eval(med_mod$subset, envir = data, enclos = parent.frame())
      cond <- cond & !is.na(cond)
    } else {
      cond <- rep(TRUE, nrow(data))
    }
    # For survival analyses, weight the pool by Sc (survival probability under a*).
    # This implements the S(1:t-1)=1 conditioning in Lin et al. (2017, Eq. 4):
    # individuals with lower cumulative survival are drawn with proportionally
    # lower probability when the Phi10 arm samples from this pool.
    surv_present <- any(sapply(models, function(m) m$mod_type == "survival"))
    if (surv_present && "Sc" %in% names(data)) {
      m_pool[[as.character(t_index)]] <- list(
        vals    = data[[med_var]][cond],
        weights = data$Sc[cond]
      )
    } else {
      m_pool[[as.character(t_index)]] <- list(
        vals    = data[[med_var]][cond],
        weights = NULL
      )
    }

    if (!is.null(out_recode) && t_index != min_time)
      apply_recodes(data, out_recode)
  }

  return(m_pool)
}


#' Main calculation function
#'
#' This function will receive the parameters and fit model. After the model is fitted,
#' random samples will be drawn from the data and apply the intervention.
#'
#' @inheritParams gformula
#' @param mediation_type Type of the mediation analysis, if the value is \code{NA}
#' no mediation analysis will be performed (default). It will be ignored if the intervention
#'  is not \code{NULL}
#' @param return_fitted Return the fitted model (default is FALSE).
#' @keywords internal

.gformula <- function(data,
                      id_var,
                      base_vars,
                      time_var,
                      exposure,
                      models,
                      intervention,
                      in_recode = NULL,
                      out_recode = NULL,
                      init_recode = NULL,
                      mediation_type = c(NA, "N", "I"),
                      mc_sample = 10000,
                      return_fitted = FALSE,
                      return_data = FALSE,
                      seed = mc_sample * 100) {
  if (length(mediation_type) > 1) {
    mediation_type <- mediation_type[1]
  } else if (!is.na(mediation_type) && !(mediation_type %in% c("N", "I"))) {
    stop("'mediation_type' must be NA, \"N\", or \"I\"")
  }

  fit_mods <- lapply(models, function(mods) {
    rsp_vars <- all.vars(formula(mods$call)[[2]])

    # Observed values range
    if (is.numeric(data[[rsp_vars]])) {
      val_ran <- range(na.omit(data[[rsp_vars]]))
    } else {
      val_ran <- unique(na.omit(data[[rsp_vars]]))
    }

    # Recode data before simulating
    if (!is.null(mods$recode)) {
      for (i in seq_along(mods$recode)) {
        data <- within(data, eval(parse(text = mods$recode[i])))
      }
    }

    mods$call$data <- substitute(data, env = parent.frame())

    fitmodel <- run_withwarning_collect(
      eval(mods$call),
      msg = sprintf("Outcome model: %s", rsp_vars)
    )

    # Pre-extract prediction components so the hot-path simulation loop can
    # skip predict() overhead (model.frame construction + na.action).
    # model.matrix(Xterms, newdt) %*% beta is a direct BLAS call.
    list(
      fitted     = fitmodel,
      Xterms     = delete.response(terms(fitmodel)),
      beta       = coef(fitmodel),
      linkinv    = if (inherits(fitmodel, "glm")) family(fitmodel)$linkinv
                   else identity,
      sigma      = tryCatch(sigma(fitmodel), error = function(e) NULL),
      recodes    = mods$recode,
      subset     = mods$subset,
      var_type   = mods$var_type,
      mod_type   = mods$mod_type,
      custom_sim = mods$custom_sim,
      rsp_vars   = rsp_vars,
      val_ran    = val_ran
    )
  })

  # only call set.seed() when a seed is explicitly supplied.
  # Omitting the seed inside bootstrap replicates lets each replicate draw
  # different Monte Carlo samples, giving valid bootstrap variance estimates.
  if (!is.null(seed)) set.seed(seed)

  # Baseline variables for Monte Carlo Sampling
  base_dat <- unique(data[, unique(c(id_var, base_vars)), with = FALSE])
  df_mc    <- data.table::as.data.table(
    base_dat[sample(1:length(base_dat[[id_var]]), mc_sample, replace = TRUE), ]
  )
  df_mc[, new_ID := seq_len(.N)]

  # Cache the time sequence once (used in each arm and in pool collection)
  time_seq <- sort(unique(data[[time_var]]))

  # Derive the reference exposure value (a*) from the non-NULL arms.
  # Used both for the interventional mediator pool pass and for the natural-effects
  # exposure swap inside simulate_data (passed as med_ref_val).
  # isTRUE() guards against mediation_type == NA returning NA instead of FALSE.
  if (!is.na(mediation_type) && any(sapply(intervention, is.null))) {
    non_null_vals    <- unlist(Filter(Negate(is.null), intervention))
    intervention_ref <- min(non_null_vals)
  } else {
    intervention_ref <- 0L
  }

  # For interventional mediation, run the reference arm (Phi00, a*=0) first
  # with collect_pool = TRUE. This simultaneously produces the Phi00 risk
  # estimate AND the mediator pool, saving a full simulation pass vs. calling
  # .collect_med_pool separately (Lin et al. 2017, Eq. 4).
  if (isTRUE(mediation_type == "I") && any(sapply(intervention, is.null))) {
    ref_arm_data <- data.table::copy(df_mc)
    ref_run <- monte_g(
      data           = ref_arm_data,
      models         = fit_mods,
      exposure       = exposure,
      time_var       = time_var,
      time_seq       = time_seq,
      intervention   = rep(intervention_ref, length(time_seq)),
      init_recode    = init_recode,
      in_recode      = in_recode,
      out_recode     = out_recode,
      mediation_type = mediation_type,
      return_data    = return_data,
      med_pool       = NULL,
      med_ref_val    = intervention_ref,
      collect_pool   = TRUE
    )
    m_pool      <- ref_run$pool
    cached_ref  <- ref_run$estimate
  } else {
    m_pool      <- NULL
    cached_ref  <- NULL
  }

  res <- sapply(intervention, function(i) {
    # The reference arm (Phi00) was already run above for pool collection;
    # return its cached estimate directly to avoid a redundant simulation.
    if (!is.null(cached_ref) && !is.null(i) && isTRUE(all(i == intervention_ref))) {
      return(cached_ref)
    }

    # copy df_mc so each arm starts with a clean Sc/Pred_Y state.
    # Without this, survival accumulation from the previous arm bleeds through.
    arm_data <- data.table::copy(df_mc)

    # pass the pre-collected mediator pool to the Phi10 arm only.
    # Other arms (Ph11, Phi00) simulate the mediator from their own trajectory.
    arm_med_pool <- if (is.null(i) && isTRUE(mediation_type == "I")) m_pool else NULL

    r <- monte_g(
      data           = arm_data,
      models         = fit_mods,
      exposure       = exposure,
      time_var       = time_var,
      time_seq       = time_seq,
      intervention   = i,
      init_recode    = init_recode,
      in_recode      = in_recode,
      out_recode     = out_recode,
      mediation_type = mediation_type,
      return_data    = return_data,
      med_pool       = arm_med_pool,
      med_ref_val    = intervention_ref
    )

    return(r)
  }, simplify = FALSE)

  if (return_fitted) {
    return(list(fitted.models = fit_mods, gform.data = res))
  } else {
    return(list(gform.data = res))
  }
}

#' Monte Carlo simulation
#'
#' @description
#'  Internal use only. Monte Carlo simulation.
#'
#' @inheritParams gformula
#' @param time_seq Time sequence vector of the data.
#' @param med_pool Named list of pre-simulated mediator values (from \code{.collect_med_pool}),
#'   keyed by time index as character. Used by the Phi10 arm for interventional mediation.
#' @param med_ref_val The reference exposure value (a*) passed through to
#'   \code{simulate_data} for the natural-effects mediator swap. Defaults to \code{0L}.
#' @param collect_pool Logical. If \code{TRUE}, collect the mediator pool during the
#'   forward simulation and return it alongside the risk estimate. The returned value is
#'   then a \code{list(estimate = ..., pool = ...)} rather than a plain scalar/data.table.
#'   Used by \code{.gformula} to avoid a separate pool-collection pass for interventional
#'   mediation: the reference arm (Phi00) doubles as the pool source.
#'
#' @keywords internal
#'
monte_g <- function(data,
                    models,
                    exposure,
                    time_var,
                    time_seq,
                    intervention = NULL,
                    init_recode  = NULL,
                    in_recode    = NULL,
                    out_recode   = NULL,
                    mediation_type = c(NA, "N", "I"),
                    return_data  = FALSE,
                    med_pool     = NULL,
                    med_ref_val  = 0L,
                    collect_pool = FALSE) {

  if (length(mediation_type) > 1) {
    mediation_type <- mediation_type[1]
  } else if (!is.na(mediation_type) && !(mediation_type %in% c("N", "I"))) {
    stop("'mediation_type' must be NA, \"N\", or \"I\"")
  }

  # Replicate static interventions to the same length as the time sequence.
  # dyn_int() objects are not replicated — the same rule applies at every time step.
  time_len <- length(time_seq)
  if (length(intervention) == 1 && !inherits(intervention, "causalMed_dynint")) {
    intervention <- rep(intervention, time_len)
  }

  # Pre-detect mediator model for optional pool collection.
  # Done once here so the time loop avoids repeated sapply() over models.
  if (collect_pool) {
    med_flag_pc <- sapply(models, function(m) m$mod_type == "mediator")
    if (!any(med_flag_pc)) {
      collect_pool <- FALSE   # no mediator model — nothing to collect
    } else {
      med_idx_pc   <- which(med_flag_pc)
      med_var_pc   <- models[[med_idx_pc]]$rsp_vars
      med_mod_pc   <- models[[med_idx_pc]]
      surv_flag_pc <- any(sapply(models, function(m) m$mod_type == "survival"))
      m_pool_out   <- list()
    }
  }

  # Get the position of the outcome
  out_flag <- sapply(models, function(mods) mods$mod_type %in% c("outcome", "survival"))
  out_flag <- which(out_flag)

  # Get the position of the censor
  cen_flag  <- sapply(models, function(mods) mods$mod_type == "censor")
  surv_flag <- sapply(models, function(mods) mods$mod_type == "survival")
  if (any(cen_flag) | any(surv_flag)) {
    is_survival <- TRUE
  } else {
    is_survival <- FALSE
  }

  if (any(cen_flag))
    cen_flag <- which(cen_flag)
  else
    cen_flag <- 0

  # Get the variable name of outcome and censor
  outcome <- all.vars(formula(models[[out_flag]]$fitted)[[2]])
  censor  <- NULL
  if (!all(cen_flag == 0)) {
    censor <- all.vars(formula(models[[cen_flag]]$fitted)[[2]])
  }

  # Get minimum and maximum time
  max_time <- max(time_seq, na.rm = TRUE)
  min_time <- min(time_seq, na.rm = TRUE)

  # Run g-formula
  for (indx in seq_along(time_seq)) {
    t_index <- time_seq[indx]

    set(data, j = time_var, value = t_index)

    # Re-code baseline variables at initiation
    if (t_index == min_time) {
      if (!is.null(init_recode)) {
        apply_recodes(data, init_recode)
      }
    }

    # Re-code data before simulating
    if (!is.null(in_recode) & t_index != min_time) {
      apply_recodes(data, in_recode)
    }

    # pass the time-specific mediator pool slice to simulate_data.
    t_med_pool <- if (!is.null(med_pool)) med_pool[[as.character(t_index)]] else NULL

    # Use the model to calculate the simulated value.
    # dyn_int() objects are passed as-is (same rule at every time step);
    # static interventions are indexed to get the value for this time step.
    current_int <- if (inherits(intervention, "causalMed_dynint")) intervention else intervention[indx]
    data <- simulate_data(
      data           = data,
      exposure       = exposure,
      models         = models,
      intervention   = current_int,
      mediation_type = mediation_type,
      med_pool       = t_med_pool,
      med_ref_val    = med_ref_val
    )

    # Collect mediator pool if requested (used by the Phi00 reference arm).
    if (collect_pool) {
      cond_pc <- if (!is.null(med_mod_pc$subset)) {
        cond_tmp <- eval(med_mod_pc$subset, envir = data, enclos = parent.frame())
        cond_tmp & !is.na(cond_tmp)
      } else {
        rep(TRUE, nrow(data))
      }
      if (surv_flag_pc && "Sc" %in% names(data)) {
        m_pool_out[[as.character(t_index)]] <- list(
          vals    = data[[med_var_pc]][cond_pc],
          weights = data$Sc[cond_pc]
        )
      } else {
        m_pool_out[[as.character(t_index)]] <- list(
          vals    = data[[med_var_pc]][cond_pc],
          weights = NULL
        )
      }
    }

    # For survival outcome: under intervention, disable censoring.
    # All individuals remain in the pool; risk is computed analytically
    # from Pred_Y = 1 - prod(1 - h_t) accumulated in simulate_data(),
    # matching the gfoRmula reference package's poprisk approach.
    if (is_survival) {
      # If the intervention is defined, no censoring applied.
      if (!is.null(intervention) & cen_flag != 0) {
        set(data, j = censor, value = 0)
      }

      # If censored outcome equals to 1, then the survival outcome is not observed.
      if (cen_flag != 0) {
        data[[outcome]] <- ifelse(data[[censor]] == 1, 0, data[[outcome]])
      }
    }

    # Recode data after simulating
    if (!is.null(out_recode) & t_index != min_time) {
      apply_recodes(data, out_recode)
    }
  }
  # loop ends here

  result <- if (return_data) {
    data
  } else {
    # Mean of analytic Pred_Y over all N individuals (no row removal).
    # For survival: Pred_Y = 1 - prod(1-h_t), the cumulative risk.
    # For non-survival: Pred_Y = predicted outcome at final time.
    sum(data[["Pred_Y"]]) / length(data[["Pred_Y"]])
  }

  if (collect_pool) return(list(estimate = result, pool = m_pool_out))
  return(result)
}
