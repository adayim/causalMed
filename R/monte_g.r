#' Main calculation function
#' 
#' This function will receive the parameters and fit model. After the model is fitted,
#' random samples will be drawn from the data and apply the intervention. 
#'
#' @inheritParams gformula
#' @param mediation_type Type of the mediation effect.
#' @param return_fitted Return the fitted model (default is FALSE).
#'
#' @importFrom pbapply pbsapply
#' 
#' @keywords internal

.gformula <- function(data,
                      id_var,
                      base_vars,
                      time_var,
                      exposure,
                      models,
                      intervention,
                      in_recode,
                      out_recode,
                      init_recode,
                      mediation_type = c(NA, "N", "I"),
                      mc_sample = 10000,
                      return_fitted = FALSE) {
  mediation_type <- match.arg(mediation_type)

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

    list(
      fitted = eval(mods$call),
      recodes = mods$recode,
      subset = mods$subset,
      var_type = mods$var_type,
      mod_type = mods$mod_type,
      custom_sim = mods$custom_sim,
      rsp_vars = rsp_vars,
      val_ran = val_ran
    )
  })

  # Setting seeds
  set.seed(12345 * mc_sample)

  # Baseline variables for Monte Carlo Sampling
  base_dat <- unique(data[, unique(c(id_var, base_vars))])
  df_mc <- data.table::as.data.table(base_dat[sample(1:nrow(base_dat), mc_sample, replace = TRUE), ])
  df_mc[, new_ID := seq_len(.N)]

  if (!is.na(mediation_type)) {
    intervention <- list(always = 1, never = 0, mediation = NULL)
  }

  #
  cat("\nLooping interventions:\n")

  res <- pbapply::pbsapply(intervention, function(i) {
    monte_g(
      data = df_mc,
      time_var = time_var,
      time_seq = unique(data[[time_var]]),
      exposure = exposure,
      models = fit_mods,
      intervention = i,
      in_recode = in_recode,
      out_recode = out_recode,
      init_recode = init_recode,
      mediation_type = mediation_type
    )
  }, simplify = FALSE)

  if (return_fitted) {
    return(list(fitted.models = fit_mods, gform.data = res))
  } else {
    return(res)
  }
}

#' Monte Carlo simulation
#'
#' @description
#'  Internal use only. Monte Carlo simulation.
#'
#' @inheritParams gformula
#' @param time_seq Time sequence vector of the data.
#' @param mediation_type Type of the mediation analysis, if the value is \code{NA} 
#' no mediation analysis will be performed (default). It will be ignored if the intervention
#'  is not \code{NULL}
#' 
#' @importFrom pbapply timerProgressBar setTimerProgressBar
#'
#' @keywords internal
#'
monte_g <- function(data,
                    exposure,
                    time_seq,
                    time_var,
                    models,
                    intervention = NULL,
                    in_recode = NULL,
                    out_recode = NULL,
                    init_recode = NULL,
                    mediation_type = c(NA, "N", "I")) {
  mediation_type <- match.arg(mediation_type)

  # Replicate to intervention to the same length of the time.
  time_len <- length(time_seq)
  if (length(intervention) == 1) {
    intervention <- rep(intervention, time_len)
  }

  # Get the position of the outcome
  out_flag <- sapply(models, function(mods) mods$mod_type %in% c("outcome", "survival"))
  out_flag <- which(out_flag)

  # Get the position of the censor
  cen_flag <- sapply(models, function(mods) mods$mod_type == "censor")
  surv_flag <- sapply(models, function(mods) mods$mod_type == "survival")
  if (any(cen_flag) | any(surv_flag)) {
    is_survival <- TRUE
  }
  cen_flag <- which(cen_flag)

  # Get the variable name of outcome and censor
  outcome <- all.vars(formula(models[[out_flag]]$fitted)[[2]])
  if (!all(cen_flag == 0)) {
    censor <- all.vars(formula(models[[cen_flag]]$fitted)[[2]])
  }

  # Simulate
  max_time <- max(time_seq, na.rm = TRUE)
  min_time <- min(time_seq, na.rm = TRUE)
  dat_y <- data

  # Run normal g-formula
  for (t_index in sort(time_seq)) {
    dat_y[[time_var]] <- t_index

    # Recode baseline variables at initiation
    if (t_index == min_time) {
      if (!is.null(init_recode)) {
        dat_y <- within(dat_y, eval(parse(text = init_recode)))
      }
    }

    # Recode data before simulating
    if (!is.null(in_recode)) {
      dat_y <- within(dat_y, eval(parse(text = in_recode)))
    }

    # Use the model to calculate the simulated value
    dat_y <- simulate_data(data = dat_y, exposure = exposure, models = models, intervention = intervention[t_index], mediation_type = mediation_type)

    if (!is.null(out_recode)) {
      dat_y <- within(dat_y, eval(parse(text = out_recode)))
    }

    # For survival outcome
    if (is_survival) {
      # If the intervention is defined, no censoring applied.
      if (!is.null(intervention) & cen_flag != 0) {
        dat_y[[censor]] <- 0
      }

      # If censored outcome equals to 1, then the survival outcome is not observed.
      if (cen_flag != 0) {
        dat_y[[outcome]] <- ifelse(dat_y[[censor]] == 1, 0, dat_y[[outcome]])
      }

      # Output data
      if (t_index == min_time) {
        out_y <- dat_y[, c("new_ID", time_var, outcome, "Pred_Y"), with = FALSE]
      } else {
        out_y <- rbind(out_y, dat_y[, c("new_ID", time_var, outcome, "Pred_Y"), with = FALSE])
      }

      # Continue for the not censored or no event observed
      if (cen_flag != 0) {
        dat_y <- dat_y[!(dat_y[[outcome]] == 1 | dat_y[[censor]] == 1), ]
      } else {
        dat_y <- dat_y[dat_y[[outcome]] != 1, ]
      }

      # If no data left
      if (nrow(dat_y) == 0) {
        break
      }

      # If all loop complete
      if (t_index == max_time) {
        out_y <- rbind(out_y, dat_y[, c("new_ID", time_var, outcome, "Pred_Y"), with = FALSE])
      }
    } else {
      # Output data
      if (t_index == min_time) {
        out_y <- dat_y[, c("new_ID", time_var, outcome, "Pred_Y"), with = FALSE]
      } else {
        out_y <- rbind(out_y, dat_y[, c("new_ID", time_var, outcome, "Pred_Y"), with = FALSE])
      }
    }
  }
  # loop ends here

  return(out_y)
}
