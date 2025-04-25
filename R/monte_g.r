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
#' @param progress_bar Show progress bar (default).
#'
#' @importFrom progressr handlers handler_progress progressor
#'
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
                      progress_bar = TRUE,
                      seed=mc_sample*100) {
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
  set.seed(seed)

  # Baseline variables for Monte Carlo Sampling
  base_dat <- unique(data[, unique(c(id_var, base_vars)), with = FALSE])
  df_mc <- data.table::as.data.table(base_dat[sample(1:length(base_dat[[id_var]]), mc_sample, replace = TRUE), ])
  df_mc[, new_ID := seq_len(.N)]

  if (!is.na(mediation_type)) {
    intervention <- list(always = 1, never = 0, mediation = NULL)
  }

  # Setup progress bar
  if (progress_bar) {
    progressr::handlers(global = TRUE)
    progressr::handlers(list(
      progressr::handler_progress(
        format   = ":spin :current/:total (:message) [:bar] :percent in :elapsedfull ETA: :eta",
        # width    = 60,
        complete = "+"
      )
    ))
    p <- progressr::progressor(steps = length(intervention))
  } else {
    p <- NULL
  }

  res <- sapply(intervention, function(i) {
    if (progress_bar) {
      p(message = "Running intervention", amount = 0)
    }

    r <- monte_g(
      data = df_mc,
      models = fit_mods,
      exposure = exposure,
      time_var = time_var,
      time_seq = sort(unique(data[[time_var]])),
      intervention = i,
      init_recode = init_recode,
      in_recode = in_recode,
      out_recode = out_recode,
      mediation_type = mediation_type,
      return_data = return_data,
      progress_bar = p
    )

    if (progress_bar) {
      p(message = "Running intervention")
    }

    return(r)
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
#' @param progress_bar Progress bar object.
#'
#'
#' @keywords internal
#'
monte_g <- function(data,
                    models,
                    exposure,
                    time_var,
                    time_seq,
                    intervention = NULL,
                    init_recode = NULL,
                    in_recode = NULL,
                    out_recode = NULL,
                    mediation_type = c(NA, "N", "I"),
                    return_data = FALSE,
                    progress_bar = NULL) {

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
  }else {
    is_survival <- FALSE
  }

  if(any(cen_flag))
    cen_flag <- which(cen_flag)
  else
    cen_flag <- 0

  # Get the variable name of outcome and censor
  outcome <- all.vars(formula(models[[out_flag]]$fitted)[[2]])
  if (!all(cen_flag == 0)) {
    censor <- all.vars(formula(models[[cen_flag]]$fitted)[[2]])
  }

  # Get minimum and maximum time
  max_time <- max(time_seq, na.rm = TRUE)
  min_time <- min(time_seq, na.rm = TRUE)

  # Run g-formula
  for (indx in seq_along(time_seq)) {
    t_index <- sort(time_seq)[indx]
    # Spin for running
    if (!is.null(progress_bar) & t_index %% 10 == 0) {
      progress_bar(amount = 0)
    }

    set(data, j = time_var, value = t_index)

    # Recode baseline variables at initiation
    if (t_index == min_time) {
      if (!is.null(init_recode)) {
        data <- within(data, eval(parse(text = init_recode)))
      }
    }

    # Recode data before simulating
    if (!is.null(in_recode)) {
      data <- within(data, eval(parse(text = in_recode)))
    }

    # Use the model to calculate the simulated value
    data <- simulate_data(data = data,
                          exposure = exposure,
                          models = models,
                          intervention = intervention[indx],
                          mediation_type = mediation_type)

    # For survival outcome
    if (is_survival) {
      # If the intervention is defined, no censoring applied.
      if (!is.null(intervention) & cen_flag != 0) {
        set(data, j = censor, value = 0)
      }

      # If censored outcome equals to 1, then the survival outcome is not observed.
      if (cen_flag != 0) {
        data[[outcome]] <- ifelse(data[[censor]] == 1, 0, data[[outcome]])
      }

      # Output data
      if (t_index == min_time) {
        out_y <- data[, c("new_ID", time_var, outcome, "Pred_Y"), with = FALSE]
      } else {
        out_y <- rbind(out_y, data[, c("new_ID", time_var, outcome, "Pred_Y"), with = FALSE])
      }

      # Continue for the not censored or no event observed
      if (cen_flag != 0) {
        data <- data[!(data[[outcome]] == 1 | data[[censor]] == 1), ]
      } else {
        data <- data[data[[outcome]] != 1, ]
      }

      # If no data left
      if (nrow(data) == 0) {
        break
      }

      # If all loop complete
      if (t_index == max_time) {
        out_y <- rbind(out_y, data[, c("new_ID", time_var, outcome, "Pred_Y"), with = FALSE])
      }
    } else {
      # Output data
      if (t_index == min_time) {
        out_y <- data[, c("new_ID", time_var, outcome, "Pred_Y"), with = FALSE]
      } else {
        out_y <- rbind(out_y, data[, c("new_ID", time_var, outcome, "Pred_Y"), with = FALSE])
      }
    }
    
    # Recode data after simulating
    if (!is.null(out_recode) & t_index != min_time) {
      data <- within(data, eval(parse(text = out_recode)))
    }
    
  }
  # loop ends here

  if(return_data){
    return(data)
  }else {
    # Calculate mean, this is faster than mean function
    sum(data[["Pred_Y"]])/length(data[["Pred_Y"]])
  }

}
