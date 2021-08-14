
#' Monte Carlo simulation
#'
#' @description
#'  Internal use only. Monte Carlo simulation.
#'
#' @param data a data frame in which to look for variables with which to predict.
#'
#' @param exposure Intervention/Exposure variable
#'
#' @param time.seq Time sequence.
#'
#' @param time.var Time variable.
#'
#' @param models fitted objects.
#'
#' @param intervention A vector, intervention treatment per time.
#'
#' @param init.recode optional, recoding of variables done at the
#' beginning of the Monte Carlo loop. Needed for operations initialize baseline variables.
#' This is executed at beginning of the Monte Carlo g-formula, executed only once at time 0.
#'
#' @param in.recode optional, On the fly recoding of variables done before the Monte
#'  Carlo loop starts. Needed to do any kind of functional forms for entry times.
#'   This is executed at each start of the Monte Carlo g-formula time steps
#'
#' @param out.recode optional, On the fly recoding of variables done at the
#' end of the Monte Carlo loop. Needed for operations like counting the number of
#' days with a treatment or creating lagged variables. This is executed at each end
#'  of the Monte Carlo g-formula time steps.
#'
#' @param verbose Print intervention information during calculation.
#'
#' @importFrom pbapply timerProgressBar setTimerProgressBar
#'
#' @keywords internal
#'
monte_g <- function(data,
                    exposure,
                    time.seq,
                    time.var,
                    models,
                    intervention = NULL,
                    in.recode = NULL,
                    out.recode = NULL,
                    init.recode = NULL,
                    mediation_type = c(NA, "N", "I"),
                    verbose = TRUE) {

  mediation_type <- match.arg(mediation_type)

  # Replicate to intervention to the same length of the time.
  time_len <- length(time.seq)
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

  if (verbose) {
    cat("\n======Running intervention, please be patient.=======\n")
  }

  # Simulate
  max_time <- max(time.seq, na.rm = TRUE)
  min_time <- min(time.seq, na.rm = TRUE)
  dat_y <- data

  # for progress bar
  if (verbose){
    pb = pbapply::timerProgressBar(min = 0, max = length(time.seq), style = 2)
    on.exit(close(pb))
  }

  # Run normal g-formula
  for (t_index in sort(time.seq)) {

    if (verbose)
      pbapply::setTimerProgressBar(pb, which(sort(time.seq) == t_index))

    dat_y[[time.var]] <- t_index

    # Recode baseline variables at initiation
    if (t_index == min_time) {
      if (!is.null(init.recode)) {
        dat_y <- within(dat_y, eval(parse(text = init.recode)))
      }
    }

    # Recode data before simulating
    if (!is.null(in.recode)) {
      dat_y <- within(dat_y, eval(parse(text = in.recode)))
    }

    # Use the model to calculate the simulated value
    dat_y <- simulate_data(data = dat_y, exposure = exposure, models = models, intervention = intervention[t_index], mediation_type = mediation_type)

    if (!is.null(out.recode)) {
      dat_y <- within(dat_y, eval(parse(text = out.recode)))
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
      if (t_index == min_time){
        out_y <- dat_y[, c("new_ID", time.var, outcome), with = FALSE]
      }else{
        out_y <- rbind(out_y, dat_y[, c("new_ID", time.var, outcome), with = FALSE])
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
        out_y <- rbind(out_y, dat_y[, c("new_ID", time.var, outcome), with = FALSE])
      }

    } else {
      # Output data
      if (t_index == min_time){
        out_y <- dat_y[, c("new_ID", time.var, outcome), with = FALSE]
      }else{
        out_y <- rbind(out_y, dat_y[, c("new_ID", time.var, outcome), with = FALSE])
      }
    }
  }
  # loop ends here

  return(out_y)
}
