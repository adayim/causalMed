#' Simulate Data for Phi10 for mediantion analysis
#'
#' Loop through the models and apply any recoding or subset.
#'
#' @param data Data to be used for the data generation
#' @param models Model list passed from \code{\link{gformula}} or \code{\link{mediation}}.
#' @param mediation_type Type of the mediation analysis.
#' @param progress_bar Progress bar object.
#' @inheritParams gformula
#'
#' @details
#' If the intervention is \code{NULL}, and the \code{mediation_type}
#'
#'
#' @keywords internal

mediantion_phi10 <- function(data,
                             exposure,
                             models,
                             mediation_type = c("N", "I"),
                             time_seq,
                             time_var,
                             init_recode = NULL,
                             in_recode = NULL,
                             out_recode = NULL,
                             return_data = FALSE,
                             progress_bar = NULL){

  mediation_type <- match.arg(mediation_type)

  # Replicate to intervention to the same length of the time.
  time_len <- length(time_seq)

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

  data0 <- copy(data)

  # Run g-formula
  for (indx in seq_along(time_seq)) {
    t_index <- sort(time_seq)[indx]
    # Spin for running
    if (!is.null(progress_bar) & t_index %% 10 == 0) {
      progress_bar(amount = 0)
    }

    set(data, j = time_var, value = t_index)
    set(data0, j = time_var, value = t_index)

    # Re-code baseline variables at initiation
    if (t_index == min_time) {
      if (!is.null(init_recode)) {
        data <- within(data, eval(parse(text = init_recode)))
        data0 <- within(data0, eval(parse(text = init_recode)))
      }
    }

    # Re-code data before simulating
    if (!is.null(in_recode) & t_index != min_time) {
      data <- within(data, eval(parse(text = in_recode)))
      data0 <- within(data0, eval(parse(text = in_recode)))
    }

    # Simulate everything under no intervention
    data0 <- simulate_data(data = data0,
                           exposure = exposure,
                           models = models,
                           intervention = 0)

    # Use the data under no intervention to calculate mediator
    data <- simulate_data(data = data,
                          exposure = exposure,
                          models = models,
                          intervention = 1,
                          mediation_type = mediation_type,
                          data0 = data0)

    # For survival outcome
    if (is_survival) {
      # If the intervention is defined, no censoring applied.
      if (!is.null(intervention) & cen_flag != 0) {
        set(data, j = censor, value = 0)
        set(data0, j = censor, value = 0)
      }

      # If censored outcome equals to 1, then the survival outcome is not observed.
      if (cen_flag != 0) {
        data[[outcome]] <- ifelse(data[[censor]] == 1, 0, data[[outcome]])
        # Subset follows normal data
        data0[[outcome]] <- ifelse(data[[censor]] == 1, 0, data[[outcome]])
      }

      # Continue for the not censored or no event observed
      if (cen_flag != 0) {
        data <- data[!(data[[outcome]] == 1 | data[[censor]] == 1), ]
        data0 <- data0[!(data[[outcome]] == 1 | data[[censor]] == 1), ]
      } else {
        data <- data[data[[outcome]] != 1, ]
        data0 <- data0[data[[outcome]] != 1, ]
      }

      # If no data left
      if (nrow(data) == 0) {
        break
      }

    }

    # Recode data after simulating
    if (!is.null(out_recode) & t_index != min_time) {
      data <- within(data, eval(parse(text = out_recode)))
      data0 <- within(data0, eval(parse(text = out_recode)))
    }

  }
  # loop ends here
  return(data)
}




