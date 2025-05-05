
#' Simulate Data
#'
#' Loop through the models and apply any recoding or subset.
#'
#' @param data Data to be used for the data generation
#' @param models Model list passed from \code{\link{gformula}} or \code{\link{mediation}}.
#' @param intervention A vector, intervention treatment, natural course will be evaluated
#' if the value is \code{NULL}. If the value contains any logical operators, the intervention
#' will be evaluated and the exposure variable will be set to 1 if it is \code{TRUE}.
#' @param mediation_type Type of the mediation analysis, if the value is \code{NA}
#' no mediation analysis will be performed (default). It will be ignored if the intervention
#'  is not \code{NULL}
#'
#' @details
#' If the intervention is \code{NULL}, and the \code{mediation_type}
#'
#'
#' @keywords internal
#'
simulate_data <- function(data,
                          exposure,
                          models,
                          intervention = NULL,
                          mediation_type = c(NA, "N", "I")) {
  mediation_type <- match.arg(mediation_type)

  # Check if the intervention is dynamic
  if (!is.null(intervention)) {
    is_dynamic <- any(grepl(">|<|=|!|%in%", intervention))
  } else {
    is_dynamic <- FALSE
  }

  # If Intervention is given, set the treatment to given value
  if (!is.null(intervention) & !is_dynamic) {
    set(data, j = exposure, value = intervention)
  }

  

  # if the mediation type is defined, than set the intervention to 1. But 0
  # for mediator. This is to calculate the phi_10
  if (!is.na(mediation_type) & is.null(intervention)) {
    set(data, j = exposure, value = 1)
    interv0 <- parse(text = paste0(exposure, " = 0"))
  }

  # Loop through models
  for (indx in seq_along(models)) {
    model <- models[[indx]]
    # Get the name of the response variable for the current model
    resp_var <- model$rsp_vars
    mod_type <- model$mod_type

    # Skip if the response variable is treatment and the intervention is defined.
    if (resp_var == exposure & !is.null(intervention)) {
      # Evaluate if the intervention is dynamic
      if (is_dynamic) {
        set(data, j = exposure, value = as.numeric(eval(parse(text = intervention), envir = data)))
        # data <- within(data, eval(parse(text = intervention)))
      }

      next
    }

    # Perform the recode in the model
    if (!is.null(model$recode)) {
      data <- within(data, eval(parse(text = model$recode)))
    }

    # If condition has been defined, apply it
    if (!is.null(model$subset)) {
      cond <- with(data, eval(model$subset))
      cond <- cond & !is.na(cond) # Avoid NA in the condition list.
    } else {
      cond <- rep(TRUE, nrow(data))
    }

    if (sum(cond) != 0) {
      if (mod_type == "mediator" & !is.na(mediation_type) & is.null(intervention)) {
        med_value <- sim_value(model = model, newdt = within(data[cond, ], eval(interv0)))
        # Set the mediator's intervention to 0
        if (mediation_type == "N") {
          data[[resp_var]][cond] <- med_value
        } else {
          data[[resp_var]][cond] <- sample(med_value, length(med_value), replace = TRUE)
        }
      } else {
        data[[resp_var]][cond] <- sim_value(model = model, newdt = data[cond, ])
      }

      # Get the predicted values for the outcome
      if (mod_type %in% c("outcome", "survival")) {
        data[["Pred_Y"]][cond] <- predict(model$fitted, newdata = data[cond, ], type = "response")
      }
    }
  }

  return(data)
}

