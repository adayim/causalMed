

#' Error Catching
#'
#' This function is used to check for errors in the \code{\link{gformula}}.
#'
#' @inheritParams gformula
#' @return No value is returned.
#' @keywords internal
#'
check_error <- function(data,
                        id_var,
                        base_vars,
                        exposure,
                        time_var,
                        models) {
  if (!(exists("data") && is.data.frame(get("data")))) {
    stop("Data does not exist or not a data.frame.", domain = "causalMed")
  }


  if (!is.list(models)) {
    stop("Models must be provided as a list", domain = "causalMed")
  }

  out_flag <- sapply(models, function(mods) mods$mod_type %in% c("outcome", "survival"))
  if (sum(out_flag) > 1) {
    stop("Only one outcome model or survival model allowed.", domain = "causalMed")
  }

  if (sum(out_flag) == 0) {
    stop("Outcome model or survival model must be defined.", domain = "causalMed")
  }

  out_flag <- which(out_flag)
  outcome <- all.vars(formula(models[[out_flag]]$call)[[2]])

  # Check variables in the data
  check_var_in(c(id_var, base_vars, exposure, outcome, time_var), data)

  if(!is.numeric(data[[time_var]])) {
    stop("The time variable must be numeric.", domain = "causalMed")
  }

  # Check if variables in the formula included in the data
  vars_models <- unlist(lapply(models, function(x) all.vars(formula(x$call))))
  check_var_in(vars_models, data)

  # Get the variable name of exposure and check if it the same as the exposure
  exp_flag <- sapply(models, function(mods) mods$mod_type == "exposure")
  if (sum(exp_flag) > 1) {
    stop("Only one exposure model is allowed.", domain = "causalMed")
  }

  if (any(exp_flag)) {
    exp_flag <- which(exp_flag)
    exposure_var <- all.vars(formula(models[[exp_flag]]$call)[[2]])
    if (exposure_var != exposure) {
      stop("The given exposure variable was different between exposure model in `models`.",
        domain = "causalMed"
      )
    }
  }

  # Check for survival models
  cen_flag <- sapply(models, function(mods) mods$mod_type == "censor")
  surv_flag <- sapply(models, function(mods) mods$mod_type == "survival")
  if (sum(cen_flag) > 1) {
    stop("Only one censor model is allowed.", domain = "causalMed")
  }
  if (sum(surv_flag) > 1) {
    stop("Only one survival model is allowed.", domain = "causalMed")
  }

  if (any(cen_flag) | any(surv_flag)) {
    is_survival <- TRUE
  }else{
    is_survival <- FALSE
  }

  if (is_survival) {
    y <- as.numeric(na.omit(data[[outcome]]))
    eps <- 1e-8
    if (!all(abs(y - 0) < eps | abs(y - 1) < eps)) {
      stop("For survival, outcome must be binary {0,1}.")
    }
  }

  # Outcome model or survival model must be defined and one only
  out_flag <- sapply(models, function(mods) mods$mod_type == "outcome")
  if (is_survival & any(out_flag)) {
    stop("You cannot define the outcome model with survival/censor model at the same time.", domain = "causalMed")
  }

  if (!any(out_flag) & !is_survival) {
    stop("Outcome or survival model must be defined in the `models`.", domain = "causalMed")
  }

  # Block user columns named "S" or "Sc" to avoid clashes with internal variables for survival outcome.
  if (isTRUE(is_survival)) {
    reserved <- c("S", "Sc")
    hit <- intersect(names(data), reserved)
    if (length(hit) > 0) {
      msg <- sprintf(
        "Column name(s) conflict with reserved internal names: %s. Please rename.",
        paste(hit, collapse = ", ")
      )
      stop(msg, call. = FALSE)
    }
  }



  # Check the models class
  cls <- vapply(models, function(x) !inherits(x, "causalMed_gmodel"), logical(1))
  if (any(cls)) {
    stop("Models in the list must be `causalMed_gmodel` object, please use spec_model to create!")
  }
}

#' Check for the intervention
#'
#' Check if the intervention is correctly defined.
#'
#' @inheritParams gformula
#' @param time_len length of the time in the data.
#'
#' @keywords internal
#'
check_intervention <- function(models, intervention, ref_int, time_len) {
  # Check if contain mediator
  med_flag <- sapply(models, function(mods) as.numeric(mods$mod_type == "mediator"))
  if (!all(med_flag == 0) & !is.null(intervention)) {
    stop("You cannot specify intervention with mediator at the same time, please set `intervention` to `NULL` or remove the mediator model in the `models`.",
      domain = "causalMed"
    )
  }

  # Checking for interventions is list
  if (!is.list(intervention)) {
    stop("Intervention must be a list object", domain = "causalMed")
  }

  # Each element must be NULL, numeric/logical (static), or a dyn_int() object (dynamic)
  bad <- sapply(intervention, function(x) {
    !is.null(x) && !is.numeric(x) && !is.logical(x) && !inherits(x, "causalMed_dynint")
  })
  if (any(bad)) {
    stop(
      "Each intervention element must be NULL, a numeric/logical value, or a dyn_int() object.",
      domain = "causalMed"
    )
  }

  # Check if all null
  interv_value <- sapply(intervention, is.null)
  if (sum(interv_value) > 1) {
    stop("All intervention are NULL.", domain = "causalMed")
  }

  # Check for intervention is a named list
  if (is.null(names(intervention)) | any(nchar(names(intervention)) == 0)) {
    stop("Intervention must be a named list object", domain = "causalMed")
  }

  # Check the length: each element must be NULL (0), scalar (1), or full-length.
  # dyn_int() objects are exempt — they apply a rule at every time step.
  intervention_len <- sapply(intervention, function(x) {
    if (inherits(x, "causalMed_dynint")) 1L else length(x)
  })
  if (!all(intervention_len %in% c(0, 1, time_len))) {
    stop("Length of the elements in the intervention must be 0 (`NULL`), 1 or the same length of time.", domain = "causalMed")
  }

  # Check for reference intervention
  if (length(ref_int) > 1) {
    stop("Length of `ref_int` must be 1.", domain = "causalMed")
  }

  if (is.numeric(ref_int) & ref_int > length(intervention)) {
    stop("`ref_int` must be less than the length of intervention.", domain = "causalMed")
  }

  if (is.character(ref_int) & ref_int != "natural" & !ref_int %in% names(intervention)) {
    stop("`ref_int` must be included in the names of intervention.", domain = "causalMed")
  }
}

#' Check temporal ordering of models for mediation analysis
#'
#' Issues a warning if covariate or exposure models appear in positions that
#' violate the assumed A(t) -> M(t) -> L(t) -> S(t) ordering.
#'
#' @param models List of model specifications from \code{\link{spec_model}}.
#' @keywords internal
check_mediation_order <- function(models) {
  mod_types <- sapply(models, function(m) m$mod_type)
  med_idx <- which(mod_types == "mediator")
  out_idx <- which(mod_types %in% c("outcome", "survival"))
  exp_idx <- which(mod_types == "exposure")

  if (length(med_idx) == 0) return(invisible(NULL))

  # Mediator must come before outcome/survival
  if (length(out_idx) > 0 && med_idx > min(out_idx)) {
    warning(sprintf(
      paste0(
        "Temporal ordering: mediator model at position [%d] appears after ",
        "the outcome/survival model at position [%d]. ",
        "The mediator must be simulated before the outcome."
      ),
      med_idx, min(out_idx)
    ), call. = FALSE)
  }

  # Exposure model must come before mediator
  if (length(exp_idx) > 0 && exp_idx > med_idx) {
    warning(sprintf(
      paste0(
        "Temporal ordering: exposure model at position [%d] appears after ",
        "the mediator model at position [%d]. ",
        "Exposure must be set before the mediator is simulated."
      ),
      exp_idx, med_idx
    ), call. = FALSE)
  }

  # Any covariate/mediator/exposure model after the outcome is almost certainly wrong
  if (length(out_idx) > 0) {
    after_out <- which(seq_along(mod_types) > max(out_idx) &
                       mod_types %in% c("covariate", "mediator", "exposure"))
    if (length(after_out) > 0) {
      warning(sprintf(
        paste0(
          "Temporal ordering: model(s) at position(s) [%s] appear after the ",
          "outcome/survival model. The outcome model should be last in the list."
        ),
        paste(after_out, collapse = ", ")
      ), call. = FALSE)
    }
  }

  invisible(NULL)
}


#' Check variables in data
#'
#' Check if the variables in the data, throw an error with variable names if not.
#'
#' @param vars Variables to check
#' @param data Data set to be checked
#'
#' @keywords internal
#'
check_var_in <- function(vars, data) {
  diff_vars <- setdiff(vars, names(data))
  if (!identical(diff_vars, character(0))) {
    stop("The following variables cannot be found in the data: ",
      paste(diff_vars, collapse = ", "),
      domain = "causalMed"
    )
  }
}


#' Validate recode parameters
#'
#' @param param_name The name of the parameter (for error messages).
#' @param param_value The object passed by the user.
#' @return TRUE if valid, stops execution otherwise.
check_recode_param <- function(param_name, param_value) {

  if (is.null(param_value)) return(TRUE)

  # Strict Check: Must be your custom class
  if (!inherits(param_value, "causalMed_recodes")) {
    stop(sprintf(
      "Invalid input for '%s'. You must use the recodes() helper function.\n  Correct: %s = recodes(x = y^2)",
      param_name, param_name
    ), call. = FALSE, domain = "causalMed")
  }

  return(TRUE)
}

