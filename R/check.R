

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
  vars_models <- unlist(sapply(models, function(x) all.vars(formula(x$call))))
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

  if (is_survival & !all.equal(sort(unique(na.omit(data[[outcome]]))), c(0, 1))) {
    # TODO: competing risk
    stop("Only binomial outcomes for survival are currently supported, and must be set to 0 and 1.", domain = "causalMed")
  }

  # Outcome model or survival model must be defined and one only
  out_flag <- sapply(models, function(mods) mods$mod_type == "outcome")
  if (is_survival & any(out_flag)) {
    stop("You cannot define the outcome model with survival/censor model at the same time.", domain = "causalMed")
  }

  if (!any(out_flag) & !is_survival) {
    stop("Outcome or survival model must be defined in the `models`.", domain = "causalMed")
  }

  # Models checking
  if (!is.list(models)) {
    stop("Models must be provided as a list", domain = "causalMed")
  }

  # Check the models class
  cls <- sapply(models, function(x) class(x) != "gmodel")
  if (any(cls)) {
    stop("Models in the list must be `gmodel` object, please use spec_model to create!")
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

  # The dynamic intervention must be the same
  is_dynamic <- sapply(intervention, function(x) any(grepl(">|<|=|!|%in%", x)))
  is_dynamic <- sapply(intervention[is_dynamic], function(x) length(unique(x)) > 1)
  if (any(is_dynamic)) {
    stop("Dynamic intervention must have the same value in the element", domain = "causalMed")
  }

  # Checking for interventions is list
  if (!is.list(intervention)) {
    stop("Intervention must be a list object", domain = "causalMed")
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

  # Check the length, each element in the list must be of length 1 or equals to the time length.
  intervention_len <- sapply(intervention, length)
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
