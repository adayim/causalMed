

#' Error Catching
#'
#' This function is used to check for errors in the \code{\link{Gformula}}.
#'
#' @inheritParams Gformula
#' @return No value is returned.
#' @keywords internal
#'
check_error <- function(data,
                        id.var,
                        base.vars,
                        exposure,
                        outcome,
                        time.var,
                        models,
                        intervention = NULL,
                        init.recode  = NULL,
                        in.recode    = NULL,
                        out.recode   = NULL,
                        mc.sample    = 10000,
                        verbose      = TRUE){

  # Get time length
  time_len <- length(unique(data[[time.var]]))

  # Check variables in the data
  check_var_in(c(id.var, base.vars, exposure, outcome, time.var), data)

  # Check if variables in the formula included in the data
  vars_models <- unlist(sapply(models, function(x)all.vars(formula(x$call))))
  check_var_in(vars_models, data)

  # Check for survival models
  cen_flag  <- sapply(models, function(mods) as.numeric(mods$type == "censor"))
  surv_flag <- sapply(models, function(mods) as.numeric(mods$type == "survival"))
  if(!all(cen_flag == 0) | !all(surv_flag == 0)){
    is.survival <- TRUE
  }

  if(is.survival & !all.equal(sort(unique(na.omit(data[[outcome]]))), c(0, 1)))
    # TODO: competing risk
    stop('Only binary outcomes for survival are currently supported, and must be set to 0 and 1.', domain = "causalMed")

  # Check if contain mediator
  med_flag <- sapply(models, function(mods) as.numeric(mods$type == "mediator"))
  if(!all(med_flag == 0) & !is.null(intervention))
    stop("You cannot specify intervention with mediator at the same time, please set `intervention` to `NULL` or remove the mediator model in the `models`.",
         domain = "causalMed")

  # Checking for interventions is list
  if(!is.list(intervention))
    stop("Intervention must be a list object", domain = "causalMed")

  # Check for intervention is a named list
  if(is.null(names(intervention)) | any(nchar(names(intervention))))
    stop("Intervention must be a named list object", domain = "causalMed")

  # Check the length, each element in the list must be of length 1 or equals to the time length.
  intervention_len <- sapply(intervention, length)
  if(!all(length %in% c(0, 1, time_len)))
    stop("Length of the elements in the intervention must be `NULL`, 1 or the same length of time.", domain = "causalMed")

  # Models checking
  if(!is.list(models))
    stop("Models must be provided as a list", domain = "causalMed")

  # Check the models class
  cls <- sapply(models, function(x)class(x) != "gmodel")
  if(any(cls))
    stop("Models in the list must be `gmodel` object, please use spec_model to create!")

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
check_var_in <- function(vars, data){
  diff_vars <- setdiff(vars, names(data))
  if(!identical(diff_vars, character(0)))
    stop("The following variables cannot be found in the data: ",
         paste(diff_vars, collapse = ", "), domain = "causalMed")

}






