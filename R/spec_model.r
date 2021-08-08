
#' Model specification for G-formula
#'
#' @description
#'  Add a specified regression model for the exposure. This is used for natural
#'  course estimation of the Monte Carlo g-formula. This must be specified before
#'  calling the fit function.
#'
#' @param formula an object of class formula: a symbolic description of the model
#'  to be fitted. Formula for specified models passed to the model. Must be
#'  contained within the input `data.frame` when initialized.
#'
#' @param subset an optional vector specifying a subset of observations to be used
#'  in the fitting process (see \code{\link[stats]{glm}}).
#'
#'  @param recode An Optional string vector. Recoding variable before including
#'  variables in the model. Strings will be evaluated before the fitting the model.
#'
#' @param var_type A character strings specifying the "type" of response variable
#'  of the model. The possible values are: \code{"binary"}, \code{"normal"},
#'  \code{"categorical"}, and \code{"custom"}.
#'
#' @param mod_type Model type. Exposure model (\code{exposure}), covariate
#'  model (\code{covariate}), mediator model (\code{mediator}), outcome
#'   model (\code{outcome}), mediator model (\code{survival}) or
#'   censoring model (\code{censor}).
#'
#' @param custom_fn Custom function to create the model. The \code{var_type} should
#' be set to \code{"custom"} to pass the `custom_fn`.
#'
#' @param ... Other parameters passed to the model.
#'
#' @details
#'
#' @import nnet multinom
#'
#' @export

spec_model <- function(formula,
                       subset   = NULL,
                       recode   = NULL,
                       var_type = c("normal", "binary", "categorical", "custom"),
                       mod_type = c("exposure", "covariate", "mediator",
                                    "outcome", "censor", "survival"),
                       custom_fn = NULL,
                       ...){

  tmpcall <- match.call(expand.dots = TRUE)

  args_list <- tmpcall

  # Remove unnecessary arguments
  args_list <- args_list[!names(args_list) %in% c("recode", "mod_type",
                                                 "custom_fn", "var_type")]

  var_type <- match.arg(var_type)

  if(var_type == "custom" & is.null(custom_fn))
    stop("`custom_fn` must be defined if `var_type` == 'custom'.",
         domain = "causalMed")

  if(var_type == "categorical"){
    args_list[[1]] <- substitute(nnet::multinom)

  }else if(var_type == "normal"){
    args_list[[1]] <- substitute(stats::glm)
    args_list$family <- substitute(gaussian)

  }else if(var_type == "binary"){
    args_list[[1]] <- substitute(stats::glm)
    args_list$family <- substitute(binomial())

  }else{
    args_list[[1]] <- substitute(custom_fn)
  }

  out <- list(call     = args_list,
              subset   = tmpcall$subset,
              recode   = tmpcall$recode,
              var_type = tmpcall$var_type,
              mod_type = tmpcall$mod_type)

  class(out) <- "gmodel"

  return(out)
}

