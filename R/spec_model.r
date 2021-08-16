
#' Model specification for G-formula
#'
#' @description
#'  Specify regression model for the time-varying variables. This is will be used for
#'  g-formula. An unevaluated model will be created and can further pass to \code{Gformula}.
#'
#' @param formula an object of class formula: a symbolic description of the model
#'  to be fitted. Formula for specified models passed to the model. Must be
#'  contained within the input `data.frame` when initialized. This will passed to
#' the model fitting function.
#'
#' @param subset an optional vector specifying a subset of observations to be used
#'  in the fitting process and passed to the model fitting function.
#'
#' @param recode An Optional string vector. Recoding variable before including
#'  variables in the model. Strings will be evaluated before the fitting the model.
#'
#' @param var_type A character strings specifying the "type" of response variable
#'  of the model. The possible values are: \code{"binary"}, \code{"normal"},
#'  \code{"categorical"}, and \code{"custom"}.
#'
#' @param mod_type Model type. Exposure model (\code{"exposure"}), covariate
#'  model (\code{"covariate"}), mediator model (\code{"mediator"}), outcome
#'   model (\code{"outcome"}), mediator model (\code{"survival"}) or
#'   censoring model (\code{"censor"}).
#'
#' @param custom_fit Custom model fitting function. The \code{var_type} should
#' be set to \code{"custom"} to pass the `custom_fit`.
#'
#' @param custom_sim Custom simulation function for the model, see \code{\link{sim_value}}.
#' Only two parameters will be pass to this function during the evaluation, the
#' fitted model object and the data to be used in the prediction. The first argument
#' of the function must be the fitted model object and the second must be the data
#' to be used in the simulation. If the custom simulation function was not provided,
#' the variable will be simulated using the normal distribution.
#'
#' @param ... Other parameters passed to the model fitting function or \code{custom_fit}.
#'
#' @seealso \code{\link{Gformula}}
#'
#' @details
#' This function will be used to create an unevaluated model for the g-formula.
#'
#' @return The list returned is:
#' \item{call}{An unevaluated model.}
#' \item{subset}{A string vector passed by the parameter \code{subset}.}
#' \item{recode}{A string vector passed by the parameter \code{recode}.}
#' \item{var_type}{A character string passed by the parameter \code{var_type}.}
#' \item{mod_type}{A character string passed by the parameter \code{mod_type}.}
#' \item{custom_sim}{A symbolic customized simulation function name, passed by the parameter \code{custom_sim}.}
#'
#' @importFrom nnet multinom
#' @importFrom stats glm
#'
#' @export

spec_model <- function(formula,
                       subset = NULL,
                       recode = NULL,
                       var_type = c("normal", "binary", "categorical", "custom"),
                       mod_type = c(
                         "exposure", "covariate", "mediator",
                         "outcome", "censor", "survival"
                       ),
                       custom_fit = NULL,
                       custom_sim = NULL,
                       ...) {
                         
  tmpcall <- match.call(expand.dots = TRUE)

  args_list <- tmpcall

  if(mod_type %in% c("censor", "survival") & var_type != "binary")
    stop("Only binary variable type is allowed for the survival and censor models.", domain = "causalMed")  

  if(mod_type == "outcome" & !var_type %in% c("binary", "normal"))
    stop("Only binary or normal variable type is allowed for the outcome model.", domain = "causalMed")
  
  # Remove unnecessary arguments
  args_list <- args_list[!names(args_list) %in% c(
    "recode", "mod_type", "custom_fit",
    "custom_sim", "var_type"
  )]

  var_type <- match.arg(var_type)

  if (var_type == "custom" & is.null(custom_fit)) {
    stop("`custom_fit` must be defined if `var_type` == 'custom'.",
      domain = "causalMed"
    )
  }

  if (var_type == "categorical") {
    args_list[[1]] <- substitute(nnet::multinom)
  } else if (var_type == "normal") {
    args_list[[1]] <- substitute(stats::glm)
    args_list$family <- substitute(gaussian)
  } else if (var_type == "binary") {
    args_list[[1]] <- substitute(stats::glm)
    args_list$family <- substitute(binomial())
  } else {
    args_list[[1]] <- substitute(custom_fit)
  }

  out <- list(
    call = args_list,
    subset = tmpcall$subset,
    recode = tmpcall$recode,
    var_type = tmpcall$var_type,
    mod_type = tmpcall$mod_type,
    custom_sim = tmpcall$custom_sim
  )

  class(out) <- "gmodel"

  return(out)
}
