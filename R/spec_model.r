
#' Model specification for G-formula
#'
#' @description
#'  Specify regression model for the time-varying variables. This is will be used for
#'  g-formula. This function will create unevaluated model and can further pass to
#' \code{\link{gformula}} or \code{\link{mediation}}.
#'
#' @param formula an object of class formula: a symbolic description of the model
#'  to be fitted. Formula for specified models passed to the model. Must be
#'  contained within the input `data.frame` when initialized. This will passed to
#' the model fitting function.
#'
#' @param subset an optional vector specifying a subset of observations to be used
#'  in the fitting process and passed to the model fitting function. This wil also
#' be used during the simulation of the response variable.
#'
#' @param recode An Optional string vector. Recoding variable before evaluating
#'  the model and before simulating the outcome variable in this model. This can
#'  be used for dynamic recoding of the variables.
#'
#' @param var_type A character strings specifying the "type" of response variable
#'  of the model. The possible values are: \code{"binomial"}, \code{"normal"},
#'  \code{"categorical"}, and \code{"custom"}. The binomial distribution, normal distribution
#' or multinomial distribution will be used to simulate the value of the responsible variable
#'  for \code{"binomial"}, \code{"normal"} and \code{"categorical"}, respectively.
#' The simulation of the response variable will be based on normal distribution if the \code{custom_sim}
#' is \code{NULL}. The simulated value of the \code{"normal"} and \code{"custom"} variable type
#' will be truncated within the observed value.
#'
#' @param mod_type Type of the model. Covariate model (\code{"covariate"}), exposure model
#' (\code{"exposure"}), mediator model (\code{"mediator"}), outcome model (\code{"outcome"}),
#'  survival model (\code{"survival"}) or censoring model (\code{"censor"}).
#'
#' @param custom_fit Custom model fitting function. The \code{var_type} should be set to
#' \code{"custom"} to pass the `custom_fit`. If the \code{var_type} is \code{"custom"} and the
#' \code{custom_fit} is node defined, the \code{\link[stats]{glm}} will be used for modelling.
#' This can be used to define the modeling fitting function other than \code{\link[stats]{glm}}
#'  and \code{\link[nnet]{multinom}}. The  function should include the package name to work in
#'  parallel computation. For example, if a truncated regression to be used for this model,
#'  one should define the parameter value with \code{truncreg::truncreg}.
#'
#' @param custom_sim Custom simulation function for the model. Only two parameters will be passed
#' to this function during the simulation. The first argument is the fitted model object
#'  and the second one is the data to be used in the prediction. The simulated value will be
#' truncated within observed value range. If the custom simulation function was not provided
#'  and the variable type is custom, the variable will be simulated using the normal distribution.
#'
#' @param ... Other parameters passed to the model fitting function, \code{\link[stats]{glm}},
#'  \code{\link[nnet]{multinom}} or \code{custom_fit}.
#'
#' @seealso \code{\link{gformula}},\code{\link{mediation}}
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
#'
#' @examples
#' mod_cov1 <- spec_model(platnorm ~ all + cmv + male + age + agecurs1 +
#'   agecurs2 + gvhdm1 + daysgvhd + daysnorelapse + wait,
#' var_type = "binomial",
#' mod_type = "covriate",
#' subset = platnormm1 == 0
#' )
#' ## For Poisson regression
#' predict_poisson <- function(fit, newdf) {
#'   theta <- stats::predict(object = fitcov, type = "response", newdata = newdf)
#'   prediction <- rpois(n = nrow(newdf), lambda = theta)
#'   return(prediction)
#' }
#' mod_cov1 <- spec_model(platnorm ~ all + cmv + male + age + agecurs1 +
#'   agecurs2 + gvhdm1 + daysgvhd + daysnorelapse + wait,
#' var_type = "custom",
#' mod_type = "covriate",
#' subset = platnormm1 == 0,
#' custom_sim = predict_poisson,
#' family = "poisson"(link = "log"),
#' y = TRUE
#' )
spec_model <- function(formula,
                       subset = NULL,
                       recode = NULL,
                       var_type = c("normal", "binomial", "categorical", "custom"),
                       mod_type = c(
                         "covariate", "exposure", "mediator",
                         "outcome", "censor", "survival"
                       ),
                       custom_fit = NULL,
                       custom_sim = NULL,
                       ...) {
  tmpcall <- match.call(expand.dots = TRUE)

  args_list <- tmpcall

  if (mod_type %in% c("censor", "survival") & var_type != "binomial") {
    stop("Only binomial variable type is allowed for the survival and censor models.", domain = "causalMed")
  }

  if (mod_type == "outcome" & !var_type %in% c("binomial", "normal")) {
    stop("Only binomial or normal variable type is allowed for the outcome model.", domain = "causalMed")
  }

  # Remove unnecessary arguments
  args_list <- args_list[!names(args_list) %in% c(
    "recode", "mod_type", "custom_fit",
    "custom_sim", "var_type"
  )]

  var_type <- match.arg(var_type)

  if (var_type == "categorical") {
    args_list[[1]] <- substitute(nnet::multinom)
  } else if (var_type == "normal") {
    args_list[[1]] <- substitute(stats::glm)
    args_list$family <- substitute(gaussian)
  } else if (var_type == "binomial") {
    args_list[[1]] <- substitute(stats::glm)
    args_list$family <- substitute(binomial())
  } else if (var_type == "custom" & is.null(custom_fit)) {
    args_list[[1]] <- substitute(stats::glm)
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
