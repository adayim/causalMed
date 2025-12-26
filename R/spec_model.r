
#' Model specification for G-formula
#'
#' @description
#'  Specify regression models for the time-varying variables to be used
#' within the g-formula simulation. This function creates unevaluated models and 
#' can further pass to \code{\link{gformula}} or \code{\link{mediation}}.
#'
#' @param formula an object of class formula: An object of class \code{formula}: symbolic model specification
#'   to be fitted (e.g., \code{Y ~ A + L + time}). The variables referenced in
#'   \code{formula} must exist in the analysis \code{data.frame} when fitting/evaluating
#'   the model during g-formula simulation.
#'
#' @param subset Optional. A logical expression or character vector defining a subset
#'   of observations to fit/simulate this model on; passed through when fitting and also
#'   respected during simulation of the response for this model.
#'
#' @param recode Optional character vector. Should be defined with \code{\link{recodes}}. One or more recoding statements (e.g.,
#'   \code{L_lag1 = L} or \code{M_lag1 = 0}) to be applied **before** evaluating
#'   the model and **before** simulating the response for this model (useful for dynamic recoding).
#'
#' @param var_type Character. The response type for simulation/prediction:
#'   \code{"binomial"}, \code{"normal"}, \code{"categorical"}, or \code{"custom"}.
#'   By default, values are simulated via:
#'   \itemize{
#'     \item \code{"binomial"}: Bernoulli draws using the fitted mean.
#'     \item \code{"normal"}: Gaussian draws using fitted mean.
#'     \item \code{"categorical"}: Multinomial draws via \code{\link[nnet]{multinom}}.
#'     \item \code{"custom"}: user-specified via \code{custom_fit} and/or \code{custom_sim}.
#'   }
#'
#' @param mod_type Character. The role of this model in the data-generating process:
#'   \code{"covariate"}, \code{"exposure"}, \code{"mediator"},
#'   \code{"outcome"}, \code{"censor"}, or \code{"survival"}.
#'
#' @param custom_fit Optional. A model fitting function used
#'   when \code{var_type = "custom"}. If \code{var_type = "custom"} and
#'   \code{custom_fit} is not provided, \code{\link[stats]{glm}} is used by default.
#' This can be used to define the modeling fitting function other than \code{\link[stats]{glm}}
#'  and \code{\link[nnet]{multinom}}. For parallel computation, provide the fully
#'   qualified function name (e.g., \code{truncreg::truncreg} for truncated regression).
#'
#' @param custom_sim Optional. A simulation function for the model. It should
#'   accept exactly two arguments: the fitted model object and a \code{data.frame} of
#'   new data to predict on, and return a vector of simulated responses of matching length.
#'   If omitted and \code{var_type = "custom"}, normal draws are used by default.
#'
#' @param ... Other parameters passed to the model fitting function, \code{\link[stats]{glm}},
#'  \code{\link[nnet]{multinom}} or \code{custom_fit}.
#'
#' @seealso \code{\link{gformula}},\code{\link{mediation}}
#'
#' @details
#' This function will be used to create an unevaluated model for the g-formula.
#' \code{spec_model()} does not fit the model immediately. It returns an unevaluated
#' call plus metadata (\code{var_type}, \code{mod_type}, \code{subset}, \code{recode},
#' \code{custom_sim}) that are used later by \code{\link{gformula}}/\code{\link{mediation}}
#' to fit in temporal order and to simulate counterfactual trajectories.
#' 
#' @return 
#' An object of class \code{"gmodel"}:
#' \describe{
#'   \item{\code{call}}{An unevaluated call (function + arguments) to fit the model.}
#'   \item{\code{subset}}{The subset expression/character vector provided via \code{subset}.}
#'   \item{\code{recode}}{The recoding statements provided via \code{recode}.}
#'   \item{\code{var_type}}{The response type, as provided.}
#'   \item{\code{mod_type}}{The model role, as provided.}
#'   \item{\code{custom_sim}}{A symbolic customized simulation function name, as provided.}
#' }
#'
#' @importFrom nnet multinom
#' @importFrom stats glm
#'
#' @examples
#' data(gvhd)
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
#' @export

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

  var_type <- match.arg(var_type)
  mod_type <- match.arg(mod_type)

  check_recode_param("recode", recode)

  is.formula <- function(x) is.call(x) && x[[1]] == quote(`~`)
  if(!is.formula(formula)) {
    stop("`formula` is not a formula object.", domain = "causalMed")
  }

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
    subset = substitute(subset),
    recode = recode,
    var_type = var_type,
    mod_type = mod_type,
    custom_sim = custom_sim
  )

  class(out) <- "gmodel"

  return(out)
}
