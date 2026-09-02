
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
#' @param subset Optional. An unquoted logical expression naming columns of the
#'   analysis data (e.g. \code{platnormm1 == 0}) that restricts this model to a
#'   subset of observations. It is passed through when fitting and re-evaluated
#'   at each time step during simulation, so only the rows satisfying it have
#'   their response drawn from this model. Rows that never satisfy it keep
#'   whatever value they already carry.
#'
#' @param recode Optional. One or more recoding statements built with
#'   \code{\link{recodes}} (e.g. \code{recodes(L_lag1 = L)} or
#'   \code{recodes(M_lag1 = 0)}), applied \strong{before} fitting the model and
#'   \strong{before} simulating its response (useful for dynamic recoding).
#'   Anything that is not a \code{recodes()} object is rejected.
#'
#' @param var_type Character. The response type for simulation/prediction:
#'   \code{"normal"} (the default), \code{"binary"}, \code{"categorical"}, or
#'   \code{"custom"}.
#'   By default, values are simulated via:
#'   \itemize{
#'     \item \code{"binary"}: Bernoulli draws using the fitted mean.
#'     \item \code{"normal"}: Gaussian draws using fitted mean, \strong{clipped to
#'           the observed range} of the response (see \code{truncate}).
#'     \item \code{"categorical"}: Multinomial draws via \code{\link[nnet]{multinom}}.
#'     \item \code{"custom"}: user-specified via \code{custom_fit} and/or \code{custom_sim}
#'           (numeric output is also clipped to the observed range unless
#'           \code{truncate = FALSE}).
#'   }
#'
#' @param mod_type Character. The role of this model in the data-generating
#'   process: \code{"covariate"} (the default), \code{"exposure"},
#'   \code{"mediator"}, \code{"outcome"}, \code{"censor"}, or
#'   \code{"survival"}.
#'
#'   \code{"censor"} declares a discrete-time censoring process. Its role is
#'   narrower than it looks: under every intervention the censoring indicator
#'   is set to zero, so the target is the risk under eliminated loss to
#'   follow-up, and with \code{estimator = "gcomp"} the censoring model does
#'   not alter the intervention-specific risks. It is used in the
#'   natural-course simulation of \code{\link{gformula}} and in the targeted
#'   estimator (\code{estimator = "tmle"}). Only \code{var_type = "binary"}
#'   is accepted for \code{"censor"} and \code{"survival"}.
#'
#' @param custom_fit Optional. A model fitting function, used \strong{only}
#'   when \code{var_type = "custom"} (it is ignored, with a warning, for the
#'   other types). If \code{var_type = "custom"} and \code{custom_fit} is not
#'   provided, \code{\link[stats]{glm}} is used by default.
#'   This can be used to define a fitting function other than
#'   \code{\link[stats]{glm}} and \code{\link[nnet]{multinom}}.
#'
#'   \strong{Scope:} \code{spec_model()} records the function \emph{name}, not
#'   the function object, and the fit is evaluated inside the package. The name
#'   must therefore resolve from the global environment or from a package
#'   namespace — a fitter defined inside another function or inside
#'   \code{local()} will not be found. Prefer a fully qualified name (e.g.
#'   \code{truncreg::truncreg} for truncated regression), which is also what
#'   makes the bootstrap work under a parallel \code{\link[future]{plan}}, where
#'   each worker is a fresh session.
#'
#'   \strong{What the fitted object must provide:} without \code{custom_sim}
#'   the simulation evaluates the fit's linear predictor, so it needs a
#'   \code{terms} component and a \code{\link[stats]{coef}} method; a fit
#'   lacking either is rejected with an explanatory error naming the variable.
#'   With \code{custom_sim} the drawing is delegated, and neither is required.
#'   Choosing a model that is appropriate for the variable, and a
#'   \code{custom_sim} that draws from the distribution that model implies,
#'   remains the analyst's responsibility.
#'
#' @param custom_sim Optional. A simulation function for the model. It must
#'   accept two arguments -- the fitted model object and a \code{data.frame} of
#'   new data to predict on -- and return a vector of simulated responses of
#'   matching length. When supplied it takes priority over the drawing rule
#'   implied by \code{var_type}. If omitted and \code{var_type = "custom"},
#'   normal draws are used by default, which requires the fitted object to have
#'   a \code{terms} component and a \code{\link[stats]{coef}} method.
#'
#'   \strong{It supplies a draw, not a fitted value.} At each time step the
#'   Monte Carlo engine assigns this variable the vector \code{custom_sim}
#'   returns, in place of the draw the built-in \code{var_type} rule would have
#'   made (a Bernoulli draw for \code{"binary"}, a Gaussian draw around the
#'   linear predictor for \code{"normal"}). Returning \code{predict()}'s fitted
#'   value therefore assigns the conditional mean rather than a draw from the
#'   conditional distribution. Whether that is appropriate for the variable
#'   being simulated is the analyst's decision.
#'
#'   \strong{It does not apply to the outcome.} For \code{mod_type = "outcome"}
#'   or \code{"survival"} the reported risk is computed from the fitted model's
#'   own linear predictor, so \code{custom_sim} affects only the simulated
#'   response value, never the estimate. Those models must therefore have
#'   extractable coefficients.
#'
#' @param truncate Logical. If \code{TRUE} (default), simulated numeric values are
#'   clipped to the range of the response observed in the data — i.e. a
#'   \code{"normal"} draw is clipped to \code{[min, max]} of the observed
#'   response, and numeric output of \code{custom_sim} is clipped as well.
#'   \code{TRUE} is the default. It corresponds to the \code{sim_trunc}
#'   argument of \pkg{gfoRmula}, documented there as "whether to truncate
#'   simulated covariates to their range in the observed data set", whose
#'   default is also \code{TRUE}. Set \code{FALSE} to draw from the untruncated
#'   fitted distribution (corresponding to \code{sim_trunc = FALSE}), or to let
#'   a \code{custom_sim} function be authoritative over its own output range.
#'   Which is appropriate depends on the variable being simulated and is the
#'   analyst's decision.
#'   Has no effect on \code{"binary"} or \code{"categorical"} responses.
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
#' An object of class \code{"causalMed_gmodel"}:
#' \describe{
#'   \item{\code{call}}{An unevaluated call (function + arguments) to fit the model.}
#'   \item{\code{subset}}{The unevaluated subset expression provided via \code{subset}.}
#'   \item{\code{recode}}{The recoding statements provided via \code{recode}.}
#'   \item{\code{var_type}}{The response type, as provided.}
#'   \item{\code{mod_type}}{The model role, as provided.}
#'   \item{\code{custom_sim}}{The simulation function provided via \code{custom_sim}.}
#'   \item{\code{truncate}}{The truncation flag, as provided.}
#' }
#'
#' @importFrom nnet multinom
#' @importFrom stats glm
#'
#' @examples
#' data(gvhd)
#' mod_cov1 <- spec_model(platnorm ~ all + cmv + male + age + agecurs1 +
#'   agecurs2 + gvhdm1 + daysgvhd + daysnorelapse + wait,
#' var_type = "binary",
#' mod_type = "covariate",
#' subset = platnormm1 == 0
#' )
#' ## For Poisson regression
#' predict_poisson <- function(fit, newdf) {
#'   theta <- stats::predict(object = fit, type = "response", newdata = newdf)
#'   prediction <- rpois(n = nrow(newdf), lambda = theta)
#'   return(prediction)
#' }
#' mod_cov1 <- spec_model(platnorm ~ all + cmv + male + age + agecurs1 +
#'   agecurs2 + gvhdm1 + daysgvhd + daysnorelapse + wait,
#' var_type = "custom",
#' mod_type = "covariate",
#' subset = platnormm1 == 0,
#' custom_sim = predict_poisson,
#' family = "poisson"(link = "log"),
#' y = TRUE
#' )
#' @export

spec_model <- function(formula,
                       subset = NULL,
                       recode = NULL,
                       var_type = c("normal", "binary", "categorical", "custom"),
                       mod_type = c(
                         "covariate", "exposure", "mediator",
                         "outcome", "censor", "survival"
                       ),
                       custom_fit = NULL,
                       custom_sim = NULL,
                       truncate = TRUE,
                       ...) {
  tmpcall <- match.call(expand.dots = TRUE)

  var_type <- match.arg(var_type)
  mod_type <- match.arg(mod_type)

  if (!is.logical(truncate) || length(truncate) != 1L || is.na(truncate)) {
    stop("`truncate` must be a single TRUE or FALSE.", domain = "causalMed")
  }

  check_recode_param("recode", recode)

  if (!inherits(formula, "formula")) {
    stop("`formula` is not a formula object.", domain = "causalMed")
  }

  # custom_sim is called as custom_sim(fitted_model, newdata) by sim_value();
  # anything else fails much later, inside the Monte Carlo loop.
  if (!is.null(custom_sim)) {
    if (!is.function(custom_sim)) {
      stop("`custom_sim` must be a function taking the fitted model and a ",
           "data.frame of new data.", domain = "causalMed")
    }
    if (!is.primitive(custom_sim) && length(formals(custom_sim)) < 2L) {
      stop("`custom_sim` must accept two arguments: the fitted model object ",
           "and a data.frame of new data.", domain = "causalMed")
    }
  }

  # custom_fit only replaces the fitting call when var_type = "custom";
  # supplying it otherwise silently had no effect.
  if (!is.null(custom_fit) && var_type != "custom") {
    warning(sprintf(paste0(
      "`custom_fit` is only used when var_type = \"custom\"; it is ignored for ",
      "var_type = \"%s\". The model will be fitted with the default function ",
      "for that type."), var_type), call. = FALSE, domain = "causalMed")
  }

  args_list <- tmpcall

  if (mod_type %in% c("censor", "survival") & var_type != "binary") {
    stop("Only binary variable type is allowed for the survival and censor models.", domain = "causalMed")
  }

  if (mod_type == "outcome" & !var_type %in% c("binary", "normal")) {
    stop("Only binary or normal variable type is allowed for the outcome model.", domain = "causalMed")
  }

  # Remove unnecessary arguments
  args_list <- args_list[!names(args_list) %in% c(
    "recode", "mod_type", "custom_fit",
    "custom_sim", "var_type", "truncate"
  )]

  if (var_type == "categorical") {
    args_list[[1]] <- substitute(nnet::multinom)
  } else if (var_type == "normal") {
    args_list[[1]] <- substitute(stats::glm)
    args_list$family <- substitute(gaussian)
  } else if (var_type == "binary") {
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
    custom_sim = custom_sim,
    truncate = truncate
  )

  class(out) <- c("causalMed_gmodel", "list")

  return(out)
}
