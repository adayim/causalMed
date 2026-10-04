
#' Specify a model for one time-varying variable
#'
#' @description
#' Describes the model for one time-varying variable: its formula, its role in
#' the data-generating process (\code{mod_type}) and how its values are
#' simulated (\code{var_type}). Nothing is fitted here; \code{\link{gformula}}
#' and \code{\link{mediation}} fit the models, in list order, and use them to
#' simulate counterfactual trajectories.
#'
#' @param formula Model formula, e.g. \code{L ~ A + lag1_L + time}. Every
#'   variable must exist in the analysis data, or be created by a recode hook.
#' @param subset Optional unquoted logical expression, e.g.
#'   \code{platnormm1 == 0}. The model is fitted on the rows where it holds and,
#'   at each simulated time step, draws only those rows; other rows keep their
#'   current value. This is how an absorbing state is modelled.
#' @param recode Optional \code{\link{recodes}} applied before this model is
#'   fitted and before its response is simulated.
#' @param var_type How values are simulated:
#'   \describe{
#'     \item{\code{"normal"}}{(default) Gaussian draws around the fitted
#'       linear model, with the residual standard deviation.}
#'     \item{\code{"binary"}}{Bernoulli draws from a logistic model.}
#'     \item{\code{"categorical"}}{Draws from a multinomial logistic model
#'       (\code{\link[nnet]{multinom}}), returned in the variable's type in the
#'       data (a factor keeps its levels, numeric codes stay numeric); the
#'       \pkg{Hmisc} package must be installed.}
#'     \item{\code{"custom"}}{A user-supplied fitting function
#'       (\code{custom_fit}) and/or simulation function (\code{custom_sim}).}
#'   }
#'   \code{"censor"} and \code{"survival"} models must be \code{"binary"};
#'   an \code{"outcome"} model must be \code{"binary"} or \code{"normal"}.
#' @param mod_type Role of the variable: \code{"covariate"} (default),
#'   \code{"exposure"}, \code{"mediator"}, \code{"outcome"} (end-of-follow-up
#'   outcome), \code{"survival"} (discrete-time event indicator) or
#'   \code{"censor"} (loss to follow-up). Under every intervention the
#'   censoring indicator is set to zero, so risks are the risks under
#'   eliminated loss to follow-up, and with \code{estimator = "gcomp"} a
#'   censoring model does not change them. It is simulated in the natural
#'   course of \code{\link{gformula}} and used by \code{estimator = "tmle"} in
#'   \code{\link{mediation}}.
#' @param custom_fit Fitting function for \code{var_type = "custom"}
#'   (ignored, with a warning, otherwise); default \code{\link[stats]{glm}}.
#'   Its \emph{name} is recorded and evaluated later, so it must be reachable
#'   from the global environment or a package namespace: a function defined
#'   inside another function or \code{local()} will not be found. A
#'   namespace-qualified name (e.g. \code{truncreg::truncreg}) also works on
#'   parallel bootstrap workers. Without \code{custom_sim}, the fitted object
#'   must have a \code{terms} component and a \code{\link[stats]{coef}}
#'   method.
#' @param custom_sim Simulation function \code{function(fit, newdata)}
#'   returning one simulated value per row of \code{newdata}. It replaces the
#'   draw implied by \code{var_type}, so it should return a draw from the
#'   model's distribution rather than a fitted mean. If it is omitted with
#'   \code{var_type = "custom"}, values are drawn from a normal distribution
#'   around the linear predictor; that is the fitted mean only under an
#'   identity link, so a fit with another link is rejected. For
#'   \code{"outcome"} and \code{"survival"} models the reported risk is always
#'   computed from the fitted model, not from \code{custom_sim}.
#' @param truncate Logical (default \code{TRUE}). Clip simulated numeric values
#'   (\code{"normal"} draws and numeric \code{custom_sim} output) to the range
#'   of the response observed in the data. This matches the default
#'   \code{sim_trunc = TRUE} of \pkg{gfoRmula}. \code{FALSE} draws from the
#'   untruncated distribution. No effect on \code{"binary"} or
#'   \code{"categorical"} variables.
#' @param ... Further arguments to the fitting function (\code{glm},
#'   \code{multinom} or \code{custom_fit}), e.g. \code{family}.
#'
#' @details
#' Choosing a model and a simulation rule that suit each variable is the
#' analyst's responsibility; the package checks only what it can verify
#' mechanically.
#'
#' @return An object of class \code{"causalMed_gmodel"}: the unevaluated
#'   fitting call plus \code{subset}, \code{recode}, \code{var_type},
#'   \code{mod_type}, \code{custom_sim} and \code{truncate}.
#'
#' @seealso \code{\link{gformula}}, \code{\link{mediation}},
#'   \code{\link{recodes}}
#'
#' @importFrom nnet multinom
#' @importFrom stats glm
#'
#' @examples
#' # A binary covariate modelled only among those not yet in the state
#' spec_model(platnorm ~ all + cmv + male + age + gvhdm1 + daysgvhd + wait,
#'            var_type = "binary", mod_type = "covariate",
#'            subset = platnormm1 == 0)
#'
#' # A count covariate: Poisson fit with a matching simulation function
#' sim_poisson <- function(fit, newdata) {
#'   rpois(nrow(newdata), predict(fit, newdata = newdata, type = "response"))
#' }
#' spec_model(daysgvhd ~ all + cmv + male + age + wait,
#'            var_type = "custom", mod_type = "covariate",
#'            custom_sim = sim_poisson, family = poisson(link = "log"))
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
