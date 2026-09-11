

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

  # Validate the model class up front: every check below reads $mod_type and
  # $call, which an element that did not come from spec_model() does not have,
  # so leaving this until the end produced a confusing downstream error first.
  cls <- vapply(models, function(x) !inherits(x, "causalMed_gmodel"), logical(1))
  if (any(cls)) {
    stop("Models in the list must be `causalMed_gmodel` object, please use spec_model to create!")
  }

  # Model roles, read once.
  mod_types <- vapply(models, function(mods) mods$mod_type, character(1))

  out_flag <- mod_types %in% c("outcome", "survival")
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

  # `base_vars` must be time-fixed. The Monte Carlo cohort is drawn from
  # unique(data[, c(id_var, base_vars)]); a column that varies within a subject
  # contributes several rows for that subject, so the sampler silently
  # over-represents subjects with more distinct values. Warn rather than stop:
  # the run is still meaningful, but the baseline distribution is not the one
  # the user intended. Implemented with base tapply to avoid adding data.table
  # NSE columns to utils::globalVariables().
  bvars <- setdiff(base_vars, c(id_var, time_var))
  if (length(bvars) > 0) {
    varying <- bvars[vapply(bvars, function(v) {
      any(tapply(data[[v]], data[[id_var]],
                 function(x) length(unique(x))) > 1L)
    }, logical(1))]
    if (length(varying) > 0) {
      warning(sprintf(paste0(
        "Baseline variable(s) {%s} are not time-fixed: they take more than one ",
        "value within at least one '%s'. The Monte Carlo cohort is sampled from ",
        "the unique (%s, baseline) rows, so those subjects are over-represented ",
        "in every intervention. Keep only time-fixed columns in `base_vars` and model ",
        "time-varying ones with spec_model(mod_type = \"covariate\")."),
        paste(varying, collapse = ", "), id_var, id_var),
        call. = FALSE, domain = "causalMed")
    }
  }

  # Check if variables in the formula included in the data
  vars_models <- unlist(lapply(models, function(x) all.vars(formula(x$call))))
  check_var_in(vars_models, data)

  # Get the variable name of exposure and check if it the same as the exposure
  exp_flag <- mod_types == "exposure"
  if (sum(exp_flag) > 1) {
    stop("Only one exposure model is allowed.", domain = "causalMed")
  }

  if (any(exp_flag)) {
    exposure_var <- all.vars(formula(models[[which(exp_flag)]]$call)[[2]])
    if (exposure_var != exposure) {
      stop("The given exposure variable was different between exposure model in `models`.",
        domain = "causalMed"
      )
    }
  }

  # Check for survival models
  if (sum(mod_types == "censor") > 1) {
    stop("Only one censor model is allowed.", domain = "causalMed")
  }
  # (A second "survival" model is already caught by the outcome/survival count
  # above; only the censor model needs its own uniqueness check.)

  is_survival <- any(mod_types %in% c("censor", "survival"))

  if (is_survival) {
    y <- as.numeric(na.omit(data[[outcome]]))
    eps <- 1e-8
    if (!all(abs(y - 0) < eps | abs(y - 1) < eps)) {
      stop("For survival, outcome must be binary {0,1}.")
    }
  }

  # A censor model implies a survival setting, so it cannot be paired with a
  # plain "outcome" model. (The reverse case -- neither an outcome nor a
  # survival model -- is already rejected above.)
  if (is_survival && any(mod_types == "outcome")) {
    stop("You cannot define the outcome model with survival/censor model at the same time.", domain = "causalMed")
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
}

#' Validate one static intervention vector
#'
#' One static exposure regime is a numeric or logical vector, with no
#' \code{NA}, of length 1 (the same value at every time point) or
#' \code{time_len} (one value per \strong{distinct} time point in the data, in
#' sorted order). Shared by
#' \code{\link{check_intervention}} (for \code{gformula()}) and
#' \code{\link{mediation}} (for its two exposure regimes).
#'
#' @param x The vector to check.
#' @param time_len Number of distinct time points in the data.
#' @param arg Argument name used in error messages.
#' @param values Optional vector of allowed values. \code{NULL} (the default)
#'   allows any numeric value.
#' @return \code{x} coerced to numeric.
#' @keywords internal
check_static_regime <- function(x, time_len, arg, values = NULL) {
  if (!is.numeric(x) && !is.logical(x)) {
    stop(sprintf("`%s` must be a numeric or logical vector.", arg),
         domain = "causalMed")
  }
  if (!(length(x) %in% c(1L, time_len))) {
    stop(sprintf(
      "`%s` must have length 1 or %d (one value per distinct time point); got length %d.",
      arg, as.integer(time_len), length(x)), domain = "causalMed")
  }
  x <- as.numeric(x)
  # Before the values check, so an NA gets its own message on both paths
  # (mediation()'s 0/1 check would otherwise report it as a bad value, and
  # gformula(), with values = NULL, would pass it into the simulation).
  if (anyNA(x)) {
    stop(sprintf("`%s` must not contain NA.", arg), domain = "causalMed")
  }
  if (!is.null(values) && !all(x %in% values)) {
    bad <- x[!x %in% values]
    stop(sprintf("`%s` must take values in {%s}; got {%s}.",
                 arg, paste(values, collapse = ", "),
                 paste(sort(unique(bad), na.last = TRUE), collapse = ", ")),
         domain = "causalMed")
  }
  x
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

  # At most one element may be NULL (the natural-course reference)
  interv_value <- sapply(intervention, is.null)
  if (sum(interv_value) > 1) {
    stop("Only one intervention element may be NULL (the natural course).",
         domain = "causalMed")
  }

  # Check for intervention is a named list
  if (is.null(names(intervention)) | any(nchar(names(intervention)) == 0)) {
    stop("Intervention must be a named list object", domain = "causalMed")
  }

  # Each static element must be numeric/logical of length 1 or time_len --
  # the same rule mediation() applies to its exposure regimes. NULL is the
  # natural course; dyn_int() objects apply a rule at every time step.
  # Positional, not by name: [[name]] returns the first match, so a
  # duplicate-named element would escape validation.
  nms <- names(intervention)
  for (i in seq_along(intervention)) {
    x <- intervention[[i]]
    if (is.null(x) || inherits(x, "causalMed_dynint")) next
    check_static_regime(x, time_len, arg = sprintf("intervention$%s", nms[i]))
  }

  # Check for reference intervention. `&&` (not `&`) throughout: with `&` both
  # sides are always evaluated, so a character `ref_int` reached the numeric
  # comparison and was silently compared as a string.
  if (length(ref_int) != 1L || is.na(ref_int)) {
    stop("Length of `ref_int` must be 1.", domain = "causalMed")
  }

  if (is.numeric(ref_int)) {
    if (ref_int < 0 || ref_int != round(ref_int)) {
      stop("Numeric `ref_int` must be a whole number: 0 for the natural course, or the position of an intervention.",
           domain = "causalMed")
    }
    if (ref_int > length(intervention)) {
      stop("`ref_int` must be less than the length of intervention.", domain = "causalMed")
    }
  } else if (is.character(ref_int)) {
    # "natural" is always allowed: gformula() creates the natural-course
    # intervention when the user asks for it and did not supply one.
    if (ref_int != "natural" && !ref_int %in% names(intervention)) {
      stop("`ref_int` must be included in the names of intervention.", domain = "causalMed")
    }
  } else {
    stop("`ref_int` must be a whole number or the name of an intervention.",
         domain = "causalMed")
  }
}


#' Check that a natural-course intervention can be simulated
#'
#' A natural-course intervention -- a \code{NULL} element of \code{intervention},
#' or \code{intervention = NULL} -- draws the exposure from its own fitted model
#' at each time step. Without a \code{mod_type = "exposure"} model there is
#' nothing to draw from, and the exposure column never enters the Monte Carlo
#' dataset, so the run fails deep inside \code{model.matrix()} with an opaque
#' "object '<exposure>' not found". Fail early with an actionable message.
#'
#' @param models List of model specifications from \code{\link{spec_model}}.
#' @param intervention The resolved (named, non-\code{NULL}) intervention list.
#' @param exposure Character scalar. Name of the exposure variable.
#' @keywords internal
check_natural_course <- function(models, intervention, exposure) {
  if (!any(vapply(intervention, is.null, logical(1)))) return(invisible(NULL))
  if (any(vapply(models, function(m) identical(m$mod_type, "exposure"), logical(1))))
    return(invisible(NULL))

  # Do NOT advise "use only static/dynamic interventions": with the default
  # ref_int = 0 a static-only list still gets a natural course added as the
  # reference, so following that advice hits this same error again. The way out
  # is to name one of the supplied interventions as the reference.
  stop(sprintf(paste0(
    "A natural-course intervention was requested but no exposure model was ",
    "supplied, so there is nothing to draw '%s' from. Either add ",
    "spec_model(%s ~ ..., mod_type = \"exposure\") to `models`, or -- if you ",
    "did not want a natural course -- set `ref_int` to the name of one of ",
    "your interventions, since the default ref_int = 0 asks for the natural ",
    "course as the reference."),
    exposure, exposure), call. = FALSE, domain = "causalMed")
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

  # Every mediator must come before the outcome/survival model.
  if (length(out_idx) > 0 && any(med_idx > min(out_idx))) {
    offenders <- med_idx[med_idx > min(out_idx)]
    warning(sprintf(
      paste0(
        "Temporal ordering: mediator model(s) at position(s) [%s] appear after ",
        "the outcome/survival model at position [%d]. ",
        "Mediators must be simulated before the outcome."
      ),
      paste(offenders, collapse = ", "), min(out_idx)
    ), call. = FALSE)
  }

  # Exposure model must come before the first mediator.
  if (length(exp_idx) > 0 && exp_idx > min(med_idx)) {
    warning(sprintf(
      paste0(
        "Temporal ordering: exposure model at position [%d] appears after ",
        "the first mediator model at position [%d]. ",
        "Exposure must be set before any mediator is simulated."
      ),
      exp_idx, min(med_idx)
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


#' Report covariates modelled as exposure-affected under natural effects
#'
#' The effects reported under \code{mediation_type = "N"} are those of Zheng &
#' van der Laan (2017): the mediator is drawn from its conditional distribution
#' under the other regime given each subject's own history, and their Lemma 1
#' identifies them under sequential randomization and positivity, whether or not
#' a covariate is exposure-affected. Reading them as \emph{individual-level}
#' natural effects, contrasts of each subject's own counterfactual mediator,
#' additionally requires a cross-world independence assumption that is not
#' expected to hold when an exposure-affected covariate also confounds the
#' mediator-outcome relationship (Avin, Shpitser & Pearl 2005; VanderWeele 2014;
#' VanderWeele & Tchetgen Tchetgen 2017). This
#' check reads the model formulas only: it reports covariate models that carry
#' the exposure on the right-hand side, i.e. covariates the user has modelled
#' as exposure-affected. It does not establish that such a covariate also
#' confounds the mediator-outcome relationship, and its silence does not
#' establish that no such confounder exists.
#'
#' @param models List of model specifications from \code{\link{spec_model}}.
#' @param exposure Character scalar. Name of the exposure variable.
#' @return Invisibly returns a character vector of the covariate (response)
#'   names whose model includes the exposure on the right-hand side; empty if
#'   none. The caller can persist this so downstream methods (e.g.
#'   \code{print.gformula}) can re-surface the identifiability caveat.
#' @keywords internal
check_natural_identifiability <- function(models, exposure) {
  cov_idx <- which(sapply(models, function(m) m$mod_type == "covariate"))
  if (length(cov_idx) == 0) return(invisible(character(0)))

  hit <- character(0)
  for (i in cov_idx) {
    f <- formula(models[[i]]$call)
    if (length(f) < 3) next  # one-sided formula, no RHS
    rhs_vars <- all.vars(f[[3]])
    if (exposure %in% rhs_vars) {
      lhs <- all.vars(f[[2]])[1]
      hit <- c(hit, lhs)
    }
  }

  if (length(hit) > 0) {
    warning(sprintf(
      paste0(
        "mediation_type = \"N\" requested, and covariate model(s) for {%s} include ",
        "the exposure '%s' on the right-hand side, i.e. they are modelled as ",
        "exposure-affected. The reported effects are those of Zheng & van der ",
        "Laan (2017), identified under sequential randomization and positivity ",
        "(their Lemma 1). Reading them as individual-level natural effects, ",
        "contrasts of each subject's own counterfactual mediator, additionally ",
        "requires a cross-world independence assumption that is not expected to ",
        "hold if such a covariate also confounds the mediator-outcome ",
        "relationship (Avin, Shpitser & Pearl 2005; VanderWeele 2014; ",
        "VanderWeele & Tchetgen Tchetgen 2017, who propose the randomized ",
        "interventional analogues, available here as mediation_type = \"I\"). ",
        "This check reads the model formulas only and cannot verify the causal ",
        "structure."
      ),
      paste(hit, collapse = ", "), exposure
    ), call. = FALSE, domain = "causalMed")
  }

  invisible(hit)
}


#' Check that the natural-effect mediator can be evaluated on the other regime
#'
#' Under \code{mediation_type = "N"} the cross-world mediator is drawn from its
#' model evaluated on the intervention's own covariate history with the
#' exposure history set to the other regime (Zheng & van der Laan 2017, Eq. 5).
#' The Monte Carlo engine sets exactly two kinds of input to that regime: the
#' exposure itself (so exposure terms written in the mediator formula, e.g.
#' \code{A:L}, are evaluated on it) and first-order exposure lags (an
#' \code{in_recode} entry that only copies the exposure, e.g.
#' \code{recodes(lag_A = A)}). Every other column a recode derives from the
#' exposure -- a chained lag, a cumulative count or other expression, an
#' \code{out_recode} copy, a column created by a model's own \code{recode} --
#' keeps the intervention's own exposure history, and the mediator
#' \code{subset} is evaluated on the intervention's own data. The check fails
#' closed: a mediator formula that reads such a column, or a mediator
#' \code{subset} that reads the exposure or anything derived from it, is
#' rejected.
#'
#' @param models List of model specifications from \code{\link{spec_model}}.
#' @param exposure Character scalar. Name of the exposure variable.
#' @param init_recode,in_recode,out_recode The recode hooks passed to
#'   \code{\link{mediation}}.
#' @return Invisibly, the names of the first-order exposure lag columns.
#' @keywords internal
check_natural_exposure_history <- function(models, exposure, init_recode = NULL,
                                           in_recode = NULL, out_recode = NULL) {
  lags <- exposure_lag_cols(in_recode, exposure)

  # Every recode assignment in the run as (target, columns read), leaving out
  # the in_recode entries that only copy the exposure: the swap sets those.
  # Skipped by POSITION, not by name, so a second entry for the same column
  # that transforms it (lag1_A = 2 * lag1_A) stays in and marks it derived.
  assignments <- function(rc, skip = NULL) {
    if (length(rc) == 0L) return(list())
    keep <- if (is.null(skip)) rep(TRUE, length(rc)) else !skip
    Map(function(nm, ex) list(target = nm, reads = all.vars(ex)),
        names(rc)[keep], as.list(rc)[keep])
  }
  is_copy <- vapply(in_recode, function(ex)
    is.name(ex) && identical(as.character(ex), exposure), logical(1))
  pairs <- c(assignments(init_recode),
             assignments(in_recode, skip = is_copy),
             assignments(out_recode),
             unlist(lapply(models, function(m) assignments(m$recode)),
                    recursive = FALSE))

  # Columns derived from the exposure history through a recode the swap does
  # not set, followed to a fixed point (a chained lag reads a lag, and so on).
  derived <- character(0)
  repeat {
    src <- c(exposure, lags, derived)
    hit <- vapply(pairs, function(p) any(p$reads %in% src), logical(1))
    new <- setdiff(vapply(pairs[hit], `[[`, character(1), "target"),
                   c(exposure, derived))
    if (length(new) == 0L) break
    derived <- c(derived, new)
  }

  med_idx <- which(vapply(models, function(m) identical(m$mod_type, "mediator"),
                          logical(1)))
  for (i in med_idx) {
    m   <- models[[i]]
    f   <- formula(m$call)
    rhs <- if (length(f) >= 3L) all.vars(f[[3L]]) else character(0)
    sub <- if (is.null(m$subset)) character(0) else all.vars(m$subset)
    # The formula is evaluated on the swapped exposure and lags, so only other
    # derived columns are a problem there. The subset is evaluated on the
    # intervention's own data, so ANY exposure dependence is.
    bad <- list(formula = intersect(rhs, derived),
                subset  = intersect(sub, c(exposure, lags, derived)))
    bad <- bad[lengths(bad) > 0L]
    if (length(bad) > 0L) {
      where <- paste(sprintf("its %s reads {%s}", names(bad),
                             vapply(bad, paste, character(1), collapse = ", ")),
                     collapse = " and ")
      stop(sprintf(paste0(
        "mediation_type = \"N\" draws the cross-world mediator '%s' from its ",
        "model with the exposure history set to the other regime (Zheng & van ",
        "der Laan 2017, Eq. 5), but %s. The simulation sets only the exposure ",
        "'%s' itself and first-order exposure lags (an in_recode entry that ",
        "only copies the exposure, e.g. recodes(lag1_%s = %s)) to that regime; ",
        "exposure terms written in the formula (e.g. %s:L) are evaluated on ",
        "them. Other columns derived from the exposure (chained lags, ",
        "cumulative or other expressions, out_recode copies, columns created by ",
        "a model's own recode) and the model's subset are evaluated on the ",
        "intervention's own exposure history."),
        all.vars(f[[2L]])[1L], where, exposure, exposure, exposure, exposure),
        call. = FALSE, domain = "causalMed")
    }
  }

  invisible(lags)
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
#' @keywords internal
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

