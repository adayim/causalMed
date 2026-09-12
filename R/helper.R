#' Define parameters for recoding
#' 
#' @description Capture expressions for variable assignment.
#' @param ... Named expressions (e.g., daysq = day^2).
#' @return A list of expressions with class 'causalMed_recodes'.
#' @export
recodes <- function(...) {
  # 1. Capture arguments as unevaluated expressions
  # alist() ensures we get 'day^2' without trying to run it
  exprs <- eval(substitute(alist(...)))

  class(exprs) <- c("causalMed_recodes", "list")
  
  return(exprs)
}

#' Define a dynamic intervention rule
#'
#' @description Captures an R expression for a dynamic (rule-based) exposure
#'   intervention. The expression is evaluated at each time step inside the
#'   simulated \code{data.table}, so any column in the current Monte Carlo
#'   dataset can be referenced by name (including the exposure itself, which
#'   holds its natural-course draw at the time of evaluation).
#'
#' @param expr An R expression that returns a numeric (or logical coerced to
#'   numeric) vector with one value per row. Column names from the current
#'   simulated dataset are in scope (e.g. \code{L1}, \code{A}).
#'
#' @return An object of class \code{"causalMed_dynint"}.
#'
#' @examples
#' # Treat only if the natural-course value of A exceeds 0:
#' dyn_int(as.numeric(A > 0))
#'
#' # Treat if L1 > 0 or L2 equals 1:
#' dyn_int(as.numeric(L1 > 0 | L2 == 1))
#'
#' @export
dyn_int <- function(expr) {
  structure(list(substitute(expr)), class = "causalMed_dynint")
}

#' @import data.table
apply_recodes <- function(data, recode_params) {
  
  # Handle Custom Class (from recodes())
  for (i in seq_along(recode_params)) {
    var_name <- names(recode_params)[i]
    expr <- recode_params[[i]]
    
    # Evaluate expression within the data.table environment
    # 'var_name' is wrapped in parentheses to use the string as the column name
    data[, (var_name) := eval(expr)]
  }
  
  invisible(data)
}

# The in_recode entries that only copy the exposure (recodes(lag_A = A)): the
# first-order exposure lags. in_recode runs at the start of every step after
# the first, before the exposure is set, so at step k >= 2 such a column holds
# the exposure of step k - 1. (An out_recode copy is NOT one: out_recode also
# skips the first step, so it still holds its init value at step 2.) Uses the
# TMLE engine's lag map so both "N" estimators share one definition.
exposure_lag_cols <- function(in_recode, exposure) {
  lag_map <- .tmle_lag_map(in_recode)
  as.character(names(lag_map)[vapply(lag_map, identical, logical(1), exposure)])
}

# Which columns carry exposure history that a regime evaluation does NOT set,
# followed to a fixed point (a chained lag reads a lag, and so on). `pairs` is
# the run's recode assignments as list(target =, reads =); `src` the columns
# the regime DOES set (the exposure and its first-order lags). The exposure
# itself is never reported: both estimators always set it. A lag column CAN
# come back, when a later recode entry transforms it rather than copying the
# exposure. Shared by the two "N" guards -- check_natural_exposure_history()
# and .tmle_check_model_recodes() -- so gcomp and TMLE cannot drift apart on
# what counts as exposure-derived.
exposure_derived_cols <- function(pairs, src, exposure) {
  derived <- character(0)
  repeat {
    hit <- vapply(pairs, function(p) any(p$reads %in% c(src, derived)),
                  logical(1))
    new <- setdiff(vapply(pairs[hit], `[[`, character(1), "target"),
                   c(exposure, derived))
    if (length(new) == 0L) break
    derived <- c(derived, new)
  }
  derived
}

# TRUE for the default regime pair, always (1) vs never (0) exposed. An object
# carrying NEITHER regime -- saved before the arguments existed -- is the
# default. Both are tested explicitly for length: all() is TRUE on a
# zero-length vector, so a half-populated argument list would otherwise pass as
# the default and gate estimator = "tmle" open on an arbitrary regime pair.
default_regime_pair <- function(exposure_regime, reference_regime) {
  if (is.null(exposure_regime) && is.null(reference_regime)) return(TRUE)
  length(exposure_regime) > 0L && length(reference_regime) > 0L &&
    all(exposure_regime == 1) && all(reference_regime == 0)
}

#' Derive parameters for a function from the current environment
#'
#' @description
#' Extracts the required parameters for a given function from the current environment.
#'
#' @param fun The function whose parameters are to be extracted.
#' @param env The environment from which to extract parameters (default is the parent environment).
#' @param ... Additional arguments passed to the function.
#' @param dots Optional argument for handling extra dots.
#'
#' @return A list of arguments for the specified function.
#' 
#' @source \url{https://stackoverflow.com/a/51002887}
#'
#' @keywords internal

get_args_for <- function(fun, env = parent.frame(), ..., dots = NULL) {
  potential <- names(formals(fun))

  if ("..." %in% potential) {
    if (missing(dots)) {
      # return everything from parent frame
      return(as.list(env))
    } else if (!is.list(dots)) {
      stop("If provided, 'dots' should be a list.")
    }

    potential <- setdiff(potential, "...")
  }

  # get all formal arguments that can be found in parent frame
  args <- mget(potential, env, ..., ifnotfound = list(NULL), inherits = FALSE)
  # remove not found
  args <- args[vapply(args, Negate(is.null), logical(1))]
  # return found args and dots
  c(args, dots)
}

# Internal function to intialise the warning
init_warn <- function(){
  causalmed_env$warning <- c()
}

# Collect warnings, muffling them so they do not also print immediately.
# All collected warnings are summarised (with repeat counts) at function exit
# via emit_warnings().
run_withwarning_collect <- function(expr, msg) {
  withCallingHandlers(
    expr,
    warning = function(w) {
      warn_messg <- sprintf("%s;\nWarning: %s",
                            msg,
                            paste(w$message, collapse = "\n"))
      causalmed_env$warning <- c(causalmed_env$warning, warn_messg)
      invokeRestart("muffleWarning")
    }
  )
}

# Emit a deduplicated summary of all collected warnings.
# Identical messages are counted and shown once with a repeat count,
# preventing hundreds of identical lines when the same warning fires
# in every bootstrap replicate.
emit_warnings <- function() {
  if (length(causalmed_env$warning) == 0) return(invisible(NULL))
  warn_tbl  <- table(causalmed_env$warning)
  warn_msgs <- mapply(function(msg, n) {
    if (n > 1L) paste0(msg, "\n  [repeated ", n, " time(s)]") else msg
  }, names(warn_tbl), as.integer(warn_tbl), SIMPLIFY = TRUE)
  message(paste(warn_msgs, collapse = "\n=============\n"), domain = "causalMed")
  invisible(NULL)
}



# Wrap fitted models for the return object of gformula()/mediation(): the
# full model objects when return_fitted = TRUE, otherwise a compact list
# (call + coefficient table). Either way, each element carries the
# recodes/subset/var_type/mod_type attributes that print.gformula(models =
# TRUE) and summary.gformula() use to label the models, and the list is named
# by each model's response variable. Every return path must use this wrapper
# so those labels stay available.
wrap_fitted_models <- function(fitted_models, return_fitted) {
  out <- lapply(fitted_models, function(x) {
    r <- if (return_fitted) {
      x$fitted
    } else {
      # The compact summary is a convenience for printing, so it must never
      # take down an otherwise successful run. A `custom_fit` model can be any
      # class: `$` fails on a non-list fit, and summary() may return something
      # with no $coefficients (summary.default gives a matrix, and `$` on a
      # matrix is an error) or a non-numeric table that print(round(.)) would
      # then reject. Degrade to NULL; print()/summary() already say
      # "rerun with return_fitted = TRUE" when the table is missing.
      list(call  = tryCatch(x$fitted$call, error = function(e) NULL),
           coeff = tryCatch({
             cf <- summary(x$fitted)$coefficients
             if (is.null(cf) || !is.numeric(cf)) NULL else cf
           }, error = function(e) NULL))
    }
    structure(r,
              recodes  = x$recodes,
              subset   = x$subset,
              var_type = x$var_type,
              mod_type = x$mod_type)
  })
  names(out) <- vapply(fitted_models, function(x) x$rsp_vars, character(1))
  out
}

# The simulation grid: the distinct observed values of `time_var`, sorted.
# gformula() and mediation() compute it ONCE from the input data and pass it
# down, so every pass -- bootstrap replicates included -- simulates the same
# steps. Every consumer goes through this one definition.
time_grid <- function(data, time_var) sort(unique(na.omit(data[[time_var]])))

# Summarize the input data for display by print.gformula():
# number of individuals, observations, and time points. `time_seq` is the
# caller's simulation grid (time_grid()), stored so print() can label a
# regime's positions.
summarize_input_data <- function(data, id_var, time_seq) {
  list(
    n_id     = data.table::uniqueN(data[[id_var]]),
    n_obs    = nrow(data),
    n_times  = length(time_seq),
    t_min    = min(time_seq),
    t_max    = max(time_seq),
    time_seq = time_seq
  )
}

# How many observed subjects follow each exposure regime: a subject follows a
# regime when its observed exposure equals the regime's value at EVERY time
# point at which that subject is observed (matched by time value, so a short
# follow-up is compared on the times it has; an NA exposure at any observed
# time means the subject follows neither; a row with an NA time has no grid
# position and is ignored). A count for the data summary -- nothing in the
# estimation uses it. Base R only (tapply), so no data.table NSE columns to
# register in zzz.R.
#
# `n_following` alone overstates support under right-censoring: a subject
# observed only at t = 0 matches the first element of every regime that
# shares that value, so it is counted as "following" regimes it was never
# followed past the first step on -- and the SAME subject can be counted for
# more than one regime this way. `n_complete` restricts the count to
# subjects who also have an observation at every time point in `time_seq`,
# so it reports how many were actually followed for the whole grid.
#
# `regimes` is a named list of full-length regime vectors aligned to
# `time_seq`. Returns a data.frame with one row per regime: `regime`,
# `values`, `n_following`, `prop_following`, `n_complete`.
regime_support <- function(data, id_var, time_var, exposure, time_seq, regimes) {
  # as.character(): tapply() groups by factor(ids), and an unused factor level
  # would give an empty group whose all() is NA, poisoning the count.
  ids_all <- as.character(data[[id_var]])
  # tapply() groups by factor(ids) and DISCARDS an NA id, so an NA-id subject
  # can never reach the numerator; counting it in the denominator would
  # understate every regime's support. Drop it from both.
  n_id    <- data.table::uniqueN(ids_all[!is.na(ids_all)])
  pos     <- match(data[[time_var]], time_seq)    # NA for an NA time
  # Compare subjects only on rows that sit on the grid. Treating an off-grid
  # row as a mismatch would disqualify the subject from every regime.
  on_grid <- !is.na(pos)
  ids <- ids_all[on_grid]
  pos <- pos[on_grid]
  aa  <- data[[exposure]][on_grid]
  # Subjects observed at every grid time. Counted once here rather than per
  # regime: with right-censored data most subjects have short follow-up, and a
  # subject observed only at t = 0 matches the first element of MANY regimes,
  # so `n_following` alone reads far larger than the number of subjects who
  # actually followed the whole trajectory. DISTINCT grid times, so a
  # duplicated (id, time) row cannot stand in for a missing time point.
  n_obs_grid <- tapply(pos, ids, function(p) length(unique(p)))
  complete   <- n_obs_grid == length(time_seq)
  rows <- lapply(regimes, function(reg) {
    follows <- !is.na(aa) & aa == reg[pos]
    ok <- tapply(follows, ids, all)
    n  <- sum(ok)
    data.frame(values         = paste(reg, collapse = " "),
               n_following    = as.integer(n),
               prop_following = n / n_id,
               n_complete     = as.integer(sum(ok & complete[names(ok)])),
               stringsAsFactors = FALSE)
  })
  out <- cbind(data.frame(regime = names(regimes), stringsAsFactors = FALSE),
               do.call(rbind, rows))
  rownames(out) <- NULL
  out
}

# Observed nonparametric benchmark of the outcome, printed alongside the
# simulated intervention means as an informal model check.
#   - Survival outcomes: product-limit (Kaplan-Meier-type) cumulative
#     incidence by the last time point. Long format with post-event rows
#     removed means the per-time at-risk set shrinks with events and
#     censoring, so hazard = events / at-risk handles right-censoring.
#   - Non-survival outcomes: observed mean of the outcome at the last
#     time point.
# Implemented with base tapply (no data.table NSE) to avoid adding columns
# to utils::globalVariables().
observed_benchmark <- function(data, outcome, time_var, is_survival) {
  tt <- data[[time_var]]
  yy <- data[[outcome]]
  ok <- !is.na(tt) & !is.na(yy)
  tt <- tt[ok]
  yy <- yy[ok]
  if (length(yy) == 0L) {
    return(list(value = NA_real_, label = sprintf("mean of %s", outcome)))
  }
  if (is_survival) {
    events  <- tapply(yy, tt, sum)
    at_risk <- tapply(yy, tt, length)
    value   <- 1 - prod(1 - events / at_risk)
    label   <- sprintf("cumulative incidence of %s by t = %s (product-limit)",
                       outcome, max(tt))
  } else {
    last  <- tt == max(tt)
    value <- mean(yy[last])
    label <- sprintf("mean of %s at t = %s (end of follow-up)",
                     outcome, max(tt))
  }
  list(value = value, label = label)
}
