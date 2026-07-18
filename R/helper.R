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
      list(call  = x$fitted$call,
           coeff = summary(x$fitted)$coefficients)
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

# Summarize the input data for display by print.gformula():
# number of individuals, observations, and time points.
summarize_input_data <- function(data, id_var, time_var) {
  times <- sort(unique(na.omit(data[[time_var]])))
  list(
    n_id    = data.table::uniqueN(data[[id_var]]),
    n_obs   = nrow(data),
    n_times = length(times),
    t_min   = min(times),
    t_max   = max(times)
  )
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
