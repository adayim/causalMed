
# Internal helper --------------------------------------------------------------
# Format a data.table for display:
#   - round numeric columns to `digits`
#   - rename bootstrap CI columns to human-readable labels
# Uses dt[] to force a visible data.table return.
.fmt_for_print <- function(x, digits) {
  dt <- data.table::copy(data.table::as.data.table(x))
  # Round numeric columns in-place
  num_cols <- names(dt)[vapply(dt, is.numeric, logical(1L))]
  if (length(num_cols))
    dt[, (num_cols) := lapply(.SD, round, digits), .SDcols = num_cols]
  # Rename RD-scale CI columns
  old_rd <- c("perct_lcl",      "perct_ucl",       "norm_lcl",       "norm_ucl")
  new_rd <- c("RD 2.5%(pct)",   "RD 97.5%(pct)",   "RD 2.5%(norm)",  "RD 97.5%(norm)")
  present_rd <- old_rd %in% names(dt)
  if (any(present_rd))
    data.table::setnames(dt, old_rd[present_rd], new_rd[present_rd])
  # Rename RR-scale CI columns (mediation output only)
  old_rr <- c("perct_lcl_RR",   "perct_ucl_RR",    "norm_lcl_RR",    "norm_ucl_RR")
  new_rr <- c("RR 2.5%(pct)",   "RR 97.5%(pct)",   "RR 2.5%(norm)",  "RR 97.5%(norm)")
  present_rr <- old_rr %in% names(dt)
  if (any(present_rr))
    data.table::setnames(dt, old_rr[present_rr], new_rr[present_rr])
  return(dt[])   # dt[] forces a visible data.table return
}


#' Print
#'
#' Print method for objects returned by \code{\link{gformula}} or \code{\link{mediation}}.
#'
#' @param x Object of class \code{"gformula"}.
#' @param models Logical. If \code{TRUE}, print fitted model details (call and
#'   coefficients, or full summary when \code{return_fitted = TRUE}).
#'   Default \code{FALSE}.
#' @param digits Integer. Number of significant digits for numeric output.
#'   Default \code{max(3, getOption("digits") - 3)}.
#' @param ... Not used.
#'
#' @export
print.gformula <- function(x,
                           models = FALSE,
                           digits = max(3, getOption("digits") - 3),
                           ...) {

  is_mediation <- !is.null(x$all.args$mediation_type)
  args         <- x$all.args

  # -- Call ------------------------------------------------------------------
  cat("Call:\n")
  print(x$call)

  # -- Analysis metadata -----------------------------------------------------
  cat("\n--- Analysis setup ---\n")
  cat(sprintf("  Exposure     : %s\n", args$exposure))
  cat(sprintf("  Time variable: %s\n", args$time_var))
  cat(sprintf("  ID variable  : %s\n", args$id_var))
  if (!is.null(args$base_vars) && length(args$base_vars) > 0)
    cat(sprintf("  Baseline vars: %s\n", paste(args$base_vars, collapse = ", ")))

  seed_str <- if (is.null(args$seed)) "none" else as.character(args$seed)
  boot_str <- if (args$R <= 1L) "none" else as.character(args$R)
  cat(sprintf("  MC sample    : %d\n", args$mc_sample))
  cat(sprintf("  Bootstrap R  : %s\n", boot_str))
  cat(sprintf("  Seed         : %s\n", seed_str))

  if (is_mediation) {
    type_label <- if (args$mediation_type == "I")
      "Interventional effects (IDE/IIE) -- Lin et al. (2017)"
    else
      "Natural effects (NDE/NIE) -- Zheng & van der Laan (2017)"
    cat(sprintf("  Mediation    : %s\n", type_label))
  } else {
    # show reference intervention for gformula
    if (!is.null(args$ref_int))
      cat(sprintf("  Reference    : %s\n", as.character(args$ref_int)))
  }

  # -- Section 1: arm-level mean outcomes ------------------------------------
  if (is_mediation) {
    s1_hdr  <- "\n--- Marginal mean outcome per arm (Q-functionals) ---"
    s1_note <- paste0(
      "  Ph11  = E[Y(a=1, M(1))]:  exposure=1, mediator under a=1\n",
      "  Phi10 = E[Y(a=1, M(0))]:  exposure=1, mediator under a=0  [cross-world]\n",
      "  Phi00 = E[Y(a=0, M(0))]:  exposure=0, mediator under a=0\n"
    )
  } else {
    s1_hdr  <- "\n--- Mean outcome by intervention ---"
    s1_note <- NULL
  }
  cat(s1_hdr, "\n")
  if (!is.null(s1_note)) cat(s1_note)
  if (!is.null(x$effect_size)) {
    print(.fmt_for_print(x$effect_size, digits))
  }

  # -- Section 2: contrasts / decomposition ----------------------------------
  if (is_mediation) {
    s2_hdr  <- "\n--- Effect decomposition ---"
    s2_note <- paste0(
      "  Total effect    = Ph11 - Phi00  =  E[Y(1,M(1))] - E[Y(0,M(0))]\n",
      "  Direct effect   = Phi10 - Phi00 =  E[Y(1,M(0))] - E[Y(0,M(0))]\n",
      "  Indirect effect = Ph11 - Phi10  =  E[Y(1,M(1))] - E[Y(1,M(0))]\n",
      "  Mediation Prop. = Indirect / Total  (as a percentage; RR not applicable)\n",
      "  RD = risk difference;  RR = risk ratio\n"
    )
  } else {
    s2_hdr  <- "\n--- Contrasts vs. reference intervention ---"
    s2_note <- NULL
  }
  cat(s2_hdr, "\n")
  if (!is.null(s2_note)) cat(s2_note)
  if (is.null(x$estimate)) {
    cat("  Note: contrasts require at least 2 interventions.\n")
  } else {
    print(.fmt_for_print(x$estimate, digits))
  }

  # -- Model details (optional) ----------------------------------------------
  if (models) {
    cat("\n--- Fitted models ---\n")
    for (nm in names(x$fitted_models)) {
      fitmod   <- x$fitted_models[[nm]]
      attr_lst <- attributes(fitmod)
      cat(sprintf("\nVariable : %s  [type: %s | role: %s]\n",
                  nm, attr_lst$var_type, attr_lst$mod_type))
      cat(paste(rep("-", 40L), collapse = ""), "\n")
      if (!is.null(attr_lst$recodes)) {
        cat("Recodes: ")
        print(attr_lst$recodes)
      }
      if (args$return_fitted) {
        print(summary(fitmod))
      } else {
        cat("Call: ")
        print(fitmod$call)
        cat("Coefficients:\n")
        print(fitmod$coeff)
      }
    }
  }

  invisible(x)
}


#' Summary
#'
#' Summary method for objects returned by \code{\link{gformula}} or
#' \code{\link{mediation}}. Shows the full estimation results followed by
#' fitted model coefficient tables.
#'
#' @param object Object of class \code{"gformula"}.
#' @param digits Integer. Number of significant digits.
#'   Default \code{max(3, getOption("digits") - 3)}.
#' @param ... Not used.
#'
#' @export
summary.gformula <- function(object,
                             digits = max(3, getOption("digits") - 3),
                             ...) {
  # Print estimation results first
  print(object, models = FALSE, digits = digits)

  # Then append full model coefficient tables
  cat("\n--- Fitted model coefficients ---\n")
  for (nm in names(object$fitted_models)) {
    fitmod   <- object$fitted_models[[nm]]
    attr_lst <- attributes(fitmod)
    cat(sprintf("\nVariable : %s  [type: %s | role: %s]\n",
                nm, attr_lst$var_type, attr_lst$mod_type))
    cat(paste(rep("-", 40L), collapse = ""), "\n")

    coef_tbl <- if (object$all.args$return_fitted) {
      tryCatch(summary(fitmod)$coefficients, error = function(e) NULL)
    } else {
      fitmod$coeff
    }

    if (!is.null(coef_tbl)) {
      print(round(coef_tbl, digits))
    } else {
      cat("  (rerun with return_fitted = TRUE for full coefficient table)\n")
    }
  }

  invisible(object)
}
