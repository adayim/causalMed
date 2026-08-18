
# Internal helper --------------------------------------------------------------
# Format a data.table for display:
#   - round numeric columns to `digits` decimal places
#   - rename bootstrap CI columns to human-readable labels
#
# `ci_labels` controls how CI columns are labelled:
#   - "plain": the CI columns describe whatever single estimate column the
#     table holds (per-intervention means, or the stacked Difference/Ratio
#     contrast table where the scale is given per row by Risk_type), so no
#     scale prefix is attached.
#   - "rd_rr": the mediation decomposition table, where perct_lcl/norm_lcl
#     belong to the RD column and *_RR to the RR column.
# Uses dt[] to force a visible data.table return.
.fmt_for_print <- function(x, digits, ci_labels = c("plain", "rd_rr")) {
  ci_labels <- match.arg(ci_labels)
  dt <- data.table::copy(data.table::as.data.table(x))
  # Round numeric columns in-place
  num_cols <- names(dt)[vapply(dt, is.numeric, logical(1L))]
  if (length(num_cols))
    dt[, (num_cols) := lapply(.SD, round, digits), .SDcols = num_cols]
  if (ci_labels == "rd_rr") {
    # RD-scale columns (Sd/CIs computed from the RD column)
    old_rd <- c("Sd",       "perct_lcl",     "perct_ucl",      "norm_lcl",      "norm_ucl")
    new_rd <- c("Sd(RD)",   "RD 2.5%(pct)",  "RD 97.5%(pct)",  "RD 2.5%(norm)", "RD 97.5%(norm)")
    # RR-scale columns
    old_rr <- c("Sd_RR",    "perct_lcl_RR",  "perct_ucl_RR",   "norm_lcl_RR",   "norm_ucl_RR")
    new_rr <- c("Sd(RR)",   "RR 2.5%(pct)",  "RR 97.5%(pct)",  "RR 2.5%(norm)", "RR 97.5%(norm)")
    old_all <- c(old_rd, old_rr)
    new_all <- c(new_rd, new_rr)
  } else {
    # Single-scale table: CIs describe the estimate column itself.
    old_all <- c("perct_lcl",  "perct_ucl",   "norm_lcl",   "norm_ucl")
    new_all <- c("2.5%(pct)",  "97.5%(pct)",  "2.5%(norm)", "97.5%(norm)")
  }
  present <- old_all %in% names(dt)
  if (any(present))
    data.table::setnames(dt, old_all[present], new_all[present])
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
#' @param digits Integer. Number of decimal places used when rounding
#'   numeric output. Default \code{max(3, getOption("digits") - 3)}.
#' @param ... Not used.
#'
#' @export
print.gformula <- function(x,
                           models = FALSE,
                           digits = max(3, getOption("digits") - 3),
                           ...) {

  is_mediation <- !is.null(x$all.args$mediation_type)
  args         <- x$all.args
  is_interv    <- is_mediation && isTRUE(args$mediation_type == "I")

  # Model roles from the fitted_models attributes (set by both estimators).
  mod_types <- vapply(x$fitted_models, function(m) {
    mt <- attr(m, "mod_type", exact = TRUE)
    if (is.null(mt)) NA_character_ else as.character(mt)
  }, character(1))
  med_vars    <- names(mod_types)[mod_types %in% "mediator"]
  outcome_nm  <- names(mod_types)[mod_types %in% c("outcome", "survival")][1]
  is_survival <- any(mod_types %in% c("survival", "censor"), na.rm = TRUE)

  # Data summary (absent on objects created by older package versions).
  ds <- x$data_summary

  # -- Call ------------------------------------------------------------------
  cat("Call:\n")
  print(x$call)

  # -- Analysis metadata -----------------------------------------------------
  cat("\n--- Analysis setup ---\n")
  cat(sprintf("  Exposure     : %s\n", args$exposure))
  if (is_mediation && length(med_vars) > 0)
    cat(sprintf("  Mediator(s)  : %s\n", paste(med_vars, collapse = ", ")))
  if (length(outcome_nm) == 1L && !is.na(outcome_nm)) {
    out_lab <- if (!is.null(ds)) {
      if (is_survival)
        sprintf("%s  [survival: cumulative incidence by t = %s]",
                outcome_nm, ds$t_max)
      else
        sprintf("%s  [mean outcome at t = %s, end of follow-up]",
                outcome_nm, ds$t_max)
    } else {
      outcome_nm
    }
    cat(sprintf("  Outcome      : %s\n", out_lab))
  }
  if (!is.null(ds)) {
    cat(sprintf("  Time variable: %s  (%d time points: %s ... %s)\n",
                args$time_var, ds$n_times, ds$t_min, ds$t_max))
  } else {
    cat(sprintf("  Time variable: %s\n", args$time_var))
  }
  cat(sprintf("  ID variable  : %s\n", args$id_var))
  if (!is.null(args$base_vars) && length(args$base_vars) > 0)
    cat(sprintf("  Baseline vars: %s\n", paste(args$base_vars, collapse = ", ")))
  if (!is.null(ds))
    cat(sprintf("  Data         : %s individuals, %s observations\n",
                format(ds$n_id, big.mark = ","),
                format(ds$n_obs, big.mark = ",")))

  is_tmle <- is_mediation && identical(args$estimator, "tmle")

  seed_str <- if (is.null(args$seed)) "none" else as.character(args$seed)
  boot_str <- if (args$R <= 1L) "none" else as.character(args$R)
  if (is_tmle) {
    cat("  Estimator    : TMLE (targeted maximum likelihood; Zheng & van der Laan 2017)\n")
  } else {
    cat(sprintf("  MC sample    : %d\n", args$mc_sample))
  }
  cat(sprintf("  Bootstrap R  : %s\n", boot_str))
  if (is_interv && !is.null(args$n_vw))
    cat(sprintf("  n_vw         : %d  (permutation draws averaged per cross-world intervention)\n",
                as.integer(args$n_vw)))
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

  # -- Section 1: intervention-level mean outcomes ------------------------------------
  if (is_mediation) {
    s1_hdr  <- "\n--- Marginal mean outcome per intervention ---"
    if (is_interv) {
      s1_note <- paste0(
        "  Under interventional effects, each intervention draws its mediators from independently-permuted pools (G):\n",
        "  Phi11 = E[Y(a=1, G1)]:  exposure=1, mediators ~ a=1 pool  [reference]\n",
        "  Phi10 = E[Y(a=1, G0)]:  exposure=1, mediators ~ a=0 pool  [cross-world]\n",
        if (length(med_vars) > 1)
          "  Phi1_k:  exposure=1, first k mediators ~ a=1 pool, rest ~ a=0 pool  [sequential]\n",
        "  Phi00 = E[Y(a=0, G0)]:  exposure=0, mediators ~ a=0 pool  [reference]\n",
        "  nat1/nat0 = E[Y(a=1)]/E[Y(a=0)]:  natural course (used for the total effect)\n"
      )
    } else {
      s1_note <- paste0(
        "  Phi11 = E[Y(a=1, M(1))]:  exposure=1, mediator under a=1\n",
        "  Phi10 = E[Y(a=1, M(0))]:  exposure=1, mediator under a=0  [cross-world]\n",
        "  Phi00 = E[Y(a=0, M(0))]:  exposure=0, mediator under a=0\n"
      )
    }
  } else {
    s1_hdr  <- "\n--- Mean outcome by intervention ---"
    s1_note <- NULL
  }
  cat(s1_hdr, "\n")
  if (!is.null(s1_note)) cat(s1_note)
  if (!is.null(x$effect_size)) {
    # CI columns here describe the per-intervention mean Est (not RD/RR).
    print(.fmt_for_print(x$effect_size, digits, ci_labels = "plain"))
  }

  # Observed nonparametric benchmark, as an informal model check.
  if (!is.null(x$observed) && is.finite(x$observed$value)) {
    cat(sprintf("  Observed (nonparametric) %s: %s\n",
                x$observed$label, format(round(x$observed$value, digits))))
    if (is_mediation) {
      # No natural-course intervention exists here (all interventions fix the
      # exposure), so this is a rough benchmark, not a direct comparison.
      cat("  (informal benchmark; interventions fix the exposure, so exact agreement is not expected)\n")
    } else {
      cat("  (informal model check: compare with the natural-course intervention)\n")
    }
  }

  # -- Section 2: contrasts / decomposition ----------------------------------
  if (is_mediation) {
    s2_hdr  <- "\n--- Effect decomposition ---"
    if (is_interv) {
      s2_note <- paste0(
        "  Direct effect (IDE)   = Phi10 - Phi00\n",
        "  Indirect effect (IIE) = Phi11 - Phi10   (sequential per mediator when N>=2)\n",
        "  IDE + IIE             = Phi11 - Phi00    (interventional overall effect)\n",
        "  Total effect (TE)     = nat1 - nat0      (natural plug-in g-formula)\n",
        "  TE - (Direct+Indirect)= mediated-interaction residual (TE - overall)\n",
        "  Mediation Prop.       = (Total - Direct) / Total  (percentage; RR not applicable)\n",
        "  RD = risk difference;  RR = risk ratio\n"
      )
    } else {
      s2_note <- paste0(
        "  Total effect    = Phi11 - Phi00 =  E[Y(1,M(1))] - E[Y(0,M(0))]\n",
        "  Direct effect   = Phi10 - Phi00 =  E[Y(1,M(0))] - E[Y(0,M(0))]\n",
        "  Indirect effect = Phi11 - Phi10 =  E[Y(1,M(1))] - E[Y(1,M(0))]\n",
        "  Mediation Prop. = Indirect / Total  (as a percentage; RR not applicable)\n",
        "  RD = risk difference;  RR = risk ratio\n"
      )
    }
  } else {
    s2_hdr  <- "\n--- Contrasts vs. reference intervention ---"
    s2_note <- NULL
  }
  cat(s2_hdr, "\n")
  if (!is.null(s2_note)) cat(s2_note)
  if (is.null(x$estimate)) {
    cat("  Note: contrasts require at least 2 interventions.\n")
  } else {
    # Mediation decomposition carries both RD and RR columns with their own
    # CI sets; the gformula contrast table has one estimate column whose
    # scale is given per row by Risk_type.
    print(.fmt_for_print(x$estimate, digits,
                         ci_labels = if (is_mediation) "rd_rr" else "plain"))
  }

  # -- Identifiability caveat (natural effects with intermediate confounding) --
  # Re-surface the short version of the runtime warning emitted by
  # check_natural_identifiability(), which is easy to miss in the warning
  # stream after a long run. Both state only what the formula scan found and
  # what the cited papers establish; neither asserts the causal structure.
  ic <- x$intermediate_confounders
  if (is_mediation && identical(args$mediation_type, "N") &&
      !is.null(ic) && length(ic) > 0) {
    cat(sprintf(
      paste0("\n  ! Identifiability: covariate model(s) for {%s} include the ",
             "exposure '%s', i.e. are modelled as exposure-affected.\n",
             "    If such a covariate also confounds the mediator-outcome ",
             "relationship, the natural direct/indirect\n",
             "    effects are not point-identified (Avin, Shpitser & Pearl 2005; ",
             "VanderWeele & Tchetgen Tchetgen 2017,\n",
             "    who propose the interventional analogues, mediation_type = ",
             "\"I\", for that setting).\n"),
      paste(ic, collapse = ", "), args$exposure))
  }

  # -- CI footnotes -----------------------------------------------------------
  if (is_tmle) {
    cat("\n  95% CIs (norm): Wald intervals from the efficient influence curve;\n")
    cat("  multiply robust under the conditions of Zheng & van der Laan (2017, Thm 1).\n")
  }
  if (!is.null(args$R) && args$R > 1) {
    cat(sprintf("\n  95%% CIs: percentile (pct) and normal approximation (norm) from %d bootstrap replicates.\n",
                as.integer(args$R)))

    # Report bootstrap replicates that were dropped (non-finite) from the
    # CI computations, per effect. Silence means all replicates contributed.
    be <- x$boot_estimates$effects
    if (is_mediation && !is.null(be)) {
      n_rep  <- length(unique(be$replicate))
      bad_rd <- tapply(!is.finite(be$RD), be$Effect, sum)
      bad_rd <- bad_rd[bad_rd > 0]
      # RR drops only matter for effects that carry an RR scale at all.
      has_rr <- tapply(is.finite(be$RR), be$Effect, any)
      bad_rr <- tapply(!is.finite(be$RR), be$Effect, sum)
      bad_rr <- bad_rr[names(bad_rr) %in% names(has_rr)[has_rr] & bad_rr > 0]
      if (length(bad_rd) > 0 || length(bad_rr) > 0) {
        cat("  Note: non-finite bootstrap draws excluded from the CIs:\n")
        for (ef in names(bad_rd))
          cat(sprintf("    %s: %d/%d replicates (RD/value scale)\n",
                      ef, as.integer(bad_rd[[ef]]), n_rep))
        for (ef in names(bad_rr))
          cat(sprintf("    %s: %d/%d replicates (RR scale)\n",
                      ef, as.integer(bad_rr[[ef]]), n_rep))
      }
    }

    # Mediation Proportion is unstable when the total effect is near zero.
    if (is_mediation && !is.null(x$estimate) &&
        all(c("perct_lcl", "perct_ucl") %in% names(x$estimate))) {
      te_row <- x$estimate[x$estimate$Effect == "Total effect", ]
      if (nrow(te_row) == 1L &&
          is.finite(te_row$perct_lcl) && is.finite(te_row$perct_ucl) &&
          te_row$perct_lcl <= 0 && te_row$perct_ucl >= 0) {
        cat("  Warning: the Total effect 95% CI includes 0; the Mediation Proportion\n")
        cat("           is numerically unstable and should be interpreted with caution.\n")
      }
    }
  }

  # -- Model details (optional) ----------------------------------------------
  if (models) {
    cat("\n--- Fitted models ---\n")
    .print_fitted_models(x$fitted_models, args$return_fitted, digits,
                         full = TRUE)
  }

  invisible(x)
}


# Internal helper ---------------------------------------------------------------
# Render the fitted-model block shared by print.gformula(models = TRUE) and
# summary.gformula(). `full = TRUE` prints each model's recodes and, when the
# object was built with return_fitted = TRUE, its complete summary();
# `full = FALSE` prints the coefficient table only, rounded to `digits`.
.print_fitted_models <- function(fitted_models, return_fitted, digits,
                                 full = TRUE) {
  for (nm in names(fitted_models)) {
    fitmod   <- fitted_models[[nm]]
    attr_lst <- attributes(fitmod)
    cat(sprintf("\nVariable : %s  [type: %s | role: %s]\n",
                nm, attr_lst$var_type, attr_lst$mod_type))
    cat(paste(rep("-", 40L), collapse = ""), "\n")

    if (full) {
      if (!is.null(attr_lst$recodes)) {
        cat("Recodes: ")
        print(attr_lst$recodes)
      }
      if (isTRUE(return_fitted)) {
        print(summary(fitmod))
      } else {
        if (!is.null(fitmod$call)) {
          cat("Call: ")
          print(fitmod$call)
        }
        # NULL whenever no numeric coefficient table could be extracted --
        # e.g. a custom_fit model of a class summary() does not describe that
        # way. Say so instead of printing "NULL".
        if (!is.null(fitmod$coeff)) {
          cat("Coefficients:\n")
          print(fitmod$coeff)
        } else {
          cat("  (no coefficient table available for this model class;",
              "rerun with return_fitted = TRUE for the full object)\n")
        }
      }
      next
    }

    coef_tbl <- if (isTRUE(return_fitted)) {
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
  invisible(NULL)
}


#' Summary
#'
#' Summary method for objects returned by \code{\link{gformula}} or
#' \code{\link{mediation}}. Shows the full estimation results followed by
#' fitted model coefficient tables.
#'
#' @param object Object of class \code{"gformula"}.
#' @param digits Integer. Number of decimal places used when rounding
#'   numeric output. Default \code{max(3, getOption("digits") - 3)}.
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
  .print_fitted_models(object$fitted_models, object$all.args$return_fitted,
                       digits, full = FALSE)

  invisible(object)
}
