#' Point-estimate contrasts against a reference intervention
#'
#' @description
#' Internal use only. Turns the per-intervention mean outcomes produced by
#' \code{.run_interventions} into risk differences and risk ratios against the
#' reference intervention. Confidence intervals are \emph{not} computed here;
#' \code{\link{gformula}} adds them from the bootstrap replicates.
#'
#' @details
#' Companion internal functions in the same file, not documented separately:
#' \code{risk_estimate_boot} does the same job for one stacked bootstrap
#' replicate, and \code{risk_estimate_mediation} builds the mediation
#' decomposition (total, direct, per-mediator indirect) by subtraction.
#'
#' @param data_list Named list of intervention results from
#'   \code{.run_interventions}: one mean outcome per intervention, or one
#'   simulated \code{data.table} per intervention when \code{return_data} is
#'   \code{TRUE}.
#' @param ref_int Character scalar naming the reference element of
#'   \code{data_list}.
#' @param intervention The named intervention list; its names determine which
#'   contrasts are formed.
#' @param return_data Logical. \code{TRUE} when \code{data_list} holds
#'   simulated datasets rather than scalar means.
#'
#' @return
#' A \code{data.table} with columns \code{Intervention}, \code{Risk_type}
#' (\code{"Difference"} or \code{"Ratio"}) and \code{Estimate}: two rows per
#' non-reference intervention.
#'
#' @seealso \code{\link{gformula}}, \code{\link{mediation}}
#'
#' @importFrom data.table data.table
#'
#' @keywords internal
risk_estimate_point <- function(data_list, ref_int, intervention, return_data) {
  ref_dat <- data_list[[ref_int]]
  vs_nam <- setdiff(names(intervention), ref_int)
  out <- sapply(vs_nam, function(x) {
    tmp_dt <- data_list[[x]]
    if (return_data) {
      ref_mean <- sum(ref_dat[["Pred_Y"]]) / length(ref_dat[["Pred_Y"]])
      vs_mean  <- sum(tmp_dt[["Pred_Y"]])  / length(tmp_dt[["Pred_Y"]])
    } else {
      ref_mean <- ref_dat
      vs_mean  <- tmp_dt
    }
    data.table(
      Intervention = c(paste(x, ref_int, sep = " - "),
                       paste(x, ref_int, sep = " / ")),
      Risk_type = c("Difference", "Ratio"),
      Estimate = c(vs_mean - ref_mean, vs_mean / ref_mean)
    )
  }, simplify = FALSE)
  data.table::rbindlist(out, use.names = TRUE)
}

# Bootstrap replicate estimates: risk difference and risk ratio from a stacked bootstrap data.table
risk_estimate_boot <- function(data_list, ref_int, intervention, return_data) {
  ref_dat <- data_list[data_list$Intervention==ref_int,]
  vs_nam <- setdiff(names(intervention), ref_int)
  out <- sapply(vs_nam, function(x) {
    tmp_dt <- data_list[data_list$Intervention==x,]
    if (return_data) {
      ref_mean <- sum(ref_dat[["Est"]]) / length(ref_dat[["Est"]])
      vs_mean  <- sum(tmp_dt[["Est"]])  / length(tmp_dt[["Est"]])
    } else {
      ref_mean <- ref_dat$Est
      vs_mean  <- tmp_dt$Est
    }
    data.table(
      Intervention = c(paste(x, ref_int, sep = " - "),
                       paste(x, ref_int, sep = " / ")),
      Risk_type = c("Difference", "Ratio"),
      Estimate = c(vs_mean - ref_mean, vs_mean / ref_mean)
    )
  }, simplify = FALSE)
  data.table::rbindlist(out, use.names = TRUE)
}

# Scalar Phi from one intervention result: the "phi_estimate" attribute when
# .run_interventions attached one (n_vw-averaged interventions), the mean of
# Pred_Y for a returned simulated dataset (return_data = TRUE), the mean of a
# Pred_Y vector (legacy path), or the scalar itself.
phi_scalar <- function(x) {
  if (is.null(x)) return(NA_real_)
  attr_phi <- attr(x, "phi_estimate", exact = TRUE)
  if (!is.null(attr_phi)) return(as.numeric(attr_phi))
  if (is.data.frame(x)) x <- x[["Pred_Y"]]
  if (length(x) == 1L) return(as.numeric(x))
  sum(x) / length(x)
}

# Proportion mediated (additive + multiplicative, in percent) from a named
# list/vector of intervention Phi values and its risk_estimate_mediation()
# table. Additive PM follows Yamamuro et al. (2021): (TE - IDE) / TE -- the
# residual (TE - OE) is folded into the mediated portion, so for
# mediation_type = "N" (zero residual) it reduces to the usual sum(IIE)/TE.
# Multiplicative PM (Lin et al. 2017 Table 2, generalised to N mediators via
# Yamamuro et al. 2021): RR_IDE * (prod RR_IIE_k - 1) / (RR_OE - 1) on the
# interventional overall scale RR_OE = Phi11/Phi00 (under "N" that ratio is
# the total-effect RR and the formula is reported as an analogue). NA where
# the denominator is degenerate or any component RR is non-finite.
pm_from_phi <- function(phi, risk) {
  total  <- risk$RD[risk$Effect == "Total effect"]
  direct <- risk$RD[risk$Effect == "Direct effect"]
  rr_ide <- risk$RR[risk$Effect == "Direct effect"]
  rr_iie <- risk$RR[grepl("^Indirect effect", risk$Effect)]
  rr_oe  <- phi_scalar(phi[["Phi11"]]) / phi_scalar(phi[["Phi00"]])

  add <- if (isTRUE(abs(total) < 1e-10)) NA_real_ else
    (total - direct) / total * 100
  mult <- if (isTRUE(abs(rr_oe - 1) < 1e-10) ||
              !all(is.finite(c(rr_oe, rr_ide, rr_iie)))) NA_real_ else
    rr_ide * (prod(rr_iie) - 1) / (rr_oe - 1) * 100
  c(add = add, mult = mult)
}

# Calculate mediation effects from a named list of intervention Phi values.
# Supports N >= 1 mediators (Yamamuro et al. 2021 multi-mediator extension).
#
# data_list must be a named list keyed by intervention name (see build_mediation_interventions):
#   Phi00, Phi10, Phi11                   (always required)
#   Phi1_k for k = 1..N-1                 (required when N >= 2)
#   nat0, nat1                            (interventional path only; natural
#                                          never-/always-treat interventions for the
#                                          plug-in total effect)
# Each element is a scalar Phi value (mean Pred_Y) or a numeric vector that
# will be averaged. The legacy `return_data = TRUE` path (data tables) is
# handled by the caller in mediation.r before calling this function.
#
# Decomposition (sequential, Phi1_0 = Phi10, Phi1_N = Phi11):
#   IDE      = Phi10  - Phi00
#   IIE(M_k) = Phi1_k - Phi1_{k-1}
#   OE       = Phi11  - Phi00 = IDE + sum_k IIE(M_k)   (overall via the interventions)
#
# INTERVENTIONAL path (nat0/nat1 present, mediation_type = "I"): Phi00 and Phi11
# are the PERMUTED-pool interventional references (E[Y_{0,G0}], E[Y_{1,G1}]), so
# OE is the interventional overall effect. The reported TOTAL effect is the
# natural plug-in contrast nat1 - nat0, and an extra row
#   "TE - (Direct + Indirect)" = TE - OE
# reports the (generally non-zero) mediated-interaction residual
# (Lin et al. 2017; VanderWeele & Tchetgen Tchetgen 2017; Yamamuro et al. 2021).
#
# NATURAL path (nat0/nat1 absent, mediation_type = "N"): Phi00/Phi11 are the
# natural reference interventions, TE = Phi11 - Phi00, the decomposition sums exactly to
# TE and no residual row is added (Zheng & van der Laan 2017).
risk_estimate_mediation <- function(data_list, med_vars = NULL, return_data = FALSE) {

  # Robust scalar extraction (see phi_scalar() above).
  to_phi <- function(name) phi_scalar(data_list[[name]])

  phi_00 <- to_phi("Phi00")
  phi_11 <- to_phi("Phi11")
  phi_10 <- to_phi("Phi10")

  # Interventional path uses separate natural-course interventions for the total effect.
  interventional <- !is.null(data_list[["nat0"]]) && !is.null(data_list[["nat1"]])
  nat0 <- if (interventional) to_phi("nat0") else phi_00
  nat1 <- if (interventional) to_phi("nat1") else phi_11

  # Number of mediators inferred from med_vars (preferred) or from the number
  # of Phi1_k interventions present in data_list (fallback).
  if (is.null(med_vars)) {
    n_med <- 1L + sum(grepl("^Phi1_\\d+$", names(data_list)))
  } else {
    n_med <- length(med_vars)
  }

  # Build the per-mediator sequence Phi1_0..Phi1_N
  phi_seq         <- numeric(n_med + 1L)
  phi_seq[1L]     <- phi_10
  phi_seq[n_med + 1L] <- phi_11
  if (n_med >= 2L) {
    for (k in seq_len(n_med - 1L)) {
      phi_seq[k + 1L] <- to_phi(sprintf("Phi1_%d", k))
    }
  }

  # Per-mediator IIE
  iie_rd <- diff(phi_seq)
  iie_rr <- phi_seq[-1L] / phi_seq[-length(phi_seq)]

  # Build the effects table. For N=1 keep the existing "Indirect effect" /
  # "Direct effect" / "Total effect" labels (back-compat with existing
  # tests and downstream code). For N>=2 label each indirect effect by its
  # mediator's response variable name.
  if (n_med == 1L) {
    iie_labels <- "Indirect effect"
  } else {
    nm <- if (!is.null(med_vars)) med_vars else paste0("M", seq_len(n_med))
    iie_labels <- sprintf("Indirect effect (%s)", nm)
  }

  te_rd <- nat1 - nat0
  te_rr <- nat1 / nat0

  out <- data.table(
    Effect = c(iie_labels, "Direct effect", "Total effect"),
    RD     = c(iie_rd,            phi_10 - phi_00, te_rd),
    RR     = c(iie_rr,            phi_10 / phi_00, te_rr)
  )

  # Interventional residual: TE (natural) minus the interventional overall
  # effect OE = Phi11 - Phi00 = IDE + sum_k IIE(M_k).
  if (interventional) {
    oe_rd <- phi_11 - phi_00
    out <- rbind(out, data.table(
      Effect = "TE - (Direct + Indirect)",
      RD     = te_rd - oe_rd,
      RR     = NA_real_
    ))
  }

  out
}

