#' Calculate Risk Difference, Risk Ratio, or Mediation Effects based on g-formula
#'
#' @description
#' These internal functions compute the Risk Difference (RD) and Risk Ratio (RR)
#' for g-formula estimates of the total effect, natural direct effect (NDE), 
#' and natural indirect effect (NIE). The functions also calculate confidence intervals 
#' (for RD and RR) and mediation effects based on g-formula methodology.
#'
#' @details
#' These functions are used internally within the g-formula framework for effect estimation.
#' They are not intended for direct user input. The functions perform:
#' - \code{risk_estimate.point}: Computes point estimates for risk difference and risk ratio.
#' - \code{risk_estimate.boot}: Computes per-bootstrap-replicate estimates for the effects (RD/RR).
#' - \code{risk_estimate.mediation}: Computes total, direct, and indirect effects by subtraction.
#'
#' These functions operate on the output of the \code{.gformula} function, estimating 
#' effects (both point estimates and confidence intervals) from simulated data.
#'
#' @param data_list A list of simulated datasets returned by \code{.gformula}.
#' @param ref_int Reference intervention to compare against.
#' @param intervention A vector specifying the interventions to compare.
#' @param return_data Logical. If \code{TRUE}, returns the simulated data along with effects.
#'
#' @return
#' A \code{data.table} with calculated effects (Risk Difference, Risk Ratio, or Mediation effects)
#' and, if applicable, confidence intervals. For mediation effects, the result will show
#' total, direct, and indirect effects.
#'
#' @seealso \code{\link{gformula}}, \code{\link{mediation}}
#'
#' @importFrom data.table data.table
#'
#' @keywords internal
#' @export


# Point estimates: risk difference and risk ratio from the main g-formula run
risk_estimate.point <- function(data_list, ref_int, intervention, return_data) {
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
risk_estimate.boot <- function(data_list, ref_int, intervention, return_data) {
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

# Calculate mediation effects from a named list of arm Phi values.
# Supports N >= 1 mediators (Yamamuro et al. 2021 multi-mediator extension).
#
# data_list must be a named list keyed by arm name (see build_mediation_arms):
#   Phi00, Phi10, Phi11                   (always required)
#   Phi1_k for k = 1..N-1                 (required when N >= 2)
#   nat0, nat1                            (interventional path only; natural
#                                          never-/always-treat arms for the
#                                          plug-in total effect)
# Each element is a scalar Phi value (mean Pred_Y) or a numeric vector that
# will be averaged. The legacy `return_data = TRUE` path (data tables) is
# handled by the caller in mediation.r before calling this function.
#
# Decomposition (sequential, Phi1_0 = Phi10, Phi1_N = Phi11):
#   IDE      = Phi10  - Phi00
#   IIE(M_k) = Phi1_k - Phi1_{k-1}
#   OE       = Phi11  - Phi00 = IDE + sum_k IIE(M_k)   (overall via the arms)
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
# natural reference arms, TE = Phi11 - Phi00, the decomposition sums exactly to
# TE and no residual row is added (Zheng & van der Laan 2017).
risk_estimate.mediation <- function(data_list, med_vars = NULL, return_data = FALSE) {

  # Robust scalar extraction (handles vector Pred_Y when return_data legacy
  # path is used, plus the attached "phi_estimate" attribute set by .gformula
  # for n_vw-averaged arms).
  to_phi <- function(name) {
    v <- data_list[[name]]
    if (is.null(v)) return(NA_real_)
    attr_phi <- attr(v, "phi_estimate", exact = TRUE)
    if (!is.null(attr_phi)) return(as.numeric(attr_phi))
    if (length(v) == 1L) return(as.numeric(v))
    return(sum(v) / length(v))
  }

  phi_00 <- to_phi("Phi00")
  phi_11 <- to_phi("Phi11")
  phi_10 <- to_phi("Phi10")

  # Interventional path uses separate natural-course arms for the total effect.
  interventional <- !is.null(data_list[["nat0"]]) && !is.null(data_list[["nat1"]])
  nat0 <- if (interventional) to_phi("nat0") else phi_00
  nat1 <- if (interventional) to_phi("nat1") else phi_11

  # Number of mediators inferred from med_vars (preferred) or from the number
  # of Phi1_k arms present in data_list (fallback).
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

