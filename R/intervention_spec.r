# Internal intervention specification for mediation analysis.
#
# Each `causalMed_intervention` describes one Monte Carlo trajectory to run:
#   - treatment:           scalar exposure value held constant across time
#   - mediator_overrides:  named list keyed by mediator response variable
#                          name. Each value is the exposure level (0 or 1)
#                          that determines how the counterfactual mediator
#                          is generated:
#                            * mediation_type = "I": the reference-pool
#                              source. The mediator is sampled jointly from
#                              the M(1:T) trajectories collected during the
#                              Phi00 (a=0) or Phi11 (a=1) reference intervention.
#                            * mediation_type = "N": the swap-target
#                              exposure for re-evaluating the mediator
#                              model on the intervention's own covariate history.
#                          An empty list means "no cross-world override" —
#                          mediators are simulated from their fitted models
#                          under the intervention's own treatment trajectory
#                          (reference interventions Phi00 and Phi11).
#
# This object is constructed internally by `mediation()` and consumed by
# `.run_interventions()` / `simulate_intervention()` / `simulate_data()`. It is not exposed in
# the public API.
intervention_spec <- function(treatment, mediator_overrides = list()) {
  if (!is.numeric(treatment) || length(treatment) != 1L) {
    stop("intervention_spec: 'treatment' must be a numeric scalar.", domain = "causalMed")
  }
  if (!is.list(mediator_overrides)) {
    stop("intervention_spec: 'mediator_overrides' must be a (possibly empty) named list.",
         domain = "causalMed")
  }
  if (length(mediator_overrides) > 0L &&
      (is.null(names(mediator_overrides)) ||
       any(nchar(names(mediator_overrides)) == 0L))) {
    stop("intervention_spec: 'mediator_overrides' must be a named list keyed by mediator response variable.",
         domain = "causalMed")
  }

  structure(
    list(treatment          = as.numeric(treatment),
         mediator_overrides = mediator_overrides),
    class = "causalMed_intervention"
  )
}

is_intervention_spec <- function(x) inherits(x, "causalMed_intervention")


# Internal helper used by mediation() to build the list of intervention_spec objects
# required for the requested mediation decomposition.
#
# ---- Interventional effects (mediation_type = "I") --------------------------
# Lin et al. (2017); VanderWeele & Tchetgen Tchetgen (2017); Yamamuro et al.
# (2021).  EVERY intervention in the decomposition draws each mediator from an
# INDEPENDENTLY-PERMUTED marginal pool -- including the reference interventions.  Writing
# G_{a*} for such a stochastic mediator draw under treatment a*:
#   Phi00  : treatment 0, all mediators from a=0 pool   = E[Y_{0,G0}]   (Ga0...)
#   Phi10  : treatment 1, all mediators from a=0 pool   = E[Y_{1,G0}]   (Ga1m0..)
#   Phi1_k : treatment 1, first k mediators from a=1 pool, rest a=0 pool (k<N)
#   Phi11  : treatment 1, all mediators from a=1 pool   = E[Y_{1,G1}]   (Ga1...mN)
# Sequential decomposition (Phi1_0 = Phi10, Phi1_N = Phi11):
#   IDE      = Phi10  - Phi00
#   IIE(M_k) = Phi1_k - Phi1_{k-1}
#   OE       = Phi11  - Phi00 = IDE + sum_k IIE(M_k)   (interventional overall)
# The natural plug-in TOTAL effect is NOT Phi11 - Phi00 here; it is obtained
# from the separate natural-course interventions
#   nat0 : treatment 0, mediators NATURAL  = E[Y_0]
#   nat1 : treatment 1, mediators NATURAL  = E[Y_1]
# so that TE = nat1 - nat0 and the mediated-interaction residual TE - OE is
# generally non-zero (matching the SAS mGFORMULA macro / Yamamuro Table 3).
#
# ---- Natural effects (mediation_type = "N") ---------------------------------
# Zheng & van der Laan (2017), single mediator only.  Here the reference interventions
# legitimately use the NATURAL mediator distribution and the decomposition sums
# exactly to the total effect, so the interventions are left as natural-reference interventions
# (override values are swap-target exposures for re-evaluating the mediator
# model, not pool sources).
build_mediation_interventions <- function(med_vars, mediation_type = "I") {
  N <- length(med_vars)
  if (N == 0L) stop("build_mediation_interventions: no mediator variables.")

  all0 <- setNames(as.list(rep(0, N)), med_vars)
  all1 <- setNames(as.list(rep(1, N)), med_vars)

  if (identical(mediation_type, "N")) {
    interventions <- list()
    interventions$Phi00 <- intervention_spec(0, list())   # natural never-treat
    interventions$Phi11 <- intervention_spec(1, list())   # natural always-treat
    interventions$Phi10 <- intervention_spec(1, all0)     # cross-world (swap mediator to a=0)
    interventions$Phi01 <- intervention_spec(0, all1)
    return(interventions)
  }

  # Interventional: permuted-pool references + separate natural-course interventions.
  interventions <- list()
  interventions$nat0  <- intervention_spec(0, list())     # E[Y_0]  natural never-treat  (-> TE)
  interventions$nat1  <- intervention_spec(1, list())     # E[Y_1]  natural always-treat (-> TE)
  interventions$Phi00 <- intervention_spec(0, all0)       # E[Y_{0,G0}]  interventional reference
  interventions$Phi10 <- intervention_spec(1, all0)       # E[Y_{1,G0}]
  if (N >= 2L) {
    for (k in seq_len(N - 1L)) {
      interventions[[sprintf("Phi1_%d", k)]] <- intervention_spec(
        1, setNames(as.list(c(rep(1, k), rep(0, N - k))), med_vars)
      )
    }
  }
  interventions$Phi11 <- intervention_spec(1, all1)       # E[Y_{1,G1}]  interventional reference
  interventions
}
