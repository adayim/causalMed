# Internal intervention specification for mediation analysis.
#
# Each `causalMed_intervention` describes one Monte Carlo trajectory to run:
#   - treatment:           the exposure regime a(1:T) -- a numeric vector with
#                          one value per time point (length 1 = the same value
#                          at every step). simulate_intervention() slices it to
#                          the current step so simulate_data() only ever sees a
#                          scalar.
#   - mediator_overrides:  named list keyed by mediator response variable
#                          name. Each value is a REGIME (same form as
#                          `treatment`) that determines how the counterfactual
#                          mediator is generated:
#                            * mediation_type = "I": the reference-pool
#                              source. The mediator is sampled jointly from
#                              the M(1:T) trajectories collected under that
#                              regime (pools are keyed by regime_key()).
#                            * mediation_type = "N": the swap-target regime
#                              for re-evaluating the mediator model on the
#                              intervention's own covariate history (current
#                              exposure, and its previous-step value for the
#                              first-order exposure lags).
#                          An empty list means "no mediator override" --
#                          mediators are simulated from their fitted models
#                          under the intervention's own treatment regime
#                          (nat0/nat1, and the pool-collection passes).
#
# The regimes are STATIC vectors, as in Lin et al. (2017, Stat Med, Section
# 2.2) and Zheng & van der Laan (2017, Section 2.2); always/never exposed
# (1 and 0, recycled) is their worked example and mediation()'s default.
#
# Constructed internally by `mediation()` and consumed by
# `.run_interventions()` / `simulate_intervention()` / `simulate_data()`. Not
# exposed in the public API.
intervention_spec <- function(treatment, mediator_overrides = list()) {
  if (!is.numeric(treatment) || length(treatment) < 1L || anyNA(treatment) ||
      !is.null(dim(treatment))) {
    stop("intervention_spec: 'treatment' must be a numeric vector of length >= 1 (an exposure regime).",
         domain = "causalMed")
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
  bad_val <- vapply(mediator_overrides, function(v) {
    !is.numeric(v) || length(v) < 1L || anyNA(v)
  }, logical(1))
  if (any(bad_val)) {
    stop("intervention_spec: every 'mediator_overrides' value must be a numeric regime vector of length >= 1.",
         domain = "causalMed")
  }
  lens <- c(length(treatment),
            vapply(mediator_overrides, length, integer(1)))
  if (!all(lens %in% c(1L, max(lens)))) {
    stop(sprintf(paste0("intervention_spec: regime lengths must all be 1 or the same ",
                        "length; got {%s} for treatment and overrides {%s}."),
                 length(treatment),
                 paste(names(mediator_overrides), vapply(mediator_overrides, length, integer(1)),
                       sep = "=", collapse = ", ")),
         domain = "causalMed")
  }

  structure(
    list(treatment          = as.numeric(treatment),
         mediator_overrides = lapply(mediator_overrides, as.numeric)),
    class = "causalMed_intervention"
  )
}

is_intervention_spec <- function(x) inherits(x, "causalMed_intervention")

# Canonical string key of a regime, used for the mediator-pool and cached-arm
# lookups in .run_interventions(). Exact: c(1, 1) and 1 get different keys.
# Callers MUST recycle regimes to full length before building interventions,
# so that a regime and the pool collected under it produce the same key;
# `mediation()` does this. A length-1 and a recycled length-T regime are
# DIFFERENT keys.
regime_key <- function(x) paste(as.numeric(x), collapse = ",")

# The scalar intervention_spec for time step `k`: `treatment[k]` and the k-th
# element of every override (a length-1 regime applies at every step). Built
# with structure() rather than intervention_spec() so the per-step slice in
# simulate_intervention()'s time loop does not re-validate. The slice is not
# a pool-key source -- regime_key() must be called on the full spec, never
# on a slice.
slice_intervention_spec <- function(spec, k) {
  pick <- function(v) if (length(v) == 1L) v else v[[k]]
  structure(
    list(treatment          = pick(spec$treatment),
         mediator_overrides = lapply(spec$mediator_overrides, pick)),
    class = "causalMed_intervention"
  )
}

# Internal helper used by mediation() to build the list of intervention_spec objects
# required for the requested mediation decomposition.
#
# ---- Interventional effects (mediation_type = "I") --------------------------
# Lin et al. (2017); VanderWeele & Tchetgen Tchetgen (2017); Yamamuro et al.
# (2021).  EVERY intervention in the decomposition draws each mediator from an
# INDEPENDENTLY-PERMUTED marginal pool -- including the reference interventions.  Writing
# G_{a*} for such a stochastic mediator draw under treatment a*:
#   (a = `regime`, a* = `ref_regime`; the defaults 1 / 0 are always/never)
#   Phi00  : treatment a*, all mediators from a* pool   = E[Y_{a*,G_a*}]
#   Phi10  : treatment a,  all mediators from a* pool   = E[Y_{a,G_a*}]
#   Phi1_k : treatment a,  first k mediators from a pool, rest a* pool (k<N)
#   Phi11  : treatment a,  all mediators from a pool    = E[Y_{a,G_a}]
# Sequential decomposition (Phi1_0 = Phi10, Phi1_N = Phi11):
#   IDE      = Phi10  - Phi00
#   IIE(M_k) = Phi1_k - Phi1_{k-1}
#   OE       = Phi11  - Phi00 = IDE + sum_k IIE(M_k)   (interventional overall)
# The natural plug-in TOTAL effect is NOT Phi11 - Phi00 here; it comes from two
# fixed-regime, natural-mediator interventions (exposure fixed to the regime,
# mediators from their own fitted models -- NOT gformula()'s natural course,
# which draws the exposure too):
#   nat0 : treatment a*, mediators NATURAL  = E[Y_{a*}]
#   nat1 : treatment a,  mediators NATURAL  = E[Y_{a}]
# so TE = nat1 - nat0 and the residual TE - OE is generally non-zero (matching
# the SAS mGFORMULA macro / Yamamuro Table 3).
#
# ---- Natural effects (mediation_type = "N") ---------------------------------
# Zheng & van der Laan (2017), single mediator only.  The references legitimately
# use the NATURAL mediator distribution and the decomposition sums exactly to TE,
# so they are left as natural references (override values are swap-target
# regimes for re-evaluating the mediator model, not pool sources).
build_mediation_interventions <- function(med_vars, mediation_type = "I",
                                          regime = 1, ref_regime = 0) {
  N <- length(med_vars)
  if (N == 0L) {
    stop("build_mediation_interventions: no mediator variables.", domain = "causalMed")
  }

  all0 <- setNames(rep(list(ref_regime), N), med_vars)
  all1 <- setNames(rep(list(regime),     N), med_vars)

  if (identical(mediation_type, "N")) {
    interventions <- list()
    interventions$Phi00 <- intervention_spec(ref_regime, list())   # natural, reference regime
    interventions$Phi11 <- intervention_spec(regime,     list())   # natural, exposure regime
    interventions$Phi10 <- intervention_spec(regime,     all0)     # cross-regime
    interventions$Phi01 <- intervention_spec(ref_regime, all1)
    return(interventions)
  }

  # Interventional: permuted-pool references + separate natural-course interventions.
  interventions <- list()
  interventions$nat0  <- intervention_spec(ref_regime, list())     # E[Y_{a*}]  (-> TE)
  interventions$nat1  <- intervention_spec(regime,     list())     # E[Y_{a}]   (-> TE)
  interventions$Phi00 <- intervention_spec(ref_regime, all0)       # E[Y_{a*,G_a*}]  interventional reference
  interventions$Phi10 <- intervention_spec(regime,     all0)       # E[Y_{a,G_a*}]
  if (N >= 2L) {
    for (k in seq_len(N - 1L)) {
      interventions[[sprintf("Phi1_%d", k)]] <- intervention_spec(
        regime,
        setNames(c(rep(list(regime), k), rep(list(ref_regime), N - k)), med_vars)
      )
    }
  }
  interventions$Phi11 <- intervention_spec(regime, all1)           # E[Y_{a,G_a}]  interventional reference
  interventions
}
