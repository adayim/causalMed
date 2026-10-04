#' Mediation analysis with the parametric g-formula
#'
#' @description
#' Decomposes the effect of a time-varying binary exposure into a direct effect
#' and an indirect effect through one or more time-varying mediators, using the
#' parametric mediational g-formula. Two estimands are available:
#' interventional effects (\code{mediation_type = "I"}, the default; Lin et al.
#' 2017, VanderWeele & Tchetgen Tchetgen 2017) and natural effects
#' (\code{"N"}; Zheng & van der Laan 2017).
#'
#' @details
#' \strong{Data.} Long format: one row per subject per time point. For a
#' survival outcome, remove every row after the event and after loss to
#' follow-up. The exposure must be coded 0/1; the methods implemented here are
#' defined for binary exposures only.
#'
#' \strong{Models.} \code{models} is a list of \code{\link{spec_model}}
#' objects in the order the variables are generated within a time point, for
#' example exposure, confounders, mediator, outcome
#' (\strong{A(t) -> L(t) -> M(t) -> Y(t)}), or \strong{A -> M -> L -> Y} when
#' the confounders respond to the mediator. At least one \code{"mediator"}
#' model is required. Several mediators, listed in temporal order, are
#' supported under \code{"I"} only (Yamamuro et al. 2021). A warning is given
#' if the exposure is listed after the mediator or the mediator after the
#' outcome.
#'
#' \strong{Exposure regimes.} The effects contrast two \emph{static} regimes,
#' \eqn{a} (\code{exposure_regime}) and \eqn{a^*} (\code{reference_regime}).
#' The defaults, always versus never exposed, are the example Lin et al.
#' (2017, \emph{Stat Med}, Section 2.2) and Zheng & van der Laan (2017,
#' Section 2.2) give when defining the effects for arbitrary regimes. Element
#' \eqn{k} of a regime is the exposure at the \eqn{k}-th distinct value of
#' \code{time_var}: for times \{1, 3, 4, 6\} a regime has four elements. A
#' regime in which the exposure follows its fitted model outside a window is a
#' dynamic regime, a different estimand, and is not available. Which pair
#' answers a given question is the analyst's judgement; only \eqn{a = a^*} is
#' rejected.
#'
#' \strong{How the mediator is set under \code{"I"}.} The model is first
#' simulated once under each regime with the mediator following its fitted
#' model, and every simulated subject's mediator trajectory \eqn{M(1{:}T)} is
#' stored in a pool. Each decomposition intervention (\code{Phi00},
#' \code{Phi10}, \code{Phi11}, and \code{Phi1_k} with several mediators) then
#' gives each subject the whole trajectory of a randomly permuted pool member,
#' each mediator permuted independently, and averages \code{n_vw} such
#' permutations. The permutation runs over the whole simulated cohort, so the
#' draw is marginal over the baseline covariates as well as the confounder
#' history: the whole-population draw \eqn{G_a} of VanderWeele & Tchetgen
#' Tchetgen (2017), as in the algorithm of Lin et al. (2017,
#' \emph{Epidemiology}). Yamamuro et al. (2021, Eq. 2) and Lin et al. (2017,
#' \emph{Stat Med}, Eq. 4) are written for the draw conditional on baseline
#' covariates, but the algorithms of both papers, and the SAS macros
#' distributed with them, permute the whole cohort. The two forms agree when
#' the counterfactual mediator distribution does not depend on the baseline
#' covariates, or when the outcome mean is additively separable in them and the
#' mediator; otherwise they are different functionals.
#'
#' \strong{How the mediator is set under \code{"N"}.} The mediator model is
#' re-evaluated on each subject's own simulated history with the exposure
#' history set to the other regime (Zheng & van der Laan 2017, Eq. 5). Only the
#' exposure and its first-order lags (\code{in_recode} entries that copy the
#' exposure, e.g. \code{recodes(lag1_A = A)}) are switched, so a mediator model
#' whose formula or \code{subset} reads any other exposure-derived column (a
#' longer lag, a cumulative count, an \code{out_recode} copy, a column created
#' by a model's own \code{recode}) is rejected with an error.
#'
#' \strong{Total effect and residual.} The total effect is the natural plug-in
#' contrast \code{nat1 - nat0}, from runs in which the exposure follows each
#' regime and the mediator follows its fitted model. Under \code{"I"} the
#' direct and indirect effects sum to the interventional overall effect
#' \code{Phi11 - Phi00} rather than to the total effect, and the difference is
#' reported as the decomposition residual. Under \code{"N"} they sum to the
#' total effect and no residual is reported.
#'
#' \strong{Survival outcomes and censoring.} The algorithm is that of Lin et
#' al. (2017, \emph{Stat Med}, Section 4): models are fitted among those at
#' risk, no simulated subject is removed, and the cumulative risk
#' \eqn{1-\prod_t(1-\hat p_t)} is accumulated from the predicted hazards. A
#' \code{"censor"} model sets the censoring indicator to zero under every
#' intervention, so every risk is the risk under eliminated loss to follow-up;
#' with \code{estimator = "gcomp"} it does not change any reported risk. It
#' marks the analysis as a survival setting and is used by
#' \code{estimator = "tmle"}.
#'
#' \strong{Warnings.} Model-fitting warnings (e.g. non-convergence) are
#' collected and printed once, with a repeat count, when the function returns.
#'
#' @inheritParams gformula
#'
#' @param outcome Character. Name of the outcome variable; must be the
#'   response of the \code{"outcome"} or \code{"survival"} model.
#' @param seed Integer random seed (default \code{12345L}) for the Monte Carlo
#'   simulation and the bootstrap replicates; the global RNG state is restored
#'   on exit. \code{NULL} disables seeding, so repeated calls differ. Give
#'   concurrent analyses different seeds.
#' @param mediation_type \code{"I"} (default) for interventional effects or
#'   \code{"N"} for natural effects.
#'
#'   Zheng & van der Laan (2017, Lemma 1) identify the \code{"N"} effects,
#'   which draw the mediator from its conditional distribution under the
#'   other regime given each subject's own history, under sequential
#'   randomization and positivity. Reading them as \emph{individual-level}
#'   natural effects additionally requires a cross-world independence
#'   assumption that is not expected to hold when a confounder of the
#'   mediator-outcome relationship is affected by prior exposure (Avin,
#'   Shpitser & Pearl 2005; VanderWeele & Tchetgen Tchetgen 2017), the setting
#'   for which VanderWeele & Tchetgen Tchetgen (2017) propose the
#'   interventional effects of \code{"I"}. With \code{"N"}, a warning (repeated
#'   by \code{print()}) names any covariate whose model includes the exposure.
#'   The check reads formulas only: it does not establish that such a
#'   covariate confounds the mediator-outcome relationship, and its silence
#'   does not establish that none does.
#'
#'   Interventional indirect effects do not, absent stronger assumptions,
#'   satisfy the sharp null criterion of Miles (2023), which the natural
#'   indirect effect does: a non-zero interventional indirect effect does not
#'   by itself show that the mediator transmits the effect for any individual.
#' @param exposure_regime Exposure regime \eqn{a(1{:}T)} whose effect is
#'   decomposed: 0/1 values, either one value for every time point or one per
#'   distinct value of \code{time_var}, in sorted order. Default \code{1}
#'   (always exposed).
#' @param reference_regime Reference regime \eqn{a^*(1{:}T)}, in the same
#'   form. Default \code{0} (never exposed). Must differ from
#'   \code{exposure_regime}.
#' @param n_vw Integer. Number of mediator permutations averaged for each
#'   pool-drawing intervention under \code{"I"} (\code{Phi00}, \code{Phi10},
#'   \code{Phi11}, \code{Phi1_k}; not \code{nat0}/\code{nat1}), within every
#'   bootstrap replicate too. Default \code{2L}, as in the SAS
#'   \code{mGFORMULA} macro; \code{1L} is faster with more Monte Carlo noise.
#'   Ignored under \code{"N"}.
#' @param estimator \code{"gcomp"} (default): the parametric g-formula plug-in
#'   with bootstrap intervals. \code{"tmle"}: the targeted minimum loss-based
#'   estimator of Zheng & van der Laan (2017, Section 4.3), with Wald
#'   intervals from the efficient influence curve (\code{R} and
#'   \code{mc_sample} are ignored).
#'
#'   \code{"tmle"} requires \code{mediation_type = "N"}, the default regimes, a
#'   single mediator, a binary outcome or survival indicator, and an exposure
#'   model. It accepts only lag-style recodes: \code{in_recode} entries that
#'   copy one column (an exposure lag must copy the exposure itself),
#'   \code{init_recode} entries that are a constant or a column name, no
#'   \code{out_recode}, no \code{subset}, no \code{var_type = "custom"}, and no
#'   model \code{recode} that reads the exposure or its lags. Other inputs
#'   raise an error. Its targeted regressions are additive main-effects working
#'   models in the variables your formulas name (transformations and
#'   interactions such as \code{poly()} or \code{A:M} are not carried over), so
#'   the multiple robustness Zheng & van der Laan establish for correctly
#'   specified nuisance models is not claimed for this implementation. Where
#'   few subjects follow a regime, the affected fluctuation steps are skipped
#'   and weights truncated, with a warning.
#' @param tmle_weight_trunc Quantile in (0, 1] at which the TMLE
#'   clever-covariate weights are truncated. Default \code{0.995}; \code{1}
#'   disables truncation. Ignored unless \code{estimator = "tmle"}.
#'
#' @return
#' An object of class \code{"gformula"}, printed by
#' \code{\link{print.gformula}}, with components:
#' \itemize{
#'   \item \code{effect_size}: mean outcome (risk, for survival) under each
#'     intervention (\code{Intervention}, \code{Est}). Under \code{"I"}:
#'     \code{Phi00} \eqn{= E[Y_{0,G_0}]}, \code{Phi10} \eqn{= E[Y_{1,G_0}]},
#'     \code{Phi11} \eqn{= E[Y_{1,G_1}]} (mediator from a permuted pool), and
#'     \code{Phi1_k} (k = 1, ..., N-1) with N mediators; under \code{"N"},
#'     \code{Phi00} and \code{Phi11} use each subject's own mediator. Always
#'     \code{nat0} \eqn{= E[Y_0]} and \code{nat1} \eqn{= E[Y_1]}. With
#'     non-default regimes, 0 and 1 in these labels stand for \eqn{a^*} and
#'     \eqn{a}.
#'   \item \code{estimate}: the decomposition, on the risk-difference
#'     (\code{RD}) and risk-ratio (\code{RR}) scales. Rows:
#'     \code{"Indirect effect"} (\code{Phi11 - Phi10}; one row per mediator,
#'     \code{"Indirect effect (<mediator>)"}, with several),
#'     \code{"Direct effect"} (\code{Phi10 - Phi00}), \code{"Total effect"}
#'     (\code{nat1 - nat0}), \code{"TE - (Direct + Indirect)"} (\code{"I"}
#'     only), and \code{"Mediation Proportion"}: the summed indirect effects as
#'     a percentage of the interventional overall effect \code{Phi11 - Phi00},
#'     the quantity Lin et al. (2017, \emph{Stat Med}, Table 2) report. It is
#'     not a share of the total effect, from which it differs by the residual.
#'     No separate risk-ratio proportion is given:
#'     \eqn{RR_{IDE}(\prod_k RR_{IIE_k}-1)/(RR_{OE}-1)} simplifies to the same
#'     number. \code{RR} is \code{NA} for the residual and proportion rows.
#'   \item With \code{R > 1}, both tables gain \code{Sd} and percentile
#'     (\code{perct_lcl}, \code{perct_ucl}) and normal-approximation
#'     (\code{norm_lcl}, \code{norm_ucl}) limits, with an \code{_RR} suffix on
#'     the risk-ratio scale; \code{boot_estimates} holds the per-replicate
#'     estimates (\code{$interventions}, \code{$effects}).
#'   \item \code{sim_data}: with \code{return_data = TRUE}, the simulated data
#'     as an end-of-follow-up snapshot: one row per Monte Carlo subject per
#'     intervention, each variable at its last simulated time step, with the
#'     accumulated \code{Pred_Y}. With \code{n_vw > 1} only the last
#'     permutation is kept while \code{effect_size$Est} averages all of them,
#'     so \code{mean(Pred_Y)} reproduces \code{Est} for the pool-drawing
#'     interventions only when \code{n_vw = 1} (always for
#'     \code{nat0}/\code{nat1}).
#'   \item \code{fitted_models}: the fitted models, as full objects when
#'     \code{return_fitted = TRUE}, otherwise their calls and coefficients.
#'   \item \code{data_summary}: numbers of subjects, rows and time points in
#'     \code{data}, and \code{regime_support}, the number of observed subjects
#'     whose exposure matches each regime at every time they were observed
#'     (\code{n_following}) and among those observed at every time point
#'     (\code{n_complete}). A subject with a missing exposure counts for
#'     neither regime. The estimation does not use these counts.
#'   \item \code{observed}: a nonparametric benchmark printed beside the
#'     estimates: the observed mean outcome at the last time point, or the
#'     product-limit cumulative incidence for a survival outcome.
#'   \item \code{intermediate_confounders}: under \code{"N"}, covariates whose
#'     model includes the exposure (see \code{mediation_type}).
#'   \item \code{tmle_diag}: with \code{estimator = "tmle"}, the subject-level
#'     efficient influence curve (\code{eic}), its column means
#'     (\code{eic_mean}) and the number of subjects (\code{n}).
#'   \item \code{call}, \code{all.args}: the matched call and the evaluated
#'     arguments.
#' }
#'
#' @references
#' Avin, C., Shpitser, I., & Pearl, J. (2005). Identifiability of path-specific
#' effects. \emph{Proceedings of the 19th International Joint Conference on
#' Artificial Intelligence}, 357-363.
#'
#' Lin, S. H., Young, J. G., Logan, R., & VanderWeele, T. J. (2017).
#' Mediation analysis for a survival outcome with time-varying exposures, mediators, and confounders.
#' \emph{Statistics in Medicine}, 36(26), 4153–4166. \doi{10.1002/sim.7426}
#'
#' Lin, S. H., Young, J., Logan, R., Tchetgen Tchetgen, E. J., & VanderWeele,
#' T. J. (2017). Parametric mediational g-formula approach to mediation
#' analysis with time-varying exposures, mediators, and confounders.
#' \emph{Epidemiology}, 28(2), 266–274. \doi{10.1097/EDE.0000000000000609}
#'
#' Miles, C. H. (2023). On the causal interpretation of randomised interventional
#' indirect effects. \emph{Journal of the Royal Statistical Society: Series B},
#' 85(4), 1154–1172. \doi{10.1093/jrsssb/qkad066}
#'
#' VanderWeele, T. J., & Tchetgen Tchetgen, E. J. (2017).
#' Mediation analysis with time varying exposures and mediators.
#' \emph{Journal of the Royal Statistical Society: Series B}, 79(3), 917–938.
#' \doi{10.1111/rssb.12194}
#'
#' Yamamuro, S., Shinozaki, T., Iimuro, S., & Matsuyama, Y. (2021).
#' Mediational g-formula for time-varying treatment and repeated-measured
#' multiple mediators: Application to atorvastatin's effect on cardiovascular
#' disease via cholesterol lowering and anti-inflammatory actions in elderly
#' type 2 diabetics. \emph{Statistical Methods in Medical Research}, 30(8),
#' 1782–1799. \doi{10.1177/09622802211025988}
#'
#' Zheng, W., & van der Laan, M. (2017).
#' Longitudinal mediation analysis with time-varying mediators and exposures, with application to survival outcomes.
#' \emph{Journal of Causal Inference}, 5(2). \doi{10.1515/jci-2016-0006}
#'
#' @seealso \code{\link{spec_model}}, \code{\link{recodes}},
#'   \code{\link{gformula}} for total effects, and
#'   \code{vignette("causalMed-02-mediation")}.
#'
#' @examples
#' data(nonsurvivaldata)
#'
#' # Models in the order the variables are generated: A -> L1 -> L2 -> M -> Y
#' models <- list(
#'   spec_model(A ~ V + lag1_A + lag1_L1 + lag1_L2 + time,
#'              var_type = "binary", mod_type = "exposure"),
#'   spec_model(L1 ~ V + A + lag1_L1 + time,
#'              var_type = "normal", mod_type = "covariate"),
#'   spec_model(L2 ~ V + A + lag1_L2 + time,
#'              var_type = "binary", mod_type = "covariate"),
#'   spec_model(M ~ V + A + L1 + L2 + lag1_M + time,
#'              var_type = "normal", mod_type = "mediator"),
#'   spec_model(Y_bin ~ V + A + M + L1 + L2,
#'              var_type = "binary", mod_type = "outcome")
#' )
#'
#' fit <- mediation(
#'   data = nonsurvivaldata, id_var = "id", time_var = "time",
#'   base_vars = "V", exposure = "A", outcome = "Y_bin", models = models,
#'   init_recode = recodes(lag1_A = 0, lag1_L1 = 0, lag1_L2 = 0, lag1_M = 0),
#'   in_recode   = recodes(lag1_A = A, lag1_L1 = L1, lag1_L2 = L2, lag1_M = M),
#'   mc_sample = 2000,
#'   R = 1,          # use R > 1 (e.g. 500) for bootstrap confidence intervals
#'   quiet = TRUE
#' )
#' fit
#'
#' @export


mediation <- function(data,
                      id_var,
                      base_vars,
                      exposure,
                      outcome,
                      time_var,
                      models,
                      init_recode = NULL,
                      in_recode = NULL,
                      out_recode = NULL,
                      mc_sample = NULL,
                      mediation_type = c("I", "N"),
                      exposure_regime = 1,
                      reference_regime = 0,
                      n_vw = 2L,
                      estimator = c("gcomp", "tmle"),
                      tmle_weight_trunc = 0.995,
                      return_fitted = FALSE,
                      return_data = FALSE,
                      R = 500,
                      quiet = FALSE,
                      seed = 12345L) {

  tpcall <- match.call()
  all.args <- mget(names(formals()),sys.frame(sys.nframe()))

  # Initilise warning
  init_warn()

  mediation_type <- match.arg(mediation_type)
  estimator      <- match.arg(estimator)
  # all.args was captured before match.arg(), so every match.arg-processed
  # argument must be written back or a defaulted call stores the full
  # candidate vector (print() then errors on the length-2 condition).
  all.args[c("mediation_type", "estimator")] <- list(mediation_type, estimator)

  data <- as.data.table(data)

  if (!is.null(seed)) {
    seed <- as.integer(seed)
    if (exists(".Random.seed", envir = .GlobalEnv)) {
      old_seed <- get(".Random.seed", envir = .GlobalEnv)
      on.exit(assign(".Random.seed", old_seed, envir = .GlobalEnv), add = TRUE)
    } else {
      on.exit(rm(".Random.seed", envir = .GlobalEnv), add = TRUE)
    }
    set.seed(seed)
  }
  boot_seed <- if (!is.null(seed)) seed else TRUE

  # Check for error
  check_error(data, id_var, base_vars, exposure, time_var, models)

  # The simulation grid: the distinct observed time values, sorted. The local
  # `time_seq` (not `time_len`) reaches .run_interventions() and
  # bootstrap_helper() through get_args_for(), so every pass simulates the
  # same steps.
  time_seq <- time_grid(data, time_var)
  time_len <- length(time_seq)

  # Resolved here rather than in the signature so it can count subjects rather
  # than rows. TMLE does no Monte Carlo simulation, so it takes the value
  # silently and never uses it.
  if (is.null(mc_sample)) {
    mc_sample <- 50L * data.table::uniqueN(data[[id_var]])
    if (!quiet && !identical(estimator, "tmle")) {
      message(sprintf(
        "mc_sample not supplied: using %d (50 per subject, %d subjects).",
        mc_sample, data.table::uniqueN(data[[id_var]])))
    }
  }
  mc_sample <- as.integer(mc_sample)
  all.args$mc_sample <- mc_sample

  check_recode_param("in_recode", in_recode)
  check_recode_param("out_recode", out_recode)
  check_recode_param("init_recode", init_recode)

  # `outcome` is supplied separately from `models`, but it must name the same
  # variable the outcome/survival model predicts: it drives the observed
  # benchmark and, for estimator = "tmle", the whole targeted engine. A name
  # that matches no column silently produced an NA benchmark, and made the
  # TMLE binary-outcome guard pass vacuously (all(logical(0)) is TRUE). Pin the
  # two together. check_error() has already established that exactly one
  # outcome/survival model exists.
  out_idx       <- which(sapply(models, function(m)
    m$mod_type %in% c("outcome", "survival")))
  model_outcome <- all.vars(formula(models[[out_idx]]$call)[[2]])
  if (!identical(as.character(outcome), model_outcome)) {
    stop(sprintf(
      "`outcome` (\"%s\") must name the response of the %s model (\"%s\").",
      paste(outcome, collapse = ", "), models[[out_idx]]$mod_type, model_outcome
    ), domain = "causalMed")
  }

  # Validate that the exposure variable is binary {0, 1}
  exp_vals <- unique(na.omit(data[[exposure]]))
  if (!all(exp_vals %in% c(0, 1))) {
    stop(sprintf(
      "Exposure variable '%s' must be binary with values in {0, 1}. Found: {%s}.",
      exposure, paste(sort(exp_vals), collapse = ", ")
    ), domain = "causalMed")
  }

  # ---- Exposure regimes a(1:T) and a*(1:T) ----------------------------------
  # Lin et al. (2017, Stat Med, Section 2.2; Definitions 5-1 to 5-3) and Zheng
  # & van der Laan (2017, Section 2.2) define the effects for two arbitrary
  # STATIC regimes; always/never exposed is their worked example and the
  # default. One value per distinct time point, in sorted order; a scalar is
  # recycled. Recycled here so that regime_key() (pool lookups) sees full
  # vectors and all.args records what was actually simulated.
  exposure_regime  <- check_static_regime(exposure_regime,  time_len,
                                          "exposure_regime",  values = c(0, 1))
  reference_regime <- check_static_regime(reference_regime, time_len,
                                          "reference_regime", values = c(0, 1))
  if (length(exposure_regime)  == 1L) exposure_regime  <- rep(exposure_regime,  time_len)
  if (length(reference_regime) == 1L) reference_regime <- rep(reference_regime, time_len)
  # `check_static_regime()` has already coerced both to double, so `identical()`
  # is a safe equality test here (1L, TRUE and 1 all arrive as 1).
  if (identical(exposure_regime, reference_regime)) {
    stop("`exposure_regime` and `reference_regime` are identical; the contrast is empty.",
         domain = "causalMed")
  }
  all.args[c("exposure_regime", "reference_regime")] <-
    list(exposure_regime, reference_regime)
  is_default_regime <- default_regime_pair(exposure_regime, reference_regime)

  # Identify mediator response variables in temporal (list) order.
  med_idx <- which(sapply(models, function(mods) mods$mod_type == "mediator"))
  if (length(med_idx) == 0L) {
    stop("Mediator model was not defined.", domain = "causalMed")
  }
  med_vars <- vapply(models[med_idx],
                     function(m) all.vars(formula(m$call)[[2]]),
                     character(1))
  N_med <- length(med_vars)

  # Data summary and observed nonparametric benchmark (displayed by print
  # as an informal benchmark next to the simulated intervention means).
  is_survival  <- any(sapply(models, function(mods)
    mods$mod_type %in% c("survival", "censor")))
  data_summary <- summarize_input_data(data, id_var, time_seq)
  observed     <- observed_benchmark(data, outcome, time_var, is_survival)

  # Observed support for the two regimes: a count, shown by print().
  data_summary$regime_support <- regime_support(
    data, id_var, time_var, exposure, time_seq,
    regimes = list("a" = exposure_regime, "a*" = reference_regime))

  # Multi-mediator IDE/IIE is the Yamamuro et al. (2021) extension of the
  # Lin et al. (2017) interventional g-formula; it is defined for
  # mediation_type = "I" only. The natural-effects path (Zheng &
  # van der Laan 2017) is single-mediator.
  if (N_med > 1L && identical(mediation_type, "N")) {
    stop(sprintf(
      "Multiple mediators (%d) are only supported for mediation_type = 'I' (Yamamuro et al. 2021).",
      N_med
    ), domain = "causalMed")
  }

  # Build the intervention list (see build_mediation_interventions()).
  #   mediation_type = "I": fixed-exposure, natural-mediator interventions nat0, nat1 (plug-in TE)
  #     plus permuted-pool interventions Phi00, Phi10, Phi1_k (k=1..N-1),
  #     Phi11  ->  4 + N interventions.
  #   mediation_type = "N" (single mediator): natural-reference interventions
  #     Phi00, Phi11, Phi10, Phi01  ->  4 interventions.
  intervention <- build_mediation_interventions(med_vars, mediation_type,
                                                regime     = exposure_regime,
                                                ref_regime = reference_regime)

  # Warn if the model list order violates the assumed A(t) -> M(t) -> L(t) -> S(t) ordering.
  check_mediation_order(models)

  # For mediation_type = "N", an intermediate (exposure-affected) confounder
  # rules out reading the reported effects as INDIVIDUAL-LEVEL natural effects
  # (the Zheng & van der Laan estimand itself is identified under sequential
  # randomization, their Lemma 1). Detect this by
  # scanning covariate model RHS for the exposure variable. The offending
  # confounder names are retained on the return object so print()/summary()
  # can re-surface a short version of this caveat (it is easy to miss in the
  # runtime warning stream).
  intermediate_confounders <- character(0)
  if (identical(mediation_type, "N")) {
    intermediate_confounders <- check_natural_identifiability(models, exposure)
  }

  # The Monte Carlo "N" engine sets only some exposure-derived columns to the
  # other regime when it draws the cross-world mediator; refuse a mediator
  # model that reads any other (see check_natural_exposure_history()). The
  # TMLE engine applies its own recode rules (.tmle_check_recodes()).
  if (identical(mediation_type, "N") && !identical(estimator, "tmle")) {
    check_natural_exposure_history(models, exposure, init_recode, in_recode,
                                   out_recode)
  }

  # ---- TMLE-specific validation (Zheng & van der Laan 2017, Section 4.3) ----
  if (identical(estimator, "tmle")) {
    if (!identical(mediation_type, "N")) {
      stop("estimator = 'tmle' is only available for mediation_type = 'N' ",
           "(Zheng & van der Laan 2017). The interventional path has no TMLE yet.",
           domain = "causalMed")
    }
    if (!is_default_regime) {
      stop("estimator = 'tmle' supports only the default regimes ",
           "(exposure_regime = 1, reference_regime = 0); use estimator = 'gcomp' ",
           "for other exposure regimes.", domain = "causalMed")
    }
    if (!any(sapply(models, function(m) m$mod_type == "exposure"))) {
      stop("estimator = 'tmle' requires an exposure model (mod_type = 'exposure') ",
           "for the clever-covariate weights.", domain = "causalMed")
    }
    if (any(sapply(models, function(m) !is.null(m$custom_sim)))) {
      stop("estimator = 'tmle' does not support var_type = 'custom' models ",
           "(their conditional density cannot be evaluated).", domain = "causalMed")
    }
    if (any(sapply(models, function(m) !is.null(m$subset)))) {
      stop("estimator = 'tmle' does not support spec_model(subset = ...); ",
           "use estimator = 'gcomp' for subset-restricted models.",
           domain = "causalMed")
    }
    out_vals <- unique(na.omit(data[[outcome]]))
    if (!all(out_vals %in% c(0, 1))) {
      stop("estimator = 'tmle' requires a binary {0, 1} outcome ",
           "(binary endpoint or survival event indicator).", domain = "causalMed")
    }
    # The targeted engine can only honour lag-style recodes; anything else
    # would be silently dropped and distort the estimate. See .tmle_check_recodes().
    .tmle_check_recodes(in_recode, init_recode, out_recode, exposure)
    .tmle_check_model_recodes(models, exposure, in_recode)
    if (R > 1) {
      if (!quiet)
        message("estimator = 'tmle': bootstrap skipped; 95% CIs are Wald ",
                "intervals from the efficient influence curve.")
      R <- 1L
      all.args$R <- 1L
    }
  }

  # Run original estimate
  if (identical(estimator, "tmle")) {
    tmle_res <- tmle_natural_mediation(
      data         = data,
      id_var       = id_var,
      base_vars    = base_vars,
      exposure     = exposure,
      outcome      = outcome,
      time_var     = time_var,
      models       = models,
      in_recode    = in_recode,
      init_recode  = init_recode,
      weight_trunc = tmle_weight_trunc
    )
    est_ori <- list(fitted.models = tmle_res$fit_mods,
                    gform.data    = as.list(tmle_res$psi))
  } else {
    # NOTE: get_args_for() DROPS NULL-valued arguments, so `seed = NULL` would
    # never reach .run_interventions() and its default would silently take over.
    # Force the element through explicitly; `arg_est["seed"] <- list(seed)`
    # keeps a NULL element (plain `$seed <- NULL` would delete it again).
    arg_est <- get_args_for(.run_interventions)
    arg_est["seed"] <- list(seed)
    arg_est$return_fitted <- TRUE
    est_ori <- do.call(.run_interventions, arg_est)
  }

  # Convert each intervention result (scalar Phi, or a data.table when
  # return_data is TRUE) into a plain Phi value for the effect calculations.
  phi_values <- sapply(est_ori$gform.data, phi_scalar)

  # Effect-size table: one row per intervention.
  est_out <- data.table::data.table(
    Intervention = names(phi_values),
    Est          = unname(phi_values)
  )

  # Compute additive and multiplicative effects (TE, IDE, IIE per mediator).
  risk_est <- risk_estimate_mediation(as.list(phi_values), med_vars = med_vars)

  # Point estimates for proportion mediated (additive + multiplicative);
  # the formulas and their references live with the decomposition in
  # pm_from_phi().
  # One proportion row: sum_k IIE(M_k) / OE. See pm_from_phi() for why the
  # former "multiplicative" row was removed (it was the same number).
  pm_point  <- pm_from_phi(phi_values, risk_est)
  pm_keys   <- "pm"
  pm_labels <- "Mediation Proportion"
  pm_est    <- unname(pm_point[pm_keys])

  # Per-replicate bootstrap estimates, retained on the returned object
  # whenever R > 1 (scalar summaries only; size independent of the data).
  boot_interventions <- NULL
  boot_effects       <- NULL

  # Get the mean of bootstrap results
  if (R > 1) {
    # Run bootstrap
    arg_pools <- get_args_for(bootstrap_helper)
    arg_pools$progress_bar <- !quiet
    arg_pools$future_seed <- boot_seed
    pools <- do.call(bootstrap_helper, arg_pools)

    pools_res <- lapply(pools, function(bt) {
      out <- utils::stack(bt$gform.data)
      colnames(out) <- c("Est", "Intervention")
      return(out)
    })
    # Retain the per-replicate intervention means on the returned object.
    # These are scalar summaries only (size independent of the input data).
    boot_interventions <- data.table::rbindlist(pools_res, idcol = "replicate")
    data.table::setcolorder(boot_interventions,
                            c("replicate", "Intervention", "Est"))
    pools_res <- data.table::rbindlist(pools_res)

    # Calculate Sd and percentile confidence interval, both from the finite
    # replicates (print() reports any that were not)
    pools_res <- pools_res[, .(
      Sd = sd(Est, na.rm = TRUE),
      perct_lcl = quantile(Est, 0.025, na.rm = TRUE),
      perct_ucl = quantile(Est, 0.975, na.rm = TRUE)
    ),
    by = c("Intervention")
    ]

    # Merge all and calculate the normal confidence interval
    est_out <- merge(est_out, pools_res, by = c("Intervention"), sort = FALSE)
    est_out <- est_out[, `:=`(
      norm_lcl = Est - stats::qnorm(0.975) * Sd,
      norm_ucl = Est + stats::qnorm(0.975) * Sd
    )]

    # Per-bootstrap Phi values
    boot_phi <- lapply(pools, function(bt) sapply(bt$gform.data, phi_scalar))

    # Per-bootstrap effects tables, each reused for that replicate's
    # proportion-mediated draws (pm_from_phi(); one decomposition per
    # replicate instead of three).
    res_list <- lapply(boot_phi, function(p)
      risk_estimate_mediation(as.list(p), med_vars = med_vars))
    # Always a length(pm_keys) x R matrix. vapply() simplifies a single-row
    # result to a plain vector, which the row indexing below cannot use, so
    # restore the matrix shape explicitly.
    boot_pm <- vapply(seq_along(boot_phi), function(i)
      pm_from_phi(boot_phi[[i]], res_list[[i]]),
      numeric(length(pm_keys)))
    if (!is.matrix(boot_pm)) {
      boot_pm <- matrix(boot_pm, nrow = length(pm_keys),
                        dimnames = list(pm_keys, NULL))
    }

    res_pools <- data.table::rbindlist(res_list, idcol = "replicate")
    # Retain a per-replicate copy before aggregation (PM rows appended below).
    boot_effects <- data.table::copy(res_pools)

    # Append the per-replicate proportion-mediated draws to the retained
    # effects table so users (and print) can see how many replicates were
    # non-finite for each effect.
    boot_effects <- rbind(
      boot_effects,
      data.table::data.table(
        replicate = rep(seq_len(ncol(boot_pm)), times = length(pm_keys)),
        Effect    = rep(pm_labels, each = ncol(boot_pm)),
        RD        = as.vector(t(boot_pm[pm_keys, , drop = FALSE])),
        RR        = NA_real_
      )
    )
    data.table::setorderv(boot_effects, "replicate")

    # Calculate Sd and percentile confidence interval for RD and RR scales
    res_pools <- res_pools[, .(
      Sd           = sd(RD, na.rm = TRUE),
      perct_lcl    = quantile(RD, 0.025, na.rm = TRUE),
      perct_ucl    = quantile(RD, 0.975, na.rm = TRUE),
      Sd_RR        = sd(RR, na.rm = TRUE),
      perct_lcl_RR = quantile(RR, 0.025, na.rm = TRUE),
      perct_ucl_RR = quantile(RR, 0.975, na.rm = TRUE)
    ),
    by = c("Effect")
    ]

    # Merge all and calculate the normal confidence interval (RD and RR)
    risk_est <- merge(risk_est, res_pools, by = c("Effect"), sort = FALSE)
    risk_est <- risk_est[, `:=`(
      norm_lcl    = RD - stats::qnorm(0.975) * Sd,
      norm_ucl    = RD + stats::qnorm(0.975) * Sd,
      norm_lcl_RR = RR - stats::qnorm(0.975) * Sd_RR,
      norm_ucl_RR = RR + stats::qnorm(0.975) * Sd_RR
    )]

    # Proportion rows (additive, multiplicative, and the residual share on the
    # interventional path), each with bootstrap CIs from its own finite draws.
    pm_draws  <- lapply(pm_keys, function(k) {
      v <- boot_pm[k, ]
      v[is.finite(v)]
    })
    pm_sd     <- vapply(pm_draws, function(v)
      if (length(v) > 1) sd(v) else NA_real_, numeric(1))
    pm_q      <- function(p) vapply(pm_draws, function(v)
      if (length(v) > 0) unname(quantile(v, p)) else NA_real_, numeric(1))

    pm_rows <- data.frame(
      Effect       = pm_labels,
      RD           = pm_est,
      RR           = NA_real_,
      Sd           = pm_sd,
      perct_lcl    = pm_q(0.025),
      perct_ucl    = pm_q(0.975),
      norm_lcl     = pm_est - stats::qnorm(0.975) * pm_sd,
      norm_ucl     = pm_est + stats::qnorm(0.975) * pm_sd,
      Sd_RR        = NA_real_,
      perct_lcl_RR = NA_real_,
      perct_ucl_RR = NA_real_,
      norm_lcl_RR  = NA_real_,
      norm_ucl_RR  = NA_real_,
      stringsAsFactors = FALSE
    )
    risk_est <- rbind(risk_est, pm_rows)
  } else {
    # No bootstrap: point estimates only
    risk_est <- rbind(
      risk_est,
      data.frame(
        Effect = pm_labels,
        RD     = pm_est,
        RR     = NA_real_,
        stringsAsFactors = FALSE
      )
    )
  }

  # ---- EIC-based Wald CIs for the TMLE (no bootstrap) ----------------------
  # The influence-curve of each functional gives SE = sd(EIC)/sqrt(n); effect
  # SEs follow by the delta method on per-subject EIC contrasts (Zheng &
  # van der Laan 2017, Corollary 1).
  if (identical(estimator, "tmle")) {
    eic  <- tmle_res$eic
    nsub <- tmle_res$n
    p    <- phi_values
    zq   <- stats::qnorm(0.975)
    se_of <- function(v) stats::sd(v) / sqrt(nsub)

    # Per-intervention SEs
    int_sd <- apply(eic, 2, se_of)
    est_out[, Sd := unname(int_sd[Intervention])]
    est_out[, `:=`(norm_lcl = Est - zq * Sd, norm_ucl = Est + zq * Sd)]

    # RD-scale effect EICs
    eic_nie <- eic[, "Phi11"] - eic[, "Phi10"]
    eic_nde <- eic[, "Phi10"] - eic[, "Phi00"]
    eic_te  <- eic[, "Phi11"] - eic[, "Phi00"]
    # RR-scale delta method: r = pa/pb  =>  EIC_r = (EIC_a - r * EIC_b) / pb
    eic_ratio <- function(nam_a, nam_b) {
      r <- p[[nam_a]] / p[[nam_b]]
      (eic[, nam_a] - r * eic[, nam_b]) / p[[nam_b]]
    }
    # PM = 100 * NIE / OE, delta method. The TMLE path is "N"-only, where there
    # are no nat0/nat1 arms, so OE = Phi11 - Phi00 IS the total effect and the
    # two denominators coincide.
    te  <- p[["Phi11"]] - p[["Phi00"]]
    nie <- p[["Phi11"]] - p[["Phi10"]]
    eic_pm <- if (abs(te) < 1e-10) rep(NA_real_, nsub) else
      100 * (eic_nie * te - nie * eic_te) / te^2

    eff_tbl <- data.table::data.table(
      Effect = c("Indirect effect", "Direct effect", "Total effect", pm_labels),
      Sd     = c(se_of(eic_nie), se_of(eic_nde), se_of(eic_te),
                 if (all(is.na(eic_pm))) NA_real_ else se_of(eic_pm),
                 rep(NA_real_, length(pm_labels) - 1L)),
      Sd_RR  = c(se_of(eic_ratio("Phi11", "Phi10")),
                 se_of(eic_ratio("Phi10", "Phi00")),
                 se_of(eic_ratio("Phi11", "Phi00")),
                 rep(NA_real_, length(pm_labels)))
    )
    risk_est <- data.table::as.data.table(risk_est)
    risk_est <- merge(risk_est, eff_tbl, by = "Effect", sort = FALSE)
    risk_est[, `:=`(
      norm_lcl    = RD - zq * Sd,
      norm_ucl    = RD + zq * Sd,
      norm_lcl_RR = RR - zq * Sd_RR,
      norm_ucl_RR = RR + zq * Sd_RR
    )]
  }

  # Extract fitted model information (attribute-wrapped so print/summary can
  # label each model; shared with gformula()).
  fitted_mods <- wrap_fitted_models(est_ori$fitted.models, return_fitted)

  # Return data
  if (return_data && identical(estimator, "tmle")) {
    warning("return_data is not available for estimator = 'tmle' ",
            "(no Monte Carlo dataset is simulated); returning NULL sim_data.",
            call. = FALSE)
    return_data <- FALSE
  }
  if(return_data){
    dat_out <- data.table::rbindlist(est_ori$gform.data, idcol = "Intervention")
  }else{
    dat_out <- NULL
  }

  emit_warnings()

  boot_estimates <- if (R > 1) {
    list(interventions = boot_interventions, effects = boot_effects)
  } else {
    NULL
  }

  # TMLE diagnostics: per-functional EIC means (should be ~0 at the targeted
  # fit whenever the fluctuations were not skipped for sparse support) and
  # the subject-level influence-curve matrix for custom contrasts.
  tmle_diag <- if (identical(estimator, "tmle")) {
    list(eic_mean = colMeans(tmle_res$eic),
         eic      = tmle_res$eic,
         n        = tmle_res$n)
  } else {
    NULL
  }

  # risk_est[] resets data.table's print flag: the TMLE branch ends with a
  # `:=`, after which the first auto-print of fit$estimate would show nothing.
  y <- list(call = tpcall,
            all.args = all.args,
            estimate = risk_est[],
            effect_size = est_out,
            sim_data = dat_out,
            fitted_models = fitted_mods,
            boot_estimates = boot_estimates,
            data_summary = data_summary,
            observed = observed,
            tmle_diag = tmle_diag,
            intermediate_confounders = intermediate_confounders
          )
  class(y) <- c("gformula", class(y))
  return(y)
}
