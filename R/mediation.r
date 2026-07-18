#' G-formula based mediation analysis
#'
#' @description
#' Conduct mediation analysis with time-varying mediators using the g-formula.
#' This function estimates total effect, natural direct effect, and natural indirect
#' effect for both natural mediation (\code{"N"}) and interventional mediation (\code{"I"}).
#' A \code{data.frame} summarizing the estimates will be returned.
#'
#' @details
#'
#' **Data requirements**
#' The input dataset must be in longitudinal ("long") format: one record per subject per
#' time period. In survival settings, all records after the event must be removed.
#'
#' The exposure variable must be **binary, coded as 0 (reference/untreated) and 1
#' (active/treated)**. Methods implemented here (Lin et al. 2017; Zheng & van der Laan 2017)
#' are defined for binary exposures only.
#'
#' **Model specification**
#' Each element of \code{models} is created by \code{\link{spec_model}} and must specify:
#' (i) a formula, (ii) \code{mod_type} (\code{"exposure"}, \code{"covariate"},
#' \code{"mediator"}, \code{"outcome"}, \code{"survival"}, or \code{"censor"}), and
#' (iii) \code{var_type} (\code{"binary"}, \code{"normal"}, \code{"categorical"}, or
#' \code{"custom"}). At least one \code{"mediator"} model is required.
#' Multiple mediator models are supported under \code{mediation_type = "I"}
#' (Yamamuro et al. 2021); list them in temporal order. Multi-mediator
#' analyses are not supported under \code{mediation_type = "N"}.
#'
#' The list order determines the simulation sequence at each time step and must match
#' your assumed data-generating process. A common ordering is
#' \strong{A(t) -> L(t) -> M(t) -> S(t)} (confounders not affected by mediator) or
#' \strong{A(t) -> M(t) -> L(t) -> S(t)} (confounders affected by both exposure and
#' mediator). The function checks that the mediator precedes the outcome and that exposure
#' precedes the mediator, and warns if these are violated.
#'
#' For non-standard covariate distributions (bounded, zero-inflated, truncated), use
#' \code{var_type = "custom"} with \code{custom_fit} and \code{custom_sim} arguments to
#' \code{\link{spec_model}}. See the custom covariate distributions section of
#' \code{vignette("causalMed-gformula")}.
#'
#' **Mediator pool (interventional effects)**
#' For \code{mediation_type = "I"}, a natural-course pass under each treatment
#' level \eqn{a^*} stores every simulated individual's full mediator trajectory
#' \eqn{M(1{:}T)} in a pool matrix. Each decomposition intervention that fixes a
#' mediator to its \eqn{a^*} value — \emph{including the reference interventions}
#' \code{Phi00} (\eqn{a^* = 0}) and \code{Phi11} (\eqn{a^* = 1}) — permutes the
#' rows of that matrix once and assigns subject \eqn{i} the entire trajectory of
#' pool individual \eqn{\pi(i)}: a joint, stochastic draw \eqn{G_{a^*}} from the
#' simulated distribution of \eqn{M(1{:}T)}, independent of \eqn{i}'s own
#' confounder history (and permuted independently across mediators). This makes
#' \emph{all} interventions in the decomposition — references and cross-world alike — use the
#' randomized interventional mediator distribution, as prescribed by Lin et al.
#' (2017, Eq. 4), VanderWeele & Tchetgen Tchetgen (2017), and Yamamuro et al.
#' (2021, Figure 3 step 3), and as implemented by the SAS mGFORMULA macro. The
#' natural plug-in total effect is obtained from separate natural-course interventions
#' (\code{nat0}, \code{nat1}), so \eqn{TE} need not equal the sum of the
#' interventional direct and indirect effects; their difference is reported as
#' the mediated-interaction residual. The pool is not survival-weighted; the
#' full reference cohort is used, mirroring the reference SAS implementations.
#'
#' **Warnings**
#' Warnings from model fitting (e.g., convergence, near-separation) are collected
#' during the run and printed as a deduplicated summary at function exit, including
#' a repeat count when the same message fires across bootstrap replicates.
#'
#' **Re-coding hooks**
#' \itemize{
#'   \item \code{init_recode}: executed once at time 0 before simulation (initialize baselines).
#'   \item \code{in_recode}: executed at the start of each time step (e.g., entry-time logic).
#'   \item \code{out_recode}: executed at the end of each time step (e.g., cumulative counts, lags).
#' }
#'
#' @inheritParams gformula
#'
#' @param outcome Character scalar. Name of the outcome variable in \code{data}.
#'   Must match the response variable in the outcome or survival model.
#' @param seed Integer random seed for reproducibility. Default \code{12345L}.
#'   Setting \code{seed} fixes the RNG stream used for the original-data Monte
#'   Carlo simulation \emph{and} the bootstrap replicates (the latter via the
#'   \code{future_seed} interface of \pkg{future.apply}); the caller's global
#'   RNG state is saved and restored on exit, so a seeded call does not disturb
#'   the ambient random stream. Pass \code{NULL} to disable seeding entirely:
#'   no \code{set.seed()} is called, the Monte Carlo draws consume and advance
#'   the ambient RNG stream, and repeated calls therefore give \emph{different}
#'   results (reproducible only via an outer \code{set.seed()}). Use
#'   \code{seed = NULL} inside simulation loops that manage their own seeds.
#'   Users running multiple seeded analyses in the same session should set
#'   distinct seeds to ensure independent bootstrap draws across analyses.
#' @param mediation_type Character. Type of mediation effect:
#'   \code{"I"} for interventional effect (default) or \code{"N"} for natural
#'   effect. The default was chosen because the interventional estimand
#'   remains identifiable in the general time-varying-confounder setting the
#'   package is designed for (see below).
#'   \strong{Identifiability.} Natural direct and indirect effects
#'   (\code{"N"}) are \emph{not identifiable} from observational data when
#'   there exists a time-varying confounder of the mediator-outcome
#'   relationship that is itself affected by prior exposure (Avin, Shpitser
#'   & Pearl 2005; VanderWeele 2014; VanderWeele & Tchetgen Tchetgen 2017).
#'   In that setting use \code{mediation_type = "I"}, which targets
#'   randomized interventional analogues of the direct/indirect effects and
#'   remains identifiable. A warning is emitted at runtime if \code{"N"} is
#'   chosen and any covariate model includes the exposure on its right-hand
#'   side, which is a sufficient (though not necessary) signal of intermediate
#'   confounding; the caveat is repeated by \code{print()}.
#'   \strong{Interpretation of \code{"I"}.} Randomized interventional indirect
#'   effects buy identifiability at a price: they do not satisfy a sharp
#'   mediational null criterion, so a non-zero IIE does not by itself prove
#'   that the mediator transmits the effect for any individual (Miles 2023).
#'   The choice between \code{"I"} and \code{"N"} is therefore a substantive
#'   one, not merely computational.
#'
#' @param n_vw Integer. Number of independent permutation draws averaged for
#'   each interventional cross-world intervention. The same averaging is applied
#'   within every bootstrap replicate. Default \code{2L} to match the SAS
#'   \code{mGFORMULA} macro's \code{n_vw = 2}. Set to \code{1L} to disable
#'   averaging (faster, slightly noisier Monte Carlo). Has no effect when
#'   \code{mediation_type = "N"} (the natural-effect mediator swap is not
#'   permutation-based).
#'
#' @param estimator Character. \code{"gcomp"} (default) for the parametric
#'   g-formula plug-in (Monte Carlo simulation + bootstrap CIs), or
#'   \code{"tmle"} for the targeted maximum likelihood estimator of Zheng &
#'   van der Laan (2017, Section 4.3). \code{"tmle"} is available for
#'   \code{mediation_type = "N"} only (single mediator, binary \{0, 1\}
#'   outcome or survival event indicator; an exposure model is required and
#'   \code{var_type = "custom"} / \code{spec_model(subset = )} are not
#'   supported). \strong{Recode restriction:} the targeted engine evaluates
#'   only lag-style recodes, so \code{in_recode} entries must copy a single
#'   column (e.g. \code{recodes(lag_A = A)}), lags of the exposure must copy
#'   the exposure itself (chained exposure lags such as
#'   \code{recodes(lag2_A = lag_A)} are rejected; deeper exposure history
#'   requires \code{estimator = "gcomp"}), \code{init_recode} entries must
#'   be a constant or a single column name, and \code{out_recode} is not
#'   supported; derived recodes (splines, cumulative counts, carry-forward
#'   flags) require \code{estimator = "gcomp"}. Violations are rejected with
#'   an error rather than silently ignored. The TMLE uses backward iterated
#'   regressions with logistic
#'   fluctuations instead of forward simulation: it is multiply robust
#'   (consistent when specific subsets of the nuisance models are correct,
#'   not only when all are), reports Wald CIs from the efficient influence
#'   curve without bootstrapping (\code{R} is ignored), and has no Monte
#'   Carlo simulation error (\code{mc_sample} is ignored). Practical
#'   positivity violations (few subjects following an intervened regime) are
#'   handled by skipping the affected fluctuation steps and truncating
#'   extreme weights, with collected warnings — inspect these before trusting
#'   the affected functionals.
#'
#' @param tmle_weight_trunc Numeric in (0, 1]. Quantile at which each
#'   clever-covariate weight vector is truncated in the TMLE fluctuation
#'   steps (positivity protection). Default \code{0.995}; \code{1} disables
#'   truncation. Ignored for \code{estimator = "gcomp"}.
#'
#' @return
#' An object of class \code{"gformula"} with components:
#' \itemize{
  #'   \item \code{call}: the matched function call.
  #'   \item \code{all.args}: a named list of evaluated arguments for reproducibility.
  #'   \item \code{effect_size}: a \code{data.table} with one row per
  #'         simulated intervention (columns \code{Intervention} and \code{Est}).
  #'         For \code{mediation_type = "I"} the mediator(s) in every
  #'         decomposition intervention are drawn from independently-permuted marginal
  #'         pools (stochastic draws \eqn{G}): \code{Phi00} (= \eqn{E[Y_{0,G_0}]})
  #'         and \code{Phi11} (= \eqn{E[Y_{1,G_1}]}) are the interventional
  #'         reference interventions, \code{Phi10} (= \eqn{E[Y_{1,G_0}]}) is the
  #'         cross-world intervention, and for \eqn{N \ge 2} mediators \code{Phi1_k}
  #'         (k = 1, …, N-1) capture the sequential per-mediator transition.
  #'         Two additional natural-course interventions \code{nat0} (= \eqn{E[Y_0]}) and
  #'         \code{nat1} (= \eqn{E[Y_1]}) give the plug-in total effect. For
  #'         \code{mediation_type = "N"} the interventions \code{Phi00}/\code{Phi11} are
  #'         the natural never-/always-treat interventions (no permutation). With
  #'         \code{R > 1} the table also carries \code{Sd},
  #'         \code{perct_lcl}/\code{perct_ucl}, and \code{norm_lcl}/\code{norm_ucl}.
  #'   \item \code{estimate}: a \code{data.table} summarizing the decomposition
  #'         of effects. For a single mediator under \code{mediation_type = "I"}
  #'         the rows are:
  #'         \itemize{
  #'           \item \code{"Indirect effect"} = \eqn{E[Y_{1,G_1}] - E[Y_{1,G_0}]}
  #'           \item \code{"Direct effect"}   = \eqn{E[Y_{1,G_0}] - E[Y_{0,G_0}]}
  #'           \item \code{"Total effect"}    = \eqn{E[Y_1] - E[Y_0]} (natural
  #'                 plug-in g-formula; \emph{not} the sum of the components)
  #'           \item \code{"TE - (Direct + Indirect)"} = the mediated-interaction
  #'                 residual \eqn{TE - (IDE + IIE)} = \eqn{TE} minus the
  #'                 interventional overall effect \eqn{E[Y_{1,G_1}] - E[Y_{0,G_0}]}
  #'                 (generally non-zero for interventional effects; this row is
  #'                 \emph{absent} for \code{mediation_type = "N"}, where the
  #'                 decomposition sums exactly to \eqn{TE}).
  #'           \item \code{"Mediation Proportion"}
  #'                  = \eqn{(TE - IDE) / TE \times 100\%}
  #'                  (= \eqn{(IIE + residual)/TE}; Yamamuro et al. 2021)
  #'           \item \code{"Mediation Proportion (multiplicative)"}
  #'                  = \eqn{RR_{IDE}\,(RR_{IIE}-1)/(RR_{OE}-1) \times 100\%}
  #'                  on the interventional overall scale (Lin et al. 2017, Table 2).
  #'                  Under \code{mediation_type = "N"} the same formula is
  #'                  reported on the total-effect scale as an \emph{analogue};
  #'                  Zheng & van der Laan (2017) do not define a proportion
  #'                  mediated.
  #'         }
  #'         For multiple mediators each indirect effect is labelled
  #'         \code{"Indirect effect (<mediator>)"}; the additive proportion is
  #'         \eqn{(TE - IDE)/TE} and the multiplicative is
  #'         \eqn{RR_{IDE}\,(\prod_k RR_{IIE_k} - 1)/(RR_{OE}-1)}
  #'         (Yamamuro et al. 2021). Columns: \code{RD}, \code{RR}; with
  #'         \code{R > 1} also \code{Sd}, percentile CIs
  #'         (\code{perct_lcl}/\code{perct_ucl}), and normal CIs
  #'         (\code{norm_lcl}/\code{norm_ucl}) for RD; and \code{Sd_RR},
  #'         \code{perct_lcl_RR}/\code{perct_ucl_RR},
  #'         \code{norm_lcl_RR}/\code{norm_ucl_RR} for RR. \code{RR} is
  #'         \code{NA} for the Mediation Proportion rows.
  #'   \item \code{sim_data}: if \code{return_data = TRUE}, the Monte Carlo simulated
  #'         longitudinal dataset used internally (can be large).
  #'   \item \code{fitted_models}: a named list of fitted models keyed by outcome, exposure,
  #'         and mediator variables. If \code{return_fitted = TRUE}, returns full model objects
  #'         plus attributes (\code{recodes}, \code{subset}, \code{var_type}, \code{mod_type});
  #'         otherwise, a compact list with \code{call} and \code{coeff}.
  #'   \item \code{boot_estimates}: when \code{R > 1}, a list of per-replicate
  #'         bootstrap estimates: \code{$interventions} (columns
  #'         \code{replicate}, \code{Intervention}, \code{Est}) and
  #'         \code{$effects} (columns \code{replicate}, \code{Effect},
  #'         \code{RD}, \code{RR}; includes the per-replicate Mediation
  #'         Proportion draws). Useful for diagnostics such as counting
  #'         non-finite replicates or computing custom intervals. These are
  #'         scalar summaries whose size is independent of the input data.
  #'         \code{NULL} when \code{R <= 1}.
  #'   \item \code{data_summary}: list with the number of individuals
  #'         (\code{n_id}), observations (\code{n_obs}), and time points
  #'         (\code{n_times}, \code{t_min}, \code{t_max}) of the input data.
  #'   \item \code{observed}: list with the observed nonparametric benchmark
  #'         of the outcome (\code{value}, \code{label}): the mean outcome at
  #'         the last time point, or the product-limit cumulative incidence
  #'         for survival outcomes. Printed next to the simulated
  #'         intervention means as an informal benchmark.
  #'   \item \code{tmle_diag}: for \code{estimator = "tmle"} only, a list with
  #'         the subject-level efficient-influence-curve matrix (\code{eic},
  #'         one column per functional), its column means (\code{eic_mean},
  #'         near zero for well-supported functionals at the targeted fit),
  #'         and the number of subjects (\code{n}). \code{NULL} for
  #'         \code{estimator = "gcomp"}.
  #'   \item \code{intermediate_confounders}: for \code{mediation_type = "N"},
  #'         a character vector of covariate names whose model includes the
  #'         exposure (exposure-affected/intermediate confounders that make
  #'         the natural effects non-identifiable); empty when none. The
  #'         \code{print} method re-surfaces a short identifiability caveat
  #'         when this is non-empty.
  #' }
#'
#' @references
#' Lin, S. H., Young, J. G., Logan, R., & VanderWeele, T. J. (2017).
#' Mediation analysis for a survival outcome with time-varying exposures, mediators, and confounders.
#' \emph{Statistics in Medicine}, 36(26), 4153–4166. \doi{10.1002/sim.7426}
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
#' Miles, C. H. (2023). On the causal interpretation of randomised interventional
#' indirect effects. \emph{Journal of the Royal Statistical Society: Series B},
#' 85(4), 1154–1172. \doi{10.1093/jrsssb/qkad066}
#'
#' @examples
#' \dontrun{
#' data(nonsurvivaldata)
#'
#' models <- list(
#'   cov_model = spec_model(L ~ V+L_lag1+A+time,var_type= "normal",mod_type = "covariate"),
#'   mediator_model = spec_model(M ~ V + A + L + M_lag1 + time,
#'                               var_type = "normal", mod_type = "mediator"),
#'   outcome_model = spec_model(Y ~ V+A+M+A * M+L,  var_type= "binary",mod_type ="outcome")
#'   )
#'
#' fit <- mediation(
#'  data = nonsurvivaldata,
#'  id_var = "id",
#'  base_vars = "V",
#'  exposure = "A",
#'  outcome = "Y",
#'  time_var = "time",
#'  models = models,
#'  init_recode = recodes(M_lag1=0,L_lag1=0),
#'  in_recode = recodes(M_lag1=M,L_lag1=L),
#'  mediation_type = "I",
#'  mc_sample = 100000,
#'  R = 500,
#'  return_data = FALSE,
#'  return_fitted = FALSE
#'  )
#'
#' print(fit)
#' }
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
                      mc_sample = nrow(data)*100,
                      mediation_type = c("I", "N"),
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

  check_recode_param("in_recode", in_recode)
  check_recode_param("out_recode", out_recode)
  check_recode_param("init_recode", init_recode)

  # Validate that the exposure variable is binary {0, 1}
  exp_vals <- unique(na.omit(data[[exposure]]))
  if (!all(exp_vals %in% c(0, 1))) {
    stop(sprintf(
      "Exposure variable '%s' must be binary with values in {0, 1}. Found: {%s}.",
      exposure, paste(sort(exp_vals), collapse = ", ")
    ), domain = "causalMed")
  }

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
  data_summary <- summarize_input_data(data, id_var, time_var)
  observed     <- observed_benchmark(data, outcome, time_var, is_survival)

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
  #   mediation_type = "I": natural-course interventions nat0, nat1 (for the plug-in TE)
  #     plus permuted-pool interventions Phi00, Phi10, Phi1_k (k=1..N-1),
  #     Phi11  ->  4 + N interventions.
  #   mediation_type = "N" (single mediator): natural-reference interventions
  #     Phi00, Phi11, Phi10, Phi01  ->  4 interventions.
  intervention <- build_mediation_interventions(med_vars, mediation_type)

  # Warn if the model list order violates the assumed A(t) -> M(t) -> L(t) -> S(t) ordering.
  check_mediation_order(models)

  # For mediation_type = "N", natural direct/indirect effects are not identifiable
  # when an intermediate (exposure-affected) confounder is present. Detect this by
  # scanning covariate model RHS for the exposure variable. The offending
  # confounder names are retained on the return object so print()/summary()
  # can re-surface a short version of this caveat (it is easy to miss in the
  # runtime warning stream).
  intermediate_confounders <- character(0)
  if (identical(mediation_type, "N")) {
    intermediate_confounders <- check_natural_identifiability(models, exposure)
  }

  # ---- TMLE-specific validation (Zheng & van der Laan 2017, Section 4.3) ----
  if (identical(estimator, "tmle")) {
    if (!identical(mediation_type, "N")) {
      stop("estimator = 'tmle' is only available for mediation_type = 'N' ",
           "(Zheng & van der Laan 2017). The interventional path has no TMLE yet.",
           domain = "causalMed")
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
  pm_point      <- pm_from_phi(phi_values, risk_est)
  med_prop_add  <- pm_point[["add"]]
  med_prop_mult <- pm_point[["mult"]]

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

    # Calculate Sd and percentile confidence interval
    pools_res <- pools_res[, .(
      Sd = sd(Est),
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
    boot_pm  <- mapply(pm_from_phi, boot_phi, res_list)   # 2 x R matrix
    boot_pm_add  <- boot_pm["add", ]
    boot_pm_mult <- boot_pm["mult", ]

    res_pools <- data.table::rbindlist(res_list, idcol = "replicate")
    # Retain a per-replicate copy before aggregation (PM rows appended below).
    boot_effects <- data.table::copy(res_pools)

    # Append the per-replicate proportion-mediated draws to the retained
    # effects table so users (and print) can see how many replicates were
    # non-finite for each effect.
    boot_effects <- rbind(
      boot_effects,
      data.table::data.table(
        replicate = rep(seq_along(boot_pm_add), times = 2L),
        Effect    = rep(c("Mediation Proportion",
                          "Mediation Proportion (multiplicative)"),
                        each = length(boot_pm_add)),
        RD        = c(boot_pm_add, boot_pm_mult),
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

    # Proportion mediated rows: additive + multiplicative, with bootstrap CIs.
    pm_finite_add  <- boot_pm_add[is.finite(boot_pm_add)]
    pm_finite_mult <- boot_pm_mult[is.finite(boot_pm_mult)]
    pm_sd_add  <- if (length(pm_finite_add)  > 1) sd(pm_finite_add)  else NA_real_
    pm_sd_mult <- if (length(pm_finite_mult) > 1) sd(pm_finite_mult) else NA_real_

    pm_rows <- data.frame(
      Effect       = c("Mediation Proportion",
                       "Mediation Proportion (multiplicative)"),
      RD           = c(med_prop_add, med_prop_mult),
      RR           = c(NA_real_, NA_real_),
      Sd           = c(pm_sd_add, pm_sd_mult),
      perct_lcl    = c(if (length(pm_finite_add)  > 0) quantile(pm_finite_add,  0.025) else NA_real_,
                       if (length(pm_finite_mult) > 0) quantile(pm_finite_mult, 0.025) else NA_real_),
      perct_ucl    = c(if (length(pm_finite_add)  > 0) quantile(pm_finite_add,  0.975) else NA_real_,
                       if (length(pm_finite_mult) > 0) quantile(pm_finite_mult, 0.975) else NA_real_),
      norm_lcl     = c(med_prop_add  - stats::qnorm(0.975) * pm_sd_add,
                       med_prop_mult - stats::qnorm(0.975) * pm_sd_mult),
      norm_ucl     = c(med_prop_add  + stats::qnorm(0.975) * pm_sd_add,
                       med_prop_mult + stats::qnorm(0.975) * pm_sd_mult),
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
        Effect = c("Mediation Proportion",
                   "Mediation Proportion (multiplicative)"),
        RD     = c(med_prop_add, med_prop_mult),
        RR     = c(NA_real_, NA_real_),
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
    # Additive PM = 100 * NIE / TE  =>  delta method
    te  <- p[["Phi11"]] - p[["Phi00"]]
    nie <- p[["Phi11"]] - p[["Phi10"]]
    eic_pm <- if (abs(te) < 1e-10) rep(NA_real_, nsub) else
      100 * (eic_nie * te - nie * eic_te) / te^2

    eff_tbl <- data.table::data.table(
      Effect = c("Indirect effect", "Direct effect", "Total effect",
                 "Mediation Proportion",
                 "Mediation Proportion (multiplicative)"),
      Sd     = c(se_of(eic_nie), se_of(eic_nde), se_of(eic_te),
                 if (all(is.na(eic_pm))) NA_real_ else se_of(eic_pm),
                 NA_real_),
      Sd_RR  = c(se_of(eic_ratio("Phi11", "Phi10")),
                 se_of(eic_ratio("Phi10", "Phi00")),
                 se_of(eic_ratio("Phi11", "Phi00")),
                 NA_real_, NA_real_)
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

  y <- list(call = tpcall,
            all.args = all.args,
            estimate = risk_est,
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
