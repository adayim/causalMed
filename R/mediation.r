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
#' \code{"mediator"}, \code{"outcome"}, \code{"survival"}, or \code{"censoring"}), and
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
#' \code{\link{spec_model}}. See the vignette section on custom covariate types.
#'
#' **Mediator pool (interventional effects)**
#' For \code{mediation_type = "I"}, a natural-course pass under each treatment
#' level \eqn{a^*} stores every simulated individual's full mediator trajectory
#' \eqn{M(1{:}T)} in a pool matrix. Each decomposition arm that fixes a
#' mediator to its \eqn{a^*} value — \emph{including the reference arms}
#' \code{Phi00} (\eqn{a^* = 0}) and \code{Phi11} (\eqn{a^* = 1}) — permutes the
#' rows of that matrix once and assigns subject \eqn{i} the entire trajectory of
#' pool individual \eqn{\pi(i)}: a joint, stochastic draw \eqn{G_{a^*}} from the
#' simulated distribution of \eqn{M(1{:}T)}, independent of \eqn{i}'s own
#' confounder history (and permuted independently across mediators). This makes
#' \emph{all} interventional arms — references and cross-world alike — use the
#' randomized interventional mediator distribution, as prescribed by Lin et al.
#' (2017, Eq. 4), VanderWeele & Tchetgen Tchetgen (2017), and Yamamuro et al.
#' (2021, Figure 3 step 3), and as implemented by the SAS mGFORMULA macro. The
#' natural plug-in total effect is obtained from separate natural-course arms
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
#'   \code{future_seed} interface of \pkg{future.apply}). Pass \code{NULL} to
#'   leave the global RNG state untouched. Users running multiple analyses in
#'   the same session should set distinct seeds to ensure independent bootstrap
#'   draws across analyses.
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
#'   confounding.
#'
#' @param n_vw Integer. Number of independent permutation draws averaged for
#'   each interventional cross-world arm. The same averaging is applied
#'   within every bootstrap replicate. Default \code{2L} to match the SAS
#'   \code{mGFORMULA} macro's \code{n_vw = 2}. Set to \code{1L} to disable
#'   averaging (faster, slightly noisier Monte Carlo). Has no effect when
#'   \code{mediation_type = "N"} (the natural-effect mediator swap is not
#'   permutation-based).
#'
#' @return
#' An object of class \code{"gformula"} with components:
#' \itemize{
  #'   \item \code{call}: the matched function call.
  #'   \item \code{all.args}: a named list of evaluated arguments for reproducibility.
  #'   \item \code{effect_size}: a \code{data.table} with one row per
  #'         simulated arm (columns \code{Intervention} and \code{Est}).
  #'         For \code{mediation_type = "I"} the mediator(s) in every
  #'         decomposition arm are drawn from independently-permuted marginal
  #'         pools (stochastic draws \eqn{G}): \code{Phi00} (= \eqn{E[Y_{0,G_0}]})
  #'         and \code{Phi11} (= \eqn{E[Y_{1,G_1}]}) are the interventional
  #'         reference arms, \code{Phi10} (= \eqn{E[Y_{1,G_0}]}) is the
  #'         cross-world arm, and for \eqn{N \ge 2} mediators \code{Phi1_k}
  #'         (k = 1, …, N-1) capture the sequential per-mediator transition.
  #'         Two additional natural-course arms \code{nat0} (= \eqn{E[Y_0]}) and
  #'         \code{nat1} (= \eqn{E[Y_1]}) give the plug-in total effect. For
  #'         \code{mediation_type = "N"} the arms \code{Phi00}/\code{Phi11} are
  #'         the natural never-/always-treat arms (no permutation). With
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
  #'                  on the interventional overall scale (Lin et al. 2017, Table 2)
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

  # Get time length
  time_len <- length(unique(na.omit(data[[time_var]])))

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

  # Build the arm list (see build_mediation_arms()).
  #   mediation_type = "I": natural-course arms nat0, nat1 (for the plug-in TE)
  #     plus permuted-pool interventional arms Phi00, Phi10, Phi1_k (k=1..N-1),
  #     Phi11  ->  4 + N arms.
  #   mediation_type = "N" (single mediator): natural-reference arms
  #     Phi00, Phi11, Phi10, Phi01  ->  4 arms.
  intervention <- build_mediation_arms(med_vars, mediation_type)

  # Warn if the model list order violates the assumed A(t) -> M(t) -> L(t) -> S(t) ordering.
  check_mediation_order(models)

  # For mediation_type = "N", natural direct/indirect effects are not identifiable
  # when an intermediate (exposure-affected) confounder is present. Detect this by
  # scanning covariate model RHS for the exposure variable.
  if (identical(mediation_type, "N")) {
    check_natural_identifiability(models, exposure)
  }

  # Run original estimate
  arg_est <- get_args_for(.gformula)
  arg_est$return_fitted <- TRUE
  est_ori <- do.call(.gformula, arg_est)

  # Convert each arm result (scalar Phi, or a data.table when return_data is
  # TRUE) into a plain Phi value for the effect calculations.
  arm_to_phi <- function(x) {
    if (is.null(x)) return(NA_real_)
    if (data.table::is.data.table(x) || is.data.frame(x)) {
      # Prefer the attached attribute (set by .gformula when n_vw > 1 so the
      # scalar reflects the average, not the last permutation alone).
      attr_phi <- attr(x, "phi_estimate", exact = TRUE)
      if (!is.null(attr_phi)) return(as.numeric(attr_phi))
      return(sum(x[["Pred_Y"]]) / length(x[["Pred_Y"]]))
    }
    as.numeric(x)
  }
  phi_values <- sapply(est_ori$gform.data, arm_to_phi)

  # Effect-size table: one row per arm.
  est_out <- data.table::data.table(
    Intervention = names(phi_values),
    Est          = unname(phi_values)
  )

  # Compute additive and multiplicative effects (TE, IDE, IIE per mediator).
  risk_est <- risk_estimate.mediation(as.list(phi_values), med_vars = med_vars)

  # Point estimates for proportion mediated (additive + multiplicative).
  # Additive PM follows Yamamuro et al. (2021): the proportion of the total
  # effect explained by the mediators is (sum IIE + residual) / TE, which
  # equals (TE - IDE) / TE. The residual (TE - OE) is folded into the mediated
  # portion. For mediation_type = "N" the residual is 0, so this reduces to the
  # usual sum(IIE)/TE.
  total_est  <- risk_est$RD[risk_est$Effect == "Total effect"]   # natural TE
  direct_est <- risk_est$RD[risk_est$Effect == "Direct effect"]
  rr_ide     <- risk_est$RR[risk_est$Effect == "Direct effect"]
  rr_iie_k   <- risk_est$RR[grepl("^Indirect effect", risk_est$Effect)]

  # Interventional overall-effect risk ratio (the scale on which the
  # multiplicative decomposition RR_OE = RR_IDE * prod(RR_IIE_k) is exact).
  rr_oe <- as.numeric(phi_values[["Phi11"]] / phi_values[["Phi00"]])

  med_prop_add <- if (isTRUE(abs(total_est) < 1e-10)) {
    NA_real_
  } else {
    (total_est - direct_est) / total_est * 100
  }

  # Multiplicative proportion mediated (Lin et al. 2017 Table 2, generalised to
  # N mediators via Yamamuro et al. 2021), on the interventional overall scale:
  #   PM_mult = RR_IDE * (prod(RR_IIE_k) - 1) / (RR_OE - 1) * 100
  med_prop_mult <- if (isTRUE(abs(rr_oe - 1) < 1e-10) ||
                       !all(is.finite(c(rr_oe, rr_ide, rr_iie_k)))) {
    NA_real_
  } else {
    rr_ide * (prod(rr_iie_k) - 1) / (rr_oe - 1) * 100
  }

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
    boot_phi <- lapply(pools, function(bt) sapply(bt$gform.data, arm_to_phi))

    # Per-bootstrap effects table
    res_pools <- lapply(boot_phi,
                        function(p) risk_estimate.mediation(as.list(p),
                                                            med_vars = med_vars))
    res_pools <- data.table::rbindlist(res_pools)

    # Per-bootstrap proportion mediated (additive + multiplicative).
    # Computed per replicate from boot_phi (the rbindlist above lost the
    # replicate index, so we recompute from the per-replicate Phi map).
    boot_pm_add  <- vapply(boot_phi, function(p) {
      one <- as.list(p)
      r <- risk_estimate.mediation(one, med_vars = med_vars)
      total  <- r$RD[r$Effect == "Total effect"]
      direct <- r$RD[r$Effect == "Direct effect"]
      if (isTRUE(abs(total) < 1e-10)) NA_real_ else (total - direct) / total * 100
    }, numeric(1))
    boot_pm_mult <- vapply(boot_phi, function(p) {
      one <- as.list(p)
      r <- risk_estimate.mediation(one, med_vars = med_vars)
      rr_oe <- as.numeric(p[["Phi11"]] / p[["Phi00"]])   # interventional overall RR
      rr_d  <- r$RR[r$Effect == "Direct effect"]
      rr_i  <- r$RR[grepl("^Indirect effect", r$Effect)]
      if (isTRUE(abs(rr_oe - 1) < 1e-10) ||
          !all(is.finite(c(rr_oe, rr_d, rr_i)))) NA_real_
      else rr_d * (prod(rr_i) - 1) / (rr_oe - 1) * 100
    }, numeric(1))

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
  
  # Extract fitted model information. Wrap with the same attribute set as
  # gformula() so that print.gformula(models = TRUE) and summary.gformula()
  # can display var_type / mod_type / recodes / subset labels.
  resp_vars_list <- sapply(est_ori$fitted.models, function(x) {
    x$rsp_vars
  })
  if (return_fitted) {
    fitted_mods <- lapply(est_ori$fitted.models, function(x) {
      structure(x$fitted,
                recodes = x$recodes,
                subset = x$subset,
                var_type = x$var_type,
                mod_type = x$mod_type)
    })
  } else {
    fitted_mods <- lapply(est_ori$fitted.models, function(x) {
      r <- list(call = x$fitted$call,
                coeff = summary(x$fitted)$coefficients)
      structure(r,
                recodes = x$recodes,
                subset = x$subset,
                var_type = x$var_type,
                mod_type = x$mod_type)
    })
  }
  names(fitted_mods) <- resp_vars_list

  # Return data
  if(return_data){
    dat_out <- data.table::rbindlist(est_ori$gform.data, idcol = "Intervention")
  }else{
    dat_out <- NULL
  }

  emit_warnings()

  y <- list(call = tpcall,
            all.args = all.args,
            estimate = risk_est,
            effect_size = est_out,
            sim_data = dat_out,
            fitted_models = fitted_mods
          )
  class(y) <- c("gformula", class(y))
  return(y)
}
