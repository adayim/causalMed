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
#' \code{"custom"}). Exactly one \code{"mediator"} model is required.
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
#' For \code{mediation_type = "I"}, a separate cohort is simulated under \eqn{a^* = 0}
#' to build the marginal mediator distribution. In survival analyses the pool is
#' sampled with probability proportional to each individual's cumulative survival
#' \eqn{Sc} under \eqn{a^*}, implementing the \eqn{S(1{:}t-1)=1} condition from
#' Lin et al. (2017, Eq. 4).
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
#' @param mediation_type Character. Type of mediation effect:
#'   \code{"N"} for natural effect or \code{"I"} for interventional effect.
#'
#' @return
#' An object of class \code{"gformula"} with components:
#' \itemize{
  #'   \item \code{call}: the matched function call.
  #'   \item \code{all.args}: a named list of evaluated arguments for reproducibility.
  #'   \item \code{effect_size}: a \code{data.table} with columns \code{Intervention} and \code{Est}.
  #'         If \code{R > 1}, also includes \code{Sd}, percentile CIs (\code{perct_lcl}, \code{perct_ucl}),
  #'         and normal CIs (\code{norm_lcl}, \code{norm_ucl}).Typical labels are:
  #'         \code{Phi11 = E[Y_{1,M(1)}]}, \code{Phi10 = E[Y_{1,M(0)}]},
  #'         \code{Phi00 = E[Y_{0,M(0)}]}, which serve as building blocks for
  #'         natural (or interventional analogue) direct and indirect effects.
  #'   \item \code{estimate}: a \code{data.table} summarizing the decomposition
  #'         of effects, with rows:
  #'         \itemize{
  #'           \item \code{"Indirect effect"} = \eqn{E[Y_{1,M(1)}] - E[Y_{1,M(0)}]}
  #'           \item \code{"Direct effect"}   = \eqn{E[Y_{1,M(0)}] - E[Y_{0,M(0)}]}
  #'           \item \code{"Total effect"}    = \eqn{E[Y_{1,M(1)}] - E[Y_{0,M(0)}]}
  #'           \item \code{"Mediation Proportion"} = Indirect / Total (as a percentage)
  #'         }
  #'         Columns: \code{RD} (risk difference), \code{RR} (risk ratio). When
  #'         \code{R > 1}: also \code{Sd}, \code{perct_lcl}/\code{perct_ucl}
  #'         (percentile CI), \code{norm_lcl}/\code{norm_ucl} (normal CI) for RD;
  #'         and \code{Sd_RR}, \code{perct_lcl_RR}/\code{perct_ucl_RR},
  #'         \code{norm_lcl_RR}/\code{norm_ucl_RR} for RR.
  #'         \code{RR} is \code{NA} for the Mediation Proportion row.
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
#'  data = non-survivaldata,
#'  id_var = "id",
#'  base_vars = "V",
#'  exposure = "A",
#'  outcome = "Y",
#'  time_var = "time",
#'  models = models10,
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
                      mediation_type = c("N", "I"),
                      return_fitted = FALSE,
                      return_data = FALSE,
                      R = 500,
                      quiet = FALSE,
                      seed = mc_sample*100) {

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

  intervention <- list(Phi11 = 1, Phi00 = 0, Phi10 = NULL)

  # Test if mediator model exists and that there is exactly one
  med_flag <- sapply(models, function(mods) mods$mod_type == "mediator")
  if (!any(med_flag)) {
    stop("Mediator model was not defined.", domain = "causalMed")
  }
  if (sum(med_flag) > 1L) {
    stop(sprintf(
      "Only one mediator model is allowed; %d were found. Multiple mediators are not currently supported.",
      sum(med_flag)
    ), domain = "causalMed")
  }

  # Warn if the model list order violates the assumed A(t) -> M(t) -> L(t) -> S(t) ordering.
  check_mediation_order(models)

  # Run original estimate
  arg_est <- get_args_for(.gformula)
  arg_est$return_fitted <- TRUE
  est_ori <- do.call(.gformula, arg_est)

  # Mean value of the outcome at each time point by intervention
  if(return_data){
    est_out <- data.table::rbindlist(est_ori$gform.data, idcol = "Intervention")
    est_out <- est_out[, list(Est = sum(Pred_Y) / length(Pred_Y)), by = c("Intervention")]
  }else{
    est_out <- data.table::as.data.table(utils::stack(est_ori$gform.data))
    colnames(est_out) <- c("Est", "Intervention")
  }

  if(return_data)
    estimate_extract <- sapply(est_ori$gform.data, "[[", "Pred_Y", simplify = FALSE)
  else
    estimate_extract <- est_ori$gform.data

  risk_est <- risk_estimate.mediation(estimate_extract, return_data = return_data)
  
  # Point estimate for Mediation Proportion
  # guard against division by zero when total effect is near zero.
  total_est    <- risk_est$RD[risk_est$Effect == "Total effect"]
  indirect_est <- risk_est$RD[risk_est$Effect == "Indirect effect"]
  med_prop_est <- if (isTRUE(abs(total_est) < 1e-10)) NA_real_ else indirect_est / total_est * 100

  # Get the mean of bootstrap results
  if (R > 1) {
    # Run bootstrap
    arg_pools <- get_args_for(bootstrap_helper)
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

    res_pools <- lapply(pools, function(x) x$gform.data)
    res_pools <- lapply(res_pools, risk_estimate.mediation, return_data = return_data)
    res_pools <- data.table::rbindlist(res_pools)

    # Compute mediation proportion within each bootstrap replicate, then
    # derive CIs from the bootstrap distribution (not ratio of CI endpoints).
    boot_med_prop <- res_pools[Effect == "Indirect effect", RD] /
      res_pools[Effect == "Total effect", RD] * 100

    # Calculate Sd and percentile confidence interval for RD and RR scales
    res_pools <- res_pools[, .(
      Sd        = sd(RD, na.rm = TRUE),
      perct_lcl = quantile(RD, 0.025, na.rm = TRUE),
      perct_ucl = quantile(RD, 0.975, na.rm = TRUE),
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

    # Mediation Proportion with bootstrap-based CIs
    # filter Inf/NaN (produced when bootstrap total effect is near zero)
    # before computing SD and quantiles — na.rm = TRUE does not remove Inf.
    boot_finite <- boot_med_prop[is.finite(boot_med_prop)]
    med_prop_sd <- if (length(boot_finite) > 1) sd(boot_finite) else NA_real_
    risk_est <- rbind(
      risk_est,
      data.frame(
        Effect       = "Mediation Proportion",
        RD           = med_prop_est,
        RR           = NA_real_,
        Sd           = med_prop_sd,
        perct_lcl    = if (length(boot_finite) > 0) as.numeric(quantile(boot_finite, 0.025)) else NA_real_,
        perct_ucl    = if (length(boot_finite) > 0) as.numeric(quantile(boot_finite, 0.975)) else NA_real_,
        norm_lcl     = med_prop_est - stats::qnorm(0.975) * med_prop_sd,
        norm_ucl     = med_prop_est + stats::qnorm(0.975) * med_prop_sd,
        Sd_RR        = NA_real_,
        perct_lcl_RR = NA_real_,
        perct_ucl_RR = NA_real_,
        norm_lcl_RR  = NA_real_,
        norm_ucl_RR  = NA_real_,
        stringsAsFactors = FALSE
      )
    )
  } else {
    # No bootstrap: point estimate only
    risk_est <- rbind(
      risk_est,
      data.frame(
        Effect = "Mediation Proportion",
        RD     = med_prop_est,
        RR     = NA_real_,
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
