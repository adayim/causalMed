
#' Parametric g-formula for time-varying interventions
#'
#' @description
#' Estimates the mean outcome, or the risk for a survival outcome, that would
#' be seen if everyone followed each of several exposure strategies, and the
#' contrasts between them (the total effect), using the parametric g-formula
#' (Robins 1986). The models in \code{models} are fitted to the observed data
#' and then used to simulate each subject forward in time under every
#' intervention.
#'
#' @details
#' \strong{Data.} Long format: one row per subject per time point. For a
#' survival outcome the data must be in risk-set form, with every row after the
#' event and after loss to follow-up removed. A \code{"censor"} model may be
#' declared; it is simulated in the natural course but does not change any
#' reported risk, each of which is the risk under eliminated loss to
#' follow-up.
#'
#' \strong{Models.} \code{models} is a list of \code{\link{spec_model}}
#' objects in the order the variables are generated within a time point. The
#' \code{"outcome"} or \code{"survival"} model gives the predicted outcome
#' (\code{Pred_Y}) under each intervention: the fitted mean or risk of an
#' end-of-follow-up outcome, or the cumulative risk built from the fitted
#' hazards.
#'
#' \strong{Interventions.} \code{intervention} is a named list, for example
#' \code{list(natural = NULL, always = 1, never = 0)}. \code{ref_int = 0} or
#' \code{"natural"} (the default) compares against the natural course: the
#' \code{NULL} element if there is one, otherwise a \code{natural} element is
#' added. Because the natural course draws the exposure from its fitted
#' model, \code{models} then needs an \code{"exposure"} model. An integer
#' position or a name selects one of your interventions instead, and no
#' natural course is added.
#'
#' \strong{Rank-deficient models.} Terms that cannot be estimated (a collinear
#' term, or one that is constant in the rows used for fitting, such as
#' \code{time} in an outcome recorded only at the last time point) are dropped
#' from the simulation, as \code{\link[stats]{predict}} would drop them, and
#' named in the warning summary printed on exit.
#'
#' The results rest on the correct temporal ordering and specification of the
#' models and on the usual g-formula assumptions (consistency, positivity, no
#' unmeasured confounding), which the package cannot check.
#'
#' @param data A \code{data.frame} in long format.
#' @param id_var Character. Name of the subject identifier.
#' @param base_vars Character vector of time-fixed baseline covariates (may be
#'   empty), with no missing values. Only these and \code{id_var} are carried
#'   into the simulated cohort; every other variable a model uses, including
#'   lags, must be created by the recode hooks.
#' @param exposure Character. Name of the exposure to intervene on.
#' @param time_var Character. Name of the numeric time variable. Each distinct
#'   value is one simulated step, in ascending order.
#' @param models List of \code{\link{spec_model}} objects, in the order the
#'   variables are generated.
#' @param intervention Named list of interventions. Each element is
#'   \code{NULL} (the natural course: exposure drawn from its fitted model), a
#'   0/1 value or vector with one value per distinct time point (a static
#'   intervention), or a \code{\link{dyn_int}} rule, e.g.
#'   \code{list(natural = NULL, treat_if_high = dyn_int(as.numeric(L1 > 0)))}.
#'   \code{NULL} (the default) runs the natural course only.
#' @param ref_int Reference for the contrasts: \code{0} or \code{"natural"}
#'   (default) for the natural course, or the position or name of an element
#'   of \code{intervention}. See Details.
#' @param init_recode \code{\link{recodes}} applied at the first time step,
#'   before the models are evaluated, e.g. to set lags to their baseline value.
#' @param in_recode \code{\link{recodes}} applied at the start of each later
#'   time step, before the models are evaluated, e.g. to update lags.
#' @param out_recode \code{\link{recodes}} applied at the end of each time
#'   step, the first included, e.g. to advance cumulative counts or carry an
#'   absorbing state forward.
#' @param return_fitted Logical. Return the full fitted model objects
#'   (default \code{FALSE}: calls and coefficients only).
#' @param mc_sample Number of subjects in the simulated Monte Carlo cohort.
#'   The default, \code{NULL}, uses 50 times the number of subjects in
#'   \code{data} and reports the value unless \code{quiet = TRUE}. A larger
#'   value reduces Monte Carlo error, which falls as
#'   \code{1/sqrt(mc_sample)}; it does not change the estimand.
#' @param return_data Logical. Return the simulated data (default
#'   \code{FALSE}; can be large).
#' @param R Number of bootstrap replicates (default \code{500}); \code{R = 1}
#'   skips the bootstrap. Replicates run through
#'   \code{future.apply::future_lapply}: sequentially unless a parallel plan
#'   is set with \code{future::plan()}, e.g. \code{plan(multisession)}.
#' @param quiet Logical. Suppress progress messages (default \code{FALSE}).
#' @param seed Integer random seed (default \code{12345}) for the Monte Carlo
#'   simulation and the bootstrap replicates; the global RNG state is restored
#'   on exit. \code{NULL} disables seeding, so repeated calls differ.
#'
#' @return
#' An object of class \code{"gformula"}, printed by
#' \code{\link{print.gformula}}, with components:
#' \itemize{
#'   \item \code{effect_size}: mean outcome (risk, for survival) under each
#'     intervention (\code{Intervention}, \code{Est}).
#'   \item \code{estimate}: with two or more interventions, the risk
#'     difference and risk ratio of each against \code{ref_int}
#'     (\code{Intervention}, \code{Risk_type}, \code{Estimate}).
#'   \item With \code{R > 1}, both tables gain \code{Sd} and percentile
#'     (\code{perct_lcl}, \code{perct_ucl}) and normal-approximation
#'     (\code{norm_lcl}, \code{norm_ucl}) limits, and \code{boot_estimates}
#'     holds the per-replicate estimates (\code{$interventions},
#'     \code{$contrasts}).
#'   \item \code{sim_data}: with \code{return_data = TRUE}, the simulated data
#'     as an end-of-follow-up snapshot: one row per Monte Carlo subject per
#'     intervention, each variable at its last simulated time step, with the
#'     accumulated \code{Pred_Y}. Earlier time steps are not kept.
#'   \item \code{fitted_models}: the fitted models, as full objects when
#'     \code{return_fitted = TRUE}, otherwise their calls and coefficients.
#'   \item \code{data_summary}: numbers of subjects, rows and time points in
#'     \code{data}.
#'   \item \code{observed}: a nonparametric benchmark printed beside the
#'     estimates: the observed mean outcome at the last time point, or the
#'     product-limit cumulative incidence for a survival outcome.
#'   \item \code{call}, \code{all.args}: the matched call and the evaluated
#'     arguments.
#' }
#'
#' @references
#' Robins, J. M. (1986). A new approach to causal inference in mortality studies with a sustained
#' exposure period—application to control of the healthy worker survivor effect. \emph{Mathematical Modelling}, 7(9–12), 1393–1512.
#' \doi{10.1016/0270-0255(86)90088-6}
#'
#' Keil, A. P., Edwards, J. K., Richardson, D. B., Naimi, A. I., & Cole, S. R. (2014).
#' The parametric g-formula for time-to-event data: intuition and a worked example.
#' \emph{Epidemiology}, 25(6), 889–897. \doi{10.1097/EDE.0000000000000160}
#'
#' @seealso \code{\link{spec_model}}, \code{\link{recodes}},
#'   \code{\link{dyn_int}}, \code{\link{mediation}} for direct and indirect
#'   effects, and \code{vignette("causalMed-03-gformula")}.
#'
#' @import data.table
#' @importFrom stats qnorm
#' @importFrom utils stack
#'
#' @examples
#' data(nonsurvivaldata)
#'
#' # Models in the order the variables are generated: A -> L1 -> L2 -> Y
#' models <- list(
#'   spec_model(A ~ V + lag1_A + lag1_L1 + lag1_L2 + time,
#'              var_type = "binary", mod_type = "exposure"),
#'   spec_model(L1 ~ V + A + lag1_L1 + time,
#'              var_type = "normal", mod_type = "covariate"),
#'   spec_model(L2 ~ V + A + lag1_L2 + time,
#'              var_type = "binary", mod_type = "covariate"),
#'   spec_model(Y_bin ~ V + A + L1 + L2,
#'              var_type = "binary", mod_type = "outcome")
#' )
#'
#' fit <- gformula(
#'   data = nonsurvivaldata, id_var = "id", time_var = "time",
#'   base_vars = "V", exposure = "A", models = models,
#'   intervention = list(natural = NULL, always = 1, never = 0),
#'   init_recode = recodes(lag1_A = 0, lag1_L1 = 0, lag1_L2 = 0),
#'   in_recode   = recodes(lag1_A = A, lag1_L1 = L1, lag1_L2 = L2),
#'   mc_sample = 2000,
#'   R = 1,          # use R > 1 (e.g. 500) for bootstrap confidence intervals
#'   quiet = TRUE
#' )
#' fit
#'
#' @export

gformula <- function(data,
                     id_var,
                     base_vars,
                     exposure,
                     time_var,
                     models,
                     intervention = NULL,
                     ref_int = 0,
                     init_recode = NULL,
                     in_recode = NULL,
                     out_recode = NULL,
                     return_fitted = FALSE,
                     mc_sample = NULL,
                     return_data = FALSE,
                     R = 500,
                     quiet = FALSE,
                     seed = 12345) {

  tpcall <- match.call()
  all.args <- mget(names(formals()),sys.frame(sys.nframe()))

  # Initilise warning
  init_warn()

  # Check for error
  check_error(data, id_var, base_vars, exposure, time_var, models)

  if (is.null(mc_sample)) {
    mc_sample <- 50L * data.table::uniqueN(data[[id_var]])
    if (!quiet) {
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

  # The simulation grid: the distinct observed time values, sorted. Kept as
  # `time_seq` so get_args_for() hands it to .run_interventions() and
  # bootstrap_helper(), fixing the grid for every pass.
  time_seq <- time_grid(data, time_var)
  time_len <- length(time_seq)

  if (!is.null(intervention)) {
    check_intervention(models, intervention, ref_int, time_len)

    # Resolve `ref_int` to the NAME of an element of `intervention`, adding the
    # natural course when the requested reference does not exist yet.
    # Downstream, risk_estimate_point()/risk_estimate_boot() drop the reference
    # with setdiff(names(intervention), ref_int): a numeric index never matches
    # a name (spurious self-contrast row), and a name that is absent silently
    # yields a zero-row contrast table. So every path below must end with
    # `ref_int` naming an element that really exists.
    if (is.numeric(ref_int) && ref_int >= 1) {
      # Positional reference. Safe to index the user's list directly: the
      # `natural` prepend below cannot have happened yet.
      ref_int <- names(intervention)[ref_int]
    } else if (identical(as.character(ref_int), "0") ||
               (identical(ref_int, "natural") &&
                !("natural" %in% names(intervention)))) {
      # "The natural course", asked for either positionally (ref_int = 0) or by
      # name. Both mean the NULL element, whatever the user called it -- so
      # resolve to that element's name and only fall through to creating one
      # when there is none. Matching on the NAME alone would prepend a SECOND
      # natural course next to an existing `list(nat = NULL, ...)`, breaking
      # check_intervention()'s one-NULL-element rule and emitting a nonsense
      # `nat - natural` contrast. The `!("natural" %in% names(...))` guard keeps
      # a non-NULL intervention the user chose to call "natural" as the
      # reference, and stops the prepend from colliding with its name.
      interv_value <- vapply(intervention, is.null, logical(1))
      ref_int <- if (any(interv_value)) names(which(interv_value))[1] else NA_character_
    }

    # Add the natural course when it was asked for but not supplied -- ref_int
    # = 0 or "natural" with no NULL element anywhere. The "natural" case used
    # to sail past every check and produce an empty contrast table.
    if (is.na(ref_int)) {
      intervention <- c(list(natural = NULL), intervention)
      ref_int <- "natural"
    }
    # Record the resolved reference so print() shows the name actually used
    # (e.g. "natural") rather than the raw 0 placeholder.
    all.args$ref_int <- ref_int
  }

  if (is.null(intervention)) {
    # Natural course only. Name it "natural" so the output label matches the
    # documentation and the label used whenever an explicit intervention list
    # is given.
    intervention <- list(natural = NULL)
    all.args$ref_int <- "natural"
  }

  # `intervention` is now a resolved named list in both paths, so one check
  # covers the explicit list and the intervention = NULL shorthand.
  check_natural_course(models, intervention, exposure)

  # Get the position of the outcome
  out_flag <- sapply(models, function(mods) mods$mod_type %in% c("outcome", "survival"))
  out_flag <- which(out_flag)
  outcome_var <- all.vars(formula(models[[out_flag]]$call)[[2]])

  # Data summary and observed nonparametric benchmark (displayed by print
  # as an informal model check against the simulated intervention means).
  is_survival  <- any(sapply(models, function(mods)
    mods$mod_type %in% c("survival", "censor")))
  data_summary <- summarize_input_data(data, id_var, time_seq)
  observed     <- observed_benchmark(data, outcome_var, time_var, is_survival)

  # Run original estimate.
  # NOTE: get_args_for() DROPS NULL-valued arguments, so `seed = NULL` would
  # never reach .run_interventions() and its default would silently take over.
  # Force the element through explicitly; `arg_est["seed"] <- list(seed)` keeps
  # a NULL element (plain `$seed <- NULL` would delete it again).
  arg_est <- get_args_for(.run_interventions)
  arg_est["seed"] <- list(seed)
  arg_est$return_fitted <- TRUE
  est_ori <- do.call(.run_interventions, arg_est)

  # Mean value of the outcome at each time point by intervention
  if(return_data){
    est_out <- data.table::rbindlist(est_ori$gform.data,
                                     idcol = "Intervention",
                                     use.names = TRUE)
    est_out <- est_out[, list(Est = sum(Pred_Y) / length(Pred_Y)),
                       by = c("Intervention")]
  }else{
    est_out <- data.table::as.data.table(utils::stack(est_ori$gform.data))
    colnames(est_out) <- c("Est", "Intervention")
  }
  setcolorder(est_out, c("Intervention", "Est"))

  # Per-replicate bootstrap estimates, retained on the returned object.
  # These are scalar summaries only (their size is independent of the input
  # data), so they are always kept when R > 1.
  boot_interventions <- NULL
  boot_contrasts     <- NULL

  # Run bootstrap
  if (R > 1) {
    arg_pools <- get_args_for(bootstrap_helper)
    arg_pools$progress_bar <- !quiet
    arg_pools$future_seed <- boot_seed
    pools <- do.call(bootstrap_helper, arg_pools)

    # Get the mean of bootstrap results.
    # pools_list is kept (before rbindlist) so it can be reused for contrasts
    # below, avoiding a redundant second lapply over all R replicates.
    pools_list <- lapply(pools, function(bt) {
      out <- utils::stack(bt$gform.data)
      colnames(out) <- c("Est", "Intervention")
      return(out)
    })
    boot_interventions <- data.table::rbindlist(pools_list, use.names = TRUE,
                                                idcol = "replicate")
    data.table::setcolorder(boot_interventions,
                            c("replicate", "Intervention", "Est"))
    pools_res <- data.table::rbindlist(pools_list, use.names = TRUE)

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
  }

  # Calculate the difference and ratio
  if (length(intervention) > 1) {
    risk_est <- risk_estimate_point(est_ori$gform.data,
                               ref_int = ref_int,
                               intervention = intervention,
                               return_data = return_data)

    if (R > 1) {
      res_pools <- lapply(pools_list,
                          risk_estimate_boot,
                          ref_int = ref_int,
                          intervention = intervention,
                          return_data = return_data)
      res_pools <- data.table::rbindlist(res_pools, idcol = "replicate")
      boot_contrasts <- data.table::copy(res_pools)
      # Calculate Sd and percentile confidence interval
      res_pools <- res_pools[, .(
        Sd = sd(Estimate, na.rm = TRUE),
        perct_lcl = quantile(Estimate, 0.025, na.rm = TRUE),
        perct_ucl = quantile(Estimate, 0.975, na.rm = TRUE)
      ),
      by = c("Intervention", "Risk_type")
      ]

      # Merge all and calculate the normal confidence interval
      risk_est <- merge(risk_est, res_pools, by = c("Intervention", "Risk_type"), sort = FALSE)
      risk_est <- risk_est[, `:=`(
        norm_lcl = Estimate - stats::qnorm(0.975) * Sd,
        norm_ucl = Estimate + stats::qnorm(0.975) * Sd
      )]
    }
  } else {
    risk_est <- NULL
  }

  # Extract fitted model information (attribute-wrapped so print/summary can
  # label each model; shared with mediation()).
  fitted_mods <- wrap_fitted_models(est_ori$fitted.models, return_fitted)

  # Return data
  if(return_data){
    dat_out <- data.table::rbindlist(est_ori$gform.data, idcol = "Intervention",use.names=T)
  }else{
    dat_out <- NULL
  }

  emit_warnings()

  boot_estimates <- if (R > 1) {
    list(interventions = boot_interventions, contrasts = boot_contrasts)
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
            observed = observed
          )
  class(y) <- c("gformula", class(y))
  return(y)
}
