
#' Parametric G-formula for Time-varying Intervention Analyses
#'
#' @description
#' Implements the **parametric g-formula** to estimate counterfactual mean outcomes
#' under one or more user-specified exposure interventions in longitudinal data.
#' Supports settings with time-fixed baselines, time-varying exposure/mediator/confounders,
#' optional survival/terminal outcomes, and nonparametric uncertainty via bootstrap.
#'
#' @details
#' The function evaluates a sequence of user-specified models (see \code{models})
#' in **temporal order** to simulate counterfactual trajectories via Monte Carlo,
#' producing: (i) intervention-specific mean outcomes, and (ii) contrasts vs. a
#' reference intervention (risk ratio/difference). When \code{R > 1}, percentile
#' and normal-approximation confidence intervals are computed from bootstrap resamples.
#'
#' **Data requirements**
#' \itemize{
#'   \item Long format: one row per subject per time point.
#'   \item \code{id_var}: unique subject identifier; \code{time_var}: ordered time index.
#'   \item Final outcome must be well-defined at the last relevant time for each subject.
#'         For survival-like settings, rows after the event of interest should be removed
#'         (or a censoring model should be included and handled in \code{models}).
#' }
#'
#' **Interventions**
#' Provide a named list \code{intervention} with exposure values per time
#' (e.g., \code{list(natural = NULL, always = c(1,1,1), never = c(0,0,0))}).
#' If \code{intervention} is \code{NULL}, the function evaluates the natural course only.
#'
#' The reference is resolved before simulation, so that \code{ref_int} always
#' names an intervention that exists:
#' \itemize{
#'   \item \code{ref_int = 0} and \code{ref_int = "natural"} both ask for the
#'         natural course. The \code{NULL} element of \code{intervention} is
#'         used if one was supplied — whatever it is named — and otherwise a
#'         \code{natural} element is added to the list.
#'   \item An integer position, or a name matching an element, uses that
#'         intervention as the reference; no natural course is added.
#' }
#' Because the default is \code{ref_int = 0}, a list holding only static or
#' dynamic interventions still gains a natural course. The natural course draws
#' the exposure from its fitted model, so \code{models} must then include a
#' \code{mod_type = "exposure"} model; supplying one is required unless
#' \code{ref_int} names one of the interventions you provided.
#'
#' **Model specification**
#' Each element of \code{models} is typically created by \code{\link{spec_model}}
#'  and must include: (i) the model formula/call, (ii) a \code{mod_type} indicating its role
#' (\code{"exposure"}, \code{"covariate"}, \code{"outcome"}, \code{"survival"},
#' or \code{"censor"}), and (iii) a \code{var_type} specifying the
#' variable type used for simulation/prediction (\code{"binary"}, \code{"normal"},
#' \code{"categorical"}, and \code{"custom"}). The list order must reflect the
#' data-generating process (temporal ordering). The outcome model is detected internally
#' and used for computing predicted outcomes (\code{Pred\_Y}) under each intervention.
#'
#' **Re-coding hooks**
#' \itemize{
#'   \item \code{init_recode}: executed once at time 0 before simulation (initialize baselines).
#'   \item \code{in_recode}: executed at the start of each time step (e.g., entry-time logic).
#'   \item \code{out_recode}: executed at the end of each time step (e.g., cumulative counts, lags).
#' }
#'
#' @param data A \code{data.frame} (long format): one row per \code{id_var} per \code{time_var}.
#' @param id_var Character scalar. Subject identifier column name.
#' @param base_vars Character vector of time-fixed baseline covariates (may be empty).
#' @param exposure Character scalar. Exposure variable to intervene on (must be in \code{data}).
#' @param time_var Character scalar. Time variable column name (ordered; integer/numeric).
#' @param models A list of model specifications evaluated in temporal order. The order
#'   appeared in the list should reflect the temporal ordering of the variables, in another
#'   way data generation process. See \code{\link{spec_model}} for a recommended constructor.
#' @param intervention A named list specifying exposure interventions. Each
#'   element is one of:
#'   \itemize{
#'     \item \code{NULL} — the natural course (exposure drawn from its fitted model).
#'     \item A numeric/logical scalar or vector (length 1 or equal to the number
#'           of time points) — a static intervention setting the exposure to that
#'           value at every (or each specific) time step.
#'     \item A \code{\link{dyn_int}} object — a dynamic (rule-based) intervention
#'           whose expression is evaluated inside the simulated dataset at each
#'           time step. Column names (including the exposure after its natural-course
#'           draw) are in scope, e.g.
#'           \code{list(natural = NULL, threshold = dyn_int(as.numeric(A > 0)))}.
#'   }
#'   If \code{intervention} is \code{NULL}, only the natural course is
#'   evaluated. A \code{natural} element is also added automatically when
#'   \code{ref_int} asks for the natural course and the list contains no
#'   \code{NULL} element; see \code{ref_int} and Details.
#' @param ref_int Reference intervention for contrasts. Either an integer index
#'   (\code{0} = natural course; \code{1}, \code{2}, … = elements of
#'   \code{intervention}) or a character name matching an element (e.g.,
#'   \code{"always"}). \code{0} and \code{"natural"} both resolve to the
#'   \code{NULL} element of \code{intervention} if there is one, and otherwise
#'   add a \code{natural} element — which requires an exposure model in
#'   \code{models}. See Details. Default: \code{0}.
#' @param init_recode Optional expression/function applied once at time 0 before the Monte Carlo loop
#'   (e.g., initializing baseline-derived variables). Should be defined with \code{\link{recodes}}. See Details.
#' @param in_recode Optional expression/function applied at the **start** of each time step
#'   (e.g., entry-time functional forms). Should be defined with \code{\link{recodes}}. See Details.
#' @param out_recode Optional expression/function applied at the **end** of each time step
#'   (e.g., create lags, cumulative counts). Should be defined with \code{\link{recodes}}. See Details.
#' @param return_fitted Logical. If \code{TRUE}, return full fitted model objects; otherwise,
#'   a light-weight summary (call and coefficients). Default \code{FALSE}.
#' @param mc_sample Integer. Monte Carlo sample size used for simulation.
#'   Default \code{nrow(data) * 100}.
#' @param return_data Logical. If \code{TRUE}, return the stacked simulated data
#'   (all interventions) including predicted outcomes; may be large. Default \code{FALSE}.
#' @param R Number of bootstrap replicates. If \code{R > 1},
#'   computation uses \code{future.apply::future_lapply} and runs sequentially
#'   unless a parallel plan is set; see \code{future::plan}. Use
#'   \code{plan(multisession)} on Windows and \code{plan(multicore)} on
#'   Unix-alikes to enable parallel bootstrap. Default \code{500}.
#' @param quiet Logical. If \code{TRUE}, suppress progress messages/bars. Default \code{FALSE}.
#' @param seed Integer random seed for reproducibility. Default \code{12345}.
#'   The seed fixes the RNG stream used for the Monte Carlo simulation
#'   \emph{and} the bootstrap replicates, and the caller's global RNG state is
#'   saved and restored on exit, so a seeded call does not disturb the ambient
#'   random stream. Pass \code{NULL} to disable seeding entirely: no
#'   \code{set.seed()} is called, the Monte Carlo draws consume and advance the
#'   ambient RNG stream, and repeated calls therefore give \emph{different}
#'   results (reproducible only via an outer \code{set.seed()}). Use
#'   \code{seed = NULL} inside simulation loops that manage their own seeds.
#'
#' @return
#' An object of class \code{"gformula"} with components:
#' \itemize{
#'   \item \code{call}: the matched call.
#'   \item \code{all.args}: a named list of evaluated arguments for reproducibility.
#'   \item \code{effect_size}: \code{data.table} with columns \code{Intervention} and
#'         \code{Est} (intervention-specific mean outcome). If \code{R > 1}, also includes
#'         \code{Sd}, percentile CIs (\code{perct_lcl}, \code{perct_ucl}) and normal CIs
#'         (\code{norm_lcl}, \code{norm_ucl}).
#'   \item \code{estimate}: if multiple interventions are provided, a \code{data.table}
#'         of contrasts vs. \code{ref_int} (columns typically include \code{Intervention},
#'         \code{Risk_type}, \code{Estimate}, and (if \code{R > 1}) CI columns).
#'   \item \code{sim_data}: if \code{return_data = TRUE}, the simulated Monte
#'         Carlo dataset stacked across interventions with an
#'         \code{Intervention} column (can be large). It is the
#'         \strong{end-of-follow-up snapshot} — one row per Monte Carlo subject
#'         per intervention, holding each variable at its final simulated time
#'         step alongside the accumulated \code{Pred_Y} — not a
#'         row-per-time-point panel: the simulation overwrites each variable in
#'         place as it steps through time.
#'   \item \code{fitted_models}: a named list of fitted models.
#'         If \code{return_fitted = TRUE}, returns full model objects plus attributes
#'         (\code{recodes}, \code{subset}, \code{var_type}, \code{mod_type});
#'         otherwise, a compact list with \code{call} and \code{coeff}.
#'   \item \code{boot_estimates}: when \code{R > 1}, a list of per-replicate
#'         bootstrap estimates: \code{$interventions} (columns
#'         \code{replicate}, \code{Intervention}, \code{Est}) and
#'         \code{$contrasts} (columns \code{replicate}, \code{Intervention},
#'         \code{Risk_type}, \code{Estimate}; \code{NULL} with a single
#'         intervention). These are scalar summaries whose size is
#'         independent of the input data. \code{NULL} when \code{R <= 1}.
#'   \item \code{data_summary}: list with the number of individuals
#'         (\code{n_id}), observations (\code{n_obs}), and time points
#'         (\code{n_times}, \code{t_min}, \code{t_max}) of the input data.
#'   \item \code{observed}: list with the observed nonparametric benchmark of
#'         the outcome (\code{value}, \code{label}): the mean outcome at the
#'         last time point, or the product-limit cumulative incidence for
#'         survival outcomes. Printed alongside the simulated means as an
#'         informal model check.
#' }
#'
#' @note
#' Final outcome should be consistently defined at the terminal time for each subject.
#' For survival-type applications, remove rows after the event of interest (or include
#' and model censoring appropriately). The function may record warnings internally and
#' print them on exit. Results depend on correct temporal ordering, model specification,
#' positivity, and no unmeasured confounding assumptions customary for g-formula.
#'
#' If a fitted model turns out to be rank deficient — a collinear term, or a
#' covariate that is constant among the rows actually used to fit it, such as
#' \code{time} in an outcome recorded only at the end of follow-up — the
#' unestimable terms are dropped from the simulation, exactly as
#' \code{\link[stats]{predict}} would, and are named in the warning summary
#' printed on exit. Treat that warning as a prompt to fix the formula.
#'
#' @references
#' Robins, J. M. (1986). A new approach to causal inference in mortality studies with a sustained
#' exposure period—application to control of the healthy worker survivor effect. \emph{Mathematical Modelling}, 7(9–12), 1393–1512.
#'
#' Keil, A. P., Edwards, J. K., Richardson, D. B., Naimi, A. I., & Cole, S. R. (2014).
#' The parametric g-formula for time-to-event data: intuition and a worked example.
#' \emph{Epidemiology}, 25(6), 889–897.
#'
#' @import data.table
#' @importFrom stats qnorm
#' @importFrom utils stack
#'
#' @examples
#' \dontrun{
#' ## Toy longitudinal data (long format)
#' data(nonsurvivaldata)
#'
#' ## Specify models in temporal order, e.g.:
#' mod1<-spec_model(A ~ A_lag1  +V,var_type= "binary",mod_type = "exposure")
#' mod2<-spec_model(L ~ A + L_lag1  +V,var_type  = "normal",mod_type = "covariate")
#' mod3<-spec_model(Y ~  A + L  + V,var_type = "binary",mod_type = "outcome")
#' models1<-list(mod1,mod2,mod3)
#'
#' ## Define interventions over T time points:
#' ints <- list(natural = NULL,
#'              always = 1
#'              )
#'
#' fit <- gformula(
#'   data = sim_data0,
#'   id_var = 'id',
#'   base_vars = "V",
#'   exposure = "A",
#'   time_var = "time",
#'   models = models1,
#'   intervention = ints,
#'   ref_int = 1,
#'   init_recode = recodes(A_lag1 = 0, L_lag1 = 0),
#'   in_recode = recodes(A_lag1 = A, L_lag1 = L),
#'   out_recode = NULL,
#'   return_fitted = TRUE,
#'   mc_sample = 100000,
#'   return_data = TRUE,
#'   R = 500,
#'   quiet = TRUE,
#'   seed = 250817
#' )
#'
#' print(fit)
#' summary(fit)
#'}
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
                     mc_sample = nrow(data)*100,
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

  # Get time length
  time_len <- length(unique(na.omit(data[[time_var]])))

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
  data_summary <- summarize_input_data(data, id_var, time_var)
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
