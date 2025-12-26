
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
#' If at least one element is \code{NULL} and \code{ref_int == 0}, that element is used
#' as the reference ("natural"). Otherwise, a \code{natural} arm is added automatically,
#' and \code{ref_int} is set to \code{"natural"}.
#'
#' **Model specification**
#' Each element of \code{models} is typically created by \code{\link{spec_model}}
#'  and must include: (i) the model formula/call, (ii) a \code{mod_type} indicating its role
#' (\code{"exposure"}, \code{"covariate"}, \code{"outcome"}, \code{"survival"},
#' or \code{"censoring"}), and (iii) a \code{var_type} specifying the
#' variable type used for simulation/prediction (\code{"binomial"}, \code{"normal"},
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
#'   element is either \code{NULL} (the natural course) or a numeric/logical
#'   vector whose length equals the number of unique time points in
#'   \code{time_var}. For example,
#'   \code{list(natural = NULL, always = c(1, 1, 1), never = c(0, 0, 0))}.
#'   If \code{intervention} is \code{NULL}, only the natural course is
#'   evaluated. If a \code{natural} arm is not provided, it is added
#'   automatically and \code{ref_int} is set to \code{"natural"}.
#' @param ref_int Reference intervention for contrasts. Either an integer index
#'   (\code{0} = natural course; \code{1}, \code{2}, … = elements of
#'   \code{intervention}) or a character name matching an element (e.g.,
#'   \code{"always"}). If no \code{natural} arm is provided, it is added and
#'   \code{ref_int} is set to \code{"natural"}. Default: \code{0}.
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
#'   \item \code{sim_data}: if \code{return_data = TRUE}, the stacked simulated
#'         Monte Carlo dataset across interventions (can be large).
#'   \item \code{fitted_models}: a named list of fitted models.
#'         If \code{return_fitted = TRUE}, returns full model objects plus attributes
#'         (\code{recodes}, \code{subset}, \code{var_type}, \code{mod_type});
#'         otherwise, a compact list with \code{call} and \code{coeff}.
#' }
#'
#' @note
#' Final outcome should be consistently defined at the terminal time for each subject.
#' For survival-type applications, remove rows after the event of interest (or include
#' and model censoring appropriately). The function may record warnings internally and
#' print them on exit. Results depend on correct temporal ordering, model specification,
#' positivity, and no unmeasured confounding assumptions customary for g-formula.
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
#' @export
#'
#' @examples
#' \dontrun{
#' ## Toy longitudinal data (long format)
#' data(nonsurvivaldata)
#'
#' ## Specify models in temporal order, e.g.:
#' mod1<-spec_model(A ~ A_lag1  +V,var_type= "binomial",mod_type = "exposure")
#' mod2<-spec_model(L ~ A + L_lag1  +V,var_type  = "normal",mod_type = "covariate")
#' mod3<-spec_model(Y ~  A + L  + V,var_type = "binomial",mod_type = "outcome")
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
#'   init_recode = c(A_lag1 = 0, L_lag1 = 0),
#'   in_recode = c(A_lag1 = A, L_lag1 = L),
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

    # If any intervention is set to NULL, but reference not defined.
    if (ref_int == 0 & length(intervention) >= 1) {
      interv_value <- sapply(intervention, is.null)
      if (any(interv_value)) {
        ref_int <- names(which(interv_value))
      } else {
        intervention <- c(list(natural = NULL), intervention)
        ref_int <- "natural"
      }
    }
  }

  data <- as.data.table(data)

  if (is.null(intervention)) {
    intervention <- list(intervention = NULL)
  }

  # Get the position of the outcome
  out_flag <- sapply(models, function(mods) mods$mod_type %in% c("outcome", "survival"))
  out_flag <- which(out_flag)
  outcome_var <- all.vars(formula(models[[out_flag]]$call)[[2]])

  # Run original estimate
  arg_est <- get_args_for(.gformula)
  arg_est$return_fitted <- TRUE
  arg_est$progress_bar <- substitute(quiet, env = parent.frame())
  est_ori <- do.call(.gformula, arg_est)

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

  # Run bootstrap
  if (R > 1) {
    arg_pools <- get_args_for(bootstrap_helper)
    arg_pools$progress_bar <- substitute(quiet, env = parent.frame())
    arg_pools$future_seed <- boot_seed
    pools <- do.call(bootstrap_helper, arg_pools)

    # Get the mean of bootstrap results
    pools_res <- lapply(pools, function(bt) {
      out <- utils::stack(bt$gform.data)
      colnames(out) <- c("Est", "Intervention")
      return(out)
    })
    pools_res <- data.table::rbindlist(pools_res, use.names = TRUE)

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
    risk_est <- risk_estimate1(est_ori$gform.data,
                               ref_int = ref_int,
                               intervention = intervention,
                               return_data = return_data)

    if (R > 1) {
      pools2 <- lapply(pools, function(bt) {
        out <- utils::stack(bt$gform.data)
        colnames(out) <- c("Est", "Intervention")
        return(out)
      })
      res_pools <- lapply(pools2,
                          risk_estimate2,
                          ref_int = ref_int,
                          intervention = intervention,
                          return_data = return_data)
      res_pools <- data.table::rbindlist(res_pools)
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

  # Extract fitted model information
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
    dat_out <- data.table::rbindlist(est_ori$gform.data, idcol = "Intervention",use.names=T)
  }else{
    dat_out <- NULL
  }

  if(length(causalmed_env$warning)>0){
    message(paste(causalmed_env$warning, collapse = "\n=============\n"),
            domain = "causalMed")
  }

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
