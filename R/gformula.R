
#' G-formula Analysis
#'
#' @description Mediation analysis for time varying mediator, estimation based
#'  on g-formula. Output contains total effect, #' natural direct effect and natural
#'   indirect effect for mediation or regular g-formula. data.frame will be returned.
#'
#' @note Not that final outcome must be the same for in all rows per subject.
#'   If the dataset is  survival settings, records after the interested outcome
#'    must be deleted. The function it self do some data manipulation internally.
#'    Please prepare the data as longitudinal format.
#'
#' @param data Data set to be sued
#'
#' @param id_var ID variable per subject.
#'
#' @param base_vars A vector of time fixed baseline variables.
#'
#' @param exposure Intervention/Exposure variable
#'
#' @param time_var Time variable.
#'
#' @param models A list of models for the G-formula, including exposure model,
#'  covariate model (if any), mediator model (if any), outcome model or
#'  censoring model (if any). See details in \code{\link{spec_model}}.
#'  The order appeared in the list should reflect the temporal ordering of the
#'  variables, in another way data generation process. The model will be evaluated
#' in this process.
#'
#' @param intervention A named list with a value of intervention on exposure.
#' if kept as NULL (default), the natural intervention course will be calculated.
#'  eg: list(natural = NULL, always = c(1, 1, 1), never = c(0, 0, 0)). Note that
#' the function will add the natural effect by default if the natural effect is
#' not define here, and the reference intervention in \code{ref_int} will be set
#' to the natural effect.
#'
#' @param ref_int Integer or character string denoting the intervention to be used as
#' the reference for calculating the risk ratio and risk difference. 0 or \code{"natural"}
#'  denotes the natural course, while subsequent integers denote user-specified
#' interventions in the order that they are named in \code{intervention}. Or if the character
#' is the name of the intervention. The default is 0.
#'
#' @param init_recode optional, recoding of variables done at the
#' beginning of the Monte Carlo loop. Needed for operations initialize baseline variables.
#' This is executed at beginning of the Monte Carlo g-formula, executed only once at time 0.
#'
#' @param in_recode optional, On the fly recoding of variables done before the Monte
#'  Carlo loop starts. Needed to do any kind of functional forms for entry times.
#'   This is executed at each start of the Monte Carlo g-formula time steps
#'
#' @param out_recode optional, On the fly recoding of variables done at the
#' end of the Monte Carlo loop. Needed for operations like counting the number of
#' days with a treatment or creating lagged variables. This is executed at each end
#'  of the Monte Carlo g-formula time steps.
#'
#' @param return_fitted Logical scalar indicating whether to return the fitted model. Default
#' is \code{FALSE}, only the summary of coefficients will be returned to reduce the memory usage.
#'
#' @param mc_sample Integer, sample size of Monte Carlo simulation. Default is 100 times of
#' original sample.
#'
#' @param return_data Logical scalar indicating whether to return the Monte Carlo simulated data set,
#'  default is \code{FALSE}. In the case of large \code{mc_sample} and long follow-up, the simulated
#' data will be very large causing problems no enough memory to allocate the objects.
#'
#' @param R The number of bootstrap replicates, default is 500. If the R is larger than 1,
#' \code{\link[future.apply]{future_lapply}} is used for the parallel computation. This will
#'  run in sequential if no parallel planed, see \code{\link[future]{plan}} for more details.
#' If the parallel plan was set to \code{plan(multisession)} in Windows or \code{plan(multicore)}
#' in other system, multiple session/core is used for bootstrap calculation.
#'
#' @param quiet if \code{TRUE} then the progress bar will be suppressed.
#'
#' @references
#' Robins, J. (1986). A new approach to causal inference in mortality studies with a sustained exposure period—application to control of the healthy worker survivor effect. Mathematical modelling, 7(9-12), 1393-1512.
#' 
#' Keil, A. P., Edwards, J. K., Richardson, D. B., Naimi, A. I., & Cole, S. R. (2014). The parametric g-formula for time-to-event data: intuition and a worked example. Epidemiology, 25(6), 889-897.
#'
#'
#' @import data.table
#'
#' @export
#'
#'

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

  set.seed(seed)

  # Get time length
  time_len <- length(unique(na.omit(data[[time_var]])))

  if (!is.null(intervention)) {
    check_intervention(models, intervention, ref_int, time_len)

    # If any intervention is set to NULL, but reference not defined.
    if (ref_int == 0 & length(intervention) >= 1) {
      interv_value <- sapply(intervention, is.null)
      if (any(interv_value)) {
        ref_int <- which(interv_value)
      } else {
        intervention <- c(list(natural = NULL), intervention)
        ref_int <- 1
      }
    }
  }

  data.table::setDT(data)

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
    pools <- do.call(bootstrap_helper, arg_pools)

    # Get the mean of bootstrap results
    pools_res <- lapply(pools, function(bt) {
      out <- utils::stack(bt)
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
    est_out <- merge(est_out, pools_res, by = c("Intervention"))
    est_out <- est_out[, `:=`(
      norm_lcl = Est - stats::qnorm(0.975) * Sd,
      norm_ucl = Est + stats::qnorm(0.975) * Sd
    )]
  }

  # Risk difference and risk ratio calculation function
  risk_calc <- function(data_list, ref_int, return_data) {
    ref_nam <- names(intervention)[ref_int]
    ref_dat <- data_list[[ref_nam]]
    vs_nam <- setdiff(names(intervention), ref_nam)
    out <- sapply(vs_nam, function(x) {
      tmp_dt <- data_list[[x]]
      if(return_data){
        ref_mean <- sum(unlist(ref_dat)) / length(ref_dat)
        vs_mean <- sum(unlist(tmp_dt)) / length(tmp_dt)
      }else{
        ref_mean <- ref_dat
        vs_mean <- tmp_dt
      }
      data.table(
        Risk_type = c("Risk difference", "Risk ratio"),
        Estimate = c(vs_mean - ref_mean, vs_mean / ref_mean)
      )
    }, simplify = FALSE)
    data.table::rbindlist(out, idcol = "Intervention", use.names = TRUE)
  }


  # Calculate the difference and ratio
  if (length(intervention) > 1) {
    risk_est <- risk_calc(est_ori$gform.data, ref_int = ref_int, return_data = return_data)

    if (R > 1) {
      res_pools <- lapply(pools, risk_calc, ref_int = ref_int, return_data = return_data)
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
      risk_est <- merge(risk_est, res_pools, by = c("Intervention", "Risk_type"))
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
