
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
#' @param id.var ID variable per subject.
#'
#' @param base.vars A vector of time fixed baseline variables.
#'
#' @param exposure Intervention/Exposure variable
#'
#' @param time.var Time variable.
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
#'  eg: list(natural = NULL, always = c(1, 1, 1), never = c(0, 0, 0))
#'
#' @param ref_int Integer or character string denoting the intervention to be used as
#' the reference for calculating the risk ratio and risk difference. 0 or \code{"natural"}
#'  denotes the natural course, while subsequent integers denote user-specified
#' interventions in the order that they are named in \code{intervention}. Or if the character
#' is the name of the intervention. The default is 0.
#'
#' @param init.recode optional, recoding of variables done at the
#' beginning of the Monte Carlo loop. Needed for operations initialize baseline variables.
#' This is executed at beginning of the Monte Carlo g-formula, executed only once at time 0.
#'
#' @param in.recode optional, On the fly recoding of variables done before the Monte
#'  Carlo loop starts. Needed to do any kind of functional forms for entry times.
#'   This is executed at each start of the Monte Carlo g-formula time steps
#'
#' @param out.recode optional, On the fly recoding of variables done at the
#' end of the Monte Carlo loop. Needed for operations like counting the number of
#' days with a treatment or creating lagged variables. This is executed at each end
#'  of the Monte Carlo g-formula time steps.
#'
#' @param is.survival Is the data survival data, default is FALSE.
#'
#' @param mc.sample Sample size of Monte Carlo simulation.
#'
#' @param R The number of bootstrap replicates, default is 500.
#'
#' @param ncores integer: number of processes to be used in parallel operation: typically one would chose this to the number of
#'  available CPUs. Parallel computation for bootstrap will be applied if \code{ncores} larger than 1.
#'
#' @references
#' Lin, S. H., Young, J. G., Logan, R., & VanderWeele, T. J. (2017). Mediation analysis for a survival outcome with time-varying exposures, mediators, and confounders. \emph{Statistics in medicine}, 36(26), 4153-4166. DOI:10.1002/sim.7426
#' Zheng, W., & van der Laan, M. (2017). Longitudinal mediation analysis with time-varying mediators and exposures, with application to survival outcomes. \emph{Journal of causal inference}, 5(2). DOI:10.1515/jci-2016-0006
#'
#'
#' @import data.table
#'
#' @export
#'
#'

Gformula <- function(data,
                     id.var,
                     base.vars,
                     exposure,
                     time.var,
                     models,
                     intervention = NULL,
                     ref_int      = 0,
                     init.recode  = NULL,
                     in.recode    = NULL,
                     out.recode   = NULL,
                     mc.sample    = 10000,
                     R            = 500,
                     ncores       = 1L) {

  tpcall <- match.call()

  orig.seed <- .Random.seed
  on.exit(.Random.seed <<- orig.seed)

  # Get time length
  time_seq <- unique(na.omit(data[[time.var]]))
  time_len <- length(time_seq)

  # Check for error
  check_error(data, id.var, base.vars, exposure, time.var, models)

  if(!is.null(intervention)){
    check_intervention(models, intervention, ref_int, time_len)

    # If any intervention is set to NULL, but reference not defined.
    if(ref_int == 0 & length(intervention) > 1){
      interv_value <- sapply(intervention, is.null)
      if(any(interv_value))
        ref_int <- which(interv_value)
    }else{
      intervention <- c(list(natural = NULL), intervention)
      ref_int <- 1
    }
  }

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
  est_ori <- do.call(.gformula, arg_est)

  # Run bootstrap
  arg_pools <- get_args_for(bootstrap_helper)
  pools <- do.call(bootstrap_helper, arg_pools)

  # Mean value of the outcome at each time point by intervention
  est_out <- data.table::rbindlist(est_ori$gform.data, idcol = "Intervention")
  est_out <- est_out[, list(Est = mean(Pred_Y)), by = c("Intervention")]

  # Get the mean of bootstrap results
  pools_res <- lapply(pools, function(bt){
    out <- sapply(bt, function(x){
      x[,list(Est = mean(Pred_Y))]
    }, simplify = FALSE)
    data.table::rbindlist(out, idcol = "Intervention")
  })
  pools_res <- data.table::rbindlist(pools_res)

  # Calculate Sd and percentile confidence interval
  pools_res <- pools_res[, .(Sd = sd(Est), perct_lcl = quantile(Est, 0.025),
                             perct_ucl = quantile(Est, 0.975)) ,
                         by = c("Intervention")]

  # Merge all and calculate the normal confidence interval
  est_out <- merge(est_out, pools_res, by = c("Intervention"))
  est_out <- est_out[, `:=` (norm_lcl = Est - stats::qnorm(0.975)*Sd,
                             norm_ucl = Est + stats::qnorm(0.975)*Sd)]

  # Risk difference and risk ratio calculation function
  risk_calc <- function(data_list, ref_int){
    # Get reference
    ref_nam <- names(intervention)[ref_int]
    ref_dat <- data_list[[ref_nam]]
    vs_nam  <- setdiff(names(intervention), ref_nam)

    # Loop through each intervention name
    out <- sapply(vs_nam, function(x){
      tmp_dt <- data_list[[x]]

      # Calculate difference and ratio of mean at each time
      ref_mean <- mean(ref_dat[["Pred_Y"]])
      vs_mean <- mean(tmp_dt[["Pred_Y"]])
      data.table(Risk_type = c("Risk difference", "Risk ratio"),
                 Estimate = c( vs_mean - ref_mean, vs_mean/ref_mean))
    }, simplify = FALSE)
    data.table::rbindlist(out, idcol = "Intervention")
  }

  # Calculate the difference and ratio
  if(length(intervention) >= 1){
    risk_ori <- risk_calc(est_ori$gform.data, ref_int)
    res_pools <- lapply(pools, risk_calc, ref_int)
    res_pools <- data.table::rbindlist(res_pools)
    # Calculate Sd and percentile confidence interval
    res_pools <- res_pools[, .(Sd = sd(Estimate, na.rm = TRUE),
                               perct_lcl = quantile(Estimate, 0.025, na.rm = TRUE),
                               perct_ucl = quantile(Estimate, 0.975, na.rm = TRUE)) ,
                           by = c("Intervention", "Risk_type")]

    # Merge all and calculate the normal confidence interval
    risk_est <- merge(risk_ori, res_pools, by = c("Intervention", "Risk_type"))
    risk_est <- risk_est[, `:=` (norm_lcl = Estimate - stats::qnorm(0.975)*Sd,
                                 norm_ucl = Estimate + stats::qnorm(0.975)*Sd)]

  }

  return(list(
    call = tpcall,
    estimate = est_out,
    risk_est = risk_est,
    gform.data = data.table::rbindlist(est_ori$gform.data, idcol = "Intervention"),
    fitted.models = est_ori$fitted.models
  ))
}




