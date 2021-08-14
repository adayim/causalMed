
#' G-formula Analysis
#'
#' @description Mediation analysis for time varying mediator, estimation based
#'  on g-formula. Output contains total effect, #' natrual direct effect and natural
#'   indirect effect for mediation or regular g-formula. data.frame will be returned.
#'
#' @note Not that final outcome must be the same for in all rows per subject.
#'   If the dataset is  survival settings, records after the interested outcome
#'    must be deleted. The funciton it self do some data manupilation internally.
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
#' @param outcome Name of the outcome variable.
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
#' @param mediation_type Type of the mediation effect, natural effect (\code{"N"}) or interventional effect (\code{"I"}).
#'
#' @param verbose Print intervention information during calculation.
#'
#' @references
#' Lin, S. H., Young, J. G., Logan, R., & VanderWeele, T. J. (2017). Mediation analysis for a survival outcome with time-varying exposures, mediators, and confounders. \emph{Statistics in medicine}, 36(26), 4153-4166. DOI:10.1002/sim.7426
#' Zheng, W., & van der Laan, M. (2017). Longitudinal mediation analysis with time-varying mediators and exposures, with application to survival outcomes. \emph{Journal of causal inference}, 5(2). DOI:10.1515/jci-2016-0006
#'
#'
#' TODO: weights, time varying intervention
#'
#' @export
#'
#'

mediation <- function(data,
                      id.var,
                      base.vars,
                      exposure,
                      outcome,
                      time.var,
                      models,
                      init.recode = NULL,
                      in.recode = NULL,
                      out.recode = NULL,
                      mc.sample = 10000,
                      mediation_type = c("N", "I"),
                      verbose = TRUE) {
  tpcall <- match.call()
  mediation_type <- match.arg(mediation_type)

  # Setting seeds
  if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) runif(1)
  seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)

  # Check for error
  do.call(check_error, tpcall)

  if (is.null(intervention)) {
    intervention <- list(intervention = NULL)
  }

  # Run along the models
  fit_mods <- lapply(models, function(mods) {
    mods$call$data <- substitute(data)
    rsp_vars <- all.vars(formula(mods$call)[[2]])
    is_numeric <- is.numeric(data[[rsp_vars]])

    list(
      fitted = eval(mods$call),
      recodes = mods$recode,
      subset = mods$subset,
      var_type = mods$var_type,
      mod_type = mods$mod_type,
      custom_sim = mods$custom_sim,
      rsp_vars = rsp_vars,
      val_ran = ifelse(is_numeric, range(na.omit(data[[rsp_vars]])),
        unique(na.omit(data[[rsp_vars]]))
      )
    ) # Observed values range
  })

  # Baseline variables for Monte Carlo Sampling
  base_dat <- unique(data[, unique(c(id.var, base.vars))])
  df_mc <- base_dat[sample(1:nrow(base_dat), mc.sample, replace = TRUE), ]

  # Mediation analysis
  if (!all(med_flag == 0)) {
    intervention <- list(
      always = rep(1, time_len),
      never = rep(0, time_len),
      mediation = rep("mediation", time_len)
    )
  }

  res <- sapply(intervention, function(i) {
    monte_g(
      data = df_mc, 
      exposure = exposure,
      time.var = time.var,
      time.seq = unique(data[[time.var]]),
      models = fit_mods, intervention = i,
      in.recode = in.recode,
      out.recode = out.recode, init.recode = init.recode,
      verbose = verbose
    )
  }, simplify = FALSE)

  return(list(
    call = tpcall,
    fitted.models = fit_mods,
    gform.data = res
  ))
}
