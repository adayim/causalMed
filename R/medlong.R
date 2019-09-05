#' Mediation analysis for time-varying mediaiton, survival and non-survival
#'
#' @description Mediation analysis for time varying mediator, estimation based
#'  on g-formula and IPTW. Output contains total effect,
#' natrual direct effect and natural indirect effect. The confidence interval will
#' be calculated using 1000 bootstrap with normal approximation.
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
#'  censoring model (if any). See details in \code{\link[causalMed]{spec_model}}.
#'
#' @param intervention A vector or a value of intervention on exposure. if kept as NULL (default),
#' the natrual intervention cousre will be calculated.
#'
#' @param init.recode optional, recoding of variables done at the
#' begaining of the Monte Carlo loop. Needed for operations initalize baseline variables.
#' This is executed at begaining of the Monte Carlo g-formula, excuted only once.
#'
#' @param out.recode optional, On the fly recoding of variables done at the
#' end of the Monte Carlo loop. Needed for operations like counting the number of
#' days with a treatment or creating lagged variables. This is executed at each end
#'  of the Monte Carlo g-formula time steps.

#' @param mc.sample Sample size of Monte Carlo simulation.
#'
#' @references
#' Lin, S. H., Young, J. G., Logan, R., & VanderWeele, T. J. (2017). Mediation analysis for a survival outcome with time‐varying exposures, mediators, and confounders. \emph{Statistics in medicine}, 36(26), 4153-4166. DOI:10.1002/sim.7426
#' Zheng, W., & van der Laan, M. (2017). Longitudinal mediation analysis with time-varying mediators and exposures, with application to survival outcomes. \emph{Journal of causal inference}, 5(2). DOI:10.1515/jci-2016-0006
#'
#'
#' TODO: weights, time varying intervention
#'
#' @export
#'
#'

Gformula <- function(data,
                     id.var,
                     base.vars,
                     exposure,
                     outcome,
                     time.var,
                     models,
                     intervention = NULL,
                     out.recode = NULL,
                     init.recode = NULL,
                     mc.sample = 10000){

  tpcall <- match.call()

  # Setting seeds
  old <- .Random.seed
  on.exit( { .Random.seed <<- old } ) # Set to it's roginal seed after exit
  set.seed(2)

  if(!all.equal(sort(unique(na.omit(data[[outcome]]))), c(0, 1)))
    stop('Only binary treatments/exposures are currently implemented, and must be set to 0 and 1.')

  if(!all.equal(sort(unique(na.omit(data[[outcome]]))), c(0, 1)))
    # TODO: competing risk
    stop('Only binary outcomes are currently supported, and must be set to 0 and 1.')

  # TODO: weights
  # if(!is.null(weights))
  #   stop('Weights currently not supported')

  # Check if contain mediator
  med_flag <- sapply(models, function(mods) as.numeric(mods$type == "mediator"))
  #med_flag <- which(med_flag == 1)
  if(!all(med_flag == 0) & !is.null(intervention))
    stop("You cannot specify intervention with mediator at the same time!")

  # Checking for interventions
  if(!is.null(intervention)){
    if(length(intervention) == 1){
      intervention <- rep(intervention, length(unique(na.omit(data[[time.var]]))))
    }

    if(is.null(ncol(intervention))){
      # If only one intervention
      if(length(unique(na.omit(data[[time.var]]))) != length(intervention)){
        stop("Length of inteverntion must be the same as time length!")
      }else{
        num.intervene <- 1
      }

    }else{
      # More than one intervention
      if(length(unique(na.omit(data[[time.var]]))) != nrow(intervention)){
        stop("Length of inteverntion must be the same as time length!")
      }else{
        num.intervene <- ncol(intervention)
      }
    }
  }

  if(is.null(intervention)) num.intervene <- 1

  # The following are model checking and model evaluations

  # Models checking
  # if(class(models) != "g.model" & !is.list(models))
  #   stop("Models must be given as a list")
  # if(class(models) != "g.model"){
  #   mods <- list()
  #   for (indx in seq_along(models)) {
  #     mods[[indx]] <- do.call(spec_model, list(models[[indx]]))
  #   }
  #   models <- mods
  # }

  # Run along the models
  # TODO: Do need to consider the same ordering??
  if(!is.null(models)){
    fit_mods <- vector("list", length = length(models))
    for(indx_mod in seq_along(models)){
      mods <- models[[indx_mod]]
      # Check for model recode
      # No need for recode during the modelling process.
      # if(!is.null(mods$recode)){
      #   for(indx_rec in seq_along(mods$recode)){
      #     data <- within(data, eval(mods$recode[[indx_rec]]))
      #   }
      # }
      ord <- mods$order
      mods$call$data <- data

      fit_mods[[ord]] <- list(mods    = eval(mods$call),
                              recodes = mods$recode,
                              subset  = mods$subset,
                              family  = mods$family,
                              type    = mods$type)

    }
  }

  # Baseline variables for Monte Carlo Sampling
  base_dat <- unique(data[, unique(c(id.var, base.vars))])
  df_mc <- base_dat[sample(1:nrow(base_dat), mc.sample, replace = TRUE), ]

  # Mediation analysis
  if(!all(med_flag == 0)){
    y00 <- monte_g(data = df_mc, time.var = time.var, time.seq = unique(data[[time.var]]),
                   models = fit_mods, intervention = 0,
                   out.recode = out.recode, init.recode = init.recode)
    y11 <- monte_g(data = df_mc, time.var = time.var, time.seq = unique(data[[time.var]]),
                   models = fit_mods, intervention = 1,
                   out.recode = out.recode, init.recode = init.recode)
    y10 <- monte_g(data = df_mc, time.var = time.var, time.seq = unique(data[[time.var]]),
                   models = fit_mods, mediation = TRUE,
                   out.recode = out.recode, init.recode = init.recode)

    nie   <- mean(y11[[outcome]], na.rm = T) - mean(y10[[outcome]], na.rm = T)
    nde   <- mean(y10[[outcome]], na.rm = T) - mean(y00[[outcome]], na.rm = T)
    total <- mean(y11[[outcome]], na.rm = T) - mean(y00[[outcome]], na.rm = T)

    res <- list(always   = y11,
                never    = y00,
                mediator = y10)

  }else{
    if(num.intervene == 1){
      res <- monte_g(data = df_mc, time.var = time.var, time.seq = unique(data[[time.var]]),
                     models = fit_mods, intervention = intervention,
                     out.recode = out.recode, init.recode = init.recode)
    }else{
      res <- lapply(1:num.intervene, function(i){
        monte_g(data = df_mc, time.var = time.var, time.seq = unique(data[[time.var]]),
                models = fit_mods, intervention = intervention[, i],
                out.recode = out.recode, init.recode = init.recode)
      })
    }
  }

  return(list(call = tpcall,
              out  = res))

}

#' Monte Carlo simulation
#'
#' @description
#'  Internal use only. Monte Carlo simulation.
#'
#' @param data a data frame in which to look for variables with which to predict.
#'
#' @param time.seq Time sequence.
#'
#' @param time.var Time varaible.
#'
#' @param models fitted objects.
#'
#' @param intervention A vector, intervention treatment per time.
#'
#' @param init.recode optional, recoding of variables done at the
#' begaining of the Monte Carlo loop. Needed for operations initalize baseline variables.
#' This is executed at begaining of the Monte Carlo g-formula, excuted only once.
#'
#' @param out.recode optional, On the fly recoding of variables done at the
#' end of the Monte Carlo loop. Needed for operations like counting the number of
#' days with a treatment or creating lagged variables. This is executed at each end
#'  of the Monte Carlo g-formula time steps.
#'
#' @export
#'
monte_g <- function(data, time.seq, time.var, models,
                    intervention = NULL, mediation = FALSE,
                    out.recode = NULL, init.recode = NULL){

  # Time varying intervention
  if(!is.null(intervention) & length(intervention) == 1){
    intervention <- rep(intervention, length(unique(na.omit(data[[time.var]]))))
  }

  # Get the position of exposure
  exp_flag <- sapply(models, function(mods) as.numeric(mods$type == "exposure"))
  exp_flag <- which(exp_flag == 1)

  if(mediation){
    # Get the position of the mediator
    med_flag <- sapply(models, function(mods) as.numeric(mods$type == "mediator"))
    med_flag <- which(med_flag == 1)
  }

  # Get the position of the outcome
  out_flag <- sapply(models, function(mods) as.numeric(mods$type == "outcome"))
  out_flag <- which(out_flag == 1)

  # Get the position of the censor
  cen_flag <- sapply(models, function(mods) as.numeric(mods$type == "censor"))
  cen_flag <- which(cen_flag == 1)

  # Get the variable name of exposure, outcome and censor
  exposure <- all.vars(formula(models[[exp_flag]]$mods)[[2]])
  outcome  <- all.vars(formula(models[[out_flag]]$mods)[[2]])
  if(!all(cen_flag == 0)){
    censor <- all.vars(formula(models[[cen_flag]]$mods)[[2]])
  }

  # Simulate
  max_time <- max(time.seq, na.rm = TRUE)
  min_time <- min(time.seq, na.rm = TRUE)
  dat_y  <- data  # Always interven

  for(t in min_time:max_time){
    dat_y[[time.var]] <- t

    for(indx in seq_along(models)){
      resp_var <- all.vars(formula(models[[indx]]$mods)[[2]])

      # Intervention
      if(!is.null(intervention)){
        interv <- parse(text = paste0(exposure, " = ", intervention[t - min_time + 1]))
      }else{
        interv <- NULL
      }

      # Recode baseline variables
      if(t == min_time){
        if(!is.null(init.recode)){
          for(i in seq_along(init.recode)){
            dat_y <- within(dat_y, eval(parse(text = init.recode[i])))
          }
        }
      }

      if(!mediation){
        dat_y[[resp_var]] <- monte_sim(dat_y, models[[indx]], interv)
      }else{
        # For mediator.
        # All the variables between exposure and mediator are set to 0, else to 1
        interv1 <- parse(text = paste0(exposure, " = 1"))
        interv0 <- parse(text = paste0(exposure, " = 0"))
        if(indx <= med_flag){
          dat_y[[resp_var]] <- monte_sim(dat_y, models[[indx]], interv0)
        }else{
          dat_y[[resp_var]] <- monte_sim(dat_y, models[[indx]], interv1)
        }
      }
      # End for mediator

      # Recode data
      for(indx_rec in seq_along(models[[indx]]$recode)){
        dat_y <- within(dat_y, eval(parse(text = models[[indx]]$recode[[indx_rec]])))
      }

    }


    if(!is.null(out.recode)){
      for(i in seq_along(out.recode)){
        dat_y <- within(dat_y, eval(parse(text = out.recode[i])))
      }
    }

    # Output censoring and death, only the not intervened have censoring
    if(t == min_time){
      if(cen_flag != 0 & !is.null(interv)){
        out_y <- dat_y[dat_y[[outcome]] == 1 | dat_ynat[[censor]] == 1, ]
      }else{
        out_y <- dat_y[dat_y[[outcome]] == 1, ]
      }

    }else{
      if(cen_flag != 0 & !is.null(interv)){
        out_y <- rbind(out_y, dat_y[dat_y[[outcome]] == 1 | dat_y[[censor]] == 1, ])
      }else{
        out_y <- rbind(out_y, dat_y[dat_y[[outcome]] == 1, ])
      }
    }


    # loop not censored and dead, only the not intervened have censoring
    if(cen_flag != 0 & !is.null(interv)){
      dat_y <- dat_y[dat_y[[outcome]] != 1 & dat_y[[censor]] != 1, ]
    }else{
      dat_y <- dat_y[dat_y[[outcome]] != 1, ]
    }

    # if loop ends
    if(nrow(dat_y) == 0)
      break
    if(t == max_time){
      out_y <- rbind(out_y, dat_y)
    }
  }

  # loop ends here
  return(out_y)
}

#' Random data simulation from predicted value.
#'
#' @description
#'  Internal use only. Predict response and generate random data.
#'
#' @param data a data frame in which to look for variables with which to predict.
#'
#' @param models fitted objects.
#'
#' @param intervention Intervention.
#'
#' @export


monte_sim <- function(data, models, intervention = NULL){

  family <- models$family
  fit <- models$mods

  # Intervein
  if(is.null(intervention)){
    newdt <- data
  }else{
    newdt <- within(data, eval(intervention))
  }

  # Randomg number generation for Monte Carlo simulation
  if(family == "multinomial"){
    pred <- predict(fit, newdata = newdt, type =  "probs")
    rMultinom(pred, 1)
  }else if(family == "binomial"){
    pred <- predict(fit, newdata = newdt, type = "response")
    rbinom(nrow(data), 1, pred)
  }else{
    pred <- predict(fit, newdata = newdt, type = "response")
    rnorm(nrow(data), pred, sd(fit$residuals, na.rm = TRUE))
  }
}

#' Model specification for G-formula
#'
#' @description
#'  Add a specified regression model for the exposure. This is used for natural course estimation
#'  of the Monte Carlo g-formula. This must be specified before calling the fit function.
#'
#' @param formula Formula for specified models passed to \code{\link[stats]{glm}}.
#' Must be contained within the input  dataframe when initialized.
#'
#' @param subset an optional vector specifying a subset of observations to be used
#'  in the fitting process (see \code{\link[stats]{glm}}).
#'
#' @param family A description of the error distribution and link function to be used in the \code{glm} model.
#'
#' @param order Numeric, temporal ordering of the Covariates.
#'
#' @param type Model type. Exposure model (\code{exposure}), covariate model (\code{covriate}), mediator model (\code{mediator}), outcome model (\code{outcome}) or censoring model (\code{censor})
#'
#' @param recode Optional, recoding expression list. Used for dynamic recoding for Monte Carlo simulation.
#'  This is executed at after model predicted simulation.
#'
#' @import nnet
#'
#' @export

spec_model <- function(formula,
                       subset   = NULL,
                       family   = "gaussian",
                       order,
                       type     = c("exposure", "covriate", "mediator", "outcome", "censor"),
                       recode   = NULL){

  tmpcall <- match.call()

  if (!tmpcall$family %in% c("binomial", "multinomial", "gaussian"))
    stop("No valid family specified (\"binomial\", \"multinomial\", \"gaussian\")")

  if(tmpcall$type %in% c("exposure", "outcome", "censor") & tmpcall$family != "binomial")
    stop("Only binomial family supported for (\"exposure\", \"outcome\", \"censor\")")

  # if(!is.null(tmpcall$recode) & !is.list(tmpcall$recode))
  #   stop("Recode must be provided as list!")

  if(!is.numeric(tmpcall$order))
    stop("Order must be numeric!")

  if(family == "multinomial"){
    out_model <- quote(multinom())
  }else{
    out_model <- quote(glm())
    out_model$family <- substitute(family)
  }
  out_model$formula <- substitute(formula)
  out_model$subset  <- substitute(subset)

  out <- list(call   = out_model,
              order  = tmpcall$order,
              type   = tmpcall$type,
              family = tmpcall$family,
              recode = tmpcall$recode)

  class(out) <- "g.model"

  return(out)
}




