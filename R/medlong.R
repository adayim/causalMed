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
#' @param id.var ID variable per subject.
#' @param base.vars A vector of time fixed baseline variables.
#' @param exposure Intervention/Exposure variable
#' @param outcome Name of the outcome variable.
#' @param time.var Time variable.
#' @param models A list of models for the G-formula, including exposure model,
#'  covariate model (if any), mediator model (if any), outcome model or
#'  censoring model (if any). See details in \code{\link[causalMed]{spec_model}}.
#' @param mc.sample Sample size of Monte Carlo simulation.
#'
#' @references
#' Lin, S. H., Young, J. G., Logan, R., & VanderWeele, T. J. (2017). Mediation analysis for a survival outcome with time‐varying exposures, mediators, and confounders. \emph{Statistics in medicine}, 36(26), 4153-4166. DOI:10.1002/sim.7426
#' Zheng, W., & van der Laan, M. (2017). Longitudinal mediation analysis with time-varying mediators and exposures, with application to survival outcomes. \emph{Journal of causal inference}, 5(2). DOI:10.1515/jci-2016-0006
#'
#'
#' TODO: weights, time varying intervention

Gformula <- function(data,
                     id.var,
                     base.vars,
                     exposure,
                     outcome,
                     time.var,
                     models,
                     mc.sample = 10000){

  tpcall <- match.call()

  if(!unique(na.omit(data[[exposure]])) %in% c(0, 1))
    stop('Only binary treatments/exposures are currently implemented, and must be set to 0 and 1.')

  if(!unique(na.omit(data[[outcome]])) %in% c(0, 1))
    # TODO: competing risk
    stop('Only binary outcomes are currently supported, and must be set to 0 and 1.')

  # TODO: weights
  # if(!is.null(weights))
  #   stop('Weights currently not supported')

  # The following are model checking and model evaluations

  # Models checking
  if(class(models) != "g.model" & !is.list(models))
    stop("Models must be given as a list")
  if(class(models) != "g.model"){
    mods <- list()
    for (indx in seq_along(models)) {
      mods[[indx]] <- do.call(spec_model, list(models[[indx]]))
    }
    models <- mods
  }

  # Run along the models
  # TODO: Do need to consider the same ordering??
  if(!is.null(models)){
    fit_covs <- vector("list", length = length(models))
    for(indx_mod in seq_along(models)){
      mods <- models[[indx_mod]]
      # Check for model recode
      if(!is.null(mods$recode)){
        for(indx_rec in seq_along(mods$recode)){
          data <- within(data, eval(mods$recode[[indx_rec]]))
        }
      }
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
  res <- monte_g(data = df_mc, time.var = time.var, models = fit_mods)

  # For mediation analysis
  if("mediator" %in% names(res)){
    nie   <- mean(res$always[[outcome]], na.rm = T) - mean(res$mediator[[outcome]], na.rm = T)
    nde   <- mean(res$mediator[[outcome]], na.rm = T) - mean(res$never[[outcome]], na.rm = T)
    total <- mean(res$always[[outcome]], na.rm = T) - mean(res$never[[outcome]], na.rm = T)
    res$mediation_effect <- c("Total effect"    = total,
                              "Indirect effect" = nie,
                              "Direct effect"   = nde)
  }

  return(res)

}

#' Monte Carlo simulation
#'
#' @description
#'  Internal use only. Monte Carlo simulation.
#'
#' @param data a data frame in which to look for variables with which to predict.
#'
#' @param time.var Time varaible.
#'
#' @param models fitted objects.
#'
#' @param intervention Intervention vector.
#'

monte_g <- function(data, time.var, models, intervention = NULL){

  # TODO: Time varying intervention

  # Get the position of exposure
  exp_flag <- sapply(models, function(mods) as.numeric(mods$type == "exposure"))
  exp_flag <- which(exp_flag == 1)

  # Get the position of the mediator
  med_flag <- sapply(models, function(mods) as.numeric(mods$type == "mediator"))
  med_flag <- which(med_flag == 1)

  # Get the position of the outcome
  out_flag <- sapply(models, function(mods) as.numeric(mods$type == "outcome"))
  out_flag <- which(out_flag == 1)

  # Get the position of the censor
  cen_flag <- sapply(models, function(mods) as.numeric(mods$type == "censor"))
  cen_flag <- which(cen_flag == 1)

  # Get the variable name of exposure, outcome and censor
  exposure <- all.vars(formula(models[[exp_flag]]$mods)[[2]])
  outcome  <- all.vars(formula(models[[out_flag]]$mods)[[2]])
  if(cen_flag != 0){
    censor <- all.vars(formula(models[[cen_flag]]$mods)[[2]])
  }

  # Simulate
  max_time <- max(data[[time.var]], na.rm = TRUE)
  min_time <- min(data[[time.var]], na.rm = TRUE)
  dat_y11  <- data  # Always interven
  dat_y10  <- data  # Mediator intervention to 0, others to 1
  dat_y00  <- data  # Never interven
  dat_ynat <- data  # Natrual

  for(t in min_time:max_time){
    dat_y11[[time.var]] <- t
    dat_y10[[time.var]] <- t
    dat_y00[[time.var]] <- t
    dat_ynat[[time.var]] <- t

    for(indx in seq_along(models)){
      resp_var <- all.vars(formula(models[[indx]]$mods)[[2]])
      if(indx >= exp_flag){
        interv1 <- parse(text = paste0(exposure, " = 1"))
        interv0 <- parse(text = paste0(exposure, " = 0"))
      }else{
        interv1 <- interv0 <- NULL
      }
      for(indx_rec in seq_along(models[[indx]]$recode)){
        dat_y11 <- within(dat_y11, eval(models[[indx]]$recode[[indx_rec]]))
        dat_y00 <- within(dat_y00, eval(models[[indx]]$recode[[indx_rec]]))
        dat_ynat <- within(dat_ynat, eval(models[[indx]]$recode[[indx_rec]]))
        dat_y10 <- within(dat_y10, eval(models[[indx]]$recode[[indx_rec]]))
      }
      dat_y11[[resp_var]] <- monte_sim(dat_y11, models[[indx]]$mods, interv1)
      dat_y00[[resp_var]] <- monte_sim(dat_y00, models[[indx]]$mods, interv0)
      dat_ynat[[resp_var]] <- monte_sim(dat_ynat, models[[indx]]$mods)

      # For mediator.
      # All the variables between exposure and mediator are set to 0, else to 1
      if(med_flag != 0){
        if(indx <= med_flag){
          dat_y10[[resp_var]] <- monte_sim(dat_y10, models[[indx]]$mods, interv0)
        }else{
          dat_y10[[resp_var]] <- monte_sim(dat_y10, models[[indx]]$mods, interv1)
        }
      }
      # End for mediator
    }

    # output censoring and death
    if(t == min_time){
      out_y11 <- dat_y11[dat_y11[[outcome]] == 1, ]
      out_y10 <- dat_y10[dat_y10[[outcome]] == 1, ]
      out_y00 <- dat_y00[dat_y00[[outcome]] == 1, ]
      if(cen_flag != 0){
        out_ynat <- dat_ynat[dat_ynat[[outcome]] == 1 | dat_ynat[[censor]] == 1, ]
      }else{
        out_ynat <- dat_ynat[dat_ynat[[outcome]] == 1, ]
      }
    }else{
      out_y11 <- rbind(out_y11, dat_y11[dat_y11[[outcome]] == 1, ])
      out_y10 <- rbind(out_y10, dat_y10[dat_y10[[outcome]] == 1, ])
      out_y00 <- rbind(out_y00, dat_y00[dat_y00[[outcome]] == 1, ])
      if(cen_flag != 0){
        out_ynat <- rbind(out_ynat, dat_ynat[dat_ynat[[outcome]] == 1 | dat_ynat[[censor]] == 1, ])
      }else{
        out_ynat <- rbind(out_ynat, dat_ynat[dat_ynat[[outcome]] == 1, ])
      }
    }

    # loop not censored and dead
    dat_y11 <- dat_y11[dat_y11[[outcome]] != 1, ]
    dat_y10 <- dat_y10[dat_y10[[outcome]] != 1, ]
    dat_y00 <- dat_y00[dat_y00[[outcome]] != 1, ]

    if(cen_flag != 0){
      dat_ynat <- dat_ynat[dat_ynat[[outcome]] != 1 & dat_ynat[[censor]] != 1, ]
    }else{
      dat_ynat <- dat_ynat[dat_ynat[[outcome]] != 1, ]
    }

    # if loop ends
    if(nrow(next_y11) == 0 & nrow(next_y10) == 0 & nrow(next_y00) == 0)
      break
    if(t == max_time){
      out_y11 <- rbind(out_y11, dat_y11)
      out_y10 <- rbind(out_y10, dat_y10)
      out_y00 <- rbind(out_y00, dat_y00)
      out_ynat <- rbind(out_ynat, dat_ynat)
    }
  }
  # loop ends here

  if(med_flag != 0){
    return(list(natural  = out_ynat,
                never    = out_y00,
                always   = out_y11,
                mediator = out_y10))
  }else{
    return(list(natural  = out_ynat,
                never    = out_y00,
                always   = out_y11))
  }

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
#' @param model Formula for specified models passed to \code{\link[stats]{glm}}.
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
#'
#' @import nnet multinom
#'
#' @export

spec_model <- function(model,
                       subset   = NULL,
                       family   = "gaussian",
                       order,
                       type     = c("exposure", "covriate", "mediator", "outcome", "censor"),
                       recode   = NULL){

  if (!family %in% c("binomial", "multinomial", "gaussian"))
    stop("No valid family specified (\"binomial\", \"multinomial\", \"gaussian\")")

  if(type %in% c("exposure", "outcome", "censor") & family != "binomial")
    stop("Only binomial family supported for (\"exposure\", \"outcome\", \"censor\")")

  if(!is.null(recode) & !is.list(recode))
    stop("Recode must be provided as list!")

  if(!is.numeric(order))
    stop("Order must be numeric!")

  if(family == "multinomial"){
    out_model <- quote(nnet::multinom())
  }else{
    out_model <- quote(glm())
    out_model$family <- substitute(family)
  }
  out_model$formula <- substitute(model)
  out_model$subset  <- substitute(subset)

  out <- list(call   = out_model,
              order  = order,
              type   = type,
              recode = recode)

  class(out) <- "g.model"

  return(out)
}

###### TODO Below ============

#' @rdname medlong-methods
#' @method summary medlong
#' @export
summary.medlong <- function (object, ...){

  nams <- names(object)

  # Extract bootstrap results
  extrcF <- function(x, ..){
    coef.table <- extract_boot(x, ...)
    coeftab <- as.matrix(coef.table[, -1])
    dimnames(coeftab) <- list(coef.table$term,
                              c("Estimate", "Bias", "Std.error", "conf.low", "conf.high"))
    prop.med <- coef.table[coef.table$term == "Natural Indirect effect", "statistic"] /
      coef.table[coef.table$term == "Total effect", "statistic"]
    list(coeftab = coeftab, prop.med = prop.med)
  }

  if("iptw" %in% nams){
    iptw.res <- extrcF(object$iptw, ...)
  }else{
    iptw.res <- list(coeftab = NULL, prop.med = NULL)
  }

  if("g-formula" %in% nams){
    gform.res <- extrcF(object$`g-formula`, ...)
  }else{
    gform.res <- list(coeftab = NULL, prop.med = NULL)
  }

  summary <- list(call         = object$call,
                  coeff.iptw   = iptw.res$coeftab,
                  coeff.gform  = gform.res$coeftab,
                  prop.iptw    = 100 * iptw.res$prop.med,
                  prop.gform   = 100 * gform.res$prop.med)

  class(summary) <- "summary.medlong"
  attr(summary, "class_object") <- class(object)
  return(summary)
}

#' @method print summary.medlong
#' @export
print.summary.medlong <- function (x, digits = max(3, getOption("digits") - 3), ...){
  cat("Call:\n")
  print(x$call)

  cat("------\n")
  cat("Exposure:", x$call$exposure, "\nMediator:", paste(x$call$mediator,
                                                         collapse = ", "))

  cat("\n---\nEstimation of standard errors based on the non-parametric bootstrap")
  if (!is.null(x$coeff.gform)){
    cat("\n------\n")
    cat("Natural effect parameter estimates with G-formula estimation:\n")
    print(round(x$coeff.gform, digits = digits))
    cat(paste0("------\nProportion Mediated: ", round(x$prop.gform, digits), "%"))
  }

  if (!is.null(x$coeff.iptw)){
    cat("\n------\n")
    cat("\nNatural effect parameter estimates with IPTW estimation:\n")
    print(round(x$coeff.iptw, digits = digits))
    cat(paste0("------\nProportion Mediated: ", round(x$prop.iptw, digits), "%"))
  }
}





