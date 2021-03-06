
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
#'  censoring model (if any). See details in \code{\link[causalMed]{spec_model}}.
#'
#' @param intervention A named list with a value of intervention on exposure.
#' if kept as NULL (default), the natrual intervention cousre will be calculated.
#'  eg: list(natural = NULL, always = 1, never = 0)
#'
#' @param init.recode optional, recoding of variables done at the
#' begaining of the Monte Carlo loop. Needed for operations initalize baseline variables.
#' This is executed at begaining of the Monte Carlo g-formula, excuted only once at time 0.
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
#'  @param is.survival Is the data survival data, defalt is FALSE.
#'
#' @param mc.sample Sample size of Monte Carlo simulation.
#'
#' @param verbose Print intervention information during calculation.
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
                     init.recode = NULL,
                     in.recode = NULL,
                     out.recode = NULL,
                     mc.sample = 10000,
                     verbose  = TRUE){

  tpcall <- match.call()

  # Setting seeds
  if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) runif(1)
  seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)

  cen_flag  <- sapply(models, function(mods) as.numeric(mods$type == "censor"))
  surv_flag <- sapply(models, function(mods) as.numeric(mods$type == "survival"))
  if(!all(cen_flag == 0) | !all(surv_flag == 0)){
    is.survival <- TRUE
  }

  if(!all.equal(sort(unique(na.omit(data[[exposure]]))), c(0, 1)))
    stop('Only binary treatments/exposures are currently implemented, and must be set to 0 and 1.')

  if(is.survival & !all.equal(sort(unique(na.omit(data[[outcome]]))), c(0, 1)))
    # TODO: competing risk
    stop('Only binary outcomes for survival are currently supported, and must be set to 0 and 1.')

  # TODO: weights
  # if(!is.null(weights))
  #   stop('Weights currently not supported')

  # Check if contain mediator
  med_flag <- sapply(models, function(mods) as.numeric(mods$type == "mediator"))
  #med_flag <- which(med_flag == 1)
  if(!all(med_flag == 0) & !is.null(intervention))
    stop("You cannot specify intervention with mediator at the same time!")

  # Checking for interventions
  if(is.null(intervention)){
    intervention <- list(intervention = NULL)
  }
  if(!is.list(intervention))
    stop("Intervention must be a named list object")


  # The following are model checking and model evaluations

  # Models checking
  if(!is.list(models))
    stop("Models must be provided as a list")
  cls <- sapply(models, function(x)class(x) != "gmodel")
  if(any(cls))
    stop("Models in the list must be gmodel object, please use spec_model to create!")

  # Run along the models
  # TODO: Do need to consider the same ordering??
  if(!is.null(models)){
    fit_mods <- vector("list", length = length(models))
    for(indx_mod in seq_along(models)){

      mods <- models[[indx_mod]]
      ord <- mods$order

      if(mods$family == "lcmm"){
        if(!"class" %in% names(data))
          stop("Variable 'class' indicating object membership must be included in the data if lcmm model provided!")

        if(!"class" %in% base.vars)
          base.vars <- c(base.vars, "class")

        rsp_vars <- all.vars(mods$call$fixed[[2]])

        # Observed values
        val_ran <- sapply(rsp_vars, function(x)unique(na.omit(data[[x]])))

        fit_mods[[ord]] <- list(mods     = mods$call,
                                # recodes  = mods$recode,
                                subset   = mods$subset,
                                family   = mods$family,
                                type     = mods$type,
                                val_ran  = val_ran)
      }else{
        mods$call$data <- substitute(data)

        rsp_vars <- all.vars(formula(mods$call)[[2]])

        fit_mods[[ord]] <- list(mods     = eval(mods$call),
                                # recodes  = mods$recode,
                                subset   = mods$subset,
                                family   = mods$family,
                                type     = mods$type,
                                val_ran  = unique(na.omit(data[[rsp_vars]]))) # Observed values
      }
    }
  }

  # Baseline variables for Monte Carlo Sampling
  base_dat <- unique(data[, unique(c(id.var, base.vars))])
  df_mc <- base_dat[sample(1:nrow(base_dat), mc.sample, replace = TRUE), ]

  # Mediation analysis
  if(!all(med_flag == 0))
    intervention <- list(always = 1, never = 0, mediation = "mediation")

  res <- sapply(intervention, function(i){
    monte_g(data = df_mc, time.var = time.var, time.seq = unique(data[[time.var]]),
            models = fit_mods, intervention = i,
            in.recode = in.recode,
            out.recode = out.recode, init.recode = init.recode,
            verbose = verbose)
  }, simplify = FALSE)

  return(list(call          = tpcall,
              fitted.models = fit_mods,
              gform.data    = res))

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
#' This is executed at begaining of the Monte Carlo g-formula, excuted only once at time 0.
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
#' @param verbose Print intervention information during calculation.
#'
#' @export
#'
monte_g <- function(data, time.seq, time.var, models,
                    intervention = NULL,
                    in.recode = NULL,
                    out.recode = NULL, init.recode = NULL,
                    verbose = TRUE){

  # Time varying intervention
  # if(!is.null(intervention) & length(intervention) == 1){
  #   intervention <- rep(intervention, length(unique(na.omit(data[[time.var]]))))
  # }

  # Get the position of the mediator
  med_flag <- sapply(models, function(mods) as.numeric(mods$type == "mediator"))
  med_flag <- which(med_flag == 1)
  med_flag <- ifelse(is.null(intervention), 0,
                    ifelse(as.character(intervention) == "mediation", med_flag, 0))

  # Get the position of exposure
  exp_flag <- sapply(models, function(mods) as.numeric(mods$type == "exposure"))
  exp_flag <- which(exp_flag == 1)

  # Get the position of the outcome
  out_flag <- sapply(models, function(mods) as.numeric(mods$type == "outcome"))
  out_flag <- which(out_flag == 1)

  # Get the position of the censor
  cen_flag <- sapply(models, function(mods) as.numeric(mods$type == "censor"))
  surv_flag <- sapply(models, function(mods) as.numeric(mods$type == "survival"))
  if(!all(cen_flag == 0) | !all(surv_flag == 0)){
    is.survival <- TRUE
  }

  cen_flag <- which(cen_flag == 1)


  # Get the variable name of exposure, outcome and censor
  exposure <- all.vars(formula(models[[exp_flag]]$mods)[[2]])
  outcome  <- all.vars(formula(models[[out_flag]]$mods)[[2]])
  if(!all(cen_flag == 0)){
    censor <- all.vars(formula(models[[cen_flag]]$mods)[[2]])
  }

  if(verbose){
    cat("\n=============\n")
    cat(paste0("Intervention: ", intervention))
  }

  # Simulate
  max_time <- max(time.seq, na.rm = TRUE)
  min_time <- min(time.seq, na.rm = TRUE)
  dat_y  <- data


  # Run normal g-formula
  for(t in min_time:max_time){
    dat_y[[time.var]] <- t

    # Recode baseline variables
    if(t == min_time){
      if(!is.null(init.recode)){
        for(i in seq_along(init.recode)){
          dat_y <- within(dat_y, eval(parse(text = init.recode[i])))
        }
      }
    }

    # Recode data before simulating
    if(!is.null(in.recode)){
      for(i in seq_along(in.recode)){
        dat_y <- within(dat_y, eval(parse(text = in.recode[i])))
      }
    }

    # Intervention
    if(!is.null(intervention) & med_flag == 0){
      dat_y[[exposure]] <- intervention
    }

    # For mediation
    if(med_flag != 0){
      dat_y[[exposure]] <- 1
      interv0 <- parse(text = paste0(exposure, " = 0"))
    }

    # Loop through variables
    for(indx in seq_along(models)){
      if(models[[indx]]$family == "lcmm"){
        resp_var <- all.vars(models[[indx]]$mods$call$fixed[[2]])

      }else{
        resp_var <- all.vars(formula(models[[indx]]$mods)[[2]])

        if(resp_var == exposure & !is.null(intervention))
          next   # Skip if variable is treatment
      }

      # If condition has been defined
      if(!is.null(models[[indx]]$subset)){
        cond <- with(dat_y, eval(models[[indx]]$subset))
      }else{
        cond <- rep(TRUE, nrow(dat_y))
      }
      if(sum(cond) != 0){
        if(med_flag == indx){
          # Set the mediator's intervention to 0
          if(models[[indx]]$family == "lcmm"){
            dat_y <- monte_sim.lcmm(within(dat_y[cond,  ], eval(interv0)),
                                    models[[indx]], time.var = time.var)

          }else{
            dat_y[[resp_var]][cond] <- monte_sim(within(dat_y[cond,  ], eval(interv0)),
                                                 models[[indx]])
          }

        }else{
          dat_y[[resp_var]][cond] <- monte_sim(dat_y[cond, ], models[[indx]])
        }
      }
    }

    if(!is.null(out.recode)){
      for(i in seq_along(out.recode)){
        dat_y <- within(dat_y, eval(parse(text = out.recode[i])))
      }
    }

    if(is.survival){
      if(!is.null(intervention) & cen_flag != 0){
        dat_y[[censor]] <- 0
      }

      # If censored outcome equals to 0
      if(cen_flag != 0){
        dat_y[[outcome]] <- ifelse(dat_y[[censor]] == 1, 0, dat_y[[outcome]])
      }

      # Output censoring and death, only the not intervened have censoring
      if(t == min_time){
        if(cen_flag != 0){
          out_y <- dat_y[dat_y[[outcome]] == 1 | dat_y[[censor]] == 1, ]
        }else{
          out_y <- dat_y[dat_y[[outcome]] == 1, ]
        }

      }else{
        if(cen_flag != 0){
          out_y <- rbind(out_y, dat_y[dat_y[[outcome]] == 1 | dat_y[[censor]] == 1, ])
        }else{
          out_y <- rbind(out_y, dat_y[dat_y[[outcome]] == 1, ])
        }
      }

      # loop not censored and dead, only the not intervened have censoring
      if(cen_flag != 0){
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
    }else{
      out_y <- dat_y
    }
  }
  # loop ends here

  return(as.data.frame(out_y))
}

#' Random data simulation from predicted value.
#'
#' @description
#'  Internal use only. Predict response and generate random data.
#'
#' @param newdt a data frame in which to look for variables with which to predict.
#'
#' @param models fitted objects.
#'
#' @export


monte_sim <- function(newdt, models){

  family <- models$family
  fit <- models$mods
  rmse <- sqrt(mean((stats::fitted(fit) - fit$y)^2, na.rm = TRUE))

  # Randomg number generation for Monte Carlo simulation
  if(family == "multinomial"){
    pred <- predict(fit, newdata = newdt, type =  "probs")
    rMultinom(pred, 1)
  }else if(family == "binomial"){
    pred <- predict(fit, newdata = newdt, type = "response")
    rbinom(nrow(newdt), 1, pred)
  }else{
    pred <- predict(fit, newdata = newdt, type = "response")
    out <- rnorm(nrow(newdt), pred, rmse)

    # Limit the predicted values within observed range.
    out <- pmax(out, min(models$val_ran))
    out <- pmin(out, max(models$val_ran))
    out
  }
}

monte_sim.lcmm <- function(newdt, models, time.var){

  fit <- models$mods
  resp_vars <- all.vars(fit$call$fixed[[2]])

  pred <- predictY(fit, newdata = newdt, draws = FALSE, var.time = time.var)$pred

  if(length(resp_vars) == 1){
    for(i in 1:ncol(pred)){
      newdt[newdt$class == i, resp_vars] <- pred[newdt$class == i, i]
    }

    newdt[[resp_vars]] <- newdt[[resp_vars]] + rnorm(length(newdt[[resp_vars]]))*fit$best["stderr"]

    # Limit the predicted values within observed range.
    newdt[[resp_vars]] <- pmax(newdt[[resp_vars]], min(models$val_ran[[resp_vars]]))
    newdt[[resp_vars]] <- pmin(newdt[[resp_vars]], max(models$val_ran[[resp_vars]]))

  }else{
    for(vars in resp_vars){
      for(i in 2:ncol(pred)){
        newdt[newdt$class == i-1, vars] <- pred[pred$Yname == vars & newdt$class == i-1, i]
      }
    }

    std_err <- fit$best[grepl("std.err", names(fit$best))]

    for(i in seq_along(resp_vars)){
      newdt[[resp_vars[i]]] <- newdt[[resp_vars[i]]] +
        rnorm(length(newdt[[resp_vars[i]]]))*std_err[i]

      # Limit the predicted values within observed range.
      newdt[[resp_vars[i]]] <- pmax(newdt[[resp_vars[i]]], min(models$val_ran[[resp_vars[i]]]))
      newdt[[resp_vars[i]]] <- pmin(newdt[[resp_vars[i]]], max(models$val_ran[[resp_vars[i]]]))
    }
  }

  return(newdt)

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
#' @param family A description of the error distribution and link function
#'  to be used in the \code{glm} model.
#'
#' @param order Numeric, temporal ordering of the Covariates.
#'
#' @param type Model type. Exposure model (\code{exposure}), covariate
#'  model (\code{covariate}), mediator model (\code{mediator}), outcome
#'   model (\code{outcome}), mediator model (\code{survival}) or censoring model (\code{censor})
#'
#' @import nnet
#'
#' @export

spec_model <- function(formula,
                       subset   = NULL,
                       family   = "gaussian",
                       order,
                       #recode   = NULL,
                       type     = c("exposure", "covariate", "mediator",
                                    "outcome", "censor", "survival")){

  tmpcall <- match.call()

  if(tmpcall$type %in% c("exposure", "outcome", "censor") & tmpcall$family != "binomial")
    stop("Only binomial family supported for (\"exposure\", \"outcome\", \"censor\")")

  if(!is.numeric(tmpcall$order))
    stop("Order must be numeric!")

  if(tmpcall$type != "mediator" & class(tmpcall$formula) %in% c("lcmm", "hlme", "multlcmm"))
    stop("LCMM object only supported for mediator model!")

  # For LCMM models
  if(class(tmpcall$formula) %in% c("lcmm", "hlme", "multlcmm")){
    out <- list(call   = tmpcall$formula,
                subset = tmpcall$formula$subset,
                order  = tmpcall$order,
                type   = tmpcall$type,
                #recode = tmpcall$recode,
                family = "lcmm")
  }else{
    # For none LCMM models
    if (!tmpcall$family %in% c("binomial", "multinomial", "gaussian"))
      stop("No valid family specified (\"binomial\", \"multinomial\", \"gaussian\")")

    # if(!is.null(tmpcall$recode) & !is.list(tmpcall$recode))
    #   stop("Recode must be provided as list!")

    if(family == "multinomial"){
      out_model <- quote(multinom())
    }else{
      out_model <- quote(glm())
      out_model$family <- substitute(family)
    }
    out_model$formula <- substitute(formula)
    out_model$subset  <- substitute(subset)

    out <- list(call   = out_model,
                subset = tmpcall$subset,
                order  = tmpcall$order,
                type   = tmpcall$type,
                #recode = tmpcall$recode,
                family = tmpcall$family)
  }

  class(out) <- "gmodel"

  return(out)
}

#' Calculate mediation analysis confidence interval
#'
#' @description Used to calculate confidence interval using non-parametric bootstrap methods.
#'
#' @param object Object returned from Gformula function (see \code{\link[calsalMed]{Gformula}}).
#'
#' @param R The number of bootstrap replicates, default is 500. Same with \code{boot}, see \code{\link[boot]{boot}} for detail.
#'
#' @param parallel If parallel operation to be used, the default is FALSE.
#'
#' @param ncores integer: number of processes to be used in parallel operation: typically one would chose this to the number of available CPUs.
#'
#' @export
#'
gformula_med_ci <- function(object, R = 500,
                            parallel = FALSE,
                            ncpus = 1L){

  if(as.character(object$call[[1]]) != "Gformula")
    stop("Object must be an Gformula!")
  if(!names(object$gform.data) %in% c("always", "never", "mediation"))
    stop("Only mediation analysis supported!")

  tmpcall <- match.call()

  # Check parallel
  if (parallel && ncpus > 1L) {
    if (.Platform$OS.type != "windows") {
      have_mc <- TRUE
    }else{
      have_mc <- FALSE
    }
    loadNamespace("parallel") # get this out of the way before recording seed
  }

  if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) runif(1)
  seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)

  obj_call <- object$call

  # Get original data from call
  data <- eval(obj_call$data)
  obj_call$verbose <- FALSE

  # Calculate proportion mediated
  out_come <- as.character(obj_call$outcome)

  if("class" %in% obj_call$base.vars){
    out.res <- lapply(unique(data[["class"]]), function(i){
      med <- mean(object$gform.data$mediation[object$gform.data$mediation$class == i, out_come])
      y11 <- mean(object$gform.data$always[object$gform.data$always$class == i, out_come])
      y00 <- mean(object$gform.data$never[object$gform.data$never$class == i, out_come])
      nde <- med - y00
      nie <- y11 - med

      c("Class"          = i,
        "Total Effect"   = nde + nie,
        "Direct Effect"  = nde,
        "Indirect Effct" = nie)
    })
    out.res <- do.call("rbind", out.res)

  }else{
    nde <- mean(object$gform.data$mediation[, out_come]) -
      mean(object$gform.data$never[, out_come])

    nie <- mean(object$gform.data$always[, out_come]) -
      mean(object$gform.data$mediation[, out_come])

    out.res <- c("Total Effect"   = nde + nie,
                 "Direct Effect"  = nde,
                 "Indirect Effct" = nie)

    # prop_med <- nie/(nie + nde)
  }

  # Create bootstrap function
  boot_func <- function(data){
    dat <- data[index, ]
    obj_call$data <- substitute(dat)
    out_come <- as.character(obj_call$outcome)
    out <- eval(obj_call)

    if("class" %in% obj_call$base.vars){
      res <- lapply(unique(data[["class"]]), function(i){
        med <- mean(object$gform.data$mediation[object$gform.data$mediation$class == i, out_come])
        y11 <- mean(object$gform.data$always[object$gform.data$always$class == i, out_come])
        y00 <- mean(object$gform.data$never[object$gform.data$never$class == i, out_come])
        nde <- med - y00
        nie <- y11 - med

        c("Class"          = i,
          "Total Effect"   = nde + nie,
          "Direct Effect"  = nde,
          "Indirect Effct" = nie)
      })
      res <- do.call("rbind", res)

    }else{
      nde <- mean(object$gform.data$mediation[, out_come]) -
        mean(object$gform.data$never[, out_come])

      nie <- mean(object$gform.data$always[, out_come]) -
        mean(object$gform.data$mediation[, out_come])

      res <- c("Total Effect"   = nde + nie,
               "Direct Effect"  = nde,
               "Indirect Effct" = nie)
    }
    return(res)
  }

  cat("Be patient, bootstrap is running...\n")
  dfm <- lapply(1:R, function(x)data[sample(1:nrow(data), nrow(data), replace = TRUE), ])

  boot_res <- if (parallel && ncpus > 1L) {
    if (have_mc) {
      parallel::mclapply(dfm, boot_func, mc.cores = ncpus)
    } else {
      list(...) # evaluate any promises
      cl <- parallel::makePSOCKcluster(rep("localhost", ncpus))
      if(RNGkind()[1L] == "L'Ecuyer-CMRG")
        parallel::clusterSetRNGStream(cl)
      res <- parallel::parLapply(cl, dfm, boot_func)
      parallel::stopCluster(cl)
      res
    }
  } else lapply(dfm, boot_func)

  conf <- do.call("rbind", boot_res)

  if("class" %in% obj_call$base.vars){
    conf <- lapply(unique(data[["class"]]), function(j){
      conf <- sapply(2:4, function(i){
        c(quantile(conf[conf$Class== j, i], prob=0.025),
          quantile(conf[conf$Class== j, i], prob=0.075))
      })
      conf <- t(conf)
      row.names(conf) <- c("Total Effect", "Direct Effect", "Indirect Effct")
      colnames(conf) <- c("LCL", "UCL")

      conf[, "Class"] <- j
      return(conf)
    })
    conf <- do.call("rbind", conf)

  }else{
    conf <- sapply(1:3, function(i){
      c(quantile(conf[, i], prob=0.025),
        quantile(conf[, i], prob=0.075))
    })

    conf <- t(conf)
    row.names(conf) <- c("Total Effect", "Direct Effect", "Indirect Effct")
    colnames(conf) <- c("LCL", "UCL")
  }

  out <- list(call      = tmpcall,
              estimated = out.res,
              confint   = conf)
  class(out) <- 'gmed_boot'
  return(out)
}


#' @method print gmed_boot
#'
#' @export
#'

print.gmed_boot <- function (x, digits = max(3, getOption("digits") - 3), ...){
  cat("Call:\n")
  print(x$call)

  cat("------\n")
  cat("Natural effect model\n")
  cat("with standard errors based on the non-parametric bootstrap\n---\n")
  cat("Parameter estimates:\n")
  out <- cbind(x$estimated, x$confint)
  colnames(out)[1] <- "Estimated"
  print(round(out, digits = digits))

  # cat(paste0("------\nProportion Mediated: ", round(x$prop_med, digits), "%"))
}






