
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
#' @param intervention A named list with a value of intervention on exposure.
#' if kept as NULL (default), the natural intervention course will be calculated.
#'  eg: list(natural = NULL, always = c(1, 1, 1), never = c(0, 0, 0))
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
                     init.recode  = NULL,
                     in.recode    = NULL,
                     out.recode   = NULL,
                     mc.sample    = 10000,
                     verbose      = TRUE){

  tpcall <- match.call()

  # Setting seeds
  if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) runif(1)
  seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)

  # Check for error
  do.call(check_error, tpcall)

  if(is.null(intervention)){
    intervention <- list(intervention = NULL)
  }

  # Run along the models
  fit_mods <- lapply(models, function(mods){
    mods$call$data <- substitute(data)
    rsp_vars <- all.vars(formula(mods$call)[[2]])
    is_numeric <- is.numeric(data[[rsp_vars]])

    list(fitted    = eval(mods$call),
         recodes  = mods$recode,
         subset   = mods$subset,
         var_type = mods$var_type,
         mod_type = mods$mod_type,
         custom_sim = mods$custom_sim,
         rsp_vars = rsp_vars,
         val_ran  = ifelse(is_numeric, range(na.omit(data[[rsp_vars]])),
                           unique(na.omit(data[[rsp_vars]])))) # Observed values range
  })

  # Baseline variables for Monte Carlo Sampling
  base_dat <- unique(data[, unique(c(id.var, base.vars))])
  df_mc <- base_dat[sample(1:nrow(base_dat), mc.sample, replace = TRUE), ]

  # Mediation analysis
  if(!all(med_flag == 0))
    intervention <- list(always    = rep(1, time_len),
                         never     = rep(0, time_len),
                         mediation = rep("mediation", time_len))

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
#' @keywords internal
#'
monte_g <- function(data,
                    time.seq,
                    time.var,
                    models,
                    intervention = NULL,
                    in.recode = NULL,
                    out.recode = NULL,
                    init.recode = NULL,
                    verbose = TRUE){

  # Replicate to intervention to the same length of the time.
  time_len <- length(time.seq)
  if(all(sapply(intervention, length) == 1))
    intervention <- lapply(intervention, function(x)rep(x, time_len))

  # Get the position of the mediator
  med_flag <- sapply(models, function(mods) as.numeric(mods$mod_type == "mediator"))
  med_flag <- which(med_flag == 1)
  med_flag <- ifelse(is.null(intervention), 0,
                    ifelse(as.character(intervention) == "mediation", med_flag, 0))

  # Get the position of the outcome
  out_flag <- sapply(models, function(mods) as.numeric(mods$mod_type == "outcome"))
  out_flag <- which(out_flag == 1)

  # Get the position of the censor
  cen_flag <- sapply(models, function(mods) as.numeric(mods$mod_type == "censor"))
  surv_flag <- sapply(models, function(mods) as.numeric(mods$mod_type == "survival"))
  if(!all(cen_flag == 0) | !all(surv_flag == 0)){
    is.survival <- TRUE
  }
  cen_flag <- which(cen_flag == 1)

  # Get the variable name of outcome and censor
  outcome  <- all.vars(formula(models[[out_flag]]$fitted)[[2]])
  if(!all(cen_flag == 0)){
    censor <- all.vars(formula(models[[cen_flag]]$fitted)[[2]])
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
  for(t_index in min_time:max_time){
    dat_y[[time.var]] <- t_index

    # Recode baseline variables at initiation
    if(t_index == min_time){
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

    # Use the model to calculate the simulated value
    dat_y <- simulate_data(data = dat_y, models = models, intervention = intervention[t_index])
    

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
      if(t_index == min_time){
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
      if(t_index == max_time){
        out_y <- rbind(out_y, dat_y)
      }
    }else{
      out_y <- dat_y
    }
  }
  # loop ends here

  return(as.data.frame(out_y))
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






