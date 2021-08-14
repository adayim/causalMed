
.gformula <- function(data,
                      id.var,
                      base.vars,
                      time.var,
                      exposure,
                      models,
                      intervention,
                      in.recode,
                      out.recode,
                      init.recode,
                      mediation_type = c(NA, "N", "I"),
                      mc.sample      = 10000,
                      verbose        = TRUE,
                      return_fitted  = FALSE){

  mediation_type <- match.arg(mediation_type)

  fit_mods <- lapply(models, function(mods) {

    rsp_vars <- all.vars(formula(mods$call)[[2]])

    # Observed values range
    if(is.numeric(data[[rsp_vars]])){
      val_ran <- range(na.omit(data[[rsp_vars]]))
    }else {
      val_ran <- unique(na.omit(data[[rsp_vars]]))
    }

    # Recode data before simulating
    if (!is.null(mods$recode)) {
      for (i in seq_along(mods$recode)) {
        data <- within(data, eval(parse(text = mods$recode[i])))
      }
    }

    mods$call$data <- substitute(data)

    list(
      fitted = eval(mods$call),
      recodes = mods$recode,
      subset = mods$subset,
      var_type = mods$var_type,
      mod_type = mods$mod_type,
      custom_sim = mods$custom_sim,
      rsp_vars = rsp_vars,
      val_ran = val_ran
    )
  })

  # Baseline variables for Monte Carlo Sampling
  base_dat <- unique(data[, unique(c(id.var, base.vars))])
  df_mc <- data.table::as.data.table(base_dat[sample(1:nrow(base_dat), mc.sample, replace = TRUE), ])
  df_mc[, new_ID:=seq_len(.N)]

  if(!is.na(mediation_type))
    intervention <- list(always = 1, never = 0, mediation = NULL)

  res <- sapply(intervention, function(i) {
    monte_g(
      data = df_mc,
      time.var = time.var,
      time.seq = unique(data[[time.var]]),
      exposure = exposure,
      models = fit_mods,
      intervention = i,
      in.recode = in.recode,
      out.recode = out.recode,
      init.recode = init.recode,
      mediation_type = mediation_type,
      verbose = verbose
    )
  }, simplify = FALSE)

  if(return_fitted){
    return(list(fitted.models = fit_mods, gform.data = res))
  } else {
    return(res)
  }
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
#' @importFrom parallel makePSOCKcluster clusterSetRNGStream stopCluster
#' @importFrom pbapply pblapply
#'
bootstrap_helper <- function(data,
                             id.var,
                             base.vars,
                             time.var,
                             exposure,
                             models,
                             intervention,
                             init.recode,
                             in.recode,
                             out.recode,
                             mc.sample      = 10000,
                             verbose        = TRUE,
                             mediation_type = c(NA, "N", "I"),
                             R              = 500,
                             parallel       = FALSE,
                             ncpus          = 1L) {

  mediation_type <- match.arg(mediation_type)

  # Check parallel
  if (parallel && ncpus > 1L) {
    if (.Platform$OS.type != "windows") {
      have_mc <- TRUE
    } else {
      have_mc <- FALSE
    }
    loadNamespace("parallel") # get this out of the way before recording seed
  }

  if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) runif(1)
  seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)

  cat("Be patient, bootstrap is running...\n")
  dfm <- lapply(1:R, function(x) data[sample(1:nrow(data), nrow(data), replace = TRUE), ])

  boot_res <- if (parallel && ncpus > 1L) {
    if (have_mc) {
      cl = ncpus
    } else {
      list(...) # evaluate any promises
      cl <- parallel::makePSOCKcluster(rep("localhost", ncpus))
      if (RNGkind()[1L] == "L'Ecuyer-CMRG") {
        parallel::clusterSetRNGStream(cl)
      }
      # clusterExport(cl, c("nobs", "iris_2species"))
      on.exit(close(parallel::stopCluster(cl)))
    }
  } else {
    cl <- NULL
  }

  boot_res <- pbapply::pblapply(dfm, .gformula,
                               id.var = id.var,
                               base.vars = base.vars,
                               time.var = time.var,
                               exposure = exposure,
                               models = models,
                               intervention = intervention,
                               in.recode = in.recode,
                               out.recode = out.recode,
                               init.recode = init.recode,
                               mc.sample      = mc.sample,
                               verbose        = FALSE,
                               mediation_type = mediation_type,
                               cl = cl)

  return(boot_res)
}
