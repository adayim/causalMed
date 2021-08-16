#' Calculate mediation analysis confidence interval
#'
#' @description Used to calculate confidence interval using non-parametric bootstrap methods.
#'
#' @param object Object returned from Gformula function (see \code{\link[calsalMed]{Gformula}}).
#'
#' @param R The number of bootstrap replicates, default is 500. Same with \code{boot}, see \code{\link[boot]{boot}} for detail.
#'
#' @param ncores integer: number of processes to be used in parallel operation: typically one would chose this to the number of available CPUs.
#'
#' @importFrom parallel makePSOCKcluster clusterSetRNGStream stopCluster clusterExport clusterEvalQ
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
                             mediation_type = c(NA, "N", "I"),
                             R              = 500,
                             ncores         = 1L) {

  mediation_type <- match.arg(mediation_type)

  # Check parallel
  if (ncores > 1L) {
    if (.Platform$OS.type != "windows") {
      have_mc <- TRUE
    } else {
      have_mc <- FALSE
    }
    loadNamespace("parallel") # get this out of the way before recording seed
  }

  set.seed(12345)

  dfm <- lapply(1:R, function(x) data[sample(1:nrow(data), nrow(data), replace = TRUE), ])

  if (ncores > 1L) {
    if (have_mc) {
      cl = ncores
    } else {
      cl <- parallel::makePSOCKcluster(rep("localhost", ncores))
      parallel::clusterSetRNGStream(cl)
      parallel::clusterExport(cl, ls(all.names = TRUE, env = environment()), envir = environment())
      parallel::clusterExport(cl, ls(all.names=TRUE, env=globalenv()), envir=globalenv())
      parallel::clusterExport(cl, c("monte_g", "simulate_data", "sim_value"))
      parallel::clusterEvalQ(cl, library(data.table))
      # clusterExport(cl, c("nobs", "iris_2species"))
      on.exit(parallel::stopCluster(cl))
    }
  } else {
    cl <- NULL
  }

  cat("\nBootstrapping:\n")

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
                               mediation_type = mediation_type,
                               cl = cl)

  return(boot_res)
}
