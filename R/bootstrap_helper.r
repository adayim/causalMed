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
                             id_var,
                             base_vars,
                             time_var,
                             exposure,
                             models,
                             intervention,
                             init_recode,
                             in_recode,
                             out_recode,
                             mc_sample = 10000,
                             mediation_type = c(NA, "N", "I"),
                             R = 500,
                             ncores = 1L) {
  mediation_type <- match.arg(mediation_type)

  set.seed(12345)

  dfm <- lapply(1:R, function(x) data[sample(1:nrow(data), nrow(data), replace = TRUE), ])

  if (ncores > 1L) {
    if (have_mc) {
      cl <- ncores
    } else {
      cl <- parallel::makePSOCKcluster(rep("localhost", ncores))
      parallel::clusterSetRNGStream(cl)
      parallel::clusterExport(cl, ls(all.names = TRUE, env = environment()), envir = environment())
      parallel::clusterExport(cl, ls(all.names = TRUE, env = globalenv()), envir = globalenv())
      parallel::clusterExport(cl, c("monte_g", "simulate_data", "sim_value"))
      parallel::clusterEvalQ(cl, library(data.table))

      on.exit(parallel::stopCluster(cl))
    }
  } else {
    cl <- NULL
  }

  cat("\nBootstrapping:\n")

  boot_res <- pbapply::pblapply(dfm, .gformula,
    id_var = id_var,
    base_vars = base_vars,
    time_var = time_var,
    exposure = exposure,
    models = models,
    intervention = intervention,
    in_recode = in_recode,
    out_recode = out_recode,
    init_recode = init_recode,
    mc_sample = mc_sample,
    mediation_type = mediation_type,
    cl = cl
  )

  return(boot_res)
}
