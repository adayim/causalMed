#' Calculate mediation analysis confidence interval
#'
#' @description Used to calculate confidence interval using non-parametric bootstrap methods.
#'
#' @param object Object returned from gformula function (see \code{\link[calsalMed]{gformula}}).
#'
#' @param R The number of bootstrap replicates, default is 500. Same with \code{boot}, see \code{\link[boot]{boot}} for detail.
#'
#' @importFrom future.apply future_lapply
#' @importFrom progressr handler_progress handlers progressor
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
                             R = 500) {
  mediation_type <- match.arg(mediation_type)

  # Progress bar
  progressr::handlers(global = TRUE)
  progressr::handlers(list(
    progressr::handler_progress(
      format   = ":spin :current/:total (:message) [:bar] :percent in :elapsedfull ETA: :eta",
      # width    = 60,
      complete = "+"
    )
  ))
  p <- progressr::progressor(steps = R)

  boot_res <- future.apply::future_lapply(1:R, function(i) {
    set.seed(12345 * i)
    indx <- sample(1:nrow(data), nrow(data), replace = TRUE)

    p(message = "Bootstrapping", amount = 0)

    res <- .gformula(data = data[indx, ],
                     id_var = id_var,
                     base_vars = base_vars,
                     time_var = time_var,
                     exposure = exposure,
                     models = models,
                     intervention = intervention,
                     init_recode = init_recode,
                     in_recode = in_recode,
                     out_recode = out_recode,
                     mc_sample = mc_sample,
                     mediation_type = mediation_type,
                     return_fitted = FALSE,
                     return_data = FALSE,
                     progress_bar = FALSE)

    p(message = "Bootstrapping")

    return(res)

  }, future.seed = TRUE)

  return(boot_res)
}
