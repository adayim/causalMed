#' Calculate mediation analysis confidence interval
#'
#' @description Used to calculate confidence interval using non-parametric bootstrap methods.
#'
#' @inheritParams gformula
#' @param future_seed Logical or integer. Seed is passed to future_lapply.
#' @importFrom future.apply future_lapply
#' @importFrom progressr handler_progress handlers progressor
#' @keywords internal
#'
bootstrap_helper <- function(data,
                             id_var,
                             base_vars,
                             time_var,
                             exposure,
                             models,
                             intervention,
                             init_recode = NULL,
                             in_recode = NULL,
                             out_recode = NULL,
                             mc_sample = 10000,
                             mediation_type = c(NA, "N", "I"),
                             R = 500,
                             progress_bar = TRUE,
                             future_seed = TRUE) {
  mediation_type <- match.arg(mediation_type)

  # Progress bar
  if (progress_bar) {
    progressr::handlers(list(
      progressr::handler_progress(
        format   = "Bootstrap [:bar] :current/:total (:percent) | Elapsed: :elapsed | ETA: :eta",
        complete = "+"
      )
    ))
    p <- progressr::progressor(steps = R)
  }

  boot_res <- future.apply::future_lapply(1:R, function(i) {
    # Resample at the individual level (not row level) to preserve
    # within-individual longitudinal structure, matching gfoRmula.
    unique_ids <- unique(data[[id_var]])
    boot_ids <- sample(unique_ids, length(unique_ids), replace = TRUE)
    id_map <- data.table::data.table(orig_id = boot_ids,
                                     new_id = seq_along(boot_ids))
    data.table::setnames(id_map, "orig_id", id_var)
    boot_data <- merge(id_map, data, by = id_var, allow.cartesian = TRUE)
    boot_data[, (id_var) := new_id]
    boot_data[, new_id := NULL]

    res <- .gformula(data = boot_data,
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
                     seed = NULL)   # let each replicate draw its own MC sample
    
    if (progress_bar)
      p()

    return(res)

  }, future.seed = future_seed)

  return(boot_res)
}
