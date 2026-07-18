#' Calculate mediation analysis confidence interval
#'
#' @description Used to calculate confidence interval using non-parametric bootstrap methods.
#'
#' @inheritParams gformula
#' @param future_seed Logical or integer. Seed is passed to future_lapply.
#' @importFrom future.apply future_lapply
#' @importFrom future nbrOfWorkers
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
                             n_vw = 1L,
                             R = 500,
                             progress_bar = TRUE,
                             future_seed = TRUE) {
  mediation_type <- match.arg(mediation_type)

  # Progress bar. progressr::handlers() sets a SESSION-WIDE option, so we save
  # the caller's current handlers and restore them on exit — otherwise a
  # bootstrap run would silently replace the user's progressr configuration for
  # the rest of their session (affecting every other progressr-aware package).
  if (progress_bar) {
    old_handlers <- progressr::handlers()
    on.exit(progressr::handlers(old_handlers), add = TRUE)
    progressr::handlers(list(
      progressr::handler_progress(
        format   = "Bootstrap [:bar] :current/:total (:percent) | Elapsed: :elapsed | ETA: :eta",
        complete = "+"
      )
    ))
    p <- progressr::progressor(steps = R)
  }

  # Pin each worker to one data.table thread (~20% faster: W workers x ~cores/2
  # OpenMP threads each would oversubscribe the machine). Must be set inside the
  # worker, and only when parallel -- sequentially the "worker" is the caller's
  # own session. 
  is_parallel <- tryCatch(future::nbrOfWorkers() > 1L, error = function(e) FALSE)

  boot_res <- future.apply::future_lapply(1:R, function(i) {
    if (is_parallel) {
      old_dt <- data.table::getDTthreads()
      data.table::setDTthreads(1L)
      on.exit(data.table::setDTthreads(old_dt), add = TRUE)
    }

    # Resample at the individual level (not row level) to preserve
    # within-individual longitudinal structure, matching gfoRmula.
    unique_ids <- unique(data[[id_var]])
    boot_ids <- sample(unique_ids, length(unique_ids), replace = TRUE)
    id_map <- data.table::data.table(orig_id = boot_ids,
                                     new_id = seq_along(boot_ids))
    data.table::setnames(id_map, "orig_id", id_var)
    boot_data <- merge(id_map, data, by = id_var, allow.cartesian = TRUE)
    # Overwrite the id with the dense bootstrap index, then drop the helper
    # column. set() (not :=) keeps `new_id` a plain string, so it needs no
    # utils::globalVariables() entry.
    data.table::set(boot_data, j = id_var, value = boot_data[["new_id"]])
    data.table::set(boot_data, j = "new_id", value = NULL)

    res <- .run_interventions(data = boot_data,
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
                     n_vw = n_vw,
                     return_fitted = FALSE,
                     return_data = FALSE,
                     seed = NULL)   # let each replicate draw its own MC sample
    
    if (progress_bar)
      p()

    return(res)

  }, future.seed = future_seed)

  return(boot_res)
}
