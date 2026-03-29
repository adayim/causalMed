
#' Simulate Data
#'
#' Loop through the models and apply any recoding or subset.
#'
#' @param data Data to be used for the data generation
#' @param models Model list passed from \code{\link{gformula}} or \code{\link{mediation}}.
#' @param intervention A vector, intervention treatment, natural course will be evaluated
#' if the value is \code{NULL}. If the value contains any logical operators, the intervention
#' will be evaluated and the exposure variable will be set to 1 if it is \code{TRUE}.
#' @param mediation_type Type of the mediation analysis, if the value is \code{NA}
#' no mediation analysis will be performed (default). It will be ignored if the intervention
#'  is not \code{NULL}
#' @param med_pool Optional named list keyed by time index (as character string). Each element
#'   is itself a list with two components: \code{vals} (simulated mediator values under a*) and
#'   \code{weights} (the corresponding survival probability \eqn{Sc} under a*, or \code{NULL}
#'   for non-survival analyses). When supplied and \code{mediation_type = "I"}, the mediator
#'   for the Phi10 (cross-world) arm is sampled from \code{vals} with probability proportional
#'   to \code{weights}, implementing the \eqn{S(1{:}t-1)=1} conditioning in Lin et al. (2017,
#'   Eq. 4).
#' @param med_ref_val The reference exposure value (a*) used when evaluating the mediator model
#'   for the Phi10 arm. Defaults to \code{0L}. Set to match the \code{ref_val} argument of
#'   \code{\link{mediation}}.
#'
#' @keywords internal
#'
simulate_data <- function(data,
                          exposure,
                          models,
                          intervention = NULL,
                          mediation_type = c(NA, "N", "I"),
                          med_pool = NULL,
                          med_ref_val = 0L) {

  # Replicate match.arg behaviour for a c(NA, "N", "I") default:
  # when the full default vector is passed (user did not specify), use first element (NA).
  if (length(mediation_type) > 1) {
    mediation_type <- mediation_type[1]
  } else if (!is.na(mediation_type) && !(mediation_type %in% c("N", "I"))) {
    stop("'mediation_type' must be NA, \"N\", or \"I\"")
  }

  is_phi10   <- !is.na(mediation_type) && is.null(intervention)
  is_dynamic <- inherits(intervention, "causalMed_dynint")

  # Ensure a valid data.table self-reference before any set() / := calls.
  data.table::setDT(data)

  # Set exposure to the intervention value for static interventions
  if (!is.null(intervention) && !is_dynamic) {
    set(data, j = exposure, value = intervention)
  }

  # Phi10 arm: fix exposure to exposure_val and cache the a*=ref_val swap expression
  if (is_phi10) {
    set(data, j = exposure, value = 1L)
    interv0 <- parse(text = paste0(exposure, " = ", med_ref_val))
  }

  # Loop through models
  for (indx in seq_along(models)) {
    model    <- models[[indx]]
    resp_var <- model$rsp_vars
    mod_type <- model$mod_type

    # Perform any model-specific recodes
    if (!is.null(model$recode)) {
      apply_recodes(data, model$recode)
    }

    # Evaluate subset condition
    if (!is.null(model$subset)) {
      cond <- eval(model$subset, envir = data, enclos = parent.frame())
      cond <- cond & !is.na(cond)
    } else {
      cond <- rep(TRUE, nrow(data))
    }

    # Skip the exposure model when:
    #   (a) a static or dynamic intervention is given (value already set above), or
    #   (b) this is the Phi10 arm (exposure already fixed to 1 above).
    if (resp_var == exposure && (!is.null(intervention) || is_phi10)) {
      if (is_dynamic) {
        # copy() so we don't modify the caller's data.table in-place;
        # copy() returns a fresh data.table with selfref = 1.
        data     <- data.table::copy(data)
        vals_nat <- sim_value(model = model, newdt = data[cond])
        data[cond, (exposure) := vals_nat]
        vals_dyn <- as.numeric(eval(intervention[[1]], envir = data[cond]))
        data[cond, (exposure) := vals_dyn]
      }
      next
    }

    if (sum(cond) != 0L) {

      if (mod_type == "mediator" && is_phi10) {
        # ── Phi10 cross-world mediator ──────────────────────────────────────────
        if (mediation_type == "N") {
          # Natural effects (Zheng et al., 2012, Eq. 5): evaluate mediator model
          # at a*=0 while keeping the individual's own covariate history under a=1.
          med_value <- sim_value(model = model,
                                 newdt = within(data[cond], eval(interv0)))
          data[cond, (resp_var) := med_value]

        } else {
          # Interventional effects (Lin et al., 2017, Eq. 4): draw from
          # a pre-collected mediator pool (correct l' trajectory).
          if (!is.null(med_pool)) {
            n_alive   <- sum(cond)
            pool_vals <- med_pool$vals
            pool_wts  <- med_pool$weights   # NULL for non-survival analyses
            pool_size <- length(pool_vals)
            if (pool_size >= n_alive) {
              if (pool_size < 10L) {
                warning(sprintf(
                  "Interventional mediation: mediator pool size (%d) is very small (< 10). Permutation may be unreliable.",
                  pool_size
                ))
              }
              # prob = pool_wts implements the S(1:t-1)=1 condition from
              # Lin et al. (2017, Eq. 4): individuals with lower cumulative
              # survival probability contribute less to the draw.
              # When pool_wts is NULL (no survival model), uniform sampling.
              data[cond, (resp_var) := sample(pool_vals, n_alive, replace = FALSE, prob = pool_wts)]
            } else if (pool_size > 0L) {
              warning(sprintf(
                "Interventional mediation: mediator pool size (%d) < alive count (%d). Sampling with replacement.",
                pool_size, n_alive
              ))
              data[cond, (resp_var) := sample(pool_vals, n_alive, replace = TRUE, prob = pool_wts)]
            }
            # If pool_size == 0, mediator is not updated (stale value retained).
          } else {
            # Fallback: draw from model at a*=0 then permute.
            med_value <- sim_value(model = model,
                                   newdt = within(data[cond], eval(interv0)))
            data[cond, (resp_var) := sample(med_value, length(med_value), replace = FALSE)]
          }
        }

      } else {
        # ── Standard model-based simulation ─────────────────────────────────────
        vals <- sim_value(model = model, newdt = data[cond])
        data[cond, (resp_var) := vals]
      }

      # ── Outcome prediction ───────────────────────────────────────────────────
      # Use pre-stored Xterms + beta (direct BLAS) instead of predict() to
      # avoid model.frame() reconstruction overhead on every time step.
      if (mod_type == "outcome") {
        mm   <- model.matrix(model$Xterms, data = data[cond])
        lp   <- drop(mm %*% model$beta)
        pred <- model$linkinv(lp)
        pred <- pmin(pmax(pred, 1e-10), 1 - 1e-10)
        data[cond, Pred_Y := pred]
      }

      if (mod_type == "survival") {
        mm <- model.matrix(model$Xterms, data = data[cond])
        lp <- drop(mm %*% model$beta)
        h  <- model$linkinv(lp)
        h  <- pmin(pmax(h, 1e-10), 1 - 1e-10)
        data[cond, S := h]

        # Initialise Sc to 1 on the first call (no Sc column yet), then
        # accumulate the product-limit survival across time steps.
        if (!("Sc" %in% names(data))) {
          data[, Sc := 1.0]
        }
        data[cond, Sc     := Sc * (1 - h)]
        data[cond, Pred_Y := 1 - Sc]
      }
    }
  }

  return(data)
}
