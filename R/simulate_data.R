
#' Simulate Data
#'
#' Loop through the models and apply any recoding or subset.
#'
#' @param data Data to be used for the data generation
#' @param models Model list passed from \code{\link{gformula}} or \code{\link{mediation}}.
#' @param intervention One of:
#'   \itemize{
#'     \item \code{NULL} — natural-course draw of the exposure (gformula) or
#'           Phi10 cross-world arm (mediation, when \code{mediation_type} is
#'           non-\code{NA}). The legacy NULL = Phi10 path is kept for
#'           backwards compatibility but \code{mediation()} now constructs
#'           \code{causalMed_arm} objects internally.
#'     \item Numeric scalar/vector or \code{causalMed_dynint} — static or
#'           dynamic exposure rule for \code{gformula()}.
#'     \item \code{causalMed_arm} — structured mediation arm carrying a
#'           fixed treatment value plus an optional named list of
#'           per-mediator overrides. See \code{arm_spec()}.
#'   }
#' @param mediation_type Type of the mediation analysis, if the value is \code{NA}
#' no mediation analysis will be performed (default).
#' @param med_pool Optional named list keyed by mediator response variable.
#'   Each element is the time-\eqn{t} slice of a pre-permuted joint
#'   mediator-trajectory matrix built by \code{.gformula} from the
#'   corresponding reference arm (Phi00 or Phi11). When the arm
#'   \code{intervention} requires a mediator override under
#'   \code{mediation_type = "I"}, the mediator is assigned directly from this
#'   vector (joint draw matching Lin et al. 2017 Eq. 4, the SAS mGFORMULA
#'   macro, and Yamamuro et al. 2021 Figure 3 step 3).
#'
#' @keywords internal
#'
simulate_data <- function(data,
                          exposure,
                          models,
                          intervention = NULL,
                          mediation_type = c(NA, "N", "I"),
                          med_pool = NULL) {

  # Replicate match.arg behaviour for a c(NA, "N", "I") default:
  # when the full default vector is passed (user did not specify), use first element (NA).
  if (length(mediation_type) > 1) {
    mediation_type <- mediation_type[1]
  } else if (!is.na(mediation_type) && !(mediation_type %in% c("N", "I"))) {
    stop("'mediation_type' must be NA, \"N\", or \"I\"")
  }

  is_arm     <- inherits(intervention, "causalMed_arm")
  is_phi10   <- !is.na(mediation_type) && is.null(intervention)
  is_dynamic <- inherits(intervention, "causalMed_dynint")

  # Ensure a valid data.table self-reference before any set() / := calls.
  data.table::setDT(data)

  # Set exposure for every flavour that fixes treatment.
  if (is_arm) {
    set(data, j = exposure, value = intervention$treatment)
  } else if (!is.null(intervention) && !is_dynamic) {
    set(data, j = exposure, value = intervention)
  } else if (is_phi10) {
    # Legacy Phi10 path (no arm_spec passed): fix exposure to 1.
    set(data, j = exposure, value = 1L)
  }

  # Per-mediator metadata for an arm-driven simulation. For "N" arms, the
  # override value is the swap-target exposure used when re-evaluating the
  # mediator model; for "I" arms the corresponding pool slice is read from
  # med_pool by mediator name.
  med_overrides <- if (is_arm) intervention$mediator_overrides else NULL

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

    # Skip the exposure model when treatment is already fixed:
    #   (a) static / dynamic intervention (gformula path),
    #   (b) arm_spec path (mediation),
    #   (c) legacy Phi10 path.
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

      # ── Mediator handling under a cross-world override ──────────────────────
      # Three triggers for cross-world handling:
      #   (1) arm_spec with this mediator in mediator_overrides,
      #   (2) legacy single-mediator Phi10 path (is_phi10),
      # Otherwise the mediator is simulated normally from its fitted model.
      mediator_override_value <- if (is_arm && mod_type == "mediator" &&
                                     resp_var %in% names(med_overrides)) {
        med_overrides[[resp_var]]
      } else {
        NULL
      }
      has_arm_override <- !is.null(mediator_override_value)

      if (mod_type == "mediator" && (has_arm_override || is_phi10)) {
        if (mediation_type == "N") {
          # Natural effects (Zheng et al., 2012, Eq. 5): evaluate the
          # mediator model on the arm's own covariate history but with the
          # exposure swapped to the cross-world value.
          swap_val <- if (has_arm_override) mediator_override_value else 0L
          interv_swap <- parse(text = paste0(exposure, " = ", swap_val))
          med_value <- sim_value(model = model,
                                 newdt = within(data[cond], eval(interv_swap)))
          data[cond, (resp_var) := med_value]

        } else {
          # Interventional effects (Lin et al. 2017 Eq. 4; Yamamuro et al.
          # 2021 Fig. 3 step 3): direct assignment from the pre-permuted
          # joint-trajectory pool slice for this mediator. arm_spec carries
          # the pool source (0 or 1) but the actual pool slice was looked up
          # by mediator name in .gformula and is delivered via med_pool.
          this_pool <- if (!is.null(med_pool)) med_pool[[resp_var]] else NULL
          if (!is.null(this_pool)) {
            data[cond, (resp_var) := this_pool[cond]]
          } else {
            # Fallback (no pool collected, e.g., reference arm had no
            # mediator model): draw from the mediator model at the
            # cross-world exposure and permute within this time step.
            swap_val <- if (has_arm_override) mediator_override_value else 0L
            interv_swap <- parse(text = paste0(exposure, " = ", swap_val))
            med_value <- sim_value(model = model,
                                   newdt = within(data[cond], eval(interv_swap)))
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
