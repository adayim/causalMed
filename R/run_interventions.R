#' Main calculation function
#'
#' This function will receive the parameters and fit model. After the model is fitted,
#' random samples will be drawn from the data and apply the intervention.
#'
#' @inheritParams gformula
#' @param mediation_type Type of the mediation analysis, if the value is \code{NA}
#' no mediation analysis will be performed (default). It will be ignored if the intervention
#'  is not \code{NULL}
#' @param n_vw Integer. Number of independent permutation draws averaged for
#'   interventional pool-drawing interventions (Vansteelandt-Williamson repetition).
#'   Reference interventions (no mediator overrides) and natural-effect interventions are
#'   unaffected. Default \code{1L}; \code{mediation()} sets this to 2 to
#'   match the SAS mGFORMULA macro.
#' @param return_fitted Return the fitted model (default is FALSE).
#' @param time_seq Numeric vector: the sorted distinct time points to simulate,
#'   as computed by \code{gformula()}/\code{mediation()} from the \emph{input}
#'   data. Passed to every bootstrap replicate so all passes simulate the same
#'   steps even when a resample lacks a rarely observed time value. When
#'   \code{NULL} the grid is derived from \code{data}.
#' @keywords internal

.run_interventions <- function(data,
                      id_var,
                      base_vars,
                      time_var,
                      exposure,
                      models,
                      intervention,
                      in_recode = NULL,
                      out_recode = NULL,
                      init_recode = NULL,
                      mediation_type = c(NA, "N", "I"),
                      mc_sample = 10000,
                      n_vw = 1L,
                      return_fitted = FALSE,
                      return_data = FALSE,
                      seed = NULL,
                      time_seq = NULL) {
  if (length(mediation_type) > 1) {
    mediation_type <- mediation_type[1]
  } else if (!is.na(mediation_type) && !(mediation_type %in% c("N", "I"))) {
    stop("'mediation_type' must be NA, \"N\", or \"I\"")
  }

  fit_mods <- fit_spec_models(models, data)

  # only call set.seed() when a seed is explicitly supplied.
  # Omitting the seed inside bootstrap replicates lets each replicate draw
  # different Monte Carlo samples, giving valid bootstrap variance estimates.
  if (!is.null(seed)) set.seed(seed)

  # Baseline variables for Monte Carlo Sampling
  base_dat <- unique(data[, unique(c(id_var, base_vars)), with = FALSE])
  df_mc    <- data.table::as.data.table(
    base_dat[sample(1:length(base_dat[[id_var]]), mc_sample, replace = TRUE), ]
  )

  # The grid comes from the caller (computed once from the INPUT data) so that
  # a bootstrap resample missing a rarely observed time value still simulates
  # every step; a resample missing an interior time value shortens the grid,
  # so `intervention[indx]` and the per-step regime slice would index the
  # wrong step. The fallback serves direct internal callers.
  if (is.null(time_seq)) time_seq <- time_grid(data, time_var)

  # ── Pool collection for interventional pool-drawing interventions ──────────────────
  # For each distinct regime appearing in any intervention's mediator_overrides,
  # run a reference intervention under that regime once with collect_pool = TRUE.
  # The resulting `pools` object is a list keyed by regime_key(); each
  # element is itself a named list keyed by mediator response variable, whose
  # value is a length-T list of `mc_sample`-long vectors — the k-th holding
  # every pool individual's simulated mediator value at `time_seq[k]`, in the
  # mediator's own type. Pool individual i's trajectory is the i-th entry of
  # every slice.
  #
  # `cached_arms` stores the scalar Phi estimate from each pool-collection
  # pass so that pure reference interventions (no mediator overrides) in the main loop
  # below can be served from the cache without a redundant simulation.
  pools       <- list()
  cached_arms <- list()

  has_mediation_spec <- any(sapply(intervention, inherits, "causalMed_intervention"))

  if (isTRUE(mediation_type == "I") && has_mediation_spec) {
    # The distinct override REGIMES across all interventions, keyed by
    # regime_key(), in order of first appearance -- first appearance across
    # interventions, then across each intervention's overrides -- because
    # pool-collection order determines RNG consumption, so it must stay
    # stable. No NA filter is needed: intervention_spec() rejects NA in
    # `treatment` and in every override, so an NA regime cannot reach here.
    pool_regimes <- list()
    for (interv in intervention) {
      if (!inherits(interv, "causalMed_intervention")) next
      for (v in interv$mediator_overrides) {
        key <- regime_key(v)
        if (is.null(pool_regimes[[key]])) pool_regimes[[key]] <- v
      }
    }

    for (key in names(pool_regimes)) {
      lvl <- pool_regimes[[key]]
      ref_interv  <- intervention_spec(treatment = lvl, mediator_overrides = list())
      ref_data <- data.table::copy(df_mc)
      ref_run  <- simulate_intervention(
        data           = ref_data,
        models         = fit_mods,
        exposure       = exposure,
        time_var       = time_var,
        time_seq       = time_seq,
        intervention   = ref_interv,
        init_recode    = init_recode,
        in_recode      = in_recode,
        out_recode     = out_recode,
        mediation_type = mediation_type,
        return_data    = return_data,
        med_pool       = NULL,
        collect_pool   = TRUE
      )
      pools[[key]]       <- ref_run$pool
      cached_arms[[key]] <- ref_run$estimate
    }
  }

  # ── Per-intervention dispatch ─────────────────────────────────────────────────────
  res <- sapply(intervention, function(interv) {

    # ----- intervention_spec path (mediation) -----
    if (inherits(interv, "causalMed_intervention")) {
      lvl_key <- regime_key(interv$treatment)

      # Reference intervention (no mediator overrides) — use cached pool-collection
      # result if available; otherwise run once.
      if (length(interv$mediator_overrides) == 0L) {
        if (!is.null(cached_arms[[lvl_key]])) return(cached_arms[[lvl_key]])
        interv_data <- data.table::copy(df_mc)
        return(simulate_intervention(
          data           = interv_data,
          models         = fit_mods,
          exposure       = exposure,
          time_var       = time_var,
          time_seq       = time_seq,
          intervention   = interv,
          init_recode    = init_recode,
          in_recode      = in_recode,
          out_recode     = out_recode,
          mediation_type = mediation_type,
          return_data    = return_data,
          med_pool       = NULL
        ))
      }

      # Intervention with mediator overrides.
      # "I": n_vw permutations averaged. "N": single pass (no permutation).
      n_reps <- if (isTRUE(mediation_type == "I")) max(1L, n_vw) else 1L

      # When n_reps > 1 we cannot meaningfully return n_reps data tables;
      # we keep the last replicate's data for inspection and average the
      # scalar Phi across replicates.
      vw_scalars <- numeric(n_reps)
      last_data  <- NULL

      for (rep_i in seq_len(n_reps)) {
        interv_data <- data.table::copy(df_mc)

        # Build per-mediator pre-permuted pool matrix for this replicate.
        # Each mediator's pool is permuted independently (Yamamuro 2021
        # Eq. 2: the pool draws are independent across mediators).
        interv_med_pool <- NULL
        if (isTRUE(mediation_type == "I")) {
          interv_med_pool <- list()
          for (med_var in names(interv$mediator_overrides)) {
            src_key      <- regime_key(interv$mediator_overrides[[med_var]])
            pool_for_med <- pools[[src_key]][[med_var]]
            # A key miss is a bug, not a fallback case: anything other than a
            # joint trajectory -- e.g. permuting the mediator WITHIN each time
            # step -- would be a different estimand, produced with no error
            # (the same failure class as the collapsed-pool regression). Fail
            # loudly here, naming the collected keys; simulate_data() also
            # stops if a slice is ever missing.
            if (is.null(pool_for_med)) {
              stop(sprintf(paste0(
                "Internal error: no mediator pool for '%s' under regime key '%s'. ",
                "Pool keys collected: {%s}. Regimes must be recycled to full ",
                "length before interventions are built."),
                med_var, src_key, paste(names(pools), collapse = "; ")),
                domain = "causalMed")
            }
            # ONE permutation applied to every time slice, so subject i
            # receives pool individual perm[i]'s whole trajectory -- a joint
            # draw of M(1:T), not an independent value per time step.
            perm <- sample.int(length(pool_for_med[[1L]]))
            interv_med_pool[[med_var]] <-
              lapply(pool_for_med, function(v) v[perm])
          }
        }

        # On the LAST replicate, ask for the simulated data when the caller
        # wants it and derive the scalar from that same run. Asking for the
        # data and then re-running would consume a second, independent RNG
        # stream, so `mean(Pred_Y)` of the returned table would not equal the
        # Phi it is reported beside -- true even at n_vw = 1.
        want_data <- return_data && rep_i == n_reps
        r <- simulate_intervention(
          data           = interv_data,
          models         = fit_mods,
          exposure       = exposure,
          time_var       = time_var,
          time_seq       = time_seq,
          intervention   = interv,
          init_recode    = init_recode,
          in_recode      = in_recode,
          out_recode     = out_recode,
          mediation_type = mediation_type,
          return_data    = want_data,
          med_pool       = interv_med_pool
        )
        if (want_data) {
          last_data <- r
          # Same definition simulate_intervention() uses for its scalar return.
          vw_scalars[rep_i] <- sum(r[["Pred_Y"]]) / length(r[["Pred_Y"]])
        } else {
          vw_scalars[rep_i] <- r
        }
      }

      avg_phi <- mean(vw_scalars)
      if (return_data) {
        # Attach the averaged Phi as an attribute so downstream code can use
        # it; the table itself reflects only the last permutation.
        data.table::setattr(last_data, "phi_estimate", avg_phi)
        return(last_data)
      }
      return(avg_phi)
    }

    # ----- legacy scalar / NULL / dyn_int path (gformula) -----
    interv_data <- data.table::copy(df_mc)
    simulate_intervention(
      data           = interv_data,
      models         = fit_mods,
      exposure       = exposure,
      time_var       = time_var,
      time_seq       = time_seq,
      intervention   = interv,
      init_recode    = init_recode,
      in_recode      = in_recode,
      out_recode     = out_recode,
      mediation_type = mediation_type,
      return_data    = return_data,
      med_pool       = NULL
    )
  }, simplify = FALSE)

  if (return_fitted) {
    return(list(fitted.models = fit_mods, gform.data = res))
  } else {
    return(list(gform.data = res))
  }
}

# Fit every spec_model() in `models` on `data` and pre-extract the
# prediction components used by both the Monte Carlo engine (sim_value,
# simulate_data) and the TMLE engine (density evaluation). Extracted from
# .run_interventions() so the two estimators share one fitting path.
fit_spec_models <- function(models, data) {
  lapply(models, function(mods) {
    rsp_vars <- all.vars(formula(mods$call)[[2]])

    # Observed values range
    if (is.numeric(data[[rsp_vars]])) {
      val_ran <- range(na.omit(data[[rsp_vars]]))
    } else {
      val_ran <- unique(na.omit(data[[rsp_vars]]))
    }

    # Recode data before fitting this model. `recode` is a causalMed_recodes
    # list of expressions, so it must go through apply_recodes() (the same
    # path used at simulation time). Applied to a copy so per-model recodes
    # do not leak across models or back to the caller's data.
    if (!is.null(mods$recode)) {
      data <- apply_recodes(data.table::copy(data), mods$recode)
    }

    mods$call$data <- substitute(data, env = parent.frame())

    fitmodel <- run_withwarning_collect(
      eval(mods$call),
      msg = sprintf("Outcome model: %s", rsp_vars)
    )

    # Aliased (NA) coefficients arise whenever the design is rank deficient --
    # a perfectly collinear term, a factor level unobserved in the fitting
    # subset, or a covariate that is constant among the complete cases (e.g.
    # `time` in an outcome recorded only at the end of follow-up). `beta` would
    # then carry an NA and `mm %*% beta` would return NA for EVERY row, turning
    # the whole analysis into a silent all-NA result. predict() drops aliased
    # terms instead, so do the same here and say which terms went: keep only the
    # estimable coefficients and record their names so lin_pred() can subset the
    # design matrix to match. Categorical (multinom) fits are exempt -- their
    # coef() is a matrix and they predict through predict(), not `beta`.
    #
    # Both terms() and coef() are extracted DEFENSIVELY. A `custom_fit` model
    # can be any class, and a data-adaptive learner (ranger, xgboost, a
    # SuperLearner ensemble) typically has neither a terms component nor a
    # coef() method. That is fine as long as nothing evaluates a linear
    # predictor for this node — see `needs_lp` below — because `custom_sim`
    # then does the drawing. Letting terms()/coef() error here would refuse a
    # perfectly usable nuisance model; letting them return NULL silently would
    # surface later as "requires numeric/complex matrix/vector arguments" out
    # of `mm %*% NULL`.
    Xterms    <- tryCatch(delete.response(terms(fitmodel)), error = function(e) NULL)
    beta      <- tryCatch(coef(fitmodel), error = function(e) NULL)
    beta_cols <- NULL
    if (!is.null(beta) && !is.matrix(beta) && anyNA(beta)) {
      # Index by POSITION when coef() returns an unnamed vector -- legal, and
      # what a hand-rolled or wrapped fit may hand back. Subsetting by name
      # then silently reduced `beta` to numeric(0) (NULL[logical] is NULL),
      # reported an empty `{}` set of dropped terms, and refused a perfectly
      # usable model with "provides no usable model terms/coefficients".
      # `lin_pred()` subsets the design matrix with `mm[, beta_cols]`, which
      # takes column indices just as happily as names.
      keep      <- !is.na(beta)
      nms       <- names(beta)
      aliased   <- if (is.null(nms)) paste("position", which(!keep)) else nms[!keep]
      beta_cols <- if (is.null(nms)) which(keep) else nms[keep]
      beta      <- beta[keep]
      causalmed_env$warning <- c(
        causalmed_env$warning,
        sprintf(paste0("Model for '%s' is rank deficient: term(s) {%s} could not ",
                       "be estimated and were dropped from the simulation (as ",
                       "predict() would). This usually means a predictor is ",
                       "collinear with another, or is constant among the rows ",
                       "used to fit the model -- check the model formula."),
                rsp_vars, paste(aliased, collapse = ", "))
      )
    }

    # Does the simulation evaluate a linear predictor for this node?
    #   - outcome/survival: always — simulate_data() computes the predicted
    #     risk/hazard from the fitted coefficients, and a custom_sim function
    #     cannot stand in for that (it supplies a DRAW, not the estimand).
    #   - other nodes: only when sim_value() falls through to lin_pred(), i.e.
    #     when there is no custom_sim and the type is not "categorical" (which
    #     predicts through predict()).
    # Refuse early, and say what to do, rather than failing inside a matrix
    # product several layers down.
    needs_lp <- mods$mod_type %in% c("outcome", "survival") ||
      (!identical(mods$var_type, "categorical") && is.null(mods$custom_sim))
    if (needs_lp && (is.null(Xterms) || is.null(beta) || length(beta) == 0L)) {
      stop(sprintf(paste0(
        "The fitted model for '%s' (mod_type = \"%s\") provides no usable model ",
        "terms/coefficients, but the simulation needs its linear predictor. %s"),
        rsp_vars, mods$mod_type,
        if (mods$mod_type %in% c("outcome", "survival"))
          paste0("The outcome/survival model supplies the predicted risk itself, ",
                 "so it must be a model with a terms component and a coef() ",
                 "method (e.g. glm). A custom_sim function cannot replace it.")
        else
          paste0("Supply a custom_sim function in spec_model() to draw '",
                 rsp_vars, "', or fit it with a model that has a coef() method.")),
        call. = FALSE, domain = "causalMed")
    }

    # sim_value() draws a node that has no custom_sim and var_type "normal" or
    # "custom" as lp + N(0, sigma), i.e. centred on the LINEAR PREDICTOR. That
    # is the fitted mean only under an identity link; with any other link (a
    # Poisson glm, say) the draws would land on the link scale without any
    # error. Refuse instead. Outcome/survival nodes are exempt: their risk is
    # computed through the inverse link, not drawn.
    draws_normal <- is.null(mods$custom_sim) &&
      mods$var_type %in% c("normal", "custom") &&
      !mods$mod_type %in% c("outcome", "survival")
    fit_link <- tryCatch(stats::family(fitmodel)$link, error = function(e) NULL)
    if (draws_normal && is.character(fit_link) && !identical(fit_link, "identity")) {
      stop(sprintf(paste0(
        "The model for '%s' (var_type = \"%s\") was fitted with a '%s' link, but ",
        "without custom_sim its values are simulated as normal draws around the ",
        "linear predictor, which is the fitted mean only under an identity link. ",
        "Supply a custom_sim function in spec_model() that draws from the ",
        "distribution the fitted model implies."),
        rsp_vars, mods$var_type, fit_link), call. = FALSE, domain = "causalMed")
    }

    # Pre-extract prediction components so the hot-path simulation loop can
    # skip predict() overhead (model.frame construction + na.action).
    # model.matrix(Xterms, newdt) %*% beta is a direct BLAS call.
    list(
      fitted     = fitmodel,
      Xterms     = Xterms,
      beta       = beta,
      # NULL for a full-rank fit, so the hot path pays nothing; otherwise the
      # estimable coefficients' names, or their column positions when coef()
      # came back unnamed.
      beta_cols  = beta_cols,
      linkinv    = if (inherits(fitmodel, "glm")) family(fitmodel)$linkinv
                   else identity,
      # suppressWarnings: sigma.default probes nobs(), which emits a
      # "no 'nobs' method is available" warning for non-standard fit classes.
      # That is noise about a value we are only opportunistically extracting.
      sigma      = suppressWarnings(
                     tryCatch(sigma(fitmodel), error = function(e) NULL)),
      recodes    = mods$recode,
      subset     = mods$subset,
      var_type   = mods$var_type,
      mod_type   = mods$mod_type,
      custom_sim = mods$custom_sim,
      # NULL for model objects built before spec_model() gained `truncate`;
      # sim_value() treats NULL as TRUE (the historical behaviour).
      truncate   = mods$truncate,
      rsp_vars   = rsp_vars,
      val_ran    = val_ran
    )
  })
}

# Linear predictor from the pre-extracted design terms and coefficients.
# Shared by sim_value() (drawing), simulate_data() (outcome/hazard prediction)
# and the TMLE engine's dens_value(), so all three treat a rank-deficient fit
# the same way. `beta_cols` is NULL for the usual full-rank case, in which case
# this is exactly the original `model.matrix() %*% beta` BLAS call.
lin_pred <- function(model, newdt) {
  mm <- model.matrix(model$Xterms, data = newdt)
  if (!is.null(model$beta_cols)) mm <- mm[, model$beta_cols, drop = FALSE]
  drop(mm %*% model$beta)
}

#' Monte Carlo simulation
#'
#' @description
#'  Internal use only. Monte Carlo simulation.
#'
#' @inheritParams gformula
#' @param time_seq Time sequence vector of the data.
#' @param med_pool Optional named list keyed by mediator response variable.
#'   Each element is itself a pre-permuted list of \code{T} vectors (one per
#'   time point, each of length \code{nrow(data)}, in the mediator's own type)
#'   supplying the cross-regime joint trajectory for that mediator. The per-time
#'   slice is delivered to \code{simulate_data}.
#' @param collect_pool Logical. If \code{TRUE}, capture the simulated
#'   trajectory of every mediator into a named list of per-time vector lists
#'   and return it alongside the risk estimate.
#'
#' @keywords internal
#'
simulate_intervention <- function(data,
                    models,
                    exposure,
                    time_var,
                    time_seq,
                    intervention = NULL,
                    init_recode  = NULL,
                    in_recode    = NULL,
                    out_recode   = NULL,
                    mediation_type = c(NA, "N", "I"),
                    return_data  = FALSE,
                    med_pool     = NULL,
                    collect_pool = FALSE) {

  if (length(mediation_type) > 1) {
    mediation_type <- mediation_type[1]
  } else if (!is.na(mediation_type) && !(mediation_type %in% c("N", "I"))) {
    stop("'mediation_type' must be NA, \"N\", or \"I\"")
  }

  # Replicate static gformula() interventions to the length of the time
  # sequence. dyn_int() objects are not replicated (the same rule applies at
  # every step); causalMed_intervention objects carry their own regimes and are
  # sliced per step below.
  time_len <- length(time_seq)
  # check_intervention() accepts logical static interventions; coerce them to
  # numeric so replication and exposure assignment treat them like {0, 1}.
  if (is.logical(intervention)) {
    intervention <- as.numeric(intervention)
  }
  if (is.numeric(intervention) && length(intervention) == 1L) {
    intervention <- rep(intervention, time_len)
  }

  # Pre-detect mediator models for pool collection.
  if (collect_pool) {
    med_flag_pc <- sapply(models, function(m) m$mod_type == "mediator")
    if (!any(med_flag_pc)) {
      collect_pool <- FALSE   # no mediator model — nothing to collect
    } else {
      med_var_pc <- vapply(models[which(med_flag_pc)],
                           function(m) m$rsp_vars, character(1))
      # One pool per mediator, stored as a LIST of per-time vectors rather than
      # a numeric matrix. A matrix forces a single storage mode, which silently
      # mangled non-numeric mediators: a factor was reduced to its integer level
      # codes and a character coerced the whole matrix. The list keeps each
      # time slice in the mediator's own type (numeric, character or factor
      # with its levels), so what is drawn from the pool is what was simulated.
      # Element k holds the pool individuals' values at time_seq[k]; the joint
      # trajectory of pool individual i is the i-th entry of every element.
      m_pool_out <- setNames(
        lapply(med_var_pc, function(v) {
          setNames(vector("list", length(time_seq)), as.character(time_seq))
        }),
        med_var_pc
      )
    }
  }

  # Get the position of the outcome
  out_flag <- sapply(models, function(mods) mods$mod_type %in% c("outcome", "survival"))
  out_flag <- which(out_flag)

  # Get the position of the censor
  cen_flag  <- sapply(models, function(mods) mods$mod_type == "censor")
  surv_flag <- sapply(models, function(mods) mods$mod_type == "survival")
  is_survival <- any(cen_flag) || any(surv_flag)

  if (any(cen_flag))
    cen_flag <- which(cen_flag)
  else
    cen_flag <- 0

  # Get the variable name of outcome and censor
  outcome <- all.vars(formula(models[[out_flag]]$fitted)[[2]])
  censor  <- NULL
  if (!all(cen_flag == 0)) {
    censor <- all.vars(formula(models[[cen_flag]]$fitted)[[2]])
  }

  # Get minimum and maximum time
  max_time <- max(time_seq, na.rm = TRUE)
  min_time <- min(time_seq, na.rm = TRUE)

  is_intervention_spec <- inherits(intervention, "causalMed_intervention")

  # Natural effects draw the cross-world mediator on the exposure HISTORY set
  # to the other regime (Zheng & van der Laan 2017, Eq. 5). in_recode has
  # copied THIS intervention's step-(k-1) exposure into the first-order lag
  # columns by the time the mediator is drawn at step k, so simulate_data() is
  # handed the other regime's step-(k-1) value for them. At the first step
  # in_recode does not run and the lags hold their regime-free init value.
  swap_lag_cols <- if (is_intervention_spec && isTRUE(mediation_type == "N")) {
    exposure_lag_cols(in_recode, exposure)
  } else {
    character(0)
  }

  # Run g-formula
  for (indx in seq_along(time_seq)) {
    t_index <- time_seq[indx]

    set(data, j = time_var, value = t_index)

    # Re-code baseline variables at initiation
    if (t_index == min_time) {
      if (!is.null(init_recode)) {
        apply_recodes(data, init_recode)
      }
    }

    # Re-code data before simulating
    if (!is.null(in_recode) & t_index != min_time) {
      apply_recodes(data, in_recode)
    }

    # Per-time slice of each mediator pool (named list of length-nrow(data)
    # vectors, already in subject-i order because the pool was pre-permuted
    # in .run_interventions).
    t_med_pool <- if (!is.null(med_pool)) {
      setNames(lapply(med_pool, function(p) p[[indx]]), names(med_pool))
    } else NULL

    # A causalMed_intervention carries whole regimes; hand simulate_data() the
    # scalar slice for this step. dyn_int() objects pass through as-is; static
    # gformula() vectors are indexed to the time step.
    current_int <- if (is_intervention_spec) {
      slice_intervention_spec(intervention, indx)
    } else if (inherits(intervention, "causalMed_dynint")) {
      intervention
    } else {
      intervention[indx]
    }

    med_swap_lags <- if (length(swap_lag_cols) > 0L && indx > 1L) {
      lapply(slice_intervention_spec(intervention, indx - 1L)$mediator_overrides,
             function(v) setNames(rep(list(v), length(swap_lag_cols)), swap_lag_cols))
    } else {
      NULL
    }

    data <- simulate_data(
      data           = data,
      exposure       = exposure,
      models         = models,
      intervention   = current_int,
      mediation_type = mediation_type,
      med_pool       = t_med_pool,
      med_swap_lags  = med_swap_lags
    )

    # Collect each mediator's value into the pool (used by the reference
    # intervention pass). Rows where the mediator model's subset would not
    # have applied at time t carry whatever value was there at the start of
    # this iteration — preserving the per-individual trajectory.
    #
    # copy() IS LOAD-BEARING. `data[[mv]]` hands back the data.table's own
    # column vector, and simulate_data() writes the next time step with
    # `data[cond, (resp_var) := vals]` — a `:=` WITH an `i` argument, which
    # sub-assigns IN PLACE. Storing the bare vector therefore stores a live
    # reference: every slice mutates together and the whole pool ends up
    # holding M(T), collapsing the joint trajectory to a constant. (The
    # earlier matrix form was safe only incidentally: `mat[, indx] <- x`
    # copies the values in.) Snapshot it explicitly. 
    #
    # The NULL guard is equally load-bearing. simulate_data() skips a model
    # WHOLESALE when its `subset` selects no rows at this time step
    # (`if (sum(cond) != 0L)`), so the mediator column may not exist yet and
    # `data[[mv]]` is then NULL. `[[<-` with a NULL value DELETES the list
    # element instead of storing it: the pool silently loses a slice, every
    # later time step shifts down one position, and the run eventually dies in
    # sim_value() with "Argument eta must be a nonempty numeric vector". Store
    # an explicit all-missing slice so the pool stays aligned with `time_seq`
    # -- nothing consumes it, because the same empty `subset` also suppresses
    # the pool-draw assignment at that step.
    if (collect_pool) {
      for (mv in med_var_pc) {
        slice <- data[[mv]]
        m_pool_out[[mv]][[indx]] <- if (is.null(slice)) {
          rep(NA, nrow(data))
        } else {
          data.table::copy(slice)
        }
      }
    }

    # For survival outcome: under intervention, disable censoring. All
    # individuals remain in the pool; risk is computed analytically from
    # Pred_Y = 1 - prod(1 - h_t) accumulated in simulate_data().
    if (is_survival) {
      if (!is.null(intervention) && cen_flag != 0) {
        set(data, j = censor, value = 0)
      }
      if (cen_flag != 0) {
        # set(), not `[[<-`: base list assignment shallow-copies the
        # data.table and invalidates its internal self-reference, so the very
        # next `:=` (out_recode, or the next time step's recodes) emitted
        # data.table's "shallow copy was taken" warning to the user, once per
        # time step.
        set(data, j = outcome,
            value = ifelse(data[[censor]] == 1, 0, data[[outcome]]))
      }
    }

    # Recode data after simulating
    if (!is.null(out_recode) & t_index != min_time) {
      apply_recodes(data, out_recode)
    }
  }
  # loop ends here

  # NA predicted outcomes silently poison the intervention mean (which is a
  # plain sum/length, so a single NA makes the estimate NA). NAs get here from
  # NA predictors -- e.g. a lag column never initialised by init_recode, or a
  # recode that produced NA -- so surface the count rather than returning a
  # silently NA (or silently partial) estimate.
  n_na <- sum(is.na(data[["Pred_Y"]]))
  if (n_na > 0L) {
    causalmed_env$warning <- c(
      causalmed_env$warning,
      sprintf(paste0("Simulation produced %d NA predicted outcome(s) out of %d ",
                     "(%.1f%%). This usually means a predictor was NA at some ",
                     "time step -- check that every lagged/derived column used ",
                     "in a model is initialised by init_recode and updated by ",
                     "in_recode/out_recode."),
              n_na, nrow(data), 100 * n_na / nrow(data))
    )
  }

  result <- if (return_data) {
    data
  } else {
    # Mean of analytic Pred_Y over all N individuals (no row removal).
    sum(data[["Pred_Y"]]) / length(data[["Pred_Y"]])
  }

  if (collect_pool) return(list(estimate = result, pool = m_pool_out))
  return(result)
}
