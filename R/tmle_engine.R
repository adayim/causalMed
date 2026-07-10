# TMLE engine for longitudinal natural direct/indirect effects
# (mediation_type = "N"), implementing Zheng & van der Laan (2017,
# J. Causal Inference 5(2), doi:10.1515/jci-2016-0006), Section 4.3.
#
# Notation mapping (paper -> causalMed):
#   A_t = (AC_t, AE_t)  -> censor model (intervened to "uncensored") + exposure
#   R_t                 -> covariate models listed BEFORE the mediator model
#   Z_t                 -> the (single) mediator model
#   L_t                 -> covariate models listed AFTER the mediator model,
#                          plus the outcome/survival node
#   Y                   -> response of the outcome/survival model
#
# For each functional J^{a,a'} = E[Y_{a, Gamma_{a'}}] the estimator runs the
# backward targeted sequential-regression loop of the paper's Section 4.3:
# three regressions per time point (L-step, Z-step, R-step), each followed by
# a one-dimensional logistic fluctuation with offset logit(Qbar) and weights
# given by the clever covariates H_{L_t}, H_{Z_t}, H_{R_t} (paper eqs. 11-13).
# The estimate is the mean of the targeted Qbar_{R_1} over subjects, and the
# influence-curve variance gives Wald CIs (paper Theorem 1 / Corollary 1).
#
# All regressions are quasibinomial GLMs on the [0,1]-bounded outcome scale;
# nuisance densities come from the same spec_model() GLM fits used by the
# Monte Carlo engine (shared fit_spec_models()).

# Bound probabilities away from 0/1 for logit transforms.
.tmle_bound <- function(p, eps = 1e-6) pmin(pmax(p, eps), 1 - eps)

# Evaluate the conditional density/pmf of a fitted spec_model at the observed
# response values in `newdt`. Sibling of sim_value(): sim_value() draws,
# dens_value() evaluates. Uses the pre-extracted Xterms/beta/linkinv/sigma.
dens_value <- function(model, newdt) {
  y <- newdt[[model$rsp_vars]]

  if (!is.null(model$custom_sim)) {
    stop("var_type = 'custom' is not supported by estimator = 'tmle'.",
         domain = "causalMed")
  }

  if (model$var_type == "categorical") {
    pred <- predict(model$fitted, newdata = newdt, type = "probs")
    if (is.null(dim(pred))) {   # two-level multinom returns a vector
      lev <- model$fitted$lev
      pred <- cbind(1 - pred, pred)
      colnames(pred) <- lev
    }
    return(.tmle_bound(pred[cbind(seq_len(nrow(pred)), as.character(y))]))
  }

  mm <- model.matrix(model$Xterms, data = newdt)
  lp <- drop(mm %*% model$beta)

  if (model$var_type == "binary") {
    p <- .tmle_bound(model$linkinv(lp))
    return(y * p + (1 - y) * (1 - p))
  }

  # normal
  pmax(dnorm(y, mean = lp, sd = model$sigma), 1e-100)
}

# Extract lag-column mappings (lag_col -> source_col) from an in_recode
# specification: entries of the form  recodes(M_lag1 = M)  where the RHS is a
# bare column name. Used to (i) identify exposure-history columns that must be
# set alongside the exposure when evaluating regime-specific predictions, and
# (ii) map lag columns to current values when a time-t fit is evaluated on
# time-(t-1) rows.
.tmle_lag_map <- function(in_recode) {
  if (is.null(in_recode)) return(list())
  out <- list()
  for (i in seq_along(in_recode)) {
    ex <- in_recode[[i]]
    if (is.name(ex)) out[[names(in_recode)[i]]] <- as.character(ex)
  }
  out
}

# Constant initial values for lag columns from init_recode (e.g. M_lag1 = 0).
.tmle_init_vals <- function(init_recode) {
  if (is.null(init_recode)) return(list())
  out <- list()
  for (i in seq_along(init_recode)) {
    ex <- init_recode[[i]]
    if (is.numeric(ex) || is.logical(ex))
      out[[names(init_recode)[i]]] <- as.numeric(ex)
  }
  out
}

# Set the exposure (and its lag columns) to regime value `x` on a copy of
# `rows`. At the first time point exposure-lag columns keep their
# init_recode value (no pre-baseline exposure under the regime).
.tmle_set_regime <- function(rows, exposure, exp_lags, x, is_first,
                             init_vals) {
  nd <- data.table::copy(rows)
  data.table::set(nd, j = exposure, value = x)
  for (lg in exp_lags) {
    if (is_first) {
      if (!is.null(init_vals[[lg]]))
        data.table::set(nd, j = lg, value = init_vals[[lg]])
      # else: keep observed value at the first time point
    } else {
      data.table::set(nd, j = lg, value = x)
    }
  }
  nd
}

# Quasibinomial working regression: fit `resp ~ predictors` on the rows where
# fit_ok is TRUE, and return bounded predictions on `pred_dt` (which may be a
# regime-modified copy). Falls back to the mean of resp if the design is
# degenerate. Warnings are routed through the package warning collector.
.tmle_qreg <- function(resp, fit_dt, pred_dt, predictors, label) {
  fit_ok <- is.finite(resp)
  preds  <- predictors[predictors %in% names(fit_dt)]
  # Drop constant columns (single observed level) to avoid rank issues.
  if (length(preds) > 0L) {
    keep <- vapply(preds, function(v) {
      x <- fit_dt[[v]][fit_ok]
      length(unique(x[!is.na(x)])) > 1L
    }, logical(1))
    preds <- preds[keep]
  }
  if (sum(fit_ok) == 0L) {
    return(rep(NA_real_, nrow(pred_dt)))
  }
  if (length(preds) == 0L) {
    return(rep(.tmle_bound(mean(resp[fit_ok])), nrow(pred_dt)))
  }
  dfit <- as.data.frame(fit_dt[fit_ok, preds, with = FALSE])
  dfit$..resp <- .tmle_bound(resp[fit_ok])
  fml <- reformulate(preds, response = "..resp")
  fit <- run_withwarning_collect(
    glm(fml, data = dfit, family = quasibinomial()),
    msg = sprintf("TMLE sequential regression (%s)", label)
  )
  .tmle_bound(as.numeric(
    predict(fit, newdata = as.data.frame(pred_dt), type = "response")
  ))
}

# One-dimensional logistic fluctuation (paper eq. 20): intercept-only
# quasibinomial GLM of `resp` with offset logit(qbar) and weights H, applied
# to rows where H > 0 and resp is defined. Returns the epsilon estimate
# (0 when no rows contribute or the fit fails).
#
# Positivity guards: when fewer than `min_obs` observations carry positive
# weight the regime has essentially no support at this node — the fluctuation
# is skipped (epsilon = 0, i.e. the estimate falls back to the
# sequential-regression extrapolation for this node) and a warning is
# collected. Epsilon is also bounded to +/- `eps_bound` to prevent a single
# sparse stratum (e.g. all-zero outcomes among a handful of subjects) from
# collapsing the targeted estimate globally.
.tmle_fluctuate <- function(resp, qbar, H, label,
                            min_obs = 10L, eps_bound = 4) {
  use <- is.finite(resp) & is.finite(H) & H > 0
  if (sum(use) == 0L) return(0)
  if (sum(use) < min_obs) {
    causalmed_env$warning <- c(
      causalmed_env$warning,
      sprintf(paste0("TMLE %s: only %d observation(s) follow the intervened ",
                     "regime at this node (practical positivity violation). ",
                     "Fluctuation skipped; the estimate relies on model ",
                     "extrapolation here."),
              label, sum(use))
    )
    return(0)
  }
  eps <- tryCatch({
    dfl <- data.frame(y = .tmle_bound(resp[use]),
                      off = qlogis(.tmle_bound(qbar[use])),
                      w = H[use])
    fit <- run_withwarning_collect(
      glm(y ~ 1 + offset(off), data = dfl, weights = w,
          family = quasibinomial()),
      msg = sprintf("TMLE fluctuation (%s)", label)
    )
    unname(coef(fit)[1L])
  }, error = function(e) 0)
  if (!is.finite(eps)) eps <- 0
  if (abs(eps) > eps_bound) {
    causalmed_env$warning <- c(
      causalmed_env$warning,
      sprintf(paste0("TMLE %s: fluctuation epsilon (%.2f) exceeded the ",
                     "stability bound and was truncated to %s%s (sparse ",
                     "support for the intervened regime)."),
              label, eps, if (eps > 0) "+" else "-", eps_bound)
    )
    eps <- sign(eps) * eps_bound
  }
  eps
}

# Truncate a clever-covariate vector at its `trunc` quantile (computed over
# strictly positive weights). Positivity protection; documented in mediation().
.tmle_trunc <- function(H, trunc) {
  pos <- H[is.finite(H) & H > 0]
  if (length(pos) == 0L || is.null(trunc) || trunc >= 1) return(H)
  cap <- quantile(pos, trunc, na.rm = TRUE, names = FALSE)
  pmin(H, cap)
}

# ---------------------------------------------------------------------------
# Main engine: TMLE of J^{a,a'} for all four (a, a') combinations.
#
# Returns list(psi = named numeric(4), eic = n x 4 matrix (subject-level
# influence-curve values, columns Phi00/Phi10/Phi01/Phi11), n = n subjects).
# ---------------------------------------------------------------------------
tmle_natural_mediation <- function(data,
                                   id_var,
                                   base_vars,
                                   exposure,
                                   outcome,
                                   time_var,
                                   models,
                                   in_recode = NULL,
                                   init_recode = NULL,
                                   weight_trunc = 0.995) {

  data <- data.table::as.data.table(data)
  data.table::setorderv(data, c(id_var, time_var))

  # Apply fit-time model recodes to the working copy once, in list order,
  # so density evaluation sees the same columns the models were fit on.
  for (m in models) {
    if (!is.null(m$recode)) data <- apply_recodes(data, m$recode)
  }

  # init_recode defines the deterministic first-time values of derived
  # columns (e.g. lag columns that are NA at baseline in the observed data).
  # Apply it to the first-time rows so density evaluation and the sequential
  # regressions see the same baseline values the simulation engine would use.
  if (!is.null(init_recode)) {
    t1  <- min(na.omit(data[[time_var]]))
    ix1 <- which(data[[time_var]] == t1)
    sub1 <- data[ix1]
    apply_recodes(sub1, init_recode)
    for (nm in names(init_recode)) {
      if (!nm %in% names(data)) {
        data.table::set(data, j = nm, value = NA_real_)
      }
      data.table::set(data, i = ix1, j = nm, value = sub1[[nm]])
    }
  }

  # ---- Nuisance model fitting (censoring-aware) ---------------------------
  # The position of the censor model in `models` encodes the within-interval
  # ordering. Nodes listed AFTER the censoring indicator are only observed
  # while the subject remains uncensored, so their models must be fit
  # "among observations that remained uncensored at time t" (paper Section
  # 4.3) — on a subject's censoring row those values are unobserved (or
  # carry data-format default values that would otherwise contaminate the
  # fits and miscalibrate the clever-covariate weights). Models listed at or
  # before the censoring node (including the censor model itself) are fit on
  # all at-risk rows.
  raw_types <- vapply(models, function(m) m$mod_type, character(1))
  cen_pos   <- which(raw_types == "censor")
  if (length(cen_pos) == 1L) {
    cvar    <- all.vars(formula(models[[cen_pos]]$call)[[2]])
    dat_unc <- data[as.numeric(data[[cvar]]) == 0]
    fit_mods <- vector("list", length(models))
    for (k in seq_along(models)) {
      fit_dat <- if (k > cen_pos) dat_unc else data
      fit_mods[[k]] <- fit_spec_models(models[k], fit_dat)[[1L]]
    }
    names(fit_mods) <- names(models)
  } else {
    fit_mods <- fit_spec_models(models, data)
  }

  mod_types <- vapply(fit_mods, function(m) m$mod_type, character(1))
  med_idx   <- which(mod_types == "mediator")
  out_idx   <- which(mod_types %in% c("outcome", "survival"))
  exp_idx   <- which(mod_types == "exposure")
  cen_idx   <- which(mod_types == "censor")
  cov_idx   <- which(mod_types == "covariate")

  med_mod <- fit_mods[[med_idx]]
  out_mod <- fit_mods[[out_idx]]
  exp_mod <- if (length(exp_idx)) fit_mods[[exp_idx]] else NULL
  cen_mod <- if (length(cen_idx)) fit_mods[[cen_idx]] else NULL

  is_survival <- fit_mods[[out_idx]]$mod_type == "survival" || length(cen_idx) > 0

  med_var    <- med_mod$rsp_vars
  censor_var <- if (!is.null(cen_mod)) cen_mod$rsp_vars else NULL

  R_mods <- fit_mods[cov_idx[cov_idx < med_idx]]
  L_mods <- fit_mods[cov_idx[cov_idx > med_idx]]
  R_vars <- vapply(R_mods, function(m) m$rsp_vars, character(1))
  L_vars <- vapply(L_mods, function(m) m$rsp_vars, character(1))

  # Lag bookkeeping
  lag_map   <- .tmle_lag_map(in_recode)
  init_vals <- .tmle_init_vals(init_recode)
  exp_lags  <- names(lag_map)[vapply(lag_map, identical, logical(1), exposure)]

  # Predictor sets for the sequential regressions. Conditioning follows the
  # paper: the L-step sees (A_t, R_t, Z_t, history); the Z-step additionally
  # drops the current mediator; the R-step additionally drops the current
  # R-block covariates (history enters through lag columns).
  rhs_all <- unique(unlist(lapply(fit_mods, function(m)
    setdiff(all.vars(formula(m$fitted)), m$rsp_vars))))
  avail   <- intersect(unique(c(rhs_all, base_vars, time_var, exposure)),
                       names(data))
  pred_L <- setdiff(avail, c(L_vars, outcome, censor_var))
  pred_Z <- setdiff(pred_L, med_var)
  pred_R <- setdiff(pred_Z, R_vars)

  # Subject/time indexing
  ids      <- unique(data[[id_var]])
  n        <- length(ids)
  id_pos   <- match(data[[id_var]], ids)           # row -> subject index
  time_seq <- sort(unique(na.omit(data[[time_var]])))
  Tn       <- length(time_seq)
  t_pos    <- match(data[[time_var]], time_seq)    # row -> time index

  # Per-row nuisance quantities under regimes 0 and 1 -----------------------
  eval_dens_at <- function(model, x, tix) {
    # density of `model`'s observed response with exposure regime set to x,
    # evaluated on the rows with time index tix (logical or integer index).
    rows <- data[tix, , drop = FALSE]
    is_first <- all(rows[[time_var]] == time_seq[1L])
    nd <- .tmle_set_regime(rows, exposure, exp_lags, x, is_first, init_vals)
    dens_value(model, nd)
  }

  # Exposure propensity P(A_t = 1 | past) at observed history.
  if (!is.null(exp_mod)) {
    mmA <- model.matrix(exp_mod$Xterms, data = data)
    pA1 <- .tmle_bound(exp_mod$linkinv(drop(mmA %*% exp_mod$beta)))
  } else {
    stop("estimator = 'tmle' requires an exposure model (mod_type = 'exposure').",
         domain = "causalMed")
  }

  # Probability of remaining uncensored at each row.
  if (!is.null(cen_mod)) {
    mmC     <- model.matrix(cen_mod$Xterms, data = data)
    pUncens <- .tmle_bound(1 - cen_mod$linkinv(drop(mmC %*% cen_mod$beta)))
  } else {
    pUncens <- rep(1, nrow(data))
  }

  # Densities at regime 0 and 1 for mediator, R-block, L-block and (survival
  # only) the outcome hazard node, evaluated per time slice so first-time lag
  # handling is correct.
  dens_by_regime <- function(model) {
    d0 <- numeric(nrow(data)); d1 <- numeric(nrow(data))
    for (ti in seq_len(Tn)) {
      tix <- which(t_pos == ti)
      if (length(tix) == 0L) next
      d0[tix] <- eval_dens_at(model, 0, tix)
      d1[tix] <- eval_dens_at(model, 1, tix)
    }
    list(d0 = d0, d1 = d1)
  }

  dZ <- dens_by_regime(med_mod)
  dR_list <- lapply(R_mods, dens_by_regime)
  # The survival/outcome node is part of the L-block history: for rows where
  # the subject survives the interval the contribution is P(Y_t = 0 | ...) =
  # 1 - h_t. Event rows never enter the history ratios (contributions beyond
  # the event are zero), so using the survival density at the observed Y is
  # correct for the rows that are consumed.
  dL_list <- lapply(L_mods, dens_by_regime)
  if (is_survival) dL_list <- c(dL_list, list(dens_by_regime(out_mod)))

  # Observed outcome/event per row (0/1 scale enforced by the caller).
  y_obs <- as.numeric(data[[outcome]])

  # Per-subject cumulative products over time (data sorted id, time).
  cum_by_id <- function(x) {
    # cumulative product within subject
    ave(x, id_pos, FUN = cumprod)
  }
  cumall_by_id <- function(x) {
    # cumulative "all TRUE so far" within subject, as 0/1
    ave(as.numeric(x), id_pos, FUN = cumprod)
  }
  lag_by_id <- function(x, fill = 1) {
    # value at t-1 within subject (fill at first row)
    ave(x, id_pos, FUN = function(v) c(fill, v[-length(v)]))
  }

  a_obs <- as.numeric(data[[exposure]])

  # "Remained uncensored through t" indicator: a subject's own censoring row
  # (censor == 1) and anything after it must not contribute to fluctuations
  # or fits — the paper's I(AC_t = 1) component of the intervention-node
  # indicator.
  if (!is.null(censor_var)) {
    uncens <- cumall_by_id(as.numeric(data[[censor_var]]) == 0)
  } else {
    uncens <- rep(1, nrow(data))
  }

  # Cumulative pieces reused across functionals
  cumA <- list(
    `0` = cum_by_id(.tmle_bound(1 - pA1) * pUncens),
    `1` = cum_by_id(pA1 * pUncens)
  )
  indA <- list(
    `0` = cumall_by_id(a_obs == 0) * uncens,
    `1` = cumall_by_id(a_obs == 1) * uncens
  )
  # ratio prod_j dZ(a')/dZ(a):
  ratioZ <- function(a, ap) {
    if (a == ap) return(rep(1, nrow(data)))
    num <- if (ap == 1) dZ$d1 else dZ$d0
    den <- if (a  == 1) dZ$d1 else dZ$d0
    cum_by_id(num / den)
  }
  # ratio prod_j dX(a)/dX(a') over a list of node densities:
  ratioX <- function(dlist, a, ap) {
    if (a == ap || length(dlist) == 0L) return(rep(1, nrow(data)))
    out <- rep(1, nrow(data))
    for (d in dlist) {
      num <- if (a  == 1) d$d1 else d$d0
      den <- if (ap == 1) d$d1 else d$d0
      out <- out * (num / den)
    }
    cum_by_id(out)
  }

  # ---- run one functional J^{a, a'} ---------------------------------------
  run_one <- function(a, ap) {
    a_chr <- as.character(a); ap_chr <- as.character(ap)

    rZ  <- ratioZ(a, ap)                     # prod_{j<=t} dZ(a')/dZ(a)
    rL  <- ratioX(dL_list, a, ap)            # prod_{j<=t} dL(a)/dL(a')
    rR  <- ratioX(dR_list, a, ap)            # prod_{j<=t} dR(a)/dR(a')

    H_L <- .tmle_trunc(indA[[a_chr]]  / cumA[[a_chr]]  * rZ, weight_trunc)
    H_Z <- .tmle_trunc(indA[[ap_chr]] / cumA[[ap_chr]] *
                         lag_by_id(rL) * rR, weight_trunc)
    H_R <- .tmle_trunc(indA[[a_chr]]  / cumA[[a_chr]]  *
                         lag_by_id(rZ), weight_trunc)

    eic  <- numeric(n)          # subject-level influence-curve accumulator
    resp <- rep(NA_real_, n)    # Qbar_{R_{t+1}} carried backward, per subject

    for (ti in rev(seq_len(Tn))) {
      tix   <- which(t_pos == ti)
      rows  <- data[tix, , drop = FALSE]
      sid   <- id_pos[tix]
      is_first <- (ti == 1L)

      # Response for the L-step: Qbar_{R_{t+1}} per subject present at t.
      if (ti == Tn) {
        r_t <- y_obs[tix]
      } else {
        r_t <- resp[sid]
        if (is_survival) {
          # Subjects whose event occurs at t are dead by t+1: Qbar = 1.
          r_t[y_obs[tix] == 1] <- 1
        }
      }
      # A subject censored at t contributes no observed response at t: the
      # regression excludes them (NA), while they still receive predictions.
      if (!is.null(censor_var)) {
        r_t[uncens[tix] == 0] <- NA_real_
      }

      # (a) L-step -----------------------------------------------------------
      nd_a  <- .tmle_set_regime(rows, exposure, exp_lags, a, is_first, init_vals)
      qL    <- .tmle_qreg(r_t, rows, nd_a, pred_L,
                          label = sprintf("L-step t=%s", time_seq[ti]))
      epsL  <- .tmle_fluctuate(r_t, qL, H_L[tix],
                               label = sprintf("L-step t=%s", time_seq[ti]))
      qLs   <- .tmle_bound(plogis(qlogis(qL) + epsL))
      ok    <- is.finite(r_t)
      eic[sid[ok]] <- eic[sid[ok]] + H_L[tix][ok] * (r_t[ok] - qLs[ok])

      # (b) Z-step -----------------------------------------------------------
      nd_ap <- .tmle_set_regime(rows, exposure, exp_lags, ap, is_first, init_vals)
      qZ    <- .tmle_qreg(qLs, rows, nd_ap, pred_Z,
                          label = sprintf("Z-step t=%s", time_seq[ti]))
      epsZ  <- .tmle_fluctuate(qLs, qZ, H_Z[tix],
                               label = sprintf("Z-step t=%s", time_seq[ti]))
      qZs   <- .tmle_bound(plogis(qlogis(qZ) + epsZ))
      eic[sid] <- eic[sid] + H_Z[tix] * (qLs - qZs)

      # (c) R-step -----------------------------------------------------------
      # Fit on time-t rows (lag columns = time-(t-1) history), evaluate on
      # the same rows for the EIC term, and on time-(t-1) rows (lag -> current
      # mapping) to produce the next backward response.
      qR    <- .tmle_qreg(qZs, rows, nd_a, pred_R,
                          label = sprintf("R-step t=%s", time_seq[ti]))
      epsR  <- .tmle_fluctuate(qZs, qR, H_R[tix],
                               label = sprintf("R-step t=%s", time_seq[ti]))
      qRs   <- .tmle_bound(plogis(qlogis(qR) + epsR))
      eic[sid] <- eic[sid] + H_R[tix] * (qZs - qRs)

      if (ti == 1L) {
        # Qbar_{R_1}(L_0): everyone has a first-time row.
        psi_i <- rep(NA_real_, n)
        psi_i[sid] <- qRs
        psi <- mean(psi_i, na.rm = TRUE)
        eic <- eic + ifelse(is.na(psi_i), 0, psi_i - psi)
        return(list(psi = psi, eic = eic))
      }

      # Propagate: evaluate the (fluctuated) time-t R-step on time-(t-1)
      # rows, mapping lag columns to their current-value sources and setting
      # the exposure history to regime a.
      pix   <- which(t_pos == (ti - 1L))
      prow  <- data[pix, , drop = FALSE]
      ndp   <- data.table::copy(prow)
      for (lg in names(lag_map)) {
        src <- lag_map[[lg]]
        if (lg %in% names(ndp) && src %in% names(ndp))
          data.table::set(ndp, j = lg, value = ndp[[src]])
      }
      # The fit is a time-t function: keep the time value at t.
      data.table::set(ndp, j = time_var, value = time_seq[ti])
      ndp <- .tmle_set_regime(ndp, exposure, exp_lags, a,
                              is_first = (ti - 1L == 1L), init_vals)
      qR_prev <- .tmle_qreg(qZs, rows, ndp, pred_R,
                            label = sprintf("R-step eval t=%s", time_seq[ti]))
      qR_prev <- .tmle_bound(plogis(qlogis(qR_prev) + epsR))

      resp_new <- rep(NA_real_, n)
      resp_new[id_pos[pix]] <- qR_prev
      resp <- resp_new
    }
  }

  combos <- list(Phi00 = c(0, 0), Phi10 = c(1, 0),
                 Phi01 = c(0, 1), Phi11 = c(1, 1))
  res <- lapply(combos, function(cc) run_one(cc[1], cc[2]))

  psi <- vapply(res, function(r) r$psi, numeric(1))
  eic <- do.call(cbind, lapply(res, function(r) r$eic))
  colnames(eic) <- names(combos)

  list(psi = psi, eic = eic, n = n, fit_mods = fit_mods)
}
