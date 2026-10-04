
#' Random data simulation from predicted value.
#'
#' @description
#'  Internal use only, predict response and simulate random data. For numeric
#'  values the simulated value is restricted to the observed value range unless
#'  the model was created with \code{spec_model(truncate = FALSE)}. Binary
#'  values are drawn from the fitted probability as is, without clamping it
#'  away from 0 and 1, as in \pkg{gfoRmula}. Categorical values are returned in
#'  the variable's type in the data: a factor keeps its levels, numeric codes
#'  stay numeric.
#'
#' @param model fitted objects defined in the `spec_model`.
#' @param newdt a data frame in which to look for variables with which to predict.
#'
#' @return A simulated random vector using the predicted value from model and newdt.
#'
#' @keywords internal
#'
sim_value <- function(model, newdt) {
  var_type <- model$var_type

  # Clip simulated numeric values to the observed range of the response unless
  # the model was created with spec_model(truncate = FALSE). NULL means the
  # model object predates the `truncate` argument: keep the historical default.
  do_trunc <- is.null(model$truncate) || isTRUE(model$truncate)

  # Custom simulation function takes priority
  if (!is.null(model$custom_sim)) {
    out <- model$custom_sim(model$fitted, newdt)
    if (do_trunc && is.numeric(model$val_ran)) {
      out <- pmin(pmax(out, model$val_ran[1L]), model$val_ran[2L])
    }
    return(out)
  }

  # Categorical: multinom returns a probability matrix — predict() is required
  if (var_type == "categorical") {
    if (!requireNamespace("Hmisc", quietly = TRUE)) {
      stop("Package 'Hmisc' is required for var_type = \"categorical\". ",
           "Please install it.", call. = FALSE)
    }
    pred <- predict(model$fitted, newdata = newdt, type = "probs")
    # predict.multinom() returns a vector, not a matrix, for a two-level
    # response (the second level's probability) and for a single row.
    if (is.null(dim(pred))) {
      lev  <- model$fitted$lev
      pred <- if (length(lev) == 2L) cbind(1 - pred, pred) else matrix(pred, nrow = 1L)
      colnames(pred) <- lev
    }
    # rMultinom() returns the level labels as text. Returned as text, a
    # numeric-coded variable no longer matched the numeric predictor its
    # dependent models were fitted on.
    return(as_observed_type(Hmisc::rMultinom(pred, 1L)[, 1L], model$rsp_proto))
  }

  # Linear predictor for binary and normal types, from the terms and
  # coefficients pre-extracted after fitting (fit_spec_models()): the design
  # matrix is built as predict() builds it, then multiplied directly.
  lp <- lin_pred(model, newdt)

  if (var_type == "binary") {
    # The fitted probability is used as is, as gfoRmula draws
    # rbinom(n, 1, predict(type = "response")). A floor of 1e-5 had added
    # events wherever the fitted daily probability was smaller.
    return(rbinom(nrow(newdt), 1L, model$linkinv(lp)))
  } else {
    # Continuous/normal: lp + N(0, sigma) avoids a full predict() call.
    # model$sigma is the residual SD pre-extracted after fitting; it is NULL
    # when sigma() failed on the fitted object (e.g. a custom_fit class with
    # no sigma method and no custom_sim), where `NULL < eps` would give the
    # opaque "argument is of length zero" error.
    if (is.null(model$sigma)) {
      stop("No residual SD could be extracted from the fitted model for '",
           model$rsp_vars, "' (sigma() failed), so normal-type simulation ",
           "is not possible. Supply a custom_sim function in spec_model().",
           call. = FALSE)
    }
    if (model$sigma < .Machine$double.eps) {
      causalmed_env$warning <- c(
        causalmed_env$warning,
        sprintf("Zero residual variance in model for '%s'. Simulated values will equal the predicted mean.",
                model$rsp_vars)
      )
    }
    out <- lp + rnorm(nrow(newdt), 0, model$sigma)
    if (do_trunc) {
      out <- pmin(pmax(out, model$val_ran[1L]), model$val_ran[2L])
    }
    return(out)
  }
}

# Convert categorical draws (level labels, as text) to the type of `proto`, a
# zero-length vector of the variable as it is in the data: a factor, ordered or
# not, keeps its levels in their order; numeric codes become numbers again.
# `proto` is NULL for a model object built before it was recorded, and the
# labels are then returned as text, as before.
as_observed_type <- function(x, proto) {
  if (is.factor(proto))  return(factor(x, levels = levels(proto),
                                       ordered = is.ordered(proto)))
  if (is.integer(proto)) return(as.integer(x))
  if (is.numeric(proto)) return(as.numeric(x))
  if (is.logical(proto)) return(as.logical(x))
  x
}
